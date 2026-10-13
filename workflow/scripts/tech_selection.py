#!/usr/bin/env python3

# by Diego De Panis, 2026
# This script is part of the GABi pipeline
# note: AI tools may have been used to improve, clean and/or comment this version of the code

"""
Select ONLY ONE sequencing technology per sample for variant calling.

Technologies are never merged: exactly one is chosen,
using DATA_PRIORITY and per-technology coverage.

Decision rule
--------------
  1. Discard technologies below MIN_TECH_COV, and any the caller cannot handle
     (the eligible set is computed at DAG-build time and passed in via
     --eligible, so this script needs no knowledge of callers).
  2. If any remaining technology reaches GOOD_TECH_COV, DATA_PRIORITY decides
     between those.
  3. If none does, the technology with the highest coverage wins.

Coverage comes from MEAN_DEPTH in the per-technology QC summaries written by
D05. Note that F01 resolves its `auto` depth thresholds from WMD (length-
weighted median depth) in the same files: mean depth is the right statistic for
"which technology has more data", median for "what depth should I expect at a
site". The divergence is deliberate.
"""

import argparse
import json
import sys
import traceback
from pathlib import Path


def parse_priority(raw):
    """DATA_PRIORITY string -> ordered list of technology names."""
    import re
    return [x.strip().lower() for x in re.split(r"[>\s,]+", str(raw)) if x.strip()]


def read_mean_depth(summary_path, log):
    """Pull MEAN_DEPTH out of a D05 QC summary. Returns 0.0 if unreadable.

    Missing coverage is treated as zero rather than fatal: the technology then
    falls below MIN_TECH_COV and drops out of contention, which is the safe
    outcome. A sample whose technologies ALL fail to parse is reported loudly
    by select_tech() below.
    """
    try:
        with open(summary_path) as fh:
            for line in fh:
                if line.startswith("MEAN_DEPTH"):
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 2:
                        return float(parts[1])
        log(f"WARNING: no MEAN_DEPTH line in {summary_path}; treating as 0.0")
    except (OSError, ValueError) as e:
        log(f"WARNING: could not read {summary_path}: {e}; treating as 0.0")
    return 0.0


def select_tech(coverage, eligible, priority, good_cov, min_cov, log):
    """Pick one technology. Returns (chosen, reason, warning_or_None).

    `coverage` maps every technology the sample has -> mean depth.
    `eligible` is the subset the configured caller can actually use.
    """
    # --- step 1: eligible, and above the absolute floor ---------------------
    usable = {t: c for t, c in coverage.items() if t in eligible and c >= min_cov}

    excluded_by_caller = sorted(t for t in coverage if t not in eligible)
    excluded_by_cov = sorted(
        t for t in coverage
        if t in eligible and coverage[t] < min_cov
    )

    if not usable:
        # Nothing clears MIN_TECH_COV. Rather than abort the whole run for one
        # sample, fall back to the best eligible technology and flag it.
        pool = {t: c for t, c in coverage.items() if t in eligible}
        if not pool:
            raise SystemExit(
                f"ERROR: no technology available for this sample is usable with "
                f"the configured caller (have: {', '.join(sorted(coverage)) or 'none'}; "
                f"caller accepts: {', '.join(sorted(eligible)) or 'none'})"
            )
        chosen = max(pool, key=pool.get)
        warning = (
            f"no technology reaches MIN_TECH_COV={min_cov}x; "
            f"falling back to {chosen} at {pool[chosen]:.2f}x"
        )
        log(f"WARNING: {warning}")
        return chosen, "below MIN_TECH_COV, best available", warning

    # --- step 1b: no competition -------------------------------------------
    # With a single candidate there is nothing to rank, so we say that and byee.
    # Reporting "no technology reaches GOOD_TECH_COV; highest coverage wins"
    # for a sample that only ever had one technology reads like a degraded
    # fallback when in fact no choice was ever available.
    if len(usable) == 1:
        chosen = next(iter(usable))
        if len(coverage) == 1:
            reason = f"only technology available ({usable[chosen]:.2f}x)"
        else:
            reason = (
                f"only technology left after filtering ({usable[chosen]:.2f}x); "
                f"others dropped: "
                f"{', '.join(sorted(set(excluded_by_caller) | set(excluded_by_cov))) or 'none'}"
            )
        return chosen, reason, None

    # --- step 2: anything good enough on its own? priority decides ----------
    good = [t for t in usable if usable[t] >= good_cov]
    if good:
        for tech in priority:
            if tech in good:
                return (tech,
                        f"{usable[tech]:.2f}x >= GOOD_TECH_COV={good_cov}, "
                        f"highest DATA_PRIORITY among those",
                        None)
        # None of the good technologies appears in DATA_PRIORITY at all
        chosen = max(good, key=usable.get)
        return (chosen,
                f"{usable[chosen]:.2f}x >= GOOD_TECH_COV={good_cov}, "
                f"not listed in DATA_PRIORITY so highest coverage wins",
                None)

    # --- step 3: nothing is good enough; coverage decides -------------------
    chosen = max(usable, key=usable.get)
    return (chosen,
            f"no technology reaches GOOD_TECH_COV={good_cov}; "
            f"highest coverage wins at {usable[chosen]:.2f}x",
            None)


def main():
    ap = argparse.ArgumentParser(
        description="Select one sequencing technology per sample for calling"
    )
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--summaries", nargs="+", required=True,
                    metavar="TECH=PATH",
                    help="per-technology D05 QC summaries, e.g. hifi=/path/x.summary.txt")
    ap.add_argument("--eligible", required=True,
                    help="comma-separated technologies the configured caller can use")
    ap.add_argument("--priority", default="hifi>illumina>ont",
                    help="DATA_PRIORITY string")
    ap.add_argument("--good-cov", type=float, default=15.0)
    ap.add_argument("--min-cov", type=float, default=1.0)
    ap.add_argument("--output", required=True)
    ap.add_argument("--log", required=True)

    args = ap.parse_args()

    with open(args.log, "w") as log_fh:
        def log(msg):
            log_fh.write(f"[tech-select] {args.sample_id}: {msg}\n")
            log_fh.flush()

        try:
            # TECH=PATH pairs, so no filename parsing and no ordering assumptions
            coverage = {}
            for item in args.summaries:
                if "=" not in item:
                    raise SystemExit(
                        f"ERROR: --summaries expects TECH=PATH pairs, got: {item!r}"
                    )
                tech, path = item.split("=", 1)
                coverage[tech.strip().lower()] = round(read_mean_depth(path, log), 2)

            if not coverage:
                raise SystemExit("ERROR: no QC summaries provided")

            eligible = [t.strip().lower() for t in args.eligible.split(",") if t.strip()]
            priority = parse_priority(args.priority)

            log("coverage: " + ", ".join(
                f"{t.upper()}={c:.2f}x" for t, c in sorted(coverage.items())
            ))
            log(f"eligible for this caller: {', '.join(eligible) or 'none'}")

            chosen, reason, warning = select_tech(
                coverage, eligible, priority, args.good_cov, args.min_cov, log
            )

            decision = {
                "sample_id": args.sample_id,
                "chosen": chosen,
                "reason": reason,
                "coverage": coverage,
                "eligible": eligible,
                "excluded_by_caller": sorted(t for t in coverage if t not in eligible),
                "thresholds": {"good_cov": args.good_cov, "min_cov": args.min_cov},
                "priority": priority,
            }
            if warning:
                decision["warning"] = warning

            Path(args.output).parent.mkdir(parents=True, exist_ok=True)
            with open(args.output, "w") as fh:
                json.dump(decision, fh, indent=2)

            log(f"SELECTED {chosen} ({coverage[chosen]:.2f}x) - {reason}")
            print(f"[GABi] {args.sample_id}: selected {chosen} "
                  f"({coverage[chosen]:.2f}x) - {reason}")

        except SystemExit:
            raise
        except Exception as e:
            log(f"ERROR: {e}")
            log(traceback.format_exc())
            print(f"FATAL: {e}", file=sys.stderr)
            print(traceback.format_exc(), file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
