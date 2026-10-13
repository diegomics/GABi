#!/usr/bin/env python3

# by Diego De Panis, 2026
# This script is part of the GABi pipeline
# note: AI tools may have been used to improve, clean and/or comment this version of the code

"""
Validate a pre-aligned BAM offered as pipeline input.

GABi can now be entered from reads or from an already-aligned BAM.
The technology is DECLARED in the data table's read_type column exactly the same
as for FASTQ, we do nothing with that here.
It reads the BAM header only to CHECK the declaration and
to report what the file looks like - the same relationship D06 has with a
declared sample_sex_karyo.

Input is the text of `samtools view -H`, not the BAM itself (to avoid pysam)

Exits non-zero on anything that would produce wrong results downstream. Things
that are merely worth knowing are reported and do not stop the run.
"""

import argparse
import re
import sys

# ---------------------------------------------------------------------------
# Aligner -> read length class. Used only to contradict a declaration, never to
# supply one: a BAM whose aligner we cannot identify is accepted on trust.
# These lists needs to be trimmed!
# ---------------------------------------------------------------------------
SHORT_READ_ALIGNERS = {
    "bwa", "bwa-mem2", "bwamem2", "bowtie2", "bowtie", "novoalign",
    "snap", "snap-aligner", "dragen-os", "hisat2", "subread", "stampy",
}
LONG_READ_ALIGNERS = {"minimap2", "winnowmap", "ngmlr", "pbmm2", "lra"}

# minimap2 -x/-ax presets -> technology
MINIMAP2_PRESETS = {
    "map-ont": "ont", "lr:hq": "ont", "map-ont-hq": "ont",
    "map-hifi": "hifi", "map-pb": "hifi", "map-pacbio": "hifi",
    "sr": "short", "sr-pe": "short",
}

# @RG PL values -> technology (SAM spec controlled vocabulary)
PL_TO_TECH = {
    "ILLUMINA": "short", "DNBSEQ": "short", "ELEMENT": "short",
    "SINGULAR": "short", "ULTIMA": "short", "BGI": "short", "MGI": "short",
    "ONT": "ont", "OXFORD_NANOPORE": "ont", "NANOPORE": "ont",
    "PACBIO": "pacbio",          # covers both HiFi and CLR - see below
}

SHORT_TECHS = {"illumina", "10x"}

# Tools whose presence in @PG means duplicates have been marked
MARKDUP_TOOLS = ("markdup", "markduplicates", "picard", "sambamba", "biobambam",
                 "bamsormadup", "dedup")


def parse_header(text):
    """SAM header text -> {'HD': [dict], 'SQ': [...], 'RG': [...], 'PG': [...]}"""
    out = {"HD": [], "SQ": [], "RG": [], "PG": []}
    for line in text.splitlines():
        if not line.startswith("@"):
            continue
        parts = line.rstrip("\n").split("\t")
        tag = parts[0][1:]
        if tag not in out:
            continue
        rec = {}
        for field in parts[1:]:
            if ":" in field:
                k, v = field.split(":", 1)
                rec[k] = v
        out[tag].append(rec)
    return out


def read_fai(path):
    """[(name, length), ...] in reference order."""
    out = []
    with open(path) as fh:
        for line in fh:
            t = line.rstrip("\n").split("\t")
            if len(t) >= 2:
                out.append((t[0], int(t[1])))
    return out


def detect_from_pg(pg_records):
    """(tech_class, evidence) from @PG. tech_class is short/ont/hifi/None."""
    for rec in pg_records:
        name = (rec.get("PN") or "").strip().lower()
        cl = (rec.get("CL") or "").strip()
        cl_first = cl.split()[0].rsplit("/", 1)[-1].lower() if cl else ""
        prog = name or cl_first

        if prog in LONG_READ_ALIGNERS or cl_first in LONG_READ_ALIGNERS:
            if "pbmm2" in (prog, cl_first):
                return "hifi", f"@PG pbmm2 ({cl[:60] or prog})"
            # minimap2 and friends: the preset decides
            m = re.search(r"-a?x\s+(\S+)", cl)
            if m:
                preset = m.group(1).lower()
                tech = MINIMAP2_PRESETS.get(preset)
                if tech:
                    return tech, f"@PG {prog} -x {preset}"
                return None, f"@PG {prog} with unrecognised preset '{preset}'"
            return None, f"@PG {prog} with no -x preset in CL"

        if prog in SHORT_READ_ALIGNERS or cl_first in SHORT_READ_ALIGNERS:
            return "short", f"@PG {prog}"

    return None, "no recognised aligner in @PG"


def detect_from_pl(rg_records):
    """(tech_class, evidence) from @RG PL."""
    pls = {(r.get("PL") or "").strip().upper() for r in rg_records}
    pls.discard("")
    if not pls:
        return None, "no @RG PL"
    if len(pls) > 1:
        return None, f"conflicting @RG PL values: {', '.join(sorted(pls))}"
    pl = pls.pop()
    return PL_TO_TECH.get(pl), f"@RG PL:{pl}"


def compatible(detected, declared):
    """Is a detected class consistent with the declared read_type?"""
    if detected is None:
        return True                      # nothing to contradict
    if detected == "short":
        return declared in SHORT_TECHS
    if detected == "pacbio":             # PACBIO covers HiFi and CLR
        return declared == "hifi"
    return detected == declared          # 'ont' / 'hifi'


def main():
    ap = argparse.ArgumentParser(description="Validate a pre-aligned input BAM")
    ap.add_argument("--header", required=True, help="output of `samtools view -H`")
    ap.add_argument("--fai", required=True, help="reference .fa.fai")
    ap.add_argument("--bam", required=True, help="BAM path, for messages")
    ap.add_argument("--sample-id", required=True)
    ap.add_argument("--read-type", required=True, help="declared technology")
    ap.add_argument("--report", required=True)
    args = ap.parse_args()

    declared = args.read_type.strip().lower()
    hdr = parse_header(open(args.header).read())
    errors, notes = [], []

    # --- sort order --------------------------------------------------------
    so = (hdr["HD"][0].get("SO") if hdr["HD"] else None) or "(absent)"
    if so != "coordinate":
        errors.append(
            f"BAM is not coordinate-sorted (@HD SO:{so}).\n"
            f"    Sort it first:  samtools sort -o sorted.bam {args.bam}"
        )
    else:
        notes.append("coordinate-sorted: yes")

    # --- reference concordance --------------------------------------------
    # The single most important check: a BAM aligned to a different reference
    # produces coordinates that silently mean something else.
    fai = read_fai(args.fai)
    bam_sq = {r["SN"]: int(r["LN"]) for r in hdr["SQ"] if "SN" in r and "LN" in r}
    missing, mismatched = [], []
    for name, length in fai:
        if name not in bam_sq:
            missing.append(name)
        elif bam_sq[name] != length:
            mismatched.append(f"{name} (reference {length}, BAM {bam_sq[name]})")

    if mismatched:
        errors.append(
            "BAM and reference disagree on contig lengths - these are different "
            "assemblies:\n      " + "\n      ".join(mismatched[:5])
            + (f"\n      ... and {len(mismatched)-5} more" if len(mismatched) > 5 else "")
        )
    if missing:
        errors.append(
            f"{len(missing)} reference contig(s) absent from the BAM header, "
            f"e.g. {', '.join(missing[:5])}.\n"
            f"    The BAM must be aligned to the same assembly given in asm_file."
        )
    if not missing and not mismatched:
        notes.append(f"reference concordance: all {len(fai)} contigs match")
    extra = [n for n in bam_sq if n not in dict(fai)]
    if extra:
        notes.append(
            f"BAM has {len(extra)} contig(s) not in the reference "
            f"(e.g. {', '.join(sorted(extra)[:3])}) - they will not be called"
        )

    # --- read groups and sample name --------------------------------------
    # GATK and bcftools take the VCF sample name from @RG SM. A mismatch would
    # put the wrong name in every output, including the joint call set.
    if not hdr["RG"]:
        errors.append(
            "BAM has no @RG read group.\n"
            "    GATK HaplotypeCaller refuses to run without one, and the VCF "
            "sample name comes from @RG SM."
        )
    else:
        sms = {(r.get("SM") or "").strip() for r in hdr["RG"]}
        sms.discard("")
        if not sms:
            errors.append("BAM read groups carry no SM (sample name).")
        elif len(sms) > 1:
            errors.append(
                f"BAM contains multiple sample names: {', '.join(sorted(sms))}.\n"
                f"    GABi expects one sample per BAM."
            )
        elif sms != {args.sample_id}:
            errors.append(
                f"@RG SM is '{sms.pop()}' but the data table calls this sample "
                f"'{args.sample_id}'.\n"
                f"    Callers take the VCF sample name from the BAM, so these "
                f"must agree. Fix with samtools addreplacerg, or rename in the "
                f"data table."
            )
        else:
            notes.append(f"sample name: {args.sample_id} (matches @RG SM)")

    # --- technology cross-check -------------------------------------------
    pg_tech, pg_why = detect_from_pg(hdr["PG"])
    pl_tech, pl_why = detect_from_pl(hdr["RG"])
    notes.append(f"aligner evidence: {pg_why}")
    notes.append(f"platform evidence: {pl_why}")

    for detected, why in ((pg_tech, pg_why), (pl_tech, pl_why)):
        if not compatible(detected, declared):
            errors.append(
                f"BAM looks like '{detected}' data ({why}) but the data table "
                f"declares read_type '{declared}'.\n"
                f"    One of the two is wrong; GABi will not guess which."
            )
    if pg_tech is None and pl_tech is None:
        notes.append(
            f"technology NOT verifiable from the header - trusting the declared "
            f"read_type '{declared}'"
        )

    # --- duplicates (report only) -----------------------------------------
    pg_blob = " ".join(
        f"{r.get('PN','')} {r.get('ID','')} {r.get('CL','')}" for r in hdr["PG"]
    ).lower()
    if any(t in pg_blob for t in MARKDUP_TOOLS):
        notes.append("duplicates: appear to be marked (@PG mentions a markdup tool)")
    else:
        notes.append(
            "duplicates: no markdup tool found in @PG - if they were never marked, "
            "duplicate reads will be counted as independent evidence"
        )
    notes.append(
        "DEDUP_PE / DEDUP_LONG do not apply to imported BAMs - they only gate "
        "the mapping-time markdup step, which does not run here"
    )

    # --- write the report --------------------------------------------------
    with open(args.report, "w") as fh:
        fh.write(f"# Imported BAM: {args.sample_id}\n\n")
        fh.write(f"- **File:** `{args.bam}`\n")
        fh.write(f"- **Declared read_type:** {declared}\n\n")
        fh.write("## Checks\n\n")
        for n in notes:
            fh.write(f"- {n}\n")
        if errors:
            fh.write("\n## Failures\n\n")
            for e in errors:
                fh.write(f"- {e}\n")

    for n in notes:
        print(f"[GABi]   {n}")

    if errors:
        print(f"\n[GABi] ERROR: {args.bam} cannot be used as input:", file=sys.stderr)
        for e in errors:
            print(f"[GABi]   - {e}", file=sys.stderr)
        sys.exit(1)

    print(f"[GABi] BAM accepted for {args.sample_id}")


if __name__ == "__main__":
    main()
