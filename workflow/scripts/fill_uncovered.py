#!/usr/bin/env python3

# by Diego De Panis, 2026
# This script is part of the GABi pipeline
# note: AI tools may have been used to improve, clean and/or comment this version of the code

"""
`bcftools mpileup` only produces a record where at least one read is present, so
uncovered positions are simply missing from the output. DeepVariant and GATK
gVCFs instead span the whole reference, representing uncovered regions as blocks
with DP=0. There is a gap there between the callers we use!
This script closes that gap so all three callers describe the same genome.

The idea is to write a "no data" gVCF-style record per ref uncovered interval:
    CHROM  START  .  <REF>  .  .  .  END=<end>  GT  ./.

In basepair mode the caller expands these intervals with
`bcftools convert --gvcf2vcf`, giving one record per position.

The output reuses the called file's exact header, so `bcftools concat` sees
two files with identical headers and identical samples.
"""

import argparse
import sys


def read_fai(path):
    """[(name, length, offset, linebases, linewidth), ...] in file order."""
    out = []
    with open(path) as fh:
        for line in fh:
            t = line.rstrip("\n").split("\t")
            if len(t) >= 5:
                out.append((t[0], int(t[1]), int(t[2]), int(t[3]), int(t[4])))
    return out


class Fasta:
    """Minimal faidx-backed random access, this avoids a pysam dependency."""

    def __init__(self, fasta, fai_records):
        self.fh = open(fasta, "rb")
        self.idx = {n: (ln, off, lb, lw) for n, ln, off, lb, lw in fai_records}

    def base(self, chrom, pos1):
        """1-based single base. Returns 'N' if anything looks off."""
        try:
            length, offset, linebases, linewidth = self.idx[chrom]
        except KeyError:
            return "N"
        if pos1 < 1 or pos1 > length:
            return "N"
        p = pos1 - 1
        byte = offset + (p // linebases) * linewidth + (p % linebases)
        self.fh.seek(byte)
        b = self.fh.read(1).decode("ascii", "replace").upper()
        return b if b in "ACGTN" else "N"

    def close(self):
        self.fh.close()


def covered_intervals(path):
    """Yield (chrom, start1, end1) from `bcftools query -f '%CHROM\\t%POS\\t%END\\n'`.

    %END is INFO/END when present (block mode) and POS otherwise (basepair), so
    the same stream describes both resolutions.
    """
    with open(path) as fh:
        for line in fh:
            t = line.rstrip("\n").split("\t")
            if len(t) < 3:
                continue
            chrom, start = t[0], int(t[1])
            try:
                end = int(t[2])
            except ValueError:
                end = start
            yield chrom, start, max(start, end)


def gaps_for_contig(length, intervals):
    """Complement of `intervals` (1-based inclusive, sorted) within 1..length."""
    cursor = 1
    for start, end in intervals:
        if start > cursor:
            yield cursor, start - 1
        cursor = max(cursor, end + 1)
    if cursor <= length:
        yield cursor, length


def main():
    ap = argparse.ArgumentParser(
        description="Emit no-data records for reference positions with no coverage"
    )
    ap.add_argument("--fai", required=True, help="reference .fa.fai")
    ap.add_argument("--fasta", required=True, help="reference .fa")
    ap.add_argument("--covered", required=True,
                    help="bcftools query -f '%%CHROM\\t%%POS\\t%%END\\n' output")
    ap.add_argument("--header", required=True,
                    help="header of the called file (bcftools view -h)")
    ap.add_argument("--output", required=True)
    ap.add_argument("--stats", help="optional summary written here")
    args = ap.parse_args()

    fai = read_fai(args.fai)
    lengths = {n: ln for n, ln, *_ in fai}
    order = [n for n, *_ in fai]

    # Group covered intervals by contig. Records arrive sorted within a contig,
    # so a running merge is enough - no need to hold the whole genome in memory
    # position by position, only interval endpoints.
    merged = {c: [] for c in order}
    unknown = set()
    for chrom, start, end in covered_intervals(args.covered):
        if chrom not in merged:
            unknown.add(chrom)
            continue
        bucket = merged[chrom]
        if bucket and start <= bucket[-1][1] + 1:
            bucket[-1] = (bucket[-1][0], max(bucket[-1][1], end))
        else:
            bucket.append((start, end))

    if unknown:
        print(f"WARNING: {len(unknown)} contig(s) in the call set are absent from "
              f"the .fai and were ignored: {', '.join(sorted(unknown)[:5])}",
              file=sys.stderr)

    fasta = Fasta(args.fasta, fai)
    n_intervals = 0
    n_positions = 0

    # Reuse the called file's header so `bcftools concat` sees identical headers
    # and identical samples. TO CHECK: in basepair mode the call set has no
    # reference blocks, so INFO/END is never declared - but the gap records below
    # use it, and an undeclared tag makes bcftools abort with
    # "Tag not defined in header".
    with open(args.header) as hfh:
        header_lines = hfh.readlines()

    if not any(l.startswith("##INFO=<ID=END,") for l in header_lines):
        end_decl = ('##INFO=<ID=END,Number=1,Type=Integer,Description='
                    '"End position of the variant described in this record">\n')
        for i, l in enumerate(header_lines):
            if l.startswith("#CHROM"):
                header_lines.insert(i, end_decl)
                break
        else:
            header_lines.append(end_decl)

    with open(args.output, "w") as out:
        out.writelines(header_lines)

        for chrom in order:
            length = lengths[chrom]
            for start, end in gaps_for_contig(length, merged[chrom]):
                ref = fasta.base(chrom, start)
                # INFO/END makes this a gVCF block, so `bcftools convert
                # --gvcf2vcf` can expand it to one record per position. That
                # expansion rewrites REF from the FASTA at every position, so
                # only this first base has to be right.
                out.write(
                    f"{chrom}\t{start}\t.\t{ref}\t.\t.\t.\tEND={end}\tGT\t./.\n"
                )
                n_intervals += 1
                n_positions += end - start + 1

    fasta.close()

    genome = sum(lengths.values())
    pct = (100.0 * n_positions / genome) if genome else 0.0
    msg = (f"uncovered: {n_intervals} interval(s), {n_positions} position(s) "
           f"of {genome} ({pct:.2f}% of the reference)")
    print(msg)
    if args.stats:
        with open(args.stats, "w") as fh:
            fh.write(msg + "\n")


if __name__ == "__main__":
    main()
