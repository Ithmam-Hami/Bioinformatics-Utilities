#!/usr/bin/env python3
# =====================================================================================================================
# nuclear_coverage_map.py  -  whole-genome coverage map (SVG), one panel per chromosome, straight from a BAM
# ---------------------------------------------------------------------------------------------------------------------
# PURPOSE
#   Quick visual check of how evenly reads cover a nuclear genome: uneven libraries, coverage sags, aneuploid or
#   missing chromosomes, and mapping problems stand out at a glance. Also writes the numbers behind the picture.
#
# WHAT IT DRAWS
#   Every reference sequence longer than 100 kb, end to end, split into ~1,050 screen columns. Per column:
#     blue line  = median read depth          pale band = 10th-90th percentile of depth
#     red tick   = column where more than half the bases have no reads
#   Only reads with MAPQ >= 20 and bases with quality >= 20 are counted. Each chromosome has its own y scale; the
#   number under its name is its typical (median) depth. An even map = a flat blue line near that number.
#
# OUTPUTS
#   <out.svg>   the map
#   <out.tsv>   per-column numbers: chrom, start, end, p10, median, p90, zero_fraction
#   stdout      one summary row per chromosome: length, median depth, y-scale maximum, mostly-uncovered columns
#
# REQUIREMENTS  Python 3 (standard library only) and samtools on PATH. The BAM must be coordinate-sorted and indexed.
#
# USAGE
#   python3 nuclear_coverage_map.py <bam> <out.svg> <out.tsv> "<title>" [contig_to_skip]
#   contig_to_skip: name of a sequence to leave out, e.g. the mitochondrial contig when the BAM was mapped to a
#   whole genome including mtDNA. Use NONE (or leave it out) to draw everything over 100 kb.
#
# SETTINGS TO CHANGE IN THE CODE
#   100000            minimum sequence length drawn (bp)
#   1050              number of screen columns per chromosome
#   "-Q 20 -q 20"     minimum mapping quality and base quality passed to samtools depth
# =====================================================================================================================
"""Render a full-chromosome nuclear depth overview directly from a BAM, for any species.
Every reference sequence longer than 100 kb is drawn; one named contig (e.g. the mitochondrial genome in a
whole-genome BAM) can be left out. Each chromosome is shown end to end at screen resolution: blue = median depth
per screen column, pale = 10th-90th percentile, red tick = screen column with >50% uncovered bases.
MAPQ >= 20, base quality >= 20.
Usage: python3 nuclear_coverage_map.py <bam> <out.svg> <out.tsv> <title> [contig_to_skip]
"""
import collections
import math
import statistics
import subprocess
import sys


# quantile(): depth at the given fraction (0.1, 0.5, 0.9) of a histogram {depth: number of bases}.
def quantile(hist, fraction):
    target = sum(hist.values()) * fraction
    count = 0
    for depth, n in sorted(hist.items()):
        count += n
        if count >= target:
            return depth
    return 0


# esc(): make text safe to put inside the SVG (XML escaping).
def esc(value):
    return str(value).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


# ---- Arguments; which reference sequences to draw ----
bam, out_svg, out_tsv, title = sys.argv[1:5]
skip = sys.argv[5] if len(sys.argv) > 5 else None
# Sequence names and lengths from the BAM index; keep those over 100 kb, except the skipped mito contig.
idx = subprocess.check_output(["samtools", "idxstats", bam], universal_newlines=True)
chroms = []
for line in idx.splitlines():
    fields = line.split("\t")
    if fields[0] != "*" and fields[0] != skip and int(fields[1]) > 100000:
        chroms.append((fields[0], int(fields[1])))

# ---- Read depth per screen column, chromosome by chromosome ----
rows = []
for chrom, length in chroms:
    # Width of one screen column in bp (each chromosome is drawn over ~1,050 columns).
    width_bp = math.ceil(length / 1050)
    # Per-base depth at MAPQ>=20 (-Q) and BQ>=20 (-q), including zero-depth positions (-aa).
    proc = subprocess.Popen(["samtools", "depth", "-aa", "-Q", "20", "-q", "20", "-r", chrom, bam],
                            stdout=subprocess.PIPE, text=True)
    bins, current, hist = [], -1, collections.Counter()
    # Collect each base's depth into its column; when the column changes, store p10, median, p90 and the zero fraction.
    for line in proc.stdout:
        fields = line.split("\t")
        pos, depth = int(fields[1]), int(fields[2])
        b = (pos - 1) // width_bp
        if b != current and current >= 0:
            bins.append((current, quantile(hist, .1), quantile(hist, .5), quantile(hist, .9), hist[0] / sum(hist.values())))
            hist.clear()
        current = b
        hist[depth] += 1
    # Store the last column, and stop if samtools failed.
    if current >= 0:
        bins.append((current, quantile(hist, .1), quantile(hist, .5), quantile(hist, .9), hist[0] / sum(hist.values())))
    if proc.wait() != 0:
        raise RuntimeError("samtools depth failed for " + chrom)
    rows.append((chrom, length, width_bp, bins))

# ---- Per-column table (TSV) ----
with open(out_tsv, "w") as handle:
    handle.write("chrom\tstart\tend\tp10\tmedian\tp90\tzero_fraction\n")
    for chrom, length, bp, bins in rows:
        for i, lo, med, hi, zero in bins:
            handle.write(f"{chrom}\t{i * bp + 1}\t{min((i + 1) * bp, length)}\t{lo}\t{med}\t{hi}\t{zero:.5f}\n")

# ---- SVG: canvas size (one 94-px row per chromosome), left margin, plot width, title and two caption lines ----
W, H = 1390, 90 + 94 * len(rows)
left, plot_w = 145, 1160
svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
       '<rect width="100%" height="100%" fill="white"/>',
       '<style>text{font-family:Arial,sans-serif;fill:#25343b} .small{font-size:12px}.title{font-size:18px;font-weight:bold}</style>',
       f'<text x="145" y="27" class="title">{esc(title)}</text>',
       '<text x="145" y="47" class="small">Full reference length; MAPQ 20, base quality 20. Each screen column shows base-depth median (blue) and 10th–90th percentile (pale).</text>',
       '<text x="145" y="64" class="small">Vertical scales are independent by chromosome; the number under each name is that chromosome\'s typical depth. Red marks: screen columns with &gt;50% uncovered bases.</text>']
summary = []
# One row per chromosome: name, length and typical depth on the left; the y scale = 1.8 x the 98th-percentile median.
for row_index, (chrom, length, bp, bins) in enumerate(rows):
    top = 82 + row_index * 94
    base = top + 66
    typical = statistics.median([x[2] for x in bins])
    scale = max(1, sorted(x[2] for x in bins)[int(.98 * (len(bins) - 1))] * 1.8)
    svg.append(f'<text x="15" y="{top + 25}" class="small">{esc(chrom)}</text>')
    svg.append(f'<text x="15" y="{top + 42}" class="small">{length / 1e6:.2f} Mb | {typical:.0f}x</text>')
    svg.append(f'<line x1="{left}" y1="{base}" x2="{left + plot_w}" y2="{base}" stroke="#c8d2d6"/>')
    # Per column: pale band from p10 to p90, and a red tick under columns that are mostly uncovered.
    for i, lo, med, hi, zero in bins:
        x = left + i * bp / length * plot_w
        ylo = base - min(lo, scale) / scale * 55
        yhi = base - min(hi, scale) / scale * 55
        svg.append(f'<line x1="{x:.2f}" y1="{ylo:.2f}" x2="{x:.2f}" y2="{yhi:.2f}" stroke="#b8dce3" stroke-width="1"/>')
        if zero > .5:
            svg.append(f'<line x1="{x:.2f}" y1="{base + 3}" x2="{x:.2f}" y2="{base + 7}" stroke="#b44141" stroke-width="1"/>')
    # Blue line through the column medians; summary row for standard output.
    points = ' '.join(f'{left + i * bp / length * plot_w:.1f},{base - min(med, scale) / scale * 55:.1f}' for i, _, med, _, _ in bins)
    svg.append(f'<polyline points="{points}" fill="none" stroke="#08778b" stroke-width="1.2"/>')
    summary.append(f"{chrom}\t{length}\t{typical}\t{scale}\t{sum(z > .5 for *_, z in bins)}")
# Close the SVG, write it, and print the per-chromosome summary.
svg.append('</svg>')
with open(out_svg, "w") as handle:
    handle.write("\n".join(svg))
print("chrom\tlength\tmedian_screen_depth\tscale_max\tcolumns_mostly_uncovered")
print("\n".join(summary))
