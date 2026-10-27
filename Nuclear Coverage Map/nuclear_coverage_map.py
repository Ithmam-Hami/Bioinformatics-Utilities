#!/usr/bin/env python3
# =====================================================================================================================
# nuclear_coverage_map.py  -  whole-genome coverage map (SVG), one panel per chromosome, straight from a BAM
# ---------------------------------------------------------------------------------------------------------------------
# PURPOSE
#   Quick visual check of how evenly reads cover a nuclear genome: uneven libraries, coverage sags, aneuploid or
#   missing chromosomes and mapping problems stand out at a glance. Optionally marks genes of interest (from a BED
#   file), flags the ones with low coverage, and draws a base-resolution zoom panel for each flagged gene, so a
#   low-coverage gene can be told apart from a low region, a mapping-quality problem or ordinary unevenness.
#
# WHAT IT DRAWS
#   Chromosome panels - every reference sequence longer than --min-length, end to end, in ~1,050 screen columns:
#     blue line    median read depth per column          pale band  10th-90th percentile of depth in the column
#     red tick     column where more than half the bases have no reads
#     orange line  share of reads with mapping quality >= --min-mapq (top of its strip = 100%): flat at the top where
#                  reads map confidently, dipping in repeats and other ambiguous sequence
#     red line     a marked gene, labelled "gene depth (ratio to the median gene depth)"
#   Zoom panels (one per marked gene): the gene +/- --zoom-flank bp at base resolution -
#     blue line = depth with mapping quality >= --min-mapq;  grey line = depth with all reads (a gap between the two
#     = reads lost to low mapping quality);  red boxes = the gene's pieces from the BED (e.g. coding exons);
#     dashed line = median gene depth;  orange line = share of reads with mapping quality >= --min-mapq
#
# OUTPUTS
#   <out.svg>   the map
#   <out.tsv>   per-column numbers: chrom, start, end, p10, median, p90, zero_fraction, mapq_share
#   stdout      one summary row per chromosome (length, median depth, y-scale maximum, mostly-uncovered columns), and
#               with --genes one row per gene (depth, ratio to the median gene depth, marked yes / no)
#
# REQUIREMENTS  Python 3 (standard library only) and samtools on PATH. The BAM must be coordinate-sorted and indexed.
#
# USAGE
#   python3 nuclear_coverage_map.py <bam> <out.svg> <out.tsv> "<title>" [options]
#   --skip CONTIG          leave out one sequence, e.g. the mitochondrial contig of a whole-genome BAM
#   --genes BED            genes to measure (BED: contig, start, end, name; several lines with the same name = one gene
#                          made of pieces, e.g. exons). Gene depth = median depth over its bases.
#   --flag-depth 20        mark a gene when its depth is below this ...
#   --flag-ratio 0.6       ... and below this fraction of the median gene depth
#   --mark-all             mark every gene in the BED instead of only the flagged ones
#   --zoom-flank 5000      bp shown on each side of a marked gene in its zoom panel (0 = no zoom panels)
#   --min-mapq 40          minimum mapping quality counted (also the threshold of the orange track)
#   --min-baseq 0          minimum base quality counted (0 = every base a read covers)
#   --no-duplicates        leave out duplicate-marked reads (default: counted; at high depth most "duplicates" are
#                          independent molecules that start at the same position by chance)
#   --min-length 100000    shortest sequence drawn (bp)
# =====================================================================================================================
import argparse
import collections
import math
import statistics
import subprocess

p = argparse.ArgumentParser(description="Whole-genome coverage map (SVG) from a BAM, one panel per chromosome.")
p.add_argument("bam"); p.add_argument("out_svg"); p.add_argument("out_tsv"); p.add_argument("title")
p.add_argument("--skip", default=None, help="sequence to leave out (e.g. the mito contig)")
p.add_argument("--genes", default=None, help="BED of genes to measure and mark")
p.add_argument("--flag-depth", type=float, default=20)
p.add_argument("--flag-ratio", type=float, default=0.6)
p.add_argument("--mark-all", action="store_true")
p.add_argument("--zoom-flank", type=int, default=5000)
p.add_argument("--min-mapq", type=int, default=40)
p.add_argument("--min-baseq", type=int, default=0)
p.add_argument("--no-duplicates", action="store_true")
p.add_argument("--min-length", type=int, default=100000)
a = p.parse_args()
DUP = [] if a.no_duplicates else ["-g", "DUP"]                   # samtools depth skips duplicates unless told not to


def quantile(hist, fraction):
    """Depth at the given fraction (0.1, 0.5, 0.9) of a histogram {depth: number of bases}."""
    target, count = sum(hist.values()) * fraction, 0
    for depth, n in sorted(hist.items()):
        count += n
        if count >= target:
            return depth
    return 0


def esc(value):
    """Make text safe inside the SVG."""
    return str(value).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def base_depth(chrom, start, end, mapq):
    """Depth at every base start..end (1-based, inclusive) counting reads with mapping quality >= mapq."""
    out = subprocess.check_output(["samtools", "depth", "-aa", "-Q", str(mapq), "-q", str(a.min_baseq)] + DUP +
                                  ["-r", "%s:%d-%d" % (chrom, start, end), a.bam], universal_newlines=True)
    d = {int(l.split("\t")[1]): int(l.split("\t")[2]) for l in out.splitlines()}
    return [d.get(q, 0) for q in range(start, end + 1)]


# ---- 1. sequences to draw: longer than --min-length, except the skipped one ----
idx = subprocess.check_output(["samtools", "idxstats", a.bam], universal_newlines=True)
chroms = [(f[0], int(f[1])) for f in (l.split("\t") for l in idx.splitlines())
          if f[0] not in ("*", a.skip) and int(f[1]) > a.min_length]
lengths = dict(chroms)

# ---- 2. genes (optional): depth of each, the median gene depth, and which ones to mark ----
pieces = collections.OrderedDict()                                  # name -> (contig, [(start, end)] 0-based BED)
if a.genes:
    for line in open(a.genes):
        f = line.rstrip("\n").split("\t")
        if len(f) < 4 or f[0].startswith(("#", "track", "browser")) or f[0] not in lengths:
            continue
        pieces.setdefault(f[3], (f[0], []))[1].append((int(f[1]), int(f[2])))
gene_depth = {}
for name, (c, ps) in pieces.items():
    bases = []
    for s, e in ps:
        bases += base_depth(c, s + 1, e, a.min_mapq)
    gene_depth[name] = statistics.median(bases) if bases else 0
med_gene = statistics.median(gene_depth.values()) if gene_depth else 0
marked = [n for n, d in gene_depth.items()
          if a.mark_all or (d < a.flag_depth and d < a.flag_ratio * med_gene)]
marked.sort(key=lambda n: (pieces[n][0], min(s for s, _ in pieces[n][1])))
marks = collections.defaultdict(list)                               # contig -> [(start, end, label)]
for n in marked:
    c, ps = pieces[n]
    marks[c].append((min(s for s, _ in ps), max(e for _, e in ps),
                     "%s %.0fx (%.2f)" % (n[:22], gene_depth[n], gene_depth[n] / med_gene if med_gene else 0)))

# ---- 3. per chromosome: depth per screen column, and the share of reads with good mapping quality ----
rows = []
for chrom, length in chroms:
    bp = math.ceil(length / 1050)                                   # bases per screen column
    proc = subprocess.Popen(["samtools", "depth", "-aa", "-Q", str(a.min_mapq), "-q", str(a.min_baseq)] + DUP +
                            ["-r", chrom, a.bam], stdout=subprocess.PIPE, universal_newlines=True)
    bins, current, hist = {}, -1, collections.Counter()
    for line in proc.stdout:
        f = line.split("\t")
        b = (int(f[1]) - 1) // bp
        if b != current and current >= 0:
            bins[current] = (quantile(hist, .1), quantile(hist, .5), quantile(hist, .9), hist[0] / sum(hist.values()))
            hist.clear()
        current = b
        hist[int(f[2])] += 1
    if current >= 0:
        bins[current] = (quantile(hist, .1), quantile(hist, .5), quantile(hist, .9), hist[0] / sum(hist.values()))
    if proc.wait() != 0:
        raise RuntimeError("samtools depth failed for " + chrom)
    # reads by start column (primary alignments): share with mapping quality >= --min-mapq
    proc = subprocess.Popen(["samtools", "view", "-F", "0x904", a.bam, chrom], stdout=subprocess.PIPE,
                            universal_newlines=True)
    tot, good = collections.Counter(), collections.Counter()
    for line in proc.stdout:
        f = line.split("\t", 6)
        b = (int(f[3]) - 1) // bp
        tot[b] += 1
        if int(f[4]) >= a.min_mapq:
            good[b] += 1
    if proc.wait() != 0:
        raise RuntimeError("samtools view failed for " + chrom)
    rows.append((chrom, length, bp, bins, {b: good[b] / tot[b] for b in tot}))

# ---- 4. zoom data: each marked gene +/- flank at base resolution, with good-MAPQ reads and with all reads ----
zdata = []
if a.zoom_flank > 0:
    for n in marked:
        c, ps = pieces[n]
        s, e = min(x for x, _ in ps), max(y for _, y in ps)
        z0, z1 = max(1, s + 1 - a.zoom_flank), min(lengths[c], e + a.zoom_flank)
        zdata.append((n, c, s, e, z0, z1, base_depth(c, z0, z1, a.min_mapq), base_depth(c, z0, z1, 0)))

# ---- 5. draw ----
BLUE, BAND, RED, ORANGE, GREY, INK = "#08778b", "#b8dce3", "#b44141", "#e08a2c", "#9aa5ab", "#25343b"
GENE = "#c0392b"
PANEL, ZCOLS, ZW, ZH = 132, 3, 430, 190                             # chromosome panel height; zoom grid
ztop = 104 + PANEL * len(rows) + 10
W = 1390
H = ztop + ((36 + math.ceil(len(zdata) / ZCOLS) * ZH) if zdata else 0)
left, plot_w = 145, 1160
GE = "&#8805;"
svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">' % (W, H, W, H),
       '<rect width="100%" height="100%" fill="white"/>',
       '<style>text{font-family:Arial,sans-serif;fill:%s} .small{font-size:12px} .tiny{font-size:10px} '
       '.title{font-size:18px;font-weight:bold}</style>' % INK,
       '<text x="145" y="27" class="title">%s</text>' % esc(a.title)]


def legend(items, x, y):
    """One row of colour keys: (kind, colour, label)."""
    for kind, colour, label in items:
        if kind == "line":
            svg.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="%s" stroke-width="2"/>' % (x, y, x + 18, y, colour))
        elif kind == "band":
            svg.append('<rect x="%d" y="%d" width="18" height="12" fill="%s"/>' % (x, y - 6, colour))
        elif kind == "tick":
            svg.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="%s" stroke-width="2"/>' % (x + 9, y - 5, x + 9, y + 5, colour))
        elif kind == "gene":
            svg.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="%s" stroke-width="2"/>' % (x + 9, y - 8, x + 9, y + 8, colour))
        elif kind == "box":
            svg.append('<rect x="%d" y="%d" width="18" height="8" fill="%s"/>' % (x, y - 4, colour))
        elif kind == "dash":
            svg.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="%s" stroke-dasharray="4,3"/>' % (x, y, x + 18, y, colour))
        svg.append('<text x="%d" y="%d" class="small">%s</text>' % (x + 24, y + 4, label))
        x += 24 + int(len(label.replace("&gt;", ">").replace(GE, ">")) * 6.6) + 22


keys = [("line", BLUE, "median depth"), ("band", BAND, "10th-90th percentile depth"),
        ("tick", RED, "&gt;50% of bases at 0 reads"), ("line", ORANGE, "share of reads with MAPQ %s %d" % (GE, a.min_mapq))]
if marked:
    keys.append(("gene", GENE, ("gene" if a.mark_all else "low-coverage gene") + " (depth, ratio to %.0fx)" % med_gene))
legend(keys, 145, 52)

tsv = open(a.out_tsv, "w")
tsv.write("chrom\tstart\tend\tp10\tmedian\tp90\tzero_fraction\tmapq_share\n")
print("chrom\tlength\tmedian_depth\ty_max\tmostly_uncovered_columns")
for ri, (chrom, length, bp, bins, mq) in enumerate(rows):
    top = 82 + ri * PANEL
    base = top + 72                                                 # depth baseline
    mqtop, mqh = base + 12, 16                                      # mapping-quality strip
    meds = [v[1] for v in bins.values()] or [0]
    typical = statistics.median(meds)
    scale = max(1, sorted(meds)[int(.98 * (len(meds) - 1))] * 1.8)
    for i in sorted(bins):
        lo, med, hi, zero = bins[i]
        tsv.write("%s\t%d\t%d\t%d\t%d\t%d\t%.4f\t%s\n" % (chrom, i * bp, min(length, (i + 1) * bp), lo, med, hi, zero,
                                                          "%.4f" % mq[i] if i in mq else "NA"))
    print("%s\t%d\t%.0f\t%.0f\t%d" % (chrom, length, typical, scale, sum(1 for v in bins.values() if v[3] > .5)))
    svg.append('<text x="15" y="%d" class="small">%s</text>' % (top + 30, esc(chrom)))
    svg.append('<text x="15" y="%d" class="small">%.2f Mb | %.0fx</text>' % (top + 47, length / 1e6, typical))
    svg.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#c8d2d6"/>' % (left, base, left + plot_w, base))
    svg.append('<rect x="%d" y="%d" width="%d" height="%d" fill="none" stroke="#e1e6e8"/>' % (left, mqtop, plot_w, mqh))
    for i in sorted(bins):
        lo, med, hi, zero = bins[i]
        x = left + i * bp / length * plot_w
        svg.append('<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="%s" stroke-width="1"/>'
                   % (x, base - min(lo, scale) / scale * 58, x, base - min(hi, scale) / scale * 58, BAND))
        if zero > .5:
            svg.append('<line x1="%.2f" y1="%d" x2="%.2f" y2="%d" stroke="%s" stroke-width="1"/>' % (x, base + 3, x, base + 8, RED))
    seg = []                                                        # MAPQ line, broken where a column has no reads
    for i in range(int(length / bp) + 1):
        if i in mq:
            seg.append("%.1f,%.1f" % (left + i * bp / length * plot_w, mqtop + (1 - mq[i]) * mqh))
        elif seg:
            svg.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="1"/>' % (" ".join(seg), ORANGE))
            seg = []
    if seg:
        svg.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="1"/>' % (" ".join(seg), ORANGE))
    pts = " ".join("%.1f,%.1f" % (left + i * bp / length * plot_w, base - min(bins[i][1], scale) / scale * 58) for i in sorted(bins))
    svg.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="1.2"/>' % (pts, BLUE))
    for k, (s, e, label) in enumerate(sorted(marks.get(chrom, []))):
        x = left + (s + e) / 2 / length * plot_w
        svg.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="%s" stroke-width="1.4"/>' % (x, top + 8, x, mqtop + mqh, GENE))
        svg.append('<text x="%.1f" y="%d" class="tiny" fill="%s" text-anchor="middle">%s</text>'
                   % (x, top + 6 + (k % 2) * 11, GENE, esc(label)))
tsv.close()

if zdata:
    svg.append('<text x="15" y="%d" class="title">%s, +/- %s kb</text>'
               % (ztop + 18, "Marked genes" if a.mark_all else "Low-coverage genes", format(a.zoom_flank / 1000, "g")))
    legend([("line", BLUE, "depth, MAPQ %s %d" % (GE, a.min_mapq)), ("line", GREY, "depth, all reads"),
            ("box", GENE, "gene pieces (BED)"), ("dash", "#5f5e5a", "median gene depth (%.0fx)" % med_gene),
            ("line", ORANGE, "share of reads with MAPQ %s %d" % (GE, a.min_mapq))], 330, ztop + 13)
for n_i, (name, c, s, e, z0, z1, dq, da) in enumerate(zdata):
    x0 = 40 + (n_i % ZCOLS) * (ZW + 15)
    y0 = ztop + 36 + (n_i // ZCOLS) * ZH
    dtop, dh = y0 + 26, 100                                         # depth area
    exy, mqy, mqh = dtop + dh + 8, dtop + dh + 24, 16               # gene-piece track, MAPQ strip
    size = z1 - z0 + 1
    ymax = max(med_gene * 1.6, sorted(da)[int(.98 * (len(da) - 1))] * 1.1, 1)
    X = lambda q: x0 + (q - z0) / size * ZW
    Y = lambda v: dtop + dh - min(v, ymax) / ymax * dh
    d = gene_depth[name]
    svg.append('<text x="%d" y="%d" class="small"><tspan font-weight="bold">%s</tspan> | %s:%d-%d | %.0fx (%.2f)</text>'
               % (x0, y0 + 14, esc(name[:28]), esc(c), s + 1, e, d, d / med_gene if med_gene else 0))
    svg.append('<rect x="%d" y="%d" width="%d" height="%d" fill="none" stroke="#e1e6e8"/>' % (x0, dtop, ZW, dh))
    svg.append('<text x="%d" y="%d" class="tiny" text-anchor="end">%.0f</text>' % (x0 - 2, dtop + 8, ymax))
    svg.append('<text x="%d" y="%d" class="tiny" text-anchor="end">0</text>' % (x0 - 2, dtop + dh))
    step = max(1, size // ZW)                                       # one point per pixel column
    pq, pa, pm = [], [], []
    for i in range(0, size, step):
        al, ql = da[i:i + step], dq[i:i + step]
        x = X(z0 + i + len(al) / 2)
        pa.append("%.1f,%.1f" % (x, Y(sum(al) / len(al))))
        pq.append("%.1f,%.1f" % (x, Y(sum(ql) / len(ql))))
        if sum(al) > 0:
            pm.append("%.1f,%.1f" % (x, mqy + (1 - sum(ql) / sum(al)) * mqh))
    if med_gene:
        svg.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#5f5e5a" stroke-dasharray="4,3"/>' % (x0, Y(med_gene), x0 + ZW, Y(med_gene)))
    svg.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="1"/>' % (" ".join(pa), GREY))
    svg.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="1.3"/>' % (" ".join(pq), BLUE))
    for ps, pe in pieces[name][1]:
        svg.append('<rect x="%.1f" y="%d" width="%.1f" height="8" fill="%s"/>' % (X(ps + 1), exy, max(1, X(pe) - X(ps + 1)), GENE))
    svg.append('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="#e1e6e8"/>' % (x0, exy + 4, x0 + ZW, exy + 4))
    svg.append('<rect x="%d" y="%d" width="%d" height="%d" fill="none" stroke="#e1e6e8"/>' % (x0, mqy, ZW, mqh))
    if pm:
        svg.append('<polyline points="%s" fill="none" stroke="%s" stroke-width="1"/>' % (" ".join(pm), ORANGE))
svg.append('</svg>')
open(a.out_svg, "w").write("\n".join(svg))

if gene_depth:
    print("\ngene\tdepth\tratio_to_median_gene\tmarked")
    for n, d in gene_depth.items():
        print("%s\t%.1f\t%.2f\t%s" % (n, d, d / med_gene if med_gene else 0, "yes" if n in marked else "no"))
