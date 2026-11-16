"""Synthetic pairtools-stats / preseq / benchmark generators for the QC report tests and demo."""
import os

import numpy as np


def make_stats_text(total=100_000_000, unmapped_frac=0.08, dup_frac=0.15, cis_frac=0.85,
                    short_frac=0.25, shoulder=3.0, orient_bias=0.0, seed=0):
    """Return realistic pairtools-stats text (tab-separated key/value).

    short_frac: approx fraction of cis pairs < 1 kb; shoulder: height of the ~150 bp
    nucleosome shoulder; orient_bias: excess of inward pairs at short range (0 = balanced).
    """
    rng = np.random.default_rng(seed)
    unmapped = int(total * unmapped_frac)
    mapped = total - unmapped
    dups = int(total * dup_frac)
    nodups = mapped - dups
    edges = np.unique(np.round(np.logspace(0, 8.4, 8 * 8 + 1)).astype(np.int64))
    edges = np.r_[0, edges]
    lo, hi = edges[:-1], edges[1:]
    s = np.sqrt(np.maximum(lo, 0.5) * hi)
    w = (hi - lo).astype(float)
    dens = s ** -1.0 * (1 + shoulder * np.exp(-0.5 * (np.log(s / 150.0) / 0.7) ** 2))
    dens *= 1 + 0.6 * (s < 1000) * (short_frac / 0.25)  # tune short-range weight
    wts = dens * w
    wts[hi > 3e8] = 0
    wts /= wts.sum()
    cis = int(nodups * cis_frac)
    trans = nodups - cis
    counts = np.round(wts * cis).astype(np.int64)

    def at_least(bp):
        return int(counts[lo >= bp].sum())

    lines = [("total", total), ("total_unmapped", unmapped), ("total_single_sided_mapped", int(total * 0.02)),
             ("total_mapped", mapped), ("total_dups", dups), ("total_nodups", nodups),
             ("cis", int(counts.sum())), ("trans", trans + (cis - int(counts.sum())))]
    for k in (1, 2, 4, 10, 20, 40):
        lines.append((f"cis_{k}kb+", at_least(k * 1000)))
    uu = int(mapped * 0.9)
    for k, v in (("NU", int(total * 0.05)), ("UN", int(total * 0.05)), ("UU", uu), ("DD", dups),
                 ("UR", int(mapped * 0.04)), ("RU", int(mapped * 0.04)), ("NN", int(total * 0.01)),
                 ("MU", int(total * 0.01))):
        lines.append((f"pair_types/{k}", v))
    lines += [("summary/frac_cis", cis_frac), ("summary/frac_dups", dup_frac)]
    for j in range(len(lo)):
        short = hi[j] <= 1000
        p_in = 0.25 + (orient_bias if short else 0.0)
        probs = np.array([0.25, p_in, 0.25, 0.25])  # ++, +-, -+, --
        probs = probs / probs.sum()
        split = rng.multinomial(int(counts[j]), probs)
        hi_s = "inf" if j == len(lo) - 1 else int(hi[j])
        for o, c in zip(("--", "-+", "+-", "++"), split[::-1]):
            lines.append((f"dist_freq/{int(lo[j])}-{hi_s}/{o}", int(c)))
    lines += [("chrom_freq/chr1/chr1", int(counts.sum() * 0.1)), ("chrom_freq/chr1/chr2", 1000)]
    return "\n".join(f"{k}\t{v}" for k, v in lines) + "\n"


def make_preseq_text(depth=60e6, extrap=2e9, distinct_max=900e6, step=1e8):
    """preseq lc_extrap-style table."""
    xs = np.arange(0, extrap + step, step)
    xs[0] = 0
    ys = distinct_max * xs / (xs + distinct_max / 1.4)
    out = ["TOTAL_READS\tEXPECTED_DISTINCT\tLOWER_0.95CI\tUPPER_0.95CI"]
    for x, y in zip(xs, ys):
        out.append(f"{x:.1f}\t{y:.1f}\t{y * 0.95:.1f}\t{y * 1.05:.1f}")
    return "\n".join(out) + "\n"


def make_benchmark_text(seconds=3600.0, rss=20000.0):
    return ("s\th:m:s\tmax_rss\tmax_vms\tmax_uss\tmax_pss\tio_in\tio_out\tmean_load\tcpu_time\n"
            f"{seconds}\t1:00:00\t{rss}\t{rss * 2}\t{rss}\t{rss}\t1500.0\t900.0\t90.0\t3000.0\n")


def make_group(outdir, name, preseq=True, bench=True, **kw):
    """Create OUTDIR/<name>/ with stats (+ optional preseq and benchmarks)."""
    d = os.path.join(outdir, name)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, f"{name}_dedup_stats.txt"), "w") as fh:
        fh.write(make_stats_text(**kw))
    if preseq:
        with open(os.path.join(d, f"{name}_libraryComplexity.txt"), "w") as fh:
            fh.write(make_preseq_text())
    if bench:
        os.makedirs(os.path.join(d, "benchmarks"), exist_ok=True)
        for rule, s in (("align", 5400.0), ("dedup", 1800.0)):
            with open(os.path.join(d, "benchmarks", f"{name}_{rule}.tsv"), "w") as fh:
                fh.write(make_benchmark_text(s))
    return d


def make_demo(outdir):
    """Three groups with different quality plus decoy directories."""
    make_group(outdir, "Proj_HCT116_B1_T1", total=420_000_000, seed=1)
    make_group(outdir, "Proj_HCT116_B2_T1", total=250_000_000, dup_frac=0.45, cis_frac=0.72, seed=2)
    make_group(outdir, "Proj_K562_B1_T1", total=80_000_000, dup_frac=0.7, cis_frac=0.45, orient_bias=0.25,
               short_frac=0.6, seed=3, preseq=False)
    os.makedirs(os.path.join(outdir, "qc"), exist_ok=True)
    os.makedirs(os.path.join(outdir, "logs_slurm"), exist_ok=True)
    os.makedirs(os.path.join(outdir, "empty_dir"), exist_ok=True)
