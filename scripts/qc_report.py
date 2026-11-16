#!/usr/bin/env python3
"""Static, self-contained Micro-C QC dashboard (a small, Micro-C-specific "MultiQC").

Usage:
    python scripts/qc_report.py --outdir OUTDIR --config config/MicroCConfig.yaml \
        --out OUTDIR/qc/MicroC_qc_report.html

Discovery (one "sample" per group directory OUTDIR/<group>/):
    <group>_dedup_stats.txt       pairtools stats (required; dirs without it are skipped)
    <group>_libraryComplexity.txt preseq lc_extrap table (optional)
    benchmarks/*.tsv              Snakemake benchmark files (optional)
    <group>.mcool                 summarised only if the `cooler` module imports (optional)

The output is a single HTML file with inline CSS/JS and one embedded copy of
plotly.js, so it works offline and can be e-mailed. Dependencies: pandas,
numpy, pyyaml, plotly (cooler optional).
"""
from __future__ import annotations

import argparse
import datetime as _dt
import glob
import html
import math
import os
import re
import sys
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import yaml
from plotly.offline import get_plotlyjs

SKIP_DIRS = {"qc", "logs_slurm", "logs", ".snakemake", "benchmarks", "aligned"}

# --------------------------------------------------------------------------
# Defaults, metric definitions, palettes
# --------------------------------------------------------------------------

DEFAULT_THRESHOLDS: Dict[str, Dict[str, Any]] = {
    "mapped_pct": {"direction": "higher_is_better", "warn": 80, "fail": 60},
    "duplicate_pct": {"direction": "lower_is_better", "warn": 40, "fail": 60},
    "cis_pct": {"direction": "higher_is_better", "warn": 70, "fail": 50},
    "cis_ge1kb_pct": {"direction": "higher_is_better", "warn": 40, "fail": 25},
    "cis_lt1kb_pct": {"direction": "lower_is_better", "warn": 40, "fail": 55},
    "cis_ge10kb_pct": {"direction": "higher_is_better", "warn": 20, "fail": 10},
    "unique_pairs_million": {"direction": "higher_is_better", "warn": 100, "fail": 20},
}

#: label, unit, definition (rendered in the Methods section and table headers)
METRIC_INFO: Dict[str, Dict[str, str]] = {
    "mapped_pct": {"label": "Mapped", "unit": "%",
                   "desc": "total_mapped / total read pairs."},
    "duplicate_pct": {"label": "Duplicates", "unit": "%",
                      "desc": "total_dups / total read pairs (all pairs incl. unmapped, as in get_qc.py)."},
    "cis_pct": {"label": "Cis", "unit": "%",
                "desc": "cis / non-duplicate pairs."},
    "cis_ge1kb_pct": {"label": "Cis ≥ 1 kb", "unit": "%",
                      "desc": "cis_1kb+ / non-duplicate pairs (informative long-range cis)."},
    "cis_lt1kb_pct": {"label": "Cis < 1 kb", "unit": "%",
                      "desc": "(cis - cis_1kb+) / non-duplicate pairs (dangling-end / self-ligation-like short pairs)."},
    "cis_ge10kb_pct": {"label": "Cis ≥ 10 kb", "unit": "%",
                       "desc": "cis_10kb+ / non-duplicate pairs."},
    "unique_pairs_million": {"label": "Unique pairs", "unit": "M",
                             "desc": "total_nodups in millions (absolute depth of the final .pairs)."},
}
EXTRA_INFO = {
    "total_million": ("Total pairs", "M", "total read pairs in millions."),
    "dup_of_mapped_pct": ("Dup / mapped", "%", "total_dups / total_mapped (informational, not thresholded)."),
    "valid_pairs_pct": ("Valid pairs", "%", "(trans + cis_1kb+) / non-duplicate pairs, as in get_qc.py (informational)."),
}

STATUS_RANK = {"pass": 0, "warn": 1, "fail": 2, "na": -1}
STATUS_SYMBOL = {"pass": "✓", "warn": "!", "fail": "✕", "na": "–"}

# Okabe-Ito colour-blind-safe palette (black omitted: invisible on dark backgrounds)
PALETTE = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442", "#999999"]
DASHES = ["solid", "dash", "dot", "dashdot", "longdash", "longdashdot"]
BIN_COLOURS = ["#D55E00", "#E69F00", "#F0E442", "#009E73", "#56B4E9", "#0072B2", "#332288"]
PAIR_TYPE_COLOURS = {
    "UU": "#0072B2", "UR": "#56B4E9", "RU": "#56B4E9", "NU": "#009E73", "UN": "#009E73",
    "MU": "#CC79A7", "UM": "#CC79A7", "NM": "#999999", "MN": "#999999", "NN": "#BBBBBB",
    "MM": "#999999", "NR": "#BBBBBB", "RN": "#BBBBBB", "MR": "#BBBBBB", "RM": "#BBBBBB",
    "RR": "#BBBBBB", "WW": "#E69F00", "XX": "#E69F00", "DD": "#D55E00",
}
PAIR_TYPE_NAMES = {
    "UU": "unique-unique", "UR": "unique-rescued", "RU": "rescued-unique", "NU": "null-unique",
    "UN": "unique-null", "MU": "multi-unique", "UM": "unique-multi", "DD": "duplicates",
    "NN": "null-null", "MM": "multi-multi", "NM": "null-multi", "MN": "multi-null",
}
GREY = "#8c8c8c"

# --------------------------------------------------------------------------
# Small helpers
# --------------------------------------------------------------------------


def safe_div(num: Optional[float], den: Optional[float]) -> Optional[float]:
    """Return num/den, or None if either is missing/NaN or den is zero."""
    if num is None or den is None:
        return None
    try:
        if math.isnan(num) or math.isnan(den) or den == 0:
            return None
    except TypeError:
        return None
    return float(num) / float(den)


def pct(num: Optional[float], den: Optional[float]) -> Optional[float]:
    """Percentage num/den*100 or None (zero-safe)."""
    r = safe_div(num, den)
    return None if r is None else 100.0 * r


def fmt(value: Optional[float], digits: int = 1, unit: str = "") -> str:
    """Format a number for display; None/NaN -> 'n/a'."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "n/a"
    return f"{value:,.{digits}f}{unit}"


def esc(x: Any) -> str:
    """HTML-escape any value."""
    return html.escape(str(x), quote=True)


def _sortval(v: Any) -> str:
    """Value for a data-sort attribute ('' when missing)."""
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return ""
    return esc(v)


# --------------------------------------------------------------------------
# Parsers
# --------------------------------------------------------------------------


def parse_pairtools_stats(path: str) -> Dict[str, float]:
    """Parse a pairtools stats file (``key<TAB>value`` per line) into a dict.

    Non-numeric values and malformed lines are ignored. Keys keep their full
    path, e.g. ``pair_types/UU`` or ``dist_freq/10-12/+-``.
    """
    stats: Dict[str, float] = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t") if "\t" in line else line.split()
            if len(parts) < 2:
                continue
            try:
                stats[parts[0].strip()] = float(parts[1])
            except ValueError:
                continue
    return stats


_DIST_RE = re.compile(r"^dist_freq/(\d+)-(\d+|inf)/([+-]{2})$")


def parse_dist_freq(stats: Dict[str, float]) -> pd.DataFrame:
    """Extract ``dist_freq`` rows into a DataFrame.

    Returns columns ``lo, hi`` (bp; ``hi`` = inf for an open last bin) and the
    four strand orientations ``++, +-, -+, --`` (zero-filled), sorted by ``lo``.
    Empty DataFrame (with those columns) if the stats contain no dist_freq.
    """
    orients = ["++", "+-", "-+", "--"]
    rows: Dict[tuple, Dict[str, float]] = {}
    for key, val in stats.items():
        m = _DIST_RE.match(key)
        if not m:
            continue
        lo, hi = int(m.group(1)), (math.inf if m.group(2) == "inf" else int(m.group(2)))
        rows.setdefault((lo, hi), {})[m.group(3)] = val
    if not rows:
        return pd.DataFrame(columns=["lo", "hi"] + orients)
    recs = [{"lo": lo, "hi": hi, **{o: d.get(o, 0.0) for o in orients}} for (lo, hi), d in rows.items()]
    return pd.DataFrame(recs).sort_values("lo").reset_index(drop=True)


def parse_preseq(path: str) -> pd.DataFrame:
    """Parse preseq ``lc_extrap`` output.

    Columns: reads, expected, lower, upper (all counts). Empty DataFrame when
    the file is empty or unparsable.
    """
    cols = ["reads", "expected", "lower", "upper"]
    rows = []
    try:
        with open(path) as fh:
            for line in fh:
                p = line.split()
                if len(p) < 4:
                    continue
                try:
                    rows.append([float(x) for x in p[:4]])
                except ValueError:
                    continue  # header
    except OSError:
        pass
    return pd.DataFrame(rows, columns=cols)


BENCH_COLS = ["sample", "rule", "seconds", "max_rss_mb", "io_in_mb", "io_out_mb"]


def parse_benchmarks(group_dir: str, group: str) -> pd.DataFrame:
    """Collect ``benchmarks/*.tsv`` of one group into one row per rule.

    Columns: sample, rule, seconds, max_rss_mb, io_in_mb, io_out_mb (NaN when
    a column is absent). Empty DataFrame if there are no benchmark files.
    """
    recs = []
    for f in sorted(glob.glob(os.path.join(group_dir, "benchmarks", "*.tsv"))):
        try:
            t = pd.read_csv(f, sep="\t")
        except Exception:
            continue
        if t.empty:
            continue
        rule = os.path.splitext(os.path.basename(f))[0]
        if rule.startswith(group + "_"):
            rule = rule[len(group) + 1:]

        def col(name: str) -> float:
            if name not in t.columns:
                return float("nan")
            return float(pd.to_numeric(t[name], errors="coerce").max())

        recs.append({"sample": group, "rule": rule, "seconds": col("s"), "max_rss_mb": col("max_rss"),
                     "io_in_mb": col("io_in"), "io_out_mb": col("io_out")})
    return pd.DataFrame(recs, columns=BENCH_COLS)


def parse_mcool(path: str) -> Optional[Dict[str, Any]]:
    """Summarise an .mcool with ``cooler`` (resolutions, contacts, balancing).

    Returns None silently if cooler is not importable or the file cannot be read.
    """
    try:
        import cooler  # optional dependency
    except Exception:
        return None
    try:
        def res(u: str) -> int:
            tail = u.rsplit("/", 1)[-1]
            return int(tail) if tail.isdigit() else 0

        uris = sorted(cooler.fileops.list_coolers(path), key=res)
        if not uris:
            return None
        finest = cooler.Cooler(f"{path}::{uris[0]}")
        info = finest.info
        return {
            "n_resolutions": len(uris),
            "finest_bp": int(finest.binsize) if finest.binsize else None,
            "contacts": float(info.get("sum", float("nan"))),
            "nnz": float(info.get("nnz", float("nan"))),
            "balanced": "weight" in finest.bins().columns,
        }
    except Exception:
        return None


def load_thresholds(config_path: Optional[str]) -> Dict[str, Dict[str, Any]]:
    """Read ``qc_thresholds`` from the YAML config; built-in defaults if absent/unreadable."""
    thr = {k: dict(v) for k, v in DEFAULT_THRESHOLDS.items()}
    if config_path and os.path.exists(config_path):
        try:
            with open(config_path) as fh:
                cfg = yaml.safe_load(fh) or {}
            user = cfg.get("qc_thresholds") or {}
            thr = {k: dict(v) for k, v in user.items() if isinstance(v, dict)} or thr
        except Exception as exc:  # malformed YAML shouldn't kill the report
            print(f"qc_report: could not read thresholds from {config_path}: {exc}", file=sys.stderr)
    return thr


# --------------------------------------------------------------------------
# Metrics and classification
# --------------------------------------------------------------------------

DIST_BINS = [("<1 kb", None, "cis_1kb+"), ("1–2 kb", "cis_1kb+", "cis_2kb+"),
             ("2–4 kb", "cis_2kb+", "cis_4kb+"), ("4–10 kb", "cis_4kb+", "cis_10kb+"),
             ("10–20 kb", "cis_10kb+", "cis_20kb+"), ("20–40 kb", "cis_20kb+", "cis_40kb+"),
             ("≥40 kb", "cis_40kb+", None)]


def nodup_count(stats: Dict[str, float]) -> Optional[float]:
    """Non-duplicate pair count: total_nodups, else total_mapped - total_dups, else None."""
    if stats.get("total_nodups") is not None:
        return stats["total_nodups"]
    if stats.get("total_mapped") is not None and stats.get("total_dups") is not None:
        return stats["total_mapped"] - stats["total_dups"]
    return None


def compute_metrics(stats: Dict[str, float]) -> Dict[str, Optional[float]]:
    """Compute the dashboard metrics from a parsed pairtools stats dict.

    Percentages of non-duplicate pairs unless noted (see METRIC_INFO). Every
    value is None when it cannot be computed (missing key or zero denominator).
    """
    g = stats.get
    total, mapped, dups = g("total"), g("total_mapped"), g("total_dups")
    nodups = nodup_count(stats)
    cis, trans, cis1 = g("cis"), g("trans"), g("cis_1kb+")
    cis_lt1 = None if cis is None or cis1 is None else max(cis - cis1, 0.0)
    valid = None if trans is None or cis1 is None else trans + cis1
    return {
        "total_million": None if total is None else total / 1e6,
        "mapped_pct": pct(mapped, total),
        "duplicate_pct": pct(dups, total),
        "dup_of_mapped_pct": pct(dups, mapped),
        "cis_pct": pct(cis, nodups),
        "cis_ge1kb_pct": pct(cis1, nodups),
        "cis_lt1kb_pct": pct(cis_lt1, nodups),
        "cis_ge10kb_pct": pct(g("cis_10kb+"), nodups),
        "valid_pairs_pct": pct(valid, nodups),
        "unique_pairs_million": None if nodups is None else nodups / 1e6,
    }


def classify(value: Optional[float], rule: Dict[str, Any]) -> str:
    """Return 'pass' / 'warn' / 'fail' (or 'na' if value is missing).

    higher_is_better: value < fail -> fail; value < warn -> warn.
    lower_is_better:  value > fail -> fail; value > warn -> warn.
    A value exactly on an edge is therefore not penalised.
    """
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "na"
    warn, fail = rule.get("warn"), rule.get("fail")
    if rule.get("direction", "higher_is_better") == "higher_is_better":
        if fail is not None and value < fail:
            return "fail"
        if warn is not None and value < warn:
            return "warn"
    else:
        if fail is not None and value > fail:
            return "fail"
        if warn is not None and value > warn:
            return "warn"
    return "pass"


def worst_status(statuses) -> str:
    """Worst of several statuses (fail > warn > pass); 'na' if none are usable."""
    best = "na"
    for s in statuses:
        if STATUS_RANK[s] > STATUS_RANK[best]:
            best = s
    return best


def classify_sample(metrics: Dict[str, Optional[float]], thresholds: Dict[str, Dict[str, Any]]) -> Dict[str, str]:
    """Status per thresholded metric plus an 'overall' entry for one sample."""
    out = {m: classify(metrics.get(m), r) for m, r in thresholds.items()}
    out["overall"] = worst_status(out.values())
    return out


def cis_distance_breakdown(stats: Dict[str, float]) -> Dict[str, Optional[float]]:
    """Pair counts per cis distance bin plus trans (counts, not percentages); None if unavailable."""
    out: Dict[str, Optional[float]] = {}
    cis = stats.get("cis")
    for label, lo_key, hi_key in DIST_BINS:
        lo = cis if lo_key is None else stats.get(lo_key)
        hi = 0.0 if hi_key is None else stats.get(hi_key)
        out[label] = None if lo is None or hi is None else max(lo - hi, 0.0)
    out["trans"] = stats.get("trans")
    return out


def contact_decay(dist: pd.DataFrame) -> pd.DataFrame:
    """P(s) per distance bin from dist_freq counts.

    P(s) = (pairs in bin / bin width in bp) / (all pairs in dist_freq), at the
    geometric bin centre. The first bin (lo = 0), the open last bin and empty
    bins are dropped (cannot be placed on a log axis). Also returns the
    in/out/same-strand fractions per bin.
    """
    cols = ["s", "p", "in_frac", "out_frac", "same_frac", "count"]
    if dist.empty:
        return pd.DataFrame(columns=cols)
    d = dist.copy()
    d["count"] = d[["++", "+-", "-+", "--"]].sum(axis=1)
    grand = d["count"].sum()
    d = d[(d["lo"] > 0) & np.isfinite(d["hi"].astype(float)) & (d["count"] > 0)].copy()
    if d.empty or grand <= 0:
        return pd.DataFrame(columns=cols)
    d["s"] = np.sqrt(d["lo"].astype(float) * d["hi"].astype(float))
    d["p"] = d["count"] / (d["hi"].astype(float) - d["lo"].astype(float)) / grand
    # pairtools: "+-" = read1 forward upstream, read2 reverse downstream = inward-facing
    d["in_frac"] = d["+-"] / d["count"]
    d["out_frac"] = d["-+"] / d["count"]
    d["same_frac"] = (d["++"] + d["--"]) / d["count"]
    return d[cols].reset_index(drop=True)


# --------------------------------------------------------------------------
# Data model
# --------------------------------------------------------------------------


def discover_samples(outdir: str) -> List[Dict[str, Any]]:
    """Find group directories under ``outdir`` and load every available input.

    Each sample dict has: name, stats, metrics, dist, decay, preseq, bench, mcool.
    Directories named qc / logs_slurm / etc. or lacking ``<group>_dedup_stats.txt`` are skipped.
    """
    samples: List[Dict[str, Any]] = []
    if not os.path.isdir(outdir):
        return samples
    for name in sorted(os.listdir(outdir)):
        d = os.path.join(outdir, name)
        if not os.path.isdir(d) or name in SKIP_DIRS or name.startswith("."):
            continue
        sp = os.path.join(d, f"{name}_dedup_stats.txt")
        if not os.path.isfile(sp):
            continue
        stats = parse_pairtools_stats(sp)
        if not stats:
            continue
        pp = os.path.join(d, f"{name}_libraryComplexity.txt")
        preseq = parse_preseq(pp) if os.path.isfile(pp) else parse_preseq(os.devnull)
        mp = os.path.join(d, f"{name}.mcool")
        dist = parse_dist_freq(stats)
        samples.append({
            "name": name, "stats": stats, "metrics": compute_metrics(stats),
            "dist": dist, "decay": contact_decay(dist), "preseq": preseq,
            "bench": parse_benchmarks(d, name),
            "mcool": parse_mcool(mp) if os.path.isfile(mp) else None,
        })
    return samples


def sample_style(samples: List[Dict[str, Any]]) -> Dict[str, Dict[str, str]]:
    """Stable colour + dash per sample (dash cycles after the palette is exhausted)."""
    return {s["name"]: {"colour": PALETTE[i % len(PALETTE)], "dash": DASHES[(i // len(PALETTE)) % len(DASHES)]}
            for i, s in enumerate(samples)}


# --------------------------------------------------------------------------
# Plots (each returns a plotly Figure, or None when there is nothing to show)
# --------------------------------------------------------------------------

_AXIS = dict(gridcolor="rgba(128,128,128,0.25)", linecolor="rgba(128,128,128,0.6)",
             zerolinecolor="rgba(128,128,128,0.4)", automargin=True)


def _layout(fig: go.Figure, height: int = 420, **kw: Any) -> go.Figure:
    """Apply the shared transparent, theme-neutral layout (text colour comes from CSS)."""
    fig.update_layout(
        height=height, margin=dict(l=60, r=20, t=30, b=50), paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)", font=dict(family="system-ui, sans-serif", size=13),
        legend=dict(orientation="h", yanchor="top", y=-0.2, x=0), hoverlabel=dict(namelength=-1), **kw)
    fig.update_xaxes(**_AXIS)
    fig.update_yaxes(**_AXIS)
    return fig


def _bar_height(n: int, base: int = 150, per: int = 38) -> int:
    return max(260, base + per * n)


def plot_pair_types(samples: List[Dict[str, Any]]) -> Optional[go.Figure]:
    """Horizontal stacked bars of pairtools pair types as % of total read pairs.

    Falls back to unmapped / duplicate / non-duplicate if ``pair_types/*`` is absent.
    """
    names = [s["name"] for s in samples]
    series: Dict[str, List[Optional[float]]] = {}
    for i, s in enumerate(samples):
        st = s["stats"]
        total = st.get("total")
        types = {k.split("/", 1)[1]: v for k, v in st.items() if k.startswith("pair_types/")}
        if types:
            vals = {k: pct(v, total) for k, v in types.items()}
        else:
            vals = {"unmapped": pct(st.get("total_unmapped"), total), "duplicates": pct(st.get("total_dups"), total),
                    "non-duplicate": pct(st.get("total_nodups"), total)}
        for k, v in vals.items():
            series.setdefault(k, [None] * len(samples))[i] = v
    series = {k: v for k, v in series.items() if any(x for x in v if x)}
    if not series:
        return None
    fig = go.Figure()
    order = sorted(series, key=lambda k: -sum(x or 0 for x in series[k]))
    for i, k in enumerate(order):
        label = f"{k} ({PAIR_TYPE_NAMES[k]})" if k in PAIR_TYPE_NAMES else k
        fig.add_bar(y=names, x=[x or 0 for x in series[k]], name=label, orientation="h",
                    marker=dict(color=PAIR_TYPE_COLOURS.get(k, PALETTE[i % len(PALETTE)]),
                                line=dict(color="rgba(128,128,128,0.5)", width=0.5)),
                    hovertemplate="%{y}<br>" + k + ": %{x:.2f}% of read pairs<extra></extra>")
    _layout(fig, _bar_height(len(names)), barmode="stack")
    fig.update_xaxes(title="% of total read pairs", range=[0, 100])
    fig.update_yaxes(autorange="reversed")
    return fig


def plot_cis_trans(samples: List[Dict[str, Any]]) -> Optional[go.Figure]:
    """Stacked bars of cis vs trans as % of non-duplicate pairs."""
    names, cis, trans = [], [], []
    for s in samples:
        st = s["stats"]
        c = s["metrics"]["cis_pct"]
        t = pct(st.get("trans"), nodup_count(st))
        if c is None and t is None:
            continue
        names.append(s["name"]); cis.append(c or 0); trans.append(t or 0)
    if not names:
        return None
    fig = go.Figure()
    fig.add_bar(y=names, x=cis, name="cis", orientation="h", marker_color=PALETTE[0],
                hovertemplate="%{y}<br>cis: %{x:.2f}%<extra></extra>")
    fig.add_bar(y=names, x=trans, name="trans", orientation="h", marker_color=PALETTE[1],
                hovertemplate="%{y}<br>trans: %{x:.2f}%<extra></extra>")
    _layout(fig, _bar_height(len(names)), barmode="stack")
    fig.update_xaxes(title="% of non-duplicate pairs", range=[0, 100])
    fig.update_yaxes(autorange="reversed")
    return fig


def plot_distance_bins(samples: List[Dict[str, Any]]) -> Optional[go.Figure]:
    """Stacked bars of cis distance bins (+ trans) as % of non-duplicate pairs."""
    names = [s["name"] for s in samples]
    labels = [b[0] for b in DIST_BINS] + ["trans"]
    cols: Dict[str, List[float]] = {l: [] for l in labels}
    any_val = False
    for s in samples:
        br = cis_distance_breakdown(s["stats"])
        den = nodup_count(s["stats"])
        for l in labels:
            v = pct(br.get(l), den)
            any_val = any_val or v is not None
            cols[l].append(v or 0)
    if not any_val:
        return None
    fig = go.Figure()
    for i, l in enumerate(labels):
        fig.add_bar(y=names, x=cols[l], name=l, orientation="h",
                    marker=dict(color=GREY if l == "trans" else BIN_COLOURS[i],
                                line=dict(color="rgba(128,128,128,0.5)", width=0.5)),
                    hovertemplate="%{y}<br>" + l + ": %{x:.2f}%<extra></extra>")
    _layout(fig, _bar_height(len(names)), barmode="stack")
    fig.update_xaxes(title="% of non-duplicate pairs", range=[0, 100])
    fig.update_yaxes(autorange="reversed")
    return fig


def _add_1kb_marker(fig: go.Figure) -> None:
    fig.add_vline(x=1000, line=dict(color="rgba(128,128,128,0.9)", width=1.5, dash="dash"),
                  annotation_text="1 kb", annotation_position="top")


def plot_contact_decay(samples: List[Dict[str, Any]], style: Dict[str, Dict[str, str]]) -> Optional[go.Figure]:
    """P(s) curves (log-log), all samples overlaid, each normalised to unit pair count; 1 kb marked."""
    fig = go.Figure()
    n = 0
    for s in samples:
        d = s["decay"]
        if d.empty:
            continue
        n += 1
        st = style[s["name"]]
        fig.add_scatter(x=d["s"], y=d["p"], mode="lines+markers", name=s["name"],
                        line=dict(color=st["colour"], dash=st["dash"], width=2), marker=dict(size=4),
                        hovertemplate="%{x:,.0f} bp<br>P(s)=%{y:.3e}<extra>" + esc(s["name"]) + "</extra>")
    if not n:
        return None
    _layout(fig, 460)
    _add_1kb_marker(fig)
    fig.update_xaxes(type="log", title="Genomic separation s (bp, bin centre)")
    fig.update_yaxes(type="log", title="P(s): pairs per bp, normalised", exponentformat="power")
    return fig


def plot_orientation(samples: List[Dict[str, Any]], style: Dict[str, Dict[str, str]]) -> Optional[go.Figure]:
    """In / out / same-strand pair fractions vs distance, relative to the ideal 25 / 25 / 50 %.

    Values near 1 mean balanced orientations (as for random ligation); a
    short-range excess of inward (or outward) pairs flags digestion / ligation artefacts.
    """
    expect = {"in": 0.25, "out": 0.25, "same": 0.5}
    dashes = {"in": "solid", "out": "dash", "same": "dot"}
    fig = go.Figure()
    n = 0
    for s in samples:
        d = s["decay"]
        if d.empty:
            continue
        n += 1
        st = style[s["name"]]
        for k in ("in", "out", "same"):
            fig.add_scatter(x=d["s"], y=d[f"{k}_frac"] / expect[k], mode="lines", name=f"{s['name']} – {k}",
                            legendgroup=s["name"], line=dict(color=st["colour"], dash=dashes[k], width=2),
                            hovertemplate="%{x:,.0f} bp<br>" + k + " / ideal = %{y:.2f}<extra>" + esc(s["name"]) + "</extra>")
    if not n:
        return None
    _layout(fig, 420)
    fig.add_hline(y=1, line=dict(color="rgba(128,128,128,0.8)", width=1))
    _add_1kb_marker(fig)
    fig.update_xaxes(type="log", title="Genomic separation s (bp, bin centre)")
    fig.update_yaxes(title="Observed fraction / ideal (in 25%, out 25%, same 50%)", rangemode="tozero")
    return fig


def _rgba(hex_colour: str, alpha: float) -> str:
    r, g, b = (int(hex_colour[i:i + 2], 16) for i in (1, 3, 5))
    return f"rgba({r},{g},{b},{alpha})"


def plot_library_complexity(samples: List[Dict[str, Any]], style: Dict[str, Dict[str, str]]) -> Optional[go.Figure]:
    """preseq lc_extrap curves (expected distinct vs total reads, with 95% CI), in millions."""
    fig = go.Figure()
    n = 0
    for s in samples:
        p = s["preseq"]
        if p.empty:
            continue
        n += 1
        st = style[s["name"]]
        x = p["reads"] / 1e6
        fig.add_scatter(x=pd.concat([x, x[::-1]]), y=pd.concat([p["upper"], p["lower"][::-1]]) / 1e6, fill="toself",
                        fillcolor=_rgba(st["colour"], 0.15), line=dict(width=0), hoverinfo="skip",
                        showlegend=False, legendgroup=s["name"])
        fig.add_scatter(x=x, y=p["expected"] / 1e6, mode="lines", name=s["name"], legendgroup=s["name"],
                        line=dict(color=st["colour"], dash=st["dash"], width=2),
                        hovertemplate="%{x:,.0f} M reads<br>%{y:,.1f} M distinct<extra>" + esc(s["name"]) + "</extra>")
        tm, nd = s["stats"].get("total_mapped"), s["stats"].get("total_nodups")
        if tm and nd:
            fig.add_scatter(x=[tm / 1e6], y=[nd / 1e6], mode="markers", legendgroup=s["name"], showlegend=False,
                            marker=dict(color=st["colour"], size=11, symbol="diamond", line=dict(color="white", width=1.5)),
                            hovertemplate="observed (pairtools): %{x:,.1f} M mapped, %{y:,.1f} M unique<extra>" + esc(s["name"]) + "</extra>")
    if not n:
        return None
    _layout(fig, 440)
    fig.update_xaxes(title="Total read pairs sequenced (millions)")
    fig.update_yaxes(title="Expected distinct pairs (millions)")
    return fig


# --------------------------------------------------------------------------
# HTML assembly
# --------------------------------------------------------------------------

CSS = """
:root{--bg:#fff;--fg:#1b1f24;--muted:#5b6570;--card:#f6f8fa;--border:#d0d7de;--accent:#0072B2;
--pass:#007a5a;--warn:#8a5a00;--fail:#b3470a;--na:#5b6570;--passbg:#d9f2ea;--warnbg:#fbeccb;--failbg:#f9dccb;--nabg:#e8e8e8}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#0f1318;--fg:#e6e9ed;--muted:#9aa5b1;--card:#171c23;
--border:#2f3741;--accent:#56B4E9;--pass:#3ed6a5;--warn:#f4c04d;--fail:#ff9a5c;--na:#9aa5b1;
--passbg:#0b3a2d;--warnbg:#4a3a0a;--failbg:#4d230c;--nabg:#2a2f36}}
:root[data-theme="dark"]{--bg:#0f1318;--fg:#e6e9ed;--muted:#9aa5b1;--card:#171c23;--border:#2f3741;--accent:#56B4E9;
--pass:#3ed6a5;--warn:#f4c04d;--fail:#ff9a5c;--na:#9aa5b1;--passbg:#0b3a2d;--warnbg:#4a3a0a;--failbg:#4d230c;--nabg:#2a2f36}
*{box-sizing:border-box}html{scroll-behavior:smooth}
body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif}
.wrap{max-width:1250px;margin:0 auto;padding:0 16px 64px}
header{padding:24px 0 8px;display:flex;flex-wrap:wrap;gap:16px;align-items:center;justify-content:space-between}
h1{font-size:1.6rem;margin:0}h2{font-size:1.25rem;margin:40px 0 8px;border-bottom:1px solid var(--border);padding-bottom:6px}
h3{font-size:1.02rem;margin:22px 0 4px}
.sub{color:var(--muted);font-size:.9rem}
nav{position:sticky;top:0;z-index:5;background:var(--bg);border-bottom:1px solid var(--border);padding:8px 0;display:flex;flex-wrap:wrap;gap:4px 16px;font-size:.9rem}
nav a{color:var(--accent);text-decoration:none}nav a:hover{text-decoration:underline}
.badge{display:inline-flex;align-items:center;gap:6px;border-radius:999px;padding:2px 12px;font-weight:600;font-size:.85rem;border:1px solid currentColor;white-space:nowrap}
.badge .sym{font-weight:800}
.pass{color:var(--pass);background:var(--passbg)}.warn{color:var(--warn);background:var(--warnbg)}
.fail{color:var(--fail);background:var(--failbg)}.na{color:var(--na);background:var(--nabg)}
.overall{font-size:1.1rem;padding:6px 18px}
.cards{display:flex;flex-wrap:wrap;gap:10px;margin:8px 0}
.card{background:var(--card);border:1px solid var(--border);border-radius:8px;padding:8px 14px;min-width:110px}
.card b{display:block;font-size:1.4rem}
.tablewrap{overflow-x:auto;border:1px solid var(--border);border-radius:8px}
table{border-collapse:collapse;width:100%;font-size:.9rem}
th,td{padding:7px 10px;text-align:right;border-bottom:1px solid var(--border);white-space:nowrap}
th:first-child,td:first-child{text-align:left;position:sticky;left:0;background:var(--bg)}
thead th{background:var(--card)}
td.l{text-align:left}
th button{all:unset;cursor:pointer;font-weight:600;display:inline-flex;gap:4px;align-items:center}
th button:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
th[aria-sort="ascending"] button::after{content:"\\25B2";font-size:.7em}
th[aria-sort="descending"] button::after{content:"\\25BC";font-size:.7em}
td.cell-pass{background:var(--passbg)}td.cell-warn{background:var(--warnbg)}td.cell-fail{background:var(--failbg)}
td .tag{font-size:.75em;font-weight:700;margin-left:4px}
tr:last-child td{border-bottom:none}
.note{background:var(--card);border-left:4px solid var(--na);padding:8px 12px;border-radius:4px;color:var(--muted)}
.legend span{margin-right:10px}
code{background:var(--card);padding:1px 5px;border-radius:4px}
.plot{margin:6px 0}
dl{display:grid;grid-template-columns:max-content 1fr;gap:4px 16px}dt{font-weight:600}dd{margin:0}
@media (max-width:640px){dl{grid-template-columns:1fr}h1{font-size:1.3rem}}
.js-plotly-plot .plotly text{fill:var(--fg)!important}
.js-plotly-plot .plotly .xtick text,.js-plotly-plot .plotly .ytick text,.js-plotly-plot .plotly .annotation-text{fill:var(--muted)!important}
"""

JS = """
(function(){
  document.querySelectorAll('table.sortable').forEach(function(tb){
    var ths = tb.querySelectorAll('thead th');
    ths.forEach(function(th, idx){
      var btn = th.querySelector('button'); if(!btn) return;
      btn.addEventListener('click', function(){
        var asc = th.getAttribute('aria-sort') !== 'ascending';
        ths.forEach(function(o){o.removeAttribute('aria-sort');});
        th.setAttribute('aria-sort', asc ? 'ascending' : 'descending');
        var body = tb.tBodies[0], rows = Array.prototype.slice.call(body.rows);
        rows.sort(function(a, b){
          var x = a.cells[idx].getAttribute('data-sort'), y = b.cells[idx].getAttribute('data-sort');
          var nx = parseFloat(x), ny = parseFloat(y);
          if(isNaN(nx) && !isNaN(ny)) return 1;
          if(!isNaN(nx) && isNaN(ny)) return -1;
          var r = (!isNaN(nx) && !isNaN(ny)) ? nx - ny : String(x).localeCompare(String(y));
          return asc ? r : -r;
        });
        rows.forEach(function(r){body.appendChild(r);});
      });
    });
  });
})();
"""


def _badge(status: str, text: Optional[str] = None, cls: str = "") -> str:
    word = text or {"pass": "PASS", "warn": "WARN", "fail": "FAIL", "na": "N/A"}[status]
    return (f'<span class="badge {status} {cls}"><span class="sym" aria-hidden="true">{STATUS_SYMBOL[status]}</span>'
            f"{esc(word)}</span>")


def _threshold_text(rule: Dict[str, Any]) -> str:
    higher = rule.get("direction", "higher_is_better") == "higher_is_better"
    return f"pass {'>=' if higher else '<='} {rule.get('warn')}, fail {'<' if higher else '>'} {rule.get('fail')}"


def render_summary_table(samples, statuses, thresholds) -> str:
    """Sortable per-sample summary table with traffic-light cells (colour + text tag)."""
    metric_cols = list(thresholds)
    btn = '<th scope="col"{title}><button type="button">{label}</button></th>'
    heads = [btn.format(title="", label="Sample"), btn.format(title="", label="Status"),
             btn.format(title=f' title="{esc(EXTRA_INFO["total_million"][2])}"', label="Total pairs (M)")]
    for m in metric_cols:
        info = METRIC_INFO.get(m, {"label": m, "unit": "", "desc": ""})
        unit = f" ({info['unit']})" if info["unit"] else ""
        heads.append(btn.format(title=f' title="{esc(info["desc"])} Thresholds: {esc(_threshold_text(thresholds[m]))}"',
                                label=esc(info["label"] + unit)))
    extras = ("dup_of_mapped_pct", "valid_pairs_pct")
    for k in extras:
        lab, unit, desc = EXTRA_INFO[k]
        heads.append(btn.format(title=f' title="{esc(desc)}"', label=esc(f"{lab} ({unit})")))
    rows = []
    for s in samples:
        m, st = s["metrics"], statuses[s["name"]]
        ov = st["overall"]
        cells = [f'<td data-sort="{esc(s["name"])}"><b>{esc(s["name"])}</b></td>',
                 f'<td data-sort="{STATUS_RANK[ov]}">{_badge(ov)}</td>',
                 f'<td data-sort="{_sortval(m["total_million"])}">{fmt(m["total_million"], 1)}</td>']
        for k in metric_cols:
            v, sv = m.get(k), st[k]
            tag = f'<span class="tag">{STATUS_SYMBOL[sv]} {sv.upper()}</span>' if sv in ("warn", "fail") else ""
            cells.append(f'<td class="cell-{sv}" data-sort="{_sortval(v)}">{fmt(v, 1)}{tag}</td>')
        for k in extras:
            cells.append(f'<td data-sort="{_sortval(m.get(k))}">{fmt(m.get(k), 1)}</td>')
        rows.append("<tr>" + "".join(cells) + "</tr>")
    return ('<div class="tablewrap"><table class="sortable"><thead><tr>' + "".join(heads) + "</tr></thead><tbody>"
            + "".join(rows) + "</tbody></table></div>"
            '<p class="sub legend">Click a column header to sort; hover a header for its definition and thresholds. '
            f'{_badge("pass")}<span></span>{_badge("warn")}<span></span>{_badge("fail")}</p>')


def _hms(sec: Optional[float]) -> str:
    if sec is None or math.isnan(sec):
        return "n/a"
    m, s = divmod(int(round(sec)), 60)
    h, m = divmod(m, 60)
    return f"{h:d}:{m:02d}:{s:02d}"


def _gb(mb: Optional[float]) -> str:
    return "n/a" if mb is None or math.isnan(mb) else f"{mb / 1024:.2f}"


def render_benchmarks(samples) -> Optional[str]:
    """Runtime / peak-memory table (one row per sample and rule); None if no benchmarks."""
    frames = [s["bench"] for s in samples if not s["bench"].empty]
    if not frames:
        return None
    df = pd.concat(frames).sort_values(["sample", "seconds"], ascending=[True, False])
    head = ("<thead><tr>" + "".join(f'<th scope="col"><button type="button">{h}</button></th>' for h in
            ["Sample", "Rule", "Wall time", "Peak RSS (GB)", "I/O in (GB)", "I/O out (GB)"]) + "</tr></thead>")
    body = []
    for r in df.itertuples():
        body.append(f'<tr><td data-sort="{esc(r.sample)}">{esc(r.sample)}</td><td class="l" data-sort="{esc(r.rule)}">{esc(r.rule)}</td>'
                    f'<td data-sort="{_sortval(r.seconds)}">{_hms(r.seconds)}</td>'
                    f'<td data-sort="{_sortval(r.max_rss_mb)}">{_gb(r.max_rss_mb)}</td>'
                    f'<td data-sort="{_sortval(r.io_in_mb)}">{_gb(r.io_in_mb)}</td>'
                    f'<td data-sort="{_sortval(r.io_out_mb)}">{_gb(r.io_out_mb)}</td></tr>')
    return f'<div class="tablewrap"><table class="sortable">{head}<tbody>{"".join(body)}</tbody></table></div>'


def render_matrices(samples) -> Optional[str]:
    """Table of .mcool summaries (only samples where cooler could read the file)."""
    rows = [(s["name"], s["mcool"]) for s in samples if s["mcool"]]
    if not rows:
        return None

    def mil(v: float) -> str:
        return fmt(None if v is None or v != v else v / 1e6, 1)

    body = "".join(
        f"<tr><td>{esc(n)}</td><td>{m['n_resolutions']}</td><td>{fmt(m['finest_bp'], 0)}</td>"
        f"<td>{mil(m['contacts'])}</td><td>{mil(m['nnz'])}</td><td>{'yes' if m['balanced'] else 'no'}</td></tr>"
        for n, m in rows)
    return ('<div class="tablewrap"><table><thead><tr><th>Sample</th><th>Resolutions</th><th>Finest bin (bp)</th>'
            "<th>Contacts (M, finest)</th><th>Non-zero pixels (M)</th><th>Balanced</th></tr></thead>"
            f"<tbody>{body}</tbody></table></div>")


def _fig_html(fig: Optional[go.Figure], key: str, missing: str) -> str:
    if fig is None:
        return f'<p class="note">{esc(missing)}</p>'
    return '<div class="plot">' + fig.to_html(full_html=False, include_plotlyjs=False, div_id=f"plot-{key}",
                                              config={"responsive": True, "displaylogo": False}) + "</div>"


def render_methods(thresholds) -> str:
    """Methods section: metric definitions with the active thresholds, plus plot descriptions."""
    items = []
    for m, rule in thresholds.items():
        info = METRIC_INFO.get(m, {"desc": "(custom metric)"})
        items.append(f"<dt><code>{esc(m)}</code></dt><dd>{esc(info['desc'])} <span class=\"sub\">Thresholds: "
                     f"{esc(_threshold_text(rule))} ({esc(rule.get('direction', 'higher_is_better').replace('_', ' '))}).</span></dd>")
    for k, (_lab, _unit, desc) in EXTRA_INFO.items():
        items.append(f"<dt><code>{esc(k)}</code></dt><dd>{esc(desc)}</dd>")
    return f"""
<p>Everything is computed from the pairtools <code>dedup</code> statistics of each group
(<code>&lt;group&gt;_dedup_stats.txt</code>). &ldquo;Non-duplicate pairs&rdquo; = <code>total_nodups</code>.
Status rule: for <em>higher is better</em> metrics a value below <code>fail</code> fails and below <code>warn</code> warns;
for <em>lower is better</em> metrics a value above <code>fail</code> fails and above <code>warn</code> warns; a value exactly on an
edge is not penalised. A sample's status is its worst metric; the header shows the worst sample. Thresholds come from
<code>qc_thresholds</code> in the config.</p>
<h3>Metric definitions</h3><dl>{''.join(items)}</dl>
<h3>Plots</h3><dl>
<dt>Pair types</dt><dd><code>pair_types/*</code> counts as % of all read pairs (UU unique-unique, UR/RU rescued, NU/UN null-unique, MU multi-unique, DD duplicates).</dd>
<dt>Cis distance bins</dt><dd>Differences of the cumulative <code>cis_Nkb+</code> counters (<code>cis &minus; cis_1kb+</code> for &lt; 1 kb), plus trans, as % of non-duplicate pairs.</dd>
<dt>P(s)</dt><dd>From <code>dist_freq</code>: pairs in each log-spaced distance bin (four strand orientations summed) divided by the bin width in bp and by the total pairs in <code>dist_freq</code>, plotted at the geometric bin centre, so curves are comparable regardless of depth. Micro-C maps show a nucleosome-scale shoulder below ~1 kb (marked). The 0-1 bp bin and any open last bin are not drawn.</dd>
<dt>Orientation</dt><dd>Per distance bin, the fraction of pairs that are inward (<code>+-</code>), outward (<code>-+</code>) or same-strand (<code>++</code>, <code>--</code>), divided by the ideal 25 / 25 / 50 % of random ligation. A strong short-range excess of inward or outward pairs indicates undigested, self-ligated or dangling-end artefacts.</dd>
<dt>Library complexity</dt><dd>preseq <code>lc_extrap</code> expected distinct pairs with 95 % CI. Diamonds show observed (pairtools) mapped vs non-duplicate counts for comparison; preseq counts reads in the BAM, so a small offset is expected.</dd>
</dl>"""


def build_report(samples: List[Dict[str, Any]], thresholds: Dict[str, Dict[str, Any]],
                 title: str = "Micro-C QC report", outdir: str = "") -> str:
    """Assemble the complete single-file HTML report string."""
    style = sample_style(samples)
    statuses = {s["name"]: classify_sample(s["metrics"], thresholds) for s in samples}
    overall = worst_status(st["overall"] for st in statuses.values())
    counts = {k: sum(1 for st in statuses.values() if st["overall"] == k) for k in ("pass", "warn", "fail")}
    nav = [("summary", "Summary"), ("composition", "Library composition"), ("decay", "Contact decay P(s)"),
           ("complexity", "Library complexity"), ("resources", "Runtime & memory"), ("methods", "Methods")]
    sections: List[str] = []

    if not samples:
        body = ('<p class="note">No samples found: no <code>&lt;group&gt;/&lt;group&gt;_dedup_stats.txt</code> files under '
                f"<code>{esc(outdir)}</code>.</p>")
    else:
        cards = "".join(f'<div class="card"><b>{counts[k]}</b>{_badge(k)}</div>' for k in ("pass", "warn", "fail"))
        body = (f'<div class="cards">{cards}<div class="card"><b>{len(samples)}</b>sample(s)</div></div>'
                + render_summary_table(samples, statuses, thresholds))
    sections.append(f'<section id="summary"><h2>Summary</h2>{body}</section>')

    mat = render_matrices(samples)
    sections.append(
        '<section id="composition"><h2>Library composition</h2><h3>Pair types</h3>'
        + _fig_html(plot_pair_types(samples), "pairtypes", "No pair-type information available.")
        + "<h3>Cis vs trans</h3>" + _fig_html(plot_cis_trans(samples), "cistrans", "No cis/trans counts available.")
        + "<h3>Cis by genomic distance</h3>"
        + _fig_html(plot_distance_bins(samples), "distbins", "No cis distance-bin counts available.")
        + (("<h3>Contact matrices (.mcool)</h3>" + mat) if mat else "") + "</section>")
    sections.append(
        '<section id="decay"><h2>Contact decay P(s)</h2><p class="sub">Normalised per sample. Look for the '
        "nucleosome-scale shoulder below ~1 kb (dashed line) and a smooth power-law decay beyond it.</p>"
        + _fig_html(plot_contact_decay(samples, style), "decay", "No dist_freq data in the stats files.")
        + "<h3>Strand orientation vs distance</h3>"
        + _fig_html(plot_orientation(samples, style), "orient", "No strand-resolved dist_freq data available.")
        + "</section>")
    sections.append(
        '<section id="complexity"><h2>Library complexity</h2>'
        + _fig_html(plot_library_complexity(samples, style), "preseq",
                    "No preseq output found (<group>_libraryComplexity.txt); enable library_complexity in the config.")
        + "</section>")
    bench = render_benchmarks(samples)
    sections.append('<section id="resources"><h2>Runtime &amp; memory</h2>'
                    + (bench or '<p class="note">No benchmark files found (benchmarks/*.tsv).</p>') + "</section>")
    sections.append(f'<section id="methods"><h2>Methods &amp; metric definitions</h2>{render_methods(thresholds)}</section>')

    now = _dt.datetime.now().strftime("%Y-%m-%d %H:%M")
    head_badge = _badge(overall, {"pass": "All samples PASS", "warn": "Some samples WARN",
                                  "fail": "Some samples FAIL", "na": "No data"}[overall], "overall")
    navhtml = "".join(f'<a href="#{i}">{esc(t)}</a>' for i, t in nav)
    where = f" &middot; {esc(outdir)}" if outdir else ""
    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<meta name="color-scheme" content="light dark"><title>{esc(title)}</title>
<style>{CSS}</style>
<script>{get_plotlyjs()}</script></head>
<body><div class="wrap">
<header><div><h1>{esc(title)}</h1><div class="sub">{len(samples)} sample(s) &middot; generated {now}{where}</div></div>
<div role="status" aria-label="Overall QC status">{head_badge}</div></header>
<nav aria-label="Sections">{navhtml}</nav>
{''.join(sections)}
</div><script>{JS}</script></body></html>
"""


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def main(argv: Optional[List[str]] = None) -> int:
    """Command-line entry point; always writes a report (even when no samples are found)."""
    ap = argparse.ArgumentParser(description="Generate a static Micro-C QC dashboard.")
    ap.add_argument("--outdir", required=True, help="pipeline output directory containing one folder per group")
    ap.add_argument("--config", default=None, help="MicroCConfig.yaml (for qc_thresholds)")
    ap.add_argument("--out", required=True, help="output HTML file")
    ap.add_argument("--title", default="Micro-C QC report")
    a = ap.parse_args(argv)
    thresholds = load_thresholds(a.config)
    samples = discover_samples(a.outdir)
    if not samples:
        print(f"qc_report: warning: no samples found under {a.outdir}", file=sys.stderr)
    doc = build_report(samples, thresholds, a.title, a.outdir)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write(doc)
    print(f"qc_report: wrote {a.out} ({len(samples)} sample(s))", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
