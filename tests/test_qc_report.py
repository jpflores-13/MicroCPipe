import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))

import qc_report as q  # noqa: E402
import synth  # noqa: E402

CONFIG = os.path.join(HERE, "..", "config", "MicroCConfig.yaml")


def _write_stats(tmp_path, text):
    p = tmp_path / "s.txt"
    p.write_text(text)
    return str(p)


def test_metric_computation(tmp_path):
    text = ("total\t1000\ntotal_unmapped\t100\ntotal_mapped\t900\ntotal_dups\t300\ntotal_nodups\t600\n"
            "cis\t480\ntrans\t120\ncis_1kb+\t300\ncis_10kb+\t120\n")
    m = q.compute_metrics(q.parse_pairtools_stats(_write_stats(tmp_path, text)))
    assert m["mapped_pct"] == pytest.approx(90)
    assert m["duplicate_pct"] == pytest.approx(30)
    assert m["cis_pct"] == pytest.approx(80)
    assert m["cis_ge1kb_pct"] == pytest.approx(50)
    assert m["cis_lt1kb_pct"] == pytest.approx(30)
    assert m["cis_ge10kb_pct"] == pytest.approx(20)
    assert m["valid_pairs_pct"] == pytest.approx(70)
    assert m["unique_pairs_million"] == pytest.approx(600 / 1e6)


def test_zero_and_missing_counts_do_not_raise():
    m = q.compute_metrics({"total": 0.0, "total_nodups": 0.0, "cis": 0.0})
    assert m["mapped_pct"] is None and m["cis_pct"] is None
    assert q.classify_sample(m, q.DEFAULT_THRESHOLDS)["overall"] in ("na", "fail", "warn", "pass")
    assert q.compute_metrics({})["cis_pct"] is None


@pytest.mark.parametrize("value,direction,expected", [
    (85, "higher_is_better", "pass"), (80, "higher_is_better", "pass"), (79.9, "higher_is_better", "warn"),
    (60, "higher_is_better", "warn"), (59.9, "higher_is_better", "fail"),
    (30, "lower_is_better", "pass"), (40, "lower_is_better", "pass"), (41, "lower_is_better", "warn"),
    (60, "lower_is_better", "warn"), (61, "lower_is_better", "fail"), (None, "lower_is_better", "na"),
])
def test_classify(value, direction, expected):
    warn, fail = (80, 60) if direction == "higher_is_better" else (40, 60)
    assert q.classify(value, {"direction": direction, "warn": warn, "fail": fail}) == expected


def test_overall_is_worst():
    thr = {"a": {"direction": "higher_is_better", "warn": 10, "fail": 5},
           "b": {"direction": "lower_is_better", "warn": 10, "fail": 20}}
    assert q.classify_sample({"a": 20, "b": 1}, thr)["overall"] == "pass"
    assert q.classify_sample({"a": 7, "b": 1}, thr)["overall"] == "warn"
    assert q.classify_sample({"a": 20, "b": 15}, thr)["overall"] == "warn"
    assert q.classify_sample({"a": 7, "b": 50}, thr)["overall"] == "fail"


def test_thresholds_loaded_from_repo_config():
    thr = q.load_thresholds(CONFIG)
    assert thr["cis_pct"] == {"direction": "higher_is_better", "warn": 70, "fail": 50}
    assert set(thr) == set(q.DEFAULT_THRESHOLDS)
    assert q.load_thresholds("/nonexistent.yaml") == q.DEFAULT_THRESHOLDS


def test_dist_freq_and_decay(tmp_path):
    stats = q.parse_pairtools_stats(_write_stats(tmp_path, synth.make_stats_text(seed=4)))
    dist = q.parse_dist_freq(stats)
    assert {"lo", "hi", "++", "+-", "-+", "--"} <= set(dist.columns)
    d = q.contact_decay(dist)
    assert not d.empty and (d["p"] > 0).all() and (d["s"] > 0).all()
    assert (d[["in_frac", "out_frac", "same_frac"]].sum(axis=1).round(6) == 1).all()
    # shoulder present: P(s) around 150 bp relative to 1/s trend is higher than near 3 kb
    near = d.iloc[(d["s"] - 150).abs().argsort()[:1]]
    far = d.iloc[(d["s"] - 3000).abs().argsort()[:1]]
    assert (near["p"].values * near["s"].values)[0] > (far["p"].values * far["s"].values)[0]
    assert q.contact_decay(q.parse_dist_freq({})).empty


def test_end_to_end_report(tmp_path):
    out = tmp_path / "out"
    synth.make_demo(str(out))
    html_path = out / "qc" / "MicroC_qc_report.html"
    rc = q.main(["--outdir", str(out), "--config", CONFIG, "--out", str(html_path)])
    assert rc == 0
    doc = html_path.read_text()
    for name in ("Proj_HCT116_B1_T1", "Proj_HCT116_B2_T1", "Proj_K562_B1_T1"):
        assert name in doc
    assert "empty_dir" not in doc.split("<section id=\"summary\">")[1].split("</section>")[0]
    for sec in ("summary", "composition", "decay", "complexity", "resources", "methods"):
        assert f'id="{sec}"' in doc
    assert "Plotly.newPlot" in doc and "plot-decay" in doc and "plot-orient" in doc
    assert "FAIL" in doc  # the poor K562 library fails
    assert "<script src" not in doc  # plotly is inlined, no external loads


def test_status_of_demo_samples(tmp_path):
    out = tmp_path / "out"
    synth.make_demo(str(out))
    samples = {s["name"]: s for s in q.discover_samples(str(out))}
    thr = q.load_thresholds(CONFIG)
    assert q.classify_sample(samples["Proj_K562_B1_T1"]["metrics"], thr)["overall"] == "fail"
    assert q.classify_sample(samples["Proj_HCT116_B1_T1"]["metrics"], thr)["overall"] in ("pass", "warn")
    assert set(samples) == {"Proj_HCT116_B1_T1", "Proj_HCT116_B2_T1", "Proj_K562_B1_T1"}


def test_single_sample_without_preseq_or_benchmarks(tmp_path):
    out = tmp_path / "out"
    synth.make_group(str(out), "Solo", preseq=False, bench=False, seed=9)
    html_path = tmp_path / "r.html"
    assert q.main(["--outdir", str(out), "--config", CONFIG, "--out", str(html_path)]) == 0
    doc = html_path.read_text()
    assert "Solo" in doc
    assert "No preseq output found" in doc
    assert "No benchmark files found" in doc
    assert "plot-decay" in doc


def test_missing_config_and_no_samples(tmp_path):
    out = tmp_path / "empty"
    out.mkdir()
    html_path = tmp_path / "x" / "r.html"
    assert q.main(["--outdir", str(out), "--out", str(html_path)]) == 0
    assert "No samples found" in html_path.read_text()


def test_stats_without_dist_freq_or_pair_types(tmp_path):
    out = tmp_path / "out"
    d = out / "G"
    d.mkdir(parents=True)
    (d / "G_dedup_stats.txt").write_text(
        "total\t100\ntotal_unmapped\t10\ntotal_mapped\t90\ntotal_dups\t0\ntotal_nodups\t90\ncis\t0\ntrans\t0\n")
    html_path = tmp_path / "r.html"
    assert q.main(["--outdir", str(out), "--config", CONFIG, "--out", str(html_path)]) == 0
    assert "No dist_freq data" in html_path.read_text()


def test_sample_names_are_escaped(tmp_path):
    out = tmp_path / "out"
    synth.make_group(str(out), "a<b>&c", preseq=False, bench=False, seed=1)
    html_path = tmp_path / "r.html"
    q.main(["--outdir", str(out), "--out", str(html_path)])
    doc = html_path.read_text()
    assert "a&lt;b&gt;&amp;c" in doc
