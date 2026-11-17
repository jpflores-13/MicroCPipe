## MicroCPipe core workflow: Micro-C FASTQ -> .pairs -> .cool / .mcool / .hic + QC.
##
## Per samplesheet row (one sequencing run): bwa mem -> pairtools parse -> sort
## Per group (rows sharing the mergeBy columns): pairtools merge -> dedup -> split
## Then: formats (cool / mcool / hic), QC (stats, preseq), and the static QC dashboard.

import os
import re
import pandas as pd

configfile: "config/MicroCConfig.yaml"

include: "rules/common.smk"

##### Samplesheet #####
samples = pd.read_csv(config["samplesheet"], sep="\t", dtype=str).fillna("")

_required = {"Name", "Read1", "Read2", "Sequencing_Directory", "Seq_Rep"} | set(config["mergeBy"])
_missing = _required - set(samples.columns)
if _missing:
    raise ValueError(f"{config['samplesheet']} is missing column(s): {sorted(_missing)}")

samples["R1"] = [os.path.join(d, f) for d, f in zip(samples["Sequencing_Directory"], samples["Read1"])]
samples["R2"] = [os.path.join(d, f) for d, f in zip(samples["Sequencing_Directory"], samples["Read2"])]
samples["group"] = samples[config["mergeBy"]].agg("_".join, axis=1)
## one alignment job per row; its ID must be unique
samples["run"] = samples["group"] + "_run" + samples.groupby("group").cumcount().add(1).astype(str)

GROUPS = list(dict.fromkeys(samples["group"]))
RUNS = samples.set_index("run")
RUNS_BY_GROUP = samples.groupby("group")["run"].apply(list).to_dict()

wildcard_constraints:
    group = "|".join(map(re.escape, GROUPS)),
    run = "|".join(map(re.escape, samples["run"])),

RUN_NAME = "MicroCPipe " + os.path.splitext(os.path.basename(config["samplesheet"]))[0]

onstart:
    notify("start", f"Started. Groups: {', '.join(GROUPS)}. Outputs: {', '.join(OUTPUTS)}.", RUN_NAME)

onsuccess:
    print("MicroCPipe completed successfully!")
    notify("success", f"Finished successfully.\nResults: {OUT}\n"
                      f"QC dashboard: {OUT}/qc/MicroC_qc_report.html\nTimings: {OUT}/benchmarking.tsv", RUN_NAME)

onerror:
    notify("error", f"A job failed. Per-job logs: {OUT}/<group>/logs\n"
                    f"Snakemake log: {log}\n\nLast lines:\n{log_tail(log)}", RUN_NAME)

##### Targets #####
def qc_targets():
    t = []
    for g in GROUPS:
        t.append(f"{OUT}/{g}/{g}_qc_summary.txt")
        if config.get("library_complexity", True) and config.get("make_bam", True):
            t.append(f"{OUT}/{g}/{g}_libraryComplexity.txt")
        if config.get("make_bam", True):
            t.append(f"{OUT}/{g}/{g}_mapped.PT.bam")
    if config.get("qc_report", True):
        t.append(f"{OUT}/qc/MicroC_qc_report.html")
    return t

rule all:
    input:
        format_targets(GROUPS),
        qc_targets(),
        OUT + "/benchmarking.tsv",

##### Alignment (per run) #####
rule align:
    input:
        r1 = lambda w: RUNS.loc[w.run, "R1"],
        r2 = lambda w: RUNS.loc[w.run, "R2"],
        chromsizes = config["chromSizes"],
    output:
        pairsam = temp(OUT + "/{group}/aligned/{run}.pairsam.gz"),
    log:
        OUT + "/{group}/logs/{run}_align.log"
    benchmark:
        OUT + "/{group}/benchmarks/{run}_align.tsv"
    threads: 16
    conda: ENVS + "/microc_align.yml"
    params:
        fasta = config["fasta"],
        bwa_flags = config.get("bwa_flags", "-5SP -T0"),
        mapq = config.get("min_mapq", 40),
        walks = config.get("walks_policy", "5unique"),
        gap = config.get("max_inter_align_gap", 30),
        tmpdir = OUT + "/{group}/tmp_{run}",
    shell:
        """
        : > {log}
        mkdir -p {params.tmpdir}
        bwa mem {params.bwa_flags} -t {threads} {params.fasta} {input.r1} {input.r2} 2>> {log} | \
        pairtools parse --min-mapq {params.mapq} --walks-policy {params.walks} \
            --max-inter-align-gap {params.gap} --nproc-in {threads} --nproc-out {threads} \
            --chroms-path {input.chromsizes} 2>> {log} | \
        pairtools sort --tmpdir {params.tmpdir} --nproc {threads} -o {output.pairsam} 2>> {log}
        rm -rf {params.tmpdir}
        """

## Merge the group's runs, mark/remove duplicates, collect stats
rule dedup:
    input:
        lambda w: expand(OUT + "/{{group}}/aligned/{run}.pairsam.gz", run=RUNS_BY_GROUP[w.group]),
    output:
        pairsam = temp(OUT + "/{group}/{group}_dedup.pairsam.gz"),
        stats = OUT + "/{group}/{group}_dedup_stats.txt",
    log:
        OUT + "/{group}/logs/{group}_dedup.log"
    benchmark:
        OUT + "/{group}/benchmarks/{group}_dedup.tsv"
    threads: 16
    conda: ENVS + "/microc_align.yml"
    params:
        tmpdir = OUT + "/{group}/tmp_dedup",
    shell:
        """
        : > {log}
        mkdir -p {params.tmpdir}
        pairtools merge --tmpdir {params.tmpdir} --nproc {threads} {input} 2>> {log} | \
        pairtools dedup --nproc-in {threads} --nproc-out {threads} --mark-dups \
            --output-stats {output.stats} --output {output.pairsam} 2>> {log}
        rm -rf {params.tmpdir}
        """

## Split into mapped, deduplicated .pairs.gz (+ BAM for preseq).
## `dedup --mark-dups` keeps duplicates in the stream as pair_type DD, and `split`
## passes unmapped/multimapped pairs through, so `select` keeps only unique/rescued
## pairs (UU, UR, RU) for the contact maps.
rule split:
    input:
        pairsam = rules.dedup.output.pairsam,
    output:
        pairs = OUT + "/{group}/{group}.pairs.gz",
        bam = temp(OUT + "/{group}/{group}_unsorted.bam"),
    log:
        OUT + "/{group}/logs/{group}_split.log"
    benchmark:
        OUT + "/{group}/benchmarks/{group}_split.tsv"
    threads: 16
    conda: ENVS + "/microc_align.yml"
    params:
        all_pairs = OUT + "/{group}/{group}_all.pairs.gz",
    shell:
        """
        : > {log}
        pairtools split --nproc-in {threads} --nproc-out {threads} \
            --output-pairs {params.all_pairs} --output-sam - {input.pairsam} 2>> {log} | \
        samtools view -bS -@ {threads} -o {output.bam} - 2>> {log}
        pairtools select '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU")' \
            --nproc-in {threads} --nproc-out {threads} -o {output.pairs} {params.all_pairs} 2>> {log}
        rm -f {params.all_pairs}
        """

rule sort_bam:
    input:
        bam = rules.split.output.bam,
    output:
        bam = OUT + "/{group}/{group}_mapped.PT.bam",
        bai = OUT + "/{group}/{group}_mapped.PT.bam.bai",
    log:
        OUT + "/{group}/logs/{group}_sort_bam.log"
    benchmark:
        OUT + "/{group}/benchmarks/{group}_sort_bam.tsv"
    threads: 16
    conda: ENVS + "/microc_align.yml"
    params:
        tmp = OUT + "/{group}/tmp_sortbam",
    shell:
        """
        : > {log}
        mkdir -p {params.tmp}
        samtools sort -@ {threads} -T {params.tmp}/s -o {output.bam} {input.bam} 2>> {log}
        samtools index {output.bam} 2>> {log}
        rm -rf {params.tmp}
        """

## ---- QC ----------------------------------------------------------------------
## Scripts call $CONDA_PREFIX/bin/python explicitly: on a cluster the launcher's venv bin/ is
## ahead of the conda env on PATH inside jobs, so a bare `python` is the wrong interpreter.
rule qc_summary:
    input:
        stats = rules.dedup.output.stats,
    output:
        OUT + "/{group}/{group}_qc_summary.txt",
    log:
        OUT + "/{group}/logs/{group}_qc_summary.log"
    benchmark:
        OUT + "/{group}/benchmarks/{group}_qc_summary.tsv"
    threads: 1
    conda: ENVS + "/microc_qc.yml"
    shell:
        "$CONDA_PREFIX/bin/python {workflow.basedir}/../scripts/get_qc.py -p {input.stats} > {output} 2> {log}"

rule library_complexity:
    input:
        bam = rules.sort_bam.output.bam,
    output:
        OUT + "/{group}/{group}_libraryComplexity.txt",
    log:
        OUT + "/{group}/logs/{group}_preseq.log"
    benchmark:
        OUT + "/{group}/benchmarks/{group}_preseq.tsv"
    threads: 1
    conda: ENVS + "/microc_align.yml"
    params:
        extrap = config.get("preseq_extrap", 2.1e9),
    shell:
        """
        ## preseq cannot extrapolate a very shallow library or one with too few duplicates
        ## ("sample not sufficiently deep or duplicates removed"). That is a property of the
        ## data, not a pipeline failure, so it must not block the maps, the dashboard or Merge:
        ## write an empty table (header only; the dashboard shows "no data") and keep going.
        : > {log}
        preseq lc_extrap -bam -pe -extrap {params.extrap} -step 1e7 -seg_len 1000000000 \
            -output {output} {input.bam} 2>> {log} || {{
            echo "WARNING: preseq failed for {wildcards.group} (see {log}); writing an empty table" >&2
            printf 'TOTAL_READS\tEXPECTED_DISTINCT\tLOWER_0.95CI\tUPPER_0.95CI\n' > {output}
        }}
        """

## Per-rule runtime / memory table (scripts/benchmarking.py)
rule benchmark_summary:
    input:
        map_targets(GROUPS),
        qc = [f"{OUT}/{g}/{g}_qc_summary.txt" for g in GROUPS],
        preseq = [f"{OUT}/{g}/{g}_libraryComplexity.txt" for g in GROUPS
                  if config.get("library_complexity", True) and config.get("make_bam", True)],
    output:
        OUT + "/benchmarking.tsv",
    threads: 1
    conda: ENVS + "/microc_qc.yml"
    params:
        outdir = OUT,
    shell:
        "$CONDA_PREFIX/bin/python {workflow.basedir}/../scripts/benchmarking.py {params.outdir} > {output}"

## Static, self-contained QC dashboard (scripts/qc_report.py)
rule qc_report:
    input:
        stats = expand(OUT + "/{group}/{group}_dedup_stats.txt", group=GROUPS),
        extra = [f"{OUT}/{g}/{g}_libraryComplexity.txt" for g in GROUPS
                 if config.get("library_complexity", True) and config.get("make_bam", True)],
        formats = map_targets(GROUPS),
        bench = OUT + "/benchmarking.tsv",
    output:
        OUT + "/qc/MicroC_qc_report.html",
    log:
        OUT + "/qc/qc_report.log"
    threads: 2
    conda: ENVS + "/microc_qc.yml"
    params:
        outdir = OUT,
        cfg = workflow.configfiles[0] if workflow.configfiles else "config/MicroCConfig.yaml",
    shell:
        """
        $CONDA_PREFIX/bin/python {workflow.basedir}/../scripts/qc_report.py --outdir {params.outdir} \
            --config {params.cfg} --out {output} 2> {log}
        """

##### Hooks for rules/formats.smk #####
def get_pairs(wildcards):
    return f"{OUT}/{wildcards.group}/{wildcards.group}.pairs.gz"

include: "rules/formats.smk"
