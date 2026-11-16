## MicroCMerge workflow: combine already-processed libraries (.pairs.gz from
## MicroCPipe) into one deeper "mega" map in cool / mcool / hic formats.

import os
import re
import pandas as pd

configfile: "config/MicroCConfig.yaml"

include: "rules/common.smk"

mcfg = config.get("merge", {})
sheet = mcfg.get("samplesheet", "MicroCMergeSamplesheet.txt")
merge_by = mcfg.get("mergeBy", config["mergeBy"])

samples = pd.read_csv(sheet, sep="\t", dtype=str).fillna("")
_missing = ({"pairs"} | set(merge_by)) - set(samples.columns)
if _missing:
    raise ValueError(f"{sheet} is missing column(s): {sorted(_missing)}")

## "merged_" prefix keeps outputs distinct from the core workflow's {group} files
samples["group"] = "merged_" + samples[merge_by].agg("_".join, axis=1)
PAIRS = samples.groupby("group")["pairs"].apply(lambda x: list(dict.fromkeys(x))).to_dict()
GROUPS = list(PAIRS)

wildcard_constraints:
    group = "|".join(map(re.escape, GROUPS)),

RUN_NAME = "MicroCMerge " + os.path.splitext(os.path.basename(sheet))[0]

onstart:
    notify("start", f"Started. Groups: {', '.join(GROUPS)}. Outputs: {', '.join(OUTPUTS)}.", RUN_NAME)

onsuccess:
    print("MicroCMerge completed successfully!")
    notify("success", f"Finished successfully.\nResults: {OUT}", RUN_NAME)

onerror:
    notify("error", f"A job failed. Per-job logs: {OUT}/<group>/logs\n"
                    f"Snakemake log: {log}\n\nLast lines:\n{log_tail(log)}", RUN_NAME)

rule all:
    input:
        format_targets(GROUPS),

## Merge the groups' .pairs.gz (already deduplicated per library)
rule merge_pairs:
    input:
        lambda w: PAIRS[w.group],
    output:
        pairs = OUT + "/{group}/{group}_merged.pairs.gz",
    log:
        OUT + "/{group}/logs/{group}_merge_pairs.log"
    benchmark:
        OUT + "/{group}/benchmarks/{group}_merge_pairs.tsv"
    threads: 8
    conda: ENVS + "/microc_align.yml"
    params:
        tmpdir = OUT + "/{group}/tmp_merge",
    shell:
        """
        mkdir -p {params.tmpdir}
        pairtools merge --tmpdir {params.tmpdir} --nproc {threads} -o {output.pairs} {input} 2> {log}
        rm -rf {params.tmpdir}
        """

## Hook for rules/formats.smk: formats are built from the merged pairs,
## written as output/{group}/{group}.{cool,mcool,hic}
def get_pairs(wildcards):
    return f"{OUT}/{wildcards.group}/{wildcards.group}_merged.pairs.gz"

include: "rules/formats.smk"
