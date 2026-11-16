## Contact-map formats, shared by core and merge workflows.
## The including workflow must define `get_pairs(wildcards)` -> path of the
## deduplicated, mapped .pairs.gz for {group}.

## .cool at base_resolution
rule cool:
    input:
        pairs = get_pairs,
        chromsizes = config["chromSizes"],
    output:
        cool = OUT + "/{group}/{group}.cool"
    log:
        OUT + "/{group}/logs/{group}_cool.log"
    benchmark:
        OUT + "/{group}/benchmarks/{group}_cool.tsv"
    threads: 1                         # cooler cload pairs is single-process
    conda: ENVS + "/microc_cooler.yml"
    params:
        assembly = config.get("genome", "hg38"),
        binsize = BASE_RES,
    shell:
        """
        cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 --assembly {params.assembly} \
            {input.chromsizes}:{params.binsize} {input.pairs} {output.cool} 2> {log}
        """

## Multi-resolution .mcool (optionally ICE-balanced)
rule mcool:
    input:
        cool = rules.cool.output.cool
    output:
        mcool = OUT + "/{group}/{group}.mcool"
    log:
        OUT + "/{group}/logs/{group}_mcool.log"
    benchmark:
        OUT + "/{group}/benchmarks/{group}_mcool.tsv"
    threads: 8
    conda: ENVS + "/microc_cooler.yml"
    params:
        resolutions = ",".join(str(r) for r in RESOLUTIONS),
        balance = "--balance" if config.get("balance", True) else "",
    shell:
        """
        cooler zoomify {params.balance} -p {threads} -r {params.resolutions} \
            -o {output.mcool} {input.cool} 2> {log}
        """

## Juicebox .hic via juicer_tools pre
rule hic:
    input:
        pairs = get_pairs,
        chromsizes = config["chromSizes"],
    output:
        hic = OUT + "/{group}/{group}.hic"
    log:
        OUT + "/{group}/logs/{group}_hic.log"
    benchmark:
        OUT + "/{group}/benchmarks/{group}_hic.tsv"
    threads: 16
    conda: ENVS + "/microc_hic.yml"
    params:
        jar = config.get("juicerTools", "scripts/juicer_tools_1.22.01.jar"),
        mem = lambda w, resources: config.get("javaMem") or (f"{int(resources.mem_mb * 0.85)}m" if getattr(resources, "mem_mb", None) else "32g"),
        res = ",".join(str(r) for r in HIC_RES),
        tmpdir = OUT + "/{group}/tmp_hic",
    shell:
        """
        mkdir -p {params.tmpdir}
        java -Xmx{params.mem} -Djava.awt.headless=true -Djava.io.tmpdir={params.tmpdir} \
            -jar {params.jar} pre --threads {threads} -r {params.res} \
            {input.pairs} {output.hic} {input.chromsizes} 2> {log}
        rm -rf {params.tmpdir}
        """

## ---- Optional downstream analysis on the .mcool (config['analysis']) ----------
def _mcool(wildcards):
    return f"{OUT}/{wildcards.group}/{wildcards.group}.mcool"


rule analysis_expected:
    input: mcool = _mcool
    output: tsv = OUT + "/{group}/analysis/{group}_expected.tsv", done = touch(OUT + "/{group}/analysis/{group}_expected.done")
    log: OUT + "/{group}/logs/{group}_expected.log"
    threads: 4
    conda: ENVS + "/microc_cooler.yml"
    params: res = config.get("analysis_resolution", {}).get("expected", 10000)
    shell: "cooltools expected-cis -p {threads} -o {output.tsv} {input.mcool}::resolutions/{params.res} 2> {log}"

rule analysis_compartments:
    input: mcool = _mcool
    output: vecs = OUT + "/{group}/analysis/{group}_compartments.cis.vecs.tsv", done = touch(OUT + "/{group}/analysis/{group}_compartments.done")
    log: OUT + "/{group}/logs/{group}_compartments.log"
    threads: 4
    conda: ENVS + "/microc_cooler.yml"
    params: res = config.get("analysis_resolution", {}).get("compartments", 100000),
            prefix = OUT + "/{group}/analysis/{group}_compartments"
    shell: "cooltools eigs-cis -o {params.prefix} {input.mcool}::resolutions/{params.res} 2> {log}"

rule analysis_insulation:
    input: mcool = _mcool
    output: tsv = OUT + "/{group}/analysis/{group}_insulation.tsv", done = touch(OUT + "/{group}/analysis/{group}_insulation.done")
    log: OUT + "/{group}/logs/{group}_insulation.log"
    threads: 4
    conda: ENVS + "/microc_cooler.yml"
    params: res = config.get("analysis_resolution", {}).get("insulation", 10000),
            windows = " ".join(str(w) for w in config.get("insulation_windows", [50000, 100000, 200000]))
    shell: "cooltools insulation -p {threads} -o {output.tsv} {input.mcool}::resolutions/{params.res} {params.windows} 2> {log}"

rule analysis_loops:
    input: mcool = _mcool
    output: tsv = OUT + "/{group}/analysis/{group}_loops.tsv", done = touch(OUT + "/{group}/analysis/{group}_loops.done")
    log: OUT + "/{group}/logs/{group}_loops.log"
    threads: 1
    conda: ENVS + "/microc_cooler.yml"
    params: res = config.get("analysis_resolution", {}).get("loops", 5000)
    shell: "mustache -f {input.mcool} -r {params.res} -o {output.tsv} 2> {log}"
