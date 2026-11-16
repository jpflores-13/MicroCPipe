# MicroCPipe: Micro-C Processing Pipeline

## What this pipeline does

Inside every cell, DNA is folded into 3D structures, and the folding
determines which parts of the genome touch each other. **Micro-C** measures
those contacts genome-wide. Cells are crosslinked, the chromatin is cut
into single nucleosomes with the enzyme **MNase**, and DNA ends that sit
close together in 3D are ligated and sequenced. Each sequenced read pair is
one observed contact. Because MNase cuts between nucleosomes instead of at a
restriction site, the contacts can be placed with near-nucleosome precision,
which is what makes loops, stripes and domain borders sharp at 1 kb bins and
below.

MicroCPipe turns the raw sequencing output (FASTQ files) into the files you
analyse (contact maps) and tells you whether the library was any good:

- **`.cool` / `.mcool`**: contact matrices for the open2c stack
  (cooler, cooltools) and for Python / R.
- **`.hic`**: the same contacts for Juicebox / Juicer tools.
- **A QC dashboard**: one self-contained HTML page with pass / warn / fail
  traffic lights for every library, so you can judge a run in minutes
  without opening a log.

It follows the layout of the PhanstielLab
[dietJuicer](https://github.com/PhanstielLab/dietJuicer) pipeline: a
**Core** workflow (FASTQ to maps) and a **Merge** workflow (combine
libraries into one deeper map), driven by a samplesheet and a config file,
and run with Snakemake 8 on SLURM or on one machine. The alignment and pair
steps follow the Dovetail Micro-C recipe (`bwa mem -5SP -T0`, then
`pairtools`).

A rule-by-rule description of the code is in
[`CODE_GUIDE.md`](CODE_GUIDE.md); the full documentation book is in
[`docs/`](docs/) (rendered with Quarto, published to GitHub Pages).

### How the pieces fit together

```
  MicroCSamplesheet.txt                          config/MicroCConfig.yaml
  (one row per sequencing run: metadata,         (reference, mergeBy, outputs,
   Read1, Read2, Sequencing_Directory)            resolutions, QC thresholds)
          │                                                │
          └──────────────────────┬─────────────────────────┘
                                 ▼
                    workflows/MicroCPipe.snakefile        (Core)
                                 │
        align   (per sequencing run: bwa mem -> pairtools parse -> sort)
                                 │
        dedup   (per group: pairtools merge -> dedup, writes stats)
                                 │
        split ── sort_bam           (per group: .pairs.gz, sorted BAM)
           │          │
           │          └── library_complexity  (preseq)
           │
           ├── cool ── mcool        (cooler: .cool, then multi-resolution)
           ├── hic                  (juicer_tools pre)
           │
           │   dedup stats ── qc_summary          (per-group QC table)
           │
           └── [analysis_*]         (optional, off by default: cooltools / mustache)

        benchmark_summary ── qc_report   (benchmarking.tsv, dashboard)
                                 │
                                 ▼
     output/{group}/…  +  output/qc/MicroC_qc_report.html

  MicroCMergeSamplesheet.txt  ──►  workflows/MicroCMerge.snakefile  (Merge)
  (pairs.gz from Core)             merge_pairs ── cool ── mcool / hic
```

A **group** is a set of samplesheet rows that become one library. Rows that
share the values of the `mergeBy` columns (`Project`, `Cell_Type`, `Bio_Rep`,
`Tech_Rep` by default) are one group, and all their sequencing runs are
aligned separately, then merged and deduplicated together. Every output
file is named after its group.

---

## Prerequisites

- A SLURM cluster, or one machine (see [Compute sites](#compute-sites))
- `python3` >= 3.11 and `conda` >= 24.7.1, made available by your site's
  `profiles/<site>/env.sh` (a one-time Miniforge install in
  `$HOME/miniforge3` works everywhere)
- Paired-end Micro-C FASTQ files
- A BWA-indexed reference FASTA, and a two-column chromosome sizes file for
  the chromosomes you want in the maps

Everything else (bwa, pairtools, samtools, preseq, cooler, cooltools, Java)
is installed by Snakemake into per-rule conda environments from `envs/` on
first use. The juicer_tools jar ships in `scripts/`.

## Quickstart

1. Clone the workflow:

    ```bash
    git clone https://github.com/jpflores-13/MicroCPipe.git
    cd MicroCPipe
    ```

2. Edit the tab-separated `MicroCSamplesheet.txt`. One row per sequencing
   run (a library sequenced twice has two rows):

    | Name | Project | Cell_Type | Condition | Bio_Rep | Tech_Rep | Seq_Rep | Read1 | Read2 | Sequencing_Directory |
    |------|---------|-----------|-----------|---------|----------|---------|-------|-------|----------------------|
    | GM12878_ctrl_1_1_1 | Example | GM12878 | Ctrl | 1 | 1 | 1 | MicroC_2M_R1.fastq.gz | MicroC_2M_R2.fastq.gz | /path/to/raw/ |

   `Name`, `Read1`, `Read2`, `Sequencing_Directory` and `Seq_Rep` are
   always required, plus every column listed in `mergeBy`. Other columns are
   free-form metadata.

3. Edit `config/MicroCConfig.yaml`. The lines most runs change:

    ```yaml
    samplesheet: 'MicroCSamplesheet.txt'
    mergeBy: ['Project', 'Cell_Type', 'Bio_Rep', 'Tech_Rep']
    fasta: '/path/to/hg38/BWAIndex/genome.fa'
    chromSizes: '/path/to/hg38.chrom.sizes'
    outputs: ['cool', 'mcool', 'hic']
    ```

4. Pick your compute site (see [Compute sites](#compute-sites)). The default
   is `hpcProfile: 'profiles/local'`; on a cluster, set `MICROC_PROFILE`
   once in `~/.bashrc`.

5. Dry run, then submit:

    ```bash
    bash   MicroCPipe.sh -n      # list every job; nothing runs
    sbatch MicroCPipe.sh         # whole workflow (bash MicroCPipe.sh on profiles/local)
    ```

6. Open `output/qc/MicroC_qc_report.html` in a browser, and take the maps
   from `output/<group>/`.

7. Optional: combine libraries into one deeper map. List the Core `.pairs.gz`
   files in `MicroCMergeSamplesheet.txt` and run the Merge workflow:

    ```bash
    sbatch MicroCMerge.sh
    ```

### Dry run first

```bash
bash MicroCPipe.sh -n
```

Prints every job and shell command that would be submitted; nothing runs.
Any extra arguments are passed straight to `snakemake`.

### Email notifications

Set an address and the workflow emails you once when it starts, once when it
finishes, and once if any job fails (not once per job):

```bash
export MICROC_EMAIL=you@university.edu   # e.g. in ~/.bashrc
```

The address can also go in `notify: email:` in the config; the environment
variable keeps it out of a shared repo and wins if both are set. With
neither set, nothing is sent. `notify: events:` picks which of `start`,
`success`, `error` to send, and `notify: command:` is the shell command that
sends one message (default `mail -s {subject} {email}`, with the body on
stdin; use e.g. `sendmail {email}` on clusters without `mail`, putting a
`Subject:` line in the body). The failure email includes the last 40 lines
of the Snakemake log. Dry runs (`-n`) send nothing, and a failed email never
fails the workflow. Check your cluster can send mail with
`echo test | mail -s test $MICROC_EMAIL`.

The launcher job itself also carries `#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80`;
`sbatch --mail-user=$MICROC_EMAIL MicroCPipe.sh` directs SLURM's own mail.

### After a failed run

```bash
bash unlock.sh               # release the Snakemake lock
sbatch MicroCPipe.sh         # --rerun-incomplete is on by default
```

Finished steps are kept; only the failed and missing ones are redone.

### What runs

| Rule | Jobs | What it does |
|------|------|--------------|
| `align` | 1 per samplesheet row | `bwa mem -5SP -T0` -> `pairtools parse` -> `pairtools sort` |
| `dedup` | 1 per group | `pairtools merge` the runs, `pairtools dedup --mark-dups`, writes `{group}_dedup_stats.txt` |
| `split` | 1 per group | `pairtools split`, then `pairtools select` keeps only unique (UU, UR, RU) pairs -> `{group}.pairs.gz`; also an unsorted BAM |
| `sort_bam` | 1 per group (when `make_bam` is on) | Sorted, indexed BAM |
| `qc_summary` | 1 per group | Dovetail-style QC table from the dedup stats |
| `library_complexity` | 1 per group | `preseq lc_extrap` on the BAM |
| `cool` | 1 per group | `cooler cload pairs` at `base_resolution` (1 thread; not balanced) |
| `mcool` | 1 per group | `cooler zoomify` to every resolution, each resolution ICE-balanced if `balance: true` |
| `hic` | 1 per group | `juicer_tools pre` |
| `analysis_*` | 1 per group per entry in `analysis` | Optional cooltools / mustache on the `.mcool` |
| `benchmark_summary` | 1 | Per-rule runtime / memory table |
| `qc_report` | 1 | The static HTML dashboard |

The Merge workflow has one rule of its own, `merge_pairs`, and reuses
`cool`, `mcool`, `hic` and `analysis_*` unchanged.

---

## Compute sites

Everything specific to a compute site lives in one folder under `profiles/`.
The config and samplesheet are the same on every cluster; what differs is
set once per cluster, outside git:

| Setting | How | Example |
|---------|-----|---------|
| Compute profile | `export MICROC_PROFILE=...` in `~/.bashrc` | `profiles/biowulf` |
| Where output goes | `export MICROC_OUTPUT_ROOT=...` in `~/.bashrc` (unset: the folder you launch from) | `/data/$USER/MicroCPipe` |

The profile is chosen by, in order: `--profile` on the command line,
`$MICROC_PROFILE`, then `hpcProfile` in the config (default
`profiles/local`). A relative `outputDir` goes under `$MICROC_OUTPUT_ROOT`
when it is set; an absolute one is used as written.

| Folder | For |
|--------|-----|
| `profiles/biowulf/` | NIH Biowulf (SLURM) |
| `profiles/longleaf/` | UNC Longleaf (SLURM) |
| `profiles/pioneer/` | Case Western Pioneer (SLURM) |
| `profiles/slurm-template/` | Any other SLURM cluster: copy and fill in the `<...>` values |
| `profiles/local/` | One machine with no scheduler (a workstation, or an interactive node) |

Each folder has two files:

- `env.sh`: sourced by the launchers and `unlock.sh`; makes `python3` and
  `conda` available.
- `config.yaml`: the Snakemake profile: scheduler, account, partition, and
  each rule's threads / memory / time.

To set up a new SLURM cluster:

```bash
cp -r profiles/slurm-template profiles/<cluster>
nano profiles/<cluster>/config.yaml   # account, partition
nano profiles/<cluster>/env.sh        # module loads or Miniforge path
echo 'export MICROC_PROFILE=profiles/<cluster>' >> ~/.bashrc
source ~/.bashrc
sbatch MicroCPipe.sh
```

`sinfo -s` lists partitions, `sacctmgr show assoc user=$USER` lists your
accounts. On Biowulf, deep libraries may need more than the default `hic`
memory (see [Troubleshooting](#troubleshooting)).

---

## Output structure

Everything is written under `outputDir` (default `output/`):

```
output/
├── benchmarking.tsv                    # runtime / memory per group and rule
├── qc/
│   ├── MicroC_qc_report.html           # the QC dashboard (open in a browser)
│   └── qc_report.log
├── {group}/
│   ├── {group}.pairs.gz                # unique (UU/UR/RU), deduplicated pairs (Core)
│   ├── {group}.cool                    # at base_resolution (not balanced)
│   ├── {group}.mcool                   # all resolutions (balanced if balance: true)
│   ├── {group}.hic                     # for Juicebox
│   ├── {group}_dedup_stats.txt         # pairtools stats
│   ├── {group}_qc_summary.txt          # Dovetail-style QC table
│   ├── {group}_libraryComplexity.txt   # preseq curve
│   ├── {group}_mapped.PT.bam(.bai)     # sorted BAM (when make_bam is on)
│   ├── analysis/{group}_{name}.*       # only if config 'analysis' is set
│   ├── logs/{rule}.log                 # one log per job
│   └── benchmarks/{group}_{rule}.tsv
└── logs_slurm/                         # SLURM job logs (cluster profiles)
```

Intermediate files (`aligned/{run}.pairsam.gz`, the dedup pairsam, the
unsorted BAM) are temporary and deleted once the steps that need them have
finished. The Merge workflow names its groups `merged_<mergeBy values>` and writes
`output/merged_<group>/merged_<group>.{cool,mcool,hic}` (plus
`merged_<group>_merged.pairs.gz`), so it never collides with Core outputs.

---

## The QC dashboard

`output/qc/MicroC_qc_report.html` is written by `scripts/qc_report.py` as
the last step of the Core workflow. It is a single static file (no server,
nothing to install) that you can open locally or email to a collaborator.
It shows, for every group (one colour per group in the plots):

- **Summary table** (sortable) with traffic lights for mapped %, duplicate %,
  cis %, cis < 1 kb / >= 1 kb / >= 10 kb and unique pairs (millions), and an
  overall status per group (its worst metric)
- **Library composition**: pair types, cis vs trans, cis by genomic distance
  bin, and (if `cooler` is available) a table summarising each `.mcool`
- **Contact decay P(s)** with the 1 kb mark, where the nucleosome-scale
  shoulder of Micro-C shows, and **strand orientation** (inward, outward,
  same-strand) against distance
- **Library complexity**: the preseq extrapolation, with the observed counts
- **Runtime and memory** per rule, read from the benchmark files
- **Methods**: every metric's definition and the active thresholds

Pass / warn / fail cut-offs come from `qc_thresholds` in the config. They
are starting points chosen for this pipeline, not published standards:
edit them to match your protocol. `docs/qc.qmd` explains how to read each
panel and what a failing Micro-C library looks like.

---

## Pipeline files

| File | Purpose |
|------|---------|
| `MicroCPipe.sh` | Launcher for the Core workflow; run this |
| `MicroCMerge.sh` | Launcher for the Merge workflow |
| `unlock.sh` | Release the Snakemake lock after a failed run |
| `MicroCSamplesheet.txt` | One row per sequencing run (Core input) |
| `MicroCMergeSamplesheet.txt` | One row per library to merge (Merge input) |
| `config/MicroCConfig.yaml` | Reference, grouping, outputs, resolutions, QC thresholds |
| `config/requirements.txt` | Snakemake version for the launcher venv |
| `profiles/<site>/config.yaml` | Scheduler settings and per-rule resources for one site |
| `profiles/<site>/env.sh` | Software setup for that site |
| `workflows/MicroCPipe.snakefile` | The Core workflow |
| `workflows/MicroCMerge.snakefile` | The Merge workflow |
| `workflows/rules/common.smk` | Shared helpers (output paths, format targets, validation) |
| `workflows/rules/formats.smk` | `cool`, `mcool`, `hic` and the optional analysis rules |
| `scripts/get_qc.py` | Per-group QC table from pairtools stats (from Dovetail) |
| `scripts/qc_report.py` | Builds the static QC dashboard |
| `scripts/benchmarking.py` | Collects Snakemake benchmarks into one table |
| `scripts/juicer_tools_1.22.01.jar` | juicer_tools, used by `hic` |
| `envs/*.yml` | Conda environments used by the rules |
| `tests/` | pytest and a tiny dry-run fixture (used by CI) |
| `docs/` | The Quarto documentation book |

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---------|--------------|-----|
| `... is missing column(s): [...]` | Samplesheet header typo, or spaces instead of tabs | Headers must be tab-separated and include the required columns |
| `resolution ... is not a multiple of base_resolution` | A `resolutions` entry that `base_resolution` does not divide | Use multiples of `base_resolution` |
| `config 'outputs' has unknown format(s)` | Typo in `outputs` | Use `cool`, `mcool`, `hic` |
| Snakemake lock error | Previous run crashed | `bash unlock.sh`, then resubmit |
| `hic` job killed (`OUT_OF_MEMORY` or Java heap error) | The rule's `mem_mb` is too small for the library | Raise `hic: mem_mb` in the profile; the Java heap follows it (85%) unless `javaMem` is set ([Recipes](docs/recipes.qmd)) |
| `MissingInputException` for FASTQ or `.pairs.gz` | Wrong `Sequencing_Directory` / `pairs` path | Check the paths; relative paths are relative to where you launch |
| Output appears in the folder you launched from | `MICROC_OUTPUT_ROOT` not set in the shell that ran `sbatch` | Add the export to `~/.bashrc`, `source ~/.bashrc`, resubmit |
| Jobs go to the wrong cluster's partitions or account | `MICROC_PROFILE` not set, so the config default was used | `echo $MICROC_PROFILE`; set it or pass `--profile` |
| No emails arrive | `mail` missing or blocked on the cluster | Test `echo test \| mail -s test $MICROC_EMAIL`; change `notify: command:` |
| A job failed | - | Its log is in `output/{group}/logs/{rule}.log` (QC dashboard: `output/qc/qc_report.log`) |

For a description of every rule and script, see
[`CODE_GUIDE.md`](CODE_GUIDE.md).

---

## References

- dietJuicer, <https://github.com/PhanstielLab/dietJuicer>: the pipeline
  layout MicroCPipe follows.
- Dovetail Genomics Micro-C documentation and example scripts,
  <https://github.com/dovetail-genomics/Micro-C>: alignment recipe and the
  `get_qc.py` QC table.
- Open2C pairtools, <https://github.com/open2c/pairtools>; cooler,
  <https://github.com/open2c/cooler>; cooltools,
  <https://github.com/open2c/cooltools>.
- Hsieh et al., *Cell* (2015) and *Mol Cell* (2020); Krietenstein et al.,
  *Mol Cell* (2020): the Micro-C method.
- Full list: `docs/references.qmd`.
