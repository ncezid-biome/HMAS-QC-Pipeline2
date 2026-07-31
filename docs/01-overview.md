# Step-mothur: Overview

## What It Does

**Step-mothur** (HMAS-QC-Pipeline2) is a [Nextflow](https://www.nextflow.io/docs/latest/index.html)-based bioinformatics pipeline for quality control and denoising of Highly Multiplexed Amplicon Sequencing (HMAS) data. It processes paired-end Illumina reads from highly multiplexed primer panels (thousands of primer pairs targeting a genome, e.g. *Salmonella*) and produces:

- A set of high-quality, denoised, unique representative sequences per sample (FASTA)
- Per-sample and combined run-level reports summarizing primer pair performance and read depth (CSV)
- An aggregated [MultiQC](https://multiqc.info/) HTML report for visual QC review across all samples in a run

The pipeline currently supports **paired-end Illumina data only**.

## Pipeline at a Glance

The workflow is organized into three broad phases, each made up of several Nextflow processes:

<p align="center">Reads Processing → Amplicon Filtering → Reporting</p>

1. **Reads Processing** — raw read QC, primer trimming/demultiplexing, and read-pair merging.
2. **Amplicon Filtering** — quality filtering, dereplication (collapsing to unique sequences), and denoising (removing likely sequencing-error artifacts and chimeras).
3. **Reporting** — mapping high-quality reads back to their unique representatives, building per-sample and per-run summary statistics, and generating the MultiQC report.

Each phase is documented in detail in [05a - Workflow: Reads Processing](05a-workflow-reads-processing.md), [05b - Workflow: Amplicon Filtering](05b-workflow-amplicon-filtering.md), and [05c - Workflow: Reporting](05c-workflow-reporting.md), including the specific scripts, Nextflow processes, and `nextflow.config` parameters involved at each step.

## Inputs

| Input | Description | Config parameter |
|---|---|---|
| Paired-end FASTQ files | Gzipped, one pair per sample, matching the pattern `*_R{1,2}*.fastq.gz`, searched recursively | `params.reads` |
| Primer file | Plain text, 4 columns, tab-delimited: an identifier column, forward primer sequence, reverse-complement of the reverse primer sequence, and a primer pair name | `params.primer` |
| Output directory | Where results are written (a versioned/timestamped subdirectory is created underneath) | `params.outdir` |

The primer file format looks like:

```
primer	CACGCATCATTTCGCAAAAGC	AGTACGTTCGGCCTCTTTCAG	OG0001079primerGroup1
```

## Outputs

| Output | Description |
|---|---|
| `<sample>.final.unique.fasta` | High-quality, denoised, chimera-free unique sequences per sample, with abundance embedded in the sequence ID (e.g. `size=551`) |
| `<sample>.csv` | Per-sample summary report (mean read depth, % successful primer pairs, etc.) |
| `report*.csv` | Combined summary report across all samples in the run |
| `genus_primer_stats.csv` | Read counts for a configurable subset of primers (filtered by keyword) across all samples |
| `*.html` (MultiQC) | Aggregated interactive QC report covering every stage of the pipeline |

Full details on output directory structure and file layout are covered in [06 - Output](06-output.md).

## QC Metrics

Two primary QC thresholds are used to assess sample performance:

- **Mean read depth ≥ 30×** per amplicon
- **≥ 90% primer pair success rate** (a primer pair is "successful" if it has at least 2 amplicons mapping in the sample)

These thresholds were derived empirically; see [`test_data/ROC/QC.md`](../test_data/ROC/QC.md) and [06 - Output](06-output.md) for the full rationale, including the false-negative-rate-vs-depth analysis behind the 30× cutoff.

## Technical Notes

- **Workflow engine**: [Nextflow](https://www.nextflow.io/docs/latest/index.html) using DSL2 (`nextflow.enable.dsl=2`). The pipeline was tested against Nextflow 22.10.6 (conda environment) and 25.10.4 (CI); see [02 - Installation](02-installation.md) for version guidance.
- **Execution modes**: local execution via conda-installed tools, or containerized execution via a prebuilt Docker/Singularity image (`docker://jinfinance/hmas2:v1.1`), selectable via profiles — see [04 - Profiles](04-profiles.md).
- **Sample naming and collision handling**: the pipeline groups FASTQ pairs by inferring a sample name from the filename (stripping lane identifiers like `_L001`). If two distinct file paths resolve to the same sample name (a collision), the pipeline appends a numeric suffix (`_2`, `_3`, ...) and physically renames the files on disk to disambiguate them before processing. This logic lives inline in `hmas2.nf` rather than in a separate module.
- **Undetermined reads**: any FASTQ pair whose sample name starts with `undetermined` (case-insensitive) — the standard Illumina bucket for unassigned reads — is filtered out of the channel before processing begins.
- **Logging**: the pipeline writes a running log to `step_mothur_pipeline.log` in the working directory via a custom `logMessage()` closure defined at the top of `hmas2.nf`, independent of Nextflow's own execution trace/log. On completion, an additional summary (succeeded/failed/cached process counts, duration, and any error report) is appended via a `workflow.onComplete` hook.
