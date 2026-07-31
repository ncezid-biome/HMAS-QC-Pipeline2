# 05a - Workflow: Reads Processing

This is the first of three workflow-stage documents. It covers the **Reads Processing** phase of the pipeline: raw read QC, primer trimming and per-read demultiplexing by adapter, adapter-mismatch filtering, and paired-end read merging. See [05b - Workflow: Amplicon Filtering](05b-workflow-amplicon-filtering.md) and [05c - Workflow: Reporting](05c-workflow-reporting.md) for the remaining phases.

All processes described here are orchestrated from the main `workflow` block in [`hmas2.nf`](../hmas2.nf).

---

## 1. Raw Read Quality Control (FastQC)

### Function

Each sample's raw paired-end FASTQ files are run through [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to generate baseline quality metrics before any trimming occurs. This gives a pre-processing snapshot of read quality, duplication levels, and per-base quality/N-content that later shows up in the MultiQC report's General Stats and FastQC sections.

No user-facing parameters are exposed for this step — it runs with FastQC's defaults against each sample's raw R1/R2 files.

### Technical Details

- Implemented in the `FASTQC` process, defined once in [`modules/fastqc/main.nf`](../modules/fastqc/main.nf) and imported into `hmas2.nf` under the alias `FASTQC_RAW` (`include { FASTQC as FASTQC_RAW } from './modules/fastqc/main.nf'`), since the same module is reusable at other points in a pipeline if needed.
- Input: a `tuple val(sample), path(reads)` — the per-sample channel produced by the `Channel.fromFilePairs(...)` call at the top of the workflow.
- Output: `*_fastqc.{zip,html}`, emitted as `fastqc_results` and published to `${params.final_outdir}/${sample}/fastqc`.
- Runs via a single shell invocation: `fastqc -q $reads`.
- FastQC's zipped output is later fed directly into MultiQC as one of the general log inputs (see [05c - Workflow: Reporting](05c-workflow-reporting.md)); the individual HTML reports are also published per sample for standalone inspection.

---

## 2. Primer Removal (Cutadapt)

### Function

This step trims primer sequences from both ends of each read pair and, in the same pass, tags each read with which primer pair it matched — effectively demultiplexing a highly multiplexed sample down to the primer-pair level. Reads that don't match any primer pair are discarded.

Exposed parameters (set in [`nextflow.config`](../nextflow.config)):

| Parameter | Default | Description |
|---|---|---|
| `params.cutadapt_maxerror` | `0` | Maximum error rate allowed when matching primer sequences (`0` = exact match required). Passed to cutadapt's `-e` flag. |
| `params.cutadapt_minlength` | `1` | Minimum read length retained after primer removal, discarding empty reads (`-m` flag). |
| `params.cutadapt_thread` | `4` | Number of threads used by each individual cutadapt invocation (`-j` flag). |
| `params.cutadapt_long` | `true` | Whether reads may be longer than the amplicon itself. When `true`, cutadapt uses **linked adapters** (`-a`/`-A`) to trim both the primer and any read-through into the opposing primer/adapter. When `false`, only a 5′ anchored primer is trimmed (`-g`/`-G`). Per the pipeline's own note in `nextflow.config`, leaving this `true` is safe even for reads shorter than the amplicon — see the [cutadapt linked adapter recipe](https://cutadapt.readthedocs.io/en/stable/recipes.html#trimming-amplicon-primers-from-paired-end-reads) for background. |
| `params.maxcutadapts` | `4` | Maximum number of concurrent cutadapt processes (Nextflow `maxForks`), independent of the per-process thread count above. |
| `params.save_trimmed` | `true` | If `true`, intermediate trimmed/matched/discarded FASTQ files are published to the output directory for troubleshooting, rather than being kept only in the Nextflow work directory. |

### Technical Details

- Implemented in the `cutadapt` process ([`modules/cutadapt/main.nf`](../modules/cutadapt/main.nf)), which wraps two Python scripts:
  - [`bin/run_cutadapt.py`](../bin/run_cutadapt.py) — constructs and runs the actual `cutadapt` command line.
  - [`bin/filter_mismatched_adapters.py`](../bin/filter_mismatched_adapters.py) — post-processes cutadapt's output (see step 3 below).
- `run_cutadapt.py` reads the primer file via the `utilities.Primers` class ([`bin/utilities.py`](../bin/utilities.py)), which parses the 4-column tab-delimited oligos file into a dict of `{primer_name: [forward_seq, revcomp(reverse_seq)]}`.
- For every primer pair, the script computes the reverse complement of both the forward and reverse primer (using a small IUPAC-aware `revcomp()` function in `utilities.py`) and inserts a pair of `-g`/`-G` (or `-a`/`-A` when `cutadapt_long` is true) arguments into the cutadapt command line — one command handles the entire primer panel in a single cutadapt invocation per sample, rather than looping cutadapt calls per primer.
- Cutadapt is invoked with `--rename={id}  adapter={adapter_name}={sample} {comment}`, which embeds the matched primer name and sample name directly into the read header. This header convention (`adapter=<name>=<sample>`) is relied upon by nearly every downstream script that needs to know which primer a read belongs to (e.g. `bin/utilities.py`'s `revcomp`, `bin/bin_fastq_by_adapter.py`, `bin/make_count_table.py`).
- `--discard-untrimmed` ensures reads with no primer match are dropped rather than passed through unmodified.
- Cutadapt's own JSON stats (`--json=...`) are captured for later ingestion by MultiQC's built-in Cutadapt module (see the `cutadapt_filtered_reads`, `cutadapt_trimmed_sequences_5`, and `cutadapt_trimmed_sequences_3` sections described in [06 - Output](06-output.md)).
- This process has `errorStrategy 'retry'` with `maxRetries 3`, since cutadapt can occasionally fail transiently under resource contention when run at high concurrency (`maxForks params.maxcutadapts`).
- Dependency: [cutadapt](https://cutadapt.readthedocs.io/en/stable/) (pinned to `4.8` in [`bin/hmas.yaml`](../bin/hmas.yaml)).

---

## 3. Adapter Mismatch Filtering

### Function

Because primer removal is performed independently on R1 and R2, it's possible for a read pair to end up with *different* primer assignments on each mate (e.g. R1 matches primer A's forward sequence but R2's mate unexpectedly matches primer B's reverse sequence). This step checks that R1 and R2 agree on which primer pair they belong to, keeping only concordant pairs and setting aside the rest. This runs automatically as part of primer removal — no separate parameters are exposed.

### Technical Details

- Implemented in [`bin/filter_mismatched_adapters.py`](../bin/filter_mismatched_adapters.py), invoked as a second step within the same `cutadapt` process shell block, immediately after `run_cutadapt.py`.
- Both trimmed FASTQ files are read entirely into memory in a single `read()` call each (rather than streaming line-by-line) to minimize round trips on network filesystems — a design choice the script's docstring calls out explicitly, and one that recurs across several of the `bin/` scripts (see also `bin_fastq_by_adapter.py` and `summarize_unique_fasta.py`).
- The adapter name is extracted from each read header via a compiled regex (`adapter=([^=\s]+)=`) matching the naming convention cutadapt was told to write in step 2.
- Concordant pairs (matching adapter on R1 and R2) are written to `<sample>.matched.1.fastq` / `<sample>.matched.2.fastq`; everything else goes to `<sample>.discarded.1.fastq` / `<sample>.discarded.2.fastq`.
- A per-sample summary CSV (`<sample>.adapter_filter.csv`) records total pairs, concordant pairs, discarded pairs, and discard percentage.
- If R1 and R2 have different record counts (a malformed/truncated input), the script processes only the overlapping paired records and emits a warning to stderr rather than failing outright.
- The `.matched.*.fastq` files are what feed into pair merging (step 5 below); the `.discarded.*.fastq` files and the summary CSV are published for troubleshooting when `params.save_trimmed` is `true`.

---

## 4. Split by Adapter (Optional Per-Primer Demultiplexing)

### Function

Optionally splits each sample's matched, concordant FASTQ pairs into **separate FASTQ files per primer/adapter**, rather than leaving all primers interleaved in one file. This is useful if downstream analysis (outside this pipeline) needs per-primer read sets, but is not required for the pipeline's own QC reporting.

Exposed parameter:

| Parameter | Default | Description |
|---|---|---|
| `params.split_by_adapter` | `true` | If `true`, runs this step and publishes per-adapter FASTQ files. If `false`, the process is skipped entirely (via a Nextflow `when` guard). |

### Technical Details

- Implemented in the `split_by_adapter` process ([`modules/local/split_by_adapter.nf`](../modules/local/split_by_adapter.nf)), which wraps [`bin/bin_fastq_by_adapter.py`](../bin/bin_fastq_by_adapter.py).
- The Python script uses a two-pass strategy to avoid holding many file handles open simultaneously:
  1. **Pass 1** scans headers only (from both R1 and R2, independently) to build a `read_id → adapter` index for each mate, then cross-checks concordance between the two indices (adapter set differences, read ID differences, and mismatched assignments for shared IDs are all reported as warnings). Note this is a redundant safety check — mismatches should already have been filtered out in step 3.
  2. **Pass 2** reads each file fully into memory once, bins records into per-adapter lists, then writes one output file at a time (never more than one output handle open concurrently), named `<sample>.<read_num>.<adapter>.fastq`.
- Output files are published under `${params.final_outdir}/${sample}/by_adapter/`.
- This process shares the same `errorStrategy 'retry'` / `maxRetries 3` / `maxForks params.maxcutadapts` resilience settings as the `cutadapt` process, since it operates on the same per-sample granularity.

---

## 5. Pair Merging (PEAR)

### Function

Merges (assembles) the concordant R1/R2 read pairs into a single consensus sequence per read pair using [PEAR](https://cme.h-its.org/exelixis/web/software/pear/) (Paired-End reAd mergeR), and performs an initial quality check during assembly. Pairs that can't be confidently assembled are discarded and counted separately.

Exposed parameters:

| Parameter | Default | Description |
|---|---|---|
| `params.merging_minquality` | `26` | Minimum base quality (PEAR `-q`) used during assembly quality trimming. The pipeline's own config comments recommend leaving this at the default. |
| `params.merging_maxlength` | `325` | Maximum assembled read length (PEAR `-m`). |
| `params.merging_minlength` | `100` | Minimum assembled read length (PEAR `-n`). |
| `params.merging_minoverlap` | `20` | Minimum required overlap between R1 and R2 for assembly (PEAR `-v`). |

`merging_maxlength` and `merging_minoverlap` are the two values the config comments suggest adjusting for your specific amplicon/read length combination; `merging_minquality` is recommended to stay at its default.

### Technical Details

- Implemented in the `pair_merging` process ([`modules/pair_merging/main.nf`](../modules/pair_merging/main.nf)).
- Runs PEAR with `-j !{params.medcpus}` for threading (shared with the general "medium" CPU allocation — see the CPU/memory grouping described in [04 - Profiles](04-profiles.md)).
- Uses `errorStrategy 'ignore'` — if PEAR fails outright for a given sample (e.g. because one of the input FASTQ files is empty after upstream filtering), the pipeline logs nothing further for that sample at this step but continues processing other samples rather than halting the run. The shell block itself also guards against empty inputs directly, printing a message and skipping the PEAR invocation if either `reads1` or `reads2` is empty.
- PEAR's `.assembled.fastq` output is renamed to `<sample>.fastq` for consistency with downstream naming conventions.
- PEAR's stdout log is parsed by [`bin/parse_pear_log.py`](../bin/parse_pear_log.py), which extracts assembled/discarded/unassembled read counts via simple substring matching on PEAR's fixed log format (e.g. lines starting with `"Assembled reads ...................:"`) and writes them to a small per-sample CSV. These per-sample CSVs are later concatenated across the whole run by the generic `combine_logs` process (see [05c - Workflow: Reporting](05c-workflow-reporting.md)) to populate MultiQC's `pear` custom-content section.
- A `versions.yml` snippet is emitted recording the installed PEAR version (parsed from `pear`'s own stderr banner), which is aggregated into the pipeline's overall software-versions report shown in MultiQC.
- Dependency: [PEAR](https://cme.h-its.org/exelixis/web/software/pear/) (installed via conda as `pear` in [`bin/hmas.yaml`](../bin/hmas.yaml); no version is pinned in the conda spec).

---

## Summary of Data Flow

```
raw R1/R2.fastq.gz
   ├──> FastQC (QC only, no filtering)
   └──> cutadapt (primer trim + per-read adapter tagging)
          └──> filter_mismatched_adapters (R1/R2 concordance check)
                 ├──> [optional] split_by_adapter (per-primer FASTQ files)
                 └──> pear (pair merging)
                        └──> <sample>.fastq  (passed to Amplicon Filtering, see 05b)
```

Continue to [05b - Workflow: Amplicon Filtering](05b-workflow-amplicon-filtering.md) for quality filtering, dereplication, and denoising.
