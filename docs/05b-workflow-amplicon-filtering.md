# 05b - Workflow: Amplicon Filtering

This is the second of three workflow-stage documents. It covers the **Amplicon Filtering** phase: quality filtering of merged reads, dereplication (collapsing to unique sequences), and denoising (removing likely sequencing-error artifacts and chimeras). See [05a - Workflow: Reads Processing](05a-workflow-reads-processing.md) for the preceding phase and [05c - Workflow: Reporting](05c-workflow-reporting.md) for what follows.

All processes described here are implemented in [`modules/vsearch/main.nf`](../modules/vsearch/main.nf) and orchestrated from the main `workflow` block in [`hmas2.nf`](../hmas2.nf). Every step in this stage relies on [vsearch](https://github.com/torognes/vsearch) (pinned to `2.22.1` in [`bin/hmas.yaml`](../bin/hmas.yaml)).

---

## 1. Quality Filtering

### Function

Takes the merged reads produced by PEAR (see [05a - Workflow: Reads Processing](05a-workflow-reads-processing.md)) and discards any read whose accumulated expected error, based on the combined Phred quality scores across the entire read, exceeds a fixed threshold. This is a stricter, whole-read quality check that complements PEAR's own assembly-time quality trimming.

No user-facing parameters are exposed for this step — the expected-error threshold (`--fastq_maxee 1`) is currently hardcoded rather than surfaced via `nextflow.config`.

### Technical Details

- Implemented in the `quality_filtering` process, which runs `vsearch --fastx_filter` with `--fastq_maxee 1`, converting the surviving reads from FASTQ to FASTA (`--fastaout`) in the same step — this is the point in the pipeline where quality scores are discarded and only sequence + count information persists downstream.
- vsearch's own log output is parsed by [`bin/parse_qfilter_log.py`](../bin/parse_qfilter_log.py), which scans for a line containing `"sequences kept"` and extracts the kept/discarded counts via simple whitespace splitting. This feeds the `qfilter` custom-content section in MultiQC (via `combine_logs`, described in [05c - Workflow: Reporting](05c-workflow-reporting.md)).
- Immediately after filtering, [`bin/remove_space.py`](../bin/remove_space.py) strips the space that cutadapt's `--rename` step (see [05a - Workflow: Reads Processing](05a-workflow-reads-processing.md)) inserted between the read ID and the `adapter=` field — i.e. turning:
  ```
  >M03235:53:...:16683:11154  adapter=OG0000294primerGroup8=AR_0409
  ```
  into:
  ```
  >M03235:53:...:16683:11154=OG0000294primerGroup8=AR_0409
  ```
  This is necessary because vsearch (and other tools used later, like the clustering step) treat whitespace as a sequence ID delimiter, which would otherwise silently truncate the adapter/sample metadata off of every sequence ID. The script's docstring explicitly notes that a Python regex substitution was chosen over a shell `sed` call for performance reasons on large files.
- A `versions.yml` snippet records the installed vsearch version, parsed from vsearch's own stderr banner.
- Output FASTA files are published to `${params.final_outdir}/${sample}/temp/`.

---

## 2. Dereplication

### Function

Collapses all quality-filtered reads for a sample down to their **unique sequences**, recording how many times each unique sequence occurred (its abundance). From this point forward in the pipeline, every read count is expressed in terms of unique sequences rather than raw reads — a distinction that matters when interpreting later reports and the MultiQC output.

No user-facing parameters are exposed for this step.

### Technical Details

- Implemented in the `dereplication` process, running `vsearch --derep_fulllength` with `--sizeout` (embeds abundance as `;size=N` in each sequence's ID) and `--relabel_keep` (preserves the original sequence ID as a prefix rather than replacing it with a generic label).
- vsearch's log is parsed by [`bin/parse_derep_log.py`](../bin/parse_derep_log.py), which looks for lines containing `"seqs,"` (total sequence count) and `"unique sequences,"` (unique sequence count) to populate the `dereplication` MultiQC section. Per the pipeline's own documentation comments (see [`bin/multiqc_config.yaml`](../bin/multiqc_config.yaml)), a total-to-unique ratio of 5 or larger is generally considered a good sign of adequate read depth for a sample.
- Output: `<sample>.unique.fasta`, published to `${params.final_outdir}/${sample}/temp/`.

---

## 3. Denoising

### Function

Removes unique sequences that are likely artifacts of sequencing error rather than true biological variants, using frequency- and composition-based denoising, then removes any remaining chimeric sequences. Because most true amplicon sequences in an HMAS run are expected to have high read support, it's normal and expected for the *majority* of unique sequences (by count of distinct sequences, not by total read volume) to be removed at this step — most of the actual read volume is concentrated in a smaller number of surviving true sequences.

Exposed parameters (set in [`nextflow.config`](../nextflow.config)):

| Parameter | Default | Description |
|---|---|---|
| `params.denoising_minsize` | `2` | Minimum abundance (read count) a unique sequence must have to be considered during denoising (vsearch `--minsize`). Sequences below this threshold are excluded from clustering entirely. |
| `params.denoising_alpha` | `4` | The UNOISE alpha parameter, controlling how dissimilar a low-abundance sequence must be from a higher-abundance one before it's treated as a distinct true variant rather than an error-derived artifact of it. See the [vsearch UNOISE3 documentation](https://www.drive5.com/usearch/manual/cmd_unoise3.html) for the underlying algorithm description (vsearch's `--cluster_unoise` implements the same UNOISE3 model). |

### Technical Details

- Implemented in the `denoising` process, which delegates the actual clustering work to a helper shell script, [`bin/run_cluster_unoise.sh`](../bin/run_cluster_unoise.sh), rather than calling vsearch directly in the Nextflow process body.
- **Why per-primer clustering, not one global clustering pass**: `run_cluster_unoise.sh` splits the sample's unique FASTA by primer group (using the 4th column of the oligos file — the primer pair name, e.g. `OG0003222primerGroup4` — as an `awk` pattern match against each sequence's header) and runs `vsearch --cluster_unoise` **separately for each primer group**, before concatenating and re-sorting all centroids by abundance. This matters because UNOISE's error model reasons about abundance ratios *relative to other sequences it's being compared against* — pooling all ~2,000+ primer groups into one clustering pass would let a highly abundant sequence from one primer's amplicon mask true low-abundance variants from an unrelated primer's amplicon, since they'd be judged against each other's noise floor despite representing biologically unrelated loci.
- Each primer group's centroids (`--centroids`, with `--sizein --sizeout` to carry abundance through) are clustered independently, then all centroids are concatenated into one file and finally re-sorted by abundance in descending order via `vsearch --sortbysize`.
- The script uses a `mktemp -d` scratch directory scoped to the current working directory (i.e. the Nextflow task's own work directory) with a `trap` to guarantee cleanup on exit, including on early failure.
- If a primer group's subset FASTA is empty (no reads matched that primer for this sample), that group is silently skipped rather than causing an error.
- Back in the `denoising` Nextflow process (not the shell script), a second pass removes chimeric sequences via `vsearch --uchime3_denovo`, producing the final `<sample>.final.unique.fasta`.
- The process also directly computes its own log — rather than parsing a tool's stdout, it counts `>` header lines before and after via `grep -c` and writes `<sample>_denoise_log.csv` inline in the process shell block. This is a different pattern from the other log-producing steps (`pear`, `qfilter`, `dereplication`), which all delegate parsing to a dedicated `bin/parse_*.py` script — denoising's counts are simple enough that shell arithmetic was used directly instead.
- Edge cases handled explicitly in the shell logic: if the input FASTA is missing/empty, counts are reported as `n/a`; if the post-denoising FASTA is missing/empty, the "after" count defaults to `0` rather than failing.
- Output: `<sample>.final.unique.fasta` — this is the pipeline's primary sequence-level deliverable, published directly under `${params.final_outdir}/${sample}/` (not the `temp/` subdirectory used by the two prior steps in this stage), and is the file referenced in the [Overview](01-overview.md)'s description of embedded abundance (`size=551`).

---

## Summary of Data Flow

```
<sample>.fastq  (merged reads, from 05a)
   └──> quality_filtering (vsearch --fastx_filter, maxee=1)
          └──> <sample>.fasta
                 └──> dereplication (vsearch --derep_fulllength)
                        └──> <sample>.unique.fasta
                               └──> denoising (per-primer UNOISE3 + chimera removal)
                                      └──> <sample>.final.unique.fasta  (passed to Reporting, see 05c)
```

Continue to [05c - Workflow: Reporting](05c-workflow-reporting.md) for how these outputs are mapped back to read counts and summarized into per-sample and per-run reports.
