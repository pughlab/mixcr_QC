# TCR Repertoire QC Pipeline

Barcode extraction → MiXCR alignment/assembly → merged QC metrics and plots,
runnable either as a plain batch loop or as dependency-chained Slurm jobs.

## Pipeline overview

```
raw paired FASTQs
      │
      ▼
extract_barcodes.py          (per sample)
      │  -> <sample>_barcode_R1.fastq, <sample>_barcode_R2.fastq
      ▼
MiXCR: align -> assemblePartial (x2) -> extend -> assemble -> exportClones
      │  -> LOG_ALIGN_<sample>.txt, LOG_ASSEMBLE_<sample>.txt, <sample>_clones.tsv
      ▼
parse_and_merge_qc.py         (across all samples)
      │  -> align_stats.csv, assemble_stats.csv, qc_metrics.tsv
      ▼
qc_plotter.py
      -> 6 QC plots + qc_summary_metrics.tsv
```

Each sample's outputs live under `OUT_DIR/<sample>/`:

```
OUT_DIR/
├── <sample>/
│   ├── extracted/   <sample>_barcode_R1.fastq, _R2.fastq, barcode stats
│   ├── mixcr/       .vdjca, .clns, <sample>_clones.tsv
│   └── LOG/         LOG_ALIGN_<sample>.txt, LOG_ASSEMBLE_<sample>.txt
├── align_stats.csv
├── assemble_stats.csv
├── qc_metrics.tsv
└── qc_plots/
```

## Scripts

| Script | Purpose |
|---|---|
| `run_batch_pipeline.sh` | Loops over all FASTQ pairs in a directory, running extract_barcodes + MiXCR + QC parsing sequentially on the local machine / within one job. |
| `submit_pipeline_slurm.sh` | Same pipeline, but submits each step as a Slurm job with `--dependency` chaining: extract → mixcr per sample, then one QC job queued behind every mixcr job. |
| `parse_and_merge_qc.py` | Parses `LOG_ALIGN_*`/`LOG_ASSEMBLE_*` files across all sample folders into `align_stats.csv`, `assemble_stats.csv`, and a merged `qc_metrics.tsv`. |
| `qc_plotter.py` | Reads `align_stats.csv` + `assemble_stats.csv` and produces cohort-level QC plots. |

## Requirements

- Python 3 with `pandas`, `numpy`, `matplotlib`
- Java + [MiXCR](https://github.com/milaboratory/mixcr) (tested against 3.0.12)
- `extract_barcodes.py` (your own barcode-extraction script — not included here)
- Slurm, if using `submit_pipeline_slurm.sh`

## Usage

### Option A — batch loop (single job, all samples sequential)

Edit the parameter block at the top of `run_batch_pipeline.sh`:

```bash
FASTQ_DIR=""          # directory of raw paired FASTQs
OUT_DIR=""            # base output directory
RAW_R1_PATTERN="_R1.fastq"   # e.g. "_R1.fastq.gz" if gzipped
RAW_R2_PATTERN="_R2.fastq"
EXTRACT_SCRIPT="/path/to/extract_barcodes.py"
BPATTERN=""           # or BLIST="/path/to/barcode_list.txt"
MIXCR_JAR="/path/to/mixcr.jar"
PARSER_SCRIPT="/path/to/parse_and_merge_qc.py"
```

Then:

```bash
bash run_batch_pipeline.sh
```

Samples are paired by matching `SAMPLE${RAW_R1_PATTERN}` to
`SAMPLE${RAW_R2_PATTERN}` in `FASTQ_DIR`; anything missing its mate is
skipped with a warning.

### Option B — Slurm-orchestrated (recommended for large cohorts)

Edit the parameter block at the top of `submit_pipeline_slurm.sh`, including
the Slurm section:

```bash
SBATCH_ACCOUNT="your_account"
SBATCH_PARTITION="all"
SBATCH_EXTRACT_TIME="04:00:00"
SBATCH_EXTRACT_MEM="8G"
SBATCH_MIXCR_TIME="10:00:00"
SBATCH_MIXCR_MEM="6G"
SBATCH_QC_TIME="01:00:00"
SBATCH_QC_MEM="4G"
```

Then:

```bash
bash submit_pipeline_slurm.sh
```

This submits, per sample, an `extract_${SAMPLE}` job followed by a
`mixcr_${SAMPLE}` job (`--dependency=afterok:<extract_job_id>`), and finally
one `qc_merge_plot` job that depends on **every** mixcr job finishing
(`--dependency=afterok:<id1>:<id2>:...`). Generated job scripts and their
Slurm `.out` logs are written to `OUT_DIR/jobs/`.

> **Note:** with `afterok`, if any single sample's job fails, its downstream
> job (and the final QC job) will never run. If you'd rather have QC run on
> whichever samples succeeded, switch the QC dependency to `afterany` —
> `parse_and_merge_qc.py` already skips sample folders with missing LOG files.

### Running QC parsing/plotting standalone

```bash
python3 parse_and_merge_qc.py --base_dir OUT_DIR --out_dir OUT_DIR
python3 qc_plotter.py --align_stats OUT_DIR/align_stats.csv \
                       --assemble_stats OUT_DIR/assemble_stats.csv \
                       --out_dir OUT_DIR/qc_plots
```

## QC plots produced

1. `01_sequencing_depth.png` — total sequencing reads per sample
2. `02_alignment_rate.png` — % successfully aligned reads, with threshold line (`--align_rate_threshold`, default 80%)
3. `03_reads_used_in_clonotypes.png` — % of reads landing in final clonotypes, with threshold line (`--reads_used_threshold`, default 50%)
4. `04_clonotype_count.png` — final clonotype count per sample
5. `05_chain_distribution.png` — TRA/TRB/TRG/TRD fraction of clonotype reads per sample
6. `06_clone_size_vs_count.png` — average reads per clonotype vs. final clonotype count (log-log)

Percentages for plots 2 and 3 are recomputed from raw read counts rather than
read from the log files directly — MiXCR log parsing strips the `(xx.x%)`
annotations down to the underlying count, so the percentage is derived as
`count / Total sequencing reads * 100`.

## Notes / known limitations

- `parse_and_merge_qc.py` expects one MiXCR align log and one assemble log
  per sample folder; if either is missing it logs a warning and continues.
- `qc_plotter.py` degrades gracefully (skips the affected plot/column with a
  warning) if an expected column is absent from a log, rather than failing
  the whole run.
- Sample name matching relies on filename conventions: `_barcode_R1.fastq`
  suffix (align stage) and `extended_<sample>.vdjca` (assemble stage) —
  update the regexes in the `sample_from_*` helper functions in
  `qc_plotter.py` if your naming differs.
