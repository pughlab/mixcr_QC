#!/bin/bash
# ======================================================================
# Slurm-orchestrated pipeline:
#   for each sample:  extract_barcodes job  --afterok-->  mixcr job
#   once ALL mixcr jobs finish:             QC merge + plot job
#
# Generates one sbatch script per job (written to OUT_DIR/jobs/) and
# submits them with --dependency so Slurm handles the queuing; this
# script itself just builds and submits, it does not run any analysis.
# ======================================================================
set -euo pipefail

# =========================
#   USER PARAMETERS
# =========================
FASTQ_DIR=""                       # directory containing raw paired FASTQs
OUT_DIR=""                         # base output directory (one subfolder per sample)
RAW_R1_PATTERN="_R1.fastq"         # suffix identifying Read 1 files in FASTQ_DIR
RAW_R2_PATTERN="_R2.fastq"         # matching suffix for Read 2 (same pattern, R1->R2)

# --- extract_barcodes.py ---
EXTRACT_SCRIPT="/path/to/extract_barcodes.py"
BPATTERN=""                        # e.g. "NNNNNNNNNNNNT"  (leave blank if using BLIST)
BLIST=""                           # e.g. "/path/to/barcode_list.txt" (leave blank if using BPATTERN)

# --- MiXCR ---
MIXCR_JAR="/cluster/tools/software/centos7/mixcr/3.0.12/mixcr.jar"
SPECIES="hsa"                      # -s flag for mixcr align

# --- QC ---
PARSER_SCRIPT="/path/to/parse_and_merge_qc.py"
PLOTTER_SCRIPT="/path/to/qc_plotter.py"
QC_PLOT_DIR="${OUT_DIR}/qc_plots"  # where qc_plotter.py writes its figures
# parse_and_merge_qc.py writes align_stats.csv / assemble_stats.csv straight
# into OUT_DIR, which qc_plotter.py then reads from there — no extra paths needed.

# --- Slurm ---
SBATCH_ACCOUNT="pughlab"
SBATCH_PARTITION="all"
SBATCH_EXTRACT_TIME="04:00:00"
SBATCH_EXTRACT_MEM="8G"
SBATCH_MIXCR_TIME="10:00:00"
SBATCH_MIXCR_MEM="6G"
SBATCH_QC_TIME="01:00:00"
SBATCH_QC_MEM="4G"

# =========================
#   SETUP / VALIDATION
# =========================
if [[ -z "$FASTQ_DIR" || -z "$OUT_DIR" ]]; then
    echo "ERROR: FASTQ_DIR and OUT_DIR must both be set." >&2
    exit 1
fi
if [[ -z "$BPATTERN" && -z "$BLIST" ]]; then
    echo "ERROR: one of BPATTERN or BLIST must be set." >&2
    exit 1
fi

JOBS_DIR="${OUT_DIR}/jobs"
mkdir -p "$JOBS_DIR" "$QC_PLOT_DIR"

if [[ -n "$BPATTERN" ]]; then
    BARCODE_ARG="--bpattern ${BPATTERN}"
else
    BARCODE_ARG="--blist ${BLIST}"
fi

MIXCR_JOB_IDS=()

# =========================
#   Per-sample: write + submit extract_barcodes job, then mixcr job
# =========================
shopt -s nullglob
r1_files=("${FASTQ_DIR}"/*"${RAW_R1_PATTERN}")
shopt -u nullglob

if [[ ${#r1_files[@]} -eq 0 ]]; then
    echo "ERROR: no files matching *${RAW_R1_PATTERN} found in ${FASTQ_DIR}" >&2
    exit 1
fi

for R1 in "${r1_files[@]}"; do
    SAMPLE=$(basename "$R1" "$RAW_R1_PATTERN")
    R2="${FASTQ_DIR}/${SAMPLE}${RAW_R2_PATTERN}"

    if [[ ! -f "$R2" ]]; then
        echo "WARNING: skipping ${SAMPLE} — missing matching R2 (${R2})" >&2
        continue
    fi

    EXTRACT_DIR="${OUT_DIR}/${SAMPLE}/extracted"
    MIXCR_DIR="${OUT_DIR}/${SAMPLE}/mixcr"
    LOG_DIR="${OUT_DIR}/${SAMPLE}/LOG"
    TMPDIR="${OUT_DIR}/scratch/mixcr_tmp_${SAMPLE}"
    mkdir -p "$EXTRACT_DIR" "$MIXCR_DIR" "$LOG_DIR" "$TMPDIR"

    EXTRACTED_R1="${EXTRACT_DIR}/${SAMPLE}_barcode_R1.fastq"
    EXTRACTED_R2="${EXTRACT_DIR}/${SAMPLE}_barcode_R2.fastq"

    # --- extract_barcodes job ---
    EXTRACT_JOB_FILE="${JOBS_DIR}/${SAMPLE}_extract.slurm"
    cat <<EOF > "$EXTRACT_JOB_FILE"
#!/bin/bash
#SBATCH --job-name=extract_${SAMPLE}
#SBATCH --output=${JOBS_DIR}/${SAMPLE}_extract_%j.out
#SBATCH --account=${SBATCH_ACCOUNT}
#SBATCH --partition=${SBATCH_PARTITION}
#SBATCH --mem=${SBATCH_EXTRACT_MEM}
#SBATCH --time=${SBATCH_EXTRACT_TIME}
set -euo pipefail

echo "Running extract_barcodes for ${SAMPLE}..."
python3 ${EXTRACT_SCRIPT} \\
    --read1 ${R1} \\
    --read2 ${R2} \\
    --outfile ${EXTRACT_DIR}/${SAMPLE} \\
    ${BARCODE_ARG}
echo "Done: ${SAMPLE} extract_barcodes"
EOF

    EXTRACT_JOB_ID=$(sbatch --parsable "$EXTRACT_JOB_FILE")
    echo "Submitted extract_barcodes for ${SAMPLE}: job ${EXTRACT_JOB_ID}"

    # --- mixcr job (queued behind extract_barcodes) ---
    MIXCR_JOB_FILE="${JOBS_DIR}/${SAMPLE}_mixcr.slurm"
    cat <<EOF > "$MIXCR_JOB_FILE"
#!/bin/bash
#SBATCH --job-name=mixcr_${SAMPLE}
#SBATCH --output=${JOBS_DIR}/${SAMPLE}_mixcr_%j.out
#SBATCH --account=${SBATCH_ACCOUNT}
#SBATCH --partition=${SBATCH_PARTITION}
#SBATCH --mem=${SBATCH_MIXCR_MEM}
#SBATCH --time=${SBATCH_MIXCR_TIME}
set -euo pipefail

echo "Running MiXCR align for ${SAMPLE}..."
java -Xmx${SBATCH_MIXCR_MEM} -Djava.io.tmpdir=${TMPDIR} -jar ${MIXCR_JAR} align \\
    -p rna-seq -s ${SPECIES} \\
    -OallowPartialAlignments=true \\
    -OvParameters.geneFeatureToAlign=VGeneWithP \\
    -r ${LOG_DIR}/LOG_ALIGN_${SAMPLE}.txt \\
    ${EXTRACTED_R1} ${EXTRACTED_R2} \\
    ${MIXCR_DIR}/${SAMPLE}.vdjca

echo "Running MiXCR assemblePartial (rescue 1) for ${SAMPLE}..."
java -Xmx${SBATCH_MIXCR_MEM} -Djava.io.tmpdir=${TMPDIR} -jar ${MIXCR_JAR} assemblePartial \\
    ${MIXCR_DIR}/${SAMPLE}.vdjca \\
    ${MIXCR_DIR}/rescue1_${SAMPLE}.vdjca

echo "Running MiXCR assemblePartial (rescue 2) for ${SAMPLE}..."
java -Xmx${SBATCH_MIXCR_MEM} -Djava.io.tmpdir=${TMPDIR} -jar ${MIXCR_JAR} assemblePartial \\
    ${MIXCR_DIR}/rescue1_${SAMPLE}.vdjca \\
    ${MIXCR_DIR}/rescue2_${SAMPLE}.vdjca

echo "Running MiXCR extend for ${SAMPLE}..."
java -Xmx${SBATCH_MIXCR_MEM} -Djava.io.tmpdir=${TMPDIR} -jar ${MIXCR_JAR} extend \\
    ${MIXCR_DIR}/rescue2_${SAMPLE}.vdjca \\
    ${MIXCR_DIR}/extended_${SAMPLE}.vdjca

echo "Running MiXCR assemble for ${SAMPLE}..."
java -Xmx${SBATCH_MIXCR_MEM} -Djava.io.tmpdir=${TMPDIR} -jar ${MIXCR_JAR} assemble \\
    -r ${LOG_DIR}/LOG_ASSEMBLE_${SAMPLE}.txt \\
    ${MIXCR_DIR}/extended_${SAMPLE}.vdjca \\
    ${MIXCR_DIR}/${SAMPLE}.clns

echo "Exporting clones for ${SAMPLE}..."
java -Xmx${SBATCH_MIXCR_MEM} -Djava.io.tmpdir=${TMPDIR} -jar ${MIXCR_JAR} exportClones \\
    -o -t --chains TRA,TRB,TRG,TRD \\
    ${MIXCR_DIR}/${SAMPLE}.clns \\
    -f ${MIXCR_DIR}/${SAMPLE}_clones.tsv

echo "Done: ${SAMPLE} mixcr"
EOF

    MIXCR_JOB_ID=$(sbatch --parsable --dependency=afterok:${EXTRACT_JOB_ID} "$MIXCR_JOB_FILE")
    echo "Submitted mixcr for ${SAMPLE}: job ${MIXCR_JOB_ID} (after ${EXTRACT_JOB_ID})"

    MIXCR_JOB_IDS+=("$MIXCR_JOB_ID")
done

if [[ ${#MIXCR_JOB_IDS[@]} -eq 0 ]]; then
    echo "ERROR: no samples were submitted — nothing to run QC on." >&2
    exit 1
fi

# =========================
#   QC job: parse+merge logs, then plot — queued behind every mixcr job
# =========================
MIXCR_DEPENDENCY=$(IFS=:; echo "${MIXCR_JOB_IDS[*]}")

QC_JOB_FILE="${JOBS_DIR}/qc_merge_plot.slurm"
cat <<EOF > "$QC_JOB_FILE"
#!/bin/bash
#SBATCH --job-name=qc_merge_plot
#SBATCH --output=${JOBS_DIR}/qc_merge_plot_%j.out
#SBATCH --account=${SBATCH_ACCOUNT}
#SBATCH --partition=${SBATCH_PARTITION}
#SBATCH --mem=${SBATCH_QC_MEM}
#SBATCH --time=${SBATCH_QC_TIME}
set -euo pipefail

echo "Parsing and merging LOG files across all samples..."
python3 ${PARSER_SCRIPT} --base_dir ${OUT_DIR} --out_dir ${OUT_DIR}

echo "Running QC plotter..."
python3 ${PLOTTER_SCRIPT} --align_stats ${OUT_DIR}/align_stats.csv --assemble_stats ${OUT_DIR}/assemble_stats.csv --out_dir ${QC_PLOT_DIR}

echo "Done: QC merge + plot"
EOF

QC_JOB_ID=$(sbatch --parsable --dependency=afterok:${MIXCR_DEPENDENCY} "$QC_JOB_FILE")
echo "Submitted QC merge+plot job ${QC_JOB_ID} (after all mixcr jobs: ${MIXCR_DEPENDENCY})"

echo "All jobs submitted. ${#MIXCR_JOB_IDS[@]} sample(s) queued; QC job ${QC_JOB_ID} will run once they all complete."
