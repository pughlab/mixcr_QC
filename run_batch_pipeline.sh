#!/bin/bash
# ======================================================================
# Batch pipeline: extract_barcodes -> MiXCR align/assemble -> unified QC
# table, run across every paired sample found in a FASTQ directory.
# ======================================================================
set -euo pipefail

# =========================
#   USER PARAMETERS
# =========================
FASTQ_DIR=""                       # directory containing raw paired FASTQs
OUT_DIR=""                         # base output directory (one subfolder per sample)
RAW_R1_PATTERN="_R1.fastq"         # suffix identifying Read 1 files in FASTQ_DIR
                                    # (e.g. "_R1.fastq.gz" if your raw files are gzipped)
RAW_R2_PATTERN="_R2.fastq"         # matching suffix for Read 2 (same pattern, R1->R2)

# --- extract_barcodes.py ---
EXTRACT_SCRIPT="/path/to/extract_barcodes.py"
BPATTERN=""                        # e.g. "NNNNNNNNNNNNT"  (leave blank if using BLIST)
BLIST=""                           # e.g. "/path/to/barcode_list.txt" (leave blank if using BPATTERN)

# --- MiXCR ---
MIXCR_JAR="/cluster/tools/software/centos7/mixcr/3.0.12/mixcr.jar"
SPECIES="hsa"                      # -s flag for mixcr align
MEM="6G"

# --- QC log parsing ---
PARSER_SCRIPT="/path/to/parse_and_merge_qc.py"
# writes align_stats.csv, assemble_stats.csv, and qc_metrics.tsv into OUT_DIR

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

mkdir -p "$OUT_DIR"
timestamp() { date "+%Y-%m-%d %H:%M:%S"; }

# =========================
#   STEP 1: extract_barcodes (per sample)
# =========================
run_extract_barcodes() {
    local sample="$1" r1="$2" r2="$3" extract_dir="$4"
    echo "[$(timestamp)] Running extract_barcodes for ${sample}..."

    local barcode_arg
    if [[ -n "$BPATTERN" ]]; then
        barcode_arg="--bpattern ${BPATTERN}"
    else
        barcode_arg="--blist ${BLIST}"
    fi

    python3 "$EXTRACT_SCRIPT" \
        --read1 "$r1" \
        --read2 "$r2" \
        --outfile "${extract_dir}/${sample}" \
        ${barcode_arg}

    EXTRACTED_R1="${extract_dir}/${sample}_barcode_R1.fastq"
    EXTRACTED_R2="${extract_dir}/${sample}_barcode_R2.fastq"

    if [[ ! -s "$EXTRACTED_R1" || ! -s "$EXTRACTED_R2" ]]; then
        echo "ERROR: extract_barcodes did not produce expected output fastqs for ${sample}." >&2
        exit 1
    fi
    echo "[$(timestamp)] extract_barcodes complete: ${EXTRACTED_R1}, ${EXTRACTED_R2}"
}

# =========================
#   STEP 2: MiXCR align -> assemble (per sample)
# =========================
run_mixcr() {
    local sample="$1" mixcr_dir="$2" log_dir="$3" tmpdir="$4"

    echo "[$(timestamp)] Running MiXCR align for ${sample}..."
    java -Xmx${MEM} -Djava.io.tmpdir="$tmpdir" -jar "$MIXCR_JAR" align \
        -p rna-seq -s ${SPECIES} \
        -OallowPartialAlignments=true \
        -OvParameters.geneFeatureToAlign=VGeneWithP \
        -r "${log_dir}/LOG_ALIGN_${sample}.txt" \
        "$EXTRACTED_R1" "$EXTRACTED_R2" \
        "${mixcr_dir}/${sample}.vdjca"

    echo "[$(timestamp)] Running MiXCR assemblePartial (rescue 1) for ${sample}..."
    java -Xmx${MEM} -Djava.io.tmpdir="$tmpdir" -jar "$MIXCR_JAR" assemblePartial \
        "${mixcr_dir}/${sample}.vdjca" \
        "${mixcr_dir}/rescue1_${sample}.vdjca"

    echo "[$(timestamp)] Running MiXCR assemblePartial (rescue 2) for ${sample}..."
    java -Xmx${MEM} -Djava.io.tmpdir="$tmpdir" -jar "$MIXCR_JAR" assemblePartial \
        "${mixcr_dir}/rescue1_${sample}.vdjca" \
        "${mixcr_dir}/rescue2_${sample}.vdjca"

    echo "[$(timestamp)] Running MiXCR extend for ${sample}..."
    java -Xmx${MEM} -Djava.io.tmpdir="$tmpdir" -jar "$MIXCR_JAR" extend \
        "${mixcr_dir}/rescue2_${sample}.vdjca" \
        "${mixcr_dir}/extended_${sample}.vdjca"

    echo "[$(timestamp)] Running MiXCR assemble for ${sample}..."
    java -Xmx${MEM} -Djava.io.tmpdir="$tmpdir" -jar "$MIXCR_JAR" assemble \
        -r "${log_dir}/LOG_ASSEMBLE_${sample}.txt" \
        "${mixcr_dir}/extended_${sample}.vdjca" \
        "${mixcr_dir}/${sample}.clns"

    echo "[$(timestamp)] Exporting clones for ${sample}..."
    java -Xmx${MEM} -Djava.io.tmpdir="$tmpdir" -jar "$MIXCR_JAR" exportClones \
        -o -t --chains TRA,TRB,TRG,TRD \
        "${mixcr_dir}/${sample}.clns" \
        -f "${mixcr_dir}/${sample}_clones.tsv"

    echo "[$(timestamp)] MiXCR complete for ${sample}."
}

# =========================
#   MAIN: loop over all R1/R2 pairs in FASTQ_DIR
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

    echo "[$(timestamp)] === Processing sample: ${SAMPLE} ==="

    EXTRACT_DIR="${OUT_DIR}/${SAMPLE}/extracted"
    MIXCR_DIR="${OUT_DIR}/${SAMPLE}/mixcr"
    LOG_DIR="${OUT_DIR}/${SAMPLE}/LOG"
    TMPDIR="${OUT_DIR}/scratch/mixcr_tmp_${SAMPLE}"
    mkdir -p "$EXTRACT_DIR" "$MIXCR_DIR" "$LOG_DIR" "$TMPDIR"

    run_extract_barcodes "$SAMPLE" "$R1" "$R2" "$EXTRACT_DIR"
    run_mixcr "$SAMPLE" "$MIXCR_DIR" "$LOG_DIR" "$TMPDIR"

    echo "[$(timestamp)] === Finished sample: ${SAMPLE} ==="
done

# =========================
#   STEP 3: parse + merge QC logs across all samples
# =========================
echo "[$(timestamp)] Parsing and merging LOG files across all samples..."
python3 "$PARSER_SCRIPT" --base_dir "$OUT_DIR" --out_dir "$OUT_DIR"

echo "[$(timestamp)] Batch pipeline complete. QC tables: ${OUT_DIR}/align_stats.csv, ${OUT_DIR}/assemble_stats.csv"
