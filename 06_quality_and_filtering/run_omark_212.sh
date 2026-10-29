#!/usr/bin/env bash

set -uo pipefail

# [release refactor] paths from the environment; OMAmer/OMArk options unchanged.
# historical: PROT_DIR=/home/py/Desktop/Mucoromycota_BRAKER_longest_isoform/02_longest_isoform_proteins
# historical: DB=/home/py/Desktop/omark_db/LUCA.h5
# historical: OUT_DIR=/home/py/Desktop/omark_results_212
PROT_DIR=${PROT_DIR:?set PROT_DIR (pre-filter longest-isoform proteins)}
DB=${OMAMER_DB:?set OMAMER_DB (path to LUCA.h5)}
OUT_DIR=${OUT_DIR:?set OUT_DIR (OMArk output folder)}
THREADS=16
EXPECTED_COUNT=212

OMAMER_DIR="$OUT_DIR/omamer"
OMARK_DIR="$OUT_DIR/omark"
LOG_DIR="$OUT_DIR/logs"
STATUS_DIR="$OUT_DIR/status"
FAILED_DIR="$OUT_DIR/failed"
MANIFEST="$OUT_DIR/run_manifest.tsv"
STATUS_TABLE="$OUT_DIR/batch_status.tsv"

timestamp() {
    date '+%F %T'
}

write_status() {
    local sample="$1"
    local state="$2"
    printf '%s\t%s\n' "$state" "$(timestamp)" \
        > "$STATUS_DIR/${sample}.status"
}

get_status() {
    local sample="$1"
    local status_file="$STATUS_DIR/${sample}.status"

    if [[ -s "$status_file" ]]; then
        awk 'NR==1 {print $1}' "$status_file"
    else
        printf 'NOT_STARTED\n'
    fi
}

update_status_table() {
    local tmp_file="${STATUS_TABLE}.tmp"

    printf "sample\tstatus\tomamer_file\tomark_sum\n" > "$tmp_file"

    local f sample state omamer_file sum_file
    for f in "${FILES[@]}"; do
        sample=$(basename "$f" _longest_isoform.faa)
        state=$(get_status "$sample")
        omamer_file="$OMAMER_DIR/${sample}.omamer"
        sum_file="$OMARK_DIR/${sample}/${sample}.sum"

        printf "%s\t%s\t%s\t%s\n" \
            "$sample" \
            "$state" \
            "$([[ -s "$omamer_file" ]] && echo YES || echo NO)" \
            "$([[ -s "$sum_file" ]] && echo YES || echo NO)" \
            >> "$tmp_file"
    done

    mv "$tmp_file" "$STATUS_TABLE"
}

echo "============================================================"
echo "OMArk batch workflow"
echo "Start time: $(timestamp)"
echo "Host: $(hostname)"
echo "============================================================"

if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "ERROR: Conda环境未激活"
    exit 1
fi

if [[ "${CONDA_DEFAULT_ENV:-}" != "omark" ]]; then
    echo "WARNING: 当前Conda环境是 ${CONDA_DEFAULT_ENV:-UNKNOWN}，不是 omark"
fi

for cmd in omamer omark md5sum; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: 找不到命令: $cmd"
        exit 1
    fi
done

if [[ ! -d "$PROT_DIR" ]]; then
    echo "ERROR: 蛋白目录不存在: $PROT_DIR"
    exit 1
fi

if [[ ! -r "$DB" ]]; then
    echo "ERROR: 数据库不存在或不可读: $DB"
    exit 1
fi

mkdir -p \
    "$OMAMER_DIR" \
    "$OMARK_DIR" \
    "$LOG_DIR" \
    "$STATUS_DIR" \
    "$FAILED_DIR"

exec 9>"$OUT_DIR/.omark_batch.lock"

if command -v flock >/dev/null 2>&1; then
    if ! flock -n 9; then
        echo "ERROR: 另一个OMArk批处理似乎正在运行"
        exit 1
    fi
fi

shopt -s nullglob
FILES=("$PROT_DIR"/*_longest_isoform.faa)
TOTAL=${#FILES[@]}

if [[ "$TOTAL" -ne "$EXPECTED_COUNT" ]]; then
    echo "ERROR: 预期 $EXPECTED_COUNT 个蛋白文件，实际找到 $TOTAL 个"
    exit 1
fi

echo "Protein files: $TOTAL"
echo "OMArk version: $(omark --version)"
echo "OMAmer version: $(omamer --version)"
echo "Threads: $THREADS"
echo "Database: $DB"
echo

{
    echo "Run metadata"
    echo "Generated: $(date --iso-8601=seconds)"
    echo "Host: $(hostname)"
    echo "Conda environment: ${CONDA_DEFAULT_ENV:-UNKNOWN}"
    echo "Conda prefix: ${CONDA_PREFIX:-UNKNOWN}"
    echo "OMArk version: $(omark --version)"
    echo "OMAmer version: $(omamer --version)"
    echo "Threads per OMAmer search: $THREADS"
    echo "Database path: $DB"
    echo "Database MD5: $(md5sum "$DB" | awk '{print $1}')"
    echo
    echo "OMAmer database metadata:"
    omamer info --db "$DB"
} > "$OUT_DIR/run_metadata.txt" 2>&1

printf "sample\tprotein_fasta\tprotein_count\tfile_size_bytes\n" > "$MANIFEST"

for f in "${FILES[@]}"; do
    sample=$(basename "$f" _longest_isoform.faa)
    protein_count=$(grep -c '^>' "$f")
    file_size=$(stat -c '%s' "$f")

    printf "%s\t%s\t%s\t%s\n" \
        "$sample" "$f" "$protein_count" "$file_size" \
        >> "$MANIFEST"
done

update_status_table

SUCCESS_COUNT=0
FAILED_COUNT=0
SKIPPED_COUNT=0
INDEX=0

for FAA in "${FILES[@]}"; do
    INDEX=$((INDEX + 1))

    SAMPLE=$(basename "$FAA" _longest_isoform.faa)
    OMAMER_FILE="$OMAMER_DIR/${SAMPLE}.omamer"
    OMARK_SAMPLE_DIR="$OMARK_DIR/${SAMPLE}"
    SUM_FILE="$OMARK_SAMPLE_DIR/${SAMPLE}.sum"

    OMAMER_LOG="$LOG_DIR/${SAMPLE}.omamer.log"
    OMARK_LOG="$LOG_DIR/${SAMPLE}.omark.log"

    echo
    echo "============================================================"
    echo "[$INDEX/$TOTAL] $SAMPLE"
    echo "Time: $(timestamp)"
    echo "Proteins: $(grep -c '^>' "$FAA")"
    echo "============================================================"

    if [[ -s "$OMAMER_FILE" && -s "$SUM_FILE" ]]; then
        write_status "$SAMPLE" "COMPLETE"
        echo "SKIP: OMAmer和OMArk结果已经完整存在"
        SKIPPED_COUNT=$((SKIPPED_COUNT + 1))
        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
        update_status_table
        continue
    fi

    if [[ ! -s "$FAA" ]]; then
        write_status "$SAMPLE" "INVALID_FASTA"
        echo "ERROR: FASTA不存在或为空: $FAA"
        FAILED_COUNT=$((FAILED_COUNT + 1))
        update_status_table
        continue
    fi

    if ! grep -qm1 '^>' "$FAA"; then
        write_status "$SAMPLE" "INVALID_FASTA"
        echo "ERROR: FASTA没有header: $FAA"
        FAILED_COUNT=$((FAILED_COUNT + 1))
        update_status_table
        continue
    fi

    if [[ ! -s "$OMAMER_FILE" ]]; then
        OMAMER_TMP="${OMAMER_FILE}.tmp.$$"

        write_status "$SAMPLE" "RUNNING_OMAMER"

        echo "OMAmer start: $(timestamp)"
        OMA_START=$(date +%s)

        omamer search \
            --db "$DB" \
            --query "$FAA" \
            --out "$OMAMER_TMP" \
            --nthreads "$THREADS" \
            > "$OMAMER_LOG" 2>&1

        OMAMER_RC=$?
        OMA_END=$(date +%s)
        OMA_SECONDS=$((OMA_END - OMA_START))

        if [[ "$OMAMER_RC" -eq 0 && -s "$OMAMER_TMP" ]]; then
            mv "$OMAMER_TMP" "$OMAMER_FILE"
            write_status "$SAMPLE" "OMAMER_OK"

            echo "OMAmer OK"
            echo "OMAmer seconds: $OMA_SECONDS"
            echo "OMAmer output: $OMAMER_FILE"
        else
            write_status "$SAMPLE" "OMAMER_FAILED"

            echo "ERROR: OMAmer失败"
            echo "Exit code: $OMAMER_RC"
            echo "Log: $OMAMER_LOG"

            if [[ -e "$OMAMER_TMP" ]]; then
                mv "$OMAMER_TMP" \
                    "$FAILED_DIR/${SAMPLE}.omamer.failed.$(date +%Y%m%d_%H%M%S)" \
                    2>/dev/null || true
            fi

            tail -n 30 "$OMAMER_LOG" 2>/dev/null || true

            FAILED_COUNT=$((FAILED_COUNT + 1))
            update_status_table
            continue
        fi
    else
        echo "Reuse existing OMAmer file: $OMAMER_FILE"
        write_status "$SAMPLE" "OMAMER_OK"
    fi

    if [[ -d "$OMARK_SAMPLE_DIR" && ! -s "$SUM_FILE" ]]; then
        ARCHIVE_DIR="$FAILED_DIR/${SAMPLE}_incomplete_$(date +%Y%m%d_%H%M%S)"
        echo "Archiving incomplete OMArk directory:"
        echo "$OMARK_SAMPLE_DIR -> $ARCHIVE_DIR"
        mv "$OMARK_SAMPLE_DIR" "$ARCHIVE_DIR"
    fi

    OMARK_TMP_DIR="$OMARK_DIR/${SAMPLE}.tmp.$$"
    mkdir -p "$OMARK_TMP_DIR"

    write_status "$SAMPLE" "RUNNING_OMARK"

    echo "OMArk start: $(timestamp)"
    OMK_START=$(date +%s)

    omark \
        -f "$OMAMER_FILE" \
        -d "$DB" \
        -o "$OMARK_TMP_DIR" \
        -of "$FAA" \
        -v \
        > "$OMARK_LOG" 2>&1

    OMARK_RC=$?
    OMK_END=$(date +%s)
    OMK_SECONDS=$((OMK_END - OMK_START))

    TMP_SUM="$OMARK_TMP_DIR/${SAMPLE}.sum"

    if [[ "$OMARK_RC" -eq 0 && -s "$TMP_SUM" ]]; then
        if [[ -d "$OMARK_SAMPLE_DIR" ]]; then
            ARCHIVE_DIR="$FAILED_DIR/${SAMPLE}_previous_$(date +%Y%m%d_%H%M%S)"
            mv "$OMARK_SAMPLE_DIR" "$ARCHIVE_DIR"
        fi

        mv "$OMARK_TMP_DIR" "$OMARK_SAMPLE_DIR"
        write_status "$SAMPLE" "COMPLETE"

        echo "OMArk OK"
        echo "OMArk seconds: $OMK_SECONDS"
        echo "Result: $SUM_FILE"

        SUCCESS_COUNT=$((SUCCESS_COUNT + 1))
    else
        write_status "$SAMPLE" "OMARK_FAILED"

        echo "ERROR: OMArk失败"
        echo "Exit code: $OMARK_RC"
        echo "Log: $OMARK_LOG"

        if [[ -d "$OMARK_TMP_DIR" ]]; then
            mv "$OMARK_TMP_DIR" \
                "$FAILED_DIR/${SAMPLE}_omark_failed_$(date +%Y%m%d_%H%M%S)" \
                2>/dev/null || true
        fi

        tail -n 50 "$OMARK_LOG" 2>/dev/null || true

        FAILED_COUNT=$((FAILED_COUNT + 1))
    fi

    update_status_table

    echo "Current completed: $SUCCESS_COUNT/$TOTAL"
    echo "Current failed: $FAILED_COUNT"
done

update_status_table

COMPLETE_FINAL=$(awk -F '\t' 'NR>1 && $2=="COMPLETE" {n++} END {print n+0}' "$STATUS_TABLE")
FAILED_FINAL=$(awk -F '\t' '
    NR>1 && (
        $2=="OMAMER_FAILED" ||
        $2=="OMARK_FAILED" ||
        $2=="INVALID_FASTA"
    ) {n++}
    END {print n+0}
' "$STATUS_TABLE")

echo
echo "============================================================"
echo "Batch finished: $(timestamp)"
echo "Total inputs: $TOTAL"
echo "Complete: $COMPLETE_FINAL"
echo "Failed: $FAILED_FINAL"
echo "Skipped existing: $SKIPPED_COUNT"
echo "Status table: $STATUS_TABLE"
echo "Manifest: $MANIFEST"
echo "============================================================"

if [[ "$FAILED_FINAL" -gt 0 ]]; then
    exit 1
fi

exit 0
