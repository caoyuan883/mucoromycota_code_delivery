#!/usr/bin/env bash
set -euo pipefail

BASE=/media/desk16/tly8529
PARENT="$BASE/mucoromycota_annotation_braker"
NAME=06_braker3_batch_v1

SRC="$PARENT/$NAME"
FINAL="$BASE/${NAME}.tar.zst"
TMP="$BASE/${NAME}.tar.zst.part"
CHECKSUM="$BASE/${NAME}.tar.zst.sha256"
COUNTFILE="$BASE/${NAME}.tar.zst.entries.txt"

[[ -d "$SRC" ]] || {
    echo "ERROR: source directory does not exist: $SRC"
    exit 1
}

command -v tar >/dev/null || {
    echo "ERROR: tar is unavailable"
    exit 1
}

command -v zstd >/dev/null || {
    echo "ERROR: zstd is unavailable"
    exit 1
}

[[ ! -e "$FINAL" ]] || {
    echo "ERROR: final archive already exists: $FINAL"
    exit 1
}

[[ ! -e "$TMP" ]] || {
    echo "ERROR: temporary archive already exists: $TMP"
    exit 1
}

SOURCE_ENTRIES=$(find "$SRC" -xdev | wc -l)

echo "Started: $(date --iso-8601=seconds)"
echo "Source: $SRC"
echo "Destination: $FINAL"
echo "Source entries: $SOURCE_ENTRIES"

du -xsh "$SRC"

cd "$PARENT"

tar \
    --acls \
    --xattrs \
    --warning=no-file-changed \
    -cf - "$NAME" |
zstd -T0 -3 -o "$TMP"

echo "Compression completed: $(date --iso-8601=seconds)"

echo "Testing zstd stream..."
zstd -t "$TMP"

echo "Reading tar contents..."
ARCHIVE_ENTRIES=$(
    zstd -dc "$TMP" |
    tar -tf - |
    wc -l
)

{
    echo "source_entries=$SOURCE_ENTRIES"
    echo "archive_entries=$ARCHIVE_ENTRIES"
} > "$COUNTFILE"

if [[ "$SOURCE_ENTRIES" -ne "$ARCHIVE_ENTRIES" ]]; then
    echo "ERROR: entry count mismatch"
    echo "Source entries:  $SOURCE_ENTRIES"
    echo "Archive entries: $ARCHIVE_ENTRIES"
    exit 1
fi

mv "$TMP" "$FINAL"

cd "$BASE"
sha256sum "$(basename "$FINAL")" > "$(basename "$CHECKSUM")"

echo "Archive verified successfully"
ls -lh "$FINAL" "$CHECKSUM" "$COUNTFILE"
cat "$COUNTFILE"
cat "$CHECKSUM"
echo "Finished: $(date --iso-8601=seconds)"
