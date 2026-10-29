#!/usr/bin/env bash
set -euo pipefail

BASE=/media/desk16/tly8529
PARENT="$BASE/mucoromycota_annotation_braker"

BATCHES=(
  07_braker3_batch_v2
  08_braker2_batch_v2
)

command -v tar >/dev/null || {
    echo "ERROR: tar is unavailable"
    exit 1
}

command -v zstd >/dev/null || {
    echo "ERROR: zstd is unavailable"
    exit 1
}

for name in "${BATCHES[@]}"; do
    src="$PARENT/$name"
    final="$BASE/${name}.tar.zst"
    tmp="$BASE/${name}.tar.zst.part"

    [[ -d "$src" ]] || {
        echo "ERROR: source directory is missing: $src"
        exit 1
    }

    [[ ! -e "$final" ]] || {
        echo "ERROR: final archive already exists: $final"
        exit 1
    }

    [[ ! -e "$tmp" ]] || {
        echo "ERROR: temporary archive already exists: $tmp"
        exit 1
    }
done

pack_one() {
    local name="$1"
    local src="$PARENT/$name"
    local final="$BASE/${name}.tar.zst"
    local tmp="$BASE/${name}.tar.zst.part"
    local checksum="$BASE/${name}.tar.zst.sha256"
    local countfile="$BASE/${name}.tar.zst.entries.txt"

    echo
    echo "============================================================"
    echo "Started: $(date --iso-8601=seconds)"
    echo "Batch: $name"
    echo "Source: $src"
    echo "Destination: $final"
    echo "============================================================"

    local source_entries
    source_entries=$(find "$src" -xdev | wc -l)

    echo "Source entries: $source_entries"
    du -xsh "$src"

    cd "$PARENT"

    tar \
      --acls \
      --xattrs \
      --warning=no-file-changed \
      -cf - "$name" |
    zstd -T0 -3 -o "$tmp"

    echo "Compression completed: $(date --iso-8601=seconds)"
    ls -lh "$tmp"

    echo "Testing zstd integrity..."
    zstd -t "$tmp"

    echo "Reading tar contents..."
    local archive_entries
    archive_entries=$(
        zstd -dc "$tmp" |
        tar -tf - |
        wc -l
    )

    {
        echo "source_entries=$source_entries"
        echo "archive_entries=$archive_entries"
    } > "$countfile"

    if [[ "$source_entries" -ne "$archive_entries" ]]; then
        echo "ERROR: entry count mismatch for $name"
        echo "Source entries:  $source_entries"
        echo "Archive entries: $archive_entries"
        exit 1
    fi

    mv "$tmp" "$final"

    cd "$BASE"
    sha256sum "$(basename "$final")" > "$(basename "$checksum")"

    echo "Archive verified successfully: $name"
    ls -lh "$final" "$checksum" "$countfile"
    cat "$countfile"
    cat "$checksum"
    echo "Finished: $(date --iso-8601=seconds)"
}

for name in "${BATCHES[@]}"; do
    pack_one "$name"
done

echo
echo "============================================================"
echo "ALL ARCHIVES COMPLETED SUCCESSFULLY"
echo "Finished: $(date --iso-8601=seconds)"
echo "============================================================"
