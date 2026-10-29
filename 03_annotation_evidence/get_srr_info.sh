#!/bin/bash
set -euo pipefail

INPUT_FILE="sra_availability.csv"
OUTPUT_FILE="sra_detailed_runinfo.csv"

echo "Species,Run_Accession,LibraryLayout,AvgLength,Size_MB" > "$OUTPUT_FILE"

tail -n +2 "$INPUT_FILE" | tr -d '\r' > tmp_detailed.csv

while IFS=',' read -r species count <&3; do
    if [[ "$count" -gt 0 ]]; then
        esearch -db sra -query "\"$species\"[Organism] AND \"illumina\"[Platform] AND \"rna seq\"[Strategy]" | \
        efetch -format runinfo | \
        awk -F',' -v sp="$species" 'NR>1 && $1!="" {print sp","$1","$16","$7","$8}' >> "$OUTPUT_FILE"

        echo "Fetched: $species"
        sleep 1
    fi
done 3< "tmp_detailed.csv"

rm -f tmp_detailed.csv

echo "Done. Detailed run info saved to $OUTPUT_FILE"
