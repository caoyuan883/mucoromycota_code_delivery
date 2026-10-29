#!/bin/bash
set -euo pipefail

INPUT_FILE="species_list.txt"
OUTPUT_FILE="sra_availability.csv"

echo "Species,SRA_Count" > "$OUTPUT_FILE"

tr '\r' '\n' < "$INPUT_FILE" | sed '/^$/d' > tmp_list.txt

while IFS= read -r species <&3; do
    species=$(echo "$species" | xargs)
    if [[ -z "$species" ]]; then continue; fi

    count=$(esearch -db sra -query "\"$species\"[Organism] AND \"illumina\"[Platform] AND \"rna seq\"[Strategy]" < /dev/null | xtract -pattern ENTREZ_DIRECT -element Count)

    if [[ -z "$count" ]]; then count=0; fi

    echo "$species,$count" >> "$OUTPUT_FILE"
    echo "Queried: $species -> $count RNA-seq record(s)"

    sleep 1
done 3< "tmp_list.txt"

rm -f tmp_list.txt

echo "Done. Results saved to $OUTPUT_FILE"
