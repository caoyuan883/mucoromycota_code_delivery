cd ~/Desktop

tail -n +2 12_tables_annotation_ready_genomes.csv | cut -d',' -f1 | tr -d '\r' > species_list.txt

wc -l species_list.txt
head -20 species_list.txt