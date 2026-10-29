#!/usr/bin/env python3
import csv
import re
from collections import Counter

input_file = "/Users/$USER/Desktop/JGI189.tsv"
output_file = "/Users/$USER/Desktop/JGI189_species_count.tsv"

def normalize_name(name):
    s = str(name).strip()
    s = re.sub(r'\s+', ' ', s)
    if not s:
        return None

    parts = s.split(' ')
    if len(parts) < 2:
        return s

    genus = parts[0]
    second = parts[1]

    if second == "sp.":
        if len(parts) >= 3:
            return f"{genus} sp. {parts[2]}"
        else:
            return f"{genus} sp."
    else:
        return f"{genus} {second}"

counts = Counter()
raw_total = 0

with open(f"/Users/{__import__('os').environ['USER']}/Desktop/JGI189.tsv", "r", encoding="utf-8-sig", newline="") as f:
    sample = f.read(4096)
    f.seek(0)
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
        delim = dialect.delimiter
    except:
        delim = "\t"

    reader = csv.reader(f, delimiter=delim)
    for row in reader:
        if not row:
            continue
        first_col = row[0].strip()
        if not first_col:
            continue
        raw_total += 1
        species = normalize_name(first_col)
        if species:
            counts[species] += 1

with open(f"/Users/{__import__('os').environ['USER']}/Desktop/JGI189_species_count.tsv", "w", encoding="utf-8", newline="") as out:
    out.write("species_name\tcount\n")
    for species_name, count in sorted(counts.items(), key=lambda x: (-x[1], x[0])):
        out.write(f"{species_name}\t{count}\n")

print("Done")
print("Input:", f"/Users/{__import__('os').environ['USER']}/Desktop/JGI189.tsv")
print("Output:", f"/Users/{__import__('os').environ['USER']}/Desktop/JGI189_species_count.tsv")
print("Total rows counted:", raw_total)
print("Unique species:", len(counts))
