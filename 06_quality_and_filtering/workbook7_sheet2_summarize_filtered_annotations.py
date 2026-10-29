#!/usr/bin/env python3
"""Read-only audit of paired filtered GTF/FAA files for Workbook 7.

Python >= 3.9; standard library only. No biological sequence filtering is performed.
Sample names and expected counts are read from the OMArk before/after TSV.
Counts and ID checks do not establish CDS-to-protein translation equivalence.
"""
import argparse
import csv
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

VERSION = "1.0.0"
TID = re.compile(r'transcript_id\s+"([^"]+)"')
GID = re.compile(r'gene_id\s+"([^"]+)"')


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def parse_gtf(path):
    genes, transcripts = set(), set()
    tx_genes = defaultdict(set)
    malformed = unmapped = 0
    with path.open(encoding='utf-8') as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                malformed += 1
                continue
            feature, attributes = fields[2], fields[8]
            tm, gm = TID.search(attributes), GID.search(attributes)
            tid = tm.group(1) if tm else None
            gid = gm.group(1) if gm else None
            bare = attributes.strip().strip(';').strip('"')
            if not tid and feature in ('transcript', 'mRNA'):
                if bare and not re.search(r'[;\s=]', bare):
                    tid = bare
            if not gid and feature == 'gene':
                if bare and not re.search(r'[;\s=]', bare):
                    gid = bare
            if gid:
                genes.add(gid)
            if tid:
                transcripts.add(tid)
                if gid:
                    tx_genes[tid].add(gid)
            if feature in ('exon', 'CDS') and not tid:
                unmapped += 1
    return genes, transcripts, tx_genes, malformed, unmapped


def read_faa(path):
    identifiers, lengths = [], []
    current = None
    with path.open(encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('>'):
                if current is not None:
                    lengths.append(current)
                header = line[1:].strip()
                identifiers.append(header.split()[0] if header else '')
                current = 0
            elif line.strip():
                if current is None:
                    raise ValueError('FASTA sequence before the first header')
                current += len(''.join(line.split()).replace('*', ''))
    if current is not None:
        lengths.append(current)
    return identifiers, sum(n == 0 for n in lengths)


def write_csv(path, rows):
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='', encoding='utf-8-sig') as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version', action='version', version=VERSION)
    parser.add_argument('--filtered-root', required=True, type=Path)
    parser.add_argument('--source-gtf-dir', required=True, type=Path)
    parser.add_argument('--omark-table', required=True, type=Path)
    parser.add_argument('--out-dir', required=True, type=Path,
                        help='New directory; an existing directory is never overwritten.')
    args = parser.parse_args()
    args.filtered_root = args.filtered_root.expanduser().resolve()
    args.source_gtf_dir = args.source_gtf_dir.expanduser().resolve()
    args.omark_table = args.omark_table.expanduser().resolve()
    args.out_dir = args.out_dir.expanduser().resolve()
    if not args.filtered_root.is_dir() or not args.source_gtf_dir.is_dir():
        parser.error('An input directory does not exist.')
    with args.omark_table.open(encoding='utf-8-sig', newline='') as handle:
        expected = list(csv.DictReader(handle, delimiter='\t'))
    if not expected or not {'sample','before_proteins','after_proteins'} <= set(expected[0]):
        parser.error('OMArk TSV must contain sample, before_proteins and after_proteins.')
    samples = [r['sample'] for r in expected]
    if len(samples) != len(set(samples)):
        parser.error('Duplicate sample labels in the OMArk table.')
    if any(not s or s in ('.','..') or '/' in s or '\\' in s for s in samples):
        parser.error('Unsafe or empty sample label.')
    args.out_dir.mkdir(parents=True, exist_ok=False)
    details = args.out_dir / 'id_differences'
    details.mkdir()
    rows, summary, hashes = [], [], []
    for expected_row in expected:
        sample = expected_row['sample']
        before = int(expected_row['before_proteins'])
        after_expected = int(expected_row['after_proteins'])
        folder = args.filtered_root / sample
        gtf = folder / (sample + '_longest_isoform.filtered.gtf')
        faa = folder / (sample + '_longest_isoform.filtered.faa')
        source = args.source_gtf_dir / (sample + '_longest_isoform.gtf')
        row = {'File label': sample}
        missing = [str(p) for p in (gtf, faa) if not p.is_file()]
        if missing:
            row.update({'Check status':'MISSING_FILE', 'Details':'; '.join(missing)})
            rows.append(row)
            print(sample, row['Check status'])
            continue
        try:
            tracked = [p for p in (gtf, faa, source) if p.is_file()]
            stamps = {str(p):(p.stat().st_size, p.stat().st_mtime_ns) for p in tracked}
            genes, tx, tx_genes, bad, unmapped = parse_gtf(gtf)
            ids, empty = read_faa(faa)
            proteins = set(ids)
            duplicate = len(ids)-len(proteins)
            source_genes = source_tx = None
            source_bad = source_unmapped = None
            if source.is_file():
                source_genes, source_tx, _, source_bad, source_unmapped = parse_gtf(source)
            only_gtf, only_faa = tx-proteins, proteins-tx
            for label, values in [('gtf_only_ids',only_gtf),('faa_only_ids',only_faa)]:
                (details / (sample+'.'+label+'.txt')).write_text(
                    '\n'.join(sorted(values))+ ('\n' if values else ''), encoding='utf-8')
            gene_tx = defaultdict(set)
            for tid, gids in tx_genes.items():
                for gid in gids:
                    gene_tx[gid].add(tid)
            one_per_gene = (bool(tx)
                and all(len(tx_genes.get(t,set()))==1 for t in tx)
                and set(gene_tx)==genes
                and all(len(v)==1 for v in gene_tx.values()))
            issues=[]
            if len(ids)!=after_expected: issues.append('FAA_COUNT_DIFFERS_FROM_OMARK')
            if len(genes)!=len(tx) or len(tx)!=len(ids): issues.append('COUNT_MISMATCH')
            if only_gtf or only_faa: issues.append('DIRECT_ID_MATCH_UNRESOLVED')
            if not one_per_gene: issues.append('GENE_TRANSCRIPT_MAPPING_CHECK')
            if bad or unmapped: issues.append('GTF_PARSE_CHECK')
            if source_bad or source_unmapped: issues.append('SOURCE_GTF_PARSE_CHECK')
            if duplicate or '' in proteins or empty: issues.append('FAA_HEADER_OR_SEQUENCE_CHECK')
            if source_genes is None: issues.append('SOURCE_GTF_MISSING')
            elif tx-source_tx: issues.append('FILTERED_TRANSCRIPTS_NOT_IN_SOURCE')
            for p in tracked:
                digest=sha256(p)
                stamp=(p.stat().st_size,p.stat().st_mtime_ns)
                if stamp!=stamps[str(p)]:
                    raise RuntimeError('Input changed during inspection: '+str(p))
                hashes.append({'sample':sample,'file':str(p),'sha256':digest,'bytes':stamp[0]})
            removed = len(source_genes)-len(genes) if source_genes is not None else None
            percent = (100*removed/len(source_genes)) if source_genes else None
            if removed is not None and removed<0: issues.append('NEGATIVE_GENE_REMOVAL')
            row.update({
                'Genes in source GTF (n)':len(source_genes) if source_genes is not None else '',
                'Genes in filtered GTF (n)':len(genes),
                'Transcripts in filtered GTF (n)':len(tx),
                'Proteins before filtering (OMArk)':before,
                'Proteins after filtering (FAA)':len(ids),
                'Proteins after filtering (OMArk)':after_expected,
                'Proteins removed (OMArk)':before-after_expected,
                'Proteins removed (%)':round(100*(before-after_expected)/before,2) if before else '',
                'Gene models removed (n)':removed if removed is not None else '',
                'Gene models removed (%)':round(percent,2) if percent is not None else '',
                'GTF-only transcript IDs (n)':len(only_gtf),
                'FAA-only protein IDs (n)':len(only_faa),
                'Duplicate FAA IDs (n)':duplicate,'Empty protein sequences (n)':empty,
                'One transcript per gene':one_per_gene,
                'Source GTF transcripts absent from filtered GTF (n)':len(source_tx-tx) if source_tx is not None else '',
                'Filtered GTF transcripts absent from source GTF (n)':len(tx-source_tx) if source_tx is not None else '',
                'Malformed GTF lines (n)':bad,
                'Exon/CDS records without transcript ID (n)':unmapped,
                'Check status':'; '.join(issues) if issues else 'PASS_COUNTS_AND_IDS',
                'Filtered GTF':str(gtf),'Filtered FAA':str(faa)})
            summary.append({k:row[k] for k in [
                'File label','Genes in source GTF (n)','Genes in filtered GTF (n)',
                'Transcripts in filtered GTF (n)','Proteins after filtering (FAA)',
                'Gene models removed (n)','Gene models removed (%)','Check status']})
        except Exception as exc:
            row.update({'Check status':'READ_OR_PARSE_ERROR','Details':str(exc)})
        rows.append(row)
        print(sample, row['Check status'])
    write_csv(args.out_dir/'filtered_gtf_faa_check.csv',rows)
    if summary: write_csv(args.out_dir/'filtering_summary_counts.csv',summary)
    if hashes: write_csv(args.out_dir/'input_checksums.csv',hashes)
    metadata={'script_version':VERSION,'script_sha256':sha256(Path(__file__)),
              'generated_utc':datetime.now(timezone.utc).isoformat(),
              'arguments':{k:str(v) for k,v in vars(args).items()},
              'omark_table_sha256':sha256(args.omark_table),'samples':len(expected),
              'status_counts':dict(Counter(r['Check status'] for r in rows))}
    (args.out_dir/'run_metadata.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2)+'\n')
    print('Output:', args.out_dir)
    return 1 if any(r['Check status']!='PASS_COUNTS_AND_IDS' for r in rows) else 0

if __name__=='__main__':
    sys.exit(main())
