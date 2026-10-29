#!/usr/bin/env python3
"""Check the per-strain FAA / GTF / GFF3 triplets before packaging.

Status: new_validation (read-only; complements the one-direction checks inside
make_release_202.sh).

Per label it checks
  * all three files exist (plain or .gz);
  * FASTA: no empty / duplicated IDs, no empty sequences, CRLF line endings;
  * GTF: transcript_id set; GFF3: ID of mRNA/transcript features;
  * ID sets FAA<->GTF and FAA<->GFF3 in BOTH directions;
  * CRLF in GTF/GFF3.
  The "##gff-version 3" header is only reported (gff3_header_ok), not
  required: the documented historical command (gtf2gff.pl --gff3) does not
  necessarily write it.
Exit codes: 0 all OK; 2 input/label error; 5 at least one strain failed.

Layout options: either a release folder (<label>/<label>_longest_isoform.*) via
--release-dir, or three separate folders via --faa-dir/--gtf-dir/--gff3-dir.
"""
import argparse, gzip, os, re, sys
from collections import Counter

TX = re.compile(r'transcript_id "([^"]+)"')
GID = re.compile(r'(?:^|;)ID=([^;]+)')


def opener(p):
    return gzip.open(p, "rb") if p.endswith(".gz") else open(p, "rb")


def find(folder, label, suffix):
    for p in (os.path.join(folder, label + suffix), os.path.join(folder, label + suffix + ".gz")):
        if os.path.isfile(p):
            return p
    return None


def read_faa(p):
    ids, crlf, empty_seq, cur, seqlen = [], False, 0, None, 0
    with opener(p) as fh:
        for raw in fh:
            if raw.endswith(b"\r\n"):
                crlf = True
            line = raw.decode("utf-8", "replace").rstrip("\r\n")
            if line.startswith(">"):
                if cur is not None and seqlen == 0:
                    empty_seq += 1
                cur = line[1:].split()[0] if line[1:].split() else ""
                ids.append(cur)
                seqlen = 0
            elif line.strip():
                seqlen += len(line.strip())
    if cur is not None and seqlen == 0:
        empty_seq += 1
    return ids, crlf, empty_seq


def read_ids(p, kind):
    ids, crlf, header = set(), False, None
    with opener(p) as fh:
        for raw in fh:
            if raw.endswith(b"\r\n"):
                crlf = True
            line = raw.decode("utf-8", "replace").rstrip("\r\n")
            if header is None:
                header = line
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 9:
                continue
            if kind == "gtf":
                m = TX.search(f[8])
                if m:
                    ids.add(m.group(1))
            elif f[2] in ("mRNA", "transcript"):
                m = GID.search(f[8])
                if m:
                    ids.add(m.group(1))
    return ids, crlf, header or ""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--labels", required=True)
    ap.add_argument("--release-dir")
    ap.add_argument("--faa-dir"); ap.add_argument("--gtf-dir"); ap.add_argument("--gff3-dir")
    ap.add_argument("--faa-suffix", default="_longest_isoform.faa")
    ap.add_argument("--gtf-suffix", default="_longest_isoform.gtf")
    ap.add_argument("--gff3-suffix", default="_longest_isoform.gff3")
    ap.add_argument("--expect-total", type=int, default=None, help="expected protein total, e.g. 2618164")
    ap.add_argument("--out", required=True, help="report TSV (must not exist)")
    a = ap.parse_args()

    if os.path.exists(a.out):
        sys.exit(f"ERROR: report exists: {a.out}")
    with open(a.labels, encoding="utf-8") as fh:
        labels = [l.strip() for l in fh if l.strip() and not l.startswith("#")]
    dup = [k for k, v in Counter(labels).items() if v > 1]
    if dup:
        print(f"ERROR: duplicated labels {dup}", file=sys.stderr); sys.exit(2)
    if not a.release_dir and not (a.faa_dir and a.gtf_dir and a.gff3_dir):
        print("ERROR: give --release-dir or all of --faa-dir --gtf-dir --gff3-dir", file=sys.stderr); sys.exit(2)

    cols = ["label", "n_faa", "n_gtf_tx", "n_gff3_mrna", "dup_faa_ids", "empty_seqs", "faa_not_in_gtf",
            "gtf_not_in_faa", "faa_not_in_gff3", "gff3_not_in_faa", "crlf", "gff3_header_ok", "missing", "status"]
    n_fail, total = 0, 0
    with open(a.out, "w", encoding="utf-8") as out:
        out.write("\t".join(cols) + "\n")
        for lab in labels:
            dirs = [os.path.join(a.release_dir, lab)] * 3 if a.release_dir else [a.faa_dir, a.gtf_dir, a.gff3_dir]
            paths = [find(d, lab, s) for d, s in zip(dirs, (a.faa_suffix, a.gtf_suffix, a.gff3_suffix))]
            missing = [k for k, p in zip(("faa", "gtf", "gff3"), paths) if p is None]
            r = dict.fromkeys(cols, "-"); r["label"] = lab; r["missing"] = ",".join(missing) or "none"
            ok = not missing
            if not missing:
                ids, c1, empty = read_faa(paths[0])
                g, c2, _ = read_ids(paths[1], "gtf")
                f3, c3, hdr = read_ids(paths[2], "gff3")
                s = set(ids)
                dupn = len(ids) - len(s) + sum(1 for i in ids if i == "")
                r.update(n_faa=len(ids), n_gtf_tx=len(g), n_gff3_mrna=len(f3), dup_faa_ids=dupn, empty_seqs=empty,
                         faa_not_in_gtf=len(s - g), gtf_not_in_faa=len(g - s), faa_not_in_gff3=len(s - f3),
                         gff3_not_in_faa=len(f3 - s), crlf=",".join(k for k, c in zip(("faa", "gtf", "gff3"), (c1, c2, c3)) if c) or "none",
                         gff3_header_ok=hdr.startswith("##gff-version 3"))
                total += len(ids)
                ok = (dupn == 0 and empty == 0 and not (s ^ g) and not (s ^ f3)
                      and r["crlf"] == "none")
            r["status"] = "OK" if ok else "FAIL"
            n_fail += not ok
            out.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"labels={len(labels)} failed={n_fail} protein_total={total}")
    if a.expect_total is not None and total != a.expect_total:
        print(f"ERROR: protein total {total} != expected {a.expect_total}", file=sys.stderr)
        n_fail += 1
    sys.exit(5 if n_fail else 0)


if __name__ == "__main__":
    main()
