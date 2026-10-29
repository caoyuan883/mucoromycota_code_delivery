#!/usr/bin/env python3
"""Generate SYNTHETIC test inputs (all labels start with SYN_). These files are
invented for software tests only; they are not biological data and their
numbers must never be reported as results.

Usage: python3 tests/synthetic_data.py OUTDIR
"""
import gzip, os, sys

H = "# SYNTHETIC TEST DATA - not biological data\n"


def w(path, text, mode="w"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, mode, newline="") as fh:
        fh.write(text)


def gtf_line(chrom, feat, s, e, strand, attrs):
    return f"{chrom}\tSYN\t{feat}\t{s}\t{e}\t.\t{strand}\t.\t{attrs}\n"


def main(o):
    # ---------- 05 longest isoform selection --------------------------------
    # g1: t1 CDS 300 span 300; t2 CDS 300 span 400 -> t2 (span tie-break)
    # g2: t1 and t2 identical CDS 150 / span 150 -> t1 (first appearance)
    # g3: t1 CDS 150; t2 CDS 300 -> t2 (CDS length)
    g = ""
    def tx(gene, t, segs):
        nonlocal g
        s, e = min(a for a, _ in segs), max(b for _, b in segs)
        g += gtf_line("c1", "transcript", s, e, "+", f"{gene}.{t}")
        for a, b in segs:
            g += gtf_line("c1", "CDS", a, b, "+", f'transcript_id "{gene}.{t}"; gene_id "{gene}";')
    g += gtf_line("c1", "gene", 1, 400, "+", "g1")
    tx("g1", "t1", [(1, 300)]); tx("g1", "t2", [(1, 150), (251, 400)])
    g += gtf_line("c1", "gene", 1001, 1150, "+", "g2")
    tx("g2", "t1", [(1001, 1150)]); tx("g2", "t2", [(1001, 1150)])
    g += gtf_line("c1", "gene", 2001, 2300, "+", "g3")
    tx("g3", "t1", [(2001, 2150)]); tx("g3", "t2", [(2001, 2300)])
    w(f"{o}/longest/gtf/SYN_A.gtf", g)
    faa = "".join(f">{i}\nM{'A' * 20}\n" for i in ["g1.t1", "g1.t2", "g2.t1", "g2.t2", "g3.t1", "g3.t2"])
    w(f"{o}/longest/protein/SYN_A_protein.faa", faa)

    # ---------- 10 release triplet checks -----------------------------------
    ids = ["g1.t1", "g2.t1", "g3.t1"]
    fa = "".join(f">{i}\nMAAAA\n" for i in ids)
    # BRAKER-like GTF: bare gene / transcript lines + CDS lines
    gt = "".join(gtf_line("c1", "gene", 1, 9, "+", i.split(".")[0]) + gtf_line("c1", "transcript", 1, 9, "+", i)
                 + gtf_line("c1", "CDS", 1, 9, "+", f'transcript_id "{i}"; gene_id "{i.split(".")[0]}";') for i in ids)
    gf = "##gff-version 3\n" + "".join(
        f"c1\tSYN\tgene\t1\t9\t.\t+\t.\tID={i.split('.')[0]}\n"
        f"c1\tSYN\tmRNA\t1\t9\t.\t+\t.\tID={i};Parent={i.split('.')[0]}\n"
        f"c1\tSYN\tCDS\t1\t9\t.\t+\t0\tID={i}.cds;Parent={i}\n" for i in ids)
    def trip(case, fa_, gt_, gf_, label="SYN_A", gz=False, skip=()):
        for kind, txt, suf in (("faa", fa_, "_longest_isoform.faa"), ("gtf", gt_, "_longest_isoform.gtf"), ("gff3", gf_, "_longest_isoform.gff3")):
            if kind in skip:
                continue
            p = f"{o}/release_check/{case}/{kind}/{label}{suf}"
            os.makedirs(os.path.dirname(p), exist_ok=True)
            if gz:
                with gzip.open(p + ".gz", "wt") as fh:
                    fh.write(txt)
            else:
                w(p, txt)
    trip("ok", fa, gt, gf)
    trip("ok_gz", fa, gt, gf, gz=True)
    trip("crlf", fa.replace("\n", "\r\n"), gt, gf)
    trip("dup_id", fa + ">g1.t1\nMAAAA\n", gt, gf)
    trip("gff3_extra", fa, gt, gf + "c1\tSYN\tmRNA\t1\t9\t.\t+\t.\tID=g9.t1;Parent=g9\n")
    trip("missing", fa, gt, gf, skip=("gff3",))
    w(f"{o}/release_check/labels.txt", "SYN_A\n")

    # ---------- 10 subphylum map --------------------------------------------
    hdr = "Record No.\tSource identifier\tScientific name (Workbook 1)\tName used in this study\tStrain/isolate\tSubphylum\tSource database\n"
    w(f"{o}/subphylum/wb3_ok.tsv", hdr +
      "1\tSYN1\tSynthomyces alpha\tSynthomyces alpha\tX1\tMucoromycotina\tSYN\n"
      "2\tSYN2\tSynthomyces sp. AB-1\tSynthomyces sp. AB-1\tAB-1\tMortierellomycotina\tSYN\n"
      "3\tSYN3\tSynthomyces sp. KO_9\tSynthomyces sp. KO_9\tKO 9\tGlomeromycotina\tSYN\n")
    w(f"{o}/subphylum/labels_ok.txt", "Synthomyces_alpha\nSynthomyces_sp.AB-1\nSynthomyces_sp._KO9\n")
    # record 2 names point to two different existing labels -> ambiguous
    w(f"{o}/subphylum/wb3_ambiguous.tsv", hdr +
      "1\tSYN1\tSynthomyces alpha\tSynthomyces alpha\tX1\tMucoromycotina\tSYN\n"
      "2\tSYN2\tSynthomyces beta\tSynthomyces alpha\tX2\tMucoromycotina\tSYN\n")
    w(f"{o}/subphylum/labels_ambiguous.txt", "Synthomyces_alpha\nSynthomyces_beta\n")
    w(f"{o}/subphylum/labels_dup.txt", "Synthomyces_alpha\nSynthomyces_alpha\nSynthomyces_sp.AB-1\n")

    # ---------- 05 prefilter view -------------------------------------------
    for lab in ("SYN_A", "SYN_B", "SYN_C"):
        w(f"{o}/prefilter/final/{lab}_longest_isoform.faa", f">x\nMFINAL{lab}\n")
    w(f"{o}/prefilter/pre/SYN_B_longest_isoform.faa", ">x\nMPREFILTER\n")
    w(f"{o}/prefilter/labels.txt", "SYN_A\nSYN_B\nSYN_C\n")
    w(f"{o}/prefilter/replace.txt", "SYN_B\n")
    w(f"{o}/prefilter/replace_missing.txt", "SYN_C\n")

    # ---------- 07 func_filtered10 (denominator, GO split, union/intersection)
    # kept after filtering: p1 p2 p3 p4 ; removed: p5
    ips = ("p1\tmd5\t100\tPfam\tPF1\td\t1\t50\t1e-5\tT\t01-01-2026\tIPR000001\tdesc\tGO:0000001(InterPro)\t-\n"
           "p2\tmd5\t100\tPANTHER\tPTHR1\td\t1\t50\t1e-5\tT\t01-01-2026\t-\t-\tGO:0000002(PANTHER)\t-\n"
           "p5\tmd5\t100\tPfam\tPF2\td\t1\t50\t1e-5\tT\t01-01-2026\tIPR000002\tdesc\tGO:0000003(InterPro)\t-\n")
    p = f"{o}/func/ips/SYN_A.interproscan.tsv.gz"; os.makedirs(os.path.dirname(p), exist_ok=True)
    with gzip.open(p, "wt") as fh:
        fh.write(ips)
    cols = ["query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl", "COG_category", "Description",
            "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass",
            "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs"]
    def em(q, cog="-", go="-", ko="-"):
        v = dict.fromkeys(cols, "-"); v.update(query=q, COG_category=cog, GOs=go, KEGG_ko=ko)
        return "\t".join(v[c] for c in cols) + "\n"
    w(f"{o}/func/em/SYN_A.emapper.annotations",
      "## SYNTHETIC\n#" + "\t".join(cols) + "\n" + em("p2", "S", "GO:0000002", "ko:K00001") + em("p3", "S") + em("p5", "S", "GO:0000009"))
    w(f"{o}/func/filtered/SYN_A_longest_isoform.filtered.faa", "".join(f">{i}\nMAAAA\n" for i in ("p1", "p2", "p3", "p4")))
    w(f"{o}/func/orig/SYN_A_longest_isoform.faa", "".join(f">{i}\nMAAAA\n" for i in ("p1", "p2", "p3", "p4", "p5")))

    # pre-filter base table (header = func_summary_v3 output columns); values are synthetic
    hdr = ['strain', 'n_proteins', 'ips_hit', 'ips_pct', 'ips_ipr', 'ips_ipr_pct', 'ips_Pfam', 'ips_PANTHER', 'ips_SUPERFAMILY', 'ips_CDD', 'ips_SMART', 'ips_ProSitePatterns', 'ips_ProSiteProfiles', 'ips_AntiFam', 'ips_go', 'ips_go_pct', 'ips_go_InterPro', 'ips_go_PANTHER', 'egg_hit', 'egg_pct', 'egg_go', 'egg_go_pct', 'egg_ko', 'egg_ko_pct', 'egg_ec', 'egg_cazy', 'egg_cog', 'egg_kegg_pathway', 'egg_preferred_name', 'union_hit', 'union_pct', 'go_union', 'go_union_pct', 'both_tools', 'ips_only', 'egg_only']
    row = lambda lab, n: "\t".join([lab, str(n)] + ["0"] * (len(hdr) - 2)) + "\n"
    w(f"{o}/func/01_per_strain_functional_summary.tsv", "\t".join(hdr) + "\n" + row("SYN_A", 5) + row("SYN_B", 7))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
