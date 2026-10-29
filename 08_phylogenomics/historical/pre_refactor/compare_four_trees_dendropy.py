#!/usr/bin/env python3
"""
Four-tree topology comparison for the Technical Validation (P/G x occ80/occ90).

Input : the four IQ-TREE .treefile files (ML trees, labels "SH-aLRT/UFBoot")
Output: printed summary + two TSV files in OUTDIR

  pip install dendropy==5.1.0
  python3 compare_four_trees_dendropy.py

All numbers are computed on the 202 ingroup strains after pruning the
5 outgroups, as unrooted trees. nRF = RF / (2n - 6).
Strong support = SH-aLRT >= 80 and UFBoot >= 95.
"""

import os
import statistics
import dendropy
from dendropy.calculate import treecompare

# ---------------------------------------------------------------- settings
TREE_DIR = os.path.expanduser("~/dual_tree")          # folder with the four .treefile
OUTDIR = os.path.join(TREE_DIR, "tv_comparison")
TREES = ["P_occ80_g90", "G_occ80_g90", "P_occ90_g90", "G_occ90_g90"]
PAIRS = [("P_occ80_g90", "G_occ80_g90"),
         ("P_occ90_g90", "G_occ90_g90"),
         ("P_occ80_g90", "P_occ90_g90"),
         ("G_occ80_g90", "G_occ90_g90")]
OUTGROUPS = ["Neurospora_crassa", "Aspergillus_nidulans", "Coprinopsis_cinerea",
             "Conidiobolus_coronatus", "Coemansia_reversa"]
SH_MIN, UF_MIN = 80.0, 95.0
GLO = set("Acaulospora Ambispora Appendiculispora Dentiscutata Diversispora "
          "Entrophospora Funneliformis Geosiphon Gigaspora Glomus Oehlia "
          "Paraglomus Rhizophagus".split())
MOR = set("Actinomortierella Cryptocalciella Dissophora Entomortierella Gamsiella "
          "Gryganskiella Haplosporangium Linnemannia Lobosporangium "
          "Lunasporangiospora Modicella Mortierella Podila".split())
# ------------------------------------------------------------------------

os.makedirs(OUTDIR, exist_ok=True)
tns = dendropy.TaxonNamespace()


def load(name):
    t = dendropy.Tree.get(path=os.path.join(TREE_DIR, name + ".treefile"),
                          schema="newick", taxon_namespace=tns,
                          preserve_underscores=True)
    t.prune_taxa_with_labels(OUTGROUPS)
    t.is_rooted = False
    t.encode_bipartitions()
    return t


def support(node):
    lab = node.label or ""
    if "/" in lab:
        a, b = lab.split("/")[:2]
        return float(a), float(b)
    return None, None


def splits(tree):
    """bipartition -> (leaf set of the smaller side, SH, UF, branch length)"""
    out = {}
    allleaves = frozenset(l.taxon.label for l in tree.leaf_node_iter())
    for nd in tree.postorder_internal_node_iter():
        if nd.parent_node is None:
            continue
        side = frozenset(l.taxon.label for l in nd.leaf_iter())
        if len(side) < 2 or len(side) > len(allleaves) - 2:
            continue
        other = allleaves - side
        key = min(side, other, key=lambda s: (len(s), sorted(s)[0]))
        sh, uf = support(nd)
        out[key] = (sh, uf, nd.edge.length or 0.0)
    return out, allleaves


def strong(v):
    return v[0] is not None and v[0] >= SH_MIN and v[1] >= UF_MIN


def compatible(a, b, T):
    ca = T - a
    return (not (a & b)) or a <= b or b <= a or not (ca & b) or ca <= b or b <= ca


trees = {k: load(k) for k in TREES}
sp = {k: splits(trees[k]) for k in TREES}
T = sp[TREES[0]][1]
n = len(T)

print("== per tree (202 ingroup strains, unrooted)")
for k in TREES:
    s, _ = sp[k]
    lens = [v[2] for v in s.values()]
    subs = {}
    for lab, grp in (("Glomeromycotina", GLO), ("Mortierellomycotina", MOR)):
        S = frozenset(x for x in T if x.split("_")[0] in grp)
        key = min(S, T - S, key=lambda z: (len(z), sorted(z)[0]))
        subs[lab] = "monophyletic" if key in s else "NOT monophyletic"
    print(f"{k}: internal={len(s)} strong={sum(strong(v) for v in s.values())} "
          f"median_internal_bl={statistics.median(lens):.4f} {subs}")

# genus-level monophyly (genus = first token of the tip label), genera with >= 2 strains
genera = {}
for x in T:
    genera.setdefault(x.split("_")[0], set()).add(x)
genera = {g: frozenset(v) for g, v in genera.items() if len(v) >= 2}
status = {}
for k in TREES:
    s, _ = sp[k]
    status[k] = {g: (min(v, T - v, key=lambda z: (len(z), sorted(z)[0])) in s)
                 for g, v in genera.items()}
print(f"\n== genera with >=2 strains: {len(genera)}")
for k in TREES:
    print(f"{k}: monophyletic={sum(status[k].values())} "
          f"non-monophyletic={sorted(g for g, m in status[k].items() if not m)}")
print("identical status in all four trees:",
      all(status[k] == status[TREES[0]] for k in TREES))

rows = []
detail = []
print("\n== pairwise")
for a, b in PAIRS:
    sa, sb = sp[a][0], sp[b][0]
    rf = treecompare.symmetric_difference(trees[a], trees[b])
    shared = set(sa) & set(sb)
    oa, ob = set(sa) - set(sb), set(sb) - set(sa)
    ha = [x for x in oa if strong(sa[x])]
    hb = [y for y in ob if strong(sb[y])]
    inc = [(x, y) for x in ha for y in hb if not compatible(x, y, T)]
    regions = []
    for x, y in inc:
        r = set(x) | set(y)
        for c in [c for c in regions if c & r]:
            regions.remove(c)
            r |= c
        regions.append(r)
    maxbl = max([sa[x][2] for x in oa] + [sb[y][2] for y in ob]) if rf else 0
    rows.append((a, b, len(shared), rf, rf / (2 * n - 6), len(inc), len(regions), maxbl))
    print(f"{a} vs {b}: shared={len(shared)}/{n-3} RF={rf} nRF={rf/(2*n-6):.4f} "
          f"strong-vs-strong incompatible pairs={len(inc)} regions={len(regions)} "
          f"max_bl_of_differing_branch={maxbl:.4f}")
    for r in regions:
        print("    region:", ", ".join(sorted(r)))
    for side, d, keys in ((a, sa, oa), (b, sb, ob)):
        for x in sorted(keys, key=len):
            sh, uf, bl = d[x]
            detail.append((f"{a}|{b}", side, len(x), sh, uf, "yes" if strong(d[x]) else "no",
                           f"{bl:.5f}", ";".join(sorted(x))))

with open(os.path.join(OUTDIR, "pairwise_summary.tsv"), "w") as fh:
    fh.write("tree_A\ttree_B\tshared_of_199\tRF\tnRF\tstrong_incompatible_pairs\t"
             "strong_conflict_regions\tmax_bl_differing\n")
    for r in rows:
        fh.write("%s\t%s\t%d\t%d\t%.4f\t%d\t%d\t%.4f\n" % r)
with open(os.path.join(OUTDIR, "differing_splits.tsv"), "w") as fh:
    fh.write("comparison\tsplit_only_in\tn_taxa\tSH_aLRT\tUFBoot\tstrong\tbranch_length\tmembers\n")
    for r in detail:
        fh.write("\t".join(str(x) for x in r) + "\n")
print("\nwritten:", OUTDIR)
