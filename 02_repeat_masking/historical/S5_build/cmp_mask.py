#!/usr/bin/env python3
"""比较两份软屏蔽 FASTA 的小写(屏蔽)位置重合度。
用法:
  cmp_mask.py A.fa B.fa            # 人读格式
  cmp_mask.py A.fa B.fa --tsv NAME # 输出一行 TSV
"""
import sys, os

T = bytes.maketrans(bytes(range(256)),
                    bytes(1 if 97 <= i <= 122 else 0 for i in range(256)))
CH = 1 << 22   # 4 MB 分块, 限制内存

def load(p):
    d, name, buf = {}, None, []
    with open(p, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None: d[name] = b"".join(buf)
                name = line[1:].split()[0]; buf = []
            else:
                buf.append(line.strip())
    if name is not None: d[name] = b"".join(buf)
    return d

def main():
    tsv = None
    args = sys.argv[1:]
    if "--tsv" in args:
        i = args.index("--tsv"); tsv = args[i + 1]; args = args[:i] + args[i + 2:]
    pa, pb = args[0], args[1]
    A, B = load(pa), load(pb)
    ka, kb = set(A), set(B)
    common = ka & kb
    la = lb = inter = tot = 0
    lendiff = 0
    for k in common:
        x, y = A[k], B[k]
        if len(x) != len(y):
            lendiff += 1; continue
        tot += len(x)
        for s in range(0, len(x), CH):
            mx = x[s:s + CH].translate(T)
            my = y[s:s + CH].translate(T)
            la += mx.count(1); lb += my.count(1)
            n = len(mx)
            ix = int.from_bytes(mx, "big") & int.from_bytes(my, "big")
            inter += ix.to_bytes(n, "big").count(1)
    union = la + lb - inter
    name = tsv if tsv else os.path.basename(pa)
    if tot == 0 or union == 0:
        msg = f"{name}\tNA\tNA\tNA\tNA\tNA\tNA\tNO_COMMON_SEQ(共有序列{len(common)},长度不符{lendiff})"
        (print(msg) if tsv else sys.stderr.write(msg + "\n"))
        return
    jac = inter / union * 100
    disc = (union - inter) / tot * 100
    if tsv:
        print(f"{name}\t{tot}\t{la}\t{lb}\t{inter}\t{jac:.4f}\t{disc:.4f}\tOK")
    else:
        print(f"A: {os.path.basename(pa)}\nB: {os.path.basename(pb)}")
        print(f"  共有序列 {len(common)} | 仅A {len(ka-kb)} | 仅B {len(kb-ka)} | 长度不符 {lendiff}")
        print(f"  序列总长 {tot:,}")
        print(f"  A 小写 {la:,} ({la/tot*100:.3f}%)   B 小写 {lb:,} ({lb/tot*100:.3f}%)")
        print(f"  交 {inter:,}  并 {union:,}")
        print(f"  >>> Jaccard = {jac:.3f}%   不一致占全基因组 = {disc:.3f}%")

main()
