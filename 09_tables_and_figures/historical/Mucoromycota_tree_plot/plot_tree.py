#!/usr/bin/env python3
import html
import re
from pathlib import Path

TREE_FILE = "Mucoromycota_212_300BUSCO_JTTFG4.treefile"
OUTPUT_FILE = "Mucoromycota_212_ML_raw_phylogram.svg"

class Node:
    def __init__(self, children=None, name="", length=0.0):
        self.children = children or []
        self.name = name
        self.length = length
        self.x = 0.0
        self.y = 0.0

def parse_newick(text):
    i = 0
    def token_until(stops):
        nonlocal i
        start = i
        while i < len(text) and text[i] not in stops:
            i += 1
        return text[start:i].strip()
    def parse_node():
        nonlocal i
        if text[i] == "(":
            i += 1
            children = []
            while True:
                children.append(parse_node())
                if text[i] == ",":
                    i += 1
                    continue
                if text[i] == ")":
                    i += 1
                    break
            name = token_until(":,();")
            length = 0.0
            if i < len(text) and text[i] == ":":
                i += 1
                value = token_until(",();")
                length = float(value) if value else 0.0
            return Node(children, name, length)
        name = token_until(":,();")
        length = 0.0
        if i < len(text) and text[i] == ":":
            i += 1
            value = token_until(",();")
            length = float(value) if value else 0.0
        return Node([], name, length)
    root = parse_node()
    if i < len(text) and text[i] == ";":
        i += 1
    if text[i:].strip():
        raise ValueError("Newick文件结尾存在无法识别的内容")
    return root

def get_support(label):
    match = re.fullmatch(r"([0-9.]+)/([0-9.]+)", label or "")
    if not match:
        return None
    return float(match.group(1)), float(match.group(2))

def render(root):
    tips = []
    nodes = []
    def walk(node, distance=0.0):
        node.x = distance
        nodes.append(node)
        if node.children:
            for child in node.children:
                walk(child, distance + child.length)
        else:
            tips.append(node)
    walk(root)
    for index, tip in enumerate(tips):
        tip.y = index
    def set_y(node):
        if node.children:
            for child in node.children:
                set_y(child)
            node.y = sum(child.y for child in node.children) / len(node.children)
    set_y(root)

    max_distance = max(node.x for node in nodes)
    width = 2500
    top = 95
    bottom = 110
    row_height = 20
    height = top + bottom + row_height * (len(tips) - 1)
    tree_left = 90
    tree_right = 1700
    scale = (tree_right - tree_left) / max_distance
    X = lambda value: tree_left + value * scale
    Y = lambda value: top + value * row_height

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,"Noto Sans",sans-serif}.tip{font-size:12px;fill:#202124;font-style:italic}.flag{font-size:10px;fill:#9b1c1c}.title{font-size:22px;font-weight:700}.sub{font-size:13px;fill:#555}.branch{stroke:#555;stroke-width:1;fill:none}</style>',
        '<text x="90" y="32" class="title">Mucoromycota 212 taxa — IQ-TREE ML raw phylogram</text>',
        '<text x="90" y="57" class="sub">JTT+F+G4; flagged labels are SH-aLRT/UFBoot. The input tree is unrooted.</text>'
    ]

    flagged = 0
    for node in nodes:
        if not node.children:
            continue
        child_y = [Y(child.y) for child in node.children]
        svg.append(f'<line class="branch" x1="{X(node.x):.2f}" y1="{min(child_y):.2f}" x2="{X(node.x):.2f}" y2="{max(child_y):.2f}"/>')
        for child in node.children:
            svg.append(f'<line class="branch" x1="{X(node.x):.2f}" y1="{Y(child.y):.2f}" x2="{X(child.x):.2f}" y2="{Y(child.y):.2f}"/>')
        values = get_support(node.name)
        if values and (values[0] < 80 or values[1] < 95):
            flagged += 1
            color = "#d93025" if values[0] < 80 and values[1] < 95 else "#f29900"
            svg.append(f'<circle cx="{X(node.x):.2f}" cy="{Y(node.y):.2f}" r="4" fill="{color}"/>')
            svg.append(f'<text x="{X(node.x)+6:.2f}" y="{Y(node.y)-5:.2f}" class="flag">{html.escape(node.name)}</text>')

    for tip in tips:
        label = tip.name.replace("_", " ")
        svg.append(f'<text x="{X(tip.x)+5:.2f}" y="{Y(tip.y)+4:.2f}" class="tip">{html.escape(label)}</text>')

    base_y = height - 42
    bar_width = 0.1 * scale
    svg.extend([
        f'<line x1="90" y1="{base_y}" x2="{90+bar_width:.2f}" y2="{base_y}" stroke="#111" stroke-width="2"/>',
        f'<line x1="90" y1="{base_y-5}" x2="90" y2="{base_y+5}" stroke="#111"/>',
        f'<line x1="{90+bar_width:.2f}" y1="{base_y-5}" x2="{90+bar_width:.2f}" y2="{base_y+5}" stroke="#111"/>',
        f'<text x="90" y="{base_y+22}" class="sub">0.1 substitutions/site</text>',
        f'<circle cx="480" cy="{base_y}" r="4" fill="#d93025"/><text x="491" y="{base_y+4}" class="sub">both SH-aLRT &lt; 80 and UFBoot &lt; 95</text>',
        f'<circle cx="800" cy="{base_y}" r="4" fill="#f29900"/><text x="811" y="{base_y+4}" class="sub">one value below threshold</text>',
        '</svg>'
    ])
    Path(OUTPUT_FILE).write_text("\n".join(svg), encoding="utf-8")
    return len(tips), flagged

def main():
    path = Path(TREE_FILE)
    if not path.exists():
        raise SystemExit(f"找不到树文件：{TREE_FILE}\n请确认树文件和plot_tree.py位于同一个文件夹。")
    root = parse_newick(path.read_text(encoding="utf-8").strip())
    tip_count, flagged_count = render(root)
    print("绘图成功！")
    print(f"物种数：{tip_count}")
    print(f"标记的低支持节点数：{flagged_count}")
    print(f"输出文件：{OUTPUT_FILE}")

if __name__ == "__main__":
    main()
