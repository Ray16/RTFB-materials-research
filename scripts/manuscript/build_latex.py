#!/usr/bin/env python
"""Convert manuscript_draft_v0.md into a JCIM (achemso) LaTeX submission and compile it.

    manuscript/jcim/manuscript.tex   achemso source (xelatex; DejaVu fonts for Unicode)
    manuscript/jcim/references.bib   BibTeX built from the reference list (DOIs as recorded)
    manuscript/jcim/manuscript.pdf   compiled with xelatex + bibtex (achemso.bst numbering)

Only the Markdown used in this manuscript is handled: headings, paragraphs, bullet lists,
pipe tables, **bold**, *italic*, `code`, and author-year citations [Key; Key]. Figures are
placed after the paragraph that first mentions them; their captions come from
manuscript/figures/*.caption.md (provenance and file notes are left out). Inline
"**Figure N.**" paragraphs are dropped when a caption file exists.

  python scripts/manuscript/build_latex.py [--no-compile]
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "manuscript/manuscript_draft_v0.md"
FIG = ROOT / "manuscript/figures"
OUT = ROOT / "manuscript/jcim"

FIGURES = {   # number -> (caption file, [image files])
    1: ("fig1_grafting", ["fig1_grafting.png"]),
    2: ("fig2_pipeline", ["fig2_pipeline.png"]),
    3: ("fig3_qc_funnel", ["fig3_qc_funnel.png"]),
    4: ("fig4_validation_parity", ["fig4_validation_parity.png", "fig4b_residuals.png"]),
    5: ("fig5_pareto", ["fig5a_anolyte_pareto.png", "fig5b_catholyte_pareto.png"]),
    6: ("fig6_dmpz_design_rule", ["fig6_dmpz_design_rule.png", "fig6b_dmpz_catholyte_space.png"]),
    7: ("fig7_cell_voltage", ["fig7_cell_voltage.png"]),
}
STOP_CAPTION = ("**Provenance", "**Data provenance", "**Files", "**File", "**Method note",
                "**Interpretation", "**Honesty", "Source:", "Generator:")


def key_id(k: str) -> str:
    return re.sub(r"[^a-z0-9]", "", k.lower().replace("è", "e").replace("é", "e"))


# ------------------------------------------------------------------ inline conversion
def esc(s: str) -> str:
    for a, b in (("\\", r"\textbackslash{}"), ("&", r"\&"), ("%", r"\%"), ("$", r"\$"),
                 ("#", r"\#"), ("_", r"\_"), ("{", r"\{"), ("}", r"\}"), ("~", r"\textasciitilde{}"),
                 ("^", r"\textasciicircum{}")):
        s = s.replace(a, b)
    return s


def inline(s: str) -> str:
    codes = []

    def keep_code(m):
        codes.append(r"\texttt{" + esc(m.group(1)) + "}")
        return f"\x00{len(codes) - 1}\x00"
    s = re.sub(r"`([^`]+)`", keep_code, s)
    cites = []

    def keep_cite(m):
        keys = [k.strip() for k in m.group(1).split(";")]
        cites.append(r"\cite{" + ",".join(key_id(k) for k in keys) + "}")
        return f"\x01{len(cites) - 1}\x01"
    item = r"(?:[A-Z][^\[\];]*?\d{4}|RDKit)"
    s = re.sub(r"\[(" + item + r"(?:;\s*" + item + r")*)\]", keep_cite, s)
    s = esc(s)
    s = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", s)
    s = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", r"\\emph{\1}", s)
    s = re.sub(r"\x00(\d+)\x00", lambda m: codes[int(m.group(1))], s)
    s = re.sub(r"\x01(\d+)\x01", lambda m: cites[int(m.group(1))], s)
    return s


def table(lines, caption=None, label=None):
    rows = [[c.strip() for c in l.strip().strip("|").split("|")] for l in lines if not re.match(r"^\|\s*-", l)]
    ncol = len(rows[0])
    num = re.compile(r"^[−+\-]?[\d.]+(\s*\([^)]*\))?$|^—$")
    cols = []
    for k in range(ncol):
        body_cells = [r[k] for r in rows[1:] if k < len(r)]
        numeric = body_cells and all(num.match(c) for c in body_cells)
        if numeric and len(rows[0][k]) <= 14:
            cols.append("r")
        elif numeric:                      # long header: let it wrap, keep numbers right-aligned
            cols.append(r">{\raggedleft\arraybackslash}X")
        else:
            cols.append(r">{\raggedright\arraybackslash}X")
    spec = "@{}" + "".join(cols) + "@{}"
    def cell(c):   # long identifiers and chemical names may break after _ - ,
        s = inline(c).replace(r"\_", r"\_\allowbreak{}")
        return re.sub(r"(?<=[A-Za-z0-9)\]])([-,])(?=[A-Za-z0-9(])", r"\1\\allowbreak{}", s)
    body = [" & ".join(cell(c) for c in rows[0]) + r" \\ \hline"]
    body += [" & ".join(cell(c) for c in r) + r" \\" for r in rows[1:]]
    cap = (r"\caption{" + caption + "}" + (r"\label{" + label + "}" if label else "") + "\n") if caption else ""
    size = r"\footnotesize" if ncol >= 5 else r"\small"
    return ("\\begin{table}[htbp]\n" + cap + size + "\n\\begin{tabularx}{\\textwidth}{" + spec + "}\n\\hline\n"
            + "\n".join(body) + "\n\\hline\n\\end{tabularx}\n\\end{table}\n")


def caption_text(num):
    name, _ = FIGURES[num]
    # Paragraphs of the caption file; the journal caption is the paragraph that opens with
    # "**Figure N." plus any panel paragraphs "**(a)...", "**(b)...". Discussion, caveats and
    # provenance stay in the caption files (and the SI), not in the journal caption.
    paras = [p.strip() for p in re.split(r"\n\s*\n", (FIG / f"{name}.caption.md").read_text())]
    keep, started = [], False
    for p in paras:
        if p.startswith("# ") or any(p.startswith(s) for s in STOP_CAPTION):
            if started and not p.startswith("# "):
                break
            continue
        if p.startswith(f"**Figure {num}"):
            started = True
            p = re.sub(rf"^\*\*Figure {num}\.\*\*\s*", "", p)
            p = re.sub(rf"^\*\*Figure {num}\.\s*", "**", p)
            keep.append(p)
        elif started and re.match(r"^\*\*\([a-z]\)", p):
            keep.append(p)
    text = " ".join(" ".join(l.strip() for l in k.splitlines()) for k in keep)
    text = re.sub(r"\*\*\s*\*\*", " ", text).strip()
    return inline(text)


def figure_env(num):
    _, imgs = FIGURES[num]
    w = 0.98 if len(imgs) == 1 else 0.49
    # Height cap keeps image + (long) caption on one page; keepaspectratio preserves shape.
    inc = "\n".join(rf"\includegraphics[width={w}\textwidth,height=0.42\textheight,keepaspectratio]{{figures/{i}}}"
                     for i in imgs)
    return ("\\begin{figure}[htbp]\n\\centering\n" + inc + "\n\\caption{" + caption_text(num)
            + "}\n\\label{fig:" + str(num) + "}\n\\end{figure}\n")


# ------------------------------------------------------------------ references -> bib
def bib_from_refs(ref_block: str) -> str:
    out = []
    for line in ref_block.splitlines():
        m = re.match(r"^- \*\*\[([^\]]+)\]\*\*\s+(.*)$", line)
        if not m:
            continue
        key, rest = m.groups()
        doi = re.search(r"DOI:\s*(\S+?)\.?$", rest)
        if key == "RDKit":
            out.append("@misc{rdkit,\n  title = {{RDKit}: Open-Source Cheminformatics Software, version 2026.03.5},\n"
                       "  howpublished = {\\url{https://www.rdkit.org}},\n  year = {2026},\n  note = {accessed 2026-11-18}\n}\n")
            continue
        # authors; title. *Journal* **year**, *vol* (issue), pages. DOI: x.
        mm = re.match(r"(.+?(?:\.|et al\.))\s+(.+?)\.\s+\*(.+?)\*\s+\*\*(\d{4})\*\*,\s*(.*?)\s*DOI:", rest)
        if not mm:
            raise ValueError(f"cannot parse reference {key}: {rest[:80]}")
        authors, title, journal, year, tail = mm.groups()
        names = [a.strip() for a in authors.rstrip(".").split(";") if a.strip() and a.strip() != "et al"]
        etal = "et al." in authors
        vol = re.search(r"\*([^*]+)\*", tail)
        issue = re.search(r"\(([^)]+)\)", tail)
        pages = re.sub(r".*\)\s*,\s*|^\*[^*]+\*\s*,\s*", "", tail).strip().rstrip(".").strip(",").strip()
        fields = dict(author=" and ".join(names) + (" and others" if etal else ""), title="{" + title + "}",
                      journal=journal, year=year)
        if journal == "arXiv":
            fields = dict(author=fields["author"], title=fields["title"], year=year,
                          journal="arXiv", note="arXiv:2506.23971")
        else:
            if vol:
                fields["volume"] = vol.group(1)
            if issue:
                fields["number"] = issue.group(1)
            if pages and not pages.startswith("*"):
                fields["pages"] = pages.replace("–", "--")
        if doi:
            fields["doi"] = doi.group(1)
        out.append("@article{" + key_id(key) + ",\n" + ",\n".join(f"  {k} = {{{v}}}" for k, v in fields.items()) + "\n}\n")
    return "\n".join(out)


# ------------------------------------------------------------------ document
def convert():
    md = SRC.read_text()
    body, refs = md.split("\n## References", 1)
    lines = body.splitlines()
    title = lines[0].lstrip("# ").strip()
    out, i = [], 1
    abstract, in_abstract = [], False
    placed, pending_caption = set(), None
    para = []

    def flush_para():
        nonlocal para
        if not para:
            return
        text = " ".join(x.strip() for x in para)
        para = []
        m = re.match(r"^\*\*Figure (\d+)\.\*\*", text)
        if m:                                         # caption files carry the figure captions;
            n = int(m.group(1))                       # the figure itself goes here
            if n in FIGURES and n not in placed:
                out.append(figure_env(n)); placed.add(n)
            return
        if in_abstract:
            abstract.append(inline(text))
            return
        out.append(inline(text) + "\n")
        for n in sorted(set(int(x) for x in re.findall(r"Figure\s+(\d+)", text))):
            if n in FIGURES and n not in placed:
                out.append(figure_env(n)); placed.add(n)

    while i < len(lines):
        l = lines[i]
        s = l.strip()
        if s.startswith("_Draft") or s == "---":
            flush_para(); i += 1; continue
        if s.startswith("## ") or s.startswith("### "):
            flush_para()
            level, head = (2, s[3:]) if s.startswith("## ") else (3, s[4:])
            in_abstract = head == "Abstract"
            if head in ("Abstract",):
                i += 1; continue
            if head == "TOC Graphic":
                while i + 1 < len(lines) and not lines[i + 1].startswith("---"):
                    i += 1
                i += 1; continue
            m = re.match(r"Table (\d+)\.\s*(.*)", head)
            if m:
                pending_caption = (m.group(1), inline(m.group(2))); i += 1; continue
            if head == "Associated Content":
                out.append("\\section*{Associated Content}\n"); i += 1; continue
            if head == "Supporting Information":
                out.append("\\begin{suppinfo}\n"); i += 1
                block = []
                while i < len(lines) and not lines[i].startswith("### "):
                    block.append(lines[i]); i += 1
                paras, items = [], []
                for x in block:
                    s2 = x.strip()
                    if s2.startswith("- "):
                        items.append(s2[2:])
                    elif items and x.startswith("  ") and s2:
                        items[-1] += " " + s2
                    elif s2:
                        if items:
                            paras.append("\\begin{itemize}\n" + "\n".join(r"\item " + inline(i) for i in items) + "\n\\end{itemize}")
                            items = []
                        paras.append(inline(s2))
                if items:
                    paras.append("\\begin{itemize}\n" + "\n".join(r"\item " + inline(i) for i in items) + "\n\\end{itemize}")
                out.append("\n".join(paras) + "\n\\end{suppinfo}\n")
                continue
            cmd = "\\section" if level == 2 else ("\\subsection*" if head in ("Code availability", "Data availability") else "\\subsection")
            out.append(f"{cmd}{{{inline(head)}}}\n")
            i += 1; continue
        m = re.match(r"^\*\*Table (\d+)\.\*\*\s*(.*)", s)
        if m:
            flush_para()
            cap = [m.group(2)]
            i += 1
            while i < len(lines) and lines[i].strip() and not lines[i].startswith("|"):
                cap.append(lines[i].strip()); i += 1
            pending_caption = (m.group(1), inline(" ".join(cap)))
            continue
        if s.startswith("|"):
            flush_para()
            tl = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                tl.append(lines[i]); i += 1
            cap, lab = (pending_caption[1], f"tab:{pending_caption[0]}") if pending_caption else (None, None)
            out.append(table(tl, cap, lab)); pending_caption = None
            continue
        if s.startswith("- "):
            flush_para()
            items = []
            while i < len(lines) and (lines[i].strip().startswith("- ") or (lines[i].startswith("  ") and lines[i].strip())):
                if lines[i].strip().startswith("- "):
                    items.append(lines[i].strip()[2:])
                else:
                    items[-1] += " " + lines[i].strip()
                i += 1
            out.append("\\begin{itemize}\n" + "\n".join(r"\item " + inline(x) for x in items) + "\n\\end{itemize}\n")
            continue
        if not s:
            flush_para(); i += 1; continue
        para.append(l); i += 1
    flush_para()
    for n in FIGURES:                     # any figure never mentioned goes at the end
        if n not in placed:
            out.append(figure_env(n))
    doc = r"""\documentclass[journal=jcisd8,manuscript=article]{achemso}
\usepackage{fontspec}
\setmainfont{DejaVu Serif}
\setsansfont{DejaVu Sans}
\setmonofont{DejaVu Sans Mono}[Scale=0.85]
\usepackage{graphicx}
\usepackage{tabularx}
\usepackage{caption}
\captionsetup{singlelinecheck=false}
\usepackage{url}
\usepackage[hyphens]{xurl}
\sloppy
\author{[TODO: authors]}
\affiliation{[TODO: affiliation]}
\email{[TODO: corresponding author email]}
\title{""" + inline(title) + r"""}
\begin{document}
\begin{tocentry}
\includegraphics[width=3.25in,height=1.75in,keepaspectratio]{figures/toc_graphic.png}
\end{tocentry}
\begin{abstract}
""" + " ".join(abstract) + r"""
\end{abstract}
""" + "\n".join(out) + r"""
\bibliography{references}
\end{document}
"""
    return doc, bib_from_refs(refs)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "figures").mkdir(exist_ok=True)
    for _, imgs in FIGURES.values():
        for im in imgs:
            shutil.copy(FIG / im, OUT / "figures" / im)
    shutil.copy(FIG / "toc_graphic.png", OUT / "figures/toc_graphic.png")
    doc, bib = convert()
    (OUT / "manuscript.tex").write_text(doc)
    (OUT / "references.bib").write_text(bib)
    print("wrote", OUT / "manuscript.tex", "and", bib.count("@"), "bib entries")
    if "--no-compile" in sys.argv:
        return 0
    log = []
    for cmd in (["xelatex", "-interaction=nonstopmode", "manuscript.tex"], ["bibtex", "manuscript"],
                ["xelatex", "-interaction=nonstopmode", "manuscript.tex"],
                ["xelatex", "-interaction=nonstopmode", "manuscript.tex"]):
        r = subprocess.run(cmd, cwd=OUT, capture_output=True, text=True)
        log.append(r.stdout[-3000:])
    (OUT / "build.log").write_text("\n".join(log))
    print("compiled:", (OUT / "manuscript.pdf").exists())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
