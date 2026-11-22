#!/usr/bin/env python
"""Convert manuscript/supporting_information.md into an ACS/JCIM SI LaTeX/PDF.

    manuscript/jcim/supporting_information.tex   XeLaTeX source (DejaVu fonts)
    manuscript/jcim/supporting_information.pdf   compiled with xelatex (3 passes)
    manuscript/jcim/si_build.log                 compiler log

Only the Markdown used in supporting_information.md is handled: title/subtitle,
S1..S10 headings, paragraphs, displayed equations, pipe tables (booktabs +
xltabular/longtable, with landscape pages for wide multi-column tables),
**bold**, *italic*, and `code`.

  python scripts/manuscript/build_si_pdf.py [--no-compile]
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "manuscript/supporting_information.md"
OUT = ROOT / "manuscript/jcim"
TEX_NAME = "supporting_information.tex"
PDF_NAME = "supporting_information.pdf"
LOG_NAME = "si_build.log"

# Chemical stems that benefit from discretionary hyphenation in narrow table cells.
CHEM_HYPHENS = (
    ("dimethoxyphenothiazine", r"dimethoxy\-pheno\-thiazine"),
    ("dihydrophenazine", r"dihydro\-phenazine"),
    ("phenanthrenequinone", r"phenanthrene\-quinone"),
    ("anthraquinone", r"anthra\-quinone"),
    ("naphthoquinone", r"naphtho\-quinone"),
    ("benzoquinone", r"benzo\-quinone"),
    ("pyromellitic", r"pyro\-mellitic"),
    ("methylbenzyloxy", r"methyl\-benzyl\-oxy"),
    ("methylbenzyl", r"methyl\-benzyl"),
    ("methoxyethoxy", r"methoxy\-ethoxy"),
    ("methoxyethyl", r"methoxy\-ethyl"),
    ("bipyridinium", r"bi\-pyridinium"),
    ("trifluoroethyl", r"tri\-fluoro\-ethyl"),
    ("trimethylammoniopropyl", r"trimethyl\-ammonio\-propyl"),
    ("dihydronaphthalen", r"dihydro\-naphthalen"),
    ("isoquinoline", r"iso\-quinoline"),
    ("phthalimide", r"phthal\-imide"),
    ("naphthalimide", r"naphthal\-imide"),
    ("phenothiazine", r"pheno\-thiazine"),
    ("tetrachloro", r"tetra\-chloro"),
    ("tetrafluoro", r"tetra\-fluoro"),
    ("tetramethyl", r"tetra\-methyl"),
    ("benzyloxy", r"benzyl\-oxy"),
    ("anthracene", r"anthra\-cene"),
)

# Known mathematical/physical tokens in SI prose and headers (never applied inside `code`).
MATH_TOKENS = (
    ("G(R^{q+1}) + G(R^{q−1}) − 2G(R^q)", r"$G(\mathrm{R}^{q+1}) + G(\mathrm{R}^{q-1}) - 2G(\mathrm{R}^{q})$"),
    ("ΔG_disp", r"$\Delta G_{\mathrm{disp}}$"),
    ("E_SMD", r"$E_{\mathrm{SMD}}$"),
    ("G_thermal", r"$G_{\mathrm{thermal}}$"),
    ("λ_o", r"$\lambda_{\mathrm{o}}$"),
    ("λ_i", r"$\lambda_{\mathrm{i}}$"),
    ("ε_op", r"$\varepsilon_{\mathrm{op}}$"),
    ("ε_s", r"$\varepsilon_{\mathrm{s}}$"),
    ("ε₀", r"$\varepsilon_{0}$"),
    ("n_D²", r"$n_{\mathrm{D}}^{2}$"),
    ("n_D", r"$n_{\mathrm{D}}$"),
)

# Displayed equations in S2 (matched verbatim on stripped indented lines).
DISPLAY_EQUATIONS = {
    "E_vs_Fc = −[G(R) − G(O)]/e − E_abs(Fc),  E_abs(Fc) = G(Fc⁺) − G(Fc),": (
        r"E_{\mathrm{vs\,Fc}} = -\frac{G(\mathrm{R}) - G(\mathrm{O})}{e} - E_{\mathrm{abs}}(\mathrm{Fc}),"
        r"\qquad E_{\mathrm{abs}}(\mathrm{Fc}) = G(\mathrm{Fc}^{+}) - G(\mathrm{Fc}),"
    ),
    "λ_o = z² e² (1/ε_op − 1/ε_s) / (8π ε₀ a),": (
        r"\lambda_{\mathrm{o}} = \frac{z^{2} e^{2}}{8\pi\varepsilon_{0} a}"
        r"\left(\frac{1}{\varepsilon_{\mathrm{op}}} - \frac{1}{\varepsilon_{\mathrm{s}}}\right),"
    ),
    "λ_i = [E_O(geometry of R) − E_O(geometry of O)] + [E_R(geometry of O) − E_R(geometry of R)].": (
        r"\lambda_{\mathrm{i}} = \bigl[E_{\mathrm{O}}(\text{geometry of R}) - E_{\mathrm{O}}(\text{geometry of O})\bigr]"
        r" + \bigl[E_{\mathrm{R}}(\text{geometry of O}) - E_{\mathrm{R}}(\text{geometry of R})\bigr]."
    ),
}


def esc(s: str) -> str:
    for a, b in (
        ("\\", r"\textbackslash{}"),
        ("&", r"\&"),
        ("%", r"\%"),
        ("$", r"\$"),
        ("#", r"\#"),
        ("_", r"\_"),
        ("{", r"\{"),
        ("}", r"\}"),
        ("~", r"\textasciitilde{}"),
        ("^", r"\textasciicircum{}"),
    ):
        s = s.replace(a, b)
    return s


def breakable_code(text: str) -> str:
    """Format inline code with clean line-break opportunities at path/symbol boundaries."""
    s = esc(text)
    s = s.replace("/", r"/\allowbreak{}")
    s = s.replace(r"\_", r"\_\allowbreak{}")
    s = s.replace(".", r".\allowbreak{}")
    s = s.replace("-", r"-\allowbreak{}")
    s = s.replace("(", r"(\allowbreak{}")
    s = s.replace(")", r")\allowbreak{}")
    return r"\texttt{" + s + "}"


def inline(s: str) -> str:
    codes: list[str] = []

    def keep_code(m: re.Match[str]) -> str:
        codes.append(breakable_code(m.group(1)))
        return f"\x00{len(codes) - 1}\x00"

    s = re.sub(r"`([^`]+)`", keep_code, s)

    maths: list[str] = []
    for raw, tex in MATH_TOKENS:
        if raw in s:
            idx = len(maths)
            maths.append(tex)
            s = s.replace(raw, f"\x02{idx}\x02")

    s = esc(s)
    s = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", s)
    s = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", r"\\emph{\1}", s)
    s = re.sub(r"\x00(\d+)\x00", lambda m: codes[int(m.group(1))], s)
    s = re.sub(r"\x02(\d+)\x02", lambda m: maths[int(m.group(1))], s)
    return s


def table_cell(c: str) -> str:
    """Convert a table cell, adding clean break points in long identifiers, DOIs, and names."""
    s = inline(c)
    # Allow breaks after underscores, hyphens, commas, slashes, closing parens/brackets.
    s = re.sub(r"\\_(?!\\allowbreak\{\})", r"\\_\\allowbreak{}", s)
    s = re.sub(r"(?<=[A-Za-z0-9)\]])([-,/])(?=[A-Za-z0-9(\[])", r"\1\\allowbreak{}", s)
    s = re.sub(r"(\))(?!\\allowbreak\{\})(?=[A-Za-z])", r"\1\\allowbreak{}", s)
    # Allow breaks after dots in DOIs (e.g., 10.1021/acs.chemmater.8b04770).
    if re.match(r"^10\.\d{4,9}/", c.strip()):
        s = s.replace(".", r".\allowbreak{}")
    for word, hyph in CHEM_HYPHENS:
        s = s.replace(word, hyph)
    return s


def _is_numeric_col(cells: list[str]) -> bool:
    num = re.compile(r"^[−+\-]?\d+(?:\.\d+)?(?:e[+\-]?\d+)?$|^—$|^$")
    non_empty = [c for c in cells if c]
    return bool(non_empty) and all(num.match(c) for c in non_empty)


def _is_short_col(cells: list[str], max_len: int = 12) -> bool:
    non_empty = [c for c in cells if c]
    return bool(non_empty) and all(len(c) <= max_len for c in non_empty)


def _col_spec(header: list[str], body_rows: list[list[str]], landscape: bool) -> tuple[str, bool]:
    """Choose column specs for booktabs + xltabular/longtable without overfull hboxes.

    Returns (spec_string, use_xltabular).
    """
    ncol = len(header)
    max_body = [max((len(r[k]) for r in body_rows if k < len(r)), default=0) for k in range(ncol)]
    max_any = max((max(len(header[k]), max_body[k]) for k in range(ncol)), default=0)

    # Compact narrow tables (<= 4 cols, short content) look better centered in longtable.
    if not landscape and ncol <= 4 and max_any <= 26 and sum(max(len(header[k]), max_body[k]) for k in range(ncol)) <= 62:
        cols = []
        for k in range(ncol):
            bcells = [r[k] for r in body_rows if k < len(r)]
            cols.append("r" if _is_numeric_col(bcells) else "l")
        return "@{}" + "".join(cols) + "@{}", False

    cols: list[str] = []
    x_indices: list[int] = []
    for k in range(ncol):
        h = header[k]
        bcells = [r[k] for r in body_rows if k < len(r)]
        numeric = _is_numeric_col(bcells)
        if numeric:
            if len(h) <= 12 or " " not in h:
                cols.append("r")
            elif len(h) <= 18:
                width = "1.55cm" if landscape else "1.95cm"
                cols.append(rf">{{\raggedleft\arraybackslash}}p{{{width}}}")
            else:
                width = "2.2cm" if landscape else "2.45cm"
                cols.append(rf">{{\raggedleft\arraybackslash}}p{{{width}}}")
        elif k == 0 and (max_body[k] >= 20 or ncol >= 5):
            cols.append(r">{\raggedright\arraybackslash}X")
            x_indices.append(k)
        elif h == "Source DOI":
            cols.append(r">{\raggedright\arraybackslash}p{4.9cm}")
        elif h == "QC flag":
            cols.append(r">{\raggedright\arraybackslash}p{2.9cm}")
        elif h == "Note":
            cols.append(r">{\raggedright\arraybackslash}p{1.35cm}")
        elif h == "GPU":
            cols.append(r">{\raggedright\arraybackslash}p{3.6cm}")
        elif h == "SCF guesses":
            cols.append(r">{\raggedright\arraybackslash}p{2.1cm}")
        elif h == "Front":
            cols.append("c")
        elif _is_short_col(bcells, max_len=11) and len(h) <= 11:
            cols.append("l")
        elif max_body[k] <= 24 and len(h) <= 24:
            cols.append("l")
        else:
            cols.append(r">{\raggedright\arraybackslash}X")
            x_indices.append(k)

    # xltabular requires at least one X column.
    if not x_indices:
        cols[0] = r">{\raggedright\arraybackslash}X"
    return "@{}" + "".join(cols) + "@{}", True


def render_table(lines: list[str], preamble_blocks: list[str] | None = None, table_label: str | None = None) -> str:
    """Render a Markdown pipe table as booktabs + xltabular/longtable (landscape if wide)."""
    rows = [
        [c.strip() for c in line.strip().strip("|").split("|")]
        for line in lines
        if not re.match(r"^\|\s*[-:]+", line.strip())
    ]
    header = rows[0]
    body_rows = rows[1:]
    ncol = len(header)
    landscape = ncol >= 9
    spec, use_xltabular = _col_spec(header, body_rows, landscape)

    hdr_tex = " & ".join(r"\textbf{" + table_cell(c) + "}" for c in header) + r" \\"
    body_tex = "\n".join(" & ".join(table_cell(c) for c in r) + r" \\" for r in body_rows)

    cont_label = f"{table_label} (continued)" if table_label else "Continued from previous page"
    size = r"\footnotesize" if ncol >= 5 else r"\small"
    tabcolsep = "3.5pt" if landscape else ("4pt" if ncol >= 7 else "5pt")

    env_open = (
        rf"\setlength{{\LTleft}}{{0pt}}\setlength{{\LTright}}{{0pt}}"
        "\n"
        rf"\begin{{xltabular}}{{\linewidth}}{{{spec}}}"
        if use_xltabular
        else rf"\setlength{{\LTleft}}{{\fill}}\setlength{{\LTright}}{{\fill}}"
        "\n"
        rf"\begin{{longtable}}{{{spec}}}"
    )
    env_close = r"\end{xltabular}" if use_xltabular else r"\end{longtable}"

    table_tex = (
        "{\n"
        f"{size}\n"
        f"\\setlength{{\\tabcolsep}}{{{tabcolsep}}}\n"
        "\\renewcommand{\\arraystretch}{1.12}\n"
        f"{env_open}\n"
        "\\toprule\n"
        f"{hdr_tex}\n"
        "\\midrule\n"
        "\\endfirsthead\n"
        f"\\multicolumn{{{ncol}}}{{@{{}}l@{{}}}}{{\\emph{{{cont_label}}}}} \\\\\n"
        "\\toprule\n"
        f"{hdr_tex}\n"
        "\\midrule\n"
        "\\endhead\n"
        "\\midrule\n"
        f"\\multicolumn{{{ncol}}}{{@{{}}r@{{}}}}{{\\emph{{Continued on next page}}}} \\\\\n"
        "\\endfoot\n"
        "\\bottomrule\n"
        "\\endlastfoot\n"
        f"{body_tex}\n"
        f"{env_close}\n"
        "}\n"
    )

    if landscape:
        intro = "\n".join(preamble_blocks or [])
        return "\\begin{landscape}\n" + intro + table_tex + "\\end{landscape}\n"
    return ("".join(preamble_blocks or [])) + table_tex


def convert(md_text: str | None = None) -> str:
    """Convert supporting_information.md text into a complete XeLaTeX document."""
    md = SRC.read_text() if md_text is None else md_text
    lines = md.splitlines()

    subtitle = ""
    i = 0
    if i < len(lines) and lines[i].startswith("# "):
        i += 1
    if i < len(lines) and lines[i].strip().startswith("*") and lines[i].strip().endswith("*"):
        subtitle = lines[i].strip().strip("*").rstrip(".")
        i += 1

    # Collect intro paragraph(s) before first ## section.
    intro_paras: list[str] = []
    para: list[str] = []
    while i < len(lines) and not lines[i].startswith("## "):
        s = lines[i].strip()
        if not s:
            if para:
                intro_paras.append(inline(" ".join(para)))
                para = []
        else:
            para.append(s)
        i += 1
    if para:
        intro_paras.append(inline(" ".join(para)))
        para = []

    body_out: list[str] = []
    pending_table_blocks: list[str] = []
    pending_table_label: str | None = None

    def flush_para() -> None:
        nonlocal para, pending_table_label
        if not para:
            return
        text = " ".join(x.strip() for x in para)
        para = []
        m = re.match(r"^\*\*(Table S\d+[a-z]?)\.\*\*", text)
        if m:
            # If the immediately preceding item is a bare \section{...}, keep it attached
            # so a landscape table never leaves an orphaned section heading on the prior page.
            if body_out and body_out[-1].startswith("\\section{"):
                pending_table_blocks.insert(0, body_out.pop())
            pending_table_label = m.group(1)
            pending_table_blocks.append("\\noindent " + inline(text) + "\n\n")
            return
        if pending_table_blocks and text.startswith("**Independent check.**"):
            pending_table_blocks.append("\\noindent " + inline(text) + "\n\n")
            return
        if pending_table_blocks:
            body_out.extend(pending_table_blocks)
            pending_table_blocks.clear()
            pending_table_label = None
        body_out.append(inline(text) + "\n\n")

    while i < len(lines):
        raw = lines[i]
        s = raw.strip()

        if s.startswith("## "):
            flush_para()
            if pending_table_blocks:
                body_out.extend(pending_table_blocks)
                pending_table_blocks.clear()
                pending_table_label = None
            head = s[3:].strip()
            head = re.sub(r"^S\d+\.\s*", "", head)
            body_out.append(f"\\section{{{inline(head)}}}\n")
            i += 1
            continue

        # Displayed equations (two-space indented lines in S2).
        if raw.startswith("  ") and s in DISPLAY_EQUATIONS:
            flush_para()
            body_out.append("\\[\n  " + DISPLAY_EQUATIONS[s] + "\n\\]\n")
            i += 1
            continue

        if s.startswith("|"):
            flush_para()
            tlines: list[str] = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                tlines.append(lines[i])
                i += 1
            body_out.append(
                render_table(tlines, list(pending_table_blocks), pending_table_label)
            )
            pending_table_blocks.clear()
            pending_table_label = None
            continue

        if not s:
            flush_para()
            i += 1
            continue

        para.append(raw)
        i += 1

    flush_para()
    if pending_table_blocks:
        body_out.extend(pending_table_blocks)

    title_tex = inline(subtitle) if subtitle else "Supporting Information"
    intro_tex = "\n\n".join(intro_paras)

    doc = (
        r"""\documentclass[10pt,letterpaper]{article}
\usepackage[margin=0.85in,headheight=14pt]{geometry}
\usepackage{fontspec}
\setmainfont{DejaVu Serif}
\setsansfont{DejaVu Sans}
\setmonofont{DejaVu Sans Mono}[Scale=0.82]
\usepackage{amsmath}
\usepackage{array}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{xltabular}
\usepackage{pdflscape}
\usepackage{fancyhdr}
\usepackage[hyphens]{xurl}
\usepackage[colorlinks=true,linkcolor=black,urlcolor=black,citecolor=black]{hyperref}

\renewcommand{\thepage}{S\arabic{page}}
\renewcommand{\thesection}{S\arabic{section}}
\renewcommand{\contentsname}{Contents}
\makeatletter
\renewcommand*\l@section[2]{%
  \ifnum \c@tocdepth >\z@
    \addpenalty\@secpenalty
    \addvspace{0.35em \@plus\p@}%
    \setlength\@tempdima{2.4em}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      \leavevmode \bfseries
      \advance\leftskip\@tempdima
      \hskip -\leftskip
      #1\nobreak\hfil \nobreak\hb@xt@\@pnumwidth{\hss #2}\par
    \endgroup
  \fi}
\renewcommand\section{\@startsection{section}{1}{\z@}%
  {-2.5ex \@plus -1ex \@minus -.2ex}%
  {1.2ex \@plus .2ex}%
  {\normalfont\large\bfseries\raggedright}}
\makeatother
\setlength{\parindent}{0pt}
\setlength{\parskip}{5pt}
\setlength{\LTpre}{6pt}
\setlength{\LTpost}{8pt}
\emergencystretch=1.5em

\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\footnotesize\sffamily Supporting Information}
\fancyhead[R]{\footnotesize\sffamily Graftable Redox Motifs for Non-Aqueous ORFBs}
\fancyfoot[C]{\footnotesize\sffamily\thepage}
\renewcommand{\headrulewidth}{0.4pt}
\renewcommand{\footrulewidth}{0pt}

\fancypagestyle{plain}{
  \fancyhf{}
  \fancyfoot[C]{\footnotesize\sffamily\thepage}
  \renewcommand{\headrulewidth}{0pt}
  \renewcommand{\footrulewidth}{0pt}
}

\begin{document}
\thispagestyle{plain}

\begin{center}
{\sffamily\bfseries\large SUPPORTING INFORMATION}\\[6pt]
\rule{\textwidth}{0.8pt}\\[10pt]
{\Large\bfseries """
        + title_tex
        + r"""}\\[10pt]
{\normalsize [TODO: authors]}\\[4pt]
{\small\itshape [TODO: affiliation]}\\[2pt]
{\small Email: [TODO: corresponding author email]}\\[8pt]
\rule{\textwidth}{0.4pt}
\end{center}

"""
        + intro_tex
        + r"""

\tableofcontents
\vspace{6pt}
\noindent\rule{\textwidth}{0.4pt}

"""
        + "".join(body_out)
        + r"""\end{document}
"""
    )
    return doc


def compile_document() -> None:
    """Compile all XeLaTeX passes, refusing failed commands or a stale/missing PDF."""
    pdf = OUT / PDF_NAME
    pdf.unlink(missing_ok=True)
    log_path = OUT / LOG_NAME
    cmd = ["xelatex", "-halt-on-error", "-interaction=nonstopmode", TEX_NAME]
    with log_path.open("w") as log:
        for _ in range(3):
            log.write("$ " + " ".join(cmd) + "\n")
            log.flush()
            try:
                subprocess.run(cmd, cwd=OUT, stdout=log, stderr=subprocess.STDOUT, check=True)
            except (OSError, subprocess.CalledProcessError) as exc:
                pdf.unlink(missing_ok=True)
                raise RuntimeError(f"LaTeX build failed at {cmd[0]}; see {log_path}") from exc
    if not pdf.is_file() or pdf.stat().st_size == 0:
        pdf.unlink(missing_ok=True)
        raise RuntimeError(f"LaTeX build produced no PDF; see {log_path}")
    print("compiled:", pdf)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    doc = convert()
    tex_path = OUT / TEX_NAME
    tex_path.write_text(doc)
    print("wrote", tex_path)
    if "--no-compile" in sys.argv:
        return 0
    compile_document()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
