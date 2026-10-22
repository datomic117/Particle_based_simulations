#!/usr/bin/env python3
"""Export A4.ipynb to a self-contained HTML with collapsed code and MathML formulas.

    python make_html.py A4.ipynb            # writes A4.html

Same job as the course script nb2html.py (which stays untouched and whose
functions are reused here), with one difference: the notebook is exported with
the template in html_template/collapse, which wraps every code cell tagged
`hide-input` (or with metadata.jupyter.source_hidden) in a closed <details>
element, because nbconvert 7.17 ignores source_hidden. The code is one click
away ("Show code"), the outputs stay visible.
"""
import subprocess
import sys
from pathlib import Path

from nb2html import MATH_CSS, rewrite_math, strip_mathjax

HERE = Path(__file__).resolve().parent

# details.hidecode: the collapsed code of a cell; details.oral: the Oral prep blocks
EXTRA_CSS = """
<style>
details.hidecode > summary {
  cursor: pointer;
  color: #777;
  font-size: 0.8em;
  margin: 0.2em 0 0.2em 0.4em;
}
details.hidecode[open] > summary { margin-bottom: 0.4em; }
details > summary > b { cursor: pointer; }
details:not(.hidecode) > summary { cursor: pointer; margin: 0.4em 0; }
</style>
"""


ORAL_OPEN = '<details>\n<summary><b>Oral prep</b> (click to open)</summary>\n</details>\n'


def fix_oral(html):
    """nbconvert's markdown renderer ends a raw <details> block at the first blank line and drops the closing tag, so
    the Oral prep list would end up outside its <details>. Re-open the block and close it at the end of its cell."""
    out, pos, n = [], 0, 0
    while True:
        i = html.find(ORAL_OPEN, pos)
        if i < 0:
            break
        end = html.index("\n</div>", i)              # the end of the markdown output that holds the block
        out += [html[pos:i], '<details class="oral">\n<summary><b>Oral prep</b> (click to open)</summary>\n', html[i + len(ORAL_OPEN):end], "\n</details>"]
        pos, n = end, n + 1
    out.append(html[pos:])
    return "".join(out), n


def export(notebook):
    notebook = Path(notebook)
    html_path = notebook.with_suffix(".html")
    cmd = [sys.executable, "-m", "nbconvert", "--to", "html", "--embed-images",
           "--template", "collapse",
           f"--TemplateExporter.extra_template_basedirs={HERE / 'html_template'}",
           str(notebook)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"make_html: nbconvert failed:\n{result.stderr}")

    html = html_path.read_text(encoding="utf-8")
    stats = {"converted": 0, "failed": []}
    html = strip_mathjax(html)
    html, n_oral = fix_oral(html)
    html = rewrite_math(html, stats)
    html = html.replace("</head>", MATH_CSS + EXTRA_CSS + "</head>", 1)
    html_path.write_text(html, encoding="utf-8", newline="\n")

    print(f"{notebook.name}: {html_path.name}, {stats['converted']} formulas "
          f"as MathML, {len(stats['failed'])} left as text, "
          f"{html.count('<details class=\"hidecode\">')} collapsed code cells, {n_oral} Oral prep blocks "
          f"({html.count('<details')} <details> and {html.count('</details>')} </details> in all)")
    for snippet in stats["failed"]:
        print(f"  could not convert: {snippet}")
    return 0


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    for nb in sys.argv[1:]:
        export(nb)
