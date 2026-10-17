#!/usr/bin/env python3
"""Export a notebook to a self-contained HTML whose equations always render.

    python3 nb2html.py A1.ipynb            # writes A1.html

`jupyter nbconvert --to html` leaves the LaTeX in the page and loads MathJax
from the internet to typeset it. Canvas serves uploaded HTML files sandboxed,
so that script never runs and the reader sees raw $...$ instead of formulas.

This script runs nbconvert and then rewrites every formula as MathML, which
browsers render themselves — no JavaScript, no internet, nothing to block.
Plots are embedded in the file, so the result is one self-contained page.

Requires nbconvert and latex2mathml:

    pip install nbconvert latex2mathml
"""
import re
import subprocess
import sys
from html import unescape
from pathlib import Path

try:
    from latex2mathml.converter import convert as latex_to_mathml
except ImportError:
    sys.exit("nb2html: missing dependency — run:  pip install latex2mathml")

# Regions whose text is code or markup, never prose: their '$' are not math.
SKIP_TAGS = {"script", "style", "pre", "code", "textarea"}

# $$...$$ / \[...\] display math, $...$ / \(...\) inline math. Inline math may
# not span a blank line and may not start or end on whitespace, which keeps
# stray dollar signs in prose from being swallowed.
MATH = re.compile(
    r"\$\$(?P<d1>.+?)\$\$"
    r"|\\\[(?P<d2>.+?)\\\]"
    r"|\$(?P<i1>[^\s$](?:[^$]*[^\s$])?)\$"
    r"|\\\((?P<i2>.+?)\\\)",
    re.S,
)

TAG = re.compile(r"<(/?)([a-zA-Z0-9]+)[^>]*>|<!--.*?-->", re.S)

# Bold vectors are the Unicode math-alphanumeric characters (U+1D400 block),
# so name fonts that carry them before falling back to a generic serif.
MATH_CSS = """
<style>
math {
  font-size: 1.08em;
  font-family: "Latin Modern Math", "STIX Two Math", "Cambria Math",
               "DejaVu Serif", serif;
}
/* a block formula is its own line, centred: 'margin: auto' only centres a
   box that has a width, hence fit-content */
math[display="block"] {
  display: block;
  width: fit-content;
  max-width: 100%;
  overflow-x: auto;
  margin: 0.9em auto;
}
</style>
"""


def convert_snippet(latex, display):
    """LaTeX -> MathML; returns None if latex2mathml cannot handle it."""
    try:
        mathml = latex_to_mathml(unescape(latex).strip())
    except Exception:
        return None
    if display:
        mathml = mathml.replace('display="inline"', 'display="block"', 1)
    return mathml


def convert_text(text, stats):
    """Replace every math snippet in one HTML text node."""
    def repl(m):
        display = m.group("d1") is not None or m.group("d2") is not None
        latex = next(g for g in m.groups() if g is not None)
        mathml = convert_snippet(latex, display)
        if mathml is None:
            stats["failed"].append(latex.strip()[:60])
            return m.group(0)               # leave the source untouched
        stats["converted"] += 1
        return mathml
    return MATH.sub(repl, text)


def rewrite_math(html, stats):
    """Walk the document, converting math only in prose text nodes."""
    out, pos, depth = [], 0, 0
    for tag in TAG.finditer(html):
        text = html[pos:tag.start()]
        out.append(convert_text(text, stats) if depth == 0 else text)
        out.append(tag.group(0))
        pos = tag.end()
        name = (tag.group(2) or "").lower()
        if name in SKIP_TAGS:
            depth += -1 if tag.group(1) else 1
            depth = max(depth, 0)
    tail = html[pos:]
    out.append(convert_text(tail, stats) if depth == 0 else tail)
    return "".join(out)


def strip_mathjax(html):
    """Remove the MathJax loader — the formulas are MathML now. Other scripts
    (require.js, an interactive plotting library) are left alone."""
    return re.sub(
        r"<script[^>]*\bsrc=\"https?://[^\"]*mathjax[^\"]*\"[^>]*>.*?</script>",
        "", html, flags=re.S | re.I)


def export(notebook):
    notebook = Path(notebook)
    html_path = notebook.with_suffix(".html")
    cmd = [sys.executable, "-m", "nbconvert", "--to", "html",
           "--embed-images", str(notebook)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"nb2html: nbconvert failed:\n{result.stderr}")

    html = html_path.read_text(encoding="utf-8")
    stats = {"converted": 0, "failed": []}
    html = strip_mathjax(html)
    html = rewrite_math(html, stats)
    html = html.replace("</head>", MATH_CSS + "</head>", 1)
    html_path.write_text(html, encoding="utf-8")

    print(f"{notebook.name}: {html_path.name}, {stats['converted']} formulas "
          f"as MathML, {len(stats['failed'])} left as text")
    for snippet in stats["failed"]:
        print(f"  could not convert: {snippet}")
    return 0


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    for notebook in argv[1:]:
        export(notebook)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
