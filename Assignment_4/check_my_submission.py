#!/usr/bin/env python3
"""Check your own 6EMA02 hand-in zip before you submit it.

    python3 check_my_submission.py group7_A1.zip --release A1.ipynb

`--release` is a freshly downloaded, unmodified copy of the assignment
notebook; it is what your task cells are compared against. Leave it out and
that one check is skipped.

Checked: the zip layout, that the task cells are untouched, that every code
cell has been executed, that no binaries are shipped, that your C code
compiles, that the HTML export is the one that displays formulas in Canvas,
and that every path and link in the notebook is relative and present in the
zip. Nothing here is graded: it is the same layout check the teacher runs, so
that no submission is refused for a reason you could have fixed in a minute.

Only the Python standard library is needed (and gcc for the compile check).
"""
import argparse
import json
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile

PASS, WARN, FAIL = "ok  ", "note", "FAIL"

# Strings in code cells that look like a file the notebook opens.
DATA_LIKE = re.compile(
    r"""["']([^"'\n]*\.(?:csv|dat|txt|tsv|json|npy|npz|pdb|xyz|png|jpg))["']""",
    re.I)
# Markdown links, e.g. [compute_acc](code/nbody.c)
MD_LINK = re.compile(r"\[[^\]]*\]\(([^)\s]+)")
ABSOLUTE = re.compile(r"""["'](?:[A-Za-z]:[\\/]|/(?:home|Users|mnt|tmp)/|file://)""")


class Report:
    def __init__(self):
        self.rows = []

    def add(self, status, text):
        self.rows.append((status, text))

    def failed(self):
        return any(s == FAIL for s, _ in self.rows)

    def show(self):
        for status, text in self.rows:
            print(f"  [{status}] {text}")
        if self.failed():
            print("\nNot ready to submit: fix the FAIL lines above.")
        else:
            notes = sum(s == WARN for s, _ in self.rows)
            print(f"\nLooks good{f' ({notes} notes to consider)' if notes else ''}. "
                  "You can submit this zip.")
        return 1 if self.failed() else 0


def cells_of(nb, kind=None):
    for c in nb.get("cells", []):
        if kind is None or c.get("cell_type") == kind:
            yield c


def task_cells(nb):
    return [(c["cell_type"], "".join(c.get("source", [])))
            for c in cells_of(nb)
            if "task" in c.get("metadata", {}).get("tags", [])]


def looks_binary(path):
    head = path.read_bytes()[:4]
    return head[:4] == b"\x7fELF" or head[:2] == b"MZ"


def check_layout(root, name, rep):
    nb = root / f"{name}.ipynb"
    html = root / f"{name}.html"
    for p, label in [(nb, f"{name}.ipynb"), (html, f"{name}.html"),
                     (root / "code", "code/"), (root / "data", "data/")]:
        rep.add(PASS if p.exists() else FAIL,
                f"{label} {'present' if p.exists() else 'MISSING from the zip'}")
    return nb, html


def check_notebook(nb_path, release, root, rep):
    try:
        nb = json.loads(nb_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as e:
        rep.add(FAIL, f"cannot read the notebook: {e}")
        return

    if release:
        try:
            ref = task_cells(json.loads(pathlib.Path(release).read_text(encoding="utf-8")))
        except (OSError, json.JSONDecodeError) as e:
            rep.add(WARN, f"cannot read --release notebook ({e}); task check skipped")
        else:
            mine = task_cells(nb)
            if mine == ref:
                rep.add(PASS, f"task cells unmodified ({len(ref)} cells)")
            else:
                rep.add(FAIL, f"TASK CELLS MODIFIED ({len(mine)} here vs {len(ref)} in "
                              "the release): restore them from the downloaded notebook")
    else:
        rep.add(WARN, "no --release given, so the task cells were not checked")

    code = [c for c in cells_of(nb, "code") if "".join(c.get("source", [])).strip()]
    run = [c for c in code if c.get("outputs") or c.get("execution_count")]
    rep.add(PASS if len(run) == len(code) else FAIL,
            f"{len(run)}/{len(code)} code cells executed"
            + ("" if len(run) == len(code) else "; run all cells, then save"))

    # every path the notebook mentions must be relative and inside the zip
    missing, absolute = [], []
    for c in cells_of(nb):
        src = "".join(c.get("source", []))
        if c.get("cell_type") == "code":
            if ABSOLUTE.search(src):
                absolute.append(ABSOLUTE.search(src).group(0).strip("\"'"))
            targets = DATA_LIKE.findall(src)
        else:
            targets = [t for t in MD_LINK.findall(src)
                       if not t.startswith(("http://", "https://", "mailto:", "#"))]
            absolute += [t for t in targets if t.startswith(("/", "~"))
                         or re.match(r"[A-Za-z]:[\\/]", t)]
        for t in targets:
            if t.startswith(("http", "/", "~")) or re.match(r"[A-Za-z]:[\\/]", t):
                continue
            if "{" in t or "}" in t or "%" in t:
                continue        # f-string or format template, not a literal path
            if not (root / t.split("#")[0]).exists():
                missing.append(t)

    if absolute:
        rep.add(FAIL, "absolute paths (they only work on your machine): "
                      + ", ".join(sorted(set(absolute))[:4]))
    else:
        rep.add(PASS, "no absolute paths")
    if missing:
        rep.add(FAIL, "referenced but not in the zip: " + ", ".join(sorted(set(missing))[:6]))
    else:
        rep.add(PASS, "every file and link the notebook references is in the zip")

    links = sum(len([t for t in MD_LINK.findall("".join(c.get("source", [])))
                     if not t.startswith(("http", "#"))])
                for c in cells_of(nb, "markdown"))
    rep.add(PASS if links else WARN,
            f"{links} clickable links to your files"
            + ("" if links else "; the code map is easier to use in the oral with links"))


def check_html(html_path, rep):
    if not html_path.exists():
        return
    text = html_path.read_text(encoding="utf-8", errors="replace")
    if "<math" in text:
        rep.add(PASS, "HTML formulas are MathML (they display in Canvas)")
    elif re.search(r"\$[^\s$][^$]{0,80}\$", re.sub(r"<(script|style|pre|code)\b.*?</\1>",
                                                   "", text, flags=re.S)):
        rep.add(WARN, "HTML still holds raw $...$: export it with "
                      "`python3 nb2html.py <notebook>` so formulas display in Canvas")
    if "base64" not in text:
        rep.add(WARN, "HTML seems to contain no embedded figures")


def check_code(code_dir, rep):
    if not code_dir.is_dir():
        return
    binaries = [p.name for p in code_dir.rglob("*")
                if p.is_file() and (looks_binary(p) or p.suffix.lower() == ".exe")]
    if binaries:
        rep.add(FAIL, "binaries in code/ (hand in sources only): "
                      + ", ".join(sorted(binaries)[:4]))
    if not (code_dir / "README.md").exists():
        rep.add(WARN, "code/README.md missing (build and run instructions belong there)")
    sources = sorted(str(p) for p in code_dir.glob("*.c"))
    if not sources:
        rep.add(FAIL, "no .c files in code/")
        return
    if not shutil.which("gcc"):
        rep.add(WARN, "gcc not found, so the compile check was skipped")
        return
    with tempfile.TemporaryDirectory() as tmp:
        r = subprocess.run(["gcc", "-O2", "-Wall", "-Wextra",
                            "-o", str(pathlib.Path(tmp) / "a.out")] + sources + ["-lm"],
                           capture_output=True, text=True)
    if r.returncode != 0:
        rep.add(FAIL, f"code/ does not compile ({r.stderr.count('error:')} errors); "
                      "first: " + next((l for l in r.stderr.splitlines()
                                        if "error:" in l), "")[:90])
    else:
        n = r.stderr.count("warning:")
        rep.add(PASS if n == 0 else WARN,
                f"code/ compiles with {n} warning{'s' if n != 1 else ''}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("zipfile", help="your submission zip")
    ap.add_argument("--release", help="freshly downloaded A<n>.ipynb")
    args = ap.parse_args()

    zpath = pathlib.Path(args.zipfile)
    rep = Report()
    print(f"Checking {zpath.name}\n")
    if not zpath.exists():
        sys.exit(f"no such file: {zpath}")

    with tempfile.TemporaryDirectory() as tmp:
        try:
            with zipfile.ZipFile(zpath) as z:
                z.extractall(tmp)
        except zipfile.BadZipFile:
            sys.exit("this is not a readable zip file")
        root = pathlib.Path(tmp)
        entries = [p for p in root.iterdir() if not p.name.startswith("__MACOSX")]
        if len(entries) == 1 and entries[0].is_dir():
            root = entries[0]

        names = sorted({p.stem for p in root.glob("A?.ipynb")})
        if not names:
            sys.exit("no A<n>.ipynb found in the zip: is this the right file?")
        name = names[0]

        nb_path, html_path = check_layout(root, name, rep)
        if nb_path.exists():
            check_notebook(nb_path, args.release, root, rep)
        check_html(html_path, rep)
        check_code(root / "code", rep)

    return rep.show()


if __name__ == "__main__":
    sys.exit(main())
