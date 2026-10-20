"""Fail if any task cell of A4.ipynb differs from the release.  Run from Assignment_4/:  $PY agent_prompts/verify_tasks.py"""
import json, sys
def tasks(p):
    nb = json.load(open(p, encoding="utf8"))
    return [(c["cell_type"], "".join(c["source"])) for c in nb["cells"]
            if "task" in c.get("metadata", {}).get("tags", [])]
mine, ref = tasks("A4.ipynb"), tasks("agent_prompts/A4_release.ipynb")
bad = [i for i, (a, b) in enumerate(zip(mine, ref)) if a != b]
if len(mine) != len(ref) or bad:
    print("TASK CELLS MODIFIED", len(mine), len(ref), bad); sys.exit(1)
print(f"OK: {len(ref)} task cells identical to the release")
