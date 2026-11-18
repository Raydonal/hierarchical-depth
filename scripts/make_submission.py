#!/usr/bin/env python3
"""Assemble the flat submission bundle for Elsevier's Editorial Manager.

Editorial Manager places every uploaded file in a single folder, while
``article/main.tex`` reads figures from ``../figures/`` and tables from
``generated/``. This script writes a copy of the manuscript with local figure
paths and the generated tables inlined, so the working ``main.tex`` stays the
single source. Output, in ``article/submission/``:

    main.tex, references.bib, main.bbl, main.pdf   manuscript (compiled to check)
    fig*.pdf                                       the figures main.tex uses
    highlights.txt                                 one highlight per line
    cover-letter.pdf                               from article/cover-letter.md

Run from the repo root after ``make tables figures``::

    PYTHONPATH=. python3 scripts/make_submission.py
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ARTICLE = ROOT / "article"
OUT = ARTICLE / "submission"
HIGHLIGHT_MAX = 85  # Elsevier: at most 85 characters per highlight, spaces included


def flatten(tex: str) -> tuple[str, list[str]]:
    figures = re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{\.\./figures/([^}]+)\}", tex)
    tex = tex.replace("{../figures/", "{")

    def inline(m: re.Match) -> str:
        return (ARTICLE / m.group(1)).read_text().rstrip("\n")

    tex = re.sub(r"\\input\{(generated/[^}]+\.tex)\}", inline, tex)
    return tex, figures


def highlights(tex: str) -> list[str]:
    block = re.search(r"\\begin\{highlights\}(.*?)\\end\{highlights\}", tex, re.S)
    if block is None:
        sys.exit("main.tex has no highlights environment")
    return [h.strip() for h in re.findall(r"\\item\s+(.*)", block.group(1))]


def run(cmd: list[str], cwd: Path) -> None:
    r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.stdout.write(r.stdout[-3000:] + r.stderr[-3000:])
        sys.exit(f"failed: {' '.join(cmd)}")


def main() -> None:
    src = (ARTICLE / "main.tex").read_text()
    tex, figures = flatten(src)

    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir()
    (OUT / "main.tex").write_text(tex)
    shutil.copy(ARTICLE / "references.bib", OUT)
    for f in figures:
        shutil.copy(ROOT / "figures" / f, OUT)

    hl = highlights(src)
    for h in hl:
        if len(h) > HIGHLIGHT_MAX:
            sys.exit(f"highlight longer than {HIGHLIGHT_MAX} characters: {h!r}")
    (OUT / "highlights.txt").write_text("\n".join(hl) + "\n")

    run(["latexmk", "-pdf", "-interaction=nonstopmode", "main.tex"], OUT)
    log = (OUT / "main.log").read_text(errors="replace")
    problems = [l for l in log.splitlines()
                if "undefined" in l or "Float too large" in l or "LaTeX Error" in l]
    run(["latexmk", "-c", "main.tex"], OUT)  # keeps main.pdf and main.bbl
    for leftover in OUT.glob("*.spl"):  # elsarticle's auxiliary, unknown to latexmk
        leftover.unlink()

    letter = ARTICLE / "cover-letter.md"
    if letter.exists():
        run(["pandoc", str(letter), "--pdf-engine=xelatex", "-o",
             str(OUT / "cover-letter.pdf")], ARTICLE)

    print(f"bundle: {OUT.relative_to(ROOT)}  ({len(figures)} figures, "
          f"{len(hl)} highlights)")
    for p in problems:
        print("WARNING:", p)
    if "\\todo{" in tex.split("\\begin{document}", 1)[1]:
        print("WARNING: main.tex still contains \\todo{...}; do not upload yet")


if __name__ == "__main__":
    main()
