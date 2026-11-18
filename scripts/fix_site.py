#!/usr/bin/env python3
"""Post-process the pkgdown site of htheoryR before publishing it on GitHub Pages.

1. pkgdown 2.2.1 writes ``type="”image/svg+xml”"`` (curly quotes) on the SVG
   favicon link of every page; browsers then ignore the SVG icon.
2. GitHub Pages serves the site under ``/hierarchical-depth/``, so the
   ``site.webmanifest`` icons must use relative paths (fixed in
   ``htheoryR/pkgdown/favicon/site.webmanifest``; checked here).
3. ``.nojekyll`` keeps GitHub Pages from running Jekyll over the site.

Usage: python3 scripts/fix_site.py htheoryR/docs
"""

import json
import sys
from pathlib import Path

BAD = 'type="”image/svg+xml”"'
GOOD = 'type="image/svg+xml"'


def main(site: Path) -> None:
    fixed = 0
    for html in site.rglob("*.html"):
        text = html.read_text(encoding="utf-8")
        if BAD in text:
            html.write_text(text.replace(BAD, GOOD), encoding="utf-8")
            fixed += 1
    manifest = json.loads((site / "site.webmanifest").read_text())
    absolute = [i["src"] for i in manifest["icons"] if i["src"].startswith("/")]
    if absolute:
        sys.exit(f"site.webmanifest has absolute icon paths: {absolute}")
    (site / ".nojekyll").touch()
    print(f"{fixed} pages fixed; manifest icons relative; .nojekyll written")


if __name__ == "__main__":
    main(Path(sys.argv[1] if len(sys.argv) > 1 else "htheoryR/docs"))
