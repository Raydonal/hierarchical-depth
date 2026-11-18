#!/usr/bin/env bash
# Publish htheoryR/docs (built by `make site`) to the gh-pages branch, which
# GitHub Pages serves at https://raydonal.github.io/hierarchical-depth/.
# The branch holds only the built site and is replaced on every publish.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SITE="$ROOT/htheoryR/docs"
REMOTE="$(git -C "$ROOT" remote get-url origin)"

[ -f "$SITE/index.html" ] || { echo "no site in $SITE: run make site first" >&2; exit 1; }
[ -f "$SITE/.nojekyll" ] || { echo "site not post-processed: run make site" >&2; exit 1; }

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
cp -a "$SITE/." "$TMP/"
git -C "$TMP" init -q -b gh-pages
git -C "$TMP" add -A
git -C "$TMP" commit -q -m "pkgdown site for htheoryR ($(git -C "$ROOT" rev-parse --short HEAD))"
git -C "$TMP" push -f "$REMOTE" gh-pages
