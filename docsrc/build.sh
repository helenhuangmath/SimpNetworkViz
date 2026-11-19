#!/usr/bin/env bash
# Build the documentation site from docsrc/ into docs/, which GitHub Pages serves
# from the gh-pages branch. Run from anywhere:  bash docsrc/build.sh
#
# docs/images/ is left alone: the README and the pages link to those figures.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(dirname "$here")"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

sphinx-build -b html -W --keep-going -d "$tmp/doctrees" "$here" "$tmp/html"
rm -rf "$tmp/html/.buildinfo" "$tmp/html/.doctrees"

mkdir -p "$root/docs"
rsync -a --delete --exclude 'images/' "$tmp/html/" "$root/docs/"
echo "Site written to $root/docs/ (open docs/index.html)"
