#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="$ROOT_DIR/build"
mkdir -p "$BUILD_DIR"

TMP_MD="$BUILD_DIR/dopamine_manual_tmp.md"

# Compose manual from sections in deterministic order.
cat > "$TMP_MD" <<'EOF'
---
title: "DOPAMINE: User and Developer Manual"
author: "DOPAMINE Project"
link-citations: true
numbersections: true
reference-section-title: "References"
geometry: margin=1in
---
EOF

# Strip each page's Jekyll front matter (--- title: ... ---) before
# appending: pandoc would otherwise parse it as a second YAML metadata
# block and overwrite the manual's own title/author/date above.
strip_front_matter() {
  awk 'BEGIN{n=0} /^---$/{n++; if(n<=2) next} {if(n>=2) print}' "$1"
}

# Jekyll uses $$..$$ for inline math (pandoc reads it as display math) and
# {{ site.baseurl }} in asset paths; rewrite both for the PDF/HTML build.
jekyll_to_pandoc() {
  perl -pe '
    s/\{\{ site\.baseurl \}\}\///g;
    unless (/^\s*\$\$.*\$\$\s*$/ && (() = /\$\$/g) == 2) { s/\$\$(.+?)\$\$/\$$1\$/g }
  '
}

for page in installation usage preprocessing input-reference theory validation examples tools development; do
  strip_front_matter "$ROOT_DIR/$page.md" | jekyll_to_pandoc >> "$TMP_MD"
  printf "\n\n" >> "$TMP_MD"
done

CSL_FILE=""
if [[ -f /usr/share/pandoc/data/csl/ieee.csl ]]; then
  CSL_FILE="/usr/share/pandoc/data/csl/ieee.csl"
elif [[ -f /usr/local/share/pandoc/data/csl/ieee.csl ]]; then
  CSL_FILE="/usr/local/share/pandoc/data/csl/ieee.csl"
fi

CSL_ARGS=()
if [[ -n "$CSL_FILE" ]]; then
  CSL_ARGS=(--csl "$CSL_FILE")
fi

pandoc "$TMP_MD" \
  --from markdown+tex_math_dollars \
  --to html5 \
  --standalone \
  --mathjax \
  --toc \
  --toc-depth=3 \
  --resource-path "$ROOT_DIR" \
  --bibliography "$ROOT_DIR/references.bib" \
  --citeproc \
  "${CSL_ARGS[@]}" \
  --metadata title="DOPAMINE: User and Developer Manual" \
  --output "$BUILD_DIR/dopamine_manual.html"

pandoc "$TMP_MD" \
  --from markdown+tex_math_dollars \
  --to pdf \
  --toc \
  --toc-depth=3 \
  --resource-path "$ROOT_DIR" \
  --bibliography "$ROOT_DIR/references.bib" \
  --citeproc \
  "${CSL_ARGS[@]}" \
  --pdf-engine=lualatex \
  --lua-filter "$ROOT_DIR/pdf-code-wrap.lua" \
  --include-in-header "$ROOT_DIR/pdf-header.tex" \
  -V mainfont="DejaVu Serif" \
  -V monofont="DejaVu Sans Mono" \
  -V monofontoptions=Scale=0.85 \
  --output "$BUILD_DIR/dopamine_manual.pdf"

echo "Generated:"
echo "  $BUILD_DIR/dopamine_manual.html"
echo "  $BUILD_DIR/dopamine_manual.pdf"
