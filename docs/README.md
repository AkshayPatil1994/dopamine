<!-- SPDX-FileCopyrightText: 2026 DOPAMINE contributors -->
<!-- SPDX-License-Identifier: AGPL-3.0-only -->

# DOPAMINE Documentation

This folder contains user and method documentation for DOPAMINE. It is
published in two independent ways from the same source files:

1. **GitHub Pages**, served directly from this folder (Settings → Pages →
   Deploy from a branch → folder `/docs`) using Jekyll + the
   [just-the-docs](https://just-the-docs.com) theme (`_config.yml`).
2. **A standalone PDF/HTML manual**, built locally with Pandoc (`build.sh`).

Both outputs are generated from the same pages, so a heading added to
one of them shows up in both the Pages sidebar and the manual's table of
contents.

## Files

- `index.md`: GitHub Pages home page and task-to-page router (not in the Pandoc manual).
- `installation.md`: Dependencies, CPU/PETSc/offline/CUDA/OpenACC builds, verification, troubleshooting.
- `input-reference.md`: Every input namelist parameter with default, units and meaning.
- `preprocessing.md`: Boundary patch names, mesh formats, direct vs. pre-decomposed meshes, `dopamine-decompose`.
- `examples.md`: Overview of the example cases and the regression suite.
- `tools.md`: Executables and helper scripts.
- `references.md`: Bibliography page (Pages only).
- `theory.md`: Governing equations, time integration, advection/gradient/
  diffusion methods, SGS/wall modeling, atmospheric boundary layer (MOST,
  SEM, Davies relaxation), mesh partitioning, and pressure backend details
  with method-specific citations.
- `usage.md`: Running, restart, GPU runs, boundary condition types, outputs,
  advection naming, recommended numerics, ABL setup (MOST wall model, SEM
  inflow, nudging), and long-time stability checks.
- `development.md`: Code layout, build system internals, testing, coding
  style, how to add a scheme/BC/ABL feature, git workflow, and how to report
  bugs.
- `manual.md`: Legacy Pandoc title-page/outline stub kept for reference; the
  actual manual content is assembled directly by `build.sh` (see below), not
  by processing this file.
- `references.bib`: DOI-backed bibliography used by Pandoc citeproc.
- `build.sh`: Concatenates installation, usage, pre-processing, input
  reference, theory, examples, tools and development pages (front matter stripped) behind a single title block, then renders HTML and
  PDF with table of contents and citation processing via Pandoc.
- `pdf-code-wrap.lua`, `pdf-header.tex`: PDF-only helpers that let inline code and code blocks wrap.
- `_config.yml`, `_includes/head_custom.html`: Jekyll/just-the-docs
  configuration for the GitHub Pages site (theme, nav, search, MathJax).

## Build the standalone manual

From repository root:

```bash
bash docs/build.sh
```

Outputs:

- `docs/build/dopamine_manual.html`
- `docs/build/dopamine_manual.pdf`

Requires `pandoc`, and `lualatex` with the `fvextra` package (e.g. from a TeX
Live install) for the PDF target. The build rewrites Jekyll-only syntax
(`$$..$$` inline math, `{{ site.baseurl }}` paths) for Pandoc; the page sources
are left untouched.

## Preview the GitHub Pages site locally

```bash
cd docs
bundle exec jekyll serve
```

Requires Ruby + Bundler and a `Gemfile` with `github-pages` (or `jekyll` plus
`jekyll-remote-theme`) — not checked in, since GitHub Pages itself needs none
of this locally; add one only if you want to preview before pushing.
