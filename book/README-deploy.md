# Deploying the bookdown to GitHub Pages

Three ways to serve the rendered guide on GitHub; pick one.

## Option A — Commit `docs/`, serve from `main` (simplest)

The bookdown is already configured to output to `../docs` (see
`book/_bookdown.yml`). To deploy:

```sh
# Build the book (lands at docs/index.html and friends)
Rscript book/render.R

# Commit the rendered HTML
git add docs
git commit -m "docs: rebuild book"
git push
```

Then on GitHub:

1. **Settings → Pages**
2. **Source**: *Deploy from a branch*
3. **Branch**: `main`, **Folder**: `/docs`
4. **Save**

The site goes live at `https://<your-username>.github.io/<repo>/`
in ~30 seconds.

A `docs/.nojekyll` file is committed alongside the HTML so GitHub
Pages skips Jekyll processing (otherwise filenames starting with `_`
would be hidden).

**Rebuild on each chapter edit**: re-run `Rscript book/render.R`, then
`git add docs && git commit && git push`.

## Option B — GitHub Actions auto-deploy (recommended for active editing)

A workflow at `.github/workflows/book.yml` builds the book on every
push to `main` that touches `book/`, `scminerViewer/R/`,
`scminer_viewer/src/`, or the workflow itself. It deploys the
rendered HTML to a `gh-pages` branch using
[`peaceiris/actions-gh-pages`](https://github.com/peaceiris/actions-gh-pages).

One-time setup on GitHub:

1. **Push** the `.github/workflows/book.yml` file to `main`. The first
   run creates the `gh-pages` branch.
2. **Settings → Pages**
3. **Source**: *Deploy from a branch*
4. **Branch**: `gh-pages`, **Folder**: `/` (root)
5. **Save**

After that, every push that changes the book auto-rebuilds and
publishes. You can stop committing `docs/` (and even add it to
`.gitignore`) if you go this route.

## Option C — Build manually and push to `gh-pages`

If you don't want a workflow and don't want to commit `docs/` to
`main`:

```sh
# Build the book — output goes to docs/ per _bookdown.yml.
Rscript book/render.R

# Push docs/ contents to the gh-pages branch (orphan history).
git worktree add /tmp/scminer-viewer-gh-pages gh-pages 2>/dev/null \
  || git worktree add --orphan -B gh-pages /tmp/scminer-viewer-gh-pages
rsync -a --delete docs/ /tmp/scminer-viewer-gh-pages/
cd /tmp/scminer-viewer-gh-pages
touch .nojekyll
git add -A
git commit -m "deploy: rebuild book"
git push origin gh-pages
cd -
git worktree remove /tmp/scminer-viewer-gh-pages
```

Then enable Pages on `gh-pages` (root) as in Option B.

## Troubleshooting

* **Pages serves the README, not the book.** You forgot to point Pages
  at `/docs` (Option A) or the `gh-pages` branch (B/C). The default
  setting serves `README.md` from the root.
* **Underscored chapter files (e.g. `_book/...`) return 404.** GitHub
  Pages ignores files starting with `_` by default; the `.nojekyll`
  file disables that behaviour. The build script and the workflow
  both create `.nojekyll` for you — make sure it's committed.
* **The book builds locally but the Action fails on
  "Render bookdown".** Usually a missing system dependency. Add
  `extra-packages: any::pandoc` or pin an older R version in the
  workflow. Check the Action logs for the exact `install.packages`
  error.
* **Pages says "404 — File not found" after enabling.** Allow 30–60
  seconds; the first build is slower. Hard-refresh (Cmd+Shift+R) once
  the green deployment notice shows on the Pages settings page.
