# Publishing checklist — GitHub and Zenodo

Updated 2026-09-06. **Nothing has been committed or pushed.** Everything below is for you to run.

---

## The one thing to understand first

**You do not push this repository.** The working tree still carries the AlphaMissense columns,
because the local gates read them, and pushing it would publish values that are not ours to
redistribute. `scripts/66_release_bundle.py` builds a redacted copy at **`release/github/`**, and
that is what goes to GitHub.

The working tree is never modified by any of this.

```
this repository  ──build──▶  release/github/   25.7 MB   →  GitHub
                 ──build──▶  release/zenodo/  182.0 MB   →  Zenodo
```

## Your Word file is now protected

You said you are hand-editing `results/FBN1_manuscript_v4.docx`. `61_typeset_v3.py` regenerates
that file from the markdown, which would have destroyed your edits silently. It now records a
hash of what it wrote and **refuses to run** if the file has changed since. A modification time
could not detect this, because the build writes the DOCX a second after the markdown and so it is
always the newer file.

If you rebuild and see the refusal, that is the guard working. Your options are

```bash
cp results/FBN1_manuscript_v4.docx results/FBN1_manuscript_v4.docx.handedit   # keep a copy
# then either put the edits into results/manuscript_v4.md and rebuild normally,
.venv/bin/python scripts/61_typeset_v3.py
# or discard them deliberately
.venv/bin/python scripts/61_typeset_v3.py --force
```

**The markdown is the source of truth.** Edits made only in Word are lost at the next rebuild and
will not appear in the PDF, so anything you want to keep has to go back into
`results/manuscript_v4.md`.

---

## What was excluded, and why

| excluded | size | reason |
|---|---|---|
| `.claude/`, `CLAUDE.md` | — | assistant working files, not part of the scientific record |
| `am_pathogenicity`, `am_class` | 12 columns | CC BY-NC-SA, restored in one command by `scripts/67_restore_alphamissense.py` |
| `resources/papers/` | 49 MB | publisher-copyright PDFs |
| `*.a3m` alignments | 525 MB | regenerable, and no script opens one |
| `*full_data*.json` | 109 MB | per-atom confidence matrices → Zenodo |
| model `*.cif` | 48 MB | coordinates → Zenodo, where the paper points |
| `data/raw/` | 447 MB | database responses, regenerable from `manifest/` |
| `resources/databases/` | 1.2 GB | the AlphaMissense subset itself |
| `molecules/`, `tools/` | — | FoldX 5, licensed per user |
| superseded figures and reports | ~25 MB | kept on disk so the v1 and v2 gates still run |

`.env` was never staged. The audit re-checks on every run.

---

## Step 1 — build and audit

```bash
.venv/bin/python scripts/66_release_bundle.py --build
```

It refuses to build if anything forbidden is staged, and refuses again if anything forbidden or
any AlphaMissense column survives into either built tree. Then check by hand:

```bash
find release/github -path '*claude*' -o -name 'CLAUDE.md' | wc -l   # expect 0
grep -rl am_pathogenicity release/github --include='*.tsv' | wc -l  # expect 0
du -sh release/github release/zenodo                                # expect 32M and 190M
```

## Step 2 — make the public repository

```bash
cd release/github
git init -b main
git config user.name  "[Your name]"
git config user.email "[your email]"
git add -A
git commit -m "Analysis code and derived data for the FBN1 cbEGF calcium-coordination study

All 43 cbEGF domains modeled with calcium as 21 tandem constructs plus 60
mutants. Pipeline, validation gates, processed tables, figures and manuscript.

Model coordinates are in the Zenodo archive. AlphaMissense columns are not
redistributed and are restored by scripts/67_restore_alphamissense.py. Raw
database dumps, journal PDFs and licensed binaries are excluded, as recorded
in README.md and NOTICE."
```

## Step 3 — push

Create an **empty** repository on GitHub, with no README, no .gitignore and no licence, since all
three already exist in the tree. Suggested name `fbn1-cbegf-calcium`.

```bash
git remote add origin https://github.com/[USER]/[REPO].git
git push -u origin main
```

Add a description and the topics `fibrillin-1`, `marfan-syndrome`, `alphafold`,
`variant-interpretation`, `structural-bioinformatics`.

## Step 4 — Zenodo

1. At [zenodo.org](https://zenodo.org), **Settings → GitHub**, and switch the repository on.
   Zenodo then mints a DOI for every release you tag.
2. Tag and release, which produces the DOI for the code:
   ```bash
   git tag -a v1.0.0 -m "Version accompanying the manuscript"
   git push origin v1.0.0
   ```
   Publish a release from that tag on GitHub.
3. The models are too large for that route, so upload them as a **separate Zenodo deposition**:
   ```bash
   cd release && zip -r fbn1-cbegf-archive-v1.0.0.zip zenodo MANIFEST.tsv README_ARCHIVE.md
   ```
   Zenodo allows 50 GB per record, so 182 MB is comfortable.
4. On the form set **Upload type** Dataset, add the title and yourself as author with your ORCID,
   and this description:

   > Analysis code, derived data and AlphaFold Server structural models for a study of calcium
   > coordination and fold stability at pathogenic FBN1 missense positions in the fibrillin-1
   > cbEGF module. Contains the complete pipeline and validation gates, the processed variant
   > tables, the manuscript figures with their value sidecars, and 480 predicted model coordinate
   > files from 96 AlphaFold Server jobs. Code is MIT and data CC BY 4.0. AlphaMissense
   > predictions are not redistributed and are rejoined from Zenodo record 8208688 by a script
   > included here.

   Under **Related identifiers** add the GitHub URL as *is supplemented by this upload*, and
   `10.5281/zenodo.8208688` as *references*.

## Step 5 — put the identifiers back in the paper

Five placeholders remain in `results/manuscript_v4.md`.

```
[Author name]                                    ×3
[Department, Institution, City, Postcode, Country]
[email address]
[GitHub repository URL]
[DOI]
```

Edit the **markdown**, not the Word file, then rebuild and re-gate.

```bash
.venv/bin/python scripts/61_typeset_v3.py     # will refuse if your docx has hand edits
.venv/bin/python tests/test_manuscript_v4.py
```

The gate strips URLs and DOIs before applying the no-colon house rule, so a real `https://` link
will not fail it.

---

## After anything changes

Re-run the build so the two trees track the repository, and re-run the gates.

```bash
.venv/bin/python scripts/66_release_bundle.py --build
.venv/bin/python tests/test_v2.py
.venv/bin/python tests/test_manuscript_v4.py
```

If you have already pushed, commit and push again from `release/github/` — it is a normal git
repository and its history is independent of this one.
