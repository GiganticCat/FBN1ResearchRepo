# Publishing — every step, every field

Updated 2026-09-06. Follow top to bottom. Every value you need to type is written out. The only thing still
unknown is the **dataset DOI**, which step 5 creates.

---

## Facts you will be asked for, in one place

| field | value |
|---|---|
| Paper title | Distinct mechanisms of calcium-ligand and cysteine substitutions in the fibrillin-1 cbEGF module |
| GitHub repository | `https://github.com/GiganticCat/FBN1ResearchRepo` |
| ORCID iD | `0009-0004-3974-571X` |
| Email | `harveyzhang612@gmail.com` |
| Author name | `Harvey Zhang` |
| Affiliation | `Independent Researcher` (PLOS will also want a country) |
| Archive file | `~/research/fbn1-marfan/release/fbn1-cbegf-archive-v1.0.0.zip` |
| Archive size | 43.8 MB, 2,559 files |
| Code licence | MIT |
| Data licence | CC BY 4.0 |
| Model licence | AlphaFold Server Output Terms — **non-commercial only**, not CC BY |

---

## Step 1 — push the licensing fixes

The repo you pushed predates the AlphaFold Server terms notice. Send it before publishing the
release, so the snapshot Zenodo archives is the correct one.

```bash
cd ~/research/fbn1-marfan/release/github
git add -A
git commit -m "Add AlphaFold Server Output Terms notice and the licensing split"
git push
git tag -f v1.0.0
git push -f origin v1.0.0
```

Force-moving the tag is safe **only because the release is not published yet**. If you have
already published it, stop and use `v1.0.1` instead — never rewrite a published release.

## Step 2 — connect Zenodo to GitHub, before tagging anything else

1. https://zenodo.org → **Log in with GitHub**
2. Top-right menu → **Settings** → **GitHub**
3. Find `GiganticCat/FBN1ResearchRepo` and switch the toggle **ON**

Zenodo only sees releases published *after* this switch is on. Getting the order wrong is the one
mistake that costs you a version number.

## Step 3 — publish the GitHub release

https://github.com/GiganticCat/FBN1ResearchRepo/releases/new

| field | what to enter |
|---|---|
| Choose a tag | `v1.0.0` (it already exists — pick it, do not create a new one) |
| Release title | `v1.0.0` |
| Describe this release | `Version accompanying the manuscript "Distinct mechanisms of calcium-ligand and cysteine substitutions in the fibrillin-1 cbEGF module". Analysis pipeline, validation gates, processed variant tables, manuscript figures and the manuscript. Structural models are in the companion Zenodo dataset.` |
| Set as pre-release | leave **unticked** |
| Set as latest release | leave **ticked** |

Publish. Within a minute or two a DOI appears on your Zenodo GitHub settings page. **That is the
code DOI.** You upload nothing for it — Zenodo takes GitHub's own source archive, which contains
no models, which is exactly right.

## Step 4 — get the archive out of WSL

```bash
cp ~/research/fbn1-marfan/release/fbn1-cbegf-archive-v1.0.0.zip /mnt/c/Users/$USER/Downloads/
```

If your Windows username differs from your Linux one that will fail — use Explorer instead and
paste this into the address bar:

```
\\wsl.localhost\Ubuntu\home\harvey\research\fbn1-marfan\release
```

## Step 5 — the dataset upload

https://zenodo.org/uploads/new — drag the zip in, then fill exactly this.

**Basic information**

| field | value |
|---|---|
| Resource type | **Dataset** |
| Title | `Analysis code, data and AlphaFold Server structural models for calcium coordination and fold stability of pathogenic FBN1 missense variants` |
| Publication date | today |
| Creators → Family name / Given names | `Zhang` / `Harvey` |
| Creators → Identifier | `0009-0004-3974-571X` |
| Creators → Affiliation | `Independent Researcher` |
| Contributors | leave empty |

There is **no Role field for Creators**. If you find yourself looking at one you are in the
Contributors box by mistake — go back up. Contributors is optional and you do not need it.

**Description** — paste verbatim

> Analysis code, derived data and AlphaFold Server structural models for a study of calcium
> coordination and fold stability at pathogenic FBN1 missense positions in the fibrillin-1 cbEGF
> module. All 43 calcium-binding EGF-like domains were modeled with their calcium ions as 21
> overlapping tandem constructs together with 60 mutants, and calcium ligands were identified
> from observed contact with the modeled ion instead of from the sequence consensus.
>
> Contains the complete numbered pipeline and its validation gates, the processed variant tables,
> the manuscript figures with sidecars listing every plotted value, and 480 predicted model
> coordinate files from 96 AlphaFold Server jobs.
>
> Licensing is mixed and is set out in full in NOTICE. Code is MIT and the derived data tables
> are CC BY 4.0. The 480 model coordinate files, their confidence matrices and the figures
> rendered from them are AlphaFold Server Output, subject to the AlphaFold Server Output Terms of
> Use at alphafoldserver.com/output-terms, and are for non-commercial use only. AlphaMissense
> predictions are not redistributed; scripts/67_restore_alphamissense.py rejoins them from Zenodo
> record 8208688.

**Licence** — do **not** choose CC BY 4.0 here. The models forbid commercial use and CC BY
permits it, so CC BY would grant rights you do not hold.

| field | value |
|---|---|
| Licence | **Other (Non-Commercial)**, or if the form asks for a custom licence, name it `AlphaFold Server Output Terms of Use` with URL `https://alphafoldserver.com/output-terms` |

**Related works** — add two

| relation | identifier | type |
|---|---|---|
| *is supplemented by* | `https://github.com/GiganticCat/FBN1ResearchRepo` | URL |
| *references* | `10.5281/zenodo.8208688` | DOI |

Then **Publish**. You now have a second DOI — **the dataset DOI**, which is the one the paper
cites.

## Step 6 — put the identifiers into the paper

One placeholder remains. Edit **`results/manuscript_v4.md`**, never the Word file.

| placeholder | replace with |
|---|---|
| `[DOI]` | the **dataset** DOI from step 5 |

Everything else is already filled in — name, affiliation, ORCID, email and the repository URL
went in on 2026-09-06. `[DOI]` is the only placeholder left in the manuscript.

Then rebuild and re-gate:

```bash
cd ~/research/fbn1-marfan
.venv/bin/python scripts/61_typeset_v3.py
.venv/bin/python tests/test_manuscript_v4.py
```

If the typesetter refuses, that is the guard protecting hand edits you made in Word. Copy the file
aside, move the edits into the markdown, then rerun. `--force` discards them.

## Step 7 — push the finished manuscript

```bash
cd ~/research/fbn1-marfan
.venv/bin/python scripts/66_release_bundle.py --build
cd release/github && git add -A && git commit -m "Add author details and archive DOI" && git push
```

No new tag needed unless you want a citable v1.0.1.

---

## PLOS ONE submission

| field | value |
|---|---|
| Corresponding author ORCID | `0009-0004-3974-571X` — **mandatory**, link the account when prompted |
| Data availability | choose *All data are in the manuscript and/or supporting information*, then paste the Data and code availability section from the manuscript |
| Financial disclosure | `The author received no specific funding for this work.` |
| Competing interests | `The authors have declared that no competing interests exist.` |
| Ethics statement | `This study used only variant records and allele frequencies already in the public domain. No identifiable individual-level data were accessed and no human participants were recruited, so the work required neither ethical approval nor informed consent.` |
| CRediT roles | Conceptualization, Methodology, Software, Formal analysis, Investigation, Data curation, Visualization, Writing – original draft |

Leave Funding acquisition, Resources and Supervision unticked unless someone else supplied them.
If your mentor shaped the study design or interpretation, that is normally Supervision and
Conceptualization and an authorship conversation, better settled before submission.

---

## Rebuilding, any time after this

```bash
cd ~/research/fbn1-marfan
.venv/bin/python scripts/66_release_bundle.py --build   # audits, then rebuilds both trees
.venv/bin/python tests/test_v2.py
.venv/bin/python tests/test_manuscript_v4.py
cd release/github && git add -A && git commit -m "..." && git push
```

`release/github/` keeps its own git history and its own remote across rebuilds. The bundler used
to delete the directory outright, which destroyed that history once; it now clears the contents
and leaves `.git` alone.

To regenerate the archive zip after a rebuild:

```bash
cd ~/research/fbn1-marfan/release && ~/research/fbn1-marfan/.venv/bin/python - <<'EOF'
import zipfile, hashlib
from pathlib import Path
out = Path("fbn1-cbegf-archive-v1.0.0.zip")
out.unlink(missing_ok=True)
items = [Path("MANIFEST.tsv"), Path("README_ARCHIVE.md")] + sorted(Path("zenodo").rglob("*"))
with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    for p in items:
        if p.is_file():
            z.write(p, p.as_posix())
h = hashlib.sha256(out.read_bytes()).hexdigest()
print(f"{out} {out.stat().st_size/1048576:.1f} MB\nsha256 {h}")
EOF
```

`zip` is not installed on this machine, which is why this uses Python.
