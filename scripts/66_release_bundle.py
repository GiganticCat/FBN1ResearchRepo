#!/usr/bin/env python3
"""
66_release_bundle.py — assemble the Zenodo archive, and check the repository is safe to publish.

Two jobs, and the second matters more than the first.

  AUDIT     refuse to build if anything that must not be published is staged for commit --
            publisher-copyright PDFs, the licensed FoldX binary or its molecules/ directory,
            the NCBI API key, or the AlphaMissense subset, which is CC BY-NC-SA and is not ours
            to redistribute. Also refuse on any staged file over GitHub's 100 MB hard limit.

  BUNDLE    write release/zenodo/ containing everything in the git tree PLUS the 480 model
            coordinate files and the AlphaFold Server confidence matrices, which are deliberately
            kept out of git and which the paper's availability statement places in the archive.

  REDACT    strip the two AlphaMissense-derived columns from the released tables. Those values
            are Cheng et al.'s under CC BY-NC-SA 4.0, and not redistributing them keeps the whole
            release cleanly CC BY 4.0. The analysis is unaffected -- the columns are rebuilt in
            one command by scripts/67_restore_alphamissense.py from the Zenodo record, and the
            summary values that appear in figures and in statistics_final.tsv stay, since four
            medians and a rank correlation are a scientific result and not a redistributed
            dataset.

Nothing here writes to the working tree. The manuscript, the tables and the figures are read and
copied, never modified, so a Word file open for hand editing is never at risk.

The split is the one the availability statement describes. GitHub carries the code and the derived
tables, roughly 26 MB. Zenodo carries the same plus the structures, roughly 185 MB. Nothing in
either is redistributed without the right to redistribute it.

**The public repository is BUILT, not pushed from here.** The working tree still carries the
AlphaMissense columns, because the local gates read them, so pushing this repository directly
would publish them. `release/github/` is the redacted tree to push. The working tree is never
modified by this script.

Writes: release/github/, release/zenodo/, release/MANIFEST.tsv, release/README_ARCHIVE.md
Run:    .venv/bin/python scripts/66_release_bundle.py [--build]
"""
from __future__ import annotations

import hashlib
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
REL = ROOT / "release"
OUT = REL / "zenodo"
GH = REL / "github"

GITHUB_FILE_LIMIT = 100 * 1024 * 1024

# Anything matching these must never reach a public repository.
FORBIDDEN = [
    (".claude/", "assistant working files, internal notes rather than scientific record"),
    ("CLAUDE.md", "the project brief, internal"),
    ("resources/papers/", "publisher-copyright journal PDFs"),
    ("resources/databases/", "the AlphaMissense subset, CC BY-NC-SA, not ours to redistribute"),
    ("molecules/", "FoldX 5 molecule data, licensed per user"),
    ("tools/", "licensed third-party binaries"),
    (".env", "the NCBI API key"),
    ("data/raw/", "raw database dumps, regenerable from manifest/"),
]

# Kept out of git, added to the Zenodo archive.
ARCHIVE_EXTRA = [
    ("structures/alphafold", "*model_*.cif", "predicted model coordinates"),
    ("structures/alphafold", "*full_data*.json", "per-atom confidence matrices"),
]

# Cheng et al.'s values, CC BY-NC-SA 4.0. Removed from every released table.
AM_COLUMNS = ["am_pathogenicity", "am_class"]


def git(*args: str) -> list[str]:
    r = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True)
    if r.returncode:
        sys.exit(f"git {' '.join(args)} failed\n{r.stderr}")
    return [x for x in r.stdout.splitlines() if x]


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def redact_tsv(src: Path, dst: Path) -> int:
    """Copy a TSV without the AlphaMissense columns. Returns how many were dropped."""
    lines = src.read_text(encoding="utf-8").splitlines()
    if not lines:
        shutil.copy2(src, dst)
        return 0
    head = lines[0].split("\t")
    drop = [i for i, c in enumerate(head) if c in AM_COLUMNS]
    if not drop:
        shutil.copy2(src, dst)
        return 0
    keep = [i for i in range(len(head)) if i not in drop]
    out = []
    for ln in lines:
        f = ln.split("\t")
        out.append("\t".join(f[i] if i < len(f) else "" for i in keep))
    dst.write_text("\n".join(out) + "\n", encoding="utf-8")
    return len(drop)


def redact_xlsx(src: Path, dst: Path) -> int:
    """Same, for the supplementary workbook, sheet by sheet."""
    from openpyxl import load_workbook
    wb = load_workbook(src)
    n = 0
    for ws in wb.worksheets:
        head = [c.value for c in ws[1]]
        for idx in sorted((i + 1 for i, c in enumerate(head) if c in AM_COLUMNS), reverse=True):
            ws.delete_cols(idx)
            n += 1
    wb.save(dst)
    return n


def audit(tracked: list[str]) -> int:
    print("AUDIT — what git is about to publish\n")
    bad = []
    for pat, why in FORBIDDEN:
        hits = [f for f in tracked if f.startswith(pat) or f"/{pat}" in f or f == pat.rstrip("/")]
        status = "clear" if not hits else f"*** {len(hits)} FILES *** {why}"
        print(f"  {pat:26} {status}")
        bad += [(f, why) for f in hits]

    big = [(f, (ROOT / f).stat().st_size) for f in tracked
           if (ROOT / f).is_file() and (ROOT / f).stat().st_size > GITHUB_FILE_LIMIT]
    print(f"  {'files over 100 MB':26} " + ("clear" if not big else f"*** {len(big)} ***"))

    total = sum((ROOT / f).stat().st_size for f in tracked if (ROOT / f).is_file())
    print(f"\n  {len(tracked)} tracked files, {total / 1048576:.1f} MB")
    if bad or big:
        for f, why in bad[:10]:
            print(f"    REFUSE {f}  ({why})")
        for f, n in big[:5]:
            print(f"    REFUSE {f}  ({n / 1048576:.0f} MB, over GitHub's hard limit)")
        return 1
    print("  nothing forbidden is staged\n")
    return 0


def copy_tree(tracked: list[str], dest: Path) -> tuple[list, int]:
    """Copy the tracked files into `dest`, stripping the AlphaMissense columns on the way."""
    if dest.exists():
        shutil.rmtree(dest)
    rows, n_redacted = [], 0
    for rel in tracked:
        src = ROOT / rel
        if not src.is_file():
            continue
        dst = dest / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        if src.suffix == ".tsv":
            n_redacted += redact_tsv(src, dst)
        elif src.suffix == ".xlsx":
            n_redacted += redact_xlsx(src, dst)
        else:
            shutil.copy2(src, dst)
        rows.append((rel, dst.stat().st_size, sha256(dst), "git"))
    return rows, n_redacted


def verify(dest: Path) -> None:
    """Nothing forbidden, and no AlphaMissense value, may survive into a built tree."""
    leaks = []
    for f in dest.rglob("*"):
        if not f.is_file():
            continue
        rel = str(f.relative_to(dest))
        for pat, why in FORBIDDEN:
            if rel.startswith(pat.rstrip("/")):
                leaks.append(f"{rel} ({why})")
        if f.suffix == ".tsv":
            head = f.open(encoding="utf-8", errors="ignore").readline()
            if any(c in head.split("\t") for c in AM_COLUMNS):
                leaks.append(f"{rel} (AlphaMissense column survived redaction)")
    if leaks:
        sys.exit("STOP: built tree is not publishable\n  " + "\n  ".join(leaks[:10]))


def build(tracked: list[str]) -> int:
    # ---- the public repository
    gh_rows, n_redacted = copy_tree(tracked, GH)
    verify(GH)
    gh_bytes = sum(n for _, n, _, _ in gh_rows)

    # ---- the archive, which is that plus the structures
    rows, n_extra = copy_tree(tracked, OUT)[0], 0

    for base, glob, what in ARCHIVE_EXTRA:
        for src in sorted((ROOT / base).rglob(glob)):
            rel = str(src.relative_to(ROOT))
            dst = OUT / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
            rows.append((rel, src.stat().st_size, sha256(src), what))
            n_extra += 1

    (REL / "MANIFEST.tsv").write_text(
        "path\tbytes\tsha256\tsource\n"
        + "\n".join(f"{r}\t{n}\t{h}\t{s}" for r, n, h, s in rows) + "\n", encoding="utf-8")

    total = sum(n for _, n, _, _ in rows)
    shutil.copy2(REL / "MANIFEST.tsv", GH / "MANIFEST.tsv") if (REL / "MANIFEST.tsv").is_file() \
        else None
    readme = f"""# Archive contents

{len(rows)} files, {total / 1048576:.0f} MB. This archive is the GitHub repository plus the
structural models, which are too large and too numerous for git.

  {len(rows) - n_extra} files  the repository, identical to the tagged GitHub release
  {n_extra} files  model coordinates and AlphaFold Server confidence matrices

`MANIFEST.tsv` gives a SHA-256 for every file. The `source` column says whether a file came from
the git tree or was added for the archive.

The model coordinates are AlphaFold Server predictions. Their terms of use are reproduced beside
them in `structures/alphafold/`. Licensing of everything else is in `NOTICE`.
"""
    (REL / "README_ARCHIVE.md").write_text(readme, encoding="utf-8")

    verify(OUT)
    print("BUILD\n")
    print(f"  release/github/  {len(gh_rows)} files, {gh_bytes / 1048576:.1f} MB"
          f"  <- push this, not the working repository")
    print(f"  release/zenodo/  {len(rows)} files, {total / 1048576:.1f} MB")
    print(f"\n  {n_redacted} AlphaMissense columns stripped from each tree")
    print(f"  {n_extra} model and confidence files added to the archive")
    print(f"  wrote release/MANIFEST.tsv and release/README_ARCHIVE.md")
    return 0


def main() -> int:
    tracked = git("ls-files")
    rc = audit(tracked)
    if rc:
        print("\nSTOP. Fix .gitignore and re-stage before publishing anything.")
        return rc
    if "--build" not in sys.argv:
        print("audit only. Re-run with --build to assemble release/zenodo/.")
        return 0
    return build(tracked)


if __name__ == "__main__":
    sys.exit(main())
