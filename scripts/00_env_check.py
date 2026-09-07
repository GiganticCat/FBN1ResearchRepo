#!/usr/bin/env python3
"""
00_env_check.py — Phase 0 environment verification for the FBN1/Marfan pipeline.

Verifies, without ever printing secrets:
  1. Directory scaffold matches CLAUDE.md section 2.
  2. Python interpreter + required package imports and versions.
  3. Native tools (pymol, mkdssp, foldx) are on PATH and callable.
  4. .env supplies NCBI_API_KEY / NCBI_EMAIL (presence only).
  5. The AlphaMissense P35555 subset is present, parseable and sane.
  6. One small live query per database, confirming response *shape*, not just HTTP 200.

Writes:
  logs/00_env_check_<UTC timestamp>.log   human-readable transcript
  logs/00_env_check_latest.json           machine-readable results (consumed by tests/test_env.py)
  manifest/tools.md                       tool versions + access dates for provenance

Exit code 0 if no check FAILed, 1 otherwise. WARN does not fail the run.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from importlib import import_module
from importlib.metadata import PackageNotFoundError, version as pkg_version
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
UTC_NOW = datetime.now(timezone.utc)
STAMP = UTC_NOW.strftime("%Y%m%dT%H%M%SZ")

# Identity facts pinned by CLAUDE.md section 1 — checks assert against these, never re-derive them.
UNIPROT_ACC = "P35555"
UNIPROT_LEN = 2871
NCBI_GENE_ID = "2200"

results: list[dict] = []
_transcript: list[str] = []


def log(line: str = "") -> None:
    print(line)
    _transcript.append(line)


def record(name: str, status: str, detail: str, **extra) -> dict:
    """status is one of PASS / WARN / FAIL."""
    entry = {"check": name, "status": status, "detail": detail, **extra}
    results.append(entry)
    log(f"  [{status:4}] {name}: {detail}")
    return entry


def section(title: str) -> None:
    log()
    log(title)
    log("-" * len(title))


def run(cmd: list[str], timeout: int = 60) -> tuple[int, str]:
    """Run a command, returning (returncode, combined stdout+stderr)."""
    try:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout, cwd=ROOT
        )
        return proc.returncode, (proc.stdout or "") + (proc.stderr or "")
    except FileNotFoundError:
        return 127, "executable not found"
    except subprocess.TimeoutExpired:
        return 124, f"timed out after {timeout}s"


# --------------------------------------------------------------------------
# 1. Directory scaffold
# --------------------------------------------------------------------------
EXPECTED_DIRS = [
    "resources/papers", "resources/databases", "resources/reference",
    "data/raw", "data/interim", "data/processed",
    "structures/pdb", "structures/alphafold", "structures/ca_transplanted",
    "scripts", "tests", "figures", "results", "logs", "manifest",
    "handoff/cmm/inbox", "handoff/cmm/outbox", "handoff/structures",
]


def check_scaffold() -> None:
    section("1. Directory scaffold")
    missing = [d for d in EXPECTED_DIRS if not (ROOT / d).is_dir()]
    created = []
    for d in missing:
        (ROOT / d).mkdir(parents=True, exist_ok=True)
        gitkeep = ROOT / d / ".gitkeep"
        if not any(p for p in (ROOT / d).iterdir()):
            gitkeep.touch()
        created.append(d)
    if created:
        record("scaffold", "WARN",
               f"{len(created)} directory(ies) were missing and have been created: "
               + ", ".join(created), created=created)
    else:
        record("scaffold", "PASS",
               f"all {len(EXPECTED_DIRS)} expected directories present", created=[])


# --------------------------------------------------------------------------
# 2. Python + packages
# --------------------------------------------------------------------------
# module name -> distribution name on PyPI (differ for several of these)
REQUIRED_PACKAGES = {
    "pandas": "pandas",
    "numpy": "numpy",
    "scipy": "scipy",
    "requests": "requests",
    "Bio": "biopython",
    "biotite": "biotite",
    "prody": "ProDy",
    "freesasa": "freesasa",
    "matplotlib": "matplotlib",
    "seaborn": "seaborn",
    "statsmodels": "statsmodels",
    "dotenv": "python-dotenv",
    "pypdf": "pypdf",  # Phase 0 literature intake (PDF text extraction)
}


def check_python() -> None:
    section("2. Python interpreter and packages")
    record("python.version", "PASS",
           f"{platform.python_version()} at {sys.executable}",
           executable=sys.executable, version=platform.python_version())

    venv = os.environ.get("VIRTUAL_ENV")
    if venv:
        record("python.venv", "PASS", venv, path=venv)
    else:
        record("python.venv", "WARN",
               "VIRTUAL_ENV is unset — is the .venv activated?")

    versions = {}
    for module, dist in REQUIRED_PACKAGES.items():
        try:
            import_module(module)
        except Exception as exc:  # noqa: BLE001 — any import failure is a gate failure
            record(f"import.{module}", "FAIL", f"import failed: {exc}")
            continue
        try:
            ver = pkg_version(dist)
        except PackageNotFoundError:
            ver = "unknown"
        versions[dist] = ver
        record(f"import.{module}", "PASS", f"{dist} {ver}", version=ver)
    return versions


# --------------------------------------------------------------------------
# 3. Native tools
# --------------------------------------------------------------------------
def check_native_tools() -> dict:
    section("3. Native structural tools")
    tools: dict[str, dict] = {}

    # --- PyMOL: get_version() (NOT the nonexistent get_version_message) -----
    path = shutil.which("pymol")
    if path is None:
        record("tool.pymol", "FAIL", "not on PATH (shutil.which('pymol') is None)")
        tools["pymol"] = {"path": None, "version": None}
    else:
        rc, out = run(["pymol", "-cq", "-d", "print(cmd.get_version()[0])"])
        ver = next(
            (ln.strip().split(">")[-1].strip() for ln in out.splitlines()
             if ln.strip() and ln.strip().split(">")[-1].strip()[:1].isdigit()),
            None,
        )
        if rc == 0 and ver:
            record("tool.pymol", "PASS", f"{ver} at {path}", path=path, version=ver)
        else:
            record("tool.pymol", "FAIL", f"callable but version not parsed (rc={rc}): {out.strip()[:200]}")
        tools["pymol"] = {"path": path, "version": ver}

    # --- DSSP ---------------------------------------------------------------
    path = shutil.which("mkdssp")
    if path is None:
        record("tool.mkdssp", "FAIL", "not on PATH")
        tools["mkdssp"] = {"path": None, "version": None}
    else:
        rc, out = run(["mkdssp", "--version"])
        ver = out.strip().splitlines()[0] if out.strip() else None
        if ver:
            record("tool.mkdssp", "PASS", f"{ver} at {path}", path=path, version=ver)
        else:
            record("tool.mkdssp", "FAIL", f"no version output (rc={rc})")
        tools["mkdssp"] = {"path": path, "version": ver}

    # --- FoldX --------------------------------------------------------------
    # FoldX 5 has no --version flag; it prints a banner then complains. The banner
    # line carrying "FoldX 5.x" is the version signal, and a non-zero exit here is
    # expected, so parse the banner rather than trusting the return code.
    path = shutil.which("foldx")
    if path is None:
        record("tool.foldx", "FAIL", "not on PATH — ask before omitting FoldX ddG (CLAUDE.md 5)")
        tools["foldx"] = {"path": None, "version": None, "molecules": None}
    else:
        rc, out = run(["foldx", "--version"])
        ver = next((ln.strip().strip("*").strip()
                    for ln in out.splitlines() if "FoldX" in ln and any(c.isdigit() for c in ln)),
                   None)
        real = os.path.realpath(path)
        # FoldX 5 replaced rotabase.txt with a molecules/ folder beside the binary.
        molecules = Path(real).parent / "molecules"
        if ver:
            record("tool.foldx", "PASS", f"{ver} at {path} -> {real}",
                   path=path, realpath=real, version=ver)
        else:
            record("tool.foldx", "FAIL",
                   f"on PATH but banner not parsed (rc={rc}): {out.strip()[:200]}")
        if molecules.is_dir():
            n = len(list(molecules.iterdir()))
            record("tool.foldx.molecules", "PASS",
                   f"{molecules} present ({n} files)", path=str(molecules))
        else:
            record("tool.foldx.molecules", "FAIL",
                   f"missing {molecules} — pass --rotabaseLocation= when running from elsewhere")
        tools["foldx"] = {"path": path, "version": ver, "molecules": str(molecules)}

    return tools


# --------------------------------------------------------------------------
# 4. Secrets (presence only — values are never read into the log)
# --------------------------------------------------------------------------
def check_env_secrets() -> None:
    section("4. Credentials (.env) — presence only, values never printed")
    env_path = ROOT / ".env"
    if not env_path.is_file():
        record("env.file", "FAIL", ".env not found at project root")
        return
    record("env.file", "PASS", ".env present")

    try:
        from dotenv import load_dotenv
        load_dotenv(env_path, override=False)
    except Exception as exc:  # noqa: BLE001
        record("env.load", "FAIL", f"python-dotenv failed to load .env: {exc}")
        return

    for key, min_len in (("NCBI_API_KEY", 20), ("NCBI_EMAIL", 5)):
        val = os.environ.get(key)
        if not val:
            record(f"env.{key}", "FAIL", "absent or empty")
        elif len(val) < min_len:
            record(f"env.{key}", "WARN", f"present but suspiciously short ({len(val)} chars)")
        else:
            record(f"env.{key}", "PASS", f"present ({len(val)} chars, value not shown)")

    if ".env" not in (ROOT / ".gitignore").read_text(encoding="utf-8"):
        record("env.gitignored", "FAIL", ".env is NOT in .gitignore — secret could be committed")
    else:
        record("env.gitignored", "PASS", ".env is git-ignored")


# --------------------------------------------------------------------------
# 5. AlphaMissense subset
# --------------------------------------------------------------------------
def check_alphamissense() -> dict:
    section("5. AlphaMissense P35555 subset")
    info: dict = {}
    subset = ROOT / "resources/databases/AlphaMissense_FBN1_P35555.tsv"
    if not subset.is_file():
        record("am.file", "FAIL", f"missing {subset}")
        return info

    sha = hashlib.sha256(subset.read_bytes()).hexdigest()
    record("am.file", "PASS",
           f"{subset.name} ({subset.stat().st_size:,} bytes, sha256 {sha[:16]}…)",
           path=str(subset), sha256=sha, bytes=subset.stat().st_size)
    info["subset_sha256"] = sha

    import pandas as pd

    # The subset was made as "header + rows matching ^P35555", which kept only the
    # copyright comment, not the real column header. Name the columns explicitly.
    df = pd.read_csv(subset, sep="\t", comment="#", header=None,
                     names=["uniprot_id", "protein_variant", "am_pathogenicity", "am_class"])
    info["rows"] = len(df)

    off_target = int((df["uniprot_id"] != UNIPROT_ACC).sum())
    if off_target:
        record("am.accession", "FAIL", f"{off_target} rows are not {UNIPROT_ACC}")
    else:
        record("am.accession", "PASS", f"all {len(df):,} rows are {UNIPROT_ACC}")

    wt = df["protein_variant"].str[0]
    pos = pd.to_numeric(df["protein_variant"].str[1:-1], errors="coerce")
    if pos.isna().any():
        record("am.variant_parse", "FAIL", f"{int(pos.isna().sum())} unparseable protein_variant values")
    else:
        record("am.variant_parse", "PASS", "all protein_variant strings parse as <WT><pos><ALT>")

    lo, hi = int(pos.min()), int(pos.max())
    if lo == 1 and hi == UNIPROT_LEN:
        record("am.positions", "PASS", f"positions span 1–{hi}, matching P35555 length")
    else:
        record("am.positions", "FAIL", f"positions span {lo}–{hi}, expected 1–{UNIPROT_LEN}")
    info["pos_min"], info["pos_max"] = lo, hi

    dups = int(df["protein_variant"].duplicated().sum())
    record("am.duplicates", "PASS" if dups == 0 else "FAIL",
           f"{dups} duplicate protein_variant strings")

    # Every position should carry exactly one WT identity. AlphaMissense sometimes
    # scores two reference alleles at a polymorphic position; that would silently
    # double-count in a Phase 2 join, so surface it here rather than downstream.
    multi = (df.assign(pos=pos, wt=wt)
               .groupby("pos")["wt"].nunique()
               .loc[lambda s: s > 1])
    if len(multi):
        detail = ", ".join(
            f"pos {p} has WT in {sorted(df.assign(pos=pos, wt=wt).query('pos == @p')['wt'].unique())}"
            for p in multi.index
        )
        record("am.wt_consistency", "WARN",
               f"{len(multi)} position(s) carry >1 wild-type residue — {detail}. "
               "Phase 2 must resolve these by the WT-identity check against the P35555 sequence.",
               positions=[int(p) for p in multi.index])
    else:
        record("am.wt_consistency", "PASS", "every position has exactly one wild-type residue")
    info["multi_wt_positions"] = [int(p) for p in multi.index]

    expected = UNIPROT_LEN * 19
    if len(df) == expected:
        record("am.row_count", "PASS", f"{len(df):,} rows = 2871 x 19 substitutions")
    else:
        record("am.row_count", "WARN",
               f"{len(df):,} rows vs {expected:,} expected (2871 x 19); "
               f"difference {len(df) - expected:+,} explained by the multi-WT position(s) above"
               if len(multi) else
               f"{len(df):,} rows vs {expected:,} expected (2871 x 19) — unexplained")

    full = ROOT / "resources/databases/AlphaMissense_aa_substitutions.tsv.gz"
    if full.is_file():
        record("am.source_archive", "PASS",
               f"full Zenodo archive present ({full.stat().st_size / 1e9:.2f} GB) "
               "— subset is reproducible from it")
    else:
        record("am.source_archive", "WARN", "full archive absent; subset cannot be re-derived locally")

    return info


# --------------------------------------------------------------------------
# 6. Database connectivity — shape-checked, not just HTTP 200
# --------------------------------------------------------------------------
def check_databases() -> dict:
    section("6. Database connectivity (small live test queries)")
    import requests

    ua = {"User-Agent": "fbn1-marfan-pipeline/0.1 (Phase 0 env check)"}
    reach: dict = {}
    timeout = 30

    def fail(name, exc):
        record(name, "FAIL", f"request failed: {type(exc).__name__}: {exc}")
        reach[name] = False

    # --- UniProt: the source of truth for sequence + domain features --------
    try:
        r = requests.get(
            f"https://rest.uniprot.org/uniprotkb/{UNIPROT_ACC}.json",
            params={"fields": "accession,sequence,protein_name"},
            headers=ua, timeout=timeout)
        r.raise_for_status()
        js = r.json()
        seq = js["sequence"]["value"]
        ok = js["primaryAccession"] == UNIPROT_ACC and len(seq) == UNIPROT_LEN
        record("db.uniprot", "PASS" if ok else "FAIL",
               f"{js['primaryAccession']} sequence length {len(seq)} "
               f"(expected {UNIPROT_LEN})", length=len(seq))
        reach["uniprot"] = ok
        # Cache the reference sequence checksum now; Phase 2's WT-identity check
        # must run against exactly this sequence.
        reach["uniprot_seq_md5"] = hashlib.md5(seq.encode()).hexdigest()
        record("db.uniprot.seq_md5", "PASS", reach["uniprot_seq_md5"])
        if reach.get("_am_multi_wt"):
            pass
        reach["_seq"] = seq
    except Exception as exc:  # noqa: BLE001
        fail("db.uniprot", exc)

    # --- NCBI E-utilities (ClinVar). API key read from env, never logged. ---
    # Use FBN1[gene]. NOT "<id>[gene_id]": that field does not exist in the ClinVar index and
    # NCBI silently degrades it to a free-text [All Fields] match, which returns records for
    # completely unrelated genes (a "2200[gene_id]" search matches COL2A1 c.2200T>C). The
    # querytranslation is checked below so a future silent re-translation cannot slip through.
    try:
        params = {"db": "clinvar", "term": "FBN1[gene]",
                  "retmode": "json", "retmax": "1"}
        key, email = os.environ.get("NCBI_API_KEY"), os.environ.get("NCBI_EMAIL")
        if key:
            params["api_key"] = key
        if email:
            params["email"] = email
        r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                         params=params, headers=ua, timeout=timeout)
        r.raise_for_status()
        js = r.json()["esearchresult"]
        count = int(js["count"])
        translation = js.get("querytranslation", "")
        if "[All Fields]" in translation:
            record("db.ncbi_clinvar", "FAIL",
                   f"query silently degraded to a free-text search ({translation!r}) — "
                   "the count would be meaningless")
            reach["ncbi_clinvar"] = False
        else:
            record("db.ncbi_clinvar", "PASS" if count > 0 else "FAIL",
                   f"esearch FBN1[gene] returned {count:,} ClinVar records "
                   f"(translated as {translation!r})", count=count)
            reach["ncbi_clinvar"] = count > 0
            reach["clinvar_count"] = count
    except Exception as exc:  # noqa: BLE001
        fail("db.ncbi_clinvar", exc)

    # --- PDBe SIFTS: UniProt <-> PDB residue mapping ------------------------
    try:
        r = requests.get(
            f"https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{UNIPROT_ACC}",
            headers=ua, timeout=timeout)
        r.raise_for_status()
        entries = r.json().get(UNIPROT_ACC, [])
        pdbs = sorted({e["pdb_id"].upper() for e in entries})
        record("db.pdbe_sifts", "PASS" if pdbs else "FAIL",
               f"{len(pdbs)} PDB entries mapped to {UNIPROT_ACC}: {', '.join(pdbs)}",
               pdb_ids=pdbs)
        reach["pdbe_sifts"] = bool(pdbs)
        reach["pdb_ids"] = pdbs
    except Exception as exc:  # noqa: BLE001
        fail("db.pdbe_sifts", exc)

    # --- AlphaFold DB -------------------------------------------------------
    # FBN1 has NO AlphaFold DB entry: AFDB omits sequences longer than ~2700 aa and
    # P35555 is 2871 aa. A bare 404 here is ambiguous between "endpoint moved" and
    # "entry genuinely absent", so probe a known-present control to tell them apart
    # and assert the absence rather than merely failing to fetch.
    AF_CONTROL = "P02751"  # fibronectin, 2477 aa — under the cap, so must be present
    try:
        r = requests.get(f"https://alphafold.ebi.ac.uk/api/prediction/{UNIPROT_ACC}",
                         headers=ua, timeout=timeout)
        if r.status_code == 200:
            rec = r.json()[0]
            record("db.alphafold", "PASS",
                   f"model {rec.get('entryId')} v{rec.get('latestVersion')} "
                   f"({rec.get('uniprotEnd')} residues)",
                   entry=rec.get("entryId"), model_version=rec.get("latestVersion"))
            reach["alphafold"] = True
            reach["af_entry"] = rec.get("entryId")
            reach["af_version"] = rec.get("latestVersion")
        elif r.status_code == 404:
            c = requests.get(f"https://alphafold.ebi.ac.uk/api/prediction/{AF_CONTROL}",
                             headers=ua, timeout=timeout)
            if c.status_code == 200:
                record("db.alphafold", "WARN",
                       f"endpoint is live (control {AF_CONTROL}, 2477 aa -> HTTP 200) but "
                       f"{UNIPROT_ACC} returns 404: FBN1 has NO AlphaFold DB model, because "
                       "AFDB omits sequences >2700 aa and FBN1 is 2871 aa. UniProt carries no "
                       "AlphaFoldDB cross-reference for P35555 either. This contradicts "
                       "CLAUDE.md 4 and removes the basis for the Phase 4 Tier-3 AlphaFill "
                       "cross-check — see results/phase0_report.md.",
                       entry=None, absent=True, control=AF_CONTROL)
                reach["alphafold"] = False
                reach["alphafold_absent_confirmed"] = True
            else:
                record("db.alphafold", "FAIL",
                       f"{UNIPROT_ACC} 404 AND control {AF_CONTROL} returned "
                       f"HTTP {c.status_code} — endpoint itself looks broken/moved")
                reach["alphafold"] = False
        else:
            record("db.alphafold", "FAIL", f"unexpected HTTP {r.status_code}")
            reach["alphafold"] = False
    except Exception as exc:  # noqa: BLE001
        fail("db.alphafold", exc)

    # --- gnomAD GraphQL -----------------------------------------------------
    try:
        query = """
        query ($geneSymbol: String!, $refGenome: ReferenceGenomeId!) {
          gene(gene_symbol: $geneSymbol, reference_genome: $refGenome) {
            gene_id
            symbol
          }
        }"""
        r = requests.post("https://gnomad.broadinstitute.org/api",
                          json={"query": query,
                                "variables": {"geneSymbol": "FBN1", "refGenome": "GRCh38"}},
                          headers={**ua, "Content-Type": "application/json"}, timeout=timeout)
        r.raise_for_status()
        js = r.json()
        gene = (js.get("data") or {}).get("gene")
        if gene and gene.get("symbol") == "FBN1":
            record("db.gnomad", "PASS",
                   f"GraphQL reachable; FBN1 = {gene.get('gene_id')}", gene_id=gene.get("gene_id"))
            reach["gnomad"] = True
        else:
            record("db.gnomad", "WARN",
                   f"reachable but unexpected payload: {json.dumps(js)[:200]}")
            reach["gnomad"] = False
    except Exception as exc:  # noqa: BLE001
        fail("db.gnomad", exc)

    # --- HGVS validation: Mutalyzer primary, VariantValidator fallback ------
    # CLAUDE.md 4 names either. VariantValidator's REST host is unreachable from
    # here, so Phase 2 leans on Mutalyzer; both are probed so the report is factual.
    try:
        r = requests.get("https://mutalyzer.nl/api/normalize/NM_000138.5:c.364C>T",
                         headers={**ua, "Accept": "application/json"}, timeout=60)
        r.raise_for_status()
        js = r.json()
        prot = (js.get("protein") or {}).get("description") if isinstance(js, dict) else None
        record("db.mutalyzer", "PASS",
               f"NM_000138.5:c.364C>T normalized"
               + (f" -> {prot}" if prot else ""),
               protein=prot)
        reach["mutalyzer"] = True
    except Exception as exc:  # noqa: BLE001
        record("db.mutalyzer", "FAIL",
               f"unreachable ({type(exc).__name__}) — Phase 2 HGVS validation needs "
               "Mutalyzer or VariantValidator; both down blocks Phase 2")
        reach["mutalyzer"] = False

    try:
        r = requests.get(
            "https://rest.variantvalidator.org/VariantValidator/variantvalidator/"
            "GRCh38/NM_000138.5:c.364C>T/NM_000138.5",
            headers={**ua, "Accept": "application/json"}, timeout=45)
        r.raise_for_status()
        js = r.json()
        hits = [k for k in js if k not in ("flag", "metadata")]
        record("db.variantvalidator", "PASS" if hits else "WARN",
               f"NM_000138.5:c.364C>T validated -> {hits[0] if hits else 'no hit'}",
               keys=hits[:3])
        reach["variantvalidator"] = bool(hits)
    except Exception as exc:  # noqa: BLE001
        record("db.variantvalidator", "WARN",
               f"unreachable ({type(exc).__name__}) — not blocking while Mutalyzer works")
        reach["variantvalidator"] = False

    # --- AlphaFill (Phase 4 Tier-3 cross-check) -----------------------------
    # AlphaFill is built ON TOP of AlphaFold DB models, so with no AFDB entry for
    # P35555 there is nothing for it to fill. Check the status *body*, not just the
    # HTTP code: a missing entry still answers 200 with status "unknown".
    try:
        r = requests.get(f"https://alphafill.eu/v1/aff/AF-{UNIPROT_ACC}-F1/status",
                         headers=ua, timeout=timeout)
        status = (r.json() or {}).get("status") if r.status_code == 200 else None
        if status == "finished":
            record("db.alphafill", "PASS", f"AF-{UNIPROT_ACC}-F1 available at alphafill.eu")
            reach["alphafill"] = True
        else:
            record("db.alphafill", "WARN",
                   f"service reachable (HTTP {r.status_code}) but AF-{UNIPROT_ACC}-F1 "
                   f"status is '{status}' — no AlphaFill entry for FBN1, as expected given "
                   "the absent AlphaFold DB model. Phase 4 Tier-3 cross-check is unavailable "
                   "for FBN1 unless run on a per-domain model.",
                   alphafill_status=status)
            reach["alphafill"] = False
    except Exception as exc:  # noqa: BLE001
        record("db.alphafill", "WARN",
               f"unreachable ({type(exc).__name__}) — Phase 4 cross-check only, not blocking")
        reach["alphafill"] = False

    return reach


# --------------------------------------------------------------------------
# 7. Literature intake inventory (content review happens separately)
# --------------------------------------------------------------------------
def check_papers() -> dict:
    section("7. Literature folder inventory")
    papers = sorted((ROOT / "resources/papers").glob("*.pdf"))
    record("papers.present", "PASS" if papers else "FAIL",
           f"{len(papers)} PDF(s) in resources/papers/")
    return {"pdfs": [p.name for p in papers]}


# --------------------------------------------------------------------------
# Provenance output
# --------------------------------------------------------------------------
def write_tools_manifest(tools: dict, versions: dict, reach: dict) -> None:
    lines = [
        "# manifest/tools.md — non-pip tool provenance",
        "",
        f"Generated by `scripts/00_env_check.py` at {UTC_NOW.isoformat()}.",
        "Regenerate by re-running that script; do not hand-edit.",
        "",
        "## Platform",
        "",
        f"- OS: {platform.platform()}",
        f"- Python: {platform.python_version()} ({sys.executable})",
        f"- venv: {os.environ.get('VIRTUAL_ENV', 'not activated')}",
        "",
        "## Native tools (not captured by requirements.txt)",
        "",
        "| tool | version | path |",
        "|------|---------|------|",
    ]
    for name, info in tools.items():
        lines.append(f"| {name} | {info.get('version') or 'MISSING'} | {info.get('path') or '—'} |")
    lines += [
        "",
        "## Key Python packages",
        "",
        "| package | version |",
        "|---------|---------|",
    ]
    for dist, ver in sorted(versions.items()):
        lines.append(f"| {dist} | {ver} |")
    lines += [
        "",
        "## Database access (this run)",
        "",
        f"Access date: {UTC_NOW.date().isoformat()}",
        "",
        "| source | reachable | note |",
        "|--------|-----------|------|",
        f"| UniProt {UNIPROT_ACC} | {reach.get('uniprot', False)} | sequence md5 `{reach.get('uniprot_seq_md5', 'n/a')}` |",
        f"| ClinVar (E-utilities) | {reach.get('ncbi_clinvar', False)} | {reach.get('clinvar_count', 'n/a')} records for Gene ID {NCBI_GENE_ID} |",
        f"| PDBe SIFTS | {reach.get('pdbe_sifts', False)} | {len(reach.get('pdb_ids', []))} PDB entries |",
        f"| AlphaFold DB | {reach.get('alphafold', False)} | "
        + ("NO model for P35555 (>2700 aa cap); endpoint verified live via control P02751"
           if reach.get("alphafold_absent_confirmed")
           else f"{reach.get('af_entry', 'n/a')} v{reach.get('af_version', 'n/a')}") + " |",
        f"| gnomAD GraphQL | {reach.get('gnomad', False)} | v4 endpoint |",
        f"| Mutalyzer | {reach.get('mutalyzer', False)} | primary HGVS validator (NM_000138.5) |",
        f"| VariantValidator | {reach.get('variantvalidator', False)} | REST host unreachable; fallback |",
        f"| AlphaFill | {reach.get('alphafill', False)} | no entry for FBN1 (depends on AFDB model) |",
        "",
    ]
    (ROOT / "manifest/tools.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    log(f"\nWrote manifest/tools.md")


def main() -> int:
    log(f"FBN1/Marfan Phase 0 environment check")
    log(f"UTC: {UTC_NOW.isoformat()}")
    log(f"Root: {ROOT}")

    check_scaffold()
    versions = check_python() or {}
    tools = check_native_tools()
    check_env_secrets()
    am = check_alphamissense()
    reach = check_databases()
    papers = check_papers()

    # The AlphaMissense multi-WT anomaly can only be adjudicated against the real
    # sequence, so resolve it here now that UniProt has been fetched.
    seq = reach.pop("_seq", None)
    if seq and am.get("multi_wt_positions"):
        section("8. Cross-check: AlphaMissense multi-WT positions vs P35555 sequence")
        for p in am["multi_wt_positions"]:
            record(f"am.wt_vs_uniprot.{p}", "WARN",
                   f"P35555 residue {p} is {seq[p - 1]} — Phase 2 must keep only "
                   f"the AlphaMissense rows whose WT is {seq[p - 1]}",
                   position=p, uniprot_residue=seq[p - 1])

    n_fail = sum(1 for r in results if r["status"] == "FAIL")
    n_warn = sum(1 for r in results if r["status"] == "WARN")
    n_pass = sum(1 for r in results if r["status"] == "PASS")

    section("Summary")
    log(f"  PASS {n_pass}   WARN {n_warn}   FAIL {n_fail}")
    for r in results:
        if r["status"] != "PASS":
            log(f"    [{r['status']}] {r['check']}: {r['detail']}")

    write_tools_manifest(tools, versions, reach)

    payload = {
        "timestamp_utc": UTC_NOW.isoformat(),
        "root": str(ROOT),
        "summary": {"pass": n_pass, "warn": n_warn, "fail": n_fail},
        "results": results,
        "tools": tools,
        "packages": versions,
        "alphamissense": am,
        "databases": reach,
        "papers": papers,
    }
    (ROOT / "logs/00_env_check_latest.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8")
    logfile = ROOT / f"logs/00_env_check_{STAMP}.log"
    logfile.write_text("\n".join(_transcript) + "\n", encoding="utf-8")
    log(f"Wrote {logfile.relative_to(ROOT)} and logs/00_env_check_latest.json")

    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
