# REQUIREMENTS & SETUP — FBN1 / Marfan pipeline (Windows + WSL2)

Complete these top to bottom before running the Claude Code prompt. This build targets
**Windows via WSL2** (a real Ubuntu Linux on your Windows machine) — the smoothest path,
because PyMOL, DSSP, FoldX, and Claude Code all run Linux-first. Native Windows is possible but
DSSP is a genuine obstacle there, so WSL2 is the recommended default.

Legend: **[you]** = you must do it · **[claude-ok]** = Claude Code can do it once running.

---

## 1. Install WSL2 (one time)  **[you]**
- [ ] From an **admin PowerShell**: `wsl --install -d Ubuntu`, then reboot.
- [ ] On first Ubuntu launch, create a **UNIX username + password** (separate from Windows;
      this password is what `sudo` asks for).
- [ ] Confirm you're on v2 (PowerShell): `wsl -l -v` → STATE shows `VERSION 2`.
- [ ] **From here on, work in the Ubuntu terminal, not PowerShell.**

## 2. Filesystem rule (the one that matters most)  **[you]**
- [ ] Keep the whole project in the **Linux** home (`~/research/fbn1-marfan/`), **not** under
      `/mnt/c/...`. Working across the Windows boundary is slow and breaks permissions/git.
- [ ] You reach Windows files only to copy things in/out — your C: drive is mounted at `/mnt/c`.

## 3. Create the project & bring files across  **[you]**
```bash
mkdir -p ~/research/fbn1-marfan
cd ~/research/fbn1-marfan
```
- [ ] Copy the delivered files from Windows (adjust the Windows path). Note the source path and
      `*.pdf` must be joined as one argument — `.../papers/*.pdf`:
```bash
mkdir -p resources/papers
cp /mnt/c/Users/harve/Downloads/FBN1/resources/papers/*.pdf resources/papers/
cp /mnt/c/Users/harve/Downloads/FBN1/CLAUDE.md .
ls resources/papers/        # confirm the 13 PDFs landed
```
  (To copy an entire existing folder instead: `cp -r /mnt/c/Users/harve/Downloads/FBN1/. .`
  — the trailing `/.` copies contents into the current dir. `-r` is required for directories.)
- [ ] Add the ~6 paywalled papers (see the manual-download DOI/rename table) so
      `resources/papers/` matches Appendix A of `CLAUDE.md`.

## 4. Build the directory tree  **[you — run the scaffold script]**
- [ ] Copy `scaffold.sh` into the project root and run it:
```bash
bash scaffold.sh
```
  It creates `data/ structures/ scripts/ tests/ results/ figures/ logs/ manifest/ handoff/
  resources/{databases,reference} .claude/skills`, adds `.gitkeep` + a sensible `.gitignore`,
  writes the CheckMyMetal hand-off README, and inits git. Safe to re-run (never overwrites).

## 5. Python environment (venv)  **[you]**
```bash
sudo apt update && sudo apt install -y python3-venv python3-pip
python3 -m venv .venv
source .venv/bin/activate                     # prompt now shows (.venv)
pip install --upgrade pip
pip install pandas numpy scipy requests biopython biotite prody freesasa \
            matplotlib seaborn statsmodels
pip freeze > requirements.txt
```
- [ ] **Every new session** re-activate with: `cd ~/research/fbn1-marfan && source .venv/bin/activate`
      (`deactivate` to exit). Add `.venv/` to `.gitignore` (the scaffold already does).

## 6. Native structural tools  **[you]**
These are compiled programs, not pip packages: install at system level; your venv's Python
calls them.
- [ ] **PyMOL (open-source):** `sudo apt install -y pymol`
      Verify: `pymol -cq -d "print(cmd.get_version_message())"`
      Drive it headless by shelling out to `.pml` scripts (`pymol -cq figure.pml`) — reproducible
      and avoids venv/binding issues.
- [ ] **DSSP (`mkdssp`):** `sudo apt install -y dssp` then `mkdssp --version`.
      If your distro's build is old/missing, install a modern one without full Anaconda via
      **micromamba** (a single static binary):
```bash
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
./bin/micromamba create -y -p ~/tools/dssp-env -c conda-forge dssp
mkdir -p ~/.local/bin && ln -sf ~/tools/dssp-env/bin/mkdssp ~/.local/bin/mkdssp
```
      (`freesasa` in the venv covers SASA independently; DSSP adds secondary structure.)

## 7. FoldX (licensed ΔΔG — you have a license)  **[you]**
Standalone binary, venv-agnostic. Use the **Linux** build inside WSL (not `foldx.exe`).
```bash
mkdir -p ~/tools/foldx && cd ~/tools/foldx
tar -xzf foldx*.tar.gz                         # copy your downloaded tarball here first
chmod +x foldx*                                # binary may be foldx / foldx_5 / foldx5Linux64
mkdir -p ~/.local/bin && ln -sf ~/tools/foldx/foldx5Linux64 ~/.local/bin/foldx
# ensure ~/.local/bin is on PATH:
grep -q '.local/bin' ~/.bashrc || echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
foldx --version
cd ~/research/fbn1-marfan
```
- [ ] FoldX 4 needs `rotabase.txt` in the run directory; FoldX 5 bundles it.
- [ ] Optional wrapper (pure Python, installs into the venv): `pip install pyfoldx`.
- [ ] **ΔΔG workflow** (Phase 5): `RepairPDB` once, then `BuildModel` with an
      `individual_list.txt` of mutations. **Mutation strings use the PDB's author numbering, not
      UniProt** — map P35555 positions through SIFTS first (same caution as the Phase 4 gate).
      If you skip ΔΔG, tell Claude Code to omit it.

## 8. Install Node + Claude Code (in WSL)  **[you]**
```bash
curl -o- https://raw.githubusercontent.com/nvm-sh/nvm/v0.40.1/install.sh | bash
source ~/.bashrc
nvm install --lts
npm install -g @anthropic-ai/claude-code        # verify current install cmd at docs.claude.com
```
- [ ] (Optional, recommended) VS Code with the **WSL extension**: from the project run `code .`
      — editor on Windows, everything executing inside Linux (status bar shows "WSL: Ubuntu").

## 9. Database access  **[mixed]**
- [ ] **NCBI E-utilities API key** **[you]** — register at NCBI, generate a key, store it as an
      env var (e.g. in a `.env` that git ignores), **not** in a committed file.
- [ ] **AlphaMissense scores** **[you]** — pre-download the release and put the FBN1/P35555
      subset (or whole file) in `resources/databases/`; record the version. Don't rely on Claude
      fetching a multi-GB file.
- [ ] **gnomAD (v4)** and **UniProt / PDBe / SIFTS / AlphaFold DB** — public, no key. **[claude-ok]**
- [ ] Confirm outbound HTTPS to: ncbi.nlm.nih.gov, rest.uniprot.org, ebi.ac.uk,
      gnomad.broadinstitute.org, alphafill.eu. **[you to check network]**

## 10. Structure collection  **[decide; claude-ok to fetch]**
- [ ] Note the **AlphaFold model** for P35555 (large multidomain — use per-domain fragments;
      inter-domain geometry is low-confidence). AlphaFill (alphafill.eu) can supply Ca²⁺-added
      models — a shortcut for the Phase 4 transplant.
- [ ] List trusted **experimental cbEGF structures** as Ca²⁺ donors/anchors (e.g. **2W86**) in
      `resources/reference/structures.md`. Decide: Claude fetches via SIFTS, or you place PDBs
      in `structures/pdb/`.

## 11. Scientific decisions to fix up front  **[you — Claude will ask otherwise]**
- [ ] Confirm canonical transcript **NM_000138.5** (MANE) → UniProt **P35555**.
- [ ] Write inclusion thresholds to `resources/reference/inclusion_criteria.md`: ClinVar review
      stars, gnomAD (v4) allele-frequency cutoff for benign, handling of conflicts (label, don't
      drop).
- [ ] Align the **"critical residue"** definition to the ClinGen FBN1 VCEP rules (cbEGF
      cysteines, Ca-binding residues) so enrichment tests match the field standard.

## 12. The CheckMyMetal manual step  **[you]**
- [ ] Read §3b of `CLAUDE.md`. Claude Code stops after Phase 4 with prepared PDBs in
      `handoff/cmm/inbox/` + an `INSTRUCTIONS.md`. You run CheckMyMetal
      (https://csgid.org/csgid/metal_sites), select **calcium**, and save outputs to
      `handoff/cmm/outbox/` with the specified filenames, then tell Claude to resume.

## 13. Skills to author for robustness  **[you, once]**
- [ ] Confirm built-ins are available: **pdf/pdf-reading**, **xlsx**, **docx**, data-analysis/viz.
- [ ] Author project skills in `.claude/skills/` (via `skill-creator`): `variant-collection`,
      `hgvs-normalization`, `structure-mapping`, `metal-site-analysis`, `pymol-figures`,
      `data-validation`, `reproducibility`.

## 14. Reproducibility  **[you]**
`pip freeze` captures Python only; PyMOL/DSSP/FoldX live outside the venv, so record them:
```bash
{ echo "python venv: $(python --version)  ($VIRTUAL_ENV)";
  echo "pymol:  $(pymol -cq -d 'print(cmd.get_version_message())' 2>&1 | tail -1)";
  echo "mkdssp: $(mkdssp --version 2>&1 | head -1)  path=$(command -v mkdssp)";
  echo "foldx:  $(foldx --version 2>&1 | head -1)  path=$(command -v foldx)"; } \
  >> manifest/tools.md
```

---

## Then start
```bash
cd ~/research/fbn1-marfan
source .venv/bin/activate
claude
```
Send:
> Read CLAUDE.md. Start Phase 0: verify the scaffold and environment, read the PDFs in
> `resources/papers/` into `resources/reference/lit_notes.md` and cross-check them against
> Appendix A, run the env-check and its gate, and report tool versions + which databases you can
> reach. My AlphaMissense subset is in `resources/databases/` and my chosen PDB IDs are in
> `resources/reference/structures.md`. Do not start Phase 1 until I approve the Phase 0 report.

### Normal working session (quick reference)
```bash
cd ~/research/fbn1-marfan     # go to the project (Linux home, not /mnt/c)
source .venv/bin/activate     # turn on the venv  →  prompt shows (.venv)
claude                        # ...work...
deactivate                    # when done
```
