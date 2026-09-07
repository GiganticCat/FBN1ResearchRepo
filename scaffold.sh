#!/usr/bin/env bash
# scaffold.sh — build the FBN1/Marfan project directory tree.
# Safe to re-run: uses `mkdir -p`, never overwrites existing files.
# Usage:  bash scaffold.sh            (scaffolds the current directory)
#         bash scaffold.sh <path>     (scaffolds <path>, creating it if needed)

set -euo pipefail

ROOT="${1:-$PWD}"
mkdir -p "$ROOT"
cd "$ROOT"
echo "Scaffolding project in: $(pwd)"

# --- Directory tree (mirrors the layout in CLAUDE.md) ------------------------
DIRS=(
  resources/papers
  resources/databases
  resources/reference
  data/raw
  data/interim
  data/processed
  structures/pdb
  structures/alphafold
  structures/ca_transplanted
  scripts
  tests
  figures
  results
  logs
  manifest
  handoff/cmm/inbox
  handoff/cmm/outbox
  .claude/skills
)

for d in "${DIRS[@]}"; do
  mkdir -p "$d"
  # keep otherwise-empty dirs under version control
  [ -z "$(ls -A "$d" 2>/dev/null)" ] && touch "$d/.gitkeep"
done
echo "Created ${#DIRS[@]} directories."

# --- .gitignore (only if absent) ---------------------------------------------
if [ ! -f .gitignore ]; then
cat > .gitignore <<'EOF'
# Python venv & caches
.venv/
__pycache__/
*.pyc

# Large / regenerable data — keep out of git
data/raw/
resources/databases/
structures/pdb/*.cif
structures/pdb/*.pdb
structures/alphafold/*.cif
structures/alphafold/*.pdb

# Secrets — never commit API keys
.env
*.key
secrets*

# OS / editor cruft
.DS_Store
*.swp
EOF
echo "Wrote .gitignore"
else
echo ".gitignore already exists — left untouched."
fi

# --- README stub for the handoff folder --------------------------------------
if [ ! -f handoff/cmm/README.md ]; then
cat > handoff/cmm/README.md <<'EOF'
# CheckMyMetal hand-off

Claude Code writes prepared PDBs + `manifest.csv` + `INSTRUCTIONS.md` into `inbox/` and stops.
You run CheckMyMetal (https://csgid.org/csgid/metal_sites), select calcium, and save each
result into `outbox/` using the filenames in the instruction sheet. Then tell Claude to resume.
EOF
echo "Wrote handoff/cmm/README.md"
fi

# --- Optional: initialize git ------------------------------------------------
if command -v git >/dev/null 2>&1 && [ ! -d .git ]; then
  git init -q
  echo "Initialized empty git repository."
else
  echo "Git init skipped (git missing or repo already initialized)."
fi

# --- Summary -----------------------------------------------------------------
echo
echo "Done. Directory tree:"
if command -v tree >/dev/null 2>&1; then
  tree -a -L 2 --dirsfirst -I '.git'
else
  find . -type d -not -path './.git*' | sort | sed 's|[^/]*/|  |g'
fi

echo
echo "Next: place CLAUDE.md in this root, copy your PDFs into resources/papers/,"
echo "then activate your venv and launch Claude Code from here."
