# CheckMyMetal hand-off

Claude Code writes prepared PDBs + `manifest.csv` + `INSTRUCTIONS.md` into `inbox/` and stops.
You run CheckMyMetal (https://csgid.org/csgid/metal_sites), select calcium, and save each
result into `outbox/` using the filenames in the instruction sheet. Then tell Claude to resume.
