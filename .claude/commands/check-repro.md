Run the reproducibility evaluation script on the given path.

Usage: /check-repro [PATH]

Steps:
1. Find the repository root (where tools/check_reproducibility.py lives).
2. Run: `pixi run python tools/check_reproducibility.py $ARGUMENTS`
   - If $ARGUMENTS is empty, use `.` (current working directory).
   - Always run from the repository root so relative paths resolve correctly.
3. Show the full output to the user without truncation.
4. If the user asks to save the score, add `--json` to the command.

This script checks:
- Step ordering (numeric prefixes on scripts/folders)
- Documentation (header blocks, READMEs)
- Path hygiene (no hardcoded local paths)
- GCS data handling (.gitignore, copy instructions)
- Naming conventions (no spaces or special chars)
- PHI / credential safety

It does NOT fix code — it reports what to improve.
