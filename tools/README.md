# tools/

Utility scripts for the compbio-commons repository.

## check_reproducibility.py

Evaluates the reproducibility of a computational biology project directory.
Checks architecture, naming, documentation, and GCP data hygiene.
Produces a 0–100 score suitable for CI monitoring.

### What it checks

| Category | Weight | Description |
|----------|--------|-------------|
| Step Ordering | 20 | Scripts/folders use numeric prefixes (`01_`, `02_`) |
| Documentation | 25 | Header blocks in scripts; README in each folder |
| Path Hygiene | 20 | No hardcoded local paths; use relative or `gs://` paths |
| GCS Data Handling | 15 | `.gitignore` covers data files; READMEs mention copy/mount instructions |
| Naming Conventions | 10 | No spaces or special characters in file/folder names |
| PHI / Credential Safety | 10 | No SSNs, emails, passwords, or API keys in code |

### Usage

`check_reproducibility.py` requires Python 3.10+ and uses only the standard library — no extra packages needed.

```bash
# Evaluate current directory
python tools/check_reproducibility.py

# Evaluate a specific project directory
python tools/check_reproducibility.py /path/to/project

# Write a JSON score file (for CI / monitoring)
python tools/check_reproducibility.py --json

# Fail (exit 1) if score is below 80
python tools/check_reproducibility.py --min-score 80
```

### Score interpretation

| Score | Meaning |
|-------|---------|
| 80–100 | Good reproducibility |
| 60–79 | Needs improvement — fix top issues and re-run |
| 0–59 | Critical issues — must fix before sharing analysis |

### JSON output format

When run with `--json`, writes `reproducibility_score.json` in the target directory:

```json
{
  "score": 73,
  "timestamp": "2026-03-25",
  "path": "/path/to/project",
  "categories": {
    "step_ordering": { "score": 16, "max": 20, "issues": [] },
    "documentation": { "score": 12, "max": 25, "issues": ["..."] }
  }
}
```

This file can be used with GCP Monitoring, GitHub Actions, or any CI pipeline to track reproducibility over time.

### GCP / GCS best practices checked

- `.gitignore` must cover common data extensions (`.bam`, `.vcf`, `.rds`, `.h5`, `.csv`, `.fastq`)
- READMEs should mention how to pull data from GCS buckets (`gsutil cp`, `gcloud storage cp`) and how to copy results back
- No large data files (>10 MB) should be committed to the repo — store in `gs://` buckets
- No hardcoded local paths (`/Users/`, `/home/`, `/mnt/`) — use `gs://` URIs or relative paths
