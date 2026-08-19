# tools/

Helper scripts that are not figure panels.

| Script | Output |
|--------|--------|
| `build_supp_table_1_demographics.py` | `tables/Supplementary_Table_1_Demographics.xlsx`, the demographics of the 731-participant serology cohort by disease group, with each group tested against HD |
| `build_supp_table_il1b_genes.py` | `tables/Supplementary_Table_7_IL1B-response-genes.csv`, the IL-1B response gene list behind Figure 5C, flagged for presence in the object and for survival into the analysis universe |
| `check_reproducibility.py` | A 0-100 repository structure score, run weekly by CI |

```bash
python3 tools/build_supp_table_1_demographics.py
python3 tools/build_supp_table_il1b_genes.py
```

## check_reproducibility.py

Scores the repository on structure and hygiene, not on scientific correctness. Standard library
only, Python 3.10+.

| Category | Weight | Checks |
|----------|--------|--------|
| Step ordering | 20 | Scripts and folders carry numeric prefixes |
| Documentation | 25 | Header blocks in scripts, a README in each folder |
| Path hygiene | 20 | No hardcoded local paths |
| Data handling | 15 | `.gitignore` covers data files, READMEs say where the data live |
| Naming | 10 | No spaces or special characters in file names |
| PHI and credentials | 10 | No identifiers, passwords or API keys in code |

```bash
python3 tools/check_reproducibility.py            # current directory
python3 tools/check_reproducibility.py /path      # a specific directory
python3 tools/check_reproducibility.py --json     # write reproducibility_score.json
python3 tools/check_reproducibility.py --min-score 80   # exit 1 below the threshold
```

A score of 80 or above is good; below 60 means something structural is wrong.
