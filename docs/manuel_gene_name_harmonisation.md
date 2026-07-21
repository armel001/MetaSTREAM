# RGI BWT — Manual Gene Name Harmonization

**Addendum to:** RGI BWT Processing Memo
**Scope:** Manual cleanup of species-prefixed ARO terms, and a pipeline bug exposed and fixed as a result.

> Gnimadi TAC et al. · MetagenAMR Pipeline · 2026

---

## 1. Rationale

`arg_counts.tsv` / `arg_matrix_normalized.tsv` (RGI BWT, variant-level resolution) retained the raw CARD `aro_term` as gene identifier. A subset of ARO terms are **species-prefixed** (e.g. `Klebsiella pneumoniae KpnE`, `Escherichia coli emrE`, `Campylobacter coli chloramphenicol acetyltransferase`), which produced unreadable axis labels in the Top 50 ARG heatmap and the drug-class panel figure.

Two automated cleanup strategies were considered and rejected in favor of manual review:

| Approach | Rejected because |
|---|---|
| Edit `arg_matrix_normalized.tsv` directly (wide format) | Renaming columns by hand risks silent column collisions (e.g. three different species' `acrA` entries colliding into one column) — Excel/pandas would drop or overwrite data without warning |
| Edit `arg_counts.tsv` only | Long format avoids the collision risk, but the cleanup would only apply to one of 12 output files, creating inconsistency across the pipeline |
| Fully automated `clean_gene_label()` (regex + species-abbreviation table) | Deterministic and reusable, but the user explicitly requested manual control over naming decisions rather than a scripted heuristic |

**Decision:** review and rename gene identifiers in the pipeline's **input** file (`results/summary/rgi_bwt_all_samples.tsv`), upstream of `normalize_and_matrices.R`. Because this file is in long format (one row per sample × gene) and is the single source feeding all 12 downstream outputs, a one-time manual edit at this stage propagates consistently everywhere.

---

## 2. Workflow

### Step 1 — Extract a compact, unique-gene review table

```python
import pandas as pd

df = pd.read_csv('results/summary/rgi_bwt_all_samples.tsv', sep='\t')

unique_genes = (df[['aro_accession', 'aro_term', 'amr_gene_family', 'drug_class']]
                .drop_duplicates(subset='aro_accession')
                .sort_values('aro_term'))
unique_genes['clean_name'] = ''

output_path = 'results/summary/unique_genes_FOR_REVIEW.xlsx'
unique_genes.to_excel(output_path, index=False)
```

One row per **unique `aro_accession`** (CARD's stable identifier) rather than per raw text string, to avoid missing any entries due to text-matching inconsistencies. `amr_gene_family` and `drug_class` are included as context columns to support the manual naming decision.

### Step 2 — Manual review

The `clean_name` column was filled manually in Excel (`unique_genes_REVIEWED.xlsx`), with species prefixes stripped and ambiguous cases disambiguated at the user's discretion.

### Step 3 — Apply the mapping to the full input file

```python
import pandas as pd

df = pd.read_csv('results/summary/rgi_bwt_all_samples.tsv', sep='\t')
mapping = pd.read_excel('results/summary/unique_genes_REVIEWED.xlsx')

# Flag unfilled clean_name cells — fall back to the original term rather
# than introducing blank gene identifiers
empty_mask = mapping['clean_name'].isna() | (mapping['clean_name'].astype(str).str.strip() == '')
if empty_mask.any():
    mapping.loc[empty_mask, 'clean_name'] = mapping.loc[empty_mask, 'aro_term']

# Join by aro_accession (stable CARD ID — avoids text-matching errors)
mapping_lookup = mapping.set_index('aro_accession')['clean_name'].to_dict()

df['aro_term_original'] = df['aro_term']                      # traceability
df['aro_term'] = df['aro_accession'].map(mapping_lookup).fillna(df['aro_term'])

df.to_csv('results/summary/rgi_bwt_all_samples_CLEAN.tsv', sep='\t', index=False)
```

The original `aro_term` is preserved in a new `aro_term_original` column for full traceability — no information is destroyed, only relabeled.

### Step 4 — Substitute into the pipeline

```bash
cp results/summary/rgi_bwt_all_samples.tsv results/summary/rgi_bwt_all_samples_BACKUP.tsv
cp results/summary/rgi_bwt_all_samples_CLEAN.tsv results/summary/rgi_bwt_all_samples.tsv
bash run_AMR_pipeline.sh -t r_analysis_bwt
```

**Result:** 1,753 input rows preserved (no data loss), reduced from the original gene count to **196 unique, harmonized gene names**.

---

## 3. Bug discovered: `drug_lookup` join fan-out

### Symptom

```
[6/8] Creating matrices...
Error in `pivot_wider()`:
! Can't convert `fill` <double> to <list>.
```

### Root cause

In `normalize_and_matrices.R`, step `[2/8]`, the drug class lookup table was built as:

```r
drug_lookup <- rgi_all %>%
    select(sample_id, best_hit_aro = !!sym(col_aro), drug_class = !!sym(col_drug)) %>%
    distinct()
```

When manual renaming merges **two or more CARD entries with different original `drug_class` values** under the same clean gene name, `distinct()` retains one row per distinct `(sample_id, best_hit_aro, drug_class)` combination — producing **duplicate `(sample_id, best_hit_aro)` keys** with different `drug_class` values attached.

The subsequent `left_join(drug_lookup, by = c("sample_id", "best_hit_aro"))` then **fans out**: every row in `arg_counts` matching a non-unique key on the right side gets duplicated, silently re-introducing the per-sample-per-gene duplication that the earlier `group_by() %>% summarise()` step had just eliminated.

This duplication only manifests downstream, at step `[6/8]`, when `pivot_wider()` encounters multiple `count` values for the same `(sample_id, best_hit_aro)` cell and falls back to a list-column — which is incompatible with the numeric `values_fill = 0`, raising the type-cast error.

### Fix

Collapse multi-valued `drug_class` entries into a single semicolon-separated string **before** the join, guaranteeing a unique key on the lookup table:

```r
if (col_drug %in% colnames(rgi_all)) {
  drug_lookup <- rgi_all %>%
    select(sample_id,
           best_hit_aro = !!sym(col_aro),
           drug_class   = !!sym(col_drug)) %>%
    distinct() %>%
    group_by(sample_id, best_hit_aro) %>%
    summarise(drug_class = paste(sort(unique(drug_class)), collapse = "; "),
              .groups = "drop")

  arg_counts <- arg_counts %>%
    left_join(drug_lookup, by = c("sample_id", "best_hit_aro"))
}
```

For a gene whose merged identity now spans two original drug classes (e.g. one ancestor entry tagged `aminoglycoside`, another tagged `fluoroquinolone`), the resulting `drug_class` field becomes `"aminoglycoside; fluoroquinolone"` — consistent with how multi-valued CARD `drug_class` strings are already handled elsewhere in the pipeline (`separate_rows(sep = "; ")` in steps `[7/8]` and `[8/8]`).

### Scope of the fix

Only the `drug_lookup` join in step `[2/8]` was affected. The functional aggregation steps (`drug_class_abundance`, `mechanism_abundance`, `family_abundance`, `family_prefix_abundance`, steps `[7/8]`–`[8/8]`) build directly from `rgi_all` via `group_by()`/`summarise()` or `count()`, with no join step — they were not at risk of the same fan-out and required no change.

---

## 4. Verification

```bash
grep -A 3 "Computing raw counts" logs/r_analysis/rgi_bwt/normalize_and_matrices.log
```
```
[2/8] Computing raw counts...
  ✓ Using all_mapped_reads as abundance metric
  ✓ drug_class joined to arg_counts
  ✓ Saved: arg_counts.tsv (992 rows)
```

```bash
cut -f2 results/r_analysis/rgi_bwt/arg_counts.tsv | sort -u | wc -l
```
```
197   # 196 unique gene names + 1 header row
```

Confirms: no duplication, harmonized gene count propagated correctly through `arg_counts.tsv` and, by construction, to all 12 downstream RGI BWT output files (`arg_relative`, `arg_normalized`, `arg_presence`, the 4 wide matrices, and the 3 functional aggregation tables).

---

## 5. Files produced

| File | Location | Purpose |
|---|---|---|
| `unique_genes_FOR_REVIEW.xlsx` | `results/summary/` | Compact gene list for manual `clean_name` assignment |
| `unique_genes_REVIEWED.xlsx` | `results/summary/` | Manually completed mapping (not auto-generated — user-curated) |
| `rgi_bwt_all_samples_CLEAN.tsv` | `results/summary/` | Pipeline input with harmonized `aro_term` + `aro_term_original` for traceability |
| `rgi_bwt_all_samples_BACKUP.tsv` | `results/summary/` | Untouched copy of the original raw aggregated file |
| `normalize_and_matrices.R` | `workflow/scripts/` | Updated with the `drug_lookup` fan-out fix |

---

## 6. Reproducibility note

This manual harmonization step is **not yet integrated into the automated Snakemake DAG** — it is a one-time manual intervention on `rgi_bwt_all_samples.tsv` performed outside the pipeline's dependency graph. If the pipeline is rerun from scratch (e.g. new samples added, RGI BWT re-executed), the gene aggregation step (`aggregate_rgi_samples`) will regenerate `rgi_bwt_all_samples.tsv` from raw `gene_mapping_data.txt` files, **overwriting the harmonized version**.

**Before any pipeline rerun affecting RGI BWT samples:**
1. Back up `unique_genes_REVIEWED.xlsx` — the manual mapping decisions are the only non-reproducible artifact in this workflow
2. Re-apply Step 3 (join) to the newly regenerated `rgi_bwt_all_samples.tsv`
3. New `aro_accession` values introduced by new samples will not appear in the existing mapping table and will need a `clean_name` entry added before re-running Step 3 — the join script's defensive fallback (unfilled → original `aro_term`) prevents pipeline failure in this case, but new gene names will appear unharmonized until manually reviewed

---

## 7. Methods — Suggested wording for manuscript

> RGI BWT gene identifiers (CARD ARO terms) were manually reviewed to remove redundant species-of-origin prefixes present in a subset of CARD entries (e.g. "*Klebsiella pneumoniae* KpnE" → "KpnE"), improving readability of gene-level figures while preserving the original CARD accession-based identity for traceability. This curation was applied upstream of all downstream quantitative analyses (normalization, matrix construction, diversity, and ordination), ensuring consistency across all derived outputs.

---

*MetagenAMR Pipeline · Snakemake · Conda · Nanopore R10.4 · CARD · 2026*
