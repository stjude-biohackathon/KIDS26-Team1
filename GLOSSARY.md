# SCRAP-AI — Parameter & Output Glossary

Every sidebar control and every output table column shown in the UI, in one
place. Organized by where you encounter it: **Sidebar parameters** (drive the
4-agent pipeline and are cached together as one unit) → **Output columns**, one
section per tab, in the order tabs appear.

Sidebar controls not listed here (Theme, Seed genes text area, Run analysis
button) are UI-only and don't feed the ranking math.

---

## 1. Sidebar parameters

### Cohort

| Parameter | Type / range | Meaning | Example |
|---|---|---|---|
| **Age group** | Pediatric / Adult | Restricts the cell-line cohort by DepMap `AgeCategory`. | `Pediatric` → only pediatric-derived lines. |
| **Include Fetus** | checkbox (Pediatric only) | Adds fetal-derived lines to the pediatric cohort. | Off by default — fetal lines are rare and can skew small cohorts. |
| **Cancer type (OncotreePrimaryDisease)** | dropdown, `n=` shown | The disease used to select the cohort; `n` is how many cell lines in that age×disease bucket have any DepMap data. | `Neuroblastoma (n=50)`. |

### Seeds

| Parameter | Type | Meaning | Example |
|---|---|---|---|
| **Seed source** | INO80 / SRCAP / Both / Custom | Which chromatin-remodeling complex defines the "seed set" everything else is measured against. `Both` de-duplicates shared subunits (RUVBL1, RUVBL2, ACTL6A, ACTB appear once). `Custom` empties the seed list below and reveals a `.txt` file uploader — paste your own genes or upload a file (newline/comma/whitespace separated) instead of using the built-in INO80/SRCAP proteins. | `Both` → 23 seed genes. |

### Method & population (collapsed expander)

| Parameter | Type | Meaning | Example |
|---|---|---|---|
| **Correlation method** | pearson / spearman | Linear (Pearson) vs rank-based (Spearman) correlation between each candidate gene and each seed gene. | `pearson` — the default, sensitive to outliers; use `spearman` if a few extreme cell lines look like they're driving a hit. |
| **Correlation population** | Cohort / Pan-cancer | Which cell lines are correlated to find co-expression **candidates**. `Cohort` = only the selected age×disease lines. `Pan-cancer` = every analyzable line regardless of cohort (candidate discovery only — essentiality metrics are always cohort-scoped). | `Cohort (age + disease)`. |

### Co-expression thresholds (collapsed expander)

| Parameter | Type / range | Meaning | Example |
|---|---|---|---|
| **min \|r\|** | slider, 0.0–0.9, default **0.30** | Minimum aggregate correlation (`mean_r` if positive-only, else `mean_abs_r`) to the seed set for a candidate to survive into ranking. Cohort-level correlations rarely exceed ~0.3, so this threshold is already fairly strict. | `0.30` → only genes whose average correlation to the seed set is ≥0.30. |
| **max FDR q** | select-slider {0.001, 0.01, 0.05, 0.1, 0.25}, default **0.05** | Benjamini-Hochberg-adjusted p-value ceiling (`min_q`) — controls the false-discovery rate across all candidate×seed tests. | `0.05` → at most 5% of surviving hits expected to be false positives. |
| **min cohort size** | number input, 3–100, default 15 | Minimum number of cell lines required before correlation is even attempted; below this, results are too noisy to trust. | `15`. |
| **variance filter percentile** | slider, 0–75, default 25 | Drops the lowest-variance candidate genes (bottom Nth percentile of expression variance across the cohort) before correlating — removes near-constant genes that can't meaningfully correlate with anything. | `25` → bottom quartile by variance excluded. |
| **positive correlations only** | checkbox, default **on** | If checked, ranks by signed `mean_r` (only positive co-expression counts); if unchecked, ranks by `mean_abs_r` (anti-correlated genes also qualify). | On → a gene that's strongly *anti*-correlated with the seed set is excluded even if `\|r\|` is large. |

### Expression & essentiality (expanded by default)

| Parameter | Type / range | Meaning | Example |
|---|---|---|---|
| **expressed if TPM >** | slider, 0–10, default 3.0 | TPM cutoff (applied to the log2(TPM+1) matrix) for calling a gene "expressed" in a given cell line. | `3.0` → a gene must exceed 3 TPM in a line to count as expressed there. |
| **essential if Chronos ≤** | slider, −2.0–0.0, default −0.5 | Chronos gene-effect cutoff for calling a gene "essential" in a line. More negative Chronos = the cell line depends more on that gene for survival. | `−0.5` — DepMap's conventional dependency threshold. |
| **seed set 'on' if ≥ this frac of seeds expressed** | slider, 0–1, default 0.5 | Fraction of the seed genes that must be individually "expressed" (per the TPM cutoff above) in a cell line for the whole seed complex to be considered active ("on") in that line. | `0.5` → at least half the 23 seed genes must be expressed in a line for INO80/SRCAP to count as "on" there. |
| **common-essential penalty (ranking)** | slider, 0–1, default 0.25 | Multiplier applied to `selectivity_score` for genes on DepMap's pan-essential list. `0` effectively excludes them from the top of the ranking; `1` applies no penalty. | `0.25` → a common-essential gene's score is cut to a quarter of what it would otherwise be. |
| **Pan-cancer baseline excludes the cohort** | checkbox, default off | Whether the pan-cancer essentiality baseline (used to compute `selectivity_delta`) excludes the selected cohort's own cell lines. On = purer "in-cohort vs strictly-outside-cohort" contrast. | Off → pan-cancer baseline includes the cohort lines too. |

### Show/hide

| Parameter | Type | Meaning |
|---|---|---|
| **Show common-essential genes** | checkbox, default **off** | Hides DepMap common-essential (pan-essential) genes from tables, plots, and the heatmap everywhere in the app — and from the auto-populated gene lists on the Tissue Specificity / Drug Targets / Insights tabs. Does not change the underlying ranking math, only what's displayed. |

### STRING interactions

| Parameter | Type / range | Meaning | Example |
|---|---|---|---|
| **Fold STRING into ranking** | checkbox, default **on** | Whether to fetch STRING physical/functional association scores for the top-ranked genes and use them to re-rank (one API call per gene). | On → the Interactions tab computes eagerly on page load. |
| **top genes to re-rank** | slider, 10–50, default 25 | How many top-ranked genes get STRING-scored (cost scales with this). | `25`. |
| **selectivity weight** | slider, 0–1, default 0.5 | Weight on the DepMap `selectivity_score` (normalized 0–1) in the combined `string_score`. | — |
| **physical weight** | slider, 0–1, default 0.35 | Weight on STRING physical-interaction evidence. | — |
| **functional weight** | slider, 0–1, default 0.15 | Weight on STRING functional-association evidence (co-expression, pathway membership, etc. as scored by STRING itself). | — |
| **edge cutoff** | slider, 0.15–0.9, default 0.4 | STRING association scores below this count as **0** in the scoring (not just hidden from the network graph) — a hard floor on how weak a link can be and still contribute. | `0.4`. |

*(The three weights above auto-normalize to sum to 1, so `string_score` always stays in [0,1] regardless of the raw slider values.)*

### Tissue Specificity tab (in-tab, not sidebar)

| Parameter | Type / range | Meaning | Example |
|---|---|---|---|
| **Gene(s)** | text input | Genes to run tissue-specificity triage on; auto-populated with the current top-10 ranked genes, editable. | `BRD7, DNAJC11, ZC3H13, ...` |
| **Normal GTEx tissue(s)** | multiselect | GTEx normal tissues the gene must be "silent" in. Auto-suggested from the cohort's lineage. | `Adrenal_Gland, Brain_Amygdala, ...` |
| **Matched tumor cohort (TCGA/TARGET)** | multiselect | The tumor type(s) the gene should be *high* in — the "target" cohort. Auto-suggested from the selected cancer type. | `Neuroblastoma`. |
| **normal/pan-cancer TPM cutoff** | slider, 1–20, default 5.0 | TPM threshold used both for the normal-tissue gate and the pan-cancer gate. | `5.0`. |
| **pan-cancer low fraction** | slider, 0.0–1.0, default **0.20** | Minimum fraction of *other* TCGA/TARGET cancer types the gene must be below-cutoff in to pass the pan-cancer gate. Also multiplies directly into `selectivity_score`. | `0.20` → must be quiet in ≥20% of other cancer types. |
| **target median TPM (positive gate)** | slider, 1–50, default 10.0 | Minimum median TPM in the target cohort to pass on the "high expression" side, independent of the fraction-expressed gate below. | `10.0`. |
| **target fraction expressed (positive gate)** | slider, 0.1–1.0, default 0.50 | Alternative way to pass the positive gate: at least this fraction of target-cohort samples must exceed the TPM cutoff. | `0.50`. |
| **Normal-tissue aggregation** | median / all | `median`: tolerant of one atypical normal-tissue subregion. `all`: every selected normal tissue must individually pass. | `median` — matches the original collaborator pipeline. |

### Drug Targets tab (in-tab, not sidebar)

| Parameter | Type / range | Meaning | Example |
|---|---|---|---|
| **Gene(s)** | text input | Genes to score for druggability; auto-populated with the current top-10 ranked genes. | — |
| **max genes to score** | slider, 1–50, default 15 | Caps how many of the entered genes are actually scored (cost control — each gene queries ChEMBL/DGIdb/Open Targets/MyGene.info). | `15`. |
| **Cancer-specific 5-tier scoring** | checkbox, default **on** | Switches from the simple 1–3 scoring to the cancer-type-aware 1–5 scoring (see `Score` below). | On. |
| **Cancer type (Open Targets-resolved)** | text input | The cancer type used to check cancer-specific drug indications; defaults to the sidebar's selected disease. | `Neuroblastoma`. |
| **min score to include in main results** | slider | Filters the results table to only rows at or above this `Score` tier. | — |
| **max drugs per gene shown in heatmap** | slider, 5–50, default 15 | Caps how many drugs per gene appear in the heatmap (readability, not a data filter). | `15`. |
| **Enable PRISM/DepMap sensitivity enrichment** | checkbox, default off | Adds an extra ~100 MB download + ~5 min to enrich results with PRISM drug-sensitivity screening data (per-lineage AUC). | Off by default — opt-in due to cost. |

---

## 2. Output columns — Overview tab

| Column / metric | Meaning | Example |
|---|---|---|
| **cohort lines (both data)** | Number of cell lines in the selected age×disease cohort that have *both* expression and CRISPR (Chronos) data — the actual denominator for every essentiality percentage downstream. | `33`. |
| **co-expr candidates** | Number of genes that passed the co-expression thresholds (before essentiality ranking). | `19` (or fewer once common-essential genes are hidden). |
| **ranked genes (shown)** | Number of genes actually displayed after essentiality ranking, respecting the "Show common-essential genes" toggle. Label changes to plain "ranked genes" if that toggle is on (nothing hidden). | `14` (of 19, with 5 common-essential hidden). |
| **lines correlated** | Number of cell lines used in the correlation step itself (may differ slightly from cohort size depending on the correlation population setting). | `33`. |

---

## 3. Output columns — Co-expression tab

One row per candidate gene (rows shown depend on the seed genes and correlation population selected).

| Column | Meaning | Example |
|---|---|---|
| **mean_r** | Average signed Pearson/Spearman correlation between this candidate gene and all seed genes it was tested against. | `0.34` — positively co-expressed with the seed set on average. |
| **mean_abs_r** | Average of the *absolute* correlation across seeds — used for ranking when "positive correlations only" is unchecked. | `0.41`. |
| **max_abs_r** | The single strongest \|r\| this candidate has to any one seed gene. | `0.62`. |
| **best_seed** | Which seed gene the candidate is most strongly correlated with. | `INO80B`. |
| **best_r** | The signed correlation value with `best_seed`. | `0.58`. |
| **min_q** | The smallest (most significant) BH-adjusted p-value across all seed comparisons for this gene. | `0.012`. |
| **n_seeds_signif** | How many seed genes this candidate is significantly correlated with (q < the max FDR q threshold). | `4`. |
| **n_seeds** | How many seed genes had non-missing correlation values for this candidate (denominator context for `n_seeds_signif`). | `21`. |

---

## 4. Output columns — Essentiality & Ranking tab

Builds on the co-expression columns above, adding per-cell-line essentiality metrics and the final ranking score.

| Column | Meaning | Example |
|---|---|---|
| **pct_coexpressed** | % of cohort cell lines where the candidate gene AND the seed set are both "expressed"/"on". | `72.7` → co-expressed in ~73% of the cohort's lines. |
| **pct_essential** | % of cohort cell lines where the gene is essential (Chronos ≤ threshold), regardless of co-expression. | `45.5`. |
| **pct_essential_pan** | % of the pan-cancer baseline (all analyzable lines, optionally excluding the cohort) where the gene is essential — the "how essential is this gene *everywhere else*" comparison. | `18.2`. |
| **selectivity_delta** | `pct_essential − pct_essential_pan`. Positive = more essential in this cohort than pan-cancer (cohort-selective); negative = more essential elsewhere (anti-selective, e.g. a pan-cancer housekeeping dependency). | `+27.3` → notably more essential in this cohort than broadly. |
| **pct_joint** | % of cohort lines where the gene is BOTH co-expressed with the seed set AND essential — the primary raw ranking signal before the selectivity weighting. | `38.2`. |
| **pct_essential_given_coexpr** | Conditional: of the lines where the gene *is* co-expressed, what fraction is it also essential in? | `52.5` — given co-expression, essential over half the time. |
| **n_coexpressed / n_essential / n_joint** | Raw cell-line counts behind the three percentages above. | `24 / 15 / 13`. |
| **n_lines** | Cohort denominator (same as the "cohort lines (both data)" Overview metric). | `33`. |
| **common_essential** | Whether the gene is on DepMap's pan-essential gene list (independent of this cohort). | `True` / `False`. |
| **selectivity_score** | The final ranking metric: `pct_joint × (1 + selectivity_delta/100) × ce_penalty_multiplier`. Higher = more cohort-selectively essential and co-expressed with the seed complex, penalized if pan-essential. | `29.4`. |

---

## 5. Output columns — Interactions (STRING) tab

Adds STRING evidence on top of the ranked table above (only for the top-N genes re-ranked).

| Column | Meaning | Example |
|---|---|---|
| **sel_norm** | The candidate's `selectivity_score`, min-max normalized to [0,1] across the top-N genes being re-ranked. | `0.83`. |
| **string_phys** | Maximum STRING **physical**-interaction confidence score (0–1) to any one seed gene — evidence of direct physical binding/complex membership. | `0.91` → strong physical-interaction evidence with at least one seed protein. |
| **string_func** | Maximum STRING **functional**-association confidence score (0–1) to any one seed gene — broader evidence (co-expression, pathway, text-mining, etc., as STRING itself defines it). | `0.67`. |
| **string_best_phys_seed** | Which seed gene gives the strongest physical-association score. | `RUVBL1`. |
| **string_best_func_seed** | Which seed gene gives the strongest functional-association score. | `ACTR8`. |
| **string_score** | Combined score: `w_sel·sel_norm + w_phys·string_phys(cutoff-thresholded) + w_func·string_func(cutoff-thresholded)`, weights auto-normalized to sum 1. This is what the table is sorted by. | `0.74`. |

*(`pct_joint`, `selectivity_score`, and `common_essential` are carried through from the Essentiality tab and shown alongside these for context.)*

---

## 6. Output columns — Tissue Specificity tab

One row per queried gene.

| Column | Meaning | Example |
|---|---|---|
| **rank** | Row rank after sorting by `n_gates_passed` (how many of the 3 gates below the gene passes, descending) then `selectivity_score` descending. A strict superset of the old pass_all-only ordering: any gene with `pass_all=True` always has `n_gates_passed=3` (tier `high`), so it still sorts first — but a gene passing 2/3 gates now correctly outranks one passing only 0 or 1, instead of being lumped together with every other `pass_all=False` gene. | `1`. |
| *(GTEx tissue columns)* | One column per selected normal tissue, showing that gene's median TPM in that tissue. | `Brain_Cortex: 2.1`. |
| **normal_max** | The gene's highest expression across all selected normal tissues — the single worst-case "how loud is it in normal tissue" number. | `4.3`. |
| **pass_normal** | Whether the gene is below the TPM cutoff across the selected normal tissues (per the aggregation mode: median or all). | `True`. |
| **pancancer_low_fraction** | Fraction of *other* TCGA/TARGET cancer types (excluding the target cohort) where this gene's median expression is below the TPM cutoff — i.e., how tumor-selective it is. `1.0` = silent in every other cancer type; near `0` = broadly expressed across cancers. | `0.75` → quiet in 75% of other cancer types. |
| **pass_pancancer** | Whether `pancancer_low_fraction` meets the "pan-cancer low fraction" threshold slider. | `True`. |
| **target_median_tpm** | Median expression of the gene across samples in the selected target tumor cohort. | `18.6`. |
| **target_fraction_expressed** | Fraction of target-cohort samples where the gene exceeds the TPM cutoff. | `1.0` (expressed in every sample). |
| **pass_target** | Whether the gene passes the "high in the target cohort" gate (via median TPM OR fraction expressed). | `True`. |
| **pass_all** | `pass_normal AND pass_pancancer AND pass_target` — the gene is silent in normal tissue, rare across other cancers, and high in the target cohort: a tissue-specific therapeutic-window candidate. Kept exactly as originally defined (it's what the collaborator-validated PHOX2B/MYCN reference result checks) — see `n_gates_passed`/`specificity_tier` below for a graded alternative. | `True`. |
| **n_gates_passed** | How many of the 3 gates above (`pass_normal`, `pass_pancancer`, `pass_target`) this gene actually passes, 0–3. Purely additive to `pass_all` — distinguishes a gene passing 2/3 gates from one passing 0/3, which `pass_all` alone collapses into the same "False" bucket. Drives the table's sort order (see `rank` above). | `2`. |
| **specificity_tier** | `n_gates_passed` mapped to a coarse label: 0→`fail`, 1→`low`, 2→`medium`, 3→`high`. `high` is exactly equivalent to `pass_all=True`. This is the value shown in the Insights tab's Comparative Summary Table's "Tissue gates" column (see section 8). | `medium`. |
| **selectivity_score** | `target_median_tpm × pancancer_low_fraction / (1 + normal_max)` — rewards high target expression and pan-cancer selectivity while penalizing normal-tissue expression. | `0.254`. |

---

## 7. Output columns — Drug Targets tab

One row per gene-drug pair.

| Column | Meaning | Example |
|---|---|---|
| **Gene** | The queried gene symbol. | `ALK`. |
| **Drug_Name** | Name of the drug/compound targeting this gene. | `CERITINIB`. |
| **Drug_ChEMBL_ID** | ChEMBL compound identifier. | `CHEMBL2403108`. |
| **Max_Phase** | Highest clinical-trial phase reached (4 = FDA-approved); `N/A (DGIdb)` if only DGIdb (not ChEMBL) confirms approval without a phase number. | `4`. |
| **FDA_Approved** | Yes/No — FDA approval status (from `max_phase ≥ 4` or a DGIdb approved flag). | `Yes`. |
| **Interaction_Type** | Nature of the gene-drug interaction (e.g. inhibitor, agonist), as reported by the source database. | `Inhibitor`. |
| **Best_Potency_nM** | Best reported potency (e.g. IC50/Ki) in nanomolar, if available. | `2.5`. |
| **Potency_Type** | Which potency metric `Best_Potency_nM` refers to. | `IC50`. |
| **Query_Cancer** *(cancer-specific mode)* | The cancer type you're checking against. | `Neuroblastoma`. |
| **Indicated_for_Query_Cancer** *(cancer-specific mode)* | Yes/No — whether this drug has a specific indication for your queried cancer type. | `Yes`. |
| **Query_Cancer_Indication_Phase** *(cancer-specific mode)* | Clinical phase of that specific indication. | `2`. |
| **Other_Cancer_Indications** | Other cancer types this drug is indicated for (excluding the query cancer). | `Melanoma; NSCLC`. |
| **Cancer_Indication** *(non-cancer-specific mode)* | Yes/No — whether the drug has *any* cancer indication at all. | `Yes`. |
| **Cancer_Source** | Which source(s) supplied the cancer-indication evidence. | `ChEMBL`. |
| **Cancer_Types / Diseases** | All cancer types / all diseases this drug is associated with. | `Neuroblastoma; Melanoma`. |
| **Gene_Cancer_Associations** *(non-cancer-specific mode)* | Cancer types this *gene itself* (not necessarily the drug) is associated with. | `Neuroblastoma`. |
| **Score** | Overall druggability tier — see the two scoring modes below. | `5`. |
| **Score_Label** *(cancer-specific mode)* | Human-readable label for the `Score` tier. | `Approved for your cancer`. |
| **PRISM_Sensitivity / PRISM_Mean_AUC / PRISM_Top_Lineages** *(PRISM enrichment only)* | DepMap PRISM drug-screening sensitivity call, mean AUC, and the cell-line lineages most sensitive to this compound. | `Sensitive`, `0.42`, `neuroblastoma, lung`. |
| **Sources** | Which upstream databases (ChEMBL, DGIdb, Open Targets, MyGene.info) contributed to this row. | `ChEMBL; Open Targets`. |

**`Score` meaning — cancer-specific 5-tier mode (default):**

| Score | Label |
|---|---|
| 5 | Approved for your cancer |
| 4 | Clinical trial for your cancer |
| 3 | Approved for other cancers |
| 2 | Approved (non-cancer) or trial (any cancer) |
| 1 | Experimental |

**`Score` meaning — simple 3-tier mode:**

| Score | Meaning |
|---|---|
| 3 | FDA-approved and has a cancer indication (or the gene itself is cancer-associated) |
| 2 | FDA-approved (any disease) OR in clinical trial with some cancer evidence |
| 1 | Experimental / preclinical |

---

## 8. Output columns — Model Predictions tab

Displays **precomputed** predictions from a protein-language-model ensemble classifier (built and
validated outside this app — see `.pi/skills/protein-complex-classifier` and the session notes on
its ProteomeLM-S batch-context sensitivity). This tab never runs the classifier live; it loads a
fixed CSV export and looks up whichever genes you enter. One row per requested gene.

| Column | Meaning | Example |
|---|---|---|
| **Prediction** | `INTERACTOR` / `NON-INTERACTOR` from the precomputed file, or **`Not generated/cached`** if the gene isn't in the loaded export — never a fabricated guess. | `INTERACTOR`. |
| **Probability** | The classifier's P(interactor), as stored in the export. Blank/`None` if not cached. | `0.9801`. |
| **Confidence** | The export's own confidence label (`High`/`Medium`/`Low`). | `High`. |
| **UniProt_ID** | UniProt accession the classifier resolved the gene symbol to. | `Q9NPI1`. |
| **cached** | Whether this gene was found in the loaded file at all — the flag driving the fallback above. | `True` / `False`. |

⚠️ **Important caveat carried over from validating this exact classifier this session:** one of its
feature blocks (ProteomeLM-S) is sensitive to which other proteins were embedded in the same
scoring batch. Probabilities are only safely comparable **within the same export/batch** — not
across different CSV exports, and not against a small ad-hoc re-run of the classifier on just a
few genes (which was empirically shown to inflate scores). This tab surfaces exactly one export at
a time for that reason.

---

## 9. Output columns — Insights tab (narrative)

### In-tab controls

| Control | Type | Meaning | Example |
|---|---|---|---|
| **Gene(s)** | text input | Genes to narrate; auto-populated with the current top-3 ranked genes, editable. | `BRD7, DNAJC11, ZC3H13`. |
| **Fetch literature + protein/drug annotations** | checkbox, default **on** | Adds Europe PMC papers + UniProt protein annotation + ChEMBL known drugs per gene, shown in-tab and (for the LLM path) fed to the model. | — |
| **papers per gene** | slider, 3–10, default 6 | How many Europe PMC papers to retrieve per gene. | `6`. |
| **Fetch STRING association with seeds** | checkbox, default **on** | Adds raw functional/physical STRING scores vs. the seed complex, fed to the LLM. Uses the sidebar's STRING edge cutoff. | — |
| **Include last Tissue Specificity triage in the comparative table (if available)** | checkbox, default **on** | Reuses the most recent Tissue Specificity tab result for any typed genes it covered. Feeds **only** the deterministic Comparative Summary Table's "Tissue gates" column below — deliberately **not** sent to the LLM's own prompt (see note below). | — |
| **Include last Drug Targets scoring as evidence (if available)** | checkbox, default **on** | Reuses the most recent Drug Targets tab result for any typed genes it covered — feeds both the LLM's prose and the comparative table's "Known drugs" column. | — |
| **Fetch BioMCP pathways + tissue atlas (if installed)** | checkbox, default **off** | Reactome/KEGG pathway membership + Human Protein Atlas (HPA) tissue expression and subcellular localization, via the optional `biomcp` CLI (biomcp.org) — install separately with `uv tool install biomcp-cli`; not a Python dependency of this app. Off by default since it needs an external binary. Feeds both the in-tab display and the LLM's prose (unlike tissue-specificity, this data doesn't have a comparable nuance/gate-mismatch risk). Shows an honest "not installed" message, never a fabricated result, if enabled without installing biomcp. | — |
| **Include plots from other tabs (…) in this report** | checkbox, default **off** | Reuses each tab's own cached last-run result to show the Co-expression volcano, Essentiality scatter, ranked-selectivity bar (all with the requested gene(s) circled/recolored crimson), STRING physical heatmap, Tissue Specificity heatmap, Drug Targets bar/heatmap, and Model Predictions classifier bar in a "📊 Plots from other tabs" expander. Co-expression/Essentiality/Model Predictions are always available once ranking exists; Interactions/Tissue Specificity/Drug Targets only if you've run that tab this session. Never triggers new computation from Insights. On-screen this is free; embedding the same plots in the PDF export additionally needs the `kaleido` package **and a real Chrome/Chromium already installed on the machine** (kaleido v1.x no longer bundles one) — if that prerequisite is missing, the PDF adds one clear note and you still get the on-screen plots. | — |
| **LLM backend / API key** (expander) | provider selector + key + base URL fields | Configure which cloud provider/model narrates (see "LLM backend" in `app/README.md`). The **base URL** field lets `anthropic` point at a Claude-compatible endpoint other than api.anthropic.com — e.g. Microsoft Azure AI Foundry. | — |
| **Generate narrative** | button | Produces the deterministic narrative (no model; always available). | — |
| **Generate with LLM** | button | Produces the LLM-generated structured report (see below); requires a configured backend. | — |
| **\U0001f4c4 Export as PDF** | button, appears after each narrative | Downloads exactly what's displayed — deterministic narrative or LLM report + comparative table (+ cross-tab plots, if enabled) — as a branded PDF (`app/pdf_export.py`). One button per narrative type. | — |

### Deterministic narrative

Plain-text summary grounded directly in the ranked-table columns (sections 2–4 above): per-gene co-expression/essentiality facts, common-essential flags, and a joint-percentage-ranked "top candidates" list. Needs no model at all and never changes structure. The headline `selectivity_delta` clause for each gene is wrapped in this app's own bounded `==highlight==` syntax (see "Emphasis / highlighting" below).

### Emphasis / highlighting

Both narrative types use a small, app-owned emphasis vocabulary beyond plain **bold**/*italic*: a bounded `==highlighted text==` syntax (not standard markdown) that this app's own code parses and renders as a soft highlight — never raw HTML/CSS from the LLM. On-screen, `ui.render_highlighted_html()`/`ui.markdown_with_highlights()` escape any stray HTML first, then substitute `==...==` for a styled `<mark>` span. In the PDF, `pdf_export.py`'s `_split_inline()` renders the same spans in crimson. The LLM is instructed to use it sparingly — at most once or twice per gene, for the single most critical takeaway — and only ever around facts already present in the JSON, exactly like every other claim in the report.

### LLM narrative — structured report

When a backend is configured, **Generate with LLM** produces a report in exactly six required `## `-headed sections, covering all requested genes together (not repeated once per gene):

1. **Executive Summary** — one short paragraph per gene, headline `selectivity_delta`.
2. **Cohort-Selective Dependency Evidence** — `pct_essential_cohort` vs `pct_essential_pan`, `selectivity_delta`, `common_essential` flag.
3. **Protein Complex / STRING Association** — physical vs functional-only link to the seed complex (or "no STRING data provided"); also notes Reactome/KEGG pathway membership and HPA subcellular localization/tissue expression when BioMCP data is present.
4. **Drug-Design & Therapeutic Relevance** — protein annotation (domains/sites/PDB) + top-scoring known drug(s) if present.
5. **Literature Evidence** — cites only PMIDs actually present in the data.
6. **Experimental Validation & Testing Guidelines** — 2–4 concrete next lab steps grounded only in the supplied data.

**Tissue specificity is deliberately excluded from this list and from the JSON the model receives.** The `pass_normal`/`pass_pancancer`/`pass_target` gates mix a lenient median-based gate with a stricter max-based ranking score (see section 6) — a nuance too easy for a model to flatten into an overconfident or misleading sentence. That evidence still reaches you, just via the deterministic table below instead of the model's prose.

The model is instructed not to add a seventh section or author its own table; strict grounding rules (only cite genes/numbers/PMIDs present in the JSON; treat `common_essential=true` as a poor selective target; no clinical-utility claims) carry over from the original single-paragraph format.

### Comparative Summary Table (appended after the LLM's text, code-computed)

One row per requested gene, built directly from the same facts fed to (or, for tissue, deliberately withheld from) the model — **never authored by the LLM**, so its numbers are exactly as trustworthy as the rest of the app's ranking. Also used as-is for the deterministic narrative's export.

| Column | Meaning | Example |
|---|---|---|
| **Gene** | The requested gene symbol, or the row shows `not ranked` (never a fabricated row) if it isn't in the ranked table. | `BRD7`. |
| **Rank** | `rank/of_total` from the Essentiality & Ranking table. | `1/19`. |
| **Selectivity score / % essential (cohort) / % essential (pan-cancer) / Selectivity delta / Common-essential?** | Carried straight from the Essentiality & Ranking columns (section 4). | — |
| **STRING physical / STRING functional** | Max STRING physical/functional score vs. the seed set, or `n/a` if STRING wasn't fetched. | `0.99` / `1.00`. |
| **Known drugs** | Prefers the Drug Targets tab's `n_drugs_found` (richer DGIdb/ChEMBL/Open Targets count) when available; falls back to the dossier's bare ChEMBL mechanism count only if Drug Targets data wasn't fetched; `n/a` if neither was fetched. | `133`. |
| **Tissue gates** | The Tissue Specificity tab's `specificity_tier` (`fail`/`low`/`medium`/`high`) for this gene if available, else `n/a`. | `fail`. |

---

## Quick cross-reference: which score is which

The app has **five independently-computed scores** across different tabs — easy to conflate, so here's the disambiguation:

| Score name | Tab | What it measures |
|---|---|---|
| `selectivity_score` (Essentiality) | Essentiality & Ranking | Cohort-selective co-expression + essentiality, from DepMap data alone. **The core ranking metric.** |
| `string_score` | Interactions | The above, re-weighted with STRING physical/functional protein-interaction evidence. |
| `selectivity_score` (Tissue Specificity) | Tissue Specificity | A *different* formula — target-tumor expression vs normal-tissue and pan-cancer expression (GTEx/TCGA), unrelated to DepMap essentiality. |
| `Score` (1–5 or 1–3) | Drug Targets | Clinical/regulatory maturity of a drug against the gene — nothing to do with co-expression or essentiality. |
| `Probability` (classifier) | Model Predictions | A precomputed protein-language-model classifier's interactor probability — batch-context-dependent, not from this app's own pipeline, and not cross-comparable to a different export. |
