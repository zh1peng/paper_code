# CMM main analysis

**Article:** Depressive Symptoms and Cardiometabolic Multimorbidity: A Multicohort Study of Adults Aged 45 Years or Older in China, South Korea, and the United States

Main analyses from preprocessed, harmonized person-wave data to Table 1 and Figures 1–6.

## Input and setup

Supply `data/master_person_wave.rds`, an R data frame with one row per `(cohort, id, wave)`. See `input_variables.csv` for required variables and coding. Participant data are not included.

Use R 4.4.3 (or compatible) and Python 3.10+. Run all commands from this folder, keeping `lib/` in place. Replace the input path below with your own; use the same output folder throughout.

```shell
Rscript install_dependencies.R
python -m pip install -r requirements.txt
```

## Run each analysis

Run preparation first, then each R analysis followed by its Python plotting script.

```shell
# Prepare analysis samples (required before analyses 1–6)
Rscript 0_prepare.R data/master_person_wave.rds results

# 1. Sample summary: Table 1 and Figure 1
Rscript 1_summary.R data/master_person_wave.rds results
python 1_summary_plot.py results

# 2. Cross-sectional associations: Figure 2
Rscript 2_crosssectional.R data/master_person_wave.rds results
python 2_crosssectional_plot.py results

# 3. AUC and discrimination: Figure 3
Rscript 3_AUC.R data/master_person_wave.rds results 5000
python 3_AUC_plot.py results

# 4. Longitudinal associations: Figure 4
Rscript 4_longitudinal.R data/master_person_wave.rds results
python 4_longitudinal_plot.py results

# 5. Cumulative incidence: Figure 5
Rscript 5_survival.R data/master_person_wave.rds results
python 5_survival_plot.py results

# 6. Joint analysis with inflammation: Figure 6
Rscript 6_inflammation.R data/master_person_wave.rds results
python 6_inflammation_plot.py results
```

The AUC analysis uses 5,000 bootstrap replicates per cohort and may take 1-2 hours.

## Outputs

- `results/tables/`: result CSV files and formatted Table 1 (`Table1.html`).
- `results/figures/`: Figures 1–6 in PNG, PDF and TIFF formats.
- `results/cache/`: intermediate analysis data; keep local because it contains participant-level records.
