# REM with Multiple Receivers: Transforming Multi-Receiver Event Data into Dyadic Form

**MSc Applied Data Science -- Utrecht University (2026)**

## About

This repository contains the code, literature notes, and project documentation for my master's thesis on Relational Event Models (REM) with multiple receivers.

Standard REMs assume each interaction event has a single sender and a single receiver. Many real-world interactions -- emails, team meetings, marketing communications -- involve multiple receivers simultaneously. This thesis investigates whether multi-receiver relational events can be transformed into multiple dyadic (single-receiver) events while preserving valid REM results, and evaluates which transformation strategies best maintain model fit and parameter stability.

## Approach

1. Generate synthetic multi-receiver event histories with known parameters using [eventnet](https://github.com/juergenlerner/eventnet) (Relational Hyperevent Model)
2. Apply different dyadic transformation strategies to convert multi-receiver events into single-receiver form
3. Fit standard dyadic REMs to the transformed data using R (`remulate`)
4. Compare recovered parameters against the RHEM ground truth
5. Validate findings on real-world data (Apollo 13 mission communications)

## Repository Structure

```
THESIS_TRACKER.md              # Project management and meeting log
literature/
  LITERATURE_NOTES.md           # Paper summaries, key concepts, glossary
  papers/                       # Reading list PDFs
code/
  README.md                     # Structure and naming conventions for code work
  thesis/                       # Thesis-relevant code branch
    01_generation_and_baselines/      # Synthetic generation + native benchmark workflows
    02_transformations_and_estimation/ # Dyadic transformation and REM estimation workflows
    03_evaluation_and_reporting/      # Threshold/diagnostic evaluation, tables, and figures
  non_thesis/                   # Non-thesis code branch
    01_exploration/             # Exploratory analyses and scratch work
    02_sandbox/                 # Non-thesis scripts and experiments
  data/                         # Datasets (synthetic + real-world)
```

## Supervision

- **Daily supervisor:** Mahdi Shafiee Kamalabad (Utrecht University)
- **Secondary supervisor:** Myrthe Prins (Utrecht University)

## Tools

- **R** (4.6+) -- REM fitting and analysis (`remulate`, `remify`, `remstats`, `remstimate`, `survival`)
- **Quarto** -- render `.qmd` notebooks to HTML
- **Java** + **eventnet 1.3** (`eventnet-1.3.jar`) -- RHEM design matrices and the cross-dataset batch harness (headless via `java -jar`)

## Reproducibility pipeline

The data-generation phase is **closed** (2026-06-03). The steps below reproduce the thesis synthetic pipeline from scratch. Methodological rationale and decision history live in `THESIS_TRACKER.md`; this section is the practical run order only.

**Prerequisites:** run all commands from the repository root unless noted. MR levels are fixed at **0 / 2 / 8 / 16 / 24%**; baseline size is **N = 2,000 events**, **20 actors**.

### Path A -- canonical single-seed workflow

Used for the 2026-05-27 presentation and the single-seed REM-vs-RHEM comparison (canonical seed `20260423`).

| Step | What to run | Output / notes |
| ---- | ----------- | -------------- |
| 1 | `quarto render code/thesis/01_generation_and_baselines/baseline_and_mr_generation_v02.qmd` | Baseline + MR-level CSVs in `code/data/synthetic/mr_levels_v02/` |
| 2 | `quarto render code/thesis/01_generation_and_baselines/rem_parameter_recovery_v01.qmd` | **Gate:** validates REM recovery on the 0% baseline before MR experiments |
| 3 | `quarto render code/thesis/01_generation_and_baselines/rhem_eventnet_prep_v01.qmd` (Sections 1–5) | eventnet input CSVs in `code/data/synthetic/eventnet_exports/` |
| 4 | Run eventnet on each export (Section 6 of the notebook, or reuse files already in `code/data/synthetic/eventnet_gui_statistics/`) | Design-matrix CSVs per MR level |
| 5 | Re-render `rhem_eventnet_prep_v01.qmd` (Section 7 runs when design matrices exist) | `code/data/synthetic/eventnet_gui_statistics/rhem_coefficient_summary.csv` |
| 6 | `quarto render code/thesis/02_transformations_and_estimation/noise_reshuffle_rem_v02.qmd` | REM-on-transformed-data estimates in `code/data/synthetic/rem_noise_estimates/` (~40 s; refits all MR × replicate cells every render) |

### Path B -- cross-dataset headline result (thesis deliverable)

Used for the dataset-averaged threshold finding (28 baseline seeds, Decision 25). **This is the primary result for the thesis write-up.**

| Step | What to run | Output / notes |
| ---- | ----------- | -------------- |
| 1 | `Rscript code/thesis/01_generation_and_baselines/cross_dataset_eventnet_batch_v01.R` | **~22 min.** Writes `code/data/synthetic/cross_dataset/{rem,rhem}_coefs_by_seed.csv`. Optional smoke test: `--n-seeds=2`. eventnet runs headless (no GUI). |
| 2 | `quarto render code/thesis/01_generation_and_baselines/cross_dataset_averaging_v01.qmd` | Aggregation, threshold, recommendation sentence; writes `cross_dataset_summary_v01.csv` (~seconds) |
| 3 | From `code/thesis/03_evaluation_and_reporting/`: `Rscript _briefing_figures.R` | PNG figures for the Meeting 11 briefing (must be run from that folder) |

Path B is self-contained: the batch script generates baselines, MR-injects, fits REM and RHEM per (seed × MR level), and does not depend on Path A.

### Reporting (optional)

| File | Command |
| ---- | ------- |
| `code/thesis/03_evaluation_and_reporting/meeting11_briefing_v01.md` | `quarto render code/thesis/03_evaluation_and_reporting/meeting11_briefing_v01.md` |
| `THESIS_TRACKER.md` | `quarto render THESIS_TRACKER.md` |

### Key data locations

```
code/data/synthetic/
  mr_levels_v02/                    # Path A: MR-level long tables
  eventnet_exports/                 # Path A: eventnet inputs
  eventnet_gui_statistics/          # Path A: eventnet design matrices + RHEM summary
  rem_noise_estimates/              # Path A: REM noise/reshuffle fits
  cross_dataset/                    # Path B: per-seed coefficients + averaged summary
```

See `code/README.md` for folder numbering conventions and file naming rules.
