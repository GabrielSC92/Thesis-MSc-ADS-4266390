# When Is Dyadic Transformation Reliable? A Threshold for Relational Event Models with Multi-Receiver Events

**MSc Applied Data Science -- Utrecht University (2026)**
Gabriel Silva Cheinquer

## About

The Relational Event Model (REM) is the standard framework for analysing relational event history data, but it is defined entirely on directed dyads (one sender, one receiver). Many real communication processes are *multi-receiver*: a single sender addresses several receivers in one event. A common workaround is to transform each multi-receiver event into several dyadic events separated by a small tie-breaking timestamp noise.

Using a controlled simulation with known dyadic generating parameters, this thesis establishes the **threshold**: the proportion of multi-receiver events up to which the dyadic transformation still recovers the underlying REM parameters reliably, and beyond which it cannot.

**Headline finding.** A three-zone operating range: recovery is reliable up to a multi-receiver share of about 8%, borderline between 8% and roughly 12%, and unreliable above 12%, with the triadic outgoing-two-paths effect binding first.

## Approach

The study is **simulation-first and uses synthetic data exclusively** (no empirical dataset is analysed). The configuration is fixed at `N = 20` actors and `M = 2,000` events, with multi-receiver (MR) levels at **0 / 2 / 8 / 16 / 24%**.

1. Generate synthetic dyadic baseline event histories with known parameters using R (`remulate`), one per independent seed
2. Construct multi-receiver versions by **additive random tagging**: add one extra (signal-free) receiver to a controlled proportion of events at each MR level
3. Apply the **tie-breaking timestamp-noise** transformation to convert multi-receiver events back to single-receiver (dyadic) form while preserving the full-likelihood exact-timing requirement
4. Fit a standard dyadic REM to the transformed data using R (`remify`, `remstats`, `remstimate`)
5. Fit a native Relational Hyperevent Model (RHEM) via [eventnet](https://github.com/juergenlerner/eventnet) + `survival` (stratified Cox) as a **contextual benchmark**
6. Compare recovered REM parameters against the **known dyadic generating truth** across 28 independent baseline datasets, and apply a per-effect recovery-noise band to determine the reliability threshold

## Repository Structure

```
thesis/                         # Thesis paper (DOCX + PDF)
code/
  thesis/
    01_generation_and_baselines/      # Synthetic generation + native benchmark workflows
    02_transformations_and_estimation/ # Dyadic transformation and REM estimation workflows
    03_evaluation_and_reporting/      # Threshold/diagnostic evaluation, tables, and figures
  data/
    synthetic/                  # Synthetic datasets (baselines, MR levels, eventnet exports, results)
```

## Tools

- **R** (4.6+) -- synthetic generation (`remulate`), dyadic REM fitting (`remify`, `remstats`, `remstimate`), and the RHEM benchmark (`survival`)
- **Quarto** -- render `.qmd` notebooks to HTML
- **Java** + **eventnet 1.3** (`eventnet-1.3.jar`) -- RHEM design matrices and the cross-dataset batch harness (headless via `java -jar`)

## Reproducibility

All commands are run from the repository root. MR levels are fixed at **0 / 2 / 8 / 16 / 24%**; the baseline is **N = 2,000 events**, **20 actors**. The full pipeline is deterministic in the baseline seed, so every reported number can be regenerated from the seed list alone.

### Headline result (cross-dataset threshold)

The primary thesis result is the dataset-averaged threshold over 28 independent baseline seeds.

| Step | What to run | Output / notes |
| ---- | ----------- | -------------- |
| 1 | `Rscript code/thesis/01_generation_and_baselines/cross_dataset_eventnet_batch_v01.R` | **~22 min.** Generates baselines, MR-injects, fits REM and RHEM per (seed × MR level). Writes `code/data/synthetic/cross_dataset/{rem,rhem}_coefs_by_seed.csv`. Optional smoke test: `--n-seeds=2`. eventnet runs headless (no GUI). |
| 2 | `quarto render code/thesis/01_generation_and_baselines/cross_dataset_averaging_v01.qmd` | Aggregation, threshold, and recommendation; writes `cross_dataset_summary_v01.csv` |
| 3 | From `code/thesis/03_evaluation_and_reporting/`: `Rscript _briefing_figures.R` | PNG figures (must be run from that folder) |

### Single-seed workflow (illustrative)

A single canonical seed (`20260423`) reproduces the REM-vs-RHEM comparison on one baseline.

| Step | What to run | Output / notes |
| ---- | ----------- | -------------- |
| 1 | `quarto render code/thesis/01_generation_and_baselines/baseline_and_mr_generation_v02.qmd` | Baseline + MR-level CSVs in `code/data/synthetic/mr_levels_v02/` |
| 2 | `quarto render code/thesis/01_generation_and_baselines/rem_parameter_recovery_v01.qmd` | Validates REM recovery on the 0% baseline |
| 3 | `quarto render code/thesis/01_generation_and_baselines/rhem_eventnet_prep_v01.qmd` | eventnet input CSVs; re-render after eventnet produces design matrices |
| 4 | `quarto render code/thesis/02_transformations_and_estimation/noise_reshuffle_rem_v02.qmd` | REM-on-transformed-data estimates in `code/data/synthetic/rem_noise_estimates/` |

## Supervision

- **Daily supervisor:** Mahdi Shafiee Kamalabad (Utrecht University)
- **Secondary supervisor:** Myrthe Prins (Utrecht University)

## License

- **Code** (everything under `code/`, excluding data and third-party tools):
  MIT License — see [LICENSE](LICENSE).
- **Thesis manuscript** (`thesis/`): licensed under
  [CC-BY-NC-ND-4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/),
  matching the Utrecht University Student Theses Repository deposit.
- **Documentation** (this README and other project notes): licensed under
  [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/).

**Not covered** by the above (retains its original license):
- **eventnet** (`eventnet-1.3.jar`) — see the
  [eventnet project](https://github.com/juergenlerner/eventnet).

The synthetic datasets under `code/data/` are generated by the code in this
repository and are released under the same terms as the code.

## How to cite

Cheinquer, G. S. (2026). *When Is Dyadic Transformation Reliable? A Threshold
for Relational Event Models with Multi-Receiver Events* [Master's thesis,
Utrecht University]. Utrecht University Student Theses Repository.
