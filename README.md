# Identifying Active Effects in Unreplicated Two-Level Designs

R code and reference output for the article

> A Z M Shafiullah and Md. Belaiet Hossain. *Identifying Active Effects in Unreplicated Two-Level Designs: A Calibrated Comparison and a Label-Invariant Qn Estimator.* Submitted to the *Journal of Statistical Computation and Simulation* (2026).

The script implements and verifies nine procedures for identifying active effects in unreplicated two-level designs:

- Daniel (1959)
- Zahn (1975a)
- Lenth (1989)
- step-down Lenth (Ye, Hamada and Wu 2001)
- Juan and Peña (1992)
- Dong (1993)
- Box and Meyer (1986)
- a median baseline
- the proposed Qn-sym estimator

It also uses the factor-level analysis of Box and Meyer (1993). Every decision rule is calibrated under the null model, and the script reproduces all simulation results, tables and figures of the article.

## Contents

| Path | Description |
|---|---|
| `Screening_Simulation.R` | Complete, self-contained script (base R) |
| `verify_reproduction.R` | Compares a new run cell by cell with the reference output |
| `reference_output/tables/` | CSV tables from the reference run |
| `reference_output/figures/` | Figures from the reference run |
| `reference_output/sessionInfo.txt` | R session of the reference run |
| `reference_output/run_log.txt` | Console log of the reference run |

The header of `Screening_Simulation.R` lists which output file corresponds to each table and figure of the article.

## Requirements

- R >= 4.0, using only the base packages `stats`, `graphics`, `grDevices` and `parallel`.
- Optional: [`BsMD`](https://CRAN.R-project.org/package=BsMD) (version 2023.920 used) for the Box and Meyer (1993) comparator in scenario S3. Without it, that comparator is skipped and all other results are unchanged.

## Usage

```sh
Rscript Screening_Simulation.R
```

Output is written to `Screening_Output/tables` and `Screening_Output/figures`. To change the location, set `OUT_DIR` in Section 0 of the script. Setting `QUICK <- TRUE` performs a reduced run for testing the installation only.

The full run took about 75 minutes on two cores (R 4.3.3, Linux).

## Reproducibility

- All random numbers derive from a fixed seed, and each simulation condition sets its own seed. Results therefore do not depend on the number of cores.
- The reference output was reproduced exactly, with all 20 tables identical, on a second machine running Windows.
- The script runs 198 automatic checks against published values. The count of passed checks is printed and saved in `tables/table_S1_unit_checks.csv`.

To compare a new run with the reference output, set the two paths at the top of `verify_reproduction.R` and run it.

## Citation

Please cite the article and this archive. Citation metadata are in `CITATION.cff`.

## License

MIT License; see `LICENSE`.

## Contact

Md. Belaiet Hossain, mdbelaiethossain@gmail.com
