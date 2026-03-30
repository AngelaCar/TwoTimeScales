## R CMD check results

0 errors | 0 warnings | 0 notes

## Resubmission

This is a resubmission of version 1.2.1. The previous version (1.2.0) failed 
the CRAN `noLD` (no long double) check on Linux due to a computationally 
singular matrix error in the `fitpgam()` example. This has been fixed by 
improving the numerical stability of the example code.

## Original resubmission

This is a new minor version (1.2.0) of TwoTimeScales, previously accepted as 1.0.0.
The main additions in this version are:

Two new model families: a log-additive model (fitpgam()) and a varying-coefficient model (fitvcm()), each with associated helper and plotting functions.
A new model comparison function select_model2ts() to compare fitted models and identify the best fit.
Functions to predict values of the hazard for specific values of the time scales and of the covariates are also added.
Minor fixes in functions and documentation.
Unit tests for prepare_data(), fit1ts(), and fit2ts().

The major changes are documented in the NEWS file.

Version 1.1.0 was never submitted to CRAN, but only published on GitHub.

── R CMD check results ──────────────── TwoTimeScales 1.2.0 ────
Duration: 3m 52.9s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

R CMD check succeeded


