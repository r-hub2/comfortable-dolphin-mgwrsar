NEWS/ChangeLog
-----------------------------
# 1.4.4 2026-10-28

Bug fixes
	•	`predict()` with `method_pred = 'TP'` (re-estimation of a GWR or MGWR at the target points) ignored the control of the training call since 1.3.1: the block that was meant to make `model@mycall$control` "a list" tested `is.list()` on the unevaluated call and always replaced it by an empty list, so the re-estimation used the defaults (`adaptive = FALSE`, `Type = 'GD'`, ...). An adaptive bandwidth (a number of neighbours) was then read as a distance: on coordinates in [0, 1] every weight was about 1 and the prediction was a global OLS; on metric coordinates the weights vanished and the predictions were zero or degenerate. The control is now evaluated (an unevaluable one falls back to an empty list), and `adaptive`, `Type` and `W` are taken from the fitted object. Version 1.1 behaved correctly; 1.3.1 to 1.4.3 are affected. For a mixed model (`Model = 'MGWR'`) the TP prediction was wrong independently of the control: `MGWR()` never excluded the target points from the sample (their rows, stacked first, were read as training rows by the neighbour index built without them), and a mixed estimator cannot be re-run at target points anyway, since its fixed part needs the local fits of every training point. The TP prediction of a mixed model now keeps the fixed coefficients of the training fit and re-estimates the varying coefficients at the target points by a GWR of the partial residual `Y - XC Betac` on the varying covariates, over the training sample. New `test-predict_oracle.R`: every prediction mode (`TP`, `tWtp_model`, `model`, `shepard`) of GWR, space-time GWR, MGWR, `tds_mgwr` and `tds_mgtwr` is compared with an explicit reference computed by hand (dense kernel weights, weighted least squares), with a relative guard between modes. The reference hashes of the TP predictions in `tests/testthat/_pred_hashes.csv`, recorded with the defect, are regenerated, and `test-predict_TP_control.R` compares the TP predictions with a hand-made local fit.

# 1.4.3 2026-10-27

Changed behaviour
	•	`TDS_MGWR()` returns the terminal sweep of the backfitting, the state that met the stopping rule, as MGWR backfitting does (Fotheringham et al. 2017; Yu et al. 2020), instead of the sweep with the smallest in-sample RMSE (versions 1.3.1 to 1.4.2). The sweeps continue from the current state whatever the change of the RMSE, as before, and the stopping rules (tolerance over `extra_iter` sweeps, patience, oscillation checks, `maxit`) are unchanged: only the returned state differs, and only when a sweep after the running minimum raised the RMSE. Selecting the returned sweep by the unpenalised in-sample RMSE favoured the most localised of the states visited, whereas the local AICc governs the bandwidth pairs only. New lever `control_tds$retained_state = "best"` restores the previous behaviour; `HBETA` and `HRMSE` hold every sweep up to the returned one.
	•	`TDS_MGWR()`: the safeguard against saturated candidate fits is now a ratio of residual variances. A candidate univariate fit whose residual variance is below `control_tds$saturation_ratio` (default `exp(-2)`) times the residual variance of the current complete fit is excluded from the search; 0 disables the rule. Until 1.4.2 a candidate was excluded when its AICc fell below half of a reference AICc (the global AICc with `get_AIC = TRUE`, the median of the candidates otherwise): the AICc shifts by `2 n log c` when the response is multiplied by `c`, so the rule depended on the unit of the response (inert on a response of unit variance, it excluded most candidates and destroyed the calibration on a response of standard deviation 0.2, and the selected bandwidths then depended on `get_AIC`). The new rule is unitless and independent of the inference option; the counts (candidates excluded, component updates with an exclusion, updates with every candidate excluded) are returned in the new slot `saturation`, with one summary warning when the rule applied, instead of one warning per update.
	•	`predict()` for `tds_mgwr` / `tds_mgtwr` models (`method_pred = 'model'`, and `'tWtp_model'` with a space-time kernel): the sharpening of the prediction weights (`exposant`) now scales each row by its largest weight before taking the power. The result is unchanged whenever `normW(W^exposant)` was finite, but a target whose weights are all small (Gaussian kernel with a small bandwidth, far target) no longer gets a row of zeros and a non-finite prediction, as a few MGWR predictions did at n = 5,000 in the Monte Carlo experiments. A target without any positive weight (kernel underflow at the kernel stage, or outside every compact support) now takes the coefficient of its nearest training point in the metric of the component (spatial distance in bandwidth units, plus temporal distance in bandwidth units for a space-time component), with a warning, instead of a zero coefficient or a NaN.
	•	`TDS_MGWR()`: the importance-driven update order ranks the covariates by the standard deviation over the observations of their fitted contribution `beta_k x_k`. Until 1.4.2 the score was the mean absolute contribution divided by the standard deviation of `x_k`: the division made the default order depend on the units of unstandardised covariates (the product `beta_k x_k` does not), and a mean absolute contribution is dominated by the level of an uncentred covariate (a time index, a day of year). The new score depends neither on the unit nor on the level of `x_k`; the order, hence the path of the descent, may change for any data set.
	•	`control_tds$ordering` now accepts `"random"` (a new permutation of the covariates at every sweep, the intercept first, reproducible through the seed) besides `"Importance"` (default) and `"cyclic"` (formula order); the value is validated and documented.
	•	New levers `control_tds$importance_order` (`"decreasing"`, the default, or `"increasing"`: direction of the importance-driven order) and `control_tds$update_last` (default `NULL`: names of varying coefficients updated last in every sweep whatever the ordering, the intercept staying first; meant for covariates that are deterministic functions of the kernel coordinates, such as a time index or a day of year).
	•	New lever `control_tds$polish` (default 3): after the search has stopped, whatever its stopping rule, the backfitting is iterated with every bandwidth fixed at the retained state until the relative change of the in-sample RMSE falls below `tol`, at most `polish` sweeps, so that the returned coefficients, projection matrices and standard errors are the fixed point of the backfitting at the selected bandwidths (to `tol`); `polish = 0` returns the state at which the search stopped. A polishing sweep costs one univariate fit per coefficient, against one per candidate pair for a search sweep.

# 1.4.2 2026-10-06

Bug fix
	•	`TDS_MGWR()` with pinned bandwidths (`control_tds$H`), an adaptive spatial kernel and a capped neighbourhood (`control$NN < n`): a coefficient pinned at the largest grid value (or `Inf`) is now fitted as a global (OLS) component, as `predict()` assumes. The global test compared the bandwidth with `n - 2` instead of the capped sentinel `min(n, NN)`, so the component was fitted locally with `NN` neighbours. New test of the out-of-sample prediction of global, spatial-only and temporal-only components with `NN < n`.
	•	`TDS_MGWR()` with `control$SE = TRUE`: the standard errors (`@sev`) and the per-coefficient effective degrees of freedom (`@edf_k`) are now those of the returned model. Since 1.4.1 the returned state is the sweep with the smallest RMSE, which need not be the last one, while `@sev` and `@edf_k` were still taken from the last sweep (relative SE errors up to about 18% on a 400-observation example; coefficients, `@Shat`, `@R_k`, `@tS` and `@AICc` were not affected). With `get_AIC = TRUE`, the standard errors are computed from the returned `@R_k`.
	•	`search_bandwidths()` keeps the best score of each grid round. A "median - 10 MAD" outlier screen discarded a genuine minimum that lay well below the flat plateau of the large bandwidths, so the next round zoomed on the wrong interval and the search returned a worse local minimum (e.g. a GWR with n = 5000 and a fixed Gaussian kernel: h = 0.077, AICc 26100, instead of h = 0.032, AICc 24216). It became visible in 1.4.1, where a zero bandwidth no longer yields a score and the screen therefore saw one value fewer. Degenerate fits are rejected at the source instead: the AICc of a saturated binomial GWR is now `Inf`, as for the Gaussian AICc since 1.4.1.

# 1.4.1 2026-10-01

New estimator
	•	New `gtwr_HWB2010()`: the GTWR of Huang, Wu and Barry (2010) in the parameterisation of the paper (a single spatio-temporal bandwidth `h_ST` and the ratio `tau`), returned as an object of the new class `gtwr` with its own `summary()` method. It is the Gaussian, fixed-bandwidth special case of the `Type = 'GDT'` space-time kernel, calibrated through `search_bandwidths()`.
	•	New helpers `bw_hwb2gdt()` and `bw_gdt2hwb()` convert between `(h_ST, tau)` and the `(h_S, h_T)` bandwidths of the package; `as_gtwr()` converts the best model of `search_bandwidths()`.

Inference
	•	The AICc is now `Inf` when the trace of the hat matrix reaches `n - 1`: a saturated fit was rewarded by a penalty term of the wrong sign instead of being rejected, which could send a bandwidth search to its smallest bandwidths.
	•	`summary()` gains `fdr_method = "spatial_BY"`, a Benjamini-Yekutieli correction of the local t-tests calibrated on the bandwidth (the default remains `"BY"`).

Performance
	•	GWR, mixed GWR and MGWR engines rewritten with a contiguous memory layout: about x20-40 for the univariate engine, x5-6.5 for the multivariate engine and x7.6-10.6 for the mixed core, with unchanged results. The mixed core no longer allocates an unused n x n x K cube.
	•	Kernel weights and their row normalisation are computed in native code (bitwise identical to the former R code).
	•	Session cache of kernel weight matrices, used by `TDS_MGWR()` and `search_bandwidths()`. Its memory budget adapts to the size of the matrices and can be set with `options(mgwrsar.weights_cache_mb = )`; 0 disables the cache.
	•	The multivariate engine no longer forms the Q factor when neither the hat matrix nor the partial hat matrices are requested.
	•	In the top-down descent, compactly supported fixed kernels are truncated to the neighbours that carry non-zero weights at the current scale; `control_tds$trunc_gauss` offers the (approximate, off by default) equivalent for a fixed Gaussian kernel.
	•	Overall, compared with 1.3.2 (n = 1000, one core, same data, same selected model): `TDS_MGWR()` (`nns = 20`, `tol = 1e-3`, OLS start) takes 2.7 s instead of 10.5 s for `tds_mgwr` with a fixed Gaussian kernel, 1.2 s instead of 6.0 s with an adaptive bisquare kernel, and 14.5 s instead of 68.4 s for `tds_mgtwr` with fixed Gaussian space-time kernels. A GWR bandwidth search with `search_bandwidths()` (n = 1000, grid of 10, 3 rounds and golden-section refinement, AICc, one core) takes 0.93 s instead of 1.96 s with an adaptive bisquare kernel, and 1.12 s instead of 1.73 s with a fixed Gaussian kernel, with the same selected bandwidth and AICc.

Dependencies
	•	The package no longer imports SMUT (which brought SKAT, SPAtest and RSpectra with it). The dense products of hat matrices use the same Eigen product, now compiled in the package, or R's `%*%` when R is linked to an optimised BLAS (Accelerate, OpenBLAS, MKL, BLIS, ATLAS, FlexiBLAS, ArmPL), which is then much faster. `options(mgwrsar.matprod = "eigen")` or `"blas"` forces one of the two.

Robustness (adversarial campaign)
	•	The local least-squares engines (GWR and mixed GWR) use a column-pivoted QR: the rank of a nearly singular local system is read reliably and the non-estimable columns are the dependent ones, whatever their position in the formula (a dependent column placed early used to make the later columns non-estimable, with non-finite coefficients in space-time adaptive fits).
	•	Leverages and hat-matrix rows are computed from the orthonormal factor only, never through the inverse of R: they stay in [0, 1] (with rows summing to 1) even with a constant or collinear predictor, where they used to be negative.
	•	Isolated points (fewer neighbours than coefficients) get the global OLS coefficients with 0, not `NA`, on an aliased column; a warning is issued when every local fit falls back to OLS (bandwidth or `NN` too small).
	•	Explicit errors replace index, dimension or NA errors for: an unknown `Model`, `control$Type`, kernel name or `criterion`; a missing `control$Z` for a space-time model; missing values in the data or the coordinates; a non-positive or missing bandwidth; an adaptive bandwidth larger than the sample; a `control$W` of the wrong dimension; `fixed_vars` absent from the model matrix; a fully collinear design in `SAR`; constant coefficients that the varying part makes non-identifiable (mixed models); a constant time variable in `TDS_MGWR()`.
	•	`search_bandwidths()` no longer closes every connection of the session on exit (`closeAllConnections()`): files open for writing, `sink()` and `capture.output()` of the caller were closed by each search.
	•	`search_bandwidths()`: the golden-section refinement could loop forever when the bracket was between one and 1.3 tolerances wide (the rounded interior point fell back on a bound); found on duplicated observations. The bracket must now shrink at every step and the loop is capped.
	•	Adaptive bandwidths are capped where the compact kernels can read them (the (H + 2)-th neighbour), instead of failing with an index error near `NN`.
	•	`fitted()` is now exported for `mgwrsar` objects.
	•	New tests `test-stress_campaign.R`: 26 adversarial data generators (repeated or duplicated locations, aligned points, clusters, isolated point, anisotropic or geographic coordinates, tiny samples, constant / collinear / locally constant / sparse / badly scaled predictors, heavy tails, outliers, heteroskedasticity, constant or perfect response, degenerate time axes) crossed with 55 estimator configurations; every fit must satisfy the invariants of `tools/stress_cases.R` or stop with an explicit message.

Bug fixes
	•	`MGWRSAR()`, `search_bandwidths()`, `golden_search_2d_bandwidth()`, `multiscale_gwr()`, `TDS_MGWR()`, `simu_multiscale()` and `mgwrsar_bootstrap_test()` no longer change the random number generator of the session: they seed their own L'Ecuyer-CMRG generator for internal draws, as before, and now restore the user's generator (kind and state) on exit. A `set.seed()` placed after a fit used to draw from L'Ecuyer-CMRG instead of the session's generator.
	•	`TDS_MGWR()`: the temporal bandwidths of a starting model taken from a nested call are named like the spatial ones, and `predict()` accepts a model whose `Ht` is a single unnamed value (no sweep kept).

TDS algorithms
	•	`TDS_MGWR()`: the control levers are repaired for `Type = 'GDT'`. `control_tds$H` and `control_tds$Ht` pin bandwidths per coefficient and per axis (one value per varying coefficient, or a named vector for a subset; `NA` leaves a bandwidth free, `Inf` makes it global); `control_tds$V` is a grid of neighbour counts; the adaptive temporal grid mirrors the spatial one.
	•	`TDS_MGWR()`: after a rejected sweep the descent continues from the rejected state (non-monotone descent) and the best sweep is the one returned. Restarting from the best state repeats the rejected sweep, since a sweep is deterministic, and freezes the descent (found by a coauthor); 1.3.2 continued from an inconsistent state (coefficients of the best sweep, residuals of the rejected one), which had the same effect by accident. Without a rejected sweep the results are unchanged; with one they may differ from 1.3.2 at the third significant digit of the RMSE. When all bandwidths are pinned, the returned model is the fixed point of the backfitting, independent of the starting model.
	•	`TDS_MGWR()`: `control_tds$TRUEBETA` (true coefficients, for simulations) no longer stops with an error; the slot `HRMSE` of the returned model then holds the RMSE of each coefficient for the starting model and every kept sweep.
	•	`TDS_MGWR()`: fixes with `fixed_vars` (dimension of the per-coefficient trace, bandwidth vectors restricted to varying coefficients).
	•	`TDS_MGWR()`: with a fixed spatial kernel and repeated locations (panel data), the spatial bandwidth grid is now extended below the distance to the first site, down to `control_tds$panel_floor` times that distance (default 0.2 for a Gaussian kernel, 0.5 otherwise; 1 restores the previous grid). The neighbour-count grid stopped at the 3rd/4th site, one lattice step on a regular panel, which censored the bandwidths of the roughest coefficients; a Gaussian kernel, whose bandwidth is a standard deviation, was censored at about three times its optimal scale.
	•	New `control_tds$V_dist`: spatial bandwidth grid given directly in distances (fixed kernel), taken as is.
	•	`control_tds$V` is read with `[[` (it was partially matched by `V_dist`).

Documentation
	•	Four of the five vignette stubs had an indented YAML header and were rendered as raw text; they now render as intended.
	•	`citation("mgwrsar")` no longer refers to another package in its header.

# 1.3.2 2026-03-02

TDS algorithms
	•	Introducing `tds_mgtwr` model for Multiscale GTWR with additive or multiplicative (cyclic or acyclic) spatio-temporal kernels.
	•	Importance-driven update schedule that prioritizes covariates according to their current scale-normalized contribution to the fitted signal.
	•	Fixed several edge cases in `TDS_MGWR()`.

Control-parameter safety
	•	Added stronger guards on neighborhood and search controls.
	•	`NN` is now capped to `n` in both `MGWRSAR()` and `TDS_MGWR()` to avoid invalid k-NN requests.
	•	In `TDS_MGWR()`, `control_tds$nns` is now capped to `min(round(n/8), round(maxit/2))`; if a user-provided value is too large, it is truncated with a warning.

Parallel search reliability
	•	Improved stability of parallel bandwidth-search workflows (foreach-based execution), with safer fallback behavior.

Testing and maintenance
	•	Added/extended tests for MGTWR/TDS-related scenarios.
	
	
# 1.3.1 2026-01-21
Major performance improvements
	•	Complete rewrite of the local fitting engine using RcppArmadillo (pivotal QR + banded optimizations).
This substantially accelerates GWR, mixed-GWR and MGWR estimations, with average speed-ups of ×4 and memory usage reduced by 30–50%.

New bandwidth selection engine
	•	New function search_bandwidth(), a unified wrapper for 1D (space) and 2D (space and time) bandwidth optimisation.
It supports multi-round grid refinement, integrates golden-section search, and relies on forking for parallel evaluation.

Visualisation
	•	New interactive plot methods based on plotly, allowing dynamic exploration of local coefficients, bandwidth paths, and model diagnostics.

TDS algorithms
	•	Significant improvements to tds_mgwr and tds_mgtwr models:
	•	smoother and more stable AICc-based decisions during bandwidth boosting,
	•	refined sequential optimisation for spatial and spatio-temporal kernels,
	•	improved handling of edge cases and isolated observations,
	•	new internal diagnostics for convergence monitoring.

Improved handling of isolated points
	•	Better detection and fallback to OLS for observations receiving zero weight under non-adaptive kernels—avoiding silent numerical instabilities.

Parallelisation robustness
	•	More reliable parallel execution:
	•	cleanup on interruptions,
	•	fallback to sequential execution when requested.

Predictive methods
	•	More stable predict_mgwrsar() logic with safer handling of model@mycall$control, avoiding previous errors when called immediately after bandwidth optimisation.

Improved numerical stability
	•	Several fixes related to:
	•	QR pivoting in local regressions,
	•	normalization of kernel weights,
	•	avoidance of underflow in Gaussian kernels for very small bandwidths.

Cross-platform build stability
	•	Fixes ensuring compatibility on macOS ARM, Linux (GCC ≥12), and Windows Rtools; better BLAS thread control (OPENBLAS, MKL, VECLIB).
	
Reproducibility across platforms
	•	Adoption of the L’Ecuyer–CMRG random number generator with inversion-based normal deviates, ensuring bitwise-stable stochastic behaviour across all platforms and parallel backends.

	
# 1.2 2025-6-24 (unreleased version)
* Introducing tds_mgtwr model allowing to estimate spatio-temporal multiscale GWR (MGTWR)
* Introducing golden_search_2d_bandwidth for automatic bandwidth selection with spatio-temporal GWR (GTWR)

# 1.1 2024-12-24
* Introducing top-down scale/multiscale GWR (tds_mgwr), adaptive top-down scale/multiscale GWR (atds_mgwr) and regular multiscale GWR (multiscale_gwr)
* Introducing generalized GWR for binomial (bionomial and quasibinomial families).
* Introducing GAM/GWR with gradient descent boosting.
* Improving kernels for spatio-temporal GWR to introduce seasonnality (GDT)
* Removing unused experimental General kernel Product functions (GDX,GDC)
* Correcting a bug for MGWRSAR_x_x_x models with spatial autocorrelation when SE=TRUE.
* Correcting a bug for GWR model with a single explanatory variable when doMC=TRUE.
* Correcting a bug for formula without explicit names of variable (like "~.").
* Correcting a bug in GWR with parallel computation (split error in gwr_beta)


# 1.0.5 2023-11-16
* Removing dependency to qlcMatrix
* Introducing experimental multiscale GWR Model
* Introducing experimental GWR with glm family
* Introducing experimental GWR with glmboost Model
* Adding the ability to calculate the trace of S without calculating the standard deviation matrix (control$get_ts parameter)
* bandwidths_mgwrsar subroutines improved
* AICc criteria added for bandwidth search
* Removing the remove_local_outlier method
* Rename all variables "coord" to "coords"
* Rename the coordinates in the data example to c('x','y')


# 1.0.4 2023-03-01
* Improved kernel_matW function by reintroducing possibility of island with no neighbours for non adaptie kernel
* Introduced an experimental version of a backfitting algorithm for the estimation of the Multiscale GWR with a selection of the bandiwth of each variable by cross validation (LOOCV). The predict_mgwrsar function allows to make predictions from a model of the 'multiscale_gwr' class.

# 1.0.3 2020-11-18
* Fixed a bug for plot_mgwrsar
* Fixed a bug for computing SE and edf when system is computationally singular
* Fixed a bug for the function plot_mgwrsar coming from the way of naming the columns of the coordinates matrix returned in the MGWRSAR function.
* Added a warning for cases where isgcv=TRUE and SE=TRUE when calling the MGWRSAR function.
* Improvement of plot_effects function to plot also effect of spatially non-varying parameters
* Improvement of plot_mgwrsar function to allow control of leaflet tile and number of quantile in legend.
* WARNING:  function to predict with target points using TP options have to be checked for model with spatial autocorrelation.


# 1.0.2 2020-06-01

* Fixed a bug in the kernel_matW function when there are duplicted spatal coordinates: duplicated coordinates are jittered with a warning if this is the case.

# 1.0.1 2020-05-04

* All models of mgwrsar package based on local linear regression can now be estimated using a target points set. Several functions that allows to choose an optimal set of target points to obtain a faster approximation of GWR coefficients has been added.

* Predictions on new data can now be done using the jacknife estimation method instead of spatial extrapolation of local coefficients from a preliminarily estimated model. Only the optimal value of the bandwidth modeled from the initial data is then used. In the function 'predict_mgwrsar', if the parameter method_pred = 'TP' (default), the prediction is done by recalculating a MGWRSAR model with the new data as target points keeping the bandwidth at the optimal value chosen with the training data, otherwise if method_pred= ('tWtp_model', 'model', 'shepard') then a matrix is used for the spatial extrapolation of the estimated coefficients, and prediction are done using these extrapolated coefficients (as in the previous version of mgwrsar package).

* 'KNN' function is deprecated and replaced by 'kernel_matW' function that allows to build spatial weight matrix and interaction matrix based on General Kernel Product. In kernel_matW function it's possible to specify the maximum number of neighbors to consider in gaussian kernel (rough gaussian kernel) to increase speed and sparsity of weights matrix.

* Fast computation of local OLS coefficients in the previous version of mgwrsar package (0.1) uses non pivotal computation that may provides undesirable results in presence of strong colinearity. The RCCP 'fastlmLLT_C' function has been replaced by R native lm.fit function in this realease.

* A new ploting function has been added: plot_effect is a function that plots the effect of a variable X_k with spatially varying coefficient, i.e X_k * Beta_k(u_i,v_i) for comparing the magnitude of effects of between variables

# 0.1 2018-05-11
* First release.

