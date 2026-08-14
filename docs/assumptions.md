# Assumptions and claim boundaries

## Peaks over threshold

`gpd_pot` assumes that positive excesses over a threshold are adequately
described by one GPD and are suitable units for the likelihood analysis.
Threshold selection, declustering, covariate selection, and measurement quality
must be justified outside the function. Sensitivity across plausible thresholds
is recommended.

A negative `xi` point estimate implies a finite endpoint only under the fitted
GPD and threshold. The property `bounded_estimate` reports that algebraic fact.
`bounded_supported` is stricter: the complete 95% profile interval must be
negative. Neither property demonstrates a physical upper bound.

Optimization is restricted to `xi > -1` because the likelihood can become
unbounded at the sample maximum at or below that boundary. The endpoint
bootstrap interval is conditional on bootstrap fits with negative
shape. Always report `bootstrap_fraction_xi_negative` and the successful
replicate count with that interval. Standard likelihood regularity is
problematic at and below `xi = -0.5`; nsevt emits a warning in this regime.

## Grouped (interval-censored) fit

`gpd_pot_grouped` treats each recorded value as censored to its rounding cell
and conditions on the per-observation selection rule (a value enters the sample
because its *recorded* mark exceeded the threshold, so the smallest true excess
that could have produced it is half a grid step). It assumes a single GPD shape
and scale over the exceedances and that the stated grid widths describe the
recording precision. The shape interval is a profile-likelihood interval under
the same interval-censored likelihood, so it refers to the estimator whose point
value is reported. The endpoint interval is a profile interval on the
reparameterised endpoint, which respects the non-linearity of `M* = u -
sigma/xi` where a percentile bootstrap does not; it remains model-conditional
and is not a physical bound. As with the continuous fit, a negative shape
supports a bounded tail only when the whole 95% profile interval lies below
zero.

## Trend permutation

The fitted alternative is linear in log scale:

    log sigma(t) = beta0 + beta1 t.

Shape is constant in this comparison. The permutation unit is the complete
block label: observations stay in their original block while block time labels
are reassigned. The randomization interpretation requires exchangeability of
those block labels under the no-trend null. Temporal autocorrelation, changes in
sampling practice, or systematic differences among blocks can invalidate that
assumption.

The LR statistic is two-sided. Report the number of permutations and the Monte
Carlo standard error with the plus-one p-value. A large p-value is a
non-rejection, not evidence that beta1 equals zero.

## Power and MDE

Power is conditional on the fitted simulation model, observed number of
exceedances, observed block layout, chosen LR critical value, and effect grid.
Positive and negative MDEs are reported separately. Simulation power does not
protect against a misspecified tail model or measurement bias.

## Multi-source robustness

Sources should be genuinely distinct products with comparable variables,
units, thresholds, support, and temporal coverage. Splitting one record into
non-overlapping eras does not create independent sources.

Cross-source agreement is not causal confirmation. Disagreement is not proof of
an instrumental artifact. The returned status distinguishes reproducibility,
opposite direction, adequate-power non-reproduction, unresolved evidence, and
absence of a reference signal.

## Experimental modules

`block_conformal` aggregates ordered scores by blocks. No general finite-sample
or beta-mixing guarantee is claimed for this implementation. `split_conformal`
has the usual marginal guarantee only when calibration and future scores are
exchangeable and the score construction, including scale fitting, is fixed
independently of calibration.

`twoscale_trend` is a residual circular moving-block bootstrap diagnostic for
empirical quantile functions. Its p-value depends on residual stationarity and
block length and is exploratory.
