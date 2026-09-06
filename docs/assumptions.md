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
and conditions on the selection rule (a value enters the sample because its
*recorded* mark is a grid point above the threshold, so the smallest true excess
that could have produced it is the lower rounding edge of the smallest such
mark: half a grid step when the threshold lies on the grid, less when it does
not). It assumes a single GPD shape and scale over the exceedances and that the
stated grid widths describe the recording precision. The shape interval is a
profile-likelihood interval under the same interval-censored likelihood, found
by searching outward from the estimate so it always contains it; a limit that
only reached the numerical search boundary is flagged rather than reported as a
statistical limit, and an endpoint that is not identified from above is reported
as infinite. The endpoint interval is a profile interval on the
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
non-rejection, not evidence that beta1 equals zero. If any optimizer refit
fails, `trend_permutation` warns that the p-value is conditional on the
successful permutations; it should not be presented as the requested complete
randomization experiment.

## Power and MDE

Power is conditional on the fitted simulation model, observed number of
exceedances, observed block layout, chosen LR critical value, and effect grid.
Positive and negative MDEs are reported separately. Simulation power does not
protect against a misspecified tail model or measurement bias. The interpolated
EMD interval perturbs the pointwise power estimates by their reported Monte
Carlo errors. It is an approximate sensitivity interval, not a joint bootstrap:
it does not model covariance across effect sizes induced by common random
numbers.

## Multi-source robustness

Sources should be genuinely distinct products with comparable variables,
units, thresholds, support, and temporal coverage. Splitting one record into
non-overlapping eras does not create independent sources.

Cross-source agreement is not causal confirmation. Disagreement is not proof of
an instrumental artifact. The returned status distinguishes reproducibility,
opposite direction, adequate-power non-reproduction, unresolved evidence, and
absence of a reference signal.

## Sequential Monte Carlo precision (`nsevt.mc`)

The MCSE formulae assume independent replicates: `mcse_proportion` and the
sequential stopping rule are valid when `draw` returns independent outcomes, and
`batch_diagnostic` is provided precisely to detect when that fails (a ratio far
from 1). The stopping rule controls Monte Carlo (simulation) error, not
statistical (sampling) error: reaching `R*` means the reported estimate is
precise for the design and data generating process supplied, not that the
underlying estimand is correct. `substream`/`block_streams` give reproducible,
independent streams from one master seed; reproducibility is exact only for the
same NumPy version and platform-independent `SeedSequence` spawning.

## Finite-sample calibration (`nsevt.calibration`)

`rejection_rate`, `coverage`, `bias_rmse` and `pseudo_true` measure the behaviour
of the estimator or test **on the data generating process you supply**. They are
diagnostics of a procedure under a stated model, not evidence about the real
world: a coverage of 0.95 against a well-specified DGP does not certify coverage
under the true data, only that the interval is calibrated for that DGP. Under
misspecification an estimator may converge to a pseudo-true value rather than
the generating parameter. `pseudo_true` approximates such a target with one
large simulated sample; it does not establish the existence or uniqueness of a
limit, and sensitivity to `R` and seed should be checked before coverage is
judged against it. Results carry Monte Carlo error (via `nsevt.mc`); read coverage and
type-I numbers to their MCSE, not beyond it.

A replicate whose test or estimator raises, or returns a non-finite p-value,
interval limit or point estimate, is a **failed replicate**: it is dropped and
counted (`n_failed`), not folded in as a non-rejection or non-coverage, and it
does not abort the campaign. `coverage` treats an interval with an infinite
limit that still brackets the target as covering (a valid but uninformative
interval) and reports how many such replicates there were in `n_infinite`.
`rejection_rate` uses one threshold, `alpha + anticonservative_margin`, for both
the stopping rule and the reported `anticonservative` verdict. Invalid controls
(`n <= 0`, `level` outside `(0, 1)`, a non-finite `target`, a margin that pushes
`alpha + margin` outside `[0, 1)`) are rejected before any simulation is drawn.

## Grouped regression and return levels (`nsevt.design`)

`fit_grouped_design` assumes the log-scale is linear in the supplied design
matrix with a common shape, and that `values`/`design` are aligned exceedances;
the covariate coding is the caller's responsibility. `profile_ci_coef` reports a
profile-likelihood interval for one coefficient and relies on the usual
chi-square calibration of the likelihood ratio, which — as for the trend test —
can be optimistic under a realistic null; calibrate it with `nsevt.calibration`
rather than assuming the asymptotic level. `return_level` and
`profile_ci_return_level` are conditional on the fitted GPD and the supplied
exceedance rate; `profile_ci_return_level` is restricted to bounded tails
(`xi < 0`) and profiles the level itself, so the interval is a model-based
summary, not a distribution-free bound, and does not propagate uncertainty in the
threshold or the exceedance rate. A return period with `m * rate <= 1` targets a
quantile at or below the threshold, which the peaks-over-threshold model does not
describe: both `nsevt.design.return_level` and `nsevt.GPDFit.return_level` return
the threshold there (with a `RuntimeWarning` below the boundary), and
`profile_ci_return_level` returns `None`.

## Experimental modules

`block_conformal` aggregates ordered scores by blocks. No general finite-sample
or beta-mixing guarantee is claimed for this implementation. `split_conformal`
has the usual marginal guarantee only when calibration and future scores are
exchangeable and the score construction, including scale fitting, is fixed
independently of calibration. When the calibration sample is too small for the
requested `alpha` (`ceil((1 - alpha)(m + 1)) > m`) the one-sided bound is `+inf`
(`ConformalBand.underpowered` is `True`); it is not clipped to the sample
maximum, whose coverage would be only `m / (m + 1)`.

`twoscale_trend` is a residual circular moving-block bootstrap diagnostic for
empirical quantile functions. Its p-value depends on residual stationarity and
block length and is exploratory.
