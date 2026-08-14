<p align="center">
  <img src="https://raw.githubusercontent.com/GaiskaSalomon/nsevt/main/assets/nsevt-logo.png" alt="nsevt" width="200">
</p>

# nsevt — non-stationary extreme-value tail inference

[![CI](https://github.com/GaiskaSalomon/nsevt/actions/workflows/ci.yml/badge.svg)](https://github.com/GaiskaSalomon/nsevt/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/nsevt.svg)](https://pypi.org/project/nsevt/)
[![Python versions](https://img.shields.io/pypi/pyversions/nsevt.svg)](https://pypi.org/project/nsevt/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![DOI](https://zenodo.org/badge/1328485209.svg)](https://zenodo.org/badge/latestdoi/1328485209)

`nsevt` is a dependency-light Python package for five connected tasks:

1. peaks-over-threshold generalized Pareto (GPD) estimation with an adaptive
   profile-likelihood interval for the shape and a bootstrap of the finite
   endpoint;
2. an interval-censored (grouped) GPD fit for discretised data, with
   profile-likelihood intervals for the shape and the endpoint, which removes
   the bias that rounding induces in both;
3. a likelihood-ratio trend test calibrated by complete-block label
   permutation, with a grid-independent minimum-detectable effect;
4. Monte Carlo power and signed minimum-detectable-effect (MDE) analysis; and
5. a pre-specified multi-source robustness analysis that distinguishes
   non-reproduction with adequate power from an unresolved comparison.

The package uses deliberately measured terminology. A negative GPD shape point
estimate gives a finite **model-conditional statistical endpoint**. `nsevt`
reports support for a negative shape only when the entire 95% profile interval
lies below zero; neither result is called a physical ceiling.

## Install

From PyPI:

```bash
pip install nsevt
```

For development or the optional Streamlit demonstration:

```bash
git clone https://github.com/GaiskaSalomon/nsevt.git
cd nsevt
pip install -e ".[dev,demo]"
```

## Quick start

```python
import numpy as np
import nsevt

rng = np.random.default_rng(7)
threshold, xi, sigma = 40.0, -0.25, 10.0
years = np.repeat(np.arange(1980, 2030), 12)
uniform = rng.uniform(size=years.size)
excess = sigma / xi * ((1 - uniform) ** (-xi) - 1)
values = threshold + excess

fit = nsevt.gpd_pot(values, threshold=threshold, n_boot=300)
print(fit.summary())
print("negative-shape estimate:", fit.bounded_estimate)
print("negative shape supported by 95% CI:", fit.bounded_supported)

trend = nsevt.trend_permutation(excess, years, n_perm=999)
print("LR permutation p:", trend["p_permutation"],
      "+/-", trend["p_permutation_mcse"])

mde = nsevt.min_detectable_effect(
    excess, years, direction="both", n_rep=200, n_perm_calibration=499
)
print("positive MDE:", mde["mde_positive"])
print("negative MDE:", mde["mde_negative"])
print("interpolated EMD:", mde["emd_positive"], mde["emd_positive_ci95"])
```

## Discretised (grouped) data

When values are recorded on a grid (for example wind speeds in 5 kt steps),
fitting the continuous GPD to the rounded values biases the shape and the finite
endpoint. `gpd_pot_grouped` fits the interval-censored likelihood instead and
reports a profile-likelihood interval for both the shape and the endpoint (the
endpoint interval profiles the reparameterised endpoint, not a percentile
bootstrap).

```python
fit = nsevt.gpd_pot_grouped(values, threshold=40, grid=5.0)
print(fit.summary())
print("shape 95% CI:", fit.xi_ci95)
print("endpoint 95% CI:", fit.endpoint_ci95)
```

## Multi-source robustness

Sources must be genuinely distinct products with comparable temporal support;
early and late halves of one record are not substitutes.

```python
result = nsevt.multisource_robustness(
    [("product A", values_a, years_a),
     ("product B", values_b, years_b)],
    threshold=40,
    reference="product A",
)
print(result.trend_status)
print(result.table())
```

Possible trend statuses include `reproduced`, `inconsistent_direction`,
`not_reproduced_with_power`, `not_resolved`, and `no_reference_signal`. Agreement
or disagreement across sources is a robustness result; it does not by itself
attribute a discrepancy to instruments, homogenization, or physical change.

## Stable and experimental functionality

| status | module | purpose |
|---|---|---|
| stable | `nsevt.gpd` | GPD fit, profile interval, conditional endpoint bootstrap, return levels |
| stable | `nsevt.grouped` | interval-censored (grouped) GPD fit; profile intervals for shape and endpoint |
| stable | `nsevt.trend` | LR block-label permutation, power/MDE, descriptive block-bootstrap interval |
| stable | `nsevt.transportability` | multi-source robustness and power-aware status |
| stable with assumptions | `split_conformal` | upper tail bound for exchangeable calibration scores |
| experimental | `block_conformal` | block-aggregate dependence sensitivity diagnostic |
| experimental | `twoscale_trend` | residual-bootstrap distribution-valued trend diagnostic |
| experimental | `wasserstein_decomposition` | numerical quantile-grid energy decomposition |

The exact assumptions and claim boundaries are documented in
[`docs/assumptions.md`](docs/assumptions.md). Experimental APIs are retained for
evaluation but are not part of the package's central inferential claim.

## Reproducibility and tests

The implementation was adapted from research pipelines and then regression-
checked; it is not represented as a verbatim copy. Randomized routines accept a
seed and report the number of successful replicates. Run the validation suite
with:

```bash
pip install -e ".[dev]"
ruff check src tests demo
pytest --cov=nsevt --cov-report=term-missing --cov-fail-under=80
python -m build
python -m twine check dist/*
```

See [`docs/validation.md`](docs/validation.md) for what each test establishes
and, equally importantly, what it does not establish.

## Citation and license

Use [`CITATION.cff`](CITATION.cff) or the Zenodo DOI shown above. `nsevt` is
released under the MIT License; see [`LICENSE`](LICENSE).
