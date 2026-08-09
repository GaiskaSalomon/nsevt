"""Streamlit demo for nsevt tail inference in the browser.

Run with:
    streamlit run demo/app.py
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

import nsevt

st.set_page_config(page_title="nsevt — tail inference", layout="wide")
st.title("nsevt — non-stationary extreme-value tail risk")
st.caption("GPD profile inference · block-label permutation · power/MDE · "
           "multi-source robustness · conformal tail bounds")

# --------------------------------------------------------------------------- #
# data input
# --------------------------------------------------------------------------- #
with st.sidebar:
    st.header("Data")
    mode = st.radio("Source", ["Simulated example", "Upload CSV"])
    threshold = st.number_input("Tail threshold u", value=40.0, step=1.0)
    if mode == "Upload CSV":
        up = st.file_uploader("CSV with columns: value, block (e.g. year)", type="csv")
        df = pd.read_csv(up) if up is not None else None
        vcol = st.text_input("value column", "value")
        bcol = st.text_input("block column", "block")
    else:
        st.markdown("A bounded tail with a mild scale trend.")
        n_year = st.slider("obs per year", 10, 60, 25)
        trend = st.slider("true trend per decade", -0.2, 0.4, 0.10, 0.02)
        seed = st.number_input("seed", value=0, step=1)
        df = None

if mode == "Simulated example":
    rng = np.random.default_rng(int(seed))
    years = np.arange(1980, 2024)
    vals, blocks = [], []
    xi_true, s0 = -0.30, 12.0
    for y in years:
        s = s0 * np.exp(trend * (y - 1980) / 10.0)
        u = rng.uniform(size=n_year)
        z = s / xi_true * ((1 - u) ** (-xi_true) - 1)
        vals.append(threshold + z)
        blocks.append(np.full(n_year, y))
    value = np.concatenate(vals)
    block = np.concatenate(blocks)
elif df is not None:
    value = df[vcol].to_numpy(float)
    block = df[bcol].to_numpy()
else:
    st.info("Upload a CSV to begin, or switch to the simulated example.")
    st.stop()

exc = value > threshold
z, blk = value[exc] - threshold, block[exc]
st.write(f"**{value.size}** observations · **{z.size}** exceedances over u={threshold:g}")

# --------------------------------------------------------------------------- #
# 1. bounded tail
# --------------------------------------------------------------------------- #
c1, c2 = st.columns(2)
with c1:
    st.subheader("1 · GPD shape and statistical endpoint")
    fit = nsevt.gpd_pot(value, threshold, n_boot=400)
    st.metric("shape ξ", f"{fit.xi:.3f}", help=f"95% profile CI {fit.xi_ci95}")
    if fit.bounded_supported:
        st.write("The 95% profile interval supports **ξ < 0**.")
    elif fit.bounded_estimate:
        st.write("The point estimate has **ξ̂ < 0**, but its 95% interval includes zero.")
    else:
        st.write("The point estimate has **ξ̂ ≥ 0**.")
    if fit.bounded_estimate:
        st.metric("model-conditional endpoint estimate", f"{fit.endpoint:.1f}",
                  help=(f"Conditional bootstrap interval {fit.endpoint_ci95}; report with "
                        "the bootstrap fraction ξ<0"))
    st.metric("bootstrap fraction ξ<0", f"{fit.bootstrap_fraction_xi_negative:.3f}")
    st.caption("A finite GPD endpoint is conditional on the threshold and model; it is not a physical ceiling.")
    rate = exc.mean()
    rps = [10, 50, 100, 500, 1000]
    rl = pd.DataFrame({"return period (obs)": rps,
                       "return level": [round(fit.return_level(r, rate=rate), 1) for r in rps]})
    st.dataframe(rl, hide_index=True)

# --------------------------------------------------------------------------- #
# 2. trend + power
# --------------------------------------------------------------------------- #
with c2:
    st.subheader("2 · Trend test (permutation) + power")
    tr = nsevt.trend_permutation(z, blk, n_perm=1000)
    st.metric("trend per decade", f"{tr['trend_per_decade']:+.3f}",
              help=f"σ change over record: {tr['sigma_change_pct']:+.0f}%")
    st.metric("permutation p-value", f"{tr['p_permutation']:.4f}",
              help=f"Monte Carlo SE {tr['p_permutation_mcse']:.4f}")
    mde = nsevt.min_detectable_effect(z, blk, n_rep=100,
                                      n_perm_calibration=299)
    st.metric("min. detectable effect (80% power)",
              "—" if mde["mde_per_decade"] is None else f"{mde['mde_per_decade']:.2f}/dec")
    pc = pd.DataFrame(mde["power_curve"])
    fig, ax = plt.subplots(figsize=(4, 2.4))
    ax.plot(pc["trend_per_decade"], pc["power"], "o-")
    ax.axhline(0.8, ls="--", color="grey")
    ax.axvline(abs(tr["trend_per_decade"]), ls=":", color="red")
    ax.set_xlabel("trend/decade")
    ax.set_ylabel("power")
    ax.set_ylim(0, 1.02)
    st.pyplot(fig)
    st.caption("A non-rejection with low power means *not resolvable*, not *no trend*.")

# --------------------------------------------------------------------------- #
# 3. conformal bound
# --------------------------------------------------------------------------- #
st.subheader("3 · Conformal upper bound for tail observations")
alpha = st.slider("miscoverage α", 0.01, 0.30, 0.10, 0.01)
band = nsevt.split_conformal(value, threshold, alpha=alpha, scale=fit.sigma)
upper = band.predict_upper(fit.sigma)
st.write(f"One-sided **{100*(1-alpha):.0f}%** upper prediction bound at σ={fit.sigma:.2f}: "
         f"**{upper:.1f}** (standardized q={band.q_standardized:.2f}, "
         f"calibration n={band.block_length}). "
         f"In-sample diagnostic coverage: {band.coverage(value[exc], fit.sigma):.3f}.")
st.caption("The finite-sample split-conformal guarantee requires exchangeable scores and a scale "
           "fixed independently of calibration. In-sample coverage is descriptive, not validation.")

# --------------------------------------------------------------------------- #
# 4. multi-source robustness
# --------------------------------------------------------------------------- #
st.subheader("4 · Multi-source robustness")
st.caption("This analysis requires genuinely distinct products covering comparable periods. "
           "Splitting one record into early and late halves is not a valid substitute.")
st.code(
    "nsevt.multisource_robustness(\n"
    "    [(\"product A\", values_a, years_a),\n"
    "     (\"product B\", values_b, years_b)],\n"
    "    threshold=u, reference=\"product A\"\n"
    ")",
    language="python",
)
