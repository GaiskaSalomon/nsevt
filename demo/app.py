"""Streamlit demo for nsevt: honest non-stationary tail-risk in the browser.

Run with:
    streamlit run demo/app.py
"""
from __future__ import annotations

import io

import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

import nsevt

st.set_page_config(page_title="nsevt — honest tail-risk", layout="wide")
st.title("nsevt — non-stationary extreme-value tail risk")
st.caption("Bounded-tail detection · permutation-calibrated trend · power/MDE · "
           "multi-source transportability · block-conformal bands")

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
        vals.append(threshold + z); blocks.append(np.full(n_year, y))
    value = np.concatenate(vals); block = np.concatenate(blocks)
elif df is not None:
    value = df[vcol].to_numpy(float); block = df[bcol].to_numpy()
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
    st.subheader("1 · Bounded-tail (finite ceiling)")
    fit = nsevt.gpd_pot(value, threshold, n_boot=400)
    st.metric("shape ξ", f"{fit.xi:.3f}", help=f"95% profile CI {fit.xi_ci95}")
    st.write("Tail is **bounded**" if fit.bounded else "Tail is **unbounded**")
    if fit.bounded:
        st.metric("finite endpoint", f"{fit.endpoint:.1f}",
                  help=f"95% bootstrap {fit.endpoint_ci95}")
    st.metric("P(ξ<0) bootstrap", f"{fit.bootstrap_fraction_xi_negative:.3f}")
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
    st.metric("permutation p-value", f"{tr['p_permutation']:.4f}")
    mde = nsevt.min_detectable_effect(z, blk, n_rep=150)
    st.metric("min. detectable effect (80% power)",
              "—" if mde["mde_per_decade"] is None else f"{mde['mde_per_decade']:.2f}/dec")
    pc = pd.DataFrame(mde["power_curve"])
    fig, ax = plt.subplots(figsize=(4, 2.4))
    ax.plot(pc["trend_per_decade"], pc["power"], "o-")
    ax.axhline(0.8, ls="--", color="grey"); ax.axvline(abs(tr["trend_per_decade"]), ls=":", color="red")
    ax.set_xlabel("trend/decade"); ax.set_ylabel("power"); ax.set_ylim(0, 1.02)
    st.pyplot(fig)
    st.caption("A non-rejection with low power means *not resolvable*, not *no trend*.")

# --------------------------------------------------------------------------- #
# 3. conformal band
# --------------------------------------------------------------------------- #
st.subheader("3 · Block-conformal prediction band (coverage under dependence)")
alpha = st.slider("miscoverage α", 0.01, 0.30, 0.10, 0.01)
band = nsevt.block_conformal(value, threshold, alpha=alpha)
upper = band.predict_upper(fit.sigma)
st.write(f"One-sided **{100*(1-alpha):.0f}%** upper prediction bound at σ={fit.sigma:.2f}: "
         f"**{upper:.1f}** (standardized q={band.q_standardized:.2f}, "
         f"{band.n_blocks} blocks of length {band.block_length}). "
         f"In-sample coverage: {band.coverage(z, fit.sigma):.3f}.")

# --------------------------------------------------------------------------- #
# 4. transportability (only for uploaded multi-source or demo split)
# --------------------------------------------------------------------------- #
st.subheader("4 · Multi-source transportability (evidence arena)")
st.caption("Split the record into two halves as a stand-in for two 'sources' and ask "
           "whether the trend reproduces. With real data, pass genuinely distinct products.")
half = len(value) // 2
arena = nsevt.transportability(
    [("first half", value[:half], block[:half]),
     ("second half", value[half:], block[half:])],
    threshold=threshold, n_perm=600, n_boot=200, check_power=True)
rows = [{"source": s.name, "n": s.n, "ξ": round(s.xi, 3), "ξ CI95": str(s.xi_ci95),
         "endpoint": None if not np.isfinite(s.endpoint) else round(s.endpoint, 1),
         "trend/dec": round(s.trend_per_decade, 3), "p_perm": round(s.p_permutation, 4)}
        for s in arena.sources]
st.dataframe(pd.DataFrame(rows), hide_index=True)
st.info(arena.verdict)
