# Measurement protocol ${VERSIONE}
Section ${SEZIONE} · ${DATA} · one page
Every number comes from the JSON files in `dati/`: no value is written by hand.
Reference: R = ${R_REF}, λ = ${LAMBDA_800}, z = ${Z_REF}, plasma ${MU_PLASMA}, n = ${N_MED}, T = ${T0}.
Goal: measure latch, thermal, and gate; do not prove the engine.

---

## Setup ${N_SETUP_1} — Latch: ratio and paired noise
**Measurement.** ${N_CICLI} fast half-turns against ${N_CICLI} slow half-turns on the same bead (r = ${R_BEAD}, linker ${D_LINKER}) at the same power. A half-turn is half of the ${GRADI_CICLO} cycle, i.e. one ${GRADI_MEZZO} turn. The fast half-turn (counterweight disengaged) is ${T_MEZZO_VELOCE}, the slow one (counterweight engaged) ${T_MEZZO_LENTO}; the sum is ${T720_S}, the ratio is ${RATIO_LATCH}, up to ${RATIO_UP} with the unvalidated hydrodynamic correction.
**Expected.** Ratio ${RATIO_LATCH} for the ideal geometry; ${RATIO_BEAD} with the translational bead term, ${RATIO_BEAD_TOT} with the full term, ${RATIO_UP} as the upper bound.
**Error.** CV = ${CV_IDEALE}; σR/R = CV·√(${DUE}/${N_CICLI}) = ${SIGMA_RATIO_REL}; per-turn ε = ${EPS_1_GEO} (geometric threshold) and ${EPS_1_OPT} (optimized threshold); ε on sums of ${N_CICLI} = ${EPS_3_GEO} (geometric).
**Falsification.** Ratio outside ${RATIO_LATCH}–${RATIO_UP} by ${SIGMA_2}σR (below ${RATIO_FALSIF_LOW}: no CW/No asymmetry); a measured CV above ${CV_MAX_1_GEO} makes the ${BUDGET} budget unreachable per turn, and above ${CV_MAX_3_GEO} makes it unreachable on sums.
**Note.** The per-turn criterion does not improve with N; the sums improve with the square root of N. At CV = ${CV_TARGET} the per-turn error is ${EPS_AT_302_N1} (outside budget) and the sums give ${EPS_AT_302_N3} (inside). The total-duration test would give ${EPS_TOTALE} and is not used.

## Setup ${N_SETUP_2} — Thermal: sink-mirror bracket
**Measurement.** Steady-state ΔT in blood, R = ${R_REF}, z = ${Z_REF}, at ${LAMBDA_800} and ${LAMBDA_1064}; incident power ${P_INC_800} and ${P_INC_1064}.
**Expected.** At ${LAMBDA_800}: sink ${DT_SINK_800}, no correction ${DT_0_800}, mirror ${DT_MIRROR_800}. At ${LAMBDA_1064}: ${DT_SINK_1064}, ${DT_0_1064}, ${DT_MIRROR_1064}.
**Error.** Model bracket ±${FORCELLA_800} at ${LAMBDA_800} and ±${FORCELLA_1064} at ${LAMBDA_1064}; declared resolution ${RISOLUZIONE_TERMICA}; the criterion is ${T_CRITERIO}.
**Falsification.** ΔT ≥ ${T_CRITERIO} at either wavelength, or ΔT above the mirror by ${SIGMA_2}σ. The ${T_CRITERIO} criterion is chosen, not a biological threshold.
**Margin.** The mirror stays below ${T_CRITERIO} by a factor ${MARGINE_800} at ${LAMBDA_800} and ${MARGINE_1064} at ${LAMBDA_1064}.

## Setup ${N_SETUP_3} — Gate: shear threshold and escape
**Measurement.** Shear ramp in plasma (μ = ${MU_PLASMA}–${MU_PLASMA_ALT}) and in blood (μ = ${MU_BLOOD}, σy = ${SIGMA_Y_MID}); escape at Q = ${Q_01}.
**Expected.** In plasma γ⋆ = ${GAMMA_LIM}, with ${GAMMA_LIM_15} at μ = ${MU_PLASMA_ALT}; in blood ${GAMMA_BLOOD_MID} with M_y = ${PREF_BINGHAM}πR³σy, and band ${GAMMA_BLOOD_LOW}–${GAMMA_BLOOD_HIGH}. At Q = ${Q_01}: U₀ = ${U0_Q01}, τesc = ${TAU_ESC_Q01}.
**Error.** In plasma ±${GAMMA_BAND_HALF} (${GAMMA_BAND_REL}), entirely from the μ band; in blood the band ${GAMMA_BLOOD_LOW}–${GAMMA_BLOOD_HIGH} covers μ = ${MU_BLOOD_LOW}–${MU_BLOOD_HIGH} and σy = ${SIGMA_Y_LOW}–${SIGMA_Y_HIGH}. The escape convention is τD = ξ/krad = ${TAU_D_Q01}, not w₀²/(${QUATTRO}D) = ${W0_SQ_4D}.
**Falsification.** Persistent rotation at ${GAMMA_VENULA} (venule, ${RAPPORTO_VENULA}× the threshold) or γ⋆ outside the band by ${SIGMA_2}σ; no escape within ${FINESTRA_FUGA} at Q = ${Q_01} falsifies the barrier.
**Note.** The Bingham prefactor is ${PREF_BINGHAM}πR³σy, conservative; the exact sphere solution would be a factor ${BINGHAM_RATIO} lower. Whole blood remains an effective rheology, not a measured curve.

## Common criteria
- **Threshold before data.** The geometric threshold is fixed before looking at the times; the optimized one is reported only as a check.
- **Pairing.** Every fast cycle is paired with its slow cycle; the total-duration test is not used.
- **Budget.** The budget is ${BUDGET} per comparison; the maximum CV is ${CV_MAX_1_GEO} per turn and ${CV_MAX_3_GEO} on sums.
- **Report.** Value, error, and falsification are reported together; an out-of-band result is declared, not discarded.

## What this protocol does not prove
- **It does not prove biological safety.** ${T_CRITERIO} is a chosen thermal criterion, not a damage threshold.
- **It does not prove the latch.** It measures ratio and noise; it does not measure actuation, energy, or fatigue.
- **It does not prove operation in flow.** It measures the stasis threshold; whole blood remains an effective rheology, not a measured curve.
