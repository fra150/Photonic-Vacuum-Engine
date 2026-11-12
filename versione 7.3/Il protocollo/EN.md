# Measurement protocol V7.3
Section 13 · 12 November 2026 · one page
Every number comes from the JSON files in `dati/`: no value is written by hand.
Reference: R = 0.5 µm, λ = 800 nm, z = 2.0 µm, plasma 1.3 mPa·s, n = 1.33, T = 310 K.
Goal: measure latch, thermal, and gate; do not prove the engine.

---

## Setup 1 — Latch: ratio and paired noise
**Measurement.** 3 fast half-turns against 3 slow half-turns on the same bead (r = 0.3969 µm, linker 1.00 µm) at the same power. A half-turn is half of the 720° cycle, i.e. one 360° turn. The fast half-turn (counterweight disengaged) is 0.5996 s, the slow one (counterweight engaged) 2.327 s; the sum is 2.926 s, the ratio is 3.881, up to 4.492 with the unvalidated hydrodynamic correction.
**Expected.** Ratio 3.881 for the ideal geometry; 4.256 with the translational bead term, 4.335 with the full term, 4.492 as the upper bound.
**Error.** CV = 17.84%; σR/R = CV·√(2/3) = 14.57%; per-turn ε = 1.07×10⁻⁴ (geometric threshold) and 1.01×10⁻⁴ (optimized threshold); ε on sums of 3 = 1.95×10⁻¹¹ (geometric).
**Falsification.** Ratio outside 3.881–4.492 by 2σR (below 2: no CW/No asymmetry); a measured CV above 17.76% makes the 1×10⁻⁴ budget unreachable per turn, and above 30.77% makes it unreachable on sums.
**Note.** The per-turn criterion does not improve with N; the sums improve with the square root of N. At CV = 30.2% the per-turn error is 2.21×10⁻² (outside budget) and the sums give 7.38×10⁻⁵ (inside). The total-duration test would give 1.14×10⁻² and is not used.

## Setup 2 — Thermal: sink-mirror bracket
**Measurement.** Steady-state ΔT in blood, R = 0.5 µm, z = 2.0 µm, at 800 nm and 1064 nm; incident power 0.6628 mW and 0.8507 mW.
**Expected.** At 800 nm: sink 0.3471 K, no correction 0.5086 K, mirror 0.6701 K. At 1064 nm: 0.2827 K, 0.3906 K, 0.4986 K.
**Error.** Model bracket ±0.161 K at 800 nm and ±0.108 K at 1064 nm; declared resolution 0.05 K; the criterion is 5.0 K.
**Falsification.** ΔT ≥ 5.0 K at either wavelength, or ΔT above the mirror by 2σ. The 5.0 K criterion is chosen, not a biological threshold.
**Margin.** The mirror stays below 5.0 K by a factor 7.46 at 800 nm and 10.03 at 1064 nm.

## Setup 3 — Gate: shear threshold and escape
**Measurement.** Shear ramp in plasma (μ = 1.3 mPa·s–1.5 mPa·s) and in blood (μ = 8.0 mPa·s, σy = 10.0 mPa); escape at Q = 0.1.
**Expected.** In plasma γ⋆ = 2.700 s⁻¹, with 2.340 s⁻¹ at μ = 1.5 mPa·s; in blood 0.2777 s⁻¹ with M_y = 4πR³σy, and band 0.1037 s⁻¹–0.4094 s⁻¹. At Q = 0.1: U₀ = 8.59 kT, τesc = 200.7 s.
**Error.** In plasma ±0.180 s⁻¹ (6.67%), entirely from the μ band; in blood the band 0.1037 s⁻¹–0.4094 s⁻¹ covers μ = 7.0 mPa·s–9.0 mPa·s and σy = 5.0 mPa–20.0 mPa. The escape convention is τD = ξ/krad = 0.0374 s, not w₀²/(4D) = 0.321 s.
**Falsification.** Persistent rotation at 267.0 s⁻¹ (venule, 98.88× the threshold) or γ⋆ outside the band by 2σ; no escape within 2007 s at Q = 0.1 falsifies the barrier.
**Note.** The Bingham prefactor is 4πR³σy, conservative; the exact sphere solution would be a factor 1.5 lower. Whole blood remains an effective rheology, not a measured curve.

## Common criteria
- **Threshold before data.** The geometric threshold is fixed before looking at the times; the optimized one is reported only as a check.
- **Pairing.** Every fast cycle is paired with its slow cycle; the total-duration test is not used.
- **Budget.** The budget is 1×10⁻⁴ per comparison; the maximum CV is 17.76% per turn and 30.77% on sums.
- **Report.** Value, error, and falsification are reported together; an out-of-band result is declared, not discarded.

## What this protocol does not prove
- **It does not prove biological safety.** 5.0 K is a chosen thermal criterion, not a damage threshold.
- **It does not prove the latch.** It measures ratio and noise; it does not measure actuation, energy, or fatigue.
- **It does not prove operation in flow.** It measures the stasis threshold; whole blood remains an effective rheology, not a measured curve.
