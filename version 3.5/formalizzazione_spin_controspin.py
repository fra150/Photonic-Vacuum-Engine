"""Punto 4 degli sviluppi futuri del paper
"Modellazione geometrica 3D del sistema Spin / Contro-spin per il Motore Fotonico v3"

FORMALIZZAZIONE MATEMATICA DEI RAPPORTI TRA SPIN E CONTRO-SPIN
===============================================================
Questo script produce i numeri e la figura che accompagnano la nota
Nota_formalizzazione_spin_controspin.md.

Contenuto:
 1. copertura doppia SU(2): U(2pi) = -1, U(4pi) = +1
 2. bilancio del momento angolare: omega_2/omega_1 = -I_1/I_2
    (momenti di inerzia calcolati dalla geometria del modello)
 3. condizione di cancellazione della coppia e della spinta laterale
 4. spinta assiale F_z = (P/c)(eta_s + eta_c cos alpha)
 5. struttura quasi-periodica del moto combinato (due frequenze)

Output: Figura_11_formalizzazione.png, parametri_formalizzazione.json
"""

import json
import math
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------------------------------
# 1. Momenti di inerzia dalla geometria del modello (1 u.a. = 100 nm)
# ----------------------------------------------------------------------------
NM_PER_UA = 100.0
UA = NM_PER_UA * 1e-9                 # m
RHO_CU = 8960.0                       # kg/m^3

# <x^2+y^2> misurato sui vertici/punti delle curve della scena (in u.a.^2)
R2_SPIN = 2.4807
R2_CONTRO = 5.7567

# lunghezze dei fili e raggi di sezione (u.a.)
ORBITE = [(2.20, 1.10), (2.95, 1.52), (3.70, 1.95)]
TAU_ORBITA = 0.032
TAU_LENTE = 0.045
TAU_DISCO = 0.075
R_DISCO = 3.20
LENTE = (1.05, 2.45)


def perimetro_ellisse(a, b):
    """Approssimazione di Ramanujan."""
    return math.pi * (3 * (a + b) - math.sqrt((3 * a + b) * (a + 3 * b)))


L_SPIN = sum(perimetro_ellisse(a, b) for a, b in ORBITE)
L_LENTI = 2 * perimetro_ellisse(*LENTE)
L_DISCO = 2 * math.pi * R_DISCO

V_SPIN = L_SPIN * math.pi * TAU_ORBITA ** 2
V_CONTRO = (L_DISCO * math.pi * TAU_DISCO ** 2
            + L_LENTI * math.pi * TAU_LENTE ** 2)

M_SPIN = RHO_CU * V_SPIN * UA ** 3
M_CONTRO = RHO_CU * V_CONTRO * UA ** 3

I_SPIN = M_SPIN * R2_SPIN * UA ** 2
I_CONTRO = M_CONTRO * R2_CONTRO * UA ** 2
RAPPORTO = I_SPIN / I_CONTRO

# ----------------------------------------------------------------------------
# 2. Copertura doppia SU(2)
# ----------------------------------------------------------------------------
theta = np.linspace(0.0, 4.0 * np.pi, 1441)
psi_re = np.cos(theta / 2.0)
psi_im = np.sin(theta / 2.0)

# ----------------------------------------------------------------------------
# 3. Spinta assiale
# ----------------------------------------------------------------------------
C_LIGHT = 299792458.0
ETA_S, ETA_C = 0.70, 0.85
ALPHA = math.radians(15.0)


def spinta(P, alpha=ALPHA, eta_s=ETA_S, eta_c=ETA_C):
    return (P / C_LIGHT) * (eta_s + eta_c * np.cos(alpha))


# ----------------------------------------------------------------------------
# 4. Figura 11
# ----------------------------------------------------------------------------
fig, axs = plt.subplots(2, 2, figsize=(13.5, 10.5))
fig.suptitle(
    "Motore Fotonico v4 - Formalizzazione dei rapporti Spin / Contro-spin\n"
    "copertura SU(2), bilancio del momento angolare, condizione di coppia nulla",
    fontsize=13, fontweight="bold",
)

# --- (a) SU(2) -----------------------------------------------------------------
ax = axs[0, 0]
ax.plot(np.degrees(theta), psi_re, color="#c0392b", lw=2.2, label="Re ψ = cos(θ/2)")
ax.plot(np.degrees(theta), psi_im, color="#2980b9", lw=2.2, label="Im ψ = sin(θ/2)")
ax.axvline(360, color="k", ls="--", lw=1.2)
ax.axvline(720, color="#27ae60", ls=":", lw=1.4)
ax.annotate("U(2π) = −1\n(ψ cambia segno)", xy=(360, -1), xytext=(120, -0.78),
            fontsize=9, arrowprops=dict(arrowstyle="->", lw=1.0))
ax.annotate("U(4π) = +1\n(ritorno allo stato iniziale)", xy=(720, 1),
            xytext=(405, 1.06), fontsize=9,
            arrowprops=dict(arrowstyle="->", lw=1.0))
ax.set_xlabel("angolo di rotazione θ (gradi)")
ax.set_ylabel("ampiezza dello spinore")
ax.set_title("(a) Lo stato di spin e' uno spinore: serve una rotazione di 720°\n"
             "per tornare allo stato iniziale", fontsize=10)
ax.grid(alpha=0.3)
ax.legend(fontsize=9, loc="lower left")
ax.set_xlim(0, 720)
ax.set_ylim(-1.35, 1.35)

# --- (b) bilancio del momento angolare -----------------------------------------
ax = axs[0, 1]
t = np.linspace(0.0, 1.0, 800)
omega1 = 2.0 * np.pi
L1 = I_SPIN * omega1 / 1e-31
L2 = -I_CONTRO * omega1 * RAPPORTO / 1e-31
L2_sbagliato = -I_CONTRO * omega1 * 0.85 / 1e-31
ax.plot(t, L1 * np.ones_like(t), color="#c0392b", lw=2.4, label="L$_1$ (Spin)")
ax.plot(t, L2 * np.ones_like(t), color="#2980b9", lw=2.4,
        label="L$_2$ = −(I$_1$/I$_2$)L$_1$  (Contro-spin)")
ax.plot(t, (L1 + L2) * np.ones_like(t), color="#27ae60", lw=3.4,
        label="L$_{tot}$ = 0  (coppia nulla)")
ax.plot(t, (L1 + L2_sbagliato) * np.ones_like(t), color="#7f8c8d", lw=2.0, ls="--",
        label="L$_{tot}$ residuo se I$_1$/I$_2$ = 0.85")
ax.axhline(0, color="k", lw=0.8)
ax.set_xlabel("frazione del ciclo")
ax.set_ylabel("momento angolare (10$^{-31}$ kg m²/s)")
ax.set_title("(b) ω$_2$/ω$_1$ = −I$_1$/I$_2$ = −%.3f (valore geometrico)\n"
             "il momento angolare totale resta nullo a ogni istante" % RAPPORTO,
             fontsize=10)
ax.grid(alpha=0.3)
ax.legend(fontsize=8, loc="center right")

# --- (c) spinta assiale ---------------------------------------------------------
ax = axs[1, 0]
P = np.linspace(0, 10.0, 200)
for a_deg, col in ((0.0, "#c0392b"), (15.0, "#e67e22"), (30.0, "#2980b9")):
    ax.plot(P, spinta(P * 1e-6, math.radians(a_deg)) * 1e15, color=col, lw=2.2,
            label="α = %.0f°" % a_deg)
ax.plot([1.0], [spinta(1e-6) * 1e15], "o", color="k", ms=9, zorder=5)
ax.annotate("1 µW → %.2f fN" % (spinta(1e-6) * 1e15), xy=(1.0, spinta(1e-6) * 1e15),
            xytext=(2.4, spinta(1e-6) * 1e15 * 1.9), fontsize=9,
            arrowprops=dict(arrowstyle="->", lw=1.0))
ax.set_xlabel("potenza ottica P (µW)")
ax.set_ylabel("spinta assiale F$_z$ (fN)")
ax.set_title("(c) F$_z$ = (P/c)(η$_s$ + η$_c$ cos α)\n"
             "η$_s$ = 0.70, η$_c$ = 0.85; dipendenza lineare dalla potenza", fontsize=10)
ax.grid(alpha=0.3)
ax.legend(fontsize=9)

# --- (d) struttura quasi-periodica del moto combinato ---------------------------
ax = axs[1, 1]
th1 = np.linspace(0.0, 4.0 * np.pi, 1000)
for rapp, col in ((RAPPORTO, "#c0392b"), (0.85, "#7f8c8d")):
    ax.plot(np.degrees(th1), np.degrees(-rapp * th1), color=col, lw=2.2,
            label="I$_1$/I$_2$ = %.3f" % rapp)
ax.axvline(360, color="k", ls="--", lw=1.0)
ax.axvline(720, color="k", ls="--", lw=1.0)
ax.axhline(-360, color="k", ls=":", lw=1.0)
ax.axhline(-720, color="k", ls=":", lw=1.0)
ax.set_xlabel("θ$_1$ (gradi), Spin")
ax.set_ylabel("θ$_2$ (gradi), Contro-spin")
ax.set_title("(d) Il moto combinato vive su un toro: le due rotazioni chiudono\n"
             "insieme solo se I$_1$/I$_2$ e' razionale  →  moto quasi-periodico",
             fontsize=10)
ax.grid(alpha=0.3)
ax.legend(fontsize=9, loc="lower left")

plt.tight_layout(rect=(0, 0, 1, 0.93))
p_fig = os.path.join(OUT, "Figura_11_formalizzazione.png")
plt.savefig(p_fig, dpi=200)
plt.close(fig)

# ----------------------------------------------------------------------------
# 5. Parametri
# ----------------------------------------------------------------------------
parametri = {
    "unita_blender_in_nm": NM_PER_UA,
    "geometria": {
        "lunghezza_fili_spin_ua": L_SPIN,
        "lunghezza_fili_contro_ua": L_DISCO + L_LENTI,
        "volume_spin_ua3": V_SPIN,
        "volume_contro_ua3": V_CONTRO,
        "massa_spin_kg": M_SPIN,
        "massa_contro_kg": M_CONTRO,
        "r2_medio_spin_ua2": R2_SPIN,
        "r2_medio_contro_ua2": R2_CONTRO,
        "I_spin_kg_m2": I_SPIN,
        "I_contro_kg_m2": I_CONTRO,
        "rapporto_I1_I2": RAPPORTO,
    },
    "relazioni": {
        "copertura_SU2": "U(2pi) = -I ; U(4pi) = +I",
        "bilancio_momento_angolare": "omega_2 = -(I_1/I_2) omega_1",
        "coppia_nulla": "tau_net = I_1 alpha_1 + I_2 alpha_2 = 0",
        "spinta_assiale": "F_z = (P/c)(eta_s + eta_c cos(alpha))",
        # Limite fisico fissato dal momento angolare per fotone (532 nm), non P/omega meccanico.
        "coppia_fotonica": "tau_max = P/omega_ottica",
    },
    "valori": {
        "eta_s": ETA_S,
        "eta_c": ETA_C,
        "alpha_deg": math.degrees(ALPHA),
        "F_z_fN_per_uW": spinta(1e-6) * 1e15,
        "rapporto_animazione": 0.85,
    },
}
p_json = os.path.join(OUT, "parametri_formalizzazione.json")
with open(p_json, "w", encoding="utf-8") as fh:
    json.dump(parametri, fh, indent=2, ensure_ascii=False)

print("Figura 11:", p_fig)
print("I_spin  = %.3e kg m2" % I_SPIN)
print("I_contro = %.3e kg m2" % I_CONTRO)
print("I1/I2   = %.4f   (valore geometrico del modello)" % RAPPORTO)
print("omega_2/omega_1 = %.4f" % -RAPPORTO)
print("F_z(1 uW) = %.3f fN" % (spinta(1e-6) * 1e15))
print("Parametri:", p_json)
