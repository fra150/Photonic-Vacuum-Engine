"""Punto 2 degli sviluppi futuri del paper
"Modellazione geometrica 3D del sistema Spin / Contro-spin per il Motore Fotonico v3"

STUDIO DELLE INTERAZIONI TRA I FLUSSI DI ROTAZIONE
==================================================
Modello: il sistema Spin e il sistema Contro-spin sono schematizzati come due
vortici di Rankine con circolazione opposta (+G, -G) e separazione d tra i
nuclei.  Il modello e' bidimensionale, incomprimibile: descrive la CINEMATICA
dei flussi, non valida il dispositivo.

Risultati analitici di riferimento:
  v_theta(r) = G r / (2 pi a^2)   per r <  a   (nucleo rigido)
  v_theta(r) = G     / (2 pi r)   per r >= a   (campo irrotazionale)
  Traslazione della coppia:  V = G / (2 pi d), perpendicolare alla linea dei centri
  Energia di interazione:    E_int = (rho G^2 / 2 pi) ln(d/a)

Risultato centrale: nel gap tra i due nuclei i due campi hanno lo STESSO verso
(getto trasversale, +y) mentre la vorticity si annulla; all'esterno i campi
hanno verso opposto e compaiono due punti di stagnazione a (0, +-d/2).
La contro-rotazione produce quindi un getto confinato, non un'attrazione:
in regime inviscido la coppia trasla; solo la viscosita' porta alla
cancellazione mutua della vorticity (Lamb-Oseen).

Output: Figura_9_interazione_flussi.png, dati_interazione_flussi.json
"""

import json
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

OUT = os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------------------------------
# 1. Parametri (coerenti con definizione_scale.py: 1 u.a. = 100 nm)
# ----------------------------------------------------------------------------
U_A_NM = 100.0                      # nm per unita' Blender
R_DISCO_NM = 320.0                  # raggio del disco Contro-spin (3.20 u.a.)
RHO = 998.2                         # kg/m^3 (acqua)
# Il disco in rotazione rigida genera un vortice di Rankine con nucleo a = R_disco
A_CORE = R_DISCO_NM * 1e-9          # m

REGIMI = {
    # omega = 2 pi / 10 s = 0.628 rad/s corrisponde a UN giro in 10 s (valore nominale;
    # il ciclo reale dura 239/24 = 9.958 s, omega esatto = 0.631 rad/s).
    "animazione (1 giro in 10 s)": 2.0 * np.pi / 10.0,
    "dispositivo 1 kHz": 2.0 * np.pi * 1.0e3,
    "dispositivo 1 MHz": 2.0 * np.pi * 1.0e6,
}


def circolazione(omega, R=R_DISCO_NM * 1e-9):
    """G = 2 pi omega R^2  (disco in rotazione rigida)."""
    return 2.0 * np.pi * omega * R ** 2


# ----------------------------------------------------------------------------
# 2. Campo della coppia contro-rotante
# ----------------------------------------------------------------------------
def campo(x, y, d, gamma, a=A_CORE):
    """Velocita' (u, v) e vorticity omega_z della coppia +G (a -d/2) / -G (a +d/2).

    psi = -(G/4pi) ln(r1^2+a^2) + (G/4pi) ln(r2^2+a^2)
    u = dpsi/dy ,  v = -dpsi/dx ,  omega_z = -laplacian(psi)
    """
    x1, x2 = -d / 2.0, +d / 2.0
    r1 = (x - x1) ** 2 + y ** 2 + a ** 2
    r2 = (x - x2) ** 2 + y ** 2 + a ** 2
    k = gamma / (4.0 * np.pi)

    u = k * (-2.0 * y / r1 + 2.0 * y / r2)
    v = k * (2.0 * (x - x1) / r1 - 2.0 * (x - x2) / r2)
    wz = k * (4.0 * a ** 2 / r1 ** 2 - 4.0 * a ** 2 / r2 ** 2)
    return u, v, wz


def energia_interazione(d, gamma, a=A_CORE):
    """E_int = (rho G^2 / 2pi) ln(d/a)  [J/m]."""
    return (RHO * gamma ** 2 / (2.0 * np.pi)) * np.log(d / a)


def velocita_coppia(d, gamma):
    """V = G / (2 pi d)."""
    return gamma / (2.0 * np.pi * d)


def punti_stagnazione(d, gamma):
    """Zeri del campo sull'asse x=0 (per simmetria le stagnazioni stanno li')."""
    yy = np.linspace(-3.0 * d, 3.0 * d, 200001)
    u, v, _ = campo(np.zeros_like(yy), yy, d, gamma)
    sp = np.hypot(u, v)
    idx = np.where((sp[1:-1] < sp[:-2]) & (sp[1:-1] < sp[2:]))[0] + 1
    return [float(yy[i]) for i in idx if abs(yy[i]) < 2.5 * d]


# ----------------------------------------------------------------------------
# 3. Figura 8
# ----------------------------------------------------------------------------
D_NM = 950.0     # separazione di riferimento = 9.5 u.a. (frame 1 della scena)
d0 = D_NM * 1e-9
GAMMA = circolazione(REGIMI["dispositivo 1 MHz"])

n = 560
L = 3.0 * d0
xs = np.linspace(-L, L, n)
ys = np.linspace(-L, L, n)
X, Y = np.meshgrid(xs, ys)
U, V, WZ = campo(X, Y, d0, GAMMA)
SPEED = np.hypot(U, V)

fig, axs = plt.subplots(2, 2, figsize=(13.5, 11))
fig.suptitle(
    "Motore Fotonico v4 - Interazione tra i flussi di rotazione Spin / Contro-spin\n"
    "coppia di vortici contro-rotanti (modello cinematico bidimensionale)",
    fontsize=13, fontweight="bold",
)

# --- (a) linee di corrente: il getto nel gap ------------------------------------
ax = axs[0, 0]
sp = ax.streamplot(xs * 1e9, ys * 1e9, U, V, density=1.5, linewidth=0.85,
                   color=SPEED / SPEED.max(), cmap="viridis", arrowsize=0.9)
fig.colorbar(sp.lines, ax=ax, label="|v| / |v|$_{max}$")
ax.add_patch(Circle((-D_NM / 2, 0), A_CORE * 1e9, color="#e67e22", zorder=5))
ax.add_patch(Circle((+D_NM / 2, 0), A_CORE * 1e9, color="#f1c40f", zorder=5))
ax.annotate("Spin  +Γ", (-D_NM / 2, 0), xytext=(-1.45 * D_NM, 0.55 * D_NM),
            color="#a04000", fontsize=11, fontweight="bold",
            arrowprops=dict(arrowstyle="->", color="#a04000"))
ax.annotate("Contro-spin  −Γ", (+D_NM / 2, 0), xytext=(0.30 * D_NM, 0.62 * D_NM),
            color="#8a6d00", fontsize=11, fontweight="bold",
            arrowprops=dict(arrowstyle="->", color="#8a6d00"))
ax.annotate("", xy=(0, 0.85 * D_NM), xytext=(0, 0.12 * D_NM),
            arrowprops=dict(arrowstyle="-|>", color="#2ecc71", lw=3.0))
ax.text(0.06 * D_NM, 0.55 * D_NM, "getto  V = Γ/(2πd)", color="#2ecc71",
        fontsize=10, fontweight="bold")
for ysp in punti_stagnazione(d0, GAMMA):
    ax.plot(0, ysp * 1e9, "X", color="white", mec="black", ms=13, zorder=6)
ax.text(-0.55 * D_NM, -1.35 * D_NM, "X = punti di stagnazione", color="white",
        fontsize=9)
ax.set_xlim(-1.6 * D_NM, 1.6 * D_NM)
ax.set_ylim(-1.6 * D_NM, 1.6 * D_NM)
ax.set_xlabel("x (nm)")
ax.set_ylabel("y (nm)")
ax.set_title("(a) Nel gap le due rotazioni hanno lo stesso verso:\n"
             "nasce un getto trasversale (d = 950 nm, ω = 1 MHz)", fontsize=10)

# --- (b) vorticity: annullamento al centro del getto ----------------------------
ax = axs[0, 1]
lim = np.abs(WZ).max()
im = ax.imshow(WZ, origin="lower", extent=[xs[0] * 1e9, xs[-1] * 1e9,
                                           ys[0] * 1e9, ys[-1] * 1e9],
               cmap="RdBu_r", vmin=-lim, vmax=lim)
cs = ax.contour(xs * 1e9, ys * 1e9, WZ, levels=[0.0], colors="k",
                linewidths=1.4, linestyles="--")
ax.clabel(cs, fmt={0: "ω$_z$ = 0"}, fontsize=9)
fig.colorbar(im, ax=ax, label="ω$_z$ (1/s)")
ax.add_patch(Circle((-D_NM / 2, 0), A_CORE * 1e9, fill=False, ec="k", lw=1.2))
ax.add_patch(Circle((+D_NM / 2, 0), A_CORE * 1e9, fill=False, ec="k", lw=1.2))
ax.plot(0, 0, "o", mfc="none", mec="w", ms=18, mew=2.0, zorder=6)
ax.text(-0.95 * D_NM, -0.95 * D_NM, "ω$_z$ = 0 al centro del getto",
        fontsize=9, color="k",
        bbox=dict(fc="white", ec="none", alpha=0.75))
ax.set_xlim(-1.6 * D_NM, 1.6 * D_NM)
ax.set_ylim(-1.6 * D_NM, 1.6 * D_NM)
ax.set_xlabel("x (nm)")
ax.set_ylabel("y (nm)")
ax.set_title("(b) Le vorticity hanno segno opposto: si annullano\n"
             "dove i flussi si sommano (annichilazione di ω, non di v)", fontsize=10)

# --- (c) profilo di velocita' e vorticity lungo la linea dei centri -------------
ax = axs[1, 0]
xline = np.linspace(-1.5 * d0, 1.5 * d0, 4000)
_, vline, wzline = campo(xline, np.zeros_like(xline), d0, GAMMA)
ax.plot(xline * 1e9, vline / np.abs(vline).max(), color="#27ae60", lw=2.4,
        label="v$_y$(x) / v$_{max}$   (getto)")
ax.plot(xline * 1e9, wzline / np.abs(wzline).max(), color="#2980b9", lw=2.0,
        ls="--", label="ω$_z$(x) / ω$_{max}$   (vorticity)")
ax.axhline(0, color="k", lw=0.8)
ax.axvline(0, color="k", ls=":", lw=1)
ax.annotate("v massima dove ω = 0", xy=(0, 0), xytext=(0.12 * D_NM, -0.62),
            fontsize=9, arrowprops=dict(arrowstyle="->", lw=1.0))
ax.set_xlabel("x lungo la linea dei centri (nm)")
ax.set_ylabel("grandezza normalizzata")
ax.set_title("(c) Profilo lungo l'asse di approccio:\n"
             "la vorticity cambia segno e si annulla nel punto di massima velocita'",
             fontsize=10)
ax.grid(alpha=0.3)
ax.legend(fontsize=9, loc="upper right")

# --- (d) traslazione della coppia ed energia di interazione ---------------------
ax = axs[1, 1]
dv = np.logspace(np.log10(400e-9), np.log10(4000e-9), 300)
l1, = ax.loglog(dv * 1e9, velocita_coppia(dv, GAMMA), color="#2ecc71", lw=2.4)
ax.axvline(D_NM, color="k", ls="--", lw=1)
ax.set_xlabel("separazione dei centri d (nm)")
ax.set_ylabel("V$_{trasl}$ (m/s)", color="#2ecc71")
ax.tick_params(axis="y", colors="#2ecc71")
ax.grid(alpha=0.3, which="both")

ax2 = ax.twinx()
l2, = ax2.loglog(dv * 1e9, np.abs(energia_interazione(dv, GAMMA)),
                 color="#8e44ad", lw=2.0, ls="-.")
ax2.set_ylabel("|E$_{int}$| (J/m)", color="#8e44ad")
ax2.tick_params(axis="y", colors="#8e44ad")

ax.legend([l1, l2], ["V$_{trasl}$ = Γ/(2πd)  [m/s]",
                     "|E$_{int}$| = (ρΓ²/2π)·ln(d/a)  [J/m]"],
          fontsize=9, loc="upper right")
ax.set_title("(d) In regime inviscido d resta costante: la coppia trasla.\n"
             "Solo la viscosita' porta alla fusione dei nuclei", fontsize=10)

plt.tight_layout(rect=(0, 0, 1, 0.94))
path = os.path.join(OUT, "Figura_9_interazione_flussi.png")
plt.savefig(path, dpi=220)
plt.close(fig)

# ----------------------------------------------------------------------------
# 4. Dati numerici
# ----------------------------------------------------------------------------
xchk = np.linspace(-4e-7, 4e-7, 8001)
dati = {
    "nota": "modello cinematico 2D; non e' una validazione fisica del dispositivo",
    "scale": {"unita_blender_in_nm": U_A_NM, "raggio_disco_nm": R_DISCO_NM,
              "rho_kg_m3": RHO, "a_core_nm": A_CORE * 1e9,
              "d_riferimento_nm": D_NM},
    "regimi": {},
}
for nome, om in REGIMI.items():
    g = circolazione(om)
    _, vchk, _ = campo(xchk, np.zeros_like(xchk), d0, g)
    dati["regimi"][nome] = {
        "omega_rad_s": om,
        "Gamma_m2_s": g,
        "V_trasl_d_riferimento_m_s": velocita_coppia(d0, g),
        "gamma_dot_max_1_s": float(np.abs(np.gradient(vchk, xchk)).max()),
        "E_int_d_riferimento_J_m": energia_interazione(d0, g),
    }

with open(os.path.join(OUT, "dati_interazione_flussi.json"), "w", encoding="utf-8") as fh:
    json.dump(dati, fh, indent=2, ensure_ascii=False)

print("Figura 9:", path)
for nome, v in dati["regimi"].items():
    print(" %-26s Gamma=%.3e m2/s | V=%.3e m/s | gammadot=%.3e 1/s | Eint=%.3e J/m"
          % (nome, v["Gamma_m2_s"], v["V_trasl_d_riferimento_m_s"],
             v["gamma_dot_max_1_s"], v["E_int_d_riferimento_J_m"]))
