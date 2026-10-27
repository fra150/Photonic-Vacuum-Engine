"""Punto 3 degli sviluppi futuri del paper
"Modellazione geometrica 3D del sistema Spin / Contro-spin per il Motore Fotonico v3"

DEFINIZIONE DELLE SCALE DIMENSIONALI
====================================
Il modello Blender e' adimensionale (unita' arbitrarie, u.a.).  Qui si fissa un
unico fattore di conversione, motivato da un vincolo fisico, e si derivano tutte
le grandezze dimensionali del sistema.

Vincolo fisico scelto: il nucleo di rame nanostrutturato deve stare nel range dei
nanoparticellati plasmonici (30-100 nm di diametro).  Il Nucleo_Cu del modello ha
raggio 0.30 u.a.; imponendo r = 30 nm si ottiene

    1 u.a. = 100 nm = 0.1 um

Conseguenze: il disco del Contro-spin ha diametro 655 nm, l'orbita esterna dello
Spin ha raggio 370 nm.  L'ingombro complessivo del sistema e' ~1.9 um, quindi
D/lambda ~ 0.42 a 1550 nm: il dispositivo resta sub-lunghezza d'onda e lavora in
regime di campo vicino.

Output: Figura_10_scale_dimensionali.png, tabella_scale_dimensionale.csv
"""

import csv
import math
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch

OUT = os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------------------------------
# 1. Fattore di scala
# ----------------------------------------------------------------------------
NM_PER_UA = 100.0           # nm per unita' Blender
UA_PER_NM = 1.0 / NM_PER_UA


def nm(v_ua):
    return v_ua * NM_PER_UA


# geometria letta dalla scena (in u.a.)
R_NUCLEO_CU = 0.30
ORBITE = [(2.20, 1.10), (2.95, 1.52), (3.70, 1.95)]   # (semi-asse a, semi-asse b)
FILO_ORBITA = 0.032                                   # bevel_depth -> raggio
R_DISCO = 3.275                                       # 3.20 + 0.075
LENTE = (1.05, 2.45)                                  # semi-assi
NUCLEO_SPIN = 0.40

D_INIZIALE_UA = 9.5        # distanza centri al frame 1
D_CONTATTO_UA = 5.1        # distanza centri al frame di massimo avvicinamento

# tempi dell'animazione
FPS = 24
F_INIZIO, F_CONTATTO = 1, 90
T_ANIM = (240 - 1) / FPS          # 9.958 s: 239 intervalli del ciclo a doppio giro
# velocita' di avvicinamento RELATIVA: entrambi i sistemi si muovono
V_APPROCCIO_UA_S = ((D_INIZIALE_UA - D_CONTATTO_UA) / (F_CONTATTO - F_INIZIO)
                    * FPS)                     # u.a./s

# ----------------------------------------------------------------------------
# 2. Tabella delle dimensioni
# ----------------------------------------------------------------------------
righe = [
    ("Nucleo Cu (nanoparticella)", "raggio", R_NUCLEO_CU),
    ("Nucleo Spin", "raggio", NUCLEO_SPIN),
    ("Spessore filo delle orbite", "diametro", 2 * FILO_ORBITA),
    ("Orbita 1", "semiasse a", ORBITE[0][0]),
    ("Orbita 1", "semiasse b", ORBITE[0][1]),
    ("Orbita 2", "semiasse a", ORBITE[1][0]),
    ("Orbita 2", "semiasse b", ORBITE[1][1]),
    ("Orbita 3", "semiasse a", ORBITE[2][0]),
    ("Orbita 3", "semiasse b", ORBITE[2][1]),
    ("Disco del Contro-spin", "raggio esterno", R_DISCO),
    ("Lente verticale", "semiasse orizzontale", LENTE[0]),
    ("Lente verticale", "semiasse verticale", LENTE[1]),
    ("Distanza dei centri (frame 1)", "d", D_INIZIALE_UA),
    ("Distanza dei centri (contatto)", "d", D_CONTATTO_UA),
]

csv_path = os.path.join(OUT, "tabella_scale_dimensionale.csv")
with open(csv_path, "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter=";")
    w.writerow(["elemento", "grandezza", "valore_u.a.", "valore_nm", "valore_um"])
    for el, gr, v in righe:
        w.writerow([el, gr, "%.3f" % v, "%.1f" % nm(v), "%.4f" % (nm(v) / 1000.0)])

# ----------------------------------------------------------------------------
# 3. Numeri adimensionali e regimi temporali
# ----------------------------------------------------------------------------
LAMBDA = {"633 nm (He-Ne)": 633.0, "1064 nm (Nd:YAG)": 1064.0, "1550 nm (telecom)": 1550.0}
D_INGOMBRO_NM = nm(2 * R_DISCO)          # diametro del disco

REGIMI_HZ = [1.0, 1.0e3, 1.0e6, 1.0e9]
t_cycle_anim = T_ANIM                 # un giro del Contro-spin per ciclo (pari alla durata reale)
v_appr_nm_s = V_APPROCCIO_UA_S * NM_PER_UA

print("Fattore di scala: 1 u.a. = %.0f nm" % NM_PER_UA)
print("Ingombro del disco: %.0f nm" % D_INGOMBRO_NM)
for k, lam in LAMBDA.items():
    print("  D/lambda a %s : %.3f" % (k, D_INGOMBRO_NM / lam))
print("Velocita' di avvicinamento nell'animazione: %.2f nm/s" % v_appr_nm_s)
for f in REGIMI_HZ:
    T = 1.0 / f
    kappa = t_cycle_anim / T
    print("  f = %-9.0f Hz | T = %.3e s | fattore di compressione = %.2e | "
          "v equivalente = %.3e m/s" % (f, T, kappa, v_appr_nm_s * 1e-9 * kappa))

# ----------------------------------------------------------------------------
# 4. Figura 9
# ----------------------------------------------------------------------------
fig = plt.figure(figsize=(14.5, 10.5))
gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.0], hspace=0.32, wspace=0.22)
fig.suptitle(
    "Motore Fotonico v4 - Definizione delle scale dimensionali\n"
    "1 unita' Blender = 100 nm  (vincolo: nanoparticella di Cu con raggio 30 nm)",
    fontsize=13, fontweight="bold",
)

ORO = "#e67e22"
GIALLO = "#d4ac0d"


def quota(ax, p0, p1, testo, colore="k", pos=0.5, dy=0.0, dx=0.0, rot=0):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="<|-|>", mutation_scale=11,
                                 color=colore, lw=1.1, shrinkA=0, shrinkB=0))
    xt = p0[0] + (p1[0] - p0[0]) * pos + dx
    yt = p0[1] + (p1[1] - p0[1]) * pos + dy
    ax.text(xt, yt, testo, color=colore, fontsize=8.5, ha="center",
            va="bottom", rotation=rot,
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=0.7))


# --- (a) Spin quotato -----------------------------------------------------------
ax = fig.add_subplot(gs[0, 0])
for i, (a, b) in enumerate(ORBITE):
    ax.add_patch(Ellipse((0, 0), 2 * a, 2 * b, fill=False, ec=ORO, lw=1.6,
                         ls=["-", "--", ":"][i], alpha=0.95))
ax.add_patch(Circle((0, 0), NUCLEO_SPIN, color="#c0392b", zorder=4))
quota(ax, (0, 0), (ORBITE[2][0], 0), "370 nm", pos=0.64, dy=0.12)
quota(ax, (0, 0), (0, ORBITE[2][1]), "195 nm", pos=0.86)
quota(ax, (-ORBITE[1][0], 0.0), (-ORBITE[1][0], ORBITE[1][1]), "152 nm", pos=0.45)
ax.annotate("nucleo\nr = 40 nm", (NUCLEO_SPIN * 0.7, -NUCLEO_SPIN * 0.7),
            xytext=(1.75, -2.45), fontsize=8.5,
            arrowprops=dict(arrowstyle="->", lw=0.9))
ax.set_xlim(-4.6, 4.6)
ax.set_ylim(-3.0, 2.6)
ax.set_aspect("equal")
ax.grid(alpha=0.25)
ax.set_xlabel("x (u.a.)  -  1 u.a. = 100 nm")
ax.set_ylabel("z (u.a.)")
ax.set_title("(a) Sistema Spin: tre orbite ellittiche (a = 220 / 295 / 370 nm)\n"
             "filo di 6.4 nm di diametro", fontsize=10)

# --- (b) Contro-spin quotato ----------------------------------------------------
ax = fig.add_subplot(gs[0, 1])
ax.add_patch(Ellipse((0, 0), 2 * R_DISCO, 2 * R_DISCO * 0.28, fill=False,
                     ec=GIALLO, lw=2.0))
for segno, ls in ((+1, "-"), (-1, "--")):
    ax.add_patch(Ellipse((0, segno * 0.62), 2 * LENTE[0], 2 * LENTE[1],
                         fill=False, ec="#b7950b", lw=1.6, ls=ls))
ax.add_patch(Circle((0, 0), R_NUCLEO_CU, color="#c0392b", zorder=4))
quota(ax, (0, 0), (R_DISCO, 0), "328 nm", pos=0.86, dy=0.10)
quota(ax, (-LENTE[0], 0), (0, 0), "105 nm", pos=0.5, dy=0.08)
quota(ax, (3.95, -0.62), (3.95, 0.62), "124 nm", pos=0.5, rot=90)
quota(ax, (0, 0), (0, LENTE[1] + 0.62), "307 nm", pos=0.62, dx=-1.55)
ax.annotate("nucleo Cu\nr = 30 nm", (-R_NUCLEO_CU * 0.7, R_NUCLEO_CU * 1.9),
            xytext=(-4.25, -2.55), fontsize=8.5,
            arrowprops=dict(arrowstyle="->", lw=0.9))
ax.set_xlim(-4.6, 4.6)
ax.set_ylim(-3.0, 2.6)
ax.set_aspect("equal")
ax.grid(alpha=0.25)
ax.set_xlabel("x (u.a.)")
ax.set_ylabel("z (u.a.)")
ax.set_title("(b) Sistema Contro-spin: disco (D = 655 nm) e doppia lente\n"
             "verticale (210 x 490 nm)", fontsize=10)

# --- (c) scala logaritmica comparata -------------------------------------------
ax = fig.add_subplot(gs[1, 0])
voci = [
    ("atomo di Cu (diametro)", 0.256, "#7f8c8d"),
    ("nanoparticella Cu (diametro)", nm(2 * R_NUCLEO_CU), "#c0392b"),
    ("persistenza ssDNA", 50.0, "#8e44ad"),
    ("filo delle orbite (diametro)", nm(2 * FILO_ORBITA), "#e67e22"),
    ("lambda/2 a 633 nm", 316.5, "#2980b9"),
    ("raggio del disco", nm(R_DISCO), "#d4ac0d"),
    ("orbita esterna (raggio)", nm(ORBITE[2][0]), "#e67e22"),
    ("lambda 633 nm", 633.0, "#2980b9"),
    ("lambda 1064 nm", 1064.0, "#2980b9"),
    ("lambda 1550 nm", 1550.0, "#2980b9"),
    ("ingombro del disco (diametro)", D_INGOMBRO_NM, "#d4ac0d"),
]
for i, (nome, val, col) in enumerate(voci):
    ax.plot([val], [i], "o", color=col, ms=8)
    ax.hlines(i, 0.2, val, color=col, lw=0.8, alpha=0.45)
    ax.text(val * 1.12, i, ("%.0f nm" % val) if val >= 10 else ("%.2f nm" % val),
            va="center", fontsize=8, color=col)
ax.set_yticks(range(len(voci)))
ax.set_yticklabels([v[0] for v in voci], fontsize=8.5)
ax.set_xscale("log")
ax.set_xlim(0.2, 6000)
ax.set_xlabel("lunghezza (nm, scala logaritmica)")
ax.axvspan(316, 1550, color="#2980b9", alpha=0.06)
ax.set_title("(c) Confronto con le lunghezze di riferimento\n"
             "il dispositivo e' sub-lunghezza d'onda a 1064 e 1550 nm", fontsize=10)
ax.grid(alpha=0.25, which="both", axis="x")

# --- (d) tabella dei regimi temporali ------------------------------------------
ax = fig.add_subplot(gs[1, 1])
ax.axis("off")
intest = ["frequenza", "periodo T", "compressione\nt$_{anim}$/T", "v equivalente", "Re"]
corpi = []
R_DISCO_M = nm(R_DISCO) * 1e-9
NU = 1.0e-6          # m^2/s, viscosita' cinematica dell'acqua
for f in REGIMI_HZ:
    T = 1.0 / f
    kappa = t_cycle_anim / T
    v = v_appr_nm_s * 1e-9 * kappa
    re = v * 2 * R_DISCO_M / NU
    if f >= 1e6:
        fs = "%.0f MHz" % (f / 1e6)
    elif f >= 1e3:
        fs = "%.0f kHz" % (f / 1e3)
    else:
        fs = "%.0f Hz" % f
    corpi.append([fs, "%.1e s" % T, "%.3g" % kappa, "%.2e m/s" % v, "%.2e" % re])

tab = ax.table(cellText=corpi, colLabels=intest, loc="center", cellLoc="center")
tab.auto_set_font_size(False)
tab.set_fontsize(9)
tab.scale(1.0, 1.8)
ax.set_title("(d) Mappatura temporale: ciclo di %.3f s;\n"
             "v fisica = v dell'animazione x compressione" % T_ANIM, fontsize=10)

plt.savefig(os.path.join(OUT, "Figura_10_scale_dimensionali.png"), dpi=210,
            bbox_inches="tight")
plt.close(fig)

print("Figura 10:", os.path.join(OUT, "Figura_10_scale_dimensionali.png"))
print("Tabella  :", csv_path)
print("v_anim = %.2f u.a./s = %.1f nm/s" % (V_APPROCCIO_UA_S, v_appr_nm_s))
