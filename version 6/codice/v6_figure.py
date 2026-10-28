"""
MOTORE FOTONICO DEL VUOTO - V6
Figure del preprint. Ogni figura legge i risultati di v6_studi.py, cosi' i
numeri dei grafici e quelli del testo sono gli stessi per costruzione.

Uso:
    python v6_figure.py
"""

import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.ticker import LogLocator

BASE = os.path.dirname(os.path.abspath(__file__))
if BASE not in sys.path:
    sys.path.insert(0, BASE)

import v6_nucleo as n      # noqa: E402
import v6_studi as s       # noqa: E402

FIG = os.path.normpath(os.path.join(BASE, "..", "figure"))
os.makedirs(FIG, exist_ok=True)

# palette
BLU = "#1f4e79"
ARANCIO = "#c55a11"
ROSSO = "#9c2f22"
VERDE = "#3c6e47"
GRIGIO = "#6f6f6f"
VIOLA = "#5b3a70"
OCRA = "#b08d2e"

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif"],
    "font.size": 9,
    "axes.titlesize": 9.5,
    "axes.labelsize": 9,
    "legend.fontsize": 7.6,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.linewidth": 0.8,
    "figure.dpi": 200,
    "savefig.dpi": 220,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.6,
})


def _pulisci(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _salva(fig, nome):
    p = os.path.join(FIG, nome)
    fig.savefig(p, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("  figura:", nome)


# ---------------------------------------------------------------------------
# 1. IL MODELLO UNICO
# ---------------------------------------------------------------------------
def fig_modello(res):
    a = res["residui"]["bilancio_potenza"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.5),
                                   gridspec_kw={"width_ratios": [1.05, 1]})
    # (a) cascata di potenza
    voci = [
        ("$P$ incidente", n.P_LASER),
        ("$P$ intercettata", a["P_intercettata_W"]),
        ("$P$ assorbita", a["P_assorbita_W"]),
        ("$P$ meccanica\n$\\tau\\omega$", a["P_meccanica_estratta_W"]),
    ]
    y = np.arange(len(voci))[::-1]
    val = [v for _, v in voci]
    ax1.barh(y, val, color=[GRIGIO, BLU, VERDE, ARANCIO], height=0.55)
    ax1.set_xscale("log")
    ax1.set_yticks(y)
    ax1.set_yticklabels([k for k, _ in voci])
    ax1.set_xlim(1e-28, 1e-4)
    ax1.set_xlabel("potenza (W)")
    ax1.set_title("(a) cascata di potenza: $1\\ \\mu$W in ingresso")
    for yy, v in zip(y, val):
        ax1.annotate(f"{v:.3g} W", (v, yy), xytext=(6, 0), textcoords="offset points",
                     va="center", fontsize=7.5)
    ax1.axvline(n.KT / 1.0, color="none")
    _pulisci(ax1)

    # (b) la coppia come limite
    P = np.logspace(-8, -2, 200)
    for fac, col, lab in ((1, BLU, "$\\ell+s=1$"), (2, VERDE, "$\\ell=1,\\ s=1$"),
                          (6, ARANCIO, "$\\ell=5,\\ s=1$"), (11, ROSSO, "$\\ell=10,\\ s=1$")):
        ax2.loglog(P, n.coppia_limite(P, l=fac - 1, s=1.0), color=col, lw=1.4, label=lab)
    ax2.axhline(n.KT, color=GRIGIO, ls="--", lw=1.0)
    ax2.annotate("$k_BT$", (2e-8, n.KT * 1.6), fontsize=8, color=GRIGIO)
    ax2.axvline(a["P_assorbita_W"], color=ARANCIO, ls=":", lw=1.0)
    ax2.annotate("$P_{ass}$ di\nriferimento", (a["P_assorbita_W"] * 1.3, 2e-27),
                 fontsize=7, color=ARANCIO)
    ax2.set_xlabel("potenza assorbita (W)")
    ax2.set_ylabel("$\\tau=(\\ell+s)\\,P_{ass}/\\omega_{ott}$  (N m)")
    ax2.set_title("(b) la coppia come limite superiore")
    ax2.legend(loc="upper left", frameon=False)
    _pulisci(ax2)
    fig.tight_layout()
    _salva(fig, "fig01_modello.png")


# ---------------------------------------------------------------------------
# 2. AUDIT OTTICO
# ---------------------------------------------------------------------------
def fig_ottica(res):
    o = res["ottica"]
    fig, axs = plt.subplots(2, 2, figsize=(9.2, 6.4))
    ax = axs[0, 0]
    c = o["sfera_curva"]
    c = [r for r in c if r["regime"] == "Mie esatta"]
    x = np.array([r["x"] for r in c])
    ax.loglog(x, [r["Q_abs"] for r in c], color=BLU, lw=1.5, label="$Q_{ass}$ (Mie)")
    ax.loglog(x, [r["Q_pr"] for r in c], color=ARANCIO, lw=1.5, label="$Q_{pr}$ (Mie)")
    ax.axhline(2.0, color=GRIGIO, ls="--", lw=0.9)
    ax.annotate("limite geometrico $Q_{ext}\\to2$", (2e-3, 2.3), fontsize=7, color=GRIGIO)
    ax.axvline(0.472, color=VERDE, ls=":", lw=1.0)
    ax.annotate("nucleo\n40 nm", (0.52, 0.02), fontsize=7, color=VERDE)
    ax.set_xlabel("parametro di dimensione $x=2\\pi a/\\lambda$")
    ax.set_ylabel("efficienza")
    ax.set_title("(a) sfera di rame: Mie esatta")
    ax.legend(loc="lower right", frameon=False)
    _pulisci(ax)

    ax = axs[0, 1]
    c = o["cilindro"]
    x = np.array([r["x"] for r in c])
    ax.loglog(x, [r["Q_esatto"] for r in c], color=BLU, lw=1.5, label="serie n$\\geq$1")
    ax.loglog(x, [r["Q_rayleigh"] for r in c], color=ROSSO, ls="--", lw=1.2, label="Rayleigh")
    ax.axvline(0.01476, color=VERDE, ls=":", lw=1.0)
    ax.annotate("filo 1,25 nm\n(x = 0,0148)", (0.016, 2e-3), fontsize=7, color=VERDE)
    ax.set_xlabel("$x=2\\pi a/\\lambda$")
    ax.set_ylabel("$Q_{ext}$")
    ax.set_title("(b) filo di rame: esatto contro Rayleigh")
    ax.legend(loc="upper left", frameon=False)
    _pulisci(ax)

    ax = axs[1, 0]
    c = o["ossidazione"]
    sp = [r["spessore_guscio_nm"] for r in c]
    ra = [r["Q_abs"] / r["Q_abs_nudo"] for r in c]
    ax.plot(sp, ra, "o-", color=VERDE, lw=1.5, ms=4)
    ax.axhline(1.0, color=GRIGIO, lw=0.8)
    ax.set_xlabel("spessore del guscio di Cu$_2$O (nm)")
    ax.set_ylabel("$Q_{ass}/Q_{ass}$(rame nudo)")
    ax.set_title("(c) ossidazione del nucleo")
    ax.set_ylim(0.7, 1.6)
    _pulisci(ax)

    ax = axs[1, 1]
    c = o["eps_dimensione"]
    rr = [r["r_nm"] for r in c]
    ax.semilogx(rr, [r["k"] for r in c], "o-", color=BLU, lw=1.4, ms=4, label="$k$ con $\\gamma(a)$")
    ax.axhline(2.53, color=BLU, ls="--", lw=0.9, label="$k$ bulk")
    ax2 = ax.twinx()
    ax2.semilogx(rr, [r["scarto_rel"] * 100 for r in c], "s-", color=ARANCIO, lw=1.2, ms=3.5,
                 label="scarto $Q_{ass}$")
    ax2.axhline(0, color=GRIGIO, lw=0.7)
    ax2.set_ylabel("$\\Delta Q_{ass}$ (%)", color=ARANCIO)
    ax2.tick_params(axis="y", colors=ARANCIO)
    ax2.grid(False)
    ax.set_xlabel("raggio (nm)")
    ax.set_ylabel("$k$ del rame", color=BLU)
    ax.set_title("(d) correzione di dimensione (Drude)")
    ax.legend(loc="lower left", frameon=False)
    _pulisci(ax)
    fig.tight_layout()
    _salva(fig, "fig02_ottica.png")


# ---------------------------------------------------------------------------
# 3. BILANCIO TERMICO
# ---------------------------------------------------------------------------
def fig_termica(res):
    sw = res["termica"]["sweep"]
    P = np.array([r["P_W"] for r in sw])
    c0 = [c for c in res["termica"]["casi"] if c["materiale"] == "rame"][0]
    em1, em2 = c0["emiss_limpido"], c0["emiss_ossidato"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.6))

    ax1.loglog(P, [r["dT_fluido_K"] for r in sw], color=BLU, lw=1.6)
    ax1.axhline(100.0, color=ROSSO, ls="--", lw=1.0)
    ax1.annotate("ebollizione ($\\Delta T=100$ K)", (3e-6, 130), fontsize=7, color=ROSSO)
    ax1.axvline(1e-6, color=GRIGIO, ls=":", lw=1)
    ax1.annotate("1 $\\mu$W", (1.2e-6, 3e-3), fontsize=7, color=GRIGIO)
    ax1.set_xlabel("potenza incidente (W)")
    ax1.set_ylabel("$\\Delta T$ in acqua (K)")
    ax1.set_title("(a) fluido: conduzione, $Nu=1$")
    _pulisci(ax1)

    ax2.loglog(P, [r["T_vuoto_limpido_K"] for r in sw], color=ROSSO, lw=1.6,
               label=f"$\\epsilon={em1:.2f}$ (superficie lucida)")
    ax2.loglog(P, [r["T_vuoto_ossidato_K"] for r in sw], color=VERDE, lw=1.6,
               label=f"$\\epsilon={em2:.2f}$ (rugosa)")
    ax2.axhline(1358, color=ROSSO, ls="--", lw=1.0)
    ax2.annotate("fusione Cu 1358 K", (1.5e-7, 1430), fontsize=7, color=ROSSO)
    ax2.axhline(n.T0, color=GRIGIO, lw=0.8)
    ax2.axvline(1e-6, color=GRIGIO, ls=":", lw=1)
    ax2.annotate(f"1 $\\mu$W: {c0['T_vuoto_emiss_limpido_K']:.0f} K\n"
                 f"contro 1358 K", (1.2e-6, 300), fontsize=7, color=ROSSO)
    ax2.set_xlabel("potenza incidente (W)")
    ax2.set_ylabel("temperatura di regime in vuoto (K)")
    ax2.set_title("(b) vuoto: equilibrio radiativo")
    ax2.legend(loc="upper left", frameon=False)
    _pulisci(ax2)
    fig.tight_layout()
    _salva(fig, "fig03_termica.png")


# ---------------------------------------------------------------------------
# 4. BASELINE A PARI MASSA
# ---------------------------------------------------------------------------
def fig_baseline(res):
    c = res["baseline"]["casi"]
    nomi = [r["nome"].replace(" (V5)", "") for r in c]
    x = np.arange(len(c))
    fig, axs = plt.subplots(1, 2, figsize=(9.4, 3.8),
                            gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axs[0]
    b = np.array([r["beta_medio"] for r in c])
    lo = b - np.array([r["beta_lo"] for r in c])
    hi = np.array([r["beta_hi"] for r in c]) - b
    ax.bar(x, b, color=[BLU] + [GRIGIO] * 4 + [VERDE], width=0.6)
    ax.errorbar(x, b, yerr=[lo, hi], fmt="none", ecolor="black", capsize=3, lw=1)
    ax.set_yscale("log")
    ax.axhline(n.BETA_V5, color=ROSSO, ls="--", lw=1.0)
    ax.annotate("V5: $2{,}89\\times10^{-23}$\n(fuori intervallo di 12-22$\\times$)",
                (3.1, n.BETA_V5 * 3.2), fontsize=7, color=ROSSO)
    ax.set_xticks(x)
    ax.set_xticklabels(nomi, rotation=22, ha="right", fontsize=7.5)
    ax.set_ylabel("$\\beta$ (N m s)")
    ax.set_title("(a) attrito rotazionale a pari massa")
    _pulisci(ax)

    ax = axs[1]
    w = np.array([r["omega_ss_lo"] for r in c])
    wh = np.array([r["omega_ss_hi_ottica"] for r in c])
    ax.bar(x, w, color=[BLU] + [GRIGIO] * 4 + [VERDE], width=0.6)
    ax.errorbar(x, w, yerr=[np.zeros_like(w), wh - w], fmt="none", ecolor="black",
                capsize=3, lw=1)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(nomi, rotation=22, ha="right", fontsize=7.5)
    ax.set_ylabel("$\\omega_{ss}$ (rad/s)")
    ax.set_title("(b) velocita' di regime, stessa massa")
    ax.axhline(1.0, color=VERDE, ls=":", lw=1)
    _pulisci(ax)
    fig.tight_layout()
    _salva(fig, "fig04_baseline.png")


# ---------------------------------------------------------------------------
# 5. MATRICE DEI REGIMI
# ---------------------------------------------------------------------------
def fig_regimi(res):
    casi = res["regimi"]["casi"]
    L = ["1 nm", "240 nm", "1 $\\mu$m", "1 m"]
    mezzi = ["fluido", "vuoto"]
    M = np.zeros((2, 4))
    T = np.zeros((2, 4))
    valid = []
    for idx, c in enumerate(casi):
        i = 0 if c["mezzo"] == "fluido" else 1
        j = idx // 2
        M[i, j] = c["Pe_r"]
        T[i, j] = c["T_K"]
        valid.append((i, j, c["validita"]))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.2))
    for ax, dati, tit, fmt in ((ax1, M, "(a) $Pe_r=\\tau/k_BT$", "{:.1e}"),
                               (ax2, T, "(b) temperatura di regime (K)", "{:.0f}")):
        im = ax.imshow(np.log10(dati) if ax is ax1 else dati, cmap="viridis",
                       aspect="auto", origin="lower")
        ax.set_xticks(range(4))
        ax.set_xticklabels(L)
        ax.set_yticks(range(2))
        ax.set_yticklabels(mezzi)
        for i in range(2):
            for j in range(4):
                v = dati[i, j]
                ax.text(j, i, fmt.format(v), ha="center", va="center", fontsize=7.5,
                        bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.82))
        for (i, j, v) in valid:
            if not v.startswith("valido"):
                ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                       hatch="///", edgecolor="white", lw=0))
        ax.set_title(tit)
        ax.grid(False)
    fig.colorbar(im, ax=ax2, fraction=0.046, pad=0.03)
    fig.tight_layout()
    _salva(fig, "fig05_regimi.png")


# ---------------------------------------------------------------------------
# 6. MATERIALI
# ---------------------------------------------------------------------------
def fig_materiali(res):
    c = res["materiali"]["casi"]
    nomi = [r["nome"] for r in c]
    fig, axs = plt.subplots(1, 3, figsize=(9.6, 3.3))
    x = np.arange(len(c))
    ax = axs[0]
    ax.bar(x, [r["v_limite_m_s"] for r in c], color=[ARANCIO, GRIGIO, BLU], width=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(nomi, rotation=18, ha="right", fontsize=7.5)
    ax.set_ylabel("$v_{lim}=\\sqrt{\\sigma_y/\\rho}$ (m/s)")
    ax.set_title("(a) velocita' limite del materiale")
    _pulisci(ax)

    ax = axs[1]
    ax.bar(x, [r["sigma_abs"] / r["massa_kg"] for r in c],
           color=[ARANCIO, GRIGIO, BLU], width=0.6)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(nomi, rotation=18, ha="right", fontsize=7.5)
    ax.set_ylabel("$\\sigma_{ass}/m$ (m$^2$/kg)")
    ax.set_title("(b) assorbimento specifico")
    _pulisci(ax)

    ax = axs[2]
    for r, col in zip(c, [ARANCIO, GRIGIO, BLU]):
        ax.scatter(r["Pe_r"], r["margine_fusione"], s=70, color=col, label=r["nome"])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.axhline(1.0, color=ROSSO, ls="--", lw=1.0)
    ax.annotate("fusione", (1e-11, 1.15), fontsize=7, color=ROSSO)
    ax.set_xlabel("$Pe_r$ a 1 $\\mu$W")
    ax.set_ylabel("$T_m/T$ in vuoto")
    ax.set_title("(c) coppia contro margine termico")
    ax.legend(frameon=False, fontsize=7)
    _pulisci(ax)
    fig.tight_layout()
    _salva(fig, "fig06_materiali.png")


# ---------------------------------------------------------------------------
# 7. LEVE DI PROGETTO
# ---------------------------------------------------------------------------
def fig_leve(res):
    lv = [l for l in res["leve"]["leve"] if l["fattore_Pe"] > 0]
    nomi = [l["nome"] for l in lv]
    f = np.array([l["fattore_Pe"] for l in lv])
    fo = np.array([l["fattore_omega"] for l in lv])
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.6, 3.5), sharey=True)
    y = np.arange(len(lv))

    def pannello(ass, dati, tit, xlab, xlim):
        ass.barh(y, np.log10(dati), color=[VERDE if v >= 1 else ARANCIO for v in dati],
                 height=0.55)
        ass.axvline(0, color="black", lw=0.8)
        for yy, v in zip(y, dati):
            ass.annotate(f"$\\times${v:.2g}", (np.log10(v), yy),
                         xytext=(5 if v >= 1 else -5, 0), textcoords="offset points",
                         va="center", ha="left" if v >= 1 else "right", fontsize=7.5)
        ass.set_xlabel(xlab)
        ass.set_title(tit)
        ass.set_xlim(*xlim)
        _pulisci(ass)

    pannello(ax, f, "(a) effetto su $Pe_r$", "$Pe_r$ (decadi)", (-1.2, 1.4))
    pannello(ax2, fo, "(b) effetto su $\\omega_{ss}$", "$\\omega_{ss}$ (decadi)", (-1.2, 11.5))
    ax.set_yticks(y)
    ax.set_yticklabels(nomi, fontsize=8)
    fig.tight_layout()
    _salva(fig, "fig07_leve.png")


# ---------------------------------------------------------------------------
# 8. CONFRONTO CON LA SILICE A GHz
# ---------------------------------------------------------------------------
def fig_ghz(res):
    rot = [c for c in res["confronto"]["casi"] if c["tipo"] == "rottura"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.4, 3.7))
    r = np.logspace(-8.4, -6.9, 200)
    for sy, col in ((5e8, GRIGIO), (1e9, BLU), (3e9, VERDE), (5e9, OCRA), (8e9, ARANCIO)):
        ax1.loglog(r * 1e9, np.sqrt(sy / 2200.0) / r / (2 * math.pi), color=col, lw=1.3,
                   label=f"{sy/1e9:g} GPa")
    oss = res["confronto"]["registrazioni"]
    ax1.plot([o["r_nm"] for o in oss], [o["f_Hz"] for o in oss], "k*", ms=11,
             label="osservazioni", zorder=5)
    for o in oss:
        ax1.annotate(f"{o['f_Hz']/1e9:.0f} GHz", (o["r_nm"], o["f_Hz"]),
                     xytext=(4, -10), textcoords="offset points", fontsize=6.5)
    ax1.set_xlabel("raggio (nm)")
    ax1.set_ylabel("$f_{lim}$ (Hz)")
    ax1.set_title("(a) tetto meccanico contro frequenze osservate")
    ax1.legend(title="$\\sigma_y$ silice", frameon=False, fontsize=6.5, ncol=2)
    ax1.set_ylim(5e8, 2e10)
    _pulisci(ax1)

    pot = [c for c in res["confronto"]["casi"] if c["tipo"] == "potenza" and c["r_nm"] == 50]
    for delta, col, lab in ((0.02, ARANCIO, "$\\delta=0{,}02$"), (0.06, BLU, "$\\delta=0{,}06$")):
        sel = [c for c in pot if c["delta"] == delta]
        ax2.loglog([c["p_mbar"] for c in sel], [c["P_sorgente_necessaria_W"] for c in sel],
                   "o-", color=col, lw=1.4, ms=4, label=lab)
    ax2.axhspan(1e-2, 1e-1, color=GRIGIO, alpha=0.18)
    ax2.annotate("potenza tipica\ndegli esperimenti", (2e-10, 2.4e-2), fontsize=7, color=GRIGIO)
    ax2.set_xlabel("pressione residua (mbar)")
    ax2.set_ylabel("potenza di sorgente necessaria (W)")
    ax2.set_title("(b) coppia necessaria a sostenere 1 GHz ($r=50$ nm)")
    ax2.legend(frameon=False, fontsize=7)
    _pulisci(ax2)
    fig.tight_layout()
    _salva(fig, "fig08_ghz.png")


# ---------------------------------------------------------------------------
# 9. OAM
# ---------------------------------------------------------------------------
def fig_oam(res):
    c = [x for x in res["oam"]["casi"] if x["s"] == 1.0]
    l = [x["l"] for x in c]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.4))
    ax1.semilogy(l, [x["tau_ideale_Nm"] for x in c], "o--", color=GRIGIO, lw=1.2, ms=4,
                 label="$\\ell+s$ senza sovrapposizione")
    ax1.semilogy(l, [x["tau_reale_Nm"] for x in c], "o-", color=BLU, lw=1.6, ms=4,
                 label="con sovrapposizione $g(\\ell)$")
    ax1.axhline(n.KT, color=ROSSO, ls="--", lw=1)
    ax1.annotate("$k_BT$", (11, n.KT * 1.6), fontsize=7.5, color=ROSSO)
    ax1.set_xlabel("momento angolare orbitale $\\ell$ ($s=1$)")
    ax1.set_ylabel("$\\tau$ (N m)")
    ax1.set_title("(a) coppia contro $\\ell$")
    ax1.legend(frameon=False, fontsize=7)
    _pulisci(ax1)

    ax2.plot(l, [x["guadagno_netto"] for x in c], "o-", color=VERDE, lw=1.6, ms=4,
             label="guadagno netto")
    ax2.plot(l, [max(x["l"] + 1, 1) for x in c], ":", color=GRIGIO, lw=1.2,
             label="$\\ell+1$ (se $\\eta$ restasse uguale)")
    ax2.set_xlabel("$\\ell$")
    ax2.set_ylabel("$\\tau(\\ell)/\\tau(0)$")
    ax2.set_title("(b) il guadagno reale e' $\\sim\\sqrt{\\ell}$")
    ax2.legend(frameon=False, fontsize=7)
    _pulisci(ax2)
    fig.tight_layout()
    _salva(fig, "fig09_oam.png")


# ---------------------------------------------------------------------------
# 11-13. RENDER BLENDER
# ---------------------------------------------------------------------------
def fig_blender():
    bdir = os.path.normpath(os.path.join(BASE, "..", "blender"))
    if not os.path.isdir(bdir):
        print("  (nessun render Blender)")
        return
    # griglia delle sei geometrie a pari massa
    nomi = [("base_1.png", "Reticolo (V5)\n740 nm"), ("base_2.png", "Disco pieno\n200 nm"),
            ("base_3.png", "Anello\n586 nm"), ("base_4.png", "Sfera\n84 nm"),
            ("base_5.png", "Nanomanubrio\n200 nm"), ("base_6.png", "Silice\n134 nm")]
    disponibili = [(f, t) for f, t in nomi if os.path.exists(os.path.join(bdir, f))]
    if len(disponibili) == 6:
        fig, axs = plt.subplots(2, 3, figsize=(9.6, 6.5))
        for ax, (f, t) in zip(axs.ravel(), disponibili):
            ax.imshow(plt.imread(os.path.join(bdir, f)))
            ax.axis("off")
        fig.tight_layout(pad=0.2)
        _salva(fig, "fig11_baseline_render.png")
    # sequenza
    seq = [os.path.join(bdir, f"rotore_seq_{k}.png") for k in range(3)]
    if all(os.path.exists(p) for p in seq):
        fig, axs = plt.subplots(1, 3, figsize=(9.4, 3.3))
        for ax, p, ang in zip(axs, seq, (0, 120, 240)):
            ax.imshow(plt.imread(p))
            ax.set_title(f"{ang}$^\\circ$", fontsize=8)
            ax.axis("off")
        fig.tight_layout()
        _salva(fig, "fig12_sequenza.png")


# ---------------------------------------------------------------------------
# 10. CHIUSURA DEI RESIDUI
# ---------------------------------------------------------------------------
def fig_residui(res):
    a = res["residui"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.5))
    # (a) catena di arrotondamento
    ax1.set_xlim(0, 10)
    ax1.set_ylim(0, 10)
    ax1.axis("off")
    righe = [
        ("valore esatto", "$376/45 = 8{,}3555556$", BLU),
        ("arrotondato nella V5", "$8{,}3565$", ARANCIO),
        ("invertito", "$1/8{,}3565 = 0{,}1196676$", ARANCIO),
        ("confrontato con", "$45/376 = 0{,}1196809$", BLU),
        ("scarto apparente", "$-1{,}10\\times10^{-4}$", ROSSO),
        ("scarto vero della scena", "$-2{,}08\\times10^{-4}$", VERDE),
    ]
    for i, (lab, val, col) in enumerate(righe):
        y = 9 - i * 1.55
        ax1.text(0.2, y, lab, fontsize=8, va="center")
        ax1.text(4.6, y, val, fontsize=9, va="center", color=col, family="serif")
        if i < len(righe) - 1:
            ax1.annotate("", (2.3, y - 0.55), (2.3, y - 1.0),
                         arrowprops=dict(arrowstyle="-|>", color=GRIGIO, lw=0.9))
    ax1.set_title("(a) il residuo $I_2/I_1$: dove nasce $-1{,}10\\times10^{-4}$")

    # (b) beta: V5 contro intervallo fisico
    ax2.barh([0], [a["massa"]["massa_V5_kg"] * 0 + 1], color="none")
    et = ["V5\n(frazione d'area)", "corpo snello\n(limite inferiore)", "disco pieno\n(limite superiore)"]
    vals = [n.BETA_V5, n.beta_anello(n.R_OUT, n.WIRE_R), n.beta_disco(n.R_OUT)]
    cols = [ROSSO, BLU, BLU]
    ax2.barh(np.arange(3), vals, color=cols, height=0.5)
    ax2.set_xscale("log")
    ax2.set_yticks(np.arange(3))
    ax2.set_yticklabels(et, fontsize=7.5)
    ax2.set_xlabel("$\\beta$ (N m s)")
    ax2.set_xlim(1e-23, 2e-21)
    for i, v in enumerate(vals):
        ax2.annotate(f"{v:.2e}", (v, i), xytext=(5, 0), textcoords="offset points",
                     va="center", fontsize=7.5)
    ax2.set_title("(b) l'attrito rotazionale: la V5 e' fuori intervallo")
    _pulisci(ax2)
    fig.tight_layout()
    _salva(fig, "fig10_residui.png")


def main():
    print("FIGURE V6")
    res = {
        "residui": s.studio_residui(),
        "ottica": s.studio_ottica(),
        "termica": s.studio_termica(),
        "baseline": s.studio_baseline(),
        "regimi": s.studio_regimi(),
        "materiali": s.studio_materiali(),
        "leve": s.studio_leve(),
        "confronto": s.studio_confronto_ghz(),
        "oam": s.studio_oam(),
    }
    fig_modello(res)
    fig_ottica(res)
    fig_termica(res)
    fig_baseline(res)
    fig_regimi(res)
    fig_materiali(res)
    fig_leve(res)
    fig_ghz(res)
    fig_oam(res)
    fig_residui(res)
    fig_blender()
    print("fatto: figure in", FIG)


if __name__ == "__main__":
    main()
