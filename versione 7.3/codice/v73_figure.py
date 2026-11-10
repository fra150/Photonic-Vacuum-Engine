# -*- coding: utf-8 -*-
"""
v73_figure.py - le figure della V7.3 (paper/it e paper/en, piu' copie PNG).

Legge dati/V73_griglia.json e dati/V73_sintesi.json.  Nessuna curva e'
disegnata a mano: quello che il testo chiama "Figura" esce da qui.

Uso:
    python codice/v73_figure.py              # italiano
    set V72_FIG_LANG=en && python codice/v73_figure.py
"""

import json
import math
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.lines import Line2D

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(BASE, ".."))
DATI = os.path.join(ROOT, "dati")
LANG = os.environ.get("V72_FIG_LANG", "it")
FIG = os.path.join(ROOT, "paper", "en" if LANG == "en" else "it", "figures")
FIG_PNG = os.path.join(ROOT, "figure", "en" if LANG == "en" else "it")
os.makedirs(FIG, exist_ok=True)
os.makedirs(FIG_PNG, exist_ok=True)


def TX(it, en):
    """Testo italiano o inglese secondo V72_FIG_LANG."""
    return en if LANG == "en" else it

with open(os.path.join(DATI, "V73_sintesi.json"), encoding="utf-8") as fh:
    S = json.load(fh)
with open(os.path.join(DATI, "V73_griglia.json"), encoding="utf-8") as fh:
    G = json.load(fh)
with open(os.path.join(DATI, "V73_taglio.json"), encoding="utf-8") as fh:
    TG = json.load(fh)
with open(os.path.join(DATI, "V73_termica.json"), encoding="utf-8") as fh:
    TH = json.load(fh)
with open(os.path.join(DATI, "V73_rumore.json"), encoding="utf-8") as fh:
    RU = json.load(fh)

CELLE = G["celle"]
LAM = (532, 800, 1064)
C = {532: "#b2182b", 800: "#2166ac", 1064: "#1b7837"}
Z0 = 2.0

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 220,
    "font.size": 9.5, "font.family": "DejaVu Sans",
    "axes.grid": True, "grid.alpha": 0.22, "grid.linewidth": 0.5,
    "axes.edgecolor": "#444444", "legend.frameon": False,
    "figure.facecolor": "white", "pdf.fonttype": 42, "ps.fonttype": 42,
})


def salva(fig, nome):
    p = os.path.join(FIG, nome)
    fig.savefig(p, bbox_inches="tight")
    png = os.path.join(FIG_PNG, os.path.splitext(nome)[0] + ".png")
    fig.savefig(png, bbox_inches="tight", dpi=220)
    plt.close(fig)
    print("  %-30s pdf %6.1f kB  png %6.1f kB"
          % (nome, os.path.getsize(p) / 1024.0,
             os.path.getsize(png) / 1024.0))


def celle(lam=None, z=Z0):
    out = [c for c in CELLE if abs(c["z_um"] - z) < 1e-9
           and (lam is None or c["lam_nm"] == lam)]
    return sorted(out, key=lambda c: c["R_um"])


def per_lam(lam, chiave, z=Z0):
    cc = celle(lam, z)
    x = np.array([c["R_um"] for c in cc])
    y = np.array([c[chiave] for c in cc], dtype=float)
    return x, y


# ---------------------------------------------------------------------------
def fig_profondita():
    z = np.logspace(-1, math.log10(2000), 400)  # micrometri
    fig, ax = plt.subplots(figsize=(7.4, 4.5))
    for lam in LAM:
        d = S["ottica"][str(lam)]
        mu_t = d["mu_a"] + d["mu_s"]
        T = np.exp(-mu_t * z * 1e-6)
        ell = 1e6 / mu_t
        ax.semilogx(z, T, color=C[lam], lw=1.9,
                    label="%d nm: $\\ell=%.1f\\,\\mu$m" % (lam, ell))
    ax.axhline(math.exp(-1), color="#777777", ls=":", lw=1.0)
    ax.axvline(2.0, color="#999999", ls="--", lw=0.9)
    ax.text(2.15, 0.055, "$z=2\\,\\mu$m", fontsize=8, color="#555555")
    ax.text(0.22, math.exp(-1) + 0.025, "$1/e$", fontsize=8, color="#555555")
    ax.set_xlim(0.1, 2000)
    ax.set_ylim(0, 1.03)
    ax.set_xlabel(TX("profondita' nel sangue $z$ ($\\mu$m)",
                    "depth in blood $z$ ($\\mu$m)"))
    ax.set_ylabel(TX("frazione di luce balistica $T_s(z)$",
                    "ballistic light fraction $T_s(z)$"))
    ax.set_title(TX("La luce balistica decade su pochi o decine di micrometri\n"
                   "(sangue intero, 150 g/L Hb, $g\\simeq0.98$)",
                   "Ballistic light decays over a few or tens of micrometres\n"
                   "(whole blood, 150 g/L Hb, $g\\simeq0.98$)"))
    ax.legend(loc="upper right")
    salva(fig, "fig_v73_profondita.pdf")


def fig_potenza():
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2), sharex=True)
    for lam in LAM:
        x, pt = per_lam(lam, "P_target_mW")
        _, pi = per_lam(lam, "P_inc_mW")
        axes[0].loglog(x, pt, color=C[lam], lw=1.8, label="%d nm" % lam)
        axes[1].loglog(x, pi, color=C[lam], lw=1.8, label="%d nm" % lam)
    axes[0].set_ylabel("$P_{\\rm target}$ (mW)")
    axes[1].set_ylabel("$P_{\\rm inc}$ (mW)")
    axes[0].set_title(TX("potenza al bersaglio, $Pe_r=10$",
                        "power at the target, $Pe_r=10$"))
    axes[1].set_title(TX("potenza incidente, $T_s(2\\,\\mu{\\rm m})$ inclusa",
                        "incident power, $T_s(2\\,\\mu{\\rm m})$ included"))
    for ax in axes:
        ax.set_xlabel("$R$ ($\\mu$m)")
        ax.axvline(0.5, color="#888888", ls=":", lw=1.0)
        ax.axvline(1.0, color="#888888", ls="--", lw=1.0)
        ax.legend(loc="upper left")
    axes[0].text(0.52, axes[0].get_ylim()[0] * 1.6, "CFL", fontsize=8,
                 color="#555555", rotation=90)
    axes[1].text(1.04, axes[1].get_ylim()[0] * 1.6,
                 TX("continuo", "continuum"), fontsize=8,
                 color="#555555", rotation=90)
    fig.suptitle(TX("A $Pe_r$ fissato la potenza non e' monotona: gli zeri di forma\n"
                   "la fanno oscillare; il sangue toglie il 7--19% a $z=2\\,\\mu$m",
                   "At fixed $Pe_r$ the power is not monotonic: shape zeros\n"
                   "make it oscillate; blood removes 7--19% at $z=2\\,\\mu$m"),
                 y=1.03)
    salva(fig, "fig_v73_potenza.pdf")


def fig_dt():
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    for lam in LAM:
        x, y = per_lam(lam, "dT")
        ax.loglog(x, y, color=C[lam], lw=1.9, label="%d nm" % lam)
    ax.axhline(20, color="#b2182b", ls="--", lw=1.1)
    ax.axhline(5, color="#e08214", ls=":", lw=1.2)
    ax.text(3.0, 23.0, TX("$20$ K: coagulazione (indicativa)",
                         "$20$ K: coagulation (indicative)"),
            color="#b2182b", fontsize=8)
    ax.text(3.0, 5.7, TX("$5$ K: criterio termico scelto, non prova di sicurezza",
                        "$5$ K: chosen thermal criterion, not a safety proof"),
            color="#e08214", fontsize=8)
    ax.set_xlim(0.15, 20)
    ax.set_ylim(0.02, 200)
    ax.set_xlabel("$R$ ($\\mu$m)")
    ax.set_ylabel("$\\Delta T$ (K)")
    ax.set_title(TX("Riscaldamento a $z=2\\,\\mu$m: a 532 nm domina il percorso,\n"
                   "a 800--1064 nm la particella; i due contributi si sommano",
                   "Heating at $z=2\\,\\mu$m: at 532 nm the path dominates,\n"
                   "at 800--1064 nm the particle; the two contributions add"))
    ax.legend(loc="lower left")
    salva(fig, "fig_v73_dt.pdf")


def fig_tempo():
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    x, y = per_lam(800, "T720_half")
    ax.loglog(x, y, color="#2166ac", lw=2.0,
              label=TX("$T_{720}$ (latch a meta')", "$T_{720}$ (half-cycle latch)"))
    _, yno = per_lam(800, "T720_no")
    ax.loglog(x, yno, color="#9ecae1", lw=1.2, ls="--",
              label=TX("$T_{720}$ senza contrappeso",
                      "$T_{720}$ without counterweight"))
    ax.axvspan(0.15, 0.5, color="#1b7837", alpha=0.10,
               label=TX("tutto nel CFL ($R\\leq0.5\\,\\mu$m)",
                       "all in CFL ($R\\leq0.5\\,\\mu$m)"))
    ax.axvspan(0.5, 0.8, color="#e08214", alpha=0.10,
               label=TX("CFL parziale, $\\tau/\\tau_B>1$",
                       "partial CFL, $\\tau/\\tau_B>1$"))
    ax.axvspan(0.8, 1.0, color="#b2182b", alpha=0.10,
               label=TX("Bingham blocca", "Bingham blocks"))
    ax.axvspan(1.0, 20.0, color="#999999", alpha=0.06,
               label=TX("continuo valido ($R\\geq1\\,\\mu$m)",
                       "continuum valid ($R\\geq1\\,\\mu$m)"))
    ax.axvline(0.5, color="#1b7837", lw=1.0)
    ax.axvline(0.8, color="#b2182b", lw=1.0)
    ax.axvline(1.0, color="#333333", lw=1.0)
    ax.set_xlim(0.15, 20)
    ax.set_ylim(0.02, 1e4)
    ax.set_xlabel("$R$ ($\\mu$m)")
    ax.set_ylabel("$T_{720}$ (s), 800 nm")
    ax.set_title(TX("Il tempo del ciclo cresce come $R^3$; le tre soglie non coincidono\n"
                   "($z=2\\,\\mu$m, $Pe_r=10$, $\\delta_{\\rm CFL}=3\\,\\mu$m)",
                   "Cycle time grows as $R^3$; the three thresholds do not coincide\n"
                   "($z=2\\,\\mu$m, $Pe_r=10$, $\\delta_{\\rm CFL}=3\\,\\mu$m)"))
    ax.legend(loc="upper left", fontsize=7.5, frameon=True, framealpha=0.9)
    salva(fig, "fig_v73_tempo.pdf")


def fig_mappa():
    R = np.array(sorted(set(round(c["R_um"], 6) for c in CELLE)))
    Z = np.array(sorted(set(round(c["z_um"], 6) for c in CELLE)))
    xe = np.exp(np.linspace(math.log(R[0] / 1.12), math.log(R[-1] * 1.12),
                            len(R) + 1))
    ze = np.concatenate([[Z[0] - (Z[1] - Z[0]) / 2],
                         (Z[:-1] + Z[1:]) / 2,
                         [Z[-1] + (Z[-1] - Z[-2]) / 2]])
    cmap = ListedColormap(["#b2182b", "#f4a582", "#fee090", "#a6d96a"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)

    fig, axes = plt.subplots(1, 3, figsize=(12.2, 3.9), sharey=True)
    for ax, lam in zip(axes, LAM):
        grid = np.zeros((len(Z), len(R)))
        for c in CELLE:
            if c["lam_nm"] != lam:
                continue
            i = int(np.argmin(np.abs(Z - c["z_um"])))
            j = int(np.argmin(np.abs(R - c["R_um"])))
            if not c["valido"]:
                v = 0
            elif not c["valido_sicuro"]:
                v = 1
            elif not c["valido_cfl"]:
                v = 2
            else:
                v = 3
            grid[i, j] = v
        ax.pcolormesh(xe, ze, grid, cmap=cmap, norm=norm,
                      shading="flat", edgecolors="none")
        ax.set_xscale("log")
        ax.set_xlim(0.15, 20)
        ax.set_ylim(0.5, 22.5)
        ax.axvline(1.0, color="#333333", ls="--", lw=0.9)
        ax.set_xlabel("$R$ ($\\mu$m)")
        ax.set_title("%d nm" % lam)
    axes[0].set_ylabel("$z$ ($\\mu$m)")
    handles = [Line2D([0], [0], marker="s", color="none",
                      markerfacecolor=col, markersize=9, label=lab)
               for col, lab in zip(["#b2182b", "#f4a582", "#fee090", "#a6d96a"],
                                   TX(["sopra 20 K", "5--20 K",
                                       "sotto il criterio 5 K", "sotto 5 K e CFL"],
                                      ["> 20 K", "5--20 K",
                                       "below the 5 K criterion", "below 5 K and CFL"]))]
    axes[2].legend(handles=handles, loc="lower right", fontsize=7.5,
                   frameon=True, framealpha=0.9)
    fig.suptitle(TX("Mappa $R$--$z$ del motore V7.2 (Bingham incluso); le soglie "
                   "termiche sono un criterio scelto, non una prova di sicurezza",
                   "$R$--$z$ map of the V7.2 engine (Bingham included); the "
                   "thermal thresholds are a chosen criterion, not a safety proof"),
                 y=1.03)
    salva(fig, "fig_v73_mappa.pdf")


def fig_cancello():
    """Shear ladder: the gate versus physiology, and the blood threshold."""
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.3))

    ax = axes[0]
    fis = TG["fisiologia"]
    ordine = ["venula 30 um", "capillare 6 um", "arteriola"]
    brevi = ["venula", "capillare", "arteriola"]
    for i, k in enumerate(ordine):
        lo, hi = fis[k]["gamma_range_s-1"]
        ax.hlines(i, lo, hi, color="#bdbdbd", lw=9, zorder=1)
        ax.plot(fis[k]["gamma_punto_s-1"], i, "o", color="#555555",
                ms=5, zorder=3)
        ax.text(hi * 1.12, i, TX("$\\times%.0f$" % TG["confronto_fisiologico"][brevi[i]]["rapporto"],
                                "$\\times%.0f$" % TG["confronto_fisiologico"][brevi[i]]["rapporto"]),
                va="center", fontsize=8, color="#555555")
    ax.axvline(TG["gamma_limite_s-1"], color="#2166ac", lw=2.0, zorder=2)
    gf = [x["gamma_s-1"] for x in TG["output_2_forza"]]
    ax.axvspan(min(gf), max(gf), color="#1b7837", alpha=0.15, zorder=0)
    ax.axvline(TG["output_3_jeffery"]["gamma_a_T720_s-1"], color="#777777",
               ls=":", lw=1.2, zorder=2)
    ax.text(TG["gamma_limite_s-1"] * 1.06, 2.42,
            TX("$\\dot\\gamma_{\\rm lim}=%.2f$" % TG["gamma_limite_s-1"],
               "$\\dot\\gamma_{\\rm lim}=%.2f$" % TG["gamma_limite_s-1"]),
            color="#2166ac", fontsize=8, rotation=90, va="top")
    ax.text(6.0, 2.30, TX("forza $Q=0.1$--$0.3$", "force $Q=0.1$--$0.3$"),
            color="#1b7837", fontsize=7.5)
    ax.set_xscale("log")
    ax.set_xlim(0.05, 4000)
    ax.set_ylim(-0.6, 2.6)
    ax.set_yticks(range(3))
    ax.set_yticklabels([TX("venula", "venule"), TX("capillare", "capillary"),
                        TX("arteriola", "arteriole")])
    ax.set_xlabel(TX("shear $\\dot\\gamma$ (s$^{-1}$)", "shear $\\dot\\gamma$ (s$^{-1}$)"))
    ax.set_title(TX("(a) il cancello e' sotto la fisiologia",
                    "(a) the gate lies below physiology"))
    ax.grid(True, axis="x", alpha=0.22, lw=0.5)

    ax = axes[1]
    sig = np.linspace(4.0, 21.0, 300)
    for mu, col in ((7.0, "#1b7837"), (8.0, "#e08214"), (9.0, "#b2182b")):
        My = 4.0 * math.pi * (0.5e-6) ** 3 * sig * 1e-3
        gam = (TG["input"]["tau_ott_Nm"] - My) / (
            TG["input"]["R_eff_factor"] * 8.0 * math.pi * mu * 1e-3 * (0.5e-6) ** 3)
        ax.plot(sig, np.maximum(gam, 0.0), color=col, lw=1.8,
                label="%s %d mPa$\\,$s" % (TX("$\\mu=$", "$\\mu=$"), mu))
    ax.axhline(TG["gamma_limite_s-1"], color="#2166ac", ls="--", lw=1.0)
    ax.text(4.2, TG["gamma_limite_s-1"] * 1.05,
            TX("plasma: 2.70 s$^{-1}$", "plasma: 2.70 s$^{-1}$"),
            fontsize=7.5, color="#2166ac")
    ax.set_yscale("log")
    ax.set_ylim(0.02, 5.0)
    ax.set_xlim(4, 21)
    ax.set_xlabel(TX("$\\sigma_y$ (mPa)", "$\\sigma_y$ (mPa)"))
    ax.set_ylabel(TX("$\\dot\\gamma_{\\rm lim}$ (s$^{-1}$)",
                     "$\\dot\\gamma_{\\rm lim}$ (s$^{-1}$)"))
    ax.set_title(TX("(b) soglia nel sangue: $0.10$--$0.41$ s$^{-1}$",
                    "(b) blood threshold: $0.10$--$0.41$ s$^{-1}$"))
    ax.legend(fontsize=7.5, loc="upper right")
    ax.grid(True, alpha=0.22, lw=0.5)

    fig.suptitle(TX("Il cancello non e' una finestra: e' un limite inferiore, "
                   "e la fisiologia e' 99--741 volte sopra",
                   "The gate is not a window: it is a lower bound, and "
                   "physiology is 99--741 times above"), y=1.02)
    salva(fig, "fig_v73_cancello.pdf")


def fig_termica():
    """Thermal bracket for the uniform medium: sink and mirror walls."""
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.3))
    lam = [532, 800, 1064]
    x = np.arange(len(lam))

    def th(env, l, R):
        for q in TH["casi"]:
            if (q["ambiente"] == env and q["lambda_nm"] == l
                    and abs(q["R_um"] - R) < 1e-6
                    and abs(q["z_um"] - 2.0) < 1e-9):
                return q
        raise KeyError((env, l, R))

    ax = axes[0]
    for R, col, off, mk in ((0.5, "#2166ac", -0.10, "o"),
                            (0.3, "#b2182b", 0.10, "s")):
        lo = [th("B", l, R)["dT_sink"] for l in lam]
        hi = [th("B", l, R)["dT_mirror"] for l in lam]
        mid = [th("B", l, R)["dT_0"] for l in lam]
        for i in range(len(lam)):
            ax.plot([x[i] + off, x[i] + off], [lo[i], hi[i]], color=col, lw=2)
            ax.plot(x[i] + off, lo[i], mk, color=col, ms=5)
            ax.plot(x[i] + off, hi[i], mk, color=col, ms=5, mfc="white")
            ax.plot(x[i] + off, mid[i], "|", color=col, ms=8)
    ax.axhline(TH["criterio"]["dT_K"], color="#333333", ls="--", lw=1.1)
    ax.text(2.42, 6.2, TX("criterio scelto 5 K", "chosen criterion 5 K"),
            fontsize=7.5, color="#333333", ha="right")
    ax.set_yscale("log")
    ax.set_ylim(0.15, 150)
    ax.set_xticks(x)
    ax.set_xticklabels(["532 nm", "800 nm", "1064 nm"])
    ax.set_ylabel("$\\Delta T$ (K)")
    ax.set_title(TX("(a) forcella parete: pozzo--specchio",
                    "(a) wall bracket: sink--mirror"))
    from matplotlib.lines import Line2D as _L
    ax.legend(handles=[
        _L([0], [0], color="#2166ac", marker="o", lw=2, label="$R=0.5\\,\\mu$m"),
        _L([0], [0], color="#b2182b", marker="s", lw=2, label="$R=0.3\\,\\mu$m")],
        fontsize=7.5, loc="upper right")

    ax = axes[1]
    q = th("B", 532, 0.5)
    lab = [TX("particella", "particle"), TX("percorso pre", "path pre"),
           TX("percorso post", "path post"), TX("diffuso", "diffuse")]
    val = [q["dT_part"], q["dT_pre"], q["dT_post"], q["dT_diff"]]
    y = np.arange(len(val))
    ax.barh(y, val, color=["#2166ac", "#e08214", "#e08214", "#1b7837"],
            height=0.55)
    for yi, v in zip(y, val):
        ax.text(v * 1.12, yi, "%.2f" % v, va="center", fontsize=8)
    ax.axvline(TH["criterio"]["dT_K"], color="#333333", ls="--", lw=1.1)
    ax.axvline(q["dT_sink"], color="#999999", ls=":", lw=1.2)
    ax.axvline(q["dT_mirror"], color="#999999", ls="-.", lw=1.2)
    from matplotlib.lines import Line2D as _L2
    ax.legend(handles=[
        _L2([0], [0], color="#333333", ls="--", lw=1.1,
            label=TX("criterio %.0f K" % TH["criterio"]["dT_K"],
                     "criterion %.0f K" % TH["criterio"]["dT_K"])),
        _L2([0], [0], color="#999999", ls=":", lw=1.2,
            label=TX("pozzo %.2f K" % q["dT_sink"],
                     "sink %.2f K" % q["dT_sink"])),
        _L2([0], [0], color="#999999", ls="-.", lw=1.2,
            label=TX("specchio %.2f K" % q["dT_mirror"],
                     "mirror %.2f K" % q["dT_mirror"]))],
        fontsize=7.2, loc="lower right")
    ax.set_xscale("log")
    ax.set_xlim(0.2, 40)
    ax.set_yticks(y)
    ax.set_yticklabels(lab)
    ax.invert_yaxis()
    ax.set_xlabel("$\\Delta T$ (K)")
    ax.set_title(TX("(b) 532 nm, $R=0.5\\,\\mu$m, sangue",
                    "(b) 532 nm, $R=0.5\\,\\mu$m, blood"))
    ax.grid(True, axis="x", alpha=0.22, lw=0.5)

    fig.suptitle(TX("532 nm supera il criterio anche nel limite ottimistico; "
                   "800 e 1064 nm restano sotto in entrambe le pareti",
                   "532 nm exceeds the criterion even in the optimistic bound; "
                   "800 and 1064 nm stay below for both walls"), y=1.02)
    salva(fig, "fig_v73_termica.pdf")


def fig_fuga():
    rr = S["fuga_bingham"]
    R = np.array([r["R_um"] for r in rr])
    P = np.array([r["P_inc_mW"] for r in rr])
    dT = np.array([r["dT"] for r in rr])
    Pe = np.array([r["Pe"] for r in rr])
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    ax = axes[0]
    ax.loglog(R, P, "o-", color="#2166ac", lw=1.8, ms=5)
    for x, y, p in zip(R, P, Pe):
        ax.annotate("$Pe=%.0f$" % p, (x, y), textcoords="offset points",
                    xytext=(4, 6), fontsize=7.5, color="#333333")
    ax.axhline(1.0, color="#b2182b", ls="--", lw=1.0)
    ax.text(1.05, 1.1, "1 W", color="#b2182b", fontsize=8)
    ax.set_xlabel("$R$ ($\\mu$m)")
    ax.set_ylabel("$P_{\\rm inc}$ (mW)")
    ax.set_title(TX("(a) potenza per superare $\\tau_B$",
                   "(a) power to overcome $\\tau_B$"))
    ax = axes[1]
    ax.loglog(R, dT, "s-", color="#b2182b", lw=1.8, ms=5)
    ax.axhline(20, color="#333333", ls="--", lw=1.0)
    ax.text(1.05, 24, "$20$ K", fontsize=8)
    ax.set_xlabel("$R$ ($\\mu$m)")
    ax.set_ylabel("$\\Delta T$ (K)")
    ax.set_title(TX("(b) riscaldamento nel punto di fuga",
                   "(b) heating at the escape point"))
    fig.suptitle(TX("Uscire dallo snervamento alzando $Pe$: la finestra si apre "
                   "a $R\\geq1\\,\\mu$m, ma il costo cresce",
                   "Escaping yield stress by raising $Pe$: the window opens "
                   "at $R\\geq1\\,\\mu$m, but the cost grows"), y=1.02)
    salva(fig, "fig_v73_fuga.pdf")


def fig_zeri():
    fig, ax = plt.subplots(figsize=(8.0, 4.6))
    dn = np.linspace(0.05, 0.12, 200)
    for lam in LAM:
        Rm = lam * 1e-3 / (2 * dn)  # micrometri
        ax.plot(dn, Rm, color=C[lam], lw=1.9, label="%d nm" % lam)
    ax.axvline(0.08, color="#333333", ls="--", lw=1.0)
    ax.axhspan(0.15, 1.0, color="#1b7837", alpha=0.10)
    ax.axhline(1.0, color="#333333", ls=":", lw=1.0)
    ax.text(0.051, 0.55, TX("finestra del continuo\n$R\\leq1\\,\\mu$m",
                           "continuum window\n$R\\leq1\\,\\mu$m"),
            fontsize=8, color="#1b7837")
    ax.set_xlim(0.05, 0.12)
    ax.set_ylim(0, 14)
    ax.set_xlabel(TX("$\\Delta n$ effettivo della vaterite",
                    "effective $\\Delta n$ of vaterite"))
    ax.set_ylabel("$R_m=m\\lambda/(2\\Delta n)$ ($\\mu$m), $m=1$")
    ax.set_title(TX("Gli zeri di forma a $\\Delta n=0.08$ cadono a $3.3$--$6.6\\,\\mu$m:\n"
                   "fuori dalla finestra in cui il rotore esiste",
                   "Shape zeros at $\\Delta n=0.08$ fall at $3.3$--$6.6\\,\\mu$m:\n"
                   "outside the window where the rotor exists"))
    ax.legend(loc="upper right")
    salva(fig, "fig_v73_zeri.pdf")


def main():
    print("Figure V7.3")
    fig_profondita()
    fig_potenza()
    fig_dt()
    fig_tempo()
    fig_mappa()
    fig_cancello()
    fig_termica()
    fig_fuga()
    fig_zeri()
    print("Fatto: figure in", FIG)


if __name__ == "__main__":
    main()
