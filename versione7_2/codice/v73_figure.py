# -*- coding: utf-8 -*-
"""
v73_figure.py - le figure della V7.2 corretta (cartella paper/V73).

Legge dati/V73_griglia.json e dati/V73_sintesi.json.  Nessuna curva e'
disegnata a mano: quello che il testo chiama "Figura" esce da qui.

Uso:
    python codice/v73_figure.py
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
FIG = os.path.join(ROOT, "paper", "V73",
                   "figures_v73_en" if LANG == "en" else "figures_v73")
os.makedirs(FIG, exist_ok=True)


def TX(it, en):
    """Testo italiano o inglese secondo V72_FIG_LANG."""
    return en if LANG == "en" else it

with open(os.path.join(DATI, "V73_sintesi.json"), encoding="utf-8") as fh:
    S = json.load(fh)
with open(os.path.join(DATI, "V73_griglia.json"), encoding="utf-8") as fh:
    G = json.load(fh)

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
    plt.close(fig)
    print("  %-30s %6.1f kB" % (nome, os.path.getsize(p) / 1024.0))


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
    ax.text(3.0, 23.0, TX("$20$ K: coagulazione", "$20$ K: coagulation"),
            color="#b2182b", fontsize=8)
    ax.text(3.0, 5.7, TX("$5$ K: soglia prudente ($42^\\circ$C)",
                        "$5$ K: prudent threshold ($42^\\circ$C)"),
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
                                   TX(["non valido", "valido ($5$--$20$ K)",
                                      "sicuro ($<5$ K)", "sicuro e CFL"],
                                     ["invalid", "valid ($5$--$20$ K)",
                                      "safe ($<5$ K)", "safe and CFL"]))]
    axes[2].legend(handles=handles, loc="lower right", fontsize=7.5,
                   frameon=True, framealpha=0.9)
    fig.suptitle(TX("Mappa $R$--$z$ della validita' (Bingham incluso); la linea "
                   "tratteggiata e' il limite del continuo $R=1\\,\\mu$m",
                   "$R$--$z$ validity map (Bingham included); the dashed line "
                   "is the continuum limit $R=1\\,\\mu$m"), y=1.03)
    salva(fig, "fig_v73_mappa.pdf")


def fig_rumore():
    c = [x for x in S["confronto_lambda"]
         if abs(x["R_um"] - 0.5) < 1e-9 and x["lam_nm"] == 800][0]
    Tno, Tcw = c["T360_no"], c["T360_cw"]
    cv = c["cv_durata"]
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))

    ax = axes[0]
    for Tv, col, lab in ((Tno, "#1b7837",
                          TX("giro senza contrappeso",
                            "turn without counterweight")),
                         (Tcw, "#b2182b",
                          TX("giro con contrappeso",
                            "turn with counterweight"))):
        mu = math.log(Tv) - 0.5 * math.log(1 + cv ** 2)
        sig = math.sqrt(math.log(1 + cv ** 2))
        tt = np.linspace(Tv * 0.35, Tv * 2.8, 500)
        pdf = np.exp(-(np.log(tt) - mu) ** 2 / (2 * sig ** 2)) / \
            (tt * sig * math.sqrt(2 * math.pi))
        ax.plot(tt, pdf, color=col, lw=1.9, label=lab)
    thr = math.sqrt(Tno * Tcw)
    ax.axvline(thr, color="#333333", ls="--", lw=1.0)
    ax.set_xscale("log")
    ax.set_xlim(0.3, 8.0)
    ax.set_xlabel(TX("durata di un giro (s)", "duration of one turn (s)"))
    ax.set_ylabel(TX("densita' di probabilita'", "probability density"))
    ax.set_title(TX("(a) durate log-normali, $R=0.5\\,\\mu$m\n"
                   "$CV=18\\%$, separazione $3.88\\times$",
                   "(a) log-normal durations, $R=0.5\\,\\mu$m\n"
                   "$CV=18\\%$, separation $3.88\\times$"))
    ax.legend(fontsize=8, loc="upper right")

    ax = axes[1]
    N = np.arange(1, 9)
    z_no = c["z_no"]
    z_cw = c["z_cw"]
    eno = 0.5 * np.array([math.erfc(z_no * math.sqrt(n) / math.sqrt(2))
                          for n in N])
    ecw = 0.5 * np.array([math.erfc(z_cw * math.sqrt(n) / math.sqrt(2))
                          for n in N])
    etot = eno + ecw
    ax.semilogy(N, etot, "o-", color="#b2182b", lw=1.8, ms=5)
    ax.axhline(1e-6, color="#777777", ls=":", lw=1.0)
    ax.text(2.05, 1.4e-6, "$10^{-6}$", fontsize=8, color="#666666")
    ax.annotate("N=2: $%.1f\\times10^{-5}$\nN=3: $%.1f\\times10^{-7}$"
                % (etot[1] * 1e5, etot[2] * 1e7),
                xy=(2, etot[1]), xytext=(2.6, 3e-9), fontsize=8,
                arrowprops=dict(arrowstyle="->", lw=0.7, color="#555555"))
    ax.set_ylim(1e-10, 1e-1)
    ax.set_xlabel(TX("numero di cicli osservati (media delle durate)",
                    "number of observed cycles (average duration)"))
    ax.set_ylabel(TX("errore totale di classificazione",
                    "total classification error"))
    ax.set_title(TX("(b) bastano 2--3 cicli: la durata classifica\n"
                   "senza ricostruire $\\theta$ (modello gaussiano ideale)",
                   "(b) 2--3 cycles suffice: duration classifies\n"
                   "without reconstructing $\\theta$ (ideal Gaussian model)"))
    ax.set_xticks(N)
    fig.suptitle(TX("Il rumore non cancella il latch: cambia la durata dei due giri",
                   "Noise does not erase the latch: it changes the duration of the two turns"),
                 y=1.02)
    salva(fig, "fig_v73_rumore.pdf")


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
    print("Figure V7.2 (cartella V73)")
    fig_profondita()
    fig_potenza()
    fig_dt()
    fig_tempo()
    fig_mappa()
    fig_rumore()
    fig_fuga()
    fig_zeri()
    print("Fatto: figure in", FIG)


if __name__ == "__main__":
    main()
