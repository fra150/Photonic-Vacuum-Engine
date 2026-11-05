"""
MOTORE FOTONICO DEL VUOTO - V7.1
v71_figure.py - le figure dell'appendice 720 gradi / contrappeso / sangue.

Tutti i numeri vengono da v71_contrappeso.py: nessuna curva e' disegnata a mano.

Uso:
    python v71_figure.py
"""

import json
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
V7 = os.path.normpath(os.path.join(BASE, "..", "..", "versione 7", "codice"))
for p in (V7, BASE):
    if p not in sys.path:
        sys.path.insert(0, p)

import v71_contrappeso as v71  # noqa: E402

DATI = os.path.normpath(os.path.join(BASE, "..", "dati"))
FIG = os.path.normpath(os.path.join(BASE, "..", "figure"))
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 190, "font.size": 9.5,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
    "axes.edgecolor": "#444444", "legend.frameon": False,
    "figure.facecolor": "white",
})

C = {532: "#c1121f", 800: "#0077b6", 650: "#e07a1f", 1064: "#588157"}


def _salva(fig, nome):
    p = os.path.join(FIG, nome)
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    print("  %-28s %6.1f kB" % (nome, os.path.getsize(p) / 1024.0))


def fig_ciclo():
    R = np.array(v71.R_GRID)
    T720, dTp, dTf, val = [], [], [], []
    for r in R:
        c = v71.ciclo(float(r), "vaterite", 800, 2e-6)
        T720.append(c["T_720_s"] if c["omega"] > 0 else np.nan)
        dTp.append(c["dT_percorso"])
        dTf.append(c["dT_particella"])
        val.append(c["valido"])
    val = np.array(val)
    r_ok = R[val] * 1e6

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.4))
    ax = axes[0]
    if len(r_ok):
        ax.axvspan(r_ok.min(), r_ok.max(), color="#0077b6", alpha=0.10,
                   label="finestra valida (800 nm)")
    ax.loglog(R * 1e6, T720, color="#0077b6", lw=1.9)
    ax.set_xlabel("R (µm)")
    ax.set_ylabel("$T_{720}$ (s)")
    ax.set_title("tempo del ciclo di 720 gradi (z = 2 µm)")
    ax.set_ylim(1e-2, 1e5)
    ax.text(0.05, 0.80, "$T_{720} = 4\\pi\\beta/(Pe_r\\,kT) \\propto R^3$\n"
                        "indipendente da $\\lambda$", transform=ax.transAxes,
            fontsize=9, color="#333333")
    ax.legend(loc="lower right", fontsize=8)

    ax = axes[1]
    ax.loglog(R * 1e6, dTp, color="#0077b6", lw=1.9,
              label="$\\Delta T$ nel sangue (segue $P_{target}$, quindi $\\delta$)")
    ax.loglog(R * 1e6, dTf, color="#e07a1f", lw=1.0, ls="--",
              label="$\\Delta T$ della particella")
    ax.axhline(20, color="#c1121f", ls="--", lw=1.2, label="limite 20 K")
    ax.set_xlabel("R (µm)")
    ax.set_ylabel("$\\Delta T$ (K)")
    ax.set_title("riscaldamento a $Pe_r = 10$: due contributi, uno decide")
    ax.legend(fontsize=8)
    fig.suptitle("V7.1 - dentro il sangue il tempo del ciclo scende con $R^3$, il "
                 "riscaldamento sale: la finestra e' stretta", y=1.02)
    _salva(fig, "fig_V71_ciclo.png")


def fig_profondita():
    z = np.logspace(-6, math.log10(2e-3), 200)
    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    for lam, d in v71.OTTICA_SANGUE.items():
        T = np.exp(-(d["mu_a"] + d["mu_s"]) * z)
        ax.semilogx(z * 1e6, T, color=C[lam], lw=1.8,
                    label="%d nm (1/e = %.0f µm)" % (lam, v71.profondita_balistica(lam) * 1e6))
    ax.axhline(math.exp(-1), color="#888888", ls=":", lw=1.0)
    ax.text(1.2, 0.42, "1/e", color="#666666", fontsize=8)
    ax.set_xlabel("profondità nel sangue z (µm)")
    ax.set_ylabel("frazione di luce coerente")
    ax.set_title("V7.1 - la luce che porta la coppia deve restare coerente:\n"
                 "nel sangue intero sono decine di micrometri, non millimetri")
    ax.set_ylim(1e-4, 1.05)
    ax.legend()
    _salva(fig, "fig_V71_profondita.png")


def fig_contrappeso():
    R = np.array([0.2e-6, 0.3e-6, 0.5e-6, 1e-6, 2e-6, 5e-6])
    senza, con = [], []
    for r in R:
        senza.append(v71.ciclo(float(r), "vaterite", 800, 2e-6, innestato=False)["T_720_s"])
        con.append(v71.ciclo(float(r), "vaterite", 800, 2e-6, innestato=True)["T_720_s"])
    x = np.arange(len(R))
    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    ax.bar(x - 0.19, senza, 0.38, color="#9ecae1", label="senza contrappeso")
    ax.bar(x + 0.19, con, 0.38, color="#0077b6", label="con contrappeso (V/2 a 1.5 R)")
    for i, (a, b) in enumerate(zip(senza, con)):
        ax.text(i + 0.19, b * 1.06, "×%.2f" % (b / a), ha="center", fontsize=8,
                color="#00456b")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(["%.1f" % (r * 1e6) for r in R])
    ax.set_xlabel("R (µm)")
    ax.set_ylabel("$T_{720}$ (s), 800 nm, z = 2 µm")
    c = v71.contrappeso(1e-6, "vaterite")
    ax.set_title("V7.1 - quanto costa il contrappeso: I × %.2f, β × %.2f, "
                 "quindi il ciclo × %.2f" % (c["costo_inerzia"], c["costo_attrito"],
                                             con[2] / senza[2]))
    ax.legend()
    _salva(fig, "fig_V71_contrappeso.png")


def fig_modulazione():
    """La firma misurabile: il subarmonico a omega/2.

    Schema, non simulazione: l'ampiezza della modulazione e' quella calcolata
    (1 - I_rot/I_tot = 76%), la forma d'onda e' un'ipotesi sinusoidale con
    periodo un giro.  L'integrale e' in forma chiusa.
    """
    rif = json.load(open(os.path.join(DATI, "V71_sintesi.json"), encoding="utf-8"))[
        "cella_riferimento_3d"]
    mod = rif["modulazione"]
    t = np.linspace(0, 4.0, 2000)

    def angolo(m):
        # omega = w0 (1 - m (0.5 + 0.5 cos 2 pi t))
        return t - m * (0.5 * t + np.sin(2 * math.pi * t) / (4 * math.pi))

    fig, ax = plt.subplots(figsize=(7.8, 4.2))
    ax.plot(t, angolo(0.0), color="#588157", lw=1.2, ls=":",
            label="senza contrappeso (ω costante)")
    ax.plot(t, angolo(mod), color="#0077b6", lw=1.8,
            label="con latch: modulazione %.0f%% a ogni giro" % (mod * 100))
    ax.plot(t, (1 - mod) * t, color="#c1121f", lw=1.0, ls="--",
            label="se il contrappeso restasse innestato")
    ax.set_xlabel("giri (360 gradi per unità)")
    ax.set_ylabel("angolo cumulato (giri)")
    ax.set_title("V7.1 - la firma del latch: la velocità cala del %.0f%% a ogni giro,\n"
                 "cioè un subarmonico a $\\omega/2$ nello spettro (predizione falsificabile)"
                 % (mod * 100))
    ax.legend()
    _salva(fig, "fig_V71_firma.png")


def main():
    print("V7.1 - figure")
    fig_profondita()
    fig_ciclo()
    fig_contrappeso()
    fig_modulazione()
    print("Fatto: figure in", FIG)


if __name__ == "__main__":
    main()
