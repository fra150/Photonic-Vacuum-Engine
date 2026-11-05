"""
MOTORE FOTONICO DEL VUOTO - V7.2
v72_figure.py - le figure della revisione: 1064 nm, costo energetico, firma.

Tutti i numeri vengono da v72_contrappeso.py o dal JSON che quello script
scrive.  Nessuna curva e' disegnata a mano.

Uso:
    python v72_figure.py
"""

import json
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

BASE = os.path.dirname(os.path.abspath(__file__))
V7 = os.path.normpath(os.path.join(BASE, "..", "..", "versione 7", "codice"))
for p in (V7, BASE):
    if p not in sys.path:
        sys.path.insert(0, p)

import v72_contrappeso as v72  # noqa: E402

DATI = os.path.normpath(os.path.join(BASE, "..", "dati"))
FIG = os.path.normpath(os.path.join(BASE, "..", "figure"))
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 190, "font.size": 9.5,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
    "axes.edgecolor": "#444444", "legend.frameon": False,
    "figure.facecolor": "white",
})

C = {532: "#c1121f", 650: "#e07a1f", 800: "#0077b6", 1064: "#588157"}


def _salva(fig, nome):
    p = os.path.join(FIG, nome)
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    print("  %-30s %6.1f kB" % (nome, os.path.getsize(p) / 1024.0))


def _curva(lam, z=2e-6, n=220):
    R = np.logspace(math.log10(0.15e-6), math.log10(20e-6), n)
    out = []
    for r in R:
        out.append(v72.ciclo(float(r), "vaterite", lam, z,
                             latch_fraction=0.5))
    return R, out


def fig_ciclo():
    R = np.logspace(math.log10(0.15e-6), math.log10(50e-6), 250)
    T720 = []
    for r in R:
        c = v72.ciclo(float(r), "vaterite", 800, 2e-6, latch_fraction=0.5)
        T720.append(c["T_720_s"] if c["omega"] > 0 else np.nan)

    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.6))
    ax = axes[0]
    ax.loglog(R * 1e6, T720, color="#333333", lw=1.9)
    ax.set_xlabel("R (µm)")
    ax.set_ylabel("$T_{720}$ (s), latch a 1/2 ciclo")
    ax.set_title("tempo del ciclo di 720 gradi (z = 2 µm)")
    ax.set_ylim(1e-1, 1e7)
    ax.text(0.05, 0.82, "$T_{720}=T_{360}^{\\rm no}+T_{360}^{\\rm cw}$\n"
                        "$=2\\pi(\\beta_{\\rm no}+\\beta_{\\rm cw})/(Pe_r kT)$\n"
                        "indipendente da $\\lambda$",
            transform=ax.transAxes, fontsize=8.8, color="#333333")

    ax = axes[1]
    for lam in (532, 800, 1064):
        R2, rows = _curva(lam, 2e-6, 260)
        dT = np.array([r["dT"] for r in rows])
        ax.loglog(R2 * 1e6, dT, color=C[lam], lw=1.7,
                  label="%d nm" % lam)
    ax.axhline(20, color="#c1121f", ls="--", lw=1.2,
               label="20 K (57 °C, danno)")
    ax.axhline(5, color="#e07a1f", ls=":", lw=1.2,
               label="5 K (42 °C, prudente)")
    ax.set_xlabel("R (µm)")
    ax.set_ylabel("$\\Delta T$ nel sangue (K)")
    ax.set_title("riscaldamento a $Pe_r=10$, z = 2 µm")
    ax.set_ylim(1e-1, 3e2)
    ax.legend(fontsize=8, loc="upper right")
    fig.suptitle("V7.2 - 1064 nm chiude il buco termico a 2.3 µm, ma non sposta "
                 "il bordo a basso R", y=1.02)
    _salva(fig, "fig_V72_ciclo.png")


def fig_profondita():
    z = np.logspace(-6, math.log10(2e-3), 240)
    fig, ax = plt.subplots(figsize=(7.8, 4.7))
    for lam, d in v72.OTTICA_SANGUE.items():
        T = np.exp(-(d["mu_a"] + d["mu_s"]) * z)
        ax.semilogx(z * 1e6, T, color=C[lam], lw=1.8,
                    label="%d nm (1/e = %.0f µm)"
                          % (lam, v72.profondita_balistica(lam) * 1e6))
    ax.axhline(math.exp(-1), color="#888888", ls=":", lw=1.0)
    ax.text(1.2, 0.42, "1/e", color="#666666", fontsize=8)
    ax.set_xlabel("profondità nel sangue z (µm)")
    ax.set_ylabel("frazione di luce coerente")
    ax.set_title("V7.2 - la luce che porta la coppia deve restare coerente:\n"
                 "a 1064 nm la profondità passa da 22 a 28 µm")
    ax.set_ylim(1e-4, 1.05)
    ax.legend()
    _salva(fig, "fig_V72_profondita.png")


def fig_contrappeso():
    d = json.load(open(os.path.join(DATI, "V72_sintesi.json"),
                       encoding="utf-8"))
    ce = d["costo_energetico"]
    R = [r["R_um"] for r in ce]
    e360 = [r["E_360_no_mJ"] for r in ce]
    e720 = [r["E_720_no_mJ"] for r in ce]
    ehalf = [r["E_720_half_mJ"] for r in ce]
    efull = [r["E_720_full_mJ"] for r in ce]
    x = np.arange(len(R))
    fig, ax = plt.subplots(figsize=(9.0, 4.8))
    ax.bar(x - 0.27, e360, 0.18, color="#c7e9c0", label="360°, senza")
    ax.bar(x - 0.09, e720, 0.18, color="#9ecae1", label="720°, senza")
    ax.bar(x + 0.09, ehalf, 0.18, color="#0077b6", label="720°, latch 1/2")
    ax.bar(x + 0.27, efull, 0.18, color="#08306b", label="720°, sempre innestato")
    for i, (a, b, c, e) in enumerate(zip(e360, e720, ehalf, efull)):
        ax.text(i - 0.09, b * 1.10, "×%.2f" % (b / a), ha="center",
                fontsize=7.5, color="#00456b")
        ax.text(i + 0.09, c * 1.10, "×%.2f" % (c / a), ha="center",
                fontsize=7.5, color="#00456b")
        ax.text(i + 0.27, e * 1.10, "×%.2f" % (e / a), ha="center",
                fontsize=7.5, color="#001f3f")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(["%.1f" % r for r in R])
    ax.set_xlabel("R (µm)")
    ax.set_ylabel("energia depositata per ciclo (mJ)")
    ax.set_title("V7.2 - il costo del contrappeso: 2.44× a metà ciclo, "
                 "3.88× nel limite conservativo\n"
                 "il fattore è rispetto allo stesso angolo; rispetto a un 360° "
                 "senza, il latch a metà costa 4.88×")
    ax.legend(fontsize=8, ncol=2)
    _salva(fig, "fig_V72_contrappeso.png")


def fig_firma():
    """La firma: omega modula a f_rot, lo stato del latch a f_rot/2."""
    t = np.linspace(0, 4.0, 4000)
    # omega modulato una volta per giro (latch che si innesta/sfila)
    m = 0.7423411936237939
    omega = 1.0 - m * 0.5 * (1.0 + np.cos(2.0 * math.pi * t))
    # stato del collegamento: A per 0-180, B per 180-360, A per 360-540...
    # periodo 720 gradi = 2 giri; e' dispari sotto 360.
    stato = np.where((np.mod(t, 2.0) < 1.0), 1.0, -1.0)

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.3))
    ax = axes[0]
    ax.plot(t, omega, color="#0077b6", lw=1.5, label="$\\omega/\\omega_0$ (latch)")
    ax.plot(t, stato, color="#c1121f", lw=1.3, ls="--",
            label="stato A/B del collegamento")
    ax.set_xlabel("giri (360° per unità)")
    ax.set_ylabel("segnale normalizzato")
    ax.set_title("ω ha periodo 360°; lo stato A/B ha periodo 720°")
    ax.set_ylim(-1.15, 1.55)
    ax.legend(fontsize=8, loc="upper right")
    ax.text(0.02, 0.04, "a 360° A→B, a 720° B→A:\nè questo che raddoppia il periodo",
            transform=ax.transAxes, fontsize=7.8, color="#444444")

    ax = axes[1]
    # FFT dei due segnali, su 8 giri, per mostrare le righe
    fs = len(t) / 4.0
    for sig, nome, col, ls in ((omega - omega.mean(), "$\\omega$", "#0077b6", "-"),
                               (stato, "stato A/B", "#c1121f", "--")):
        f = np.fft.rfftfreq(len(sig), d=1.0 / fs)
        A = np.abs(np.fft.rfft(sig)) / len(sig)
        ax.semilogy(f, A + 1e-6, color=col, lw=1.4, ls=ls, label=nome)
    ax.axvline(1.0, color="#0077b6", lw=0.8, ls=":")
    ax.axvline(0.5, color="#c1121f", lw=0.8, ls=":")
    ax.text(1.02, 1e-2, "$f_{\\rm rot}$", color="#0077b6", fontsize=8)
    ax.text(0.34, 1e-2, "$f_{\\rm rot}/2$", color="#c1121f", fontsize=8)
    ax.set_xlim(0, 3.2)
    ax.set_xlabel("frequenza (giri$^{-1}$)")
    ax.set_ylabel("ampiezza FFT")
    ax.set_title("spettro: il subarmonico a $f_{\\rm rot}/2$ è nello stato")
    ax.legend(fontsize=8, loc="upper right")
    fig.suptitle("V7.2 - la firma del 720° non è nell'orientazione: è nello "
                 "stato del latch", y=1.02)
    _salva(fig, "fig_V72_firma.png")


def fig_finestre():
    """Mappa R-z: verde = sicuro (dT<5 K), giallo = valido (5-20 K),
    rosso = non valido (>20 K o altro criterio)."""
    R = np.logspace(math.log10(0.15e-6), math.log10(50e-6), 90)
    Z = np.logspace(math.log10(2e-6), math.log10(50e-6), 24)
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.4), sharey=True)
    cmap = ListedColormap(["#c1121f", "#f2c14e", "#2a9d8f"])
    for ax, lam in zip(axes, (532, 800, 1064)):
        M = np.zeros((len(Z), len(R)))
        for i, z in enumerate(Z):
            for j, r in enumerate(R):
                c = v72.ciclo(float(r), "vaterite", lam, float(z),
                              latch_fraction=0.5)
                if not c["valido"]:
                    M[i, j] = 0
                elif not c["valido_sicuro"]:
                    M[i, j] = 1
                else:
                    M[i, j] = 2
        ax.pcolormesh(R * 1e6, Z * 1e6, M, cmap=cmap, vmin=0, vmax=2,
                      shading="auto")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("R (µm)")
        ax.set_title("%d nm" % lam)
        ax.set_ylim(2, 50)
    axes[0].set_ylabel("z (µm)")
    from matplotlib.patches import Patch
    handles = [Patch(color="#2a9d8f", label="sicuro, dT < 5 K"),
               Patch(color="#f2c14e", label="valido, 5-20 K"),
               Patch(color="#c1121f", label="non valido")]
    axes[-1].legend(handles=handles, fontsize=8, loc="lower left")
    fig.suptitle("V7.2 - finestra operativa: 1064 nm elimina il buco a 2-3 µm, "
                 "ma il bordo a basso R resta termico", y=1.02)
    _salva(fig, "fig_V72_finestre.png")


def main():
    print("V7.2 - figure")
    fig_profondita()
    fig_ciclo()
    fig_contrappeso()
    fig_firma()
    fig_finestre()
    print("Fatto: figure in", FIG)


if __name__ == "__main__":
    main()
