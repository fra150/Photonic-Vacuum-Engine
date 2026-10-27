"""
MOTORE FOTONICO V5 - GRAFICI DEL MODULO SPAZIALE
Autore: Dr. Francesco Bulla / Sistema GNRX

Legge SP1-SP5 e disegna quattro figure:

  fig_SP1_scala_determinismo.png   omega di regime e Pe_r alle tre scale
  fig_SP2_terra_vs_spazio.png      forze a terra e rapporto spinta/peso
  fig_SP3_materiali_1m.png         i tre materiali del caso a 1 m
  fig_SP4_rapporto_45_376.png      regime accoppiato e invarianza del rapporto

Uso:
    python v5_grafici_spazio.py
"""

import json
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
DATI = os.path.join(BASE, "dati")
FIG = os.path.join(BASE, "figure")
os.makedirs(FIG, exist_ok=True)

sys.path.insert(0, BASE)
import v5_spazio as s  # noqa: E402
import v5_nucleo as n  # noqa: E402

BLU = "#1f4e79"
ARANCIO = "#c55a11"
VERDE = "#2e7d32"
ROSSO = "#b00020"
GRIGIO = "#5a5a5a"

plt.rcParams.update({
    "figure.dpi": 160,
    "savefig.dpi": 160,
    "font.size": 9.5,
    "axes.titlesize": 10.5,
    "axes.labelsize": 9.5,
    "axes.edgecolor": "#404040",
    "axes.linewidth": 0.8,
    "axes.grid": True,
    "grid.color": "#d9d9d9",
    "grid.linewidth": 0.6,
    "legend.frameon": False,
    "figure.facecolor": "white",
})


def carica(nome):
    with open(os.path.join(DATI, nome), encoding="utf-8") as fh:
        return json.load(fh)


def pulisci(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def salva(fig, nome):
    percorso = os.path.join(FIG, nome)
    fig.tight_layout()
    fig.savefig(percorso, bbox_inches="tight")
    plt.close(fig)
    print("  %-42s %6.1f kB" % (nome, os.path.getsize(percorso) / 1024.0))


# ---------------------------------------------------------------------------
# Sweep continuo in k, con le stesse funzioni dei dati
# ---------------------------------------------------------------------------
K_SWEEP = np.logspace(-3.0, 6.0, 90)


def sweep():
    """omega e Pe_r per Spin e Contro-spin su tutto l'intervallo di scala."""
    omega_spin, omega_cont, pe_dich, pe_scatt = [], [], [], []
    for k in K_SWEEP:
        b1 = n.beta_telaio(3.70 * n.UA * k)
        b2 = n.beta_telaio(3.275 * n.UA * k)
        Fd = (n.P_LASER / n.C_LIGHT) * (n.ETA_S_DECL + n.ETA_C_DECL *
                                       math.cos(n.ALPHA_DECL)) * min(1.0, k ** 2)
        Fs = n.scattering_globale(k)["F_z_N"]
        omega_spin.append(Fd * 3.70 * n.UA * k / b1)
        omega_cont.append(Fd * 3.275 * n.UA * k / b2)
        pe_dich.append(Fd * 3.70 * n.UA * k / n.KT)
        pe_scatt.append(Fs * 3.70 * n.UA * k / n.KT)
    return (np.array(omega_spin), np.array(omega_cont),
            np.array(pe_dich), np.array(pe_scatt))


def ingombro(k):
    return 1.9e-6 * k


# ---------------------------------------------------------------------------
# SP1 - scala e determinismo
# ---------------------------------------------------------------------------
def fig_sp1():
    omega_spin, omega_cont, pe_dich, pe_scatt = sweep()
    ing = ingombro(K_SWEEP)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.8, 4.0))

    a1.loglog(ing, omega_spin, color=BLU, lw=2.2, label="Spin (3.70 u.a.)")
    a1.loglog(ing, omega_cont, color=ARANCIO, lw=2.2, ls="--",
              label="Contro-spin (3.275 u.a.)")
    for k, etichetta in ((1e-3, "1 nm"), (1.0, "10$^{-18}$ kg (V4)"), (1e6, "1 m")):
        a1.plot(ingombro(k), 45.0 if etichetta != "1 m" else 4.5e-11,
                "o", ms=6, mfc="white", mec=BLU, mew=1.4, zorder=5)
    a1.annotate("1 nm", xy=(ingombro(1e-3), 45.0), xytext=(1.05e-9, 0.08),
                fontsize=8.5, color=GRIGIO, ha="center",
                arrowprops=dict(arrowstyle="->", color=GRIGIO, lw=0.9))
    a1.annotate("V4 · 10$^{-18}$ kg", xy=(ingombro(1.0), 45.0),
                xytext=(4e-7, 260), fontsize=8.5, color=GRIGIO,
                arrowprops=dict(arrowstyle="->", color=GRIGIO, lw=0.9))
    a1.annotate("1 m", xy=(ingombro(1e6), 4.5e-11), xytext=(0.30, 1.2e-9),
                fontsize=8.5, color=GRIGIO,
                arrowprops=dict(arrowstyle="->", color=GRIGIO, lw=0.9))
    a1.axvspan(ing[0], ingombro(1.0), color=BLU, alpha=0.05, lw=0)
    a1.text(2.2e-9, 2.2e-3, "forza ∝ k$^2$,  ω costante", fontsize=8, color=GRIGIO)
    a1.text(3e-3, 2.5e-9, "forza satura,  ω ∝ k$^{-2}$", fontsize=8, color=GRIGIO)
    a1.set_xlabel("ingombro del dispositivo (m)")
    a1.set_ylabel("$\\omega$ di regime (rad s$^{-1}$)")
    a1.set_title("(a) Velocità angolare in microgravità (g = 0)")
    a1.legend(loc="lower left")
    a1.set_xlim(ing[0], ing[-1])
    a1.set_ylim(1e-12, 1e4)

    a2.loglog(ing, pe_dich, color=BLU, lw=2.2, label="coppia dichiarata (V4)")
    a2.loglog(ing, pe_scatt, color=ROSSO, lw=2.2, ls="--",
              label="coppia da scattering")
    a2.axhline(1.0, color=GRIGIO, lw=1.1, ls=":")
    a2.fill_between(ing, 1e-20, 1.0, color=ROSSO, alpha=0.05, lw=0)
    a2.text(1.5e-8, 2e6, "regime deterministico", fontsize=8.5, color=BLU)
    a2.text(1.5e-8, 1e-8, "regime termico", fontsize=8.5, color=ROSSO)
    for k, valore, testo, colore in ((1.0, 0.4385, "0.44", BLU),
                                     (1.0, 2.665e-3, "2.7×10$^{-3}$", ROSSO)):
        a2.plot(ingombro(k), valore, "o", ms=6, mfc="white", mec=colore,
                mew=1.4, zorder=5)
        a2.text(ingombro(k) * 1.6, valore * 1.15, testo, fontsize=8.5, color=colore)
    a2.set_xlabel("ingombro del dispositivo (m)")
    a2.set_ylabel("$Pe_r = \\omega/D_r = \\tau/kT$")
    a2.set_title("(b) Numero di Péclet rotazionale")
    a2.legend(loc="lower right")
    a2.set_xlim(ing[0], ing[-1])
    a2.set_ylim(1e-20, 1e12)

    salva(fig, "fig_SP1_scala_determinismo.png")


# ---------------------------------------------------------------------------
# SP2 - terra contro spazio
# ---------------------------------------------------------------------------
def fig_sp2():
    d = carica("SP2_terra_vs_spazio.json")
    casi = carica("SP1_spazio_scale.json")["casi"]["massa_10^-18_kg"]

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.8, 4.0))

    voci = ["peso", "peso\napparente", "spinta\nscattering", "spinta\ndichiarata"]
    chiavi = [("terra", "peso_N"), ("terra", "peso_apparente_N"),
              ("terra", "spinta_peso_scattering"), ("terra", "spinta_peso_dichiarata")]
    x = np.arange(len(voci))
    larg = 0.36
    for i, (rot, colore, off) in enumerate((("spin", BLU, -larg / 2),
                                            ("contro", ARANCIO, larg / 2))):
        r = casi[rot]
        massa = r["inerzia"]["massa_kg"]
        peso = r["terra"]["peso_N"]
        valori = [peso, r["terra"]["peso_apparente_N"],
                  r["terra"]["spinta_peso_scattering"] * peso,
                  r["terra"]["spinta_peso_dichiarata"] * peso]
        a1.bar(x + off, valori, larg, color=colore, alpha=0.85,
               label="Spin" if i == 0 else "Contro-spin", edgecolor="white", lw=0.6)
        for xx, vv in zip(x + off, valori):
            a1.text(xx, vv * 1.25, "%.2e" % vv, ha="center", fontsize=7,
                    rotation=90, color=colore)

    a1.set_yscale("log")
    a1.set_xticks(x)
    a1.set_xticklabels(voci)
    a1.set_ylabel("forza (N)")
    a1.set_ylim(1e-18, 1e-12)
    a1.set_title("(a) Bilancio delle forze sul dispositivo V4, a terra")
    a1.legend(loc="upper left", ncols=2)
    a1.text(0.02, 0.86, "in orbita: peso ≡ 0, restano solo $\\tau$ e $\\beta$",
            transform=a1.transAxes, fontsize=8, color=GRIGIO)

    omega_spin, omega_cont, pe_dich, pe_scatt = sweep()
    ing = ingombro(K_SWEEP)
    massa_spin = 8960.0 * 0.1391 * (n.UA ** 3) * K_SWEEP ** 3
    peso_spin = massa_spin * s.G_TERRA
    F_dich = (n.P_LASER / n.C_LIGHT) * (n.ETA_S_DECL + n.ETA_C_DECL *
                                        math.cos(n.ALPHA_DECL)) * np.minimum(1.0, K_SWEEP ** 2)
    F_scatt = np.array([n.scattering_globale(k)["F_z_N"] for k in K_SWEEP])
    r_dich = F_dich / peso_spin
    r_scatt = F_scatt / peso_spin

    a2.loglog(ing, r_dich, color=BLU, lw=2.2, label="forza dichiarata")
    a2.loglog(ing, r_scatt, color=ROSSO, lw=2.2, ls="--", label="forza di scattering")
    a2.axhline(1.0, color=GRIGIO, lw=1.1, ls=":")
    a2.text(2e-9, 1.4, "spinta = peso", fontsize=8.5, color=GRIGIO)

    for rapporto, colore, etichetta in ((r_dich, BLU, "dichiarata"),
                                        (r_scatt, ROSSO, "scattering")):
        idx = int(np.argmin(np.abs(np.log(rapporto))))
        k_cross = K_SWEEP[idx]
        a2.plot(ingombro(k_cross), 1.0, "o", ms=6, mfc="white", mec=colore,
                mew=1.4, zorder=6)
        a2.annotate("%.1f µm" % (ingombro(k_cross) * 1e6),
                    xy=(ingombro(k_cross), 1.0),
                    xytext=(ingombro(k_cross) * (12 if colore == ROSSO else 0.05),
                            1.0 * (10 ** (-3 if colore == ROSSO else 2))),
                    fontsize=8, color=colore,
                    arrowprops=dict(arrowstyle="->", color=colore, lw=0.9))

    for k, etichetta in ((1e-3, "1 nm"), (1.0, "V4"), (1e6, "1 m")):
        a2.axvline(ingombro(k), color=GRIGIO, lw=0.6, alpha=0.35)
        a2.text(ingombro(k) * 1.25, 3e-14, etichetta, fontsize=8, color=GRIGIO,
                rotation=90, va="bottom")

    a2.set_xlabel("ingombro del dispositivo (m)")
    a2.set_ylabel("spinta / peso (a terra)")
    a2.set_title("(b) A che scala il dispositivo si sostiene da solo")
    a2.legend(loc="lower left")
    a2.set_xlim(ing[0], ing[-1])
    a2.set_ylim(1e-18, 1e8)

    salva(fig, "fig_SP2_terra_vs_spazio.png")


# ---------------------------------------------------------------------------
# SP3 - i tre materiali a 1 m
# ---------------------------------------------------------------------------
def fig_sp3():
    d = carica("SP3_materiali_1m.json")["materiali"]
    nomi = list(d)
    colori = (BLU, ARANCIO, VERDE)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.8, 3.9))

    y = np.arange(len(nomi))
    massa = [d[m]["massa_totale_kg"] for m in nomi]
    a1.barh(y, massa, color=colori, alpha=0.85, height=0.5)
    for yy, vv, m in zip(y, massa, nomi):
        a1.text(vv * 1.03, yy, "%.3f kg   (I$_1$ = %.2e kg m$^2$)"
                % (vv, d[m]["spin"]["inerzia"]["momento_kgm2"]),
                va="center", fontsize=8)
    a1.set_yticks(y)
    a1.set_yticklabels(nomi)
    a1.set_xlim(0, max(massa) * 1.55)
    a1.set_xlabel("massa totale del sistema a 1 m (kg)")
    a1.set_title("(a) Massa inerziale")
    pulisci(a1)
    a1.grid(axis="y", visible=False)

    omega_lim = [d[m]["spin"]["meccanica"]["omega_limite_rad_s"] for m in nomi]
    margine = [d[m]["spin"]["meccanica"]["margine"] for m in nomi]
    a2.barh(y, omega_lim, color=colori, alpha=0.85, height=0.5)
    for yy, vv, mm in zip(y, omega_lim, margine):
        a2.text(vv * 1.04, yy, "%.0f rad s$^{-1}$   margine 10$^{%.0f}$"
                % (vv, math.log10(mm)), va="center", fontsize=8)
    a2.axvline(4.5e-11, color=ROSSO, lw=1.6)
    a2.text(4.5e-11 * 1.5, -0.62, "ω operativa 4.5×10$^{-11}$ rad s$^{-1}$",
            fontsize=8, color=ROSSO)
    a2.set_xscale("log")
    a2.set_yticks(y)
    a2.set_yticklabels(nomi)
    a2.set_xlim(1e-11, max(omega_lim) * 3.5)
    a2.set_xlabel("$\\omega$ che porta la tensione centrifugua allo snervamento (rad s$^{-1}$)")
    a2.set_title("(b) Limite meccanico, stesso a tutte le scale")
    pulisci(a2)
    a2.grid(axis="y", visible=False)

    salva(fig, "fig_SP3_materiali_1m.png")


# ---------------------------------------------------------------------------
# SP4 - rapporto 45/376
# ---------------------------------------------------------------------------
def fig_sp4():
    d = carica("SP4_rapporto_accoppiato.json")["risultati"]
    scale = ["10^-9_m", "massa_10^-18_kg", "1_m"]

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.8, 4.0))

    ing = np.array([ingombro(1e-3), ingombro(1.0), ingombro(1e6)])
    libero_s = [d[c]["dichiarata"]["omega_spin_rad_s"] for c in scale]
    libero_c = [d[c]["dichiarata"]["omega_contro_rad_s"] for c in scale]
    acc_s = [carica("SP1_spazio_scale.json")["casi"][c]["accoppiato"]["dichiarata"]
             ["omega_spin_rad_s"] for c in scale]
    acc_c = [carica("SP1_spazio_scale.json")["casi"][c]["accoppiato"]["dichiarata"]
             ["omega_contro_rad_s"] for c in scale]

    a1.loglog(ing, np.abs(libero_s), "o-", color=BLU, lw=2.0, ms=6,
              label="Spin libero")
    a1.loglog(ing, np.abs(libero_c), "s--", color=ARANCIO, lw=2.0, ms=6,
              label="Contro-spin libero")
    a1.loglog(ing, np.abs(acc_s), "o:", color=BLU, lw=2.0, ms=6, mfc="white",
              label="Spin accoppiato")
    a1.loglog(ing, np.abs(acc_c), "s:", color=ARANCIO, lw=2.0, ms=6, mfc="white",
              label="Contro-spin accoppiato")
    a1.set_xlabel("ingombro del dispositivo (m)")
    a1.set_ylabel("$|\\omega|$ (rad s$^{-1}$)")
    a1.set_title("(a) Il vincolo $I_1\\omega_1 + I_2\\omega_2 = 0$ raddoppia lo Spin")
    a1.legend(loc="lower left", fontsize=8)
    a1.set_ylim(1e-12, 1e4)

    target = -n.RATIO
    serie = [
        ("Cu", BLU, "o", [d[c]["dichiarata"]["omega_contro_su_spin"] for c in scale],
         [1.0, 2.0, 3.0]),
        ("Grafene/CNT", VERDE, "^",
         [d["1_m_materiali"][m]["dichiarata"]["omega_contro_su_spin"]
          for m in ("Grafene/CNT",)], [4.4]),
        ("SiC", ARANCIO, "s",
         [d["1_m_materiali"][m]["dichiarata"]["omega_contro_su_spin"]
          for m in ("SiC",)], [5.2]),
        ("Metamateriale Si/SiO2", GRIGIO, "D",
         [d["1_m_materiali"][m]["dichiarata"]["omega_contro_su_spin"]
          for m in ("Metamateriale Si/SiO2",)], [6.0]),
    ]
    for etichetta, colore, marcatore, valori, xs in serie:
        a2.plot(xs, valori, marcatore, color=colore, ms=7, label=etichetta)

    a2.axhline(target, color=ROSSO, lw=1.4, ls="--")
    a2.axhspan(target * 1.001, target * 0.999, color=ROSSO, alpha=0.10, lw=0)
    a2.text(0.55, target * 1.0016, "$-45/376 = -0.1196809$", fontsize=8.5, color=ROSSO)
    a2.text(0.55, target * 0.9973,
            "banda ±0.1 %: la geometria dà $-0.1196676$\n"
            "($-1.1\\times10^{-4}$), la massa di\nbilanciamento di P3 chiude a zero",
            fontsize=8, color=GRIGIO, va="top")
    a2.set_xticks([1, 2, 3, 4.4, 5.2, 6.0])
    a2.set_xticklabels(["1 nm", "V4", "1 m", "Grafene", "SiC", "Si/SiO$_2$"],
                       fontsize=8.5)
    a2.set_xlim(0.4, 7.0)
    a2.set_ylim(target * 1.004, target * 0.996)
    a2.set_ylabel("$\\omega_2/\\omega_1$")
    a2.set_title("(b) Il rapporto non dipende né dalla scala né dal materiale")
    a2.legend(loc="lower right", fontsize=8, ncols=2)

    salva(fig, "fig_SP4_rapporto_45_376.png")


def main():
    print("=" * 68)
    print("MOTORE FOTONICO V5 - GRAFICI DEL MODULO SPAZIALE")
    print("=" * 68)
    fig_sp2()
    fig_sp4()
    # SP1 e SP3 sono ridisegnate da v5_coppia_unificata.py sull'unico modello
    # di coppia ammesso, tau = P_ass/omega_ott.
    print("Fatto: figure in", FIG)


if __name__ == "__main__":
    main()
