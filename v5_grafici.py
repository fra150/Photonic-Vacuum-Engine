"""
MOTORE FOTONICO V5 - GRAFICI
Autore: Dr. Francesco Bulla

Legge i JSON prodotti da v5_dati.py e disegna un grafico per ciascuno dei
cinque sviluppi della sezione 5.3, piu' i grafici delle scale.

Uso:
    python v5_grafici.py
"""

import json
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
DATI = os.path.join(BASE, "dati")
FIG = os.path.join(BASE, "figure")
os.makedirs(FIG, exist_ok=True)

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
# P1 - doppio giro
# ---------------------------------------------------------------------------
def fig_p1():
    d = carica("P1_doppio_giro.json")
    camp = d["campioni"]
    frame = [c["frame"] for c in camp]
    s1 = [c["spin_deg"] for c in camp]
    s2 = [c["contro_deg"] for c in camp]

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 3.9))

    a1.plot(frame, s1, color=BLU, lw=2.4, label="Spin (asse 1)")
    a1.plot(frame, s2, color=ARANCIO, lw=2.4, label="Contro-spin (asse 2)")
    a1.axhline(720, color=BLU, ls=":", lw=1.0)
    a1.axhline(360, color=BLU, ls=":", lw=0.8, alpha=0.6)
    a1.annotate("720$^{\\circ}$ = due giri", xy=(150, 720), xytext=(96, 560),
                color=BLU, fontsize=9,
                arrowprops=dict(arrowstyle="->", color=BLU, lw=1.0))
    a1.annotate("86.17$^{\\circ}$", xy=(228, 86.17), xytext=(150, 190),
                color=ARANCIO, fontsize=9,
                arrowprops=dict(arrowstyle="->", color=ARANCIO, lw=1.0))
    a1.set_xlabel("fotogramma")
    a1.set_ylabel("angolo cumulato (gradi)")
    a1.set_title("(a) Chiavi assolute, mai ridotte modulo 360")
    a1.legend(loc="upper left")
    a1.set_xlim(1, 240)
    a1.set_ylim(0, 800)

    giri = np.array([0.0, d["giri_spin_totali"]])
    a2.plot(giri, giri * d["rapporto_razionale"]["valore"], color=ARANCIO, lw=2.2)
    a2.plot([1.0 / d["rapporto_razionale"]["valore"]], [1.0], marker="o", ms=7, color=BLU)
    a2.annotate("1 giro del Contro-spin\nogni %.2f giri dello Spin"
                % d["controgiri_per_giro_di_spin"],
                xy=(1.0 / d["rapporto_razionale"]["valore"], 1.0), xytext=(158, 3.4),
                fontsize=8.6, color=BLU,
                arrowprops=dict(arrowstyle="->", color=BLU, lw=1.0))
    a2.annotate("%d giri" % d["giri_contro_totali"],
                xy=(d["giri_spin_totali"], d["giri_contro_totali"]),
                xytext=(232, 40), fontsize=8.6, color=ARANCIO,
                arrowprops=dict(arrowstyle="->", color=ARANCIO, lw=1.0))
    a2.set_xlabel("giri compiuti dallo Spin")
    a2.set_ylabel("giri compiuti dal Contro-spin")
    a2.set_title("(b) Moto relativo: rapporto 45/376")
    a2.set_xlim(0, 400)
    a2.set_ylim(0, 52)

    fig.suptitle("P1 - Doppio giro e rapporto razionale 45/376", fontsize=11.5,
                 color=BLU, y=1.02)
    salva(fig, "fig_P1_doppio_giro.png")


# ---------------------------------------------------------------------------
# P2 - dissipazione
# ---------------------------------------------------------------------------
def fig_p2():
    d = carica("P2_dissipazione.json")
    ang = [c["theta_deg"] for c in d["coppia_angolo"]]
    t_rel = [c["tau_relativo"] for c in d["coppia_angolo"]]

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 3.9))

    a1.plot(ang, t_rel, color=BLU, lw=2.2)
    a1.axhline(1.0, color=GRIGIO, ls=":", lw=0.9)
    a1.axhline(0.0, color="#888888", lw=0.8)
    a1.fill_between(ang, 0, t_rel, color=BLU, alpha=0.10)
    a1.set_xlabel("angolo del rotore $\\theta$ (gradi)")
    a1.set_ylabel("$\\tau(\\theta)/\\tau_{max}$")
    a1.set_title("(a) Coppia a due lobi, dipendente dall'angolo")
    a1.set_xlim(0, 360)
    a1.annotate("la coppia cambia segno:\ne' il termine che rende\nil ciclo non banale",
                xy=(90, -0.096), xytext=(112, 0.30), fontsize=8.2, color=ROSSO,
                arrowprops=dict(arrowstyle="->", color=ROSSO, lw=1.0))

    tt = [c["t_s"] for c in d["transitorio"]]
    om = [c["omega_normalizzato"] for c in d["transitorio"]]
    a2.semilogx(tt, om, color=ARANCIO, lw=2.2)
    a2.axhline(1.0, color=GRIGIO, ls=":", lw=0.9)
    a2.set_xlabel("tempo (s)")
    a2.set_ylabel("$\\omega/\\omega_{ss}$")
    a2.set_title("(b) Transitorio, $\\tau_{in}$ = %.2e s" % d["tempo_inerziale_s"])
    a2.set_ylim(0, 1.12)

    txt = ("$\\beta$ telaio = %.2e N m s\n$\\beta$ disco pieno = %.2e N m s\n"
           "$\\omega_{ss}$ = %.1f rad/s\n$Pe_r$ = %.2f"
           % (d["beta_telaio_Nms"], d["beta_disco_pieno_Nms"],
              d["omega_ss_dichiarata_rad_s"], d["Pe_rot"]))
    a2.text(0.42, 0.05, txt, transform=a2.transAxes, fontsize=8.4, color=GRIGIO,
            bbox=dict(boxstyle="round,pad=0.4", fc="#f7f7f7", ec="#dddddd"))

    fig.suptitle("P2 - Termine dissipativo $\\tau_{diss}=-\\beta\\omega$", fontsize=11.5,
                 color=BLU, y=1.02)
    salva(fig, "fig_P2_dissipazione.png")


# ---------------------------------------------------------------------------
# P3 - chiusura esatta
# ---------------------------------------------------------------------------
def fig_p3():
    d = carica("P3_chiusura.json")
    st = d["stato_attuale"]

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 3.9))

    ver = d["verifica"]
    valori = [ver["errore_chiusura_gradi_prima"], 0.004]
    etichette = ["senza\nbilanciamento", "con massa\ndi progetto"]
    barre = a1.bar(etichette, valori, color=[ROSSO, VERDE], width=0.5)
    a1.text(barre[0].get_x() + 0.25, ver["errore_chiusura_gradi_prima"] * 0.97,
            "%.2f$^{\\circ}$" % ver["errore_chiusura_gradi_prima"],
            ha="center", va="top", fontsize=10, color="white")
    a1.text(barre[1].get_x() + 0.25, 0.004 * 0.9, "0$^{\\circ}$",
            ha="center", va="bottom", fontsize=10, color=VERDE)
    a1.set_ylim(0, ver["errore_chiusura_gradi_prima"] * 1.25)
    a1.set_ylabel("scarto dopo %d cicli (gradi)" % ver["cicli"])
    a1.set_title("(a) Errore di chiusura: %.2f$^{\\circ}$ da annullare"
                 % ver["errore_chiusura_gradi_prima"])

    r = [c["r_bilanciamento_nm"] for c in d["massa_bilanciamento_vs_raggio"]]
    m = [c["massa_kg"] for c in d["massa_bilanciamento_vs_raggio"]]
    a2.loglog(r, m, color=BLU, lw=2.2)
    a2.plot([42], [st["massa_bilanciamento_kg"]], marker="o", ms=6, color=ARANCIO)
    a2.annotate("42 nm\n%.2e kg" % st["massa_bilanciamento_kg"], xy=(42, st["massa_bilanciamento_kg"]),
                xytext=(52, st["massa_bilanciamento_kg"] * 1.8), fontsize=8.6,
                color=ARANCIO,
                arrowprops=dict(arrowstyle="->", color=ARANCIO, lw=1.0))
    a2.set_xlabel("raggio di bilanciamento (nm)")
    a2.set_ylabel("massa di rame (kg)")
    a2.set_title("(b) Massa necessaria al variare del raggio")

    fig.suptitle("P3 - Chiusura esatta del ciclo in 188 iterazioni", fontsize=11.5,
                 color=BLU, y=1.02)
    salva(fig, "fig_P3_chiusura.png")


# ---------------------------------------------------------------------------
# P4 - scattering
# ---------------------------------------------------------------------------
def fig_p4():
    d = carica("P4_scattering.json")
    sf = d["sfera_Q_vs_x"]
    cl = d["cilindro_Q_vs_x"]
    rif = d["caso_riferimento"]

    fig, assi = plt.subplots(2, 2, figsize=(9.6, 7.0))

    a = assi[0][0]
    a.loglog([c["x"] for c in sf], [max(c["Q_pr"], 1e-4) for c in sf], color=BLU, lw=2.0)
    a.axvline(rif["x_nucleo"], color=ARANCIO, ls="--", lw=1.0)
    diam = rif["x_nucleo"] * 532 / math.pi
    a.annotate("nucleo\nd = %.0f nm" % diam,
               xy=(rif["x_nucleo"], 0.2), fontsize=8.2, color=ARANCIO)
    a.set_xlabel("parametro di dimensione $x = 2\\pi a/\\lambda$")
    a.set_ylabel("$Q_{pr}$")
    a.set_title("(a) Sfera di rame: Mie esatta")

    a = assi[0][1]
    val = [(c["x"], max(c["Q_pr"], 1e-9)) for c in cl if c["x"] <= 0.35]
    a.loglog([v[0] for v in val], [v[1] for v in val], color=VERDE, lw=2.0)
    a.loglog([0.35, 2.0], [1.1, 1.1], color=GRIGIO, ls="--", lw=1.2)
    a.axvline(rif["x_filo"], color=ARANCIO, ls="--", lw=1.0)
    a.annotate("filo\n2.5 nm", xy=(rif["x_filo"], 2e-4), fontsize=8.2, color=ARANCIO)
    a.annotate("limite geometrico\n$Q_{pr}\\to 1.1$", xy=(0.9, 1.1), xytext=(0.12, 0.35),
               fontsize=8.0, color=GRIGIO,
               arrowprops=dict(arrowstyle="->", color=GRIGIO, lw=0.9))
    a.set_xlabel("parametro di dimensione $x = 2\\pi a/\\lambda$")
    a.set_ylabel("$Q_{pr}$")
    a.set_title("(b) Cilindro sottile: Rayleigh fino a $x\\approx0.35$")

    a = assi[1][0]
    P = [c["P_W"] for c in d["forza_vs_potenza"]]
    F = [c["F_dichiarata_N"] for c in d["forza_vs_potenza"]]
    a.loglog(P, F, color=GRIGIO, lw=1.8, label="stima V4 (intercetta tutto il fascio)")
    a.loglog(P, [p * rif["sigma_effettiva_m2"] / rif["area_fascio_m2"] / 2.998e8 for p in P],
             color=BLU, lw=2.0, label="Mie, waist 500 nm")
    a.loglog(P, [p * rif["sigma_effettiva_m2"] / rif["area_min_diffrazione_m2"] / 2.998e8 for p in P],
             color=ARANCIO, lw=2.0, ls="--", label="Mie, fuoco limitato")
    a.set_xlabel("potenza ottica (W)")
    a.set_ylabel("$F_z$ (N)")
    a.set_title("(c) La forza diventa una previsione")
    a.legend(fontsize=7.4, loc="upper left")

    a = assi[1][1]
    lam = [c["lambda_nm"] for c in d["sezione_per_lambda"]]
    sig = [c["sigma_pr_m2"] for c in d["sezione_per_lambda"]]
    a.semilogy(lam, sig, marker="o", ms=5, color=BLU, lw=1.8)
    a.axvline(532, color=ARANCIO, ls=":", lw=1.0)
    a.set_xlabel("lunghezza d'onda (nm)")
    a.set_ylabel("$\\sigma_{pr}$ (m$^2$)")
    a.set_title("(d) Sezione di pressione di radiazione")

    fig.suptitle("P4 - Scattering completo al posto di $\\eta_s=0.70$, $\\eta_c=0.85$",
                 fontsize=11.5, color=BLU, y=1.00)
    salva(fig, "fig_P4_scattering.png")


# ---------------------------------------------------------------------------
# P5 - protocollo
# ---------------------------------------------------------------------------
def fig_p5():
    d = carica("P5_protocollo.json")
    e1 = d["esperimento_1_coppia"]

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 3.9))

    P = [c["P_W"] for c in e1["coppia_vs_potenza"]]
    om = [c["omega_rad_s"] for c in e1["coppia_vs_potenza"]]
    crit = d["criterio"]
    a1.loglog(P, om, color=BLU, lw=2.2, label="$\\omega_{ss}=\\tau/\\beta$")
    a1.axhline(crit["D_rot_1_s"], color=ROSSO, ls="--", lw=1.2,
               label="soglia di coerenza $D_r$")
    a1.axvline(1e-6, color=GRIGIO, ls=":", lw=1.0)
    a1.axvline(crit["potenza_necessaria_W"], color=VERDE, ls="-.", lw=1.2)
    a1.annotate("incrocio\n%.0f nW" % (crit["potenza_necessaria_W"] * 1e9),
                xy=(crit["potenza_necessaria_W"], crit["D_rot_1_s"]),
                xytext=(2e-8, 3e2), fontsize=8.2, color=VERDE,
                arrowprops=dict(arrowstyle="->", color=VERDE, lw=1.0))
    a1.annotate("1 $\\mu$W di progetto\n$Pe_r$ = %.1f" % crit["Pe_rot"],
                xy=(1e-6, 57.4), xytext=(3e-6, 3e-3), fontsize=8.2, color=GRIGIO,
                arrowprops=dict(arrowstyle="->", color=GRIGIO, lw=1.0))
    a1.set_xlabel("potenza ottica (W)")
    a1.set_ylabel("$\\omega_{ss}$ (rad s$^{-1}$)")
    a1.set_title("(a) Esperimento 1: coppia ottica diretta")
    a1.legend(fontsize=7.6, loc="upper left")

    th = np.linspace(0, 900, 2000)
    segnale = 1 + np.cos(np.radians(th) / 2)
    a2.plot(th, segnale, color=BLU, lw=2.2)
    a2.axvline(360, color=GRIGIO, ls=":", lw=1.0)
    a2.axvline(720, color=GRIGIO, ls=":", lw=1.0)
    a2.annotate("inversione di segno\na 360$^{\\circ}$", xy=(360, 0.04), xytext=(402, 0.30),
                fontsize=8.4, color=ROSSO,
                arrowprops=dict(arrowstyle="->", color=ROSSO, lw=1.0))
    a2.annotate("ritorno identico\na 720$^{\\circ}$", xy=(720, 2.0), xytext=(560, 1.55),
                fontsize=8.4, color=VERDE,
                arrowprops=dict(arrowstyle="->", color=VERDE, lw=1.0))
    a2.set_xlabel("angolo del rotore (gradi)")
    a2.set_ylabel("intensita' normalizzata")
    a2.set_title("(b) Esperimento 2: periodicita' a $4\\pi$")
    a2.set_xlim(0, 900)
    a2.set_ylim(-0.05, 2.15)

    fig.suptitle("P5 - I due esperimenti che discriminano il modello",
                 fontsize=11.5, color=BLU, y=1.02)
    salva(fig, "fig_P5_protocollo.png")


# ---------------------------------------------------------------------------
# Scale
# ---------------------------------------------------------------------------
def fig_sweep():
    d = carica("S5_sweep_scala.json")
    s = d["serie"]
    ing = [c["ingombro_nm"] for c in s]

    fig, assi = plt.subplots(2, 2, figsize=(9.6, 7.0))
    marcatori = [(1900.0, "240 nm", ARANCIO), (1.9, "1 nm", VERDE),
                 (1.9e9, "1 m", ROSSO)]

    def segna(ax, chiave):
        ax.plot(ing, [max(c[chiave], 1e-30) for c in s], color=BLU, lw=2.0)

    def riferimenti(ax):
        for x, et, col in marcatori:
            ax.axvline(x, color=col, ls="--", lw=0.9, alpha=0.8)
        ax.set_xscale("log")

    a = assi[0][0]
    segna(a, "Pe")
    riferimenti(a)
    a.set_yscale("log")
    a.axhline(1.0, color=GRIGIO, ls=":", lw=1.0)
    a.set_ylabel("$Pe = F_z r/kT$")
    a.set_title("(a) Trasporto deterministico sopra $Pe=1$")

    a = assi[0][1]
    segna(a, "Pe_rot")
    riferimenti(a)
    a.set_yscale("log")
    a.axhline(1.0, color=GRIGIO, ls=":", lw=1.0)
    a.set_ylabel("$Pe_r = \\omega_{ss}/D_r$")
    a.set_title("(b) Coerenza rotazionale")

    a = assi[1][0]
    segna(a, "omega_ss_rad_s")
    riferimenti(a)
    a.set_yscale("log")
    a.set_ylabel("$\\omega_{ss}$ (rad s$^{-1}$)")
    a.set_xlabel("ingombro (nm)")
    a.set_title("(c) Velocita' di regime")

    a = assi[1][1]
    segna(a, "energia_su_kT")
    riferimenti(a)
    a.set_yscale("log")
    a.axhline(1.0, color=GRIGIO, ls=":", lw=1.0)
    a.set_ylabel("$\\frac{1}{2}I\\omega^2/kT$")
    a.set_xlabel("ingombro (nm)")
    a.set_title("(d) Energia rotazionale contro rumore termico")

    handle = [plt.Line2D([], [], color=c, ls="--", lw=1.2, label=t)
              for _, t, c in marcatori]
    assi[1][0].legend(handles=handle, fontsize=7.8, loc="upper right")
    for ax in assi.ravel():
        pulisci(ax)

    fig.suptitle("Scale: 1 nm, 240 nm, 1 m - leggi di scala del motore",
                 fontsize=11.5, color=BLU, y=1.00)
    salva(fig, "fig_S1_sweep_scala.png")


def fig_1ghz():
    d = carica("S4_regime_1ghz.json")
    f = np.logspace(6, 10.2, 400)
    om = 2 * math.pi * f
    v = om * 327.5e-9
    sigma = 8960.0 * v ** 2

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 3.9))

    a1.loglog(f, sigma, color=BLU, lw=2.2, label="tensione centrifuga $\\rho v^2$")
    a1.axhline(d["snervamento_Cu_Pa"], color=ROSSO, ls="--", lw=1.2,
               label="snervamento del rame")
    a1.axvline(d["frequenza_limite_snervamento_Hz"], color=VERDE, ls=":", lw=1.2)
    a1.axvline(1e9, color=ARANCIO, ls="-.", lw=1.4)
    a1.annotate("1 GHz", xy=(1e9, 1e10), fontsize=8.4, color=ARANCIO, ha="right")
    a1.annotate("%.0f MHz" % (d["frequenza_limite_snervamento_Hz"] / 1e6),
                xy=(d["frequenza_limite_snervamento_Hz"], 2e8), xytext=(2.2e7, 6e8),
                fontsize=8.4, color=VERDE,
                arrowprops=dict(arrowstyle="->", color=VERDE, lw=1.0))
    a1.set_xlabel("frequenza di rotazione (Hz)")
    a1.set_ylabel("tensione (Pa)")
    a1.set_title("(a) Limite strutturale")
    a1.legend(fontsize=7.8, loc="lower right")

    a2.loglog(f, v, color=BLU, lw=2.2, label="velocita' tangenziale")
    a2.axhline(1480, color=ROSSO, ls="--", lw=1.2, label="velocita' del suono in acqua")
    a2.axvline(1e9, color=ARANCIO, ls="-.", lw=1.4)
    a2.set_xlabel("frequenza di rotazione (Hz)")
    a2.set_ylabel("velocita' (m s$^{-1}$)")
    a2.set_title("(b) Limite di cavitazione (Mach 1.39 a 1 GHz)")
    a2.legend(fontsize=7.8, loc="lower right")

    fig.suptitle("Regime a 1 GHz: la scala e' chiusa dalla resistenza del rame",
                 fontsize=11.5, color=BLU, y=1.02)
    salva(fig, "fig_S2_regime_1ghz.png")


def main():
    print("=" * 66)
    print("MOTORE FOTONICO V5 - GRAFICI")
    print("=" * 66)
    fig_p1()
    fig_p3()
    fig_p4()
    # P2, P5 e lo sweep di scala seguono l'unico modello di coppia ammesso,
    # tau = P_ass/omega_ott: le ridisegna v5_coppia_unificata.py.
    import v5_coppia_unificata as unificata
    casi = unificata.tabella()
    unificata.fig_p2(casi)
    unificata.fig_p5(casi)
    unificata.fig_sweep()
    fig_1ghz()
    print("-" * 66)
    print("Fatto: %d figure in %s" % (len(os.listdir(FIG)), FIG))


if __name__ == "__main__":
    main()
