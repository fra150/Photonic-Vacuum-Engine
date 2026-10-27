"""
MOTORE FOTONICO V5.1 - GRAFICI DEL VUOTO PROFONDO
Autore: Dr. Francesco Bulla

Legge i JSON prodotti da v51_vuoto.py e disegna:
  fig_V1  spin-up nel vuoto con il muro di rottura
  fig_V2  modalita' pulsata a scalini
  fig_V3  confronto Rame vs Grafene/CNT
  fig_V4  delta-v per tre modelli di spinta e accelerazione in g

Uso:
    python v51_grafici.py
"""

import json
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

BASE = os.path.dirname(os.path.abspath(__file__))
DATI = os.path.join(BASE, "dati")
FIG = os.path.join(BASE, "figure")
os.makedirs(FIG, exist_ok=True)

RAME = "#b06a2c"
GRAF = "#1f4e79"
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
    print("  %-38s %6.1f kB" % (nome, os.path.getsize(percorso) / 1024.0))


SCALE_ETICHETTA = {"nano_240nm": "240 nm", "macro_1m": "1 m"}
MAT_COLORE = {"Rame": RAME, "Grafene/CNT": GRAF}
MAT_ETICHETTA = {"Rame": "Rame", "Grafene/CNT": "Grafene/CNT"}


# ---------------------------------------------------------------------------
# V1 - spin-up nel vuoto
# ---------------------------------------------------------------------------
def fig_spin_up(spin, sintesi):
    fig, assi = plt.subplots(1, 3, figsize=(13.4, 3.9))

    # (a) e (b): stessa fisica, due scale di tempo molto diverse
    for ax, materiale, lettera in ((assi[0], "Rame", "a"),
                                   (assi[1], "Grafene/CNT", "b")):
        d = spin["casi"]["nano_240nm|%s|1_uW" % materiale]
        t = np.array([p[0] for p in d["curva"]])
        w = np.array([p[1] for p in d["curva"]])
        w_lim = d["omega_rottura_rad_s"]
        ax.plot(t, w, color=MAT_COLORE[materiale], lw=1.9,
                label="$\\omega = \\alpha\\,t$  (accelerazione costante)")
        ax.axhline(w_lim, color=ROSSO, ls="--", lw=1.3,
                   label="muro di rottura  $\\sigma_{max}=\\rho\\,\\omega^2R^2$")
        ax.plot([d["t_critico_s"]], [w_lim], "o", color=ROSSO, ms=5.5,
                mfc="white", mew=1.3, zorder=5)
        ax.annotate("$t_{critico}$ = %.2f s" % d["t_critico_s"],
                    xy=(d["t_critico_s"], w_lim),
                    xytext=(d["t_critico_s"] * 0.55, w_lim * 0.78),
                    fontsize=8.2, color=ROSSO,
                    arrowprops=dict(arrowstyle="->", color=ROSSO, lw=0.8))
        ax.set_xlabel("tempo  [s]")
        ax.set_ylabel("$\\omega$  [rad s$^{-1}$]")
        ax.set_title("(%s) %s — 240 nm, 1 µW\n$\\beta=0$: nessun regime"
                     % (lettera, MAT_ETICHETTA[materiale]), fontsize=9.5)
        ax.set_xlim(0, d["t_critico_s"] * 1.22)
        ax.set_ylim(0, w_lim * 1.18)
        ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        ax.legend(loc="lower right", fontsize=7.6)
        pulisci(ax)

    # (c) soglia di rottura al variare della scala
    ax2 = assi[2]
    mat = carica("V1_parametri_materiali.json")
    ingombri = np.logspace(-8, 1, 160)
    raggi = ingombri / 2.0
    for materiale in ("Rame", "Grafene/CNT"):
        sigma = mat[materiale]["sigma_max_Pa"]
        rho = mat[materiale]["rho_kg_m3"]
        w_lim = np.sqrt(sigma / (rho * raggi ** 2))
        ax2.loglog(ingombri, w_lim, color=MAT_COLORE[materiale], lw=1.9,
                   label=MAT_ETICHETTA[materiale])

    for L, testo, dy in ((240e-9, "240 nm", 1.6), (1.0, "1 m", 0.06)):
        ax2.axvline(L, color=GRIGIO, ls=":", lw=1.0)
        ax2.text(L * dy, 2.0e1, " %s" % testo, fontsize=8, color=GRIGIO,
                 rotation=90, va="bottom", ha="right")
    ax2.set_xlabel("ingombro  [m]")
    ax2.set_ylabel("$\\omega_{rottura}$  [rad s$^{-1}$]")
    ax2.set_title("(c) soglia di rottura al variare della scala",
                  fontsize=9.5)
    ax2.legend(loc="upper right", fontsize=8.5)
    pulisci(ax2)

    salva(fig, "fig_V1_spin_up_vuoto.png")


# ---------------------------------------------------------------------------
# V2 - modalita' pulsata
# ---------------------------------------------------------------------------
def fig_pulsata(puls):
    fig, assi = plt.subplots(2, 2, figsize=(11.8, 6.8))

    # (a) e (b): omega a scalini a piena potenza su tutta l'ora
    for ax, materiale, lettera in ((assi[0][0], "Rame", "a"),
                                   (assi[0][1], "Grafene/CNT", "b")):
        d = puls["casi"]["nano_240nm|%s" % materiale]
        curva = np.array(d["curva"])
        t, w = curva[:, 0], curva[:, 1]
        w_lim = d["omega_rottura_rad_s"]
        ax.plot(t, w / w_lim, color=MAT_COLORE[materiale], lw=1.1)
        ax.axhline(1.0, color=ROSSO, ls="--", lw=1.3)
        ax.set_xlim(0, 3600)
        ax.set_ylim(0, 1.32)
        ax.text(0.04, 0.90,
                "$\\omega$ finale = %.1f %% del muro\n(%d impulsi, %.3f s accesi)"
                % (100.0 * d["frazione_muro_raggiunta"], d["cicli"],
                   d["tempo_on_totale_s"]),
                transform=ax.transAxes, fontsize=8.2, color=GRIGIO, va="top")
        ax.set_xlabel("tempo  [s]")
        ax.set_ylabel("$\\omega\\,/\\,\\omega_{rottura}$")
        ax.set_title("(%s) %s — $\\omega(t)$ a piena potenza"
                     % (lettera, MAT_ETICHETTA[materiale]), fontsize=9.5)
        pulisci(ax)

    # (c) potenza scelta per fermarsi a meta' del muro in un'ora
    ax = assi[1][0]
    for materiale, stile in (("Rame", "-"), ("Grafene/CNT", "--")):
        s = puls["casi"]["nano_240nm|%s" % materiale]["sicura"]
        curva = np.array(s["curva"])
        ax.plot(curva[:, 0] / 3600.0, curva[:, 1] / s["omega_rottura_rad_s"],
                color=MAT_COLORE[materiale], lw=1.3, ls=stile,
                label="%s — $P$ = %.2f µW" % (MAT_ETICHETTA[materiale],
                                              s["potenza_W"] * 1e6))
    ax.axhline(1.0, color=ROSSO, ls="--", lw=1.3, label="muro di rottura")
    ax.axhline(0.5, color=GRIGIO, ls=":", lw=1.0, label="meta' del muro")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.35)
    ax.set_xlabel("tempo  [h]")
    ax.set_ylabel("$\\omega\\,/\\,\\omega_{rottura}$")
    ax.set_title("(c) potenza con margine 2: 3596 impulsi in un'ora",
                 fontsize=9.5)
    ax.legend(loc="lower right", fontsize=7.8)

    # riquadro: i primi impulsi a 1 µW, con la scala verticale espansa
    ins = ax.inset_axes([0.09, 0.45, 0.40, 0.40])
    massimo = 0.0
    for materiale in ("Rame", "Grafene/CNT"):
        d = puls["casi"]["nano_240nm|%s" % materiale]
        curva = np.array(d["curva"][:140])
        ins.plot(np.arange(len(curva)) / 2.0,
                 curva[:, 1] / d["omega_rottura_rad_s"],
                 color=MAT_COLORE[materiale], lw=1.0)
        massimo = max(massimo, float((curva[:, 1] / d["omega_rottura_rad_s"]).max()))
    ins.set_xlim(0, 70)
    ins.set_ylim(0, massimo * 1.10)
    ins.text(0.04, 0.90, "primi 70 impulsi, 1 µW", transform=ins.transAxes,
             fontsize=6.0, color=GRIGIO, va="top")
    ins.text(0.97, 0.06, "scala espansa\n(muro = 1.00 fuori scala)",
             transform=ins.transAxes, fontsize=5.8, color=GRIGIO,
             ha="right", va="bottom")
    ins.set_xlabel("impulsi", fontsize=6.4, labelpad=1.0)
    ins.tick_params(labelsize=6.4)
    ins.set_facecolor("#fafafa")
    pulisci(ax)

    # (d) velocita' lineare accumulata nell'ora, piena potenza vs margine 2
    ax = assi[1][1]
    x = np.arange(2)
    larghezza = 0.34
    for i, (chiave, etichetta, colore) in enumerate(
            (("piena", "1 µW, un'ora di impulsi", "#9e9e9e"),
             ("sicura", "potenza con margine 2", None))):
        valori = []
        for materiale in ("Rame", "Grafene/CNT"):
            d = puls["casi"]["nano_240nm|%s" % materiale]
            valori.append(d["delta_v_m_s"] if chiave == "piena"
                          else d["sicura"]["delta_v_m_s"])
        colori = [MAT_COLORE[m] if colore is None else colore
                  for m in ("Rame", "Grafene/CNT")]
        ax.bar(x + (i - 0.5) * larghezza, valori, larghezza, color=colori,
               alpha=1.0 if colore is None else 0.55, label=etichetta)
        for xi, vi in zip(x + (i - 0.5) * larghezza, valori):
            ax.text(xi, vi * 1.25, "%.2f" % vi, ha="center", fontsize=7.2)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([MAT_ETICHETTA[m] for m in ("Rame", "Grafene/CNT")])
    ax.set_ylabel("$\\Delta v$ in un'ora  [m s$^{-1}$]")
    ax.set_title("(d) $\\Delta v$: 1 µW contro potenza con margine 2", fontsize=9.5)
    ax.legend(loc="upper left", fontsize=8.0)
    ax.set_ylim(1e0, 1e4)
    pulisci(ax)

    salva(fig, "fig_V2_pulsata.png")


# ---------------------------------------------------------------------------
# V3 - confronto Rame vs Grafene
# ---------------------------------------------------------------------------
def fig_confronto(sintesi):
    fig, assi = plt.subplots(1, 3, figsize=(11.6, 3.9))
    scale = ("nano_240nm", "macro_1m")
    x = np.arange(len(scale))
    larghezza = 0.36

    for i, materiale in enumerate(("Rame", "Grafene/CNT")):
        masse = [sintesi["confronto"][s]["%s" % materiale]["massa_totale_kg"]
                 for s in scale]
        tempi = [sintesi["confronto"][s]["%s" % materiale]["t_critico_s"]
                 for s in scale]
        wlim = [sintesi["confronto"][s]["%s" % materiale]["omega_rottura_rad_s"]
                for s in scale]
        assi[0].bar(x + (i - 0.5) * larghezza, masse, larghezza,
                    color=MAT_COLORE[materiale], label=MAT_ETICHETTA[materiale])
        assi[1].bar(x + (i - 0.5) * larghezza, wlim, larghezza,
                    color=MAT_COLORE[materiale], label=MAT_ETICHETTA[materiale])
        assi[2].bar(x + (i - 0.5) * larghezza, tempi, larghezza,
                    color=MAT_COLORE[materiale], label=MAT_ETICHETTA[materiale])
        for ax, serie in zip(assi, (masse, wlim, tempi)):
            for xi, vi in zip(x + (i - 0.5) * larghezza, serie):
                ax.text(xi, vi * 2.0, "%.2e" % vi, ha="center", fontsize=6.8,
                        rotation=90)

    assi[0].set_yscale("log")
    assi[0].set_ylabel("massa totale  [kg]")
    assi[0].set_title("(a) massa inerziale", fontsize=9.5)

    assi[1].set_yscale("log")
    assi[1].set_ylabel("$\\omega_{rottura}$  [rad s$^{-1}$]")
    assi[1].set_title("(b) velocita' di rottura", fontsize=9.5)

    assi[2].set_yscale("log")
    assi[2].set_ylabel("$t_{critico}$  [s]")
    assi[2].set_title("(c) vita a piena potenza (1 µW)", fontsize=9.5)

    for ax in assi:
        basso, alto = ax.get_ylim()
        ax.set_ylim(basso, alto * 60.0)
        ax.set_xticks(x)
        ax.set_xticklabels([SCALE_ETICHETTA[s] for s in scale])
        pulisci(ax)
    assi[0].legend(loc="upper left", fontsize=8.5)
    salva(fig, "fig_V3_confronto_rame_grafene.png")


# ---------------------------------------------------------------------------
# V4 - delta-v in un'ora
# ---------------------------------------------------------------------------
def fig_deltav(prop, sintesi):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.4, 4.0))

    scale = ("nano_240nm", "macro_1m")
    modelli = (("a_scattering_m_s2", "scattering Mie (V5)", "#1f4e79"),
               ("a_P_su_c_m_s2", "$P/c$ — limite superiore", "#c55a11"),
               ("a_2P_su_c_m_s2", "$2P/c$ — limite superiore", "#2e7d32"))

    # (a) barre orizzontali: due gruppi di scala, tre modelli di spinta
    y = []
    valori = []
    colori = []
    etichette_y = []
    for j, s in enumerate(scale):
        base = j * 4
        d = prop["casi"]["%s|Grafene/CNT|1_uW" % s]
        for i, (chiave, etichetta, colore) in enumerate(modelli):
            y.append(base + i)
            valori.append(d[chiave] * 3600.0)
            colori.append(colore)
        etichette_y.append((base + 1, SCALE_ETICHETTA[s]))

    ax1.barh(y, valori, 0.68, color=colori)
    for yi, vi in zip(y, valori):
        ax1.text(vi * 1.5, yi, "%.2e" % vi, va="center", fontsize=6.8)
    ax1.set_xscale("log")
    ax1.set_yticks([p for p, _ in etichette_y])
    ax1.set_yticklabels([t for _, t in etichette_y], fontsize=8.5)
    ax1.invert_yaxis()
    ax1.set_xlim(1e-16, 1e14)
    ax1.set_xlabel("$\\Delta v$ in un'ora  [m s$^{-1}$]")
    ax1.axvline(2.998e8, color=ROSSO, ls=":", lw=1.1)
    ax1.text(0.828, 0.55, "c", transform=ax1.transAxes, fontsize=8,
             color=ROSSO, va="center",
             bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="none"))
    ax1.set_title("(a) $\\Delta v$ in un'ora, Grafene/CNT, 1 µW", fontsize=9.5)
    ax1.legend(handles=[Patch(facecolor=c, label=l) for _, l, c in modelli],
               loc="lower right", fontsize=7.4)
    pulisci(ax1)

    # (b) accelerazione in g per le quattro combinazioni
    etichette = []
    valori = []
    colori = []
    for s in scale:
        for materiale in ("Rame", "Grafene/CNT"):
            d = prop["casi"]["%s|%s|1_uW" % (s, materiale)]
            etichette.append("%s\n%s" % (SCALE_ETICHETTA[s],
                                          MAT_ETICHETTA[materiale]))
            valori.append(d["a_scattering_g"])
            colori.append(MAT_COLORE[materiale])

    xb = np.arange(len(valori))
    ax2.bar(xb, valori, 0.62, color=colori)
    ax2.axhline(1.0, color=GRIGIO, ls=":", lw=1.0)
    ax2.text(3.42, 1.6, "1 g", fontsize=8, color=GRIGIO, ha="right")
    for xi, vi in zip(xb, valori):
        ax2.text(xi, vi * 2.0, "%.3g g" % vi, ha="center", fontsize=7.5)
    ax2.set_yscale("log")
    ax2.set_ylabel("accelerazione  [g]")
    ax2.set_xticks(xb)
    ax2.set_xticklabels(etichette, fontsize=8)
    ax2.set_title("(b) accelerazione dalla spinta di scattering, 1 µW",
                  fontsize=9.5)
    ax2.set_ylim(1e-19, 3e2)
    pulisci(ax2)

    salva(fig, "fig_V4_deltav.png")


# ---------------------------------------------------------------------------
# V5 - tabella comparativa Rame vs Grafene/CNT
# ---------------------------------------------------------------------------
def fig_tabella(sintesi, puls):
    scale = ("nano_240nm", "macro_1m")
    colonne = ["parametro"] + ["%s — %s" % (SCALE_ETICHETTA[s], m)
                               for s in scale for m in ("Rame", "Grafene/CNT")]
    chiavi = [
        ("massa totale  [kg]", lambda c: c["massa_totale_kg"]),
        ("momento d'inerzia $I_{tot}$  [kg m$^2$]", lambda c: c["I_tot_kgm2"]),
        ("$\\omega_{rottura}$  [rad s$^{-1}$]",
         lambda c: c["omega_rottura_rad_s"]),
        ("$v$ periferica di rottura  [m s$^{-1}$]",
         lambda c: c["v_periferica_rottura_m_s"]),
        ("$\\alpha = \\tau_{max}/I_{tot}$  [rad s$^{-2}$]",
         lambda c: c["alpha_rad_s2"]),
        ("vita a piena potenza (1 µW)  [s]", lambda c: c["t_critico_s"]),
        ("accelerazione di scattering  [g]",
         lambda c: c["accelerazione_scattering_g"]),
        ("$\\Delta v$ in un'ora, pulsata 1 µW  [m s$^{-1}$]",
         lambda c: c["delta_v_m_s"]),
        ("potenza sicura per un'ora  [W]",
         lambda c: c["potenza_sicura_W"]),
    ]

    righe = []
    for etichetta, estrattore in chiavi:
        riga = [etichetta]
        for s in scale:
            for m in ("Rame", "Grafene/CNT"):
                riga.append("%.3e" % estrattore(sintesi["confronto"][s][m])
                            if etichetta not in
                            ("$\\Delta v$ in un'ora, pulsata 1 µW  [m s$^{-1}$]",
                             "potenza sicura per un'ora  [W]")
                            else "%.3e" % estrattore(puls["casi"][
                                "%s|%s" % (s, m)]))
        righe.append(riga)

    fig, ax = plt.subplots(figsize=(11.4, 4.2))
    ax.axis("off")
    tab = ax.table(cellText=righe, colLabels=colonne, loc="center",
                   cellLoc="right", colLoc="center")
    tab.auto_set_font_size(False)
    tab.set_fontsize(8.2)
    tab.scale(1.0, 1.55)
    for (riga, colonna), cella in tab.get_celld().items():
        cella.set_edgecolor("#cfcfcf")
        if riga == 0:
            cella.set_facecolor("#ededed")
            cella.set_text_props(weight="bold")
        elif colonna > 0:
            colore = RAME if colonna in (1, 3) else GRAF
            cella.set_text_props(color=colore)
        if colonna == 0:
            cella.set_text_props(ha="left")
            cella._loc = "left"
    ax.set_title("Rame vs Grafene/CNT nel vuoto profondo ($g=0$, "
                 "$\\eta=0$, $\\beta=0$)", fontsize=10.5, pad=14)
    salva(fig, "fig_V5_tabella_comparativa.png")


def main():
    spin = carica("V2_spin_up.json")
    puls = carica("V3_pulsata.json")
    prop = carica("V4_propulsione.json")
    sintesi = carica("V5_sintesi.json")
    print("[V5.1] figure:")
    fig_spin_up(spin, sintesi)
    fig_pulsata(puls)
    fig_confronto(sintesi)
    fig_deltav(prop, sintesi)
    fig_tabella(sintesi, puls)


if __name__ == "__main__":
    main()
