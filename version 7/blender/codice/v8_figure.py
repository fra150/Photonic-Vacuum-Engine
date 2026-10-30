"""
MOTORE FOTONICO DEL VUOTO - V8.0
v8_figure.py - le figure della V8.0, tutte dai dati calcolati.

Nessun numero e' scritto a mano: ogni curva viene da v8_mappa.py (mappa) o
dal modello richiamato direttamente per la scala.  Le linee di frontiera
(Pe_r = 10, delta T = 20 K, omega = 0.5 omega_rottura, T = 0.7 T_lim) sono
ricalcolate con il modello, non interpolate a occhio.

La V8.0 introduce il dominio del fascio: R <= 2 um (z_R ~ 2 um).  Le figure
mostrano la banda grigia "fuori dominio del fascio (R > 2 um)" e le finestre
vengono disegnate solo dentro il dominio.

Figure
------
fig01_canali            i canali in (R, P), acqua; scattering solo manubrio
fig02_termica_vuoto     margine termico in vuoto: rame lucido e rame nero
fig03_omega_rottura     t_spin in vuoto (il vincolo operativo del volano)
fig04_finestre          la finestra valida [P_lo, P_hi] entro R <= 2 um
fig05_geometrie         le cinque forme a parita' di canale e materiale
fig06_scala             da 1 nm a 100 um: Pe_r e omega
fig07_frontiera         P(Pe_r = 10): colore = materiale, tratto = canale

Uso:
    python v8_figure.py
"""

import csv
import json
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D

BASE = os.path.dirname(os.path.abspath(__file__))
if BASE not in sys.path:
    sys.path.insert(0, BASE)

import v8_mappa as m  # noqa: E402

DATI = os.path.normpath(os.path.join(BASE, "..", "dati_v8"))
FIG = os.path.normpath(os.path.join(BASE, "..", "figure_v8"))
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 190, "font.size": 10,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
    "axes.edgecolor": "#444444", "axes.labelcolor": "#222222",
    "axes.titlesize": 11, "legend.frameon": False, "figure.facecolor": "white",
})

COL_MAT = {"cu_lucido": "#b07d3a", "cu_nero": "#3d3d3d",
           "silice": "#2e86ab", "vaterite": "#d1495b", "diamante": "#7b6cf6"}
COL_CAN = {"assorbimento": "#c1121f", "ritardo": "#0077b6",
           "scattering": "#588157"}
LBL_MAT = {"cu_lucido": "Cu lucido", "cu_nero": "Cu nero",
           "silice": "silice", "vaterite": "vaterite"}
LBL_CAN = {"assorbimento": "assorbimento", "ritardo": "ritardo di fase",
           "scattering": "scattering"}

R_MAX_UM = m.R_MAX_VALID * 1e6
R_MAX_VALID_M = m.R_MAX_VALID


# ---------------------------------------------------------------------------
# dati
# ---------------------------------------------------------------------------
def carica_csv():
    path = os.path.join(DATI, "M1_griglia.csv")
    righe = []
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            righe.append(r)
    return righe


def griglia(righe, canale, materiale, forma, mezzo, p_mbar):
    """Matrici (R, P) dei campi della griglia calcolata."""
    R = sorted({float(r["R"]) for r in righe})
    P = sorted({float(r["P"]) for r in righe})
    iR = {v: i for i, v in enumerate(R)}
    iP = {v: i for i, v in enumerate(P)}
    campi = ("Pe_r", "omega_frac", "t_spin", "t_spin_limite", "T", "dT", "tau")
    out = {k: np.full((len(P), len(R)), np.nan) for k in campi}
    for r in righe:
        if r["canale"] != canale or r["materiale"] != materiale:
            continue
        if r["forma"] != forma or r["mezzo"] != mezzo:
            continue
        if (r["p"] not in ("None", "")) and abs(float(r["p"]) - (p_mbar or 0)) > 1e-12:
            continue
        i, j = iP[float(r["P"])], iR[float(r["R"])]
        out["Pe_r"][i, j] = max(float(r["Pe_r"]), 1e-30)
        for k in campi[1:]:
            out[k][i, j] = float(r[k])
    out["R"] = np.array(R)
    out["P"] = np.array(P)
    return out


def frontiere_R(canale, materiale, forma, mezzo, p_mbar):
    """P(Pe_r = 10), P(dT = 20 K oppure T = 0.7 T_lim), P(omega = 0.5 omega_lim)."""
    R = list(m.R_GRID)
    p10, pterm, pmecc = [], [], []
    for r in R:
        c = m.cella(canale, 1e-3, r, materiale, forma, mezzo, p_mbar)
        if not c["tau"]:
            p10.append(np.nan); pterm.append(np.nan); pmecc.append(np.nan)
            continue
        tao = c["tau"] / 1e-3
        p10.append(10.0 * m.KT / tao)
        if mezzo == "acqua":
            pterm.append((20.0 / c["dT"]) * 1e-3 if c["dT"] > 0 else np.inf)
            pmecc.append(0.5 * c["omega_rottura"] * c["beta"] / tao)
        else:
            pterm.append(m._p_termico_vuoto(c, r, materiale, forma))
            pmecc.append(c["I"] * 0.5 * c["omega_rottura"] / m.CRITERI["t_spin_min_vuoto"] / tao)
    return np.array(R), np.array(p10), np.array(pterm), np.array(pmecc)


def banda_fuori_dominio(ax, testo=True, pos=(0.97, 0.95), fontsize=7.0):
    """Banda grigia R > 2 um, con etichetta in un angolo."""
    # zorder 1.5: sopra il pcolormesh (zorder 1), sotto le linee (zorder 2).
    ax.axvspan(R_MAX_UM, 30.0, color="#888888", alpha=0.08, zorder=1.5)
    if testo:
        ax.text(pos[0], pos[1],
                "fuori dominio del fascio\n($R > 2\\,\\mu$m)",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=fontsize, color="#555555",
                bbox=dict(facecolor="white", alpha=0.65, edgecolor="none",
                          pad=1.4))


def _salva(fig, nome):
    p = os.path.join(FIG, nome)
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    print("  %-26s %6.1f kB" % (nome, os.path.getsize(p) / 1024.0))


# ---------------------------------------------------------------------------
# fig01: i tre canali
# ---------------------------------------------------------------------------
def fig01(righe):
    # Lo scattering su sfera e' stato rimosso: per una forma assiale J_z si
    # conserva.  L'unico caso di scattering lecito e' il manubrio.
    casi = (("assorbimento", "cu_nero", "sfera", "acqua"),
            ("ritardo", "vaterite", "sfera", "acqua"),
            ("scattering", "silice", "manubrio", "acqua"))
    fig, axes = plt.subplots(1, 3, figsize=(13.4, 4.3), sharey=True)
    for ax, (canale, mat, forma, mezzo) in zip(axes, casi):
        d = griglia(righe, canale, mat, forma, mezzo, None)
        R, P, pe = d["R"], d["P"], d["Pe_r"]
        im = ax.pcolormesh(R * 1e6, P, np.log10(pe), cmap="magma",
                           norm=TwoSlopeNorm(vcenter=0.0, vmin=-6, vmax=6),
                           shading="nearest")
        X, Y = np.meshgrid(R * 1e6, P)
        ax.contour(X, Y, np.log10(pe), levels=[0.0, 1.0], colors=["white", "#7be0ff"],
                   linewidths=[1.0, 1.6])
        r, p10, pterm, pmecc = frontiere_R(canale, mat, forma, mezzo, None)
        ax.plot(r * 1e6, p10, color="white", lw=1.2, ls="--", label="$Pe_r=10$")
        pmax = np.minimum(pterm, pmecc)
        ax.plot(r * 1e6, pmax, color="#ffd166", lw=1.4, label="tetto (termico/meccanico)")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("R (µm)")
        ax.set_title("%s\n%s, %s, acqua" % (LBL_CAN[canale], LBL_MAT[mat], forma),
                     color=COL_CAN[canale])
        ax.set_ylim(1e-6, 1e-1)
        banda_fuori_dominio(ax, testo=(ax is axes[0]))
    axes[2].set_xlim(0.1, 3.0)
    axes[2].text(0.03, 0.06, "solo forma non assiale:\nJ$_z$ non si conserva\n"
                             "oltre 2 µm: estrapolazione",
                 fontsize=7.2, color="white", transform=axes[2].transAxes)
    axes[0].set_ylabel("potenza incidente P (W)")
    axes[0].legend(loc="lower left", fontsize=8, labelcolor="white")
    cb = fig.colorbar(im, ax=axes, pad=0.012, fraction=0.026)
    cb.set_label("log$_{10}$ $Pe_r$")
    fig.suptitle("V8.0 - i canali di coppia: lo scattering sopravvive solo per il "
                 "manubrio; R > 2 µm e' fuori dominio del fascio", y=1.02)
    _salva(fig, "fig01_canali.png")


# ---------------------------------------------------------------------------
# fig02: termica in vuoto
# ---------------------------------------------------------------------------
def fig02(righe):
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.3), sharey=True)
    for ax, mat in zip(axes, ("cu_lucido", "cu_nero")):
        d = griglia(righe, "assorbimento", mat, "sfera", "vuoto", 1e-8)
        R, P, T = d["R"], d["P"], d["T"]
        Tm = m.MATERIALI[mat]["Tm"]
        im = ax.pcolormesh(R * 1e6, P, T / Tm, cmap="inferno",
                           norm=LogNorm(vmin=0.2, vmax=20), shading="nearest")
        X, Y = np.meshgrid(R * 1e6, P)
        ax.contour(X, Y, T / Tm, levels=[0.7, 1.0], colors=["white", "#ffd166"],
                   linewidths=[1.6, 1.0])
        r, p10, pterm, pmecc = frontiere_R("assorbimento", mat, "sfera", "vuoto", 1e-8)
        ax.plot(r * 1e6, p10, color="#7be0ff", lw=1.3, ls="--", label="$Pe_r=10$")
        ax.plot(r * 1e6, pterm, color="white", lw=1.6, label="$T=0.7\\,T_{lim}$")
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("R (µm)")
        if mat == "cu_nero":
            titolo = ("%s, sfera, vuoto $10^{-8}$ mbar\n"
                      "$\\varepsilon_{IR}=0.30$ DICHIARATO (non Mie)"
                      % LBL_MAT[mat])
            ax.text(0.04, 0.80, "50 nm non basta: servono strati micron\n"
                                "bordo termico ~ 30 µm, non 21 µm",
                    fontsize=7.5, color="#333333", transform=ax.transAxes)
        else:
            e_lo = min(m.eps_eff(mat, 1e-6), m.eps_eff(mat, 3e-5))
            e_hi = max(m.eps_eff(mat, 1e-6), m.eps_eff(mat, 3e-5))
            titolo = ("%s, sfera, vuoto $10^{-8}$ mbar\n"
                      "$\\varepsilon_{eff}$ da Mie: %.3f-%.3f"
                      % (LBL_MAT[mat], e_lo, e_hi))
        ax.set_title(titolo, color=COL_MAT[mat])
        ax.legend(loc="lower left", fontsize=8)
        banda_fuori_dominio(ax, testo=(ax is axes[0]))
    axes[0].set_ylabel("potenza incidente P (W)")
    cb = fig.colorbar(im, ax=axes, pad=0.012, fraction=0.03)
    cb.set_label("$T \\,/\\, T_{lim}$")
    fig.suptitle("V8.0 - in vuoto il tetto e' radiativo: il rame lucido e' uno "
                 "specchio nell'IR; il nero dipende da un parametro; R > 2 µm fuori "
                 "dominio", y=1.02)
    _salva(fig, "fig02_termica_vuoto.png")


# ---------------------------------------------------------------------------
# fig03: omega / omega_rottura
# ---------------------------------------------------------------------------
def fig03(righe):
    """t_spin in vuoto: e' li' che il vincolo operativo e' il tempo di avvio."""
    casi = (("ritardo", "vaterite", "sfera", "vuoto", 1e-8),
            ("assorbimento", "cu_nero", "sfera", "vuoto", 1e-8),
            ("scattering", "silice", "manubrio", "vuoto", 1e-8))
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.3), sharey=True)
    for ax, (canale, mat, forma, mezzo, p) in zip(axes, casi):
        d = griglia(righe, canale, mat, forma, mezzo, p)
        R, P, ts = d["R"], d["P"], d["t_spin"]
        ts = np.where(np.isfinite(ts) & (ts > 0), ts, np.nan)
        z = np.log10(np.clip(ts, 1e-8, 1e8))
        im = ax.pcolormesh(R * 1e6, P, z, cmap="viridis",
                           norm=Normalize(vmin=-8, vmax=8), shading="nearest")
        X, Y = np.meshgrid(R * 1e6, P)
        ax.contour(X, Y, z, levels=[0.0], colors="white", linewidths=[1.6])
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("R (µm)")
        if forma == "manubrio":
            ax.set_xlim(0.1, 3.0)
            ax.text(0.04, 0.06, "oltre 2 µm: fuori dominio\n(estrapolazione)",
                    fontsize=7.2, color="white", transform=ax.transAxes)
        ax.set_title("%s, %s, %s\nvuoto $10^{-8}$ mbar"
                     % (LBL_CAN[canale], LBL_MAT[mat], forma),
                     color=COL_CAN[canale])
        banda_fuori_dominio(ax, testo=(ax is axes[0]))
    axes[0].set_ylabel("potenza incidente P (W)")
    cb = fig.colorbar(im, ax=axes, pad=0.012, fraction=0.03)
    cb.set_label("log$_{10}$ $t_{spin}$ (s)")
    fig.suptitle("V8.0 - in vuoto il vincolo e' il tempo di avvio: la linea bianca "
                 "e' $t_{spin}=1$ s; R > 2 µm fuori dominio del fascio", y=1.02)
    _salva(fig, "fig03_omega_rottura.png")


# ---------------------------------------------------------------------------
# fig04: la finestra valida entro R <= 2 um
# ---------------------------------------------------------------------------
def conta_valide(righe):
    out = {}
    for r in righe:
        if str(r["valido"]).strip().lower() != "true":
            continue
        p = None if r["mezzo"] == "acqua" else float(r["p"])
        k = (r["canale"], r["materiale"], r["forma"], r["mezzo"], p)
        out[k] = out.get(k, 0) + 1
    return out


def fig04(righe):
    with open(os.path.join(DATI, "M1_mappa.json"), encoding="utf-8") as fh:
        d = json.load(fh)
    conta = conta_valide(righe)
    ord_mezzo = {"acqua": 0, "vuoto": 1}
    ord_can = {"assorbimento": 0, "ritardo": 1, "scattering": 2}
    ord_mat = {"cu_lucido": 0, "cu_nero": 1, "silice": 2, "vaterite": 3}
    ord_for = {"sfera": 0, "disco": 1, "anello": 2, "manubrio": 3, "reticolo": 4}
    fr = sorted(d["riassunto"],
                key=lambda r: (ord_mezzo[r["mezzo"]],
                               0 if r["p_mbar"] is None else -math.log10(r["p_mbar"]),
                               ord_can[r["canale"]], ord_mat[r["materiale"]],
                               ord_for[r["forma"]]))
    colori = ["#c1121f", "#b07d3a", "#0077b6", "#d1495b", "#2e86ab", "#588157",
              "#7b6cf6", "#e07a5f"]
    fig, ax = plt.subplots(figsize=(10.0, 5.8))
    for i, r in enumerate(fr):
        fs = r["finestre"]
        R = np.array([w["R"] for w in fs]) * 1e6
        lo = np.array([w["P_lo"] for w in fs])
        hi = np.array([w["P_hi"] for w in fs])
        col = colori[i % len(colori)]
        etichetta = "%s / %s %s / %s" % (
            LBL_CAN[r["canale"]], LBL_MAT[r["materiale"]], r["forma"],
            "acqua" if r["mezzo"] == "acqua" else "vuoto %.0e mbar" % r["p_mbar"])
        if len(R) > 1:
            ax.fill_between(R, lo, hi, alpha=0.10, color=col)
        ax.plot(R, lo, color=col, lw=1.5, marker="o" if len(R) == 1 else None,
                ms=5)
        ax.plot(R, hi, color=col, lw=1.5, ls="--",
                marker="s" if len(R) == 1 else None, ms=5)
        n_cells = conta.get((r["canale"], r["materiale"], r["forma"],
                             r["mezzo"], r["p_mbar"]), 0)
        if n_cells == 0:
            xm = float(np.mean(R))
            ym = float(np.mean(np.sqrt(lo * hi))) if len(lo) else 0.0
            ax.annotate("solo frontiera:\nnessuna cella di griglia",
                        (xm, ym), xytext=(25, 12), textcoords="offset points",
                        fontsize=6.4, color=col,
                        arrowprops=dict(arrowstyle="-", color=col, lw=0.7))
        ax.plot([], [], color=col, lw=3, alpha=0.6,
                label="%s\n[%.1e – %.1e W]" % (etichetta, float(np.nanmin(lo)),
                                              float(np.nanmax(hi))))
    ax.axhline(10e-3, color="#888888", lw=1.2, ls=":",
               label="tetto dei criteri: 10 mW")
    ax.axvline(R_MAX_UM, color="#666666", lw=1.0, ls="--")
    ax.axvspan(R_MAX_UM, 3.0, color="#888888", alpha=0.08, zorder=1.5)
    ax.text(2.94, 1.15e-2, "fuori dominio\n($R > 2\\,\\mu$m)", fontsize=7,
            color="#555555", ha="right", va="top",
            bbox=dict(facecolor="white", alpha=0.65, edgecolor="none", pad=1.2))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(0.1, 3.0)
    ax.set_xlabel("R (µm)")
    ax.set_ylabel("potenza incidente P (W)")
    ax.set_ylim(3e-5, 2e-2)
    ax.set_title("V8.0 - finestra valida entro $R \\leq 2\\,\\mu$m")
    ax.legend(fontsize=6.8, loc="upper center", bbox_to_anchor=(0.5, -0.13),
              ncol=2)
    fig.subplots_adjust(bottom=0.30)
    _salva(fig, "fig04_finestre.png")


# ---------------------------------------------------------------------------
# fig05: le forme
# ---------------------------------------------------------------------------
def fig05():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.3))
    stili = {"sfera": "-", "disco": "--", "anello": ":", "manubrio": "-.",
             "reticolo": "-."}
    marker = {"reticolo": "o"}
    mecc_sfera = []
    for forma in m.FORME:
        r, p10, pterm, pmecc = frontiere_R("ritardo", "vaterite", forma, "acqua", None)
        # l'anello ha un foro centrale: per R >> w0 il fascio ci passa dentro,
        # eta_int -> 0 e la potenza richiesta diverge.  Si taglia, non si
        # aggiusta: la divergenza e' il risultato.
        p10 = np.where((p10 > 1e1) | ~np.isfinite(p10), np.nan, p10)
        pmecc = np.where((pmecc > 1e9) | ~np.isfinite(pmecc), np.nan, pmecc)
        if forma == "sfera":
            mecc_sfera = pmecc[np.isfinite(pmecc)]
        kw = {"marker": marker.get(forma), "markevery": 2 if forma == "reticolo" else None,
              "ms": 4}
        ax1.plot(r * 1e6, p10, ls=stili[forma], color=COL_CAN["ritardo"], lw=1.4,
                 label=forma, **kw)
        ax2.plot(r * 1e6, pmecc, ls=stili[forma], color=COL_CAN["ritardo"], lw=1.4,
                 label=forma, **kw)
    mecc_sfera = np.concatenate([mecc_sfera]) if len(mecc_sfera) else np.array([])
    ax1.set_xscale("log"); ax1.set_yscale("log")
    ax1.set_ylim(3e-5, 2e-2)
    ax1.set_xlabel("R (µm)"); ax1.set_ylabel("P per $Pe_r=10$ (W)")
    ax1.set_title("potenza necessaria")
    ax2.set_xscale("log"); ax2.set_yscale("log")
    ax2.set_ylim(1e-5, 1e9)
    ax2.set_xlabel("R (µm)"); ax2.set_ylabel("P per $\\omega=0.5\\,\\omega_{rottura}$ (W)")
    ax2.set_title("tetto meccanico")
    if len(mecc_sfera):
        ax2.text(0.04, 0.12, "per la sfera il tetto meccanico va da %.1e a %.1e W:\n"
                             "non si raggiunge, il rotore non si rompe"
                 % (float(np.nanmin(mecc_sfera)), float(np.nanmax(mecc_sfera))),
                 fontsize=7.5, color="#555555", transform=ax2.transAxes)
    ax1.text(0.45, 0.06, "l'anello esce dal grafico: il fascio\npassa nel foro, $\\eta_{\\rm int}\\to 0$",
             fontsize=7.5, color="#555555", transform=ax1.transAxes,
             bbox=dict(facecolor="white", alpha=0.65, edgecolor="none", pad=1.5))
    ax1.legend(fontsize=8)
    for ax in (ax1, ax2):
        banda_fuori_dominio(ax, testo=(ax is ax1))
    fig.suptitle("V8.0 - geometrie a confronto, canale a ritardo, vaterite in acqua; "
                 "R > 2 µm fuori dominio del fascio", y=1.02)
    _salva(fig, "fig05_geometrie.png")


# ---------------------------------------------------------------------------
# fig06: da 1 nm a 100 um
# ---------------------------------------------------------------------------
def fig06():
    # Il sweep e' ricalcolato con il modello V8.0: non si usa il vecchio
    # S5 della V7.1, che aveva un altro dn e un'altra validita' del fascio.
    R = np.logspace(-9, -4, 61)
    pe, w, em = [], [], []
    for r in R:
        c = m.cella("ritardo", 10e-3, float(r), "vaterite", "sfera", "acqua")
        pe.append(c["Pe_r"])
        w.append(c["omega"])
        em.append(c["eta_mecc"])
    pe = np.array(pe); w = np.array(w); em = np.array(em)
    # finestra della mappa V8.0 (vaterite, sfera, acqua)
    with open(os.path.join(DATI, "M1_mappa.json"), encoding="utf-8") as fh:
        mappa = json.load(fh)
    fr = [f for f in mappa["riassunto"]
          if (f.get("canale"), f.get("materiale"), f.get("forma"), f.get("mezzo"))
          == ("ritardo", "vaterite", "sfera", "acqua")]
    if fr:
        rmin = fr[0]["R_min_um"] * 1e-6
        rmax = fr[0]["R_max_um"] * 1e-6
    else:
        rmin, rmax = None, None
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.1))
    axes[0].loglog(R, pe, color="#0077b6", lw=1.8)
    axes[0].axhline(10, color="#c1121f", ls="--", lw=1.2)
    axes[0].text(2e-9, 13, "$Pe_r = 10$", color="#c1121f", fontsize=9)
    axes[0].set_xlabel("R (m)"); axes[0].set_ylabel("$Pe_r$")
    axes[0].set_title("determinismo (10 mW)")
    axes[1].loglog(R, np.abs(w), color="#d1495b", lw=1.8)
    axes[1].set_xlabel("R (m)"); axes[1].set_ylabel("$\\omega$ (rad/s)")
    axes[1].set_title("velocita' (10 mW)")
    axes[2].loglog(R, em, color="#588157", lw=1.8)
    axes[2].set_xlabel("R (m)"); axes[2].set_ylabel("$\\eta_{mecc}$")
    axes[2].set_title("efficienza meccanica (10 mW)")
    for ax in axes:
        ax.set_xlim(1e-9, 1e-4)
        # Asse in metri: la banda va da 2 um a 100 um, non da 2 a 30.
        ax.axvspan(R_MAX_VALID_M, 1e-4, color="#888888", alpha=0.08, zorder=1.5)
        ax.axvline(R_MAX_VALID_M, color="#666666", lw=0.9, ls="--")
        if rmin is not None:
            ax.axvspan(rmin, rmax, color="#0077b6", alpha=0.10, zorder=0)
    axes[0].text(0.97, 0.96, "fuori dominio del fascio\n($R > 2\\,\\mu$m)",
                 transform=axes[0].transAxes, ha="right", va="top",
                 fontsize=7, color="#555555",
                 bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", pad=1.4))
    if rmin is not None:
        axes[0].text(0.52, 0.06, "finestra mappa\n%.2f–%.2f µm"
                     % (rmin * 1e6, rmax * 1e6), transform=axes[0].transAxes,
                     fontsize=7, color="#0077b6")
    fig.suptitle("V8.0 - da 1 nm a 100 µm: la rotazione esiste solo nella decade "
                 "micrometrica (vaterite, sfera, acqua, 10 mW); oltre 2 µm e' "
                 "estrapolazione", y=1.02)
    fig.subplots_adjust(wspace=0.42, bottom=0.18)
    _salva(fig, "fig06_scala.png")


# ---------------------------------------------------------------------------
# fig07: la scala di potenza
# ---------------------------------------------------------------------------
def fig07():
    """P(Pe_r=10) per forma: colore = materiale, tratto = canale."""
    with open(os.path.join(DATI, "M1_mappa.json"), encoding="utf-8") as fh:
        d = json.load(fh)
    ls_can = {"assorbimento": "-", "ritardo": "--", "scattering": ":"}
    forme = ("sfera", "disco", "manubrio")
    fig, axes = plt.subplots(1, 3, figsize=(13.6, 4.6), sharey=True)
    for ax, forma in zip(axes, forme):
        for f in d["frontiere"]:
            if not f.get("esiste") or f["mezzo"] != "acqua" or f["forma"] != forma:
                continue
            canale, mat = f["canale"], f["materiale"]
            if canale == "assorbimento" and mat not in ("cu_lucido", "cu_nero"):
                continue
            if canale == "ritardo" and mat not in ("silice", "vaterite"):
                continue
            if canale == "scattering" and not (forma == "manubrio" and mat == "silice"):
                continue
            R = np.array([c["R"] for c in f["celle"]]) * 1e6
            p10 = np.array([c["p10"] for c in f["celle"]])
            p10 = np.where((p10 > 1e-1) | (p10 < 1e-7) | ~np.isfinite(p10), np.nan, p10)
            ax.plot(R, p10, color=COL_MAT[mat], ls=ls_can[canale], lw=1.5,
                    alpha=0.95)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("R (µm)")
        ax.set_title(forma)
        ax.set_ylim(1e-7, 1e-1)
        ax.set_xlim(0.1, 3.0 if forma == "manubrio" else 30.0)
        banda_fuori_dominio(ax, testo=(ax is axes[0]))
    axes[0].set_ylabel("P per $Pe_r = 10$ (W)")
    axes[2].text(0.04, 0.06, "scattering solo qui:\nforma non assiale",
                 fontsize=7.5, color="#555555", transform=axes[2].transAxes)
    leg_mat = [Line2D([0], [0], color=COL_MAT[k], lw=2.5, label=LBL_MAT[k])
               for k in ("cu_lucido", "cu_nero", "silice", "vaterite")]
    leg_can = [Line2D([0], [0], color="#444444", lw=1.8, ls=ls_can[k],
                      label=LBL_CAN[k]) for k in ("assorbimento", "ritardo", "scattering")]
    leg1 = axes[1].legend(handles=leg_mat, fontsize=8, ncol=2, loc="upper center",
                          bbox_to_anchor=(0.5, -0.16))
    axes[1].add_artist(leg1)
    axes[1].legend(handles=leg_can, fontsize=8, ncol=3, loc="upper center",
                   bbox_to_anchor=(0.5, -0.30))
    fig.suptitle("V8.0 - potenza per orientare il rotore in acqua: colore = materiale, "
                 "tratto = canale; R > 2 µm fuori dominio del fascio", y=1.02)
    _salva(fig, "fig07_frontiera.png")


def main():
    print("V8.0 - figure")
    righe = carica_csv()
    fig01(righe)
    fig02(righe)
    fig03(righe)
    fig04(righe)
    fig05()
    fig06()
    fig07()
    print("Fatto: figure in", FIG)


if __name__ == "__main__":
    main()
