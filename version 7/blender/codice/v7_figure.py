"""
MOTORE FOTONICO DEL VUOTO - V7.1
v7_figure.py - le figure della V7.1, tutte dai dati calcolati.

Nessun numero e' scritto a mano: ogni curva viene da v7_mappa.py (mappa) o da
v7_scale.py (scale).  Le linee di frontiera (Pe_r = 10, delta T = 20 K,
omega = 0.5 omega_rottura, T = 0.7 T_lim) sono ricalcolate con il modello,
non interpolate a occhio.

Figure
------
fig01_canali            i canali in (R, P), acqua; scattering solo manubrio
fig02_termica_vuoto     margine termico in vuoto: rame lucido e rame nero
fig03_omega_rottura     t_spin in vuoto (il vincolo operativo del volano)
fig04_finestre          la finestra valida [P_lo, P_hi] in funzione di R
fig05_geometrie         le cinque forme a parita' di canale e materiale
fig06_scala             da 1 nm a 1 m: Pe_r e omega
fig07_frontiera         P(Pe_r = 10): colore = materiale, tratto = canale

Uso:
    python v7_figure.py
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

import v7_mappa as m  # noqa: E402

DATI = os.path.normpath(os.path.join(BASE, "..", "dati"))
FIG = os.path.normpath(os.path.join(BASE, "..", "figure"))
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
    axes[2].set_xlim(0.1, 3.0)
    axes[2].text(0.03, 0.06, "solo forma non assiale:\nJ$_z$ non si conserva\n"
                             "oltre ~2 µm il fascio\nnon illumina il manubrio",
                 fontsize=7.2, color="white", transform=axes[2].transAxes)
    axes[0].set_ylabel("potenza incidente P (W)")
    axes[0].legend(loc="lower left", fontsize=8, labelcolor="white")
    cb = fig.colorbar(im, ax=axes, pad=0.012, fraction=0.026)
    cb.set_label("log$_{10}$ $Pe_r$")
    fig.suptitle("V7.1 - i canali di coppia dopo la correzione della simmetria: "
                 "lo scattering sopravvive solo per il manubrio", y=1.02)
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
            titolo = ("%s, sfera, vuoto $10^{-8}$ mbar\n"
                      "$\\varepsilon_{eff}$ da Mie: %.3f-%.3f"
                      % (LBL_MAT[mat], m.eps_eff(mat, 1e-6), m.eps_eff(mat, 3e-5)))
        ax.set_title(titolo, color=COL_MAT[mat])
        ax.legend(loc="lower left", fontsize=8)
    axes[0].set_ylabel("potenza incidente P (W)")
    cb = fig.colorbar(im, ax=axes, pad=0.012, fraction=0.03)
    cb.set_label("$T \\,/\\, T_{lim}$")
    fig.suptitle("V7.1 - in vuoto il tetto e' radiativo: il rame lucido e' uno specchio "
                 "nell'IR; il nero dipende da un parametro, non da Mie", y=1.02)
    _salva(fig, "fig02_termica_vuoto.png")


# ---------------------------------------------------------------------------
# fig03: omega / omega_rottura
# ---------------------------------------------------------------------------
def fig03(righe):
    """t_spin in vuoto: e' li' che il vincolo operativo e' il tempo di avvio.

    In acqua omega/omega_rottura resta sotto 1e-4: una mappa di quel rapporto
    non ha una linea bianca da mostrare.  Il vincolo vero in vuoto e' t_spin,
    il tempo per portare il rotore a meta' del limite meccanico.
    """
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
            ax.text(0.04, 0.06, "oltre ~2 µm il fascio\nnon illumina il manubrio",
                    fontsize=7.2, color="white", transform=ax.transAxes)
        ax.set_title("%s, %s, %s\nvuoto $10^{-8}$ mbar"
                     % (LBL_CAN[canale], LBL_MAT[mat], forma),
                     color=COL_CAN[canale])
    axes[0].set_ylabel("potenza incidente P (W)")
    cb = fig.colorbar(im, ax=axes, pad=0.012, fraction=0.03)
    cb.set_label("log$_{10}$ $t_{spin}$ (s)")
    fig.suptitle("V7.1 - in vuoto il vincolo e' il tempo di avvio: la linea bianca "
                 "e' $t_{spin}=1$ s", y=1.02)
    _salva(fig, "fig03_omega_rottura.png")


# ---------------------------------------------------------------------------
# fig04: la finestra valida
# ---------------------------------------------------------------------------
def fig04():
    with open(os.path.join(DATI, "M1_mappa.json"), encoding="utf-8") as fh:
        d = json.load(fh)
    fr = [f for f in d["frontiere"] if f.get("ha_finestra")]
    scelte = [
        ("ritardo", "vaterite", "sfera", "acqua", None),
        ("assorbimento", "cu_nero", "sfera", "acqua", None),
        ("assorbimento", "cu_lucido", "sfera", "acqua", None),
        ("ritardo", "vaterite", "sfera", "vuoto", 1e-8),
        ("assorbimento", "cu_nero", "sfera", "vuoto", 1e-8),
        ("scattering", "silice", "manubrio", "acqua", None),
    ]
    fig, ax = plt.subplots(figsize=(9.6, 5.2))
    colori = ["#0077b6", "#c1121f", "#b07d3a", "#588157", "#7b6cf6", "#d1495b"]
    for i, (canale, mat, forma, mezzo, p) in enumerate(scelte):
        trovato = None
        for f in fr:
            if (f["canale"], f["materiale"], f["forma"], f["mezzo"]) != (canale, mat, forma, mezzo):
                continue
            if (f["p_mbar"] is None) != (p is None):
                continue
            if p is not None and abs(f["p_mbar"] - p) > 1e-15:
                continue
            trovato = f
            break
        if not trovato:
            continue
        R = np.array([w["R"] for w in trovato["finestre"]]) * 1e6
        lo = np.array([w["P_lo"] for w in trovato["finestre"]])
        hi = np.array([w["P_hi"] for w in trovato["finestre"]])
        etichetta = "%s / %s %s / %s" % (LBL_CAN[canale], LBL_MAT[mat], forma,
                                         "acqua" if mezzo == "acqua" else "vuoto %.0e mbar" % p)
        col = colori[i % len(colori)]
        if len(R) > 1:
            ax.fill_between(R, lo, hi, alpha=0.10, color=col)
        ax.plot(R, lo, color=col, lw=1.5, marker="o" if len(R) == 1 else None,
                ms=5)
        ax.plot(R, hi, color=col, lw=1.5, ls="--",
                marker="s" if len(R) == 1 else None, ms=5)
        if len(R) == 1:
            ax.annotate("solo bordo griglia", (R[0], lo[0]), xytext=(6, 6),
                        textcoords="offset points", fontsize=6.5, color=col)
        ax.plot([], [], color=col, lw=3, alpha=0.6,
                label="%s\n[%.1e – %.1e W]" % (etichetta, float(np.nanmin(lo)),
                                              float(np.nanmax(hi))))
    ax.axhline(10e-3, color="#888888", lw=1.2, ls=":",
               label="tetto dei criteri: 10 mW")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("R (µm)")
    ax.set_ylabel("potenza incidente P (W)")
    ax.set_ylim(3e-5, 2e-2)
    ax.set_title("V7.1 - finestra valida: P continua ($Pe_r=10$) e tratteggiata "
                 "(tetto termico/meccanico); scattering solo per il manubrio")
    ax.legend(fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2)
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
    ax1.text(0.03, 0.06, "l'anello esce dal grafico: il fascio\npassa nel foro, $\\eta_{int}\\to 0$",
             fontsize=7.5, color="#555555", transform=ax1.transAxes)
    ax1.legend(fontsize=8)
    fig.suptitle("V7.1 - geometrie a confronto, canale a ritardo, vaterite in acqua; "
                 "medie su dispersione larga, niente picchi da caso monodisperso", y=1.02)
    _salva(fig, "fig05_geometrie.png")


# ---------------------------------------------------------------------------
# fig06: da 1 nm a 1 m
# ---------------------------------------------------------------------------
def fig06():
    with open(os.path.join(DATI, "S5_sweep_scala.json"), encoding="utf-8") as fh:
        d = json.load(fh)
    with open(os.path.join(DATI, "M1_mappa.json"), encoding="utf-8") as fh:
        mappa = json.load(fh)
    R = np.array([r["R_m"] for r in d["righe"]])
    pe = np.array([r["Pe_r"] for r in d["righe"]])
    w = np.array([r["omega"] for r in d["righe"]])
    em = np.array([r["eta_mecc"] for r in d["righe"]])
    # La fascia deve essere quella vera della mappa, non un intervallo a mano.
    fr = [f for f in mappa["frontiere"]
          if (f.get("canale"), f.get("materiale"), f.get("forma"), f.get("mezzo"))
          == ("ritardo", "vaterite", "sfera", "acqua") and f.get("ha_finestra")]
    if fr:
        rmin = fr[0]["R_min_finestra"]
        rmax = fr[0]["R_max_finestra"]
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
    if rmin is not None:
        for ax in axes:
            ax.axvspan(rmin, rmax, color="#0077b6", alpha=0.08)
            ax.text(rmin * 1.2, ax.get_ylim()[0] * 3,
                    "finestra mappa\n%.2f–%.0f µm" % (rmin * 1e6, rmax * 1e6),
                    fontsize=7, color="#0077b6")
    fig.suptitle("V7.1 - da 1 nm a 1 m: la rotazione esiste solo nella decade micrometrica "
                 "(vaterite, sfera, acqua, 10 mW); la fascia e' quella della mappa", y=1.02)
    fig.subplots_adjust(wspace=0.42, bottom=0.18)
    _salva(fig, "fig06_scala.png")


# ---------------------------------------------------------------------------
# fig07: la scala di potenza
# ---------------------------------------------------------------------------
def fig07():
    """P(Pe_r=10) per forma: colore = materiale, tratto = canale.

    Le combinazioni non applicabili sono rimosse: niente assorbimento su
    silice/vaterite, niente ritardo sui metalli, niente scattering su forme
    assiali.  Lo scattering compare solo nel pannello del manubrio.
    """
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
    fig.suptitle("V7.1 - potenza per orientare il rotore in acqua: colore = materiale, "
                 "tratto = canale; solo combinazioni applicabili", y=1.02)
    _salva(fig, "fig07_frontiera.png")


def main():
    print("V7.1 - figure")
    righe = carica_csv()
    fig01(righe)
    fig02(righe)
    fig03(righe)
    fig04()
    fig05()
    fig06()
    fig07()
    print("Fatto: figure in", FIG)


if __name__ == "__main__":
    main()
