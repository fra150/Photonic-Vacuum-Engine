"""
MOTORE FOTONICO DEL VUOTO - V7.3
v73_termica.py - termica a mezzo uniforme, con i tre contributi separati.

La V7.2 sommava due soli termini (particella + cammino).  La V7.3 aggiunge
i termini che il caso a mezzo uniforme rende visibili:

  dT_part  assorbimento della sfera (Mie + tetto geometrico)
  dT_pre   assorbimento del fascio dall'interfaccia alla sfera (V7.2)
  dT_post  assorbimento del fascio oltre il fuoco, fino a z_R
  dT_diff  assorbimento, nei pressi della sfera, della luce diffusa da lei

e tratta la parete con un segno dichiarato:
  sink   parete isoterma (vetro, endotelio): immagine negativa, limite
         inferiore = il termine di cammino viene annullato;
  mirror parete adiabatica: immagine positiva, limite superiore = il termine
         di cammino viene raddoppiato.

Il caso a due zone (CFL di plasma + sangue fuori) NON e' un ambiente V7.3:
il CFL esiste solo con flusso, e il flusso e' chiuso da v73_taglio.py.  Resta
solo come riferimento storico e non entra nelle tabelle.

Il criterio dT < 5 K e' un criterio termico SCELTO, non una prova di
sicurezza biologica.

Uso:
    python codice/v73_termica.py
"""

import json
import math
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(BASE, ".."))
if BASE not in sys.path:
    sys.path.insert(0, BASE)

import v7_mappa as m  # noqa: E402
import v73_contrappeso as c  # noqa: E402
import v73_taglio as t  # noqa: E402

DATI = os.path.join(ROOT, "dati")
TAB_IT = os.path.join(ROOT, "paper", "it", "tabelle")
TAB_EN = os.path.join(ROOT, "paper", "en", "tabelle")
FIG_IT = os.path.join(ROOT, "paper", "it", "figures")
FIG_EN = os.path.join(ROOT, "paper", "en", "figures")
FIG_ROOT_IT = os.path.join(ROOT, "figure", "it")
FIG_ROOT_EN = os.path.join(ROOT, "figure", "en")
for d in (DATI, TAB_IT, TAB_EN, FIG_IT, FIG_EN, FIG_ROOT_IT, FIG_ROOT_EN):
    os.makedirs(d, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. I DUE AMBIENTI VIVI
# ---------------------------------------------------------------------------
K_BLOOD = 0.55          # W/m/K
K_PLASMA = 0.60         # W/m/K, acqua/plasma
MU_A_PLASMA = {532: 0.044, 800: 2.3, 1064: 15.0}   # 1/m, acqua
MU_S_PLASMA = {532: 0.0, 800: 0.0, 1064: 0.0}      # mezzo chiaro

AMBIENTI = {
    "A": {"nome": "plasma/tampone in vitro",
          "mu_a": MU_A_PLASMA, "mu_s": MU_S_PLASMA, "kappa": K_PLASMA},
    "B": {"nome": "sangue intero in stasi",
          "mu_a": {k: c.OTTICA[k]["mu_a"] for k in (532, 800, 1064)},
          "mu_s": {k: c.OTTICA[k]["mu_s"] for k in (532, 800, 1064)},
          "kappa": K_BLOOD},
}

R_CASI = (0.3e-6, 0.4e-6, 0.5e-6)
LAM = (532, 800, 1064)
Z = 2.0e-6
W0 = 0.5e-6
D_T_CRITERIO = 5.0      # K, scelto


def trasmissione(mu_a, mu_s, z):
    return math.exp(-(mu_a + mu_s) * z)


def dT_pre(P_inc, lam, z, mu_a, mu_s, kappa, w0=W0):
    """Termine di cammino V7.2, con la media del deposito."""
    if z <= 0:
        return 0.0
    mu_t = mu_a + mu_s
    x = mu_t * z
    f = (1.0 - math.exp(-x)) / x if x > 1e-12 else 1.0
    return (mu_a * P_inc / (2.0 * math.pi * kappa)) * f * math.log1p(z / w0)


def dT_post(P_inc, lam, z, mu_a, mu_s, kappa, w0=W0):
    """Assorbimento oltre il fuoco, fino alla lunghezza di Rayleigh.

    P_focus = P_inc exp(-mu_t z); z_R = pi w0^2/lambda.
    dT_post = mu_a P_focus/(2 pi kappa) * ln(1 + z_R/w0).

    Il prefattore 2 pi (sorgente lineare) e' la scelta conservativa: con 4 pi
    (sorgente puntiforme) il termine dimezza.  La formula e' una stima di
    ordine zero, non una soluzione dell'equazione del calore.
    """
    mu_t = mu_a + mu_s
    P_focus = P_inc * math.exp(-mu_t * z)
    z_R = math.pi * w0 ** 2 / (lam * 1e-9)
    return (mu_a * P_focus / (2.0 * math.pi * kappa)) * math.log1p(z_R / w0)


def potenza_diffusa(P_inc, R, lam, w0=W0):
    """Potenza diffusa dalla sfera, limitata dalla potenza intercettata.

    P_sca = eta_int P_inc * Q_sca/Q_ext.  Il rapporto Q_sca/Q_ext e' l'albedo
    di scattering; il fattore eta_int evita di contare piu' luce di quella che
    investe la sfera (errore opposto a quello di usare I0 pi greco R^2 senza
    tetto quando R ~ w0).
    """
    q = m.q_ottico("vaterite", R, lam * 1e-9)
    eta = m.eta_int(R, "sfera", w0)
    albedo = q["Q_sca"] / q["Q_ext"] if q["Q_ext"] > 0 else 0.0
    return {"P_sca": eta * P_inc * albedo, "eta_int": eta,
            "Q_sca": q["Q_sca"], "Q_ext": q["Q_ext"], "albedo": albedo}


def dT_diff(P_sca, mu_a, kappa):
    """Assorbimento della luce diffusa nei pressi della sfera.

    dT_diff ~ mu_a P_sca/(4 pi kappa).  Vale nei due limiti mu_a R << 1
    (assorbe una frazione mu_a R in un volume ~R^3) e mu_a R >> 1 (assorbe
    entro 1/mu_a): in entrambi il prodotto si semplifica.
    """
    return mu_a * P_sca / (4.0 * math.pi * kappa)


def termica(R, lam, ambiente="B", z=Z, mu=None):
    """Tutti i termini per una cella (R, lambda), in un ambiente dichiarato."""
    env = AMBIENTI[ambiente]
    if mu is None:
        mu = c.SANGUE["mu_p"] if ambiente == "A" else 8.0e-3
    rif = c.ciclo(R, lam, z, mu=mu)
    if rif is None:
        return None
    P_inc = rif["P_inc"]
    mu_a = env["mu_a"][lam]
    mu_s = env["mu_s"][lam]
    kappa = env["kappa"]

    part = c.potenza_assorbita_particella(P_inc, R, lam, W0)
    dT_part = part["P_abs"] / (4.0 * math.pi * kappa * R)
    dT_p = dT_pre(P_inc, lam, z, mu_a, mu_s, kappa)
    dT_po = dT_post(P_inc, lam, z, mu_a, mu_s, kappa)
    ps = potenza_diffusa(P_inc, R, lam)
    dT_d = dT_diff(ps["P_sca"], mu_a, kappa)

    # Parete: immagine sulla particella e sul cammino, segno dichiarato.
    img_part = dT_part * R / (2.0 * z)
    dT_0 = dT_part + dT_p + dT_po + dT_d
    dT_sink = dT_part - img_part + dT_po + dT_d          # isoterma: annulla il cammino
    dT_mirror = dT_part + img_part + 2.0 * dT_p + dT_po + dT_d  # adiabatica: raddoppia

    return {
        "R_um": R * 1e6, "lambda_nm": lam, "ambiente": ambiente,
        "mu_mPa_s": mu * 1e3, "z_um": z * 1e6,
        "P_inc_mW": P_inc * 1e3, "T_s": trasmissione(mu_a, mu_s, z),
        "mu_a": mu_a, "mu_s": mu_s, "kappa": kappa,
        "P_abs_mW": part["P_abs"] * 1e3, "Q_abs": part["Q_abs"],
        "P_sca_mW": ps["P_sca"] * 1e3, "albedo": ps["albedo"],
        "eta_int": ps["eta_int"],
        "dT_part": dT_part, "dT_pre": dT_p, "dT_post": dT_po, "dT_diff": dT_d,
        "dT_0": dT_0, "dT_sink": dT_sink, "dT_mirror": dT_mirror,
        "sotto_5K_sink": dT_sink < D_T_CRITERIO,
        "sotto_5K_0": dT_0 < D_T_CRITERIO,
        "sotto_5K_mirror": dT_mirror < D_T_CRITERIO,
    }


def sintesi():
    righe = []
    for ambiente in ("A", "B"):
        for lam in LAM:
            for R in R_CASI:
                x = termica(R, lam, ambiente)
                if x:
                    righe.append(x)
    # riferimento V7.2: solo particle + pre, senza post/diff/wall
    rif = []
    for lam in LAM:
        for R in R_CASI:
            x = termica(R, lam, "B")
            if x:
                rif.append({"lambda_nm": lam, "R_um": x["R_um"],
                            "dT_v72": x["dT_part"] + x["dT_pre"],
                            "dT_v73_0": x["dT_0"],
                            "dT_v73_sink": x["dT_sink"],
                            "dT_v73_mirror": x["dT_mirror"]})
    return {
        "descrizione": "V7.3: termica a mezzo uniforme, tre contributi, "
                       "parete con segno esplicito",
        "criterio": {"dT_K": D_T_CRITERIO,
                     "nota": "criterio termico scelto, non prova di "
                             "sicurezza biologica"},
        "ambienti": {k: {"nome": v["nome"],
                         "mu_a": v["mu_a"], "mu_s": v["mu_s"],
                         "kappa": v["kappa"]}
                     for k, v in AMBIENTI.items()},
        "casi": righe,
        "confronto_v72": rif,
        "note": {
            "due_zone": "il modello a due zone (CFL di plasma + sangue) non "
                        "e' un ambiente V7.3: il CFL esiste solo con flusso, "
                        "e il flusso e' chiuso da v73_taglio.py",
            "parete": "sink = parete isoterma (vetro/endotelio), immagine "
                      "negativa, limite inferiore; mirror = parete "
                      "adiabatica, immagine positiva, limite superiore; "
                      "0 = nessuna correzione di parete",
            "post": "dT_post usa 2 pi (sorgente lineare) e z_R = pi w0^2/"
                    "lambda; con 4 pi il termine dimezza",
            "diff": "dT_diff = mu_a P_sca/(4 pi kappa), con P_sca limitata "
                    "da eta_int P_inc",
            "pennes": "non usato: tau_perfusione ~ 60 s contro T720 = 2.93 s, "
                      "peso 4.8%",
        },
    }


# ---------------------------------------------------------------------------
# 2. TABELLE
# ---------------------------------------------------------------------------
def tabella(lingua="it"):
    s = sintesi()
    righe = [x for x in s["casi"] if x["ambiente"] == "B"]
    if lingua == "it":
        intest = (r"$\lambda$ & $R$ & $P_{\rm inc}$ & $\Delta T_{\rm part}$ & "
                  r"$\Delta T_{\rm pre}$ & $\Delta T_{\rm post}$ & "
                  r"$\Delta T_{\rm diff}$ & $\Delta T_0$ & "
                  r"$\Delta T_{\rm sink}$ & $\Delta T_{\rm mir}$ \\"
                  "\n" r"nm & $\mu$m & mW & K & K & K & K & K & K & K \\")
        nota = (r"Nota: sangue intero in stasi (ramo B), $z=2\,\mu$m, "
                r"$w_0=0.5\,\mu$m.  $\Delta T_0$ somma i quattro contributi "
                r"senza correzione di parete; $\Delta T_{\rm sink}$ \`e il "
                r"limite isotermo (vetro/endotelio: immagine negativa, il "
                r"cammino viene annullato); $\Delta T_{\rm mir}$ (specchio) \`e il "
                r"limite adiabatico (immagine positiva, il cammino viene "
                r"raddoppiato).  Il criterio $5$~K \`e scelto, non biologico. "
                r"In plasma/tampone (ramo A) i termini di cammino, post e "
                r"diffuso sono $<0.01$~K: resta solo $\Delta T_{\rm part}$.")
    else:
        intest = (r"$\lambda$ & $R$ & $P_{\rm inc}$ & $\Delta T_{\rm part}$ & "
                  r"$\Delta T_{\rm pre}$ & $\Delta T_{\rm post}$ & "
                  r"$\Delta T_{\rm diff}$ & $\Delta T_0$ & "
                  r"$\Delta T_{\rm sink}$ & $\Delta T_{\rm mir}$ \\"
                  "\n" r"nm & $\mu$m & mW & K & K & K & K & K & K & K \\")
        nota = (r"Note: whole blood in stasis (branch B), $z=2\,\mu$m, "
                r"$w_0=0.5\,\mu$m.  $\Delta T_0$ sums the four terms without "
                r"a wall correction; $\Delta T_{\rm sink}$ is the isothermal "
                r"bound (glass/endothelium: negative image, the path term is "
                r"cancelled); $\Delta T_{\rm mir}$ (mirror) is the adiabatic bound "
                r"(positive image, the path term is doubled).  The $5$~K "
                r"criterion is chosen, not biological.  In plasma/buffer "
                r"(branch A) the path, post and diffuse terms are $<0.01$~K: "
                r"only $\Delta T_{\rm part}$ remains.")

    lines = [
        r"\begin{tabular}{@{}rrrrrrrrrr@{}}",
        r"\toprule",
        intest,
        r"\midrule",
    ]
    for lam in LAM:
        for R in R_CASI:
            x = next((q for q in righe
                      if q["lambda_nm"] == lam and abs(q["R_um"] - R * 1e6) < 1e-9),
                     None)
            if x is None:
                continue
            lines.append(
                r"%d & %.1f & %.3f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\"
                % (lam, x["R_um"], x["P_inc_mW"], x["dT_part"], x["dT_pre"],
                   x["dT_post"], x["dT_diff"], x["dT_0"], x["dT_sink"],
                   x["dT_mirror"]))
        lines.append(r"\addlinespace")
    if lines[-1] == r"\addlinespace":
        lines.pop()
    lines += [
        r"\midrule",
        r"\multicolumn{10}{p{14.6cm}}{\scriptsize " + nota + r"} \\",
        r"\bottomrule",
        r"\end{tabular}",
    ]
    return "\n".join(lines) + "\n"


def scrivi(nome, testo):
    with open(nome, "w", encoding="utf-8") as fh:
        fh.write(testo)
    print("  scritto %-52s %5.1f kB" % (nome, os.path.getsize(nome) / 1024.0))


def figura(lingua="it"):
    """ΔT per le tre pareti, a R=0.5 e R=0.3 um, sangue intero in stasi."""
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    s = sintesi()
    it = lingua == "it"
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.5), sharey=True)
    for ax, R in zip(axes, (0.5, 0.3)):
        x = [q for q in s["casi"]
             if q["ambiente"] == "B" and abs(q["R_um"] - R) < 1e-9]
        x = sorted(x, key=lambda q: q["lambda_nm"])
        pos = np.arange(len(LAM))
        w = 0.25
        ax.bar(pos - w, [q["dT_sink"] for q in x], w, color="#2166ac",
               label=("parete isoterma (sink)" if it else "isothermal wall (sink)"))
        ax.bar(pos, [q["dT_0"] for q in x], w, color="#777777",
               label=("nessuna correzione" if it else "no wall correction"))
        ax.bar(pos + w, [q["dT_mirror"] for q in x], w, color="#b2182b",
               label=("parete adiabatica (mirror)" if it else "adiabatic wall (mirror)"))
        ax.axhline(5.0, ls="--", color="#b2182b", lw=1.1)
        ax.text(2.42, 5.4, "5 K", color="#b2182b", ha="right", va="bottom",
                fontsize=8)
        ax.set_xticks(pos)
        ax.set_xticklabels(["532", "800", "1064"])
        ax.set_xlabel(("lunghezza d'onda (nm)" if it else "wavelength (nm)"))
        ax.set_title(("$R=%.1f\\,\\mu$m" % R), fontsize=9.5)
        ax.set_yscale("log")
        ax.set_ylim(0.15, 120)
        ax.grid(True, axis="y", alpha=0.22, lw=0.5)
    axes[0].set_ylabel(("$\\Delta T$ (K), $z=2\\,\\mu$m" if it else
                        "$\\Delta T$ (K), $z=2\\,\\mu$m"))
    axes[0].legend(fontsize=7.5, loc="upper right")
    fig.suptitle(("Sangue intero in stasi: effetto della parete" if it else
                  "Whole blood in stasis: wall effect"), fontsize=10)
    fig.tight_layout()
    for d in (FIG_IT, FIG_ROOT_IT) if it else (FIG_EN, FIG_ROOT_EN):
        os.makedirs(d, exist_ok=True)
        p = os.path.join(d, "v73_termica_%s.pdf" % lingua)
        fig.savefig(p, bbox_inches="tight")
        fig.savefig(os.path.join(d, "v73_termica_%s.png" % lingua),
                    bbox_inches="tight", dpi=220)
        print("  scritto %-52s %5.1f kB" % (p, os.path.getsize(p) / 1024.0))
    plt.close(fig)


def main():
    s = sintesi()
    p = os.path.join(DATI, "V73_termica.json")
    with open(p, "w", encoding="utf-8") as fh:
        json.dump(s, fh, indent=1, ensure_ascii=False, default=float)
    print("  scritto %-52s %5.1f kB" % (p, os.path.getsize(p) / 1024.0))
    scrivi(os.path.join(TAB_IT, "tab_v73_termica.tex"), tabella("it"))
    scrivi(os.path.join(TAB_EN, "tab_v73_termica.tex"), tabella("en"))
    figura("it")
    figura("en")

    print("\n-- sangue intero in stasi, z=2 um")
    for lam in LAM:
        print("   %4d nm" % lam)
        for R in R_CASI:
            x = termica(R, lam, "B")
            print("      R=%.1f um  P=%.3f mW  part=%.2f pre=%.2f post=%.2f "
                  "diff=%.2f  0=%.2f sink=%.2f mirror=%.2f"
                  % (R * 1e6, x["P_inc_mW"], x["dT_part"], x["dT_pre"],
                     x["dT_post"], x["dT_diff"], x["dT_0"], x["dT_sink"],
                     x["dT_mirror"]))
    print("\n-- plasma/tampone in vitro, z=2 um")
    for lam in LAM:
        x = termica(0.5e-6, lam, "A")
        print("   %4d nm  P=%.3f mW  part=%.2f pre=%.4f post=%.4f diff=%.4f"
              % (lam, x["P_inc_mW"], x["dT_part"], x["dT_pre"],
                 x["dT_post"], x["dT_diff"]))
    print("\n-- confronto V7.2 (solo part+pre) e V7.3, R=0.5 um")
    for x in s["confronto_v72"]:
        if abs(x["R_um"] - 0.5) < 1e-9:
            print("   %4d nm  V72=%.2f K  V73_0=%.2f  sink=%.2f  mirror=%.2f"
                  % (x["lambda_nm"], x["dT_v72"], x["dT_v73_0"],
                     x["dT_v73_sink"], x["dT_v73_mirror"]))
    print("\n   criterio: dT < %.0f K (scelto, non biologico)" % D_T_CRITERIO)
    return s


if __name__ == "__main__":
    main()
