"""
MOTORE FOTONICO DEL VUOTO - V7.3
v73_ambienti.py - la tabella degli ambienti, generata dal cancello.

La tabella non e' un elenco di possibilita': e' la conseguenza di
v73_taglio.py.  Ogni ambiente dichiara mu, reologia, CFL, cellule, parete
termica e mu_a; l'ultima colonna riporta il limite di shear e lo stato.

Regole fisse:
  - il CFL e' migrazione indotta dal flusso: in stasi non esiste;
  - "plasma + stasi" si realizza solo in vitro;
  - la parete e' un pozzo termico (vetro, endotelio) salvo dichiarazione
    esplicita di parete adiabatica (specchio, fattore 2);
  - il guadagno CFL sangue/plasma decresce con lambda ma non svanisce:
    5.3e5 (532), 190 (800), 13 (1064).

Uso:
    python codice/v73_ambienti.py
"""

import json
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(BASE, ".."))
if BASE not in sys.path:
    sys.path.insert(0, BASE)

import v73_taglio as t  # noqa: E402

DATI = os.path.join(ROOT, "dati")
TAB_IT = os.path.join(ROOT, "paper", "it", "tabelle")
TAB_EN = os.path.join(ROOT, "paper", "en", "tabelle")
for d in (DATI, TAB_IT, TAB_EN):
    os.makedirs(d, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. I QUATTRO AMBIENTI
# ---------------------------------------------------------------------------
def ambienti():
    g_plasma = t.gamma_coppia(mu=t.MU_PLASMA)
    g_blood = [t.gamma_coppia(mu=mu, sigma_y=sy)
               for mu, sy in zip(t.MU_BLOOD, t.SIGMA_Y_BLOOD)]
    return [
        {
            "id": "A",
            "nome_it": "Plasma o tampone, in vitro, stasi",
            "nome_en": "Plasma/buffer, in vitro, stasis",
            "mu_mPa_s": "1.3--1.5",
            "reologia_it": "newtoniana; $\\sigma_y=0$",
            "reologia_en": "Newtonian; $\\sigma_y=0$",
            "cfl_it": "assente; nessuna cellula",
            "cfl_en": "none; no cells",
            "parete_it": "vetro o PDMS, pozzo ($-$)",
            "parete_en": "glass or PDMS, sink ($-$)",
            "mu_a": "0.044 / 2.3 / 15",
            "mu_a_en": "0.044 / 2.3 / 15",
            "gamma_it": "$\\dot\\gamma_{\\rm lim}=2.70$ s$^{-1}$; ammesso",
            "gamma_en": "$\\dot\\gamma_{\\rm lim}=2.70$ s$^{-1}$; allowed",
            "stato": "ammesso",
        },
        {
            "id": "B",
            "nome_it": "Sangue intero, in vivo, stasi",
            "nome_en": "Whole blood, in vivo, stasis",
            "mu_mPa_s": "7--9",
            "reologia_it": "Bingham; $\\sigma_y=5$--20 mPa",
            "reologia_en": "Bingham; $\\sigma_y=5$--20 mPa",
            "cfl_it": "assente; RBC, rouleaux",
            "cfl_en": "none; RBCs, rouleaux",
            "parete_it": "endotelio, pozzo ($-$)",
            "parete_en": "endothelium, sink ($-$)",
            "mu_a": "23500 / 437 / 200",
            "mu_a_en": "23500 / 437 / 200",
            "gamma_it": "$\\dot\\gamma_{\\rm lim}=0.10$--$0.41$; ammesso con riserva",
            "gamma_en": "$\\dot\\gamma_{\\rm lim}=0.10$--$0.41$; allowed with caveats",
            "stato": "ammesso con riserva",
        },
        {
            "id": "C",
            "nome_it": "Sangue intero, in vivo, flusso",
            "nome_en": "Whole blood, in vivo, flow",
            "mu_mPa_s": "4--10",
            "reologia_it": "Bingham; $\\sigma_y=5$--10 mPa",
            "reologia_en": "Bingham; $\\sigma_y=5$--10 mPa",
            "cfl_it": "1--3 $\\mu$m; collisioni RBC",
            "cfl_en": "1--3 $\\mu$m; RBC collisions",
            "parete_it": "endotelio, pozzo ($-$)",
            "parete_en": "endothelium, sink ($-$)",
            "mu_a": "23500 / 437 / 200",
            "mu_a_en": "23500 / 437 / 200",
            "gamma_it": "$\\dot\\gamma=10^2$--$10^3$; escluso",
            "gamma_en": "$\\dot\\gamma=10^2$--$10^3$; excluded",
            "stato": "escluso",
        },
        {
            "id": "D",
            "nome_it": "Microfluidica, ex vivo, shear controllato",
            "nome_en": "Microfluidics, ex vivo, controlled shear",
            "mu_mPa_s": "1.3--10",
            "reologia_it": "newtoniana o Bingham",
            "reologia_en": "Newtonian or Bingham",
            "cfl_it": "solo se cellule e flusso",
            "cfl_en": "only with cells and flow",
            "parete_it": "vetro o PDMS, pozzo ($-$)",
            "parete_en": "glass or PDMS, sink ($-$)",
            "mu_a": "secondo mezzo",
            "mu_a_en": "depends on medium",
            "gamma_it": "$\\dot\\gamma<\\dot\\gamma_{\\rm lim}$; ammesso solo sotto soglia",
            "gamma_en": "$\\dot\\gamma<\\dot\\gamma_{\\rm lim}$; allowed only below threshold",
            "stato": "ammesso solo sotto soglia",
        },
    ]


def sintesi():
    s = t.sintesi()
    return {
        "descrizione": "V7.3: ambienti coerenti con il cancello sul flusso",
        "gate": {
            "gamma_coppia_plasma": s["output_1_coppia"]["gamma_vincolato_s-1"],
            "gamma_coppia_sangue": {
                "7mPa_5mPa": t.gamma_coppia(mu=7e-3, sigma_y=5e-3)["gamma_vincolato"],
                "8mPa_10mPa": t.gamma_coppia(mu=8e-3, sigma_y=10e-3)["gamma_vincolato"],
                "9mPa_20mPa": t.gamma_coppia(mu=9e-3, sigma_y=20e-3)["gamma_vincolato"],
            },
            "gamma_forza": [x["gamma_s-1"] for x in s["output_2_forza"]],
            "gamma_jeffery": s["output_3_jeffery"]["gamma_a_T720_s-1"],
            "gamma_limite": s["gamma_limite_s-1"],
            "rapporti_fisiologici": s["confronto_fisiologico"],
        },
        "ottica": {
            "mu_a_m-1": {"532": {"plasma": 0.044, "sangue": 23500},
                         "800": {"plasma": 2.3, "sangue": 437},
                         "1064": {"plasma": 15.0, "sangue": 200}},
            "guadagno_cfl": {"532": 5.3e5, "800": 190, "1064": 13},
            "nota": "il guadagno CFL decresce con lambda ma resta 13x a 1064; "
                    "il CFL esiste solo con flusso, quindi non entra negli "
                    "ambienti A, B, D",
        },
        "termica": {
            "parete": {
                "vetro_PDMS": "pozzo (immagine negativa)",
                "endotelio": "pozzo; L_diff = 640 um a 3 s, il calore esce dal vaso",
                "adiabatica": "specchio (fattore 2); va dichiarata",
            },
            "pennes": "non usato: tau_perfusione ~ 60 s contro T720 = 2.93 s, "
                      "peso 4.8%",
        },
        "ambienti": ambienti(),
        "note": {
            "cfl": "limite del dominio, non ambiente: esiste solo con flusso",
            "plasma_stasi": "plasma + stasi si realizza solo in vitro",
            "verdetto": "il flusso fisiologico e' 99-741 volte sopra la soglia; "
                        "il dispositivo vive in stasi ex vivo o in stasi in vivo",
        },
    }


# ---------------------------------------------------------------------------
# 2. LATEX
# ---------------------------------------------------------------------------
def tabella(lingua="it"):
    A = ambienti()
    if lingua == "it":
        intest = ("Ambiente & $\\mu$ (mPa\\,s) & Reologia & CFL / cellule & "
                  "Parete termica & $\\mu_a$ (m$^{-1}$)\\newline 532 / 800 / 1064 & "
                  "$\\dot\\gamma$ / stato \\\\")
        nome = lambda x: x["nome_it"]
        reo = lambda x: x["reologia_it"]
        cfl = lambda x: x["cfl_it"]
        par = lambda x: x["parete_it"]
        gam = lambda x: x["gamma_it"]
        nota = ("Nota: il CFL \\`e migrazione indotta dal flusso e non esiste "
                "in stasi; ``plasma + stasi'' si realizza solo in vitro.  "
                "Guadagno CFL sangue/plasma: $5.3\\times10^5$ (532~nm), 190 "
                "(800~nm), 13 (1064~nm); decresce con $\\lambda$ ma non "
                "svanisce.  Il segno della sorgente immagine dipende dalla "
                "parete: vetro/endotelio = pozzo (immagine negativa); parete "
                "adiabatica = specchio (fattore 2).")
    else:
        intest = ("Environment & $\\mu$ (mPa\\,s) & Rheology & CFL / cells & "
                  "Thermal wall & $\\mu_a$ (m$^{-1}$)\\newline 532 / 800 / 1064 & "
                  "$\\dot\\gamma$ / status \\\\")
        nome = lambda x: x["nome_en"]
        reo = lambda x: x["reologia_en"]
        cfl = lambda x: x["cfl_en"]
        par = lambda x: x["parete_en"]
        gam = lambda x: x["gamma_en"]
        nota = ("Note: the CFL is flow-induced migration and does not exist in "
                "stasis; ``plasma + stasis'' is realized only in vitro.  CFL "
                "gain blood/plasma: $5.3\\times10^5$ (532~nm), 190 (800~nm), "
                "13 (1064~nm); it decreases with $\\lambda$ but does not "
                "vanish.  The sign of the image source depends on the wall: "
                "glass/endothelium = sink (negative image); adiabatic wall = "
                "mirror (factor 2).")

    lines = [
        r"\setlength{\tabcolsep}{4pt}",
        r"\begin{tabular}{@{}>{\raggedright\arraybackslash}p{2.2cm}"
        r">{\raggedright\arraybackslash}p{1.25cm}"
        r">{\raggedright\arraybackslash}p{2.0cm}"
        r">{\raggedright\arraybackslash}p{1.95cm}"
        r">{\raggedright\arraybackslash}p{2.05cm}"
        r">{\raggedright\arraybackslash}p{2.3cm}"
        r">{\raggedright\arraybackslash}p{2.3cm}@{}}",
        r"\toprule",
        intest,
        r"\midrule",
    ]
    for x in A:
        lines.append(r"\textbf{%s.} %s & %s & %s & %s & %s & %s & %s \\"
                     % (x["id"], nome(x), x["mu_mPa_s"], reo(x), cfl(x),
                        par(x), (x["mu_a"] if lingua == "it" else x["mu_a_en"]),
                        gam(x)))
        lines.append(r"\addlinespace")
    if lines[-1] == r"\addlinespace":
        lines.pop()
    lines += [
        r"\midrule",
        r"\multicolumn{7}{p{14.6cm}}{\scriptsize " + nota + r"} \\",
        r"\bottomrule",
        r"\end{tabular}",
    ]
    return "\n".join(lines) + "\n"


def scrivi(nome, testo):
    p = os.path.join(nome)
    with open(p, "w", encoding="utf-8") as fh:
        fh.write(testo)
    print("  scritto %-46s %5.1f kB" % (p, os.path.getsize(p) / 1024.0))


def main():
    s = sintesi()
    p = os.path.join(DATI, "V73_ambienti.json")
    with open(p, "w", encoding="utf-8") as fh:
        json.dump(s, fh, indent=1, ensure_ascii=False, default=float)
    print("  scritto %-46s %5.1f kB" % (p, os.path.getsize(p) / 1024.0))
    scrivi(os.path.join(TAB_IT, "tab_v73_ambienti.tex"), tabella("it"))
    scrivi(os.path.join(TAB_EN, "tab_v73_ambienti.tex"), tabella("en"))

    print("\n-- gate")
    g = s["gate"]
    print("   plasma: gamma_coppia = %.4f s^-1" % g["gamma_coppia_plasma"])
    for k, v in g["gamma_coppia_sangue"].items():
        print("   sangue %-14s %.4f s^-1" % (k, v))
    print("   forza: %s s^-1" % ["%.2f" % x for x in g["gamma_forza"]])
    print("   Jeffery: %.2f s^-1" % g["gamma_jeffery"])
    print("   limite: %.2f s^-1" % g["gamma_limite"])
    for k, v in g["rapporti_fisiologici"].items():
        print("   %-14s x%.1f" % (k, v["rapporto"]))
    print("\n-- ambienti")
    for x in s["ambienti"]:
        print("   %s  %-38s mu=%-7s  %s"
              % (x["id"], x["nome_it"], x["mu_mPa_s"], x["stato"]))
    print("\n   %s" % s["note"]["plasma_stasi"])
    return s


if __name__ == "__main__":
    main()
