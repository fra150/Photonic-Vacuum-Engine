# -*- coding: utf-8 -*-
"""
v73_test.py - controlli di coerenza della V7.3.

Non e' una validazione fisica: e' una suite di regressione che verifica che
i numeri pubblicati siano quelli che gli script producono, che le identita'
del cancello valgano, che la forcella termica sia ordinata e che il rumore
esatto non sia stato sostituito di nuovo dall'approssimazione gaussiana.

Uso:
    python codice/v73_test.py

Esce con codice 1 al primo controllo fallito.  Non richiede Blender ne'
LaTeX: i file prodotti da Blender e dalle tabelle vengono solo ispezionati.
"""

import json
import math
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(BASE, ".."))
DATI = os.path.join(ROOT, "dati")
PAPER = os.path.join(ROOT, "paper")
if BASE not in sys.path:
    sys.path.insert(0, BASE)

import v73_contrappeso as c  # noqa: E402
import v73_rumore as n       # noqa: E402
import v73_taglio as t       # noqa: E402

OK = [0]


def check(cond, msg):
    if cond:
        OK[0] += 1
        print("  ok   %s" % msg)
    else:
        print("  FAIL %s" % msg)
        raise AssertionError(msg)


def load(nome):
    with open(os.path.join(DATI, nome), encoding="utf-8") as fh:
        return json.load(fh)


def vicino(a, b, tol):
    return abs(a - b) <= tol * max(1.0, abs(b))


TG = load("V73_taglio.json")
TH = load("V73_termica.json")
RU = load("V73_rumore.json")
AM = load("V73_ambienti.json")
SI = load("V73_sintesi.json")
GR = load("V73_griglia.json")


# ---------------------------------------------------------------------------
# 1. cancello: le quattro uscite e il minimo
# ---------------------------------------------------------------------------
def test_cancello():
    print("\n[1] cancello sul flusso")
    g = t.gamma_coppia()
    check(vicino(g["gamma_vincolato"], 2.7002050307921315, 1e-12),
          "gamma_coppia = 2.700205 s^-1")
    check(vicino(g["gamma_sfera_libera"], 20.959540011604155, 1e-12),
          "sfera libera = 20.959540 s^-1 (verifica di convenzione)")
    check(vicino(g["gamma_vincolato"],
                 t.TAU_OTT / t.beta_cw(t.R_REF, t.MU_PLASMA), 1e-12),
          "gamma_coppia = tau_opt / beta_cw")

    for x in TG["output_2_forza"]:
        f = t.gamma_forza(Q=x["Q"])
        check(vicino(f["gamma"], x["gamma_s-1"], 1e-12),
              "gamma_forza(Q=%.1f) = F/(6 pi mu R_eff z_c)" % x["Q"])
        check(vicino(f["F_pN"], x["F_pN"], 1e-12),
              "F_trap(Q=%.1f) = Q n P / c" % x["Q"])

    minimo = min([g["gamma_vincolato"]]
                 + [x["gamma_s-1"] for x in TG["output_2_forza"]])
    check(vicino(TG["gamma_limite_s-1"], minimo, 1e-12),
          "gamma_limite = min(output 1, output 2)")
    check(vicino(TG["gamma_limite_s-1"], 2.7002050307921315, 1e-12),
          "gamma_limite = 2.700 s^-1")
    check(TG["gamma_limite_s-1"] < TG["output_3_jeffery"]["gamma_a_T720_s-1"],
          "Jeffery 5.20 s^-1 e' sopra il limite, non lo definisce")

    # Bingham additivo: mu = 8 mPa s, sigma_y = 10 mPa
    gb = t.gamma_coppia(mu=8.0e-3, sigma_y=10.0e-3)
    check(vicino(gb["gamma_vincolato"], 0.27774656351859234, 1e-12),
          "Bingham sangue (8 mPa s, 10 mPa) = 0.277747 s^-1")
    check(vicino(gb["drive"], t.TAU_OTT - gb["M_y"], 1e-18),
          "M_y e' uno scostamento additivo, non un fattore")

    # criterio di parete: 4 pi R^3 sigma_y = tau
    for x in TG["criterio_parete"]:
        my = 4.0 * math.pi * (x["R_um"] * 1e-6) ** 3 * x["sigma_y_mPa"] * 1e-3
        check(vicino(my, t.TAU_OTT, 0.02),
              "parete: 4 pi R^3 sigma_y = tau a R=%.2f um" % x["R_um"])

    corr = TG["correzione_idrodinamica"]
    check(vicino(corr["termine"], (corr["r_cw"] / 2.0) ** 2, 1e-12),
          "correzione due sfere: (r_cw/d)^2 = 0.15749")
    check(vicino(corr["ratio_bead"], 4.2561015779523, 1e-12),
          "ratio bead = 4.2561 (+9.7%)")
    check(vicino(corr["ratio_bead_totale"], 4.33484664357073, 1e-12),
          "ratio bead totale = 4.3348 (+11.7%)")
    check(vicino(corr["ratio_limite_superiore"], 4.492336774807589, 1e-12),
          "limite superiore = 4.4923 (+15.7%)")


# ---------------------------------------------------------------------------
# 2. ambienti
# ---------------------------------------------------------------------------
def test_ambienti():
    print("\n[2] ambienti")
    stati = {x["id"]: x["stato"] for x in AM["ambienti"]}
    check(stati["A"] == "ammesso", "A (plasma in vitro, stasi) ammesso")
    check(stati["B"] == "ammesso con riserva",
          "B (sangue in vivo, stasi) ammesso con riserva")
    check(stati["C"] == "escluso", "C (sangue in vivo, flusso) escluso")
    check(stati["D"] == "ammesso solo sotto soglia",
          "D (microfluidica ex vivo) ammesso solo sotto soglia")
    check(vicino(AM["gate"]["gamma_limite"], TG["gamma_limite_s-1"], 1e-12),
          "la tabella ambienti usa lo stesso gamma_limite del cancello")
    check("solo in vitro" in AM["note"]["plasma_stasi"],
          "plasma + stasi dichiarato solo in vitro")


# ---------------------------------------------------------------------------
# 3. termica: forcella e criterio
# ---------------------------------------------------------------------------
def test_termica():
    print("\n[3] termica a mezzo uniforme")
    check(vicino(TH["criterio"]["dT_K"], 5.0, 1e-12),
          "criterio = 5 K")
    check("non prova di sicurezza" in TH["criterio"]["nota"],
          "il criterio 5 K e' dichiarato scelto, non prova di sicurezza")

    for q in TH["casi"]:
        somma = q["dT_part"] + q["dT_pre"] + q["dT_post"] + q["dT_diff"]
        check(vicino(somma, q["dT_0"], 1e-9),
              "%s %d nm R=%.1f: dT_0 = somma dei quattro contributi"
              % (q["ambiente"], q["lambda_nm"], q["R_um"]))
        check(q["dT_sink"] <= q["dT_0"] + 1e-9,
              "%s %d nm R=%.1f: sink <= dT_0"
              % (q["ambiente"], q["lambda_nm"], q["R_um"]))
        check(q["dT_mirror"] >= q["dT_0"] - 1e-9,
              "%s %d nm R=%.1f: mirror >= dT_0"
              % (q["ambiente"], q["lambda_nm"], q["R_um"]))
        img = q["dT_part"] * q["R_um"] / (2.0 * q["z_um"])
        check(vicino(q["dT_sink"], q["dT_0"] - q["dT_pre"] - img, 1e-9),
              "%s %d nm R=%.1f: sink annulla il cammino e sottrae l'immagine"
              % (q["ambiente"], q["lambda_nm"], q["R_um"]))
        check(vicino(q["dT_mirror"], q["dT_0"] + q["dT_pre"] + img, 1e-9),
              "%s %d nm R=%.1f: mirror raddoppia il cammino e somma l'immagine"
              % (q["ambiente"], q["lambda_nm"], q["R_um"]))
        check(q["P_sca_mW"] <= q["eta_int"] * q["P_inc_mW"] * (1.0 + 1e-12),
              "%s %d nm R=%.1f: P_sca limitata da eta_int P_inc"
              % (q["ambiente"], q["lambda_nm"], q["R_um"]))

    def th(env, lam, R):
        for q in TH["casi"]:
            if (q["ambiente"] == env and q["lambda_nm"] == lam
                    and abs(q["R_um"] - R) < 1e-9
                    and abs(q["z_um"] - 2.0) < 1e-9):
                return q
        raise KeyError((env, lam, R))

    check(th("B", 532, 0.5)["dT_sink"] > 5.0,
          "532 nm, sangue, R=0.5: fallisce anche nel pozzo (5.72 K)")
    check(th("B", 532, 0.3)["dT_sink"] > 5.0,
          "532 nm, sangue, R=0.3: fallisce anche nel pozzo (22.25 K)")
    for lam in (800, 1064):
        for R in (0.3, 0.5):
            check(th("B", lam, R)["dT_mirror"] < 5.0,
                  "%d nm, sangue, R=%.1f: passa anche nello specchio (%.2f K)"
                  % (lam, R, th("B", lam, R)["dT_mirror"]))
    check(th("A", 532, 0.5)["dT_0"] < 1.0,
          "plasma, 532 nm, R=0.5: contributo di percorso trascurabile")


# ---------------------------------------------------------------------------
# 4. rumore esatto e latch
# ---------------------------------------------------------------------------
def test_rumore():
    print("\n[4] rumore esatto e latch")
    cv = RU["input"]["CV_ideale"]
    e = n.errore_appaiato(cv)
    check(vicino(e["eps"], 1.0093642171805574e-4, 1e-9),
          "errore appaiato a CV ideale = 1.01e-4")
    check(vicino(e["eps_geo"], 1.071763598904871e-4, 1e-9),
          "errore appaiato geometrico = 1.07e-4")
    check(vicino(n.errore_totale(cv)["eps"], 1.1372352512035116e-2, 1e-9),
          "errore a durata totale = 1.14e-2")

    base = n.errore_appaiato(cv)["eps"]
    for N in (1, 2, 6):
        check(vicino(n.errore_appaiato(cv, N=N)["eps"], base, 1e-12),
              "l'errore appaiato non dipende da N (N=%d)" % N)

    lim = n.limite_cv()
    check(vicino(lim["appaiato"], 0.17830789239388847, 1e-10),
          "limite CV appaiato = 17.83%")
    check(vicino(lim["totale"], 0.1158067007739147, 1e-10),
          "limite CV totale = 11.58%")
    check(vicino(RU["limite_cv"]["appaiato"], lim["appaiato"], 1e-12),
          "il JSON riporta il limite appaiato corretto")

    eps = [n.errore_appaiato(x)["eps"] for x in (0.15, 0.178, 0.2, 0.25, 0.30)]
    check(all(eps[i] < eps[i + 1] for i in range(len(eps) - 1)),
          "l'errore cresce con CV")
    check(vicino(n.errore_appaiato(0.302)["eps"], 2.0708246218485975e-2, 1e-9),
          "a CV=30.2% l'errore esatto e' 2.07e-2, non 1e-4")

    sp = RU["specifica_latch"]
    check(vicino(sp["cv"]["limite_10e4_appaiato"], lim["appaiato"], 1e-12),
          "la specifica del latch usa il limite esatto")
    check(sp["cv"]["target"] == 0.3 and "NON e' supportato" in sp["cv"]["stato"],
          "il target CV 30% e' dichiarato non supportato")
    check(sp["ratio"]["valore"] < sp["ratio"]["corretto_bead"]
          < sp["ratio"]["corretto_bead_tot"]
          < sp["ratio"]["limite_superiore"],
          "i quattro rapporti sono ordinati (3.881 < 4.256 < 4.335 < 4.492)")
    check("non e' dimostrato" in RU["verdetto"],
          "il verdetto del latch non lo dichiara dimostrato")

    for k in ("eps_appaiato", "eps_totale"):
        check(RU["ritrattazione_30_2"][k] > 1e-3,
              "la ritrattazione del 30.2%% riporta un errore > 1e-3 (%s)" % k)


# ---------------------------------------------------------------------------
# 5. coerenza del modello e dei JSON
# ---------------------------------------------------------------------------
def test_modello():
    print("\n[5] modello e JSON")
    vs = c.verifica_scaling()
    check(vs["lambda_ok"], "scaling lambda verificato")
    check(vs["dn_ok"], "scaling Delta_n verificato")
    check(vicino(vs["lambda_1064_su_800"], 0.762765646549912, 1e-12),
          "lambda 1064/800 = 0.7628")
    check(vicino(vs["dn_010_su_006"], 2.7631782596346244, 1e-12),
          "Delta_n 0.10/0.06 = 2.763")
    check(vicino(SI["beta_ratio"], 3.8811015779522995, 1e-12),
          "beta_cw/beta_no = 3.8811")
    check("4 pi R^3" in TG["convenzioni"]["bingham"]
          and "(8/3)" in TG["convenzioni"]["bingham"],
          "la convenzione Bingham dichiara il fattore 1.5 della soluzione esatta")
    check("solo con flusso" in TG["convenzioni"]["cfl"],
          "il CFL e' dichiarato limite del dominio, non ambiente")
    check("rumore_gaussiano_ritrattato" in SI,
          "il JSON segnala i campi gaussiani ritrattati")
    check(vicino(SI["riassunto"]["celle"], 540, 0),
          "la griglia ha 540 celle")
    check(vicino(SI["riassunto"]["valide"], 127, 0),
          "127 celle valide (dT < 20 K)")
    check(vicino(SI["riassunto"]["valide_sicure"], 93, 0),
          "93 celle sotto il criterio 5 K")
    check(vicino(SI["riassunto"]["valide_cfl"], 133, 0),
          "133 celle con CFL soddisfatto")
    check(len(GR["celle"]) == 540,
          "V73_griglia.json contiene 540 celle")


# ---------------------------------------------------------------------------
# 6. tabelle, macro e figure
# ---------------------------------------------------------------------------
def test_prodotti():
    print("\n[6] tabelle, macro e figure")
    macro = os.path.join(PAPER, "numeri_v73.tex")
    check(os.path.isfile(macro), "paper/numeri_v73.tex esiste")
    with open(macro, encoding="utf-8") as fh:
        testo = fh.read()
    for nome, valore in (("VSeventyThreeGammaCrit", "2.700"),
                         ("VSeventyThreeCVLimPair", "17.83"),
                         ("VSeventyThreeCVLimBeadTot", "19.4"),
                         ("VSeventyThreeTCriterio", "5.0")):
        check(("\\newcommand{\\%s}{%s}" % (nome, valore)) in testo,
              "macro %s = %s" % (nome, valore))

    tabelle = [
        "tab_v73_ambienti.tex", "tab_v73_termica.tex", "tab_v73_rumore.tex",
        "tab_v73_tempi.tex", "tab_v73_latch.tex", "tab_v73_finestre.tex",
    ]
    for lang in ("it", "en"):
        for nome in tabelle:
            p = os.path.join(PAPER, lang, "tabelle", nome)
            check(os.path.isfile(p) and os.path.getsize(p) > 100,
                  "tabella %s/%s non vuota" % (lang, nome))

    with open(os.path.join(PAPER, "it", "tabelle",
                           "tab_v73_finestre.tex"), encoding="utf-8") as fh:
        fin = fh.read()
    check("sicuro" not in fin, "la tabella finestre non usa piu' 'sicuro'")
    check("sotto 5 K" in fin, "la tabella finestre dice 'sotto 5 K'")

    figure = ["fig_v73_cancello.pdf", "fig_v73_termica.pdf",
              "fig_v73_profondita.pdf", "fig_v73_potenza.pdf",
              "fig_v73_dt.pdf", "fig_v73_tempo.pdf", "fig_v73_mappa.pdf",
              "fig_v73_fuga.pdf", "fig_v73_zeri.pdf"]
    for lang in ("it", "en"):
        for nome in figure:
            p = os.path.join(PAPER, lang, "figures", nome)
            check(os.path.isfile(p) and os.path.getsize(p) > 1000,
                  "figura %s/%s non vuota" % (lang, nome))
        for nome in ("v73_rumore_%s.pdf" % lang, "v73_termica_%s.pdf" % lang):
            p = os.path.join(PAPER, lang, "figures", nome)
            check(os.path.isfile(p) and os.path.getsize(p) > 1000,
                  "figura %s/%s non vuota" % (lang, nome))

    blender = ["v73_latch_000.png", "v73_latch_360.png",
               "v73_latch_720.png", "v73_latch_dettaglio.png",
               "v73_cancello.png", "motore_fotonico_v73.blend"]
    for nome in blender:
        p = os.path.join(ROOT, "blender", nome)
        check(os.path.isfile(p) and os.path.getsize(p) > 10000,
              "blender/%s non vuoto" % nome)


def main():
    print("V7.3 - controlli di coerenza")
    test_cancello()
    test_ambienti()
    test_termica()
    test_rumore()
    test_modello()
    test_prodotti()
    print("\nTutti i %d controlli sono passati." % OK[0])


if __name__ == "__main__":
    main()
