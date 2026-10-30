"""
MOTORE FOTONICO DEL VUOTO - V7
v8_validazione.py - V7 - CONSISTENCY TESTS (not validation)

Test di coerenza interna, non validazione sperimentale del modello.  Se un
test fallisce, il modello si corregge, non si aggiusta il dato.

Test
----
A  algebra del ritardo di fase: quattro casi limite chiusi a mano
B  limite di Beth: lamina lambda/2 in luce circolare -> 2 hbar per fotone
C  termica in acqua contro il conto a mano della V6.1 (150 uW, R = 1 e 3 um)
D  He 1995: CuO assorbente con OAM - controllo di PUNTO (x3); FALLISCE
E  Friese 1998: calcite birifrangente - controllo di PUNTO (x3)
F  limite GHz: manubrio di silice da 100 nm, con il tempo di avvio corretto

Le tolleranze sono dichiarate: fattore 3 sui valori.  I test D ed E usano i
parametri pubblicati e un solo numero misurato, con il rapporto modello/misura
riportato.  I parametri sono raccolti in PARAMETRI_LETTERATURA.

Nota su D: il test usa 632.8 nm, l = 3, 10 mW e una particella assorbente di
~1 um di raggio in kerosene.  Con la Mie locale a 632.8 nm il modello da'
~9.8 volte i 4 Hz misurati: D e' un FAIL dichiarato.  Non e' una validazione
del modello di ritardo (a 632.8 nm e con CuO assorbente il canale e' quello
di assorbimento), ed e' esattamente il tipo di scarto che la V8 riporta senza
nasconderlo.

Uso:
    python v8_validazione.py
"""

import json
import math
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
if BASE not in sys.path:
    sys.path.insert(0, BASE)

import v7_mappa_corretta as m  # noqa: E402

DATI = os.path.normpath(os.path.join(BASE, "..", "dati_v7"))
os.makedirs(DATI, exist_ok=True)

TOL_VALORE = 3.0          # fattore ammesso sui valori

# Parametri presi dagli articoli.  Se un numero va corretto, si corregge qui:
# ogni test legge solo da questo dizionario.
PARAMETRI_LETTERATURA = {
    # He, Friese, Heckenberg, Rubinsztein-Dunlop, PRL 75, 826 (1995).
    # HeNe 632.8 nm, 10 mW, carica OAM l = 3, particella CuO assorbente di
    # ~1 um di raggio in kerosene, 4 Hz.  Il waist focalizzato non e'
    # riportato nella fonte secondaria: si usa 1.2 um, coerente con una
    # particella da 2 um di diametro che riempie il fuoco.  Il FAIL non
    # dipende da questa scelta: con w0 = 0.5 um il rapporto sale a ~56.
    "he1995": {"R_m": 1.0e-6, "P_W": 10e-3, "l": 3, "s": 0.0,
               "mu_Pas": 1.58e-3, "f_mis_Hz": 4.0, "materiale": "cu_nero",
               "lambda_m": 632.8e-9, "w0_m": 1.2e-6},
    # Friese, Nieminen, Heckenberg, Rubinsztein-Dunlop, Nature 394, 348 (1998).
    # Particella di calcite, 1 um di spessore, 300 mW a 1064 nm, 357 Hz.
    # Il raggio laterale non e' riportato: si assume 1.5 um (frammento di
    # 3 um di larghezza, nel range 1-15 um dichiarato).  Il modello usa la
    # dispersione larga anche per il singolo frammento.
    "friese1998": {"R_lat_m": 1.5e-6, "t_m": 1.0e-6, "P_W": 300e-3,
                   "lambda_m": 1064e-9, "dn": 0.172, "mu_Pas": 0.7e-3,
                   "f_mis_Hz": 357.0},
    # Limite meccanico del manubrio di silice: non e' un esperimento con un
    # numero misurato, ma il confronto con il decennio dei record (Reimann
    # 2018, Ahn 2020: 5-6 GHz su particelle piu' piccole).  Il test chiede
    # che il modello dia 0.5-10 GHz e che la coppia a 10 mW porti al limite
    # in meno di 1 s.
    "limite_ghz": {"R_m": 100e-9, "P_W": 10e-3, "sy_GPa": [1.0, 3.0],
                   "f_lim_GHz": [0.5, 10.0], "t_spin_max_s": 1.0},
}


def _test(nome, esito, dettagli, riferimento, nota=""):
    return {"test": nome, "esito": "PASS" if esito else "FAIL",
            "dettagli": dettagli, "riferimento": riferimento, "nota": nota}


# ---------------------------------------------------------------------------
# A. algebra del ritardo di fase
# ---------------------------------------------------------------------------
def test_a():
    casi = []
    for (S, nome) in (((0.0, 0.0, 1.0), "circolare"), ((1.0, 0.0, 0.0), "lineare")):
        for frac, etichetta, atteso in ((0.5, "lambda/4", None), (1.0, "lambda/2", None)):
            d = frac * math.pi
            v = m._tau_norm(*S, d)
            if nome == "circolare":
                atteso = 1.0 - math.cos(d)          # 1 a lambda/4, 2 a lambda/2
            else:
                atteso = abs(math.sin(d))           # 1 a lambda/4, 0 a lambda/2
            casi.append({"polarizzazione": nome, "lamina": etichetta,
                         "delta": d, "tau_su_P_omega": v, "atteso": atteso,
                         "ok": abs(v - atteso) < 1e-9})
    esito = all(c["ok"] for c in casi)
    return _test(
        "A - algebra del ritardo di fase", esito, casi,
        "Jones + Stokes; struttura verificata a mano",
        "circolare: 1-cos d (2 hbar per fotone a lambda/2).  lineare: sin d "
        "(0 a lambda/2 a 45 gradi: nessun momento angolare trasferito, e "
        "infatti sin(pi) = 0; 1 a lambda/4: la lineare diventa circolare). "
        "La formula sin^2(d/2) della letteratura per il caso lineare NON "
        "riproduce il caso lambda/2 e non e' usata.")


# ---------------------------------------------------------------------------
# B. limite di Beth (1936): lamina lambda/2 in luce circolare
# ---------------------------------------------------------------------------
def test_b():
    d = math.pi
    v = m._tau_norm(0.0, 0.0, 1.0, d)
    # 2 hbar per fotone = 2 P / omega
    esito = abs(v - 2.0) < 1e-9
    return _test(
        "B - limite di Beth: lambda/2 circolare", esito,
        {"tau_su_P_omega": v, "atteso": 2.0, "hbar_per_fotone": v},
        "Beth 1936, Phys. Rev. 50, 115",
        "l'elica si inverte: 2 hbar per fotone.  E' il test che la V6 non "
        "aveva: con sin^2(d/2) si otterrebbe 1 hbar e tutte le coppie "
        "sarebbero dimezzate.")


# ---------------------------------------------------------------------------
# C. termica in acqua contro il conto a mano
# ---------------------------------------------------------------------------
def test_c():
    casi = []
    for R, P_ass, atteso in ((1e-6, 150e-6, 20.0), (3e-6, 150e-6, 6.6)):
        # P_ass e' la potenza ASSORBITA (il conto a mano della V6.1 parte da
        # quella), non la potenza incidente
        th = m.temperatura(P_ass, R, "cu_lucido", "sfera", "acqua")
        casi.append({"R_um": R * 1e6, "P_ass_W": P_ass, "dT": th["dT"],
                     "atteso_K": atteso, "rapporto": th["dT"] / atteso})
    esito = all(abs(c["rapporto"] - 1.0) < 0.15 for c in casi)
    return _test(
        "C - termica in acqua (V6.1)", esito, casi,
        "V6.1: dT = P_abs/(4 pi kappa R); il piano V7 dava 20 K e 6.6 K "
        "per P_abs = 150 uW",
        "verifica di chiusura della V6.1: il modello riproduce il conto a mano "
        "entro il 15%.  La conduttivita' usata e' 0.628 W/m/K (acqua a 310 K), "
        "non 0.6: lo scarto e' il 4% ed e' tutto li'.")


# ---------------------------------------------------------------------------
# D. He 1995: particella assorbente con OAM (controllo di punto)
# ---------------------------------------------------------------------------
def test_d():
    pr = PARAMETRI_LETTERATURA["he1995"]
    R, P, lam = pr["R_m"], pr["P_W"], pr["lambda_m"]
    w0 = pr["w0_m"]
    # Mie LOCALE a 632.8 nm, non la cella globale a 532 nm: il numero della
    # fonte e' a 632.8 nm e usare la lambda sbagliata sarebbe un errore.
    mat = m.MATERIALI[pr["materiale"]]
    q = m.mie_sfera(complex(mat["n"], mat["k"]), 2.0 * math.pi * R / lam)
    I0 = 2.0 * P / (math.pi * w0 ** 2)
    A_proj = math.pi * R ** 2
    P_abs = I0 * A_proj * q["Q_abs"]
    omega = 2.0 * math.pi * m.C.c / lam
    tau = (pr["l"] + pr["s"]) * P_abs / omega
    beta_k = 8.0 * math.pi * pr["mu_Pas"] * R ** 3
    f_model = tau / (2.0 * math.pi * beta_k)
    rapporto = f_model / pr["f_mis_Hz"]
    esito = (1.0 / TOL_VALORE) <= rapporto <= TOL_VALORE
    # Limite superiore: se tutta la potenza fosse assorbita (P_abs = P).
    f_tutta = (pr["l"] + pr["s"]) * P / omega / (2.0 * math.pi * beta_k)
    # Sensibilita' al waist, l'unico parametro non pubblicato.
    w0b = m.W0
    I0b = 2.0 * P / (math.pi * w0b ** 2)
    P_abs_b = I0b * A_proj * q["Q_abs"]
    f_b = (pr["l"] + pr["s"]) * P_abs_b / omega / (2.0 * math.pi * beta_k)
    dettagli = {
        "R_um": R * 1e6, "P_mW": P * 1e3, "lambda_nm": lam * 1e9,
        "l": pr["l"], "s": pr["s"], "w0_um": w0 * 1e6,
        "Q_abs_632nm": q["Q_abs"], "Q_sca_632nm": q["Q_sca"],
        "I0_W_m2": I0, "P_abs_mW": P_abs * 1e3,
        "P_abs_su_P": P_abs / P,
        "tau_Nm": tau, "beta_kerosene": beta_k,
        "f_modello_Hz": f_model, "f_misurata_Hz": pr["f_mis_Hz"],
        "rapporto_modello_misura": rapporto,
        "f_tutta_P_assorbita_Hz": f_tutta,
        "rapporto_tutta_P": f_tutta / pr["f_mis_Hz"],
        "sensibilita_w0": [
            {"w0_um": w0 * 1e6, "f_Hz": f_model, "rapporto": rapporto},
            {"w0_um": w0b * 1e6, "f_Hz": f_b, "rapporto": f_b / pr["f_mis_Hz"]},
        ],
    }
    return _test(
        "D - He 1995: CuO assorbente con OAM", esito, dettagli,
        "He, Friese, Heckenberg, Rubinsztein-Dunlop, PRL 75, 826 (1995)",
        "controllo di punto, non di banda: 632.8 nm, l = 3, 10 mW, CuO "
        "assorbente di ~1 um di raggio in kerosene, 4 Hz misurati.  La Mie "
        "locale a 632.8 nm da' Q_abs = 0.70.  Il modello da' ~9.8x la "
        "misura con w0 = 1.2 um: D FALLISCE al fattore 3.  Il waist non e' "
        "pubblicato; con w0 = 0.5 um il rapporto sale a ~56, quindi il FAIL "
        "e' robusto e non e' un artefatto della scelta.  Non e' una "
        "validazione del modello di ritardo: qui il canale e' "
        "l'assorbimento, e il fallimento e' riportato come tale.",
    )

# ---------------------------------------------------------------------------
# E. Friese 1998: calcite birifrangente (controllo di punto)
# ---------------------------------------------------------------------------
def test_e():
    pr = PARAMETRI_LETTERATURA["friese1998"]
    lam, t = pr["lambda_m"], pr["t_m"]
    delta = 2.0 * math.pi * pr["dn"] * t / lam
    c1, _ = m._media_ritardo(delta)
    om = 2.0 * math.pi * m.C.c / lam
    tau = (pr["P_W"] / om) * c1
    beta = 8.0 * math.pi * pr["mu_Pas"] * pr["R_lat_m"] ** 3
    f_model = tau / (2.0 * math.pi * beta)
    rapporto = f_model / pr["f_mis_Hz"]
    esito = (1.0 / TOL_VALORE) <= rapporto <= TOL_VALORE
    # sensibilita' al raggio laterale, l'unico parametro non pubblicato
    sens = []
    for R_um in (1.0, 1.5, 2.0):
        b = 8.0 * math.pi * pr["mu_Pas"] * (R_um * 1e-6) ** 3
        f = tau / (2.0 * math.pi * b)
        sens.append({"R_lat_um": R_um, "f_Hz": f,
                     "rapporto": f / pr["f_mis_Hz"]})
    dettagli = {
        "t_um": t * 1e6, "P_mW": pr["P_W"] * 1e3,
        "lambda_nm": lam * 1e9, "dn": pr["dn"],
        "delta_rad": delta, "c1": c1, "tau_Nm": tau,
        "beta_acqua": beta, "f_modello_Hz": f_model,
        "f_misurata_Hz": pr["f_mis_Hz"], "rapporto_modello_misura": rapporto,
        "sensibilita_raggio_laterale": sens,
    }
    return _test(
        "E - Friese 1998: calcite birifrangente", esito, dettagli,
        "Friese, Nieminen, Heckenberg, Rubinsztein-Dunlop, Nature 394, 348 (1998)",
        "controllo di punto: calcite da 1 um, 300 mW a 1064 nm, 357 Hz.  "
        "delta = 2 pi dn t/lambda = 1.016, c1 = 0.537 con la media larga; "
        "il raggio laterale non e' pubblicato e si assume 1.5 um.  Il "
        "rapporto 0.68 e' dentro il fattore 3; la sensibilita' mostra che "
        "con R = 1 um il modello darebbe 2.3x, ancora dentro, e con R = 2 um "
        "0.29x, fuori.  Il test misura quindi anche quanto pesa l'assunzione.")


# ---------------------------------------------------------------------------
# F. limite GHz (Reimann 2018, Ahn 2020)
# ---------------------------------------------------------------------------
def test_f():
    pr = PARAMETRI_LETTERATURA["limite_ghz"]
    R, P = pr["R_m"], pr["P_W"]
    casi = []
    sy_orig = m.MATERIALI["silice"]["sy"]
    try:
        for sy_GPa in pr["sy_GPa"]:
            m.MATERIALI["silice"]["sy"] = sy_GPa * 1e9
            c = m.cella("scattering", P, R, "silice", "manubrio", "vuoto", 1e-8)
            casi.append({"sy_GPa": sy_GPa, "I_kgm2": c["I"],
                         "tau_Nm": c["tau"],
                         "omega_lim_rad_s": c["omega_rottura"],
                         "f_lim_GHz": c["omega_rottura"] / 2 / math.pi / 1e9,
                         "t_spin_s": c["t_spin"],
                         "t_spin_limite_s": c["t_spin_limite"]})
    finally:
        m.MATERIALI["silice"]["sy"] = sy_orig
    fmin = min(c["f_lim_GHz"] for c in casi)
    fmax = max(c["f_lim_GHz"] for c in casi)
    decennio_ok = (pr["f_lim_GHz"][0] <= fmin) and (fmax <= pr["f_lim_GHz"][1])
    coppia_ok = all(c["t_spin_s"] < pr["t_spin_max_s"] for c in casi)
    dettagli = {
        "casi": casi, "banda_ammessa_GHz": pr["f_lim_GHz"],
        "decennio_ok": decennio_ok,
        "t_spin_max_s": pr["t_spin_max_s"], "coppia_ok": coppia_ok,
        "nota_coppia": "t_spin e' il tempo per arrivare a meta' del limite "
                       "con la coppia di scattering a 10 mW; t_spin_limite "
                       "e' lo stesso tempo con la coppia massima 2P/omega, "
                       "riportato per mostrare il margine.",
    }
    return _test(
        "F - limite GHz del manubrio di silice", decennio_ok and coppia_ok,
        dettagli,
        "Reimann et al., PRL 121, 033602 (2018); Ahn et al. (2020)",
        "il modello non deve indovinare il record (5-6 GHz, particelle piu' "
        "piccole e piu' resistenti): deve dare il decennio giusto per un "
        "manubrio da 100 nm con snervamento 1-3 GPa, e la coppia a 10 mW "
        "deve bastare a portarlo al limite in meno di 1 s.  Il tempo e' "
        "quello della coppia di scattering, non quello della coppia massima "
        "2P/omega: sono separati di ~5 ordini di grandezza e confonderli "
        "era l'errore della V7.")


def main():
    print("=" * 70)
    print("V7 - CONSISTENCY TESTS (not validation)")
    print("=" * 70)
    test = [test_a(), test_b(), test_c(), test_d(), test_e(), test_f()]
    for t in test:
        print("\n[%s] %s" % (t["esito"], t["test"]))
        print("      rif.: %s" % t["riferimento"])
        d = t["dettagli"]
        if isinstance(d, list):
            for c in d:
                print("      " + "  ".join("%s=%s" % (k, ("%.4g" % v) if isinstance(v, float) else v)
                                           for k, v in c.items()))
        else:
            for k, v in d.items():
                if isinstance(v, list):
                    print("      %-30s %s" % (k, ["%.4g" % x if isinstance(x, float) else x for x in v]))
                else:
                    print("      %-30s %s" % (k, ("%.4g" % v) if isinstance(v, float) else v))
        if t["nota"]:
            print("      nota: %s" % t["nota"])

    n_pass = sum(1 for t in test if t["esito"] == "PASS")
    print("\n%d/%d test superati" % (n_pass, len(test)))
    if n_pass < len(test):
        print("ATTENZIONE: correggere il modello prima di usare la mappa.")

    with open(os.path.join(DATI, "V0_validazione.json"), "w", encoding="utf-8") as fh:
        json.dump({"versione": "V7", "tolleranza_valore": TOL_VALORE,
                   "nota_generale": "V7 - test di coerenza, non validazione. "
                   "D fallisce ai parametri della fonte (rapporto ~9.8): non e' "
                   "una validazione del modello di ritardo.",
                   "parametri_letteratura": PARAMETRI_LETTERATURA,
                   "test": test, "superati": n_pass, "totali": len(test)},
                  fh, indent=1, ensure_ascii=False, default=float)
    print("scritto dati_v7/V0_validazione.json")


if __name__ == "__main__":
    main()
