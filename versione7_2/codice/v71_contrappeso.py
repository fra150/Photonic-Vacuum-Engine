"""
MOTORE FOTONICO DEL VUOTO - V7.1
v71_contrappeso.py - il ciclo di 720 gradi con contrappeso, dentro il sangue.

Appendice dichiarata: il piano V7 teneva fuori 720 gradi, SU(2) e sangue.  Qui
l'idea si porta fino ai numeri.  La regola resta quella della V7: se non regge,
il risultato negativo e' il risultato.

L'IDEA (formulata dall'autore)
------------------------------
"Ogni 360 gradi si attiva il principio del contrappeso, avendo per asse
adiacente il sangue sottostante; solo allora, dopo che il motore e' interno al
sangue, compie un ciclo completo di 720 gradi usando un contrappeso di meta'
del volume del motore."

INTERPRETAZIONE OPERATIVA (dichiarata, da confermare)
-----------------------------------------------------
1. contrappeso = corpo di volume V/2, stessa densita', su un braccio rigido a
   distanza d = 2 R.
2. "si attiva ogni 360 gradi" = latch che innesta e disinnesta una volta per
   giro.  Cio' che ne segue, e che si misura, e' una modulazione periodica di
   omega con periodo 360 gradi: un subarmonico a omega/2.
3. "asse adiacente il sangue" = la reazione del contrappeso va nel fluido, non
   in un cuscinetto: il contrappeso porta il proprio attrito viscoso e il
   mezzo e' sangue intero a 310 K.
4. "ciclo di 720" = lo stato comprende orientazione E configurazione del
   collegamento; dopo 360 il collegamento e' invertito, dopo 720 torna.  E' il
   raddoppio del rivestimento (SO(3) -> SU(2)) del trucco della cintura di
   Dirac: proprieta' TOPOLOGICA del collegamento, non sorgente di momento
   angolare.

COSA LA V7.1 NON DICE
---------------------
- Non produce spinta: F <= 2P/c, come qualunque fascio.
- Il ciclo di 720 non aumenta la coppia: la coppia dipende da P e dal ritardo,
  non dal numero di giri.  Cio' che il 720 cambia e' il TEMPO del ciclo e la
  modulazione del moto, non l'ampiezza.
- Non viola il secondo principio: il latch dissipa e il conto e' qui sotto.

LE TRE GRANDEZZE CHE DECIDONO
-----------------------------
a) profondita' balistica 1/(mu_a+mu_s): la coppia a ritardo richiede luce
   coerente.  Sangue intero: 11 um a 532 nm, 22 um a 800 nm.
b) contrappeso: I e beta crescono di ~4 e ~2.8, quindi il ciclo rallenta.
c) riscaldamento: due contributi, la particella (sfera) e il SANGUE LUNGO IL
   CAMMINO (sorgente lineare, non sferica: e' il termine che decide).

Uso:
    python v71_contrappeso.py
"""

import json
import math
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
V7 = os.path.normpath(os.path.join(BASE, "..", "..", "versione 7", "codice"))
for p in (V7, BASE):
    if p not in sys.path:
        sys.path.insert(0, p)

import v7_mappa as m  # noqa: E402

DATI = os.path.normpath(os.path.join(BASE, "..", "dati"))
os.makedirs(DATI, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. IL SANGUE
# ---------------------------------------------------------------------------
SANGUE = {
    "mu": 4.0e-3,
    "rho": 1060.0,
    "kappa": 0.55,
    "T": 310.0,
    "nota_mu": "sangue intero non newtoniano: 3.5 mPa s ad alto shear, "
               "~10 mPa s a 10 1/s.  A omega ~ 1-1000 1/s il valore 4 mPa s "
               "e' quello giusto entro +-30%.",
}

# attenuazione del sangue intero (m^-1), emoglobina 150 g/L.
OTTICA_SANGUE = {
    532: {"mu_a": 7500.0, "mu_s": 80000.0,
          "nota": "HbO2: eps ~ 13900 M^-1 cm^-1, c = 2.33 mM"},
    650: {"mu_a": 2500.0, "mu_s": 60000.0, "nota": "finestra rossa"},
    800: {"mu_a": 430.0, "mu_s": 45000.0, "nota": "finestra tissutale"},
    1064: {"mu_a": 200.0, "mu_s": 35000.0, "nota": "Hb quasi trasparente"},
}


def profondita_balistica(lam_nm):
    d = OTTICA_SANGUE[lam_nm]
    return 1.0 / (d["mu_a"] + d["mu_s"])


def trasmissione_sangue(lam_nm, z):
    d = OTTICA_SANGUE[lam_nm]
    return math.exp(-(d["mu_a"] + d["mu_s"]) * z)


# ---------------------------------------------------------------------------
# 2. IL CONTRAPPESO
# ---------------------------------------------------------------------------
def contrappeso(R, materiale, forma="sfera", d_rel=2.0, frazione_volume=0.5):
    """Inerzia e attrito del rotore CON il contrappeso innestato.

    I    = I_rot + m_cw (d^2 + 2/5 r_cw^2)
    beta = beta_rot + 8 pi mu r_cw^3 + 6 pi mu r_cw d^2
    """
    mat = m.MATERIALI[materiale]
    g = m.geometria(forma, R)
    V_cw = frazione_volume * g["V"]
    r_cw = (3.0 * V_cw / (4.0 * math.pi)) ** (1.0 / 3.0)
    m_rot = mat["rho"] * g["V"]
    m_cw = mat["rho"] * V_cw
    d = d_rel * R

    I_rot = g["kI"] * m_rot * R ** 2
    I_cw = m_cw * (d ** 2 + 0.4 * r_cw ** 2)
    beta_rot = m.beta_acqua(R, forma) * (SANGUE["mu"] / m.MU_ACQUA)
    beta_cw = 8.0 * math.pi * SANGUE["mu"] * r_cw ** 3 + \
              6.0 * math.pi * SANGUE["mu"] * r_cw * d ** 2

    return {"R": R, "forma": forma, "materiale": materiale,
            "frazione_volume": frazione_volume, "d_rel": d_rel,
            "r_cw": r_cw, "d": d, "m_rot": m_rot, "m_cw": m_cw,
            "I_rot": I_rot, "I_cw": I_cw, "I_tot": I_rot + I_cw,
            "beta_rot": beta_rot, "beta_cw": beta_cw, "beta_tot": beta_rot + beta_cw,
            "costo_inerzia": (I_rot + I_cw) / I_rot,
            "costo_attrito": (beta_rot + beta_cw) / beta_rot}


def contrappeso_assente(R, materiale, forma="sfera"):
    mat = m.MATERIALI[materiale]
    g = m.geometria(forma, R)
    I = g["kI"] * mat["rho"] * g["V"] * R ** 2
    b = m.beta_acqua(R, forma) * (SANGUE["mu"] / m.MU_ACQUA)
    return {"R": R, "I_rot": I, "I_tot": I, "beta_rot": b, "beta_tot": b,
            "m_cw": 0.0, "d": 0.0, "r_cw": 0.0,
            "costo_inerzia": 1.0, "costo_attrito": 1.0}


def omega_latch(R, materiale, forma="sfera", d_rel=2.0, contatto_rel=0.2):
    """Velocita' oltre cui il perno non tiene: m_cw omega^2 d <= sigma_y A."""
    c = contrappeso(R, materiale, forma, d_rel)
    A = (contatto_rel * R) ** 2
    return math.sqrt(m.MATERIALI[materiale]["sy"] * A / (c["m_cw"] * c["d"]))


# ---------------------------------------------------------------------------
# 3. IL CICLO
# ---------------------------------------------------------------------------
def ciclo(R, materiale, lam_nm, z, Pe_target=10.0, forma="sfera", d_rel=2.0,
          innestato=True, w0=None):
    """Ciclo di 720 gradi alla profondita' z nel sangue, tarato su Pe_r target.

    La potenza al bersaglio non e' un parametro libero: si fissa Pe_r e si
    ricava la potenza necessaria.  Poi si guarda quanto costa portarla li'
    (P_inc) e quanto scalda il sangue lungo il cammino.
    """
    w0 = w0 or m.W0
    T_s = trasmissione_sangue(lam_nm, z)
    omega_lam = 2.0 * math.pi * m.C.c / (lam_nm * 1e-9)

    co1 = m.coppia("ritardo", 1.0, R, materiale, forma, lam=lam_nm * 1e-9, w0=w0)
    tau_per_W = co1["tau"] * (m.OMEGA_OTT / omega_lam)
    if not co1["applicabile"] or tau_per_W <= 0.0:
        return {"R": R, "materiale": materiale, "lam_nm": lam_nm, "z": z,
                "valido": False, "motivi": [co1.get("motivo", "nessuna coppia")],
                "tau": 0.0, "omega": 0.0, "Pe_r": 0.0, "T_720_s": float("inf"),
                "E_720_J": float("inf"), "dT": float("inf"), "P_inc": float("inf"),
                "P_target": float("inf"), "T_sangue": T_s, "delta": 0.0,
                "modulazione_omega_su_2": 0.0, "costo_inerzia": 1.0,
                "costo_attrito": 1.0, "omega_rottura": 0.0, "omega_latch": 0.0}

    P_target = Pe_target * m.KT / tau_per_W
    P_inc = P_target / T_s
    tau = tau_per_W * P_target

    c = contrappeso(R, materiale, forma, d_rel) if innestato \
        else contrappeso_assente(R, materiale, forma)
    g = m.geometria(forma, R)
    omega = tau / c["beta_tot"]

    P_abs_part = co1["P_ass"] * P_target
    dT_part = P_abs_part / (4.0 * math.pi * SANGUE["kappa"] * g["Req"])
    # sangue lungo il cammino: sorgente LINEARE, non sferica
    P_path = max(P_inc - P_target, 0.0)
    dT_path = (P_path / (2.0 * math.pi * SANGUE["kappa"] * z)) * \
              math.log(max(z, w0) / w0) if z > 0 else 0.0
    dT = max(dT_part, dT_path)

    w_lim, _ = m.omega_rottura(R, materiale, forma)
    w_latch = omega_latch(R, materiale, forma, d_rel)
    T_360 = 2.0 * math.pi / omega if omega > 0 else float("inf")

    motivi = []
    if dT > m.CRITERI["dT_max_acqua"]:
        motivi.append("dT > 20 K")
    if omega > 0.5 * w_lim:
        motivi.append("omega > 0.5 omega_rottura")
    if omega > w_latch:
        motivi.append("omega > omega_latch")
    if T_s < math.exp(-1.0):
        motivi.append("z oltre la profondita' balistica")
    if P_inc > 1.0:
        motivi.append("P_inc > 1 W (sangue ustionato in superficie)")

    return {
        "R": R, "materiale": materiale, "lam_nm": lam_nm, "z": z,
        "innestato": innestato, "T_sangue": T_s,
        "P_target": P_target, "P_inc": P_inc, "P_ratio": P_inc / P_target,
        "delta": co1["delta"], "tau": tau, "Pe_r": Pe_target,
        "I_tot": c["I_tot"], "beta_tot": c["beta_tot"],
        "omega": omega, "omega_rottura": w_lim, "omega_latch": w_latch,
        "T_360_s": T_360, "T_720_s": 2.0 * T_360,
        "E_720_J": P_inc * 2.0 * T_360,
        "dT": dT, "dT_particella": dT_part, "dT_percorso": dT_path,
        "modulazione_omega_su_2": 1.0 - c["I_rot"] / c["I_tot"],
        "costo_inerzia": c["costo_inerzia"], "costo_attrito": c["costo_attrito"],
        "valido": len(motivi) == 0, "motivi": motivi,
    }


# ---------------------------------------------------------------------------
# 4. GRIGLIA
# ---------------------------------------------------------------------------
R_GRID = tuple(float(x) for x in np.logspace(math.log10(0.15e-6), math.log10(50e-6), 16))
Z_GRID = (2e-6, 5e-6, 10e-6, 15e-6, 20e-6, 30e-6, 50e-6, 100e-6, 300e-6, 1e-3)
LAM = (532, 800)
MAT = ("vaterite", "silice")


def griglia():
    return [ciclo(R, mat, lam, z) for mat in MAT for lam in LAM
            for z in Z_GRID for R in R_GRID]


def riassunto(righe):
    out = {"celle": len(righe),
           "valide": sum(1 for r in righe if r["valido"])}
    out["profondita"] = {
        str(lam): {"profondita_balistica_um": profondita_balistica(lam) * 1e6,
                   "mu_a_cm": OTTICA_SANGUE[lam]["mu_a"] / 100.0,
                   "mu_s_cm": OTTICA_SANGUE[lam]["mu_s"] / 100.0}
        for lam in LAM}
    for lam in LAM:
        cand = [r for r in righe if r["lam_nm"] == lam and r["valido"]]
        if cand:
            b = min(cand, key=lambda r: r["T_720_s"])
            out["migliore_%d_nm" % lam] = {
                "R_um": b["R"] * 1e6, "materiale": b["materiale"],
                "z_um": b["z"] * 1e6, "P_target_uW": b["P_target"] * 1e6,
                "P_inc_mW": b["P_inc"] * 1e3, "T_720_s": b["T_720_s"],
                "E_720_mJ": b["E_720_J"] * 1e3,
                "dT_K": b["dT"], "dT_percorso_K": b["dT_percorso"],
                "tau": b["tau"], "modulazione": b["modulazione_omega_su_2"],
                "costo_inerzia": b["costo_inerzia"]}
        else:
            # nessuna cella valida: riporto la meno peggiore e perche'
            cand = [r for r in righe if r["lam_nm"] == lam and r["omega"] > 0]
            if cand:
                b = min(cand, key=lambda r: max(r["dT"] / 20.0,
                                                math.exp(-1.0) / r["T_sangue"]))
                out["migliore_%d_nm" % lam] = {
                    "nessuna_valida": True, "R_um": b["R"] * 1e6,
                    "z_um": b["z"] * 1e6, "dT_K": b["dT"],
                    "dT_percorso_K": b["dT_percorso"], "T_720_s": b["T_720_s"],
                    "motivi": b["motivi"]}
            else:
                out["migliore_%d_nm" % lam] = None
    return out


def confronto_senza_contrappeso(lam_nm=800, z=10e-6):
    righe = []
    for R in (1e-6, 2e-6, 5e-6, 10e-6):
        a = ciclo(R, "vaterite", lam_nm, z, innestato=False)
        b = ciclo(R, "vaterite", lam_nm, z, innestato=True)
        righe.append({"R_um": R * 1e6,
                      "T_720_senza_s": a["T_720_s"], "T_720_con_s": b["T_720_s"],
                      "omega_senza": a["omega"], "omega_con": b["omega"],
                      "costo_inerzia": b["costo_inerzia"],
                      "costo_attrito": b["costo_attrito"],
                      "rallentamento": a["omega"] / b["omega"] if b["omega"] > 0 else None,
                      "modulazione": b["modulazione_omega_su_2"]})
    return righe


def _scrivi(nome, payload):
    p = os.path.join(DATI, nome)
    with open(p, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1, ensure_ascii=False, default=float)
    print("  scritto %-28s %6.1f kB" % (nome, os.path.getsize(p) / 1024.0))


def main():
    print("=" * 74)
    print("V7.1 - ciclo di 720 gradi con contrappeso, dentro il sangue")
    print("=" * 74)

    print("\n-- profondita' balistica nel sangue intero (serve luce coerente)")
    for lam in LAM:
        d = OTTICA_SANGUE[lam]
        print("   %4d nm: mu_a = %6.0f 1/m, mu_s = %6.0f 1/m -> 1/e = %5.1f um"
              % (lam, d["mu_a"], d["mu_s"], profondita_balistica(lam) * 1e6))

    print("\n-- il prezzo del contrappeso (vaterite, sfera, d = 1.5 R)")
    print("   %-8s %11s %14s %12s %12s" % ("R (um)", "I_tot/I_rot", "beta_tot/beta_rot",
                                           "omega_rott", "omega_latch"))
    for R in (1e-6, 2e-6, 5e-6, 10e-6):
        c = contrappeso(R, "vaterite")
        print("   %-8.1f %11.2f %14.2f %12.2e %12.2e"
              % (R * 1e6, c["costo_inerzia"], c["costo_attrito"],
                 m.omega_rottura(R, "vaterite")[0], omega_latch(R, "vaterite")))

    righe = griglia()
    _scrivi("V71_griglia.json", {"sangue": SANGUE, "ottica": OTTICA_SANGUE,
                                 "R_um": [r * 1e6 for r in R_GRID],
                                 "z_um": [z * 1e6 for z in Z_GRID],
                                 "lambda_nm": list(LAM), "Pe_r_target": 10.0,
                                 "celle": righe})
    rs = riassunto(righe)
    cf = confronto_senza_contrappeso()
    _scrivi("V71_sintesi.json", {
        "descrizione": "V7.1: ciclo di 720 gradi con contrappeso di mezzo volume "
                       "dentro sangue intero a 310 K, tarato su Pe_r = 10",
        "interpretazione": "contrappeso = corpo di volume V/2 a distanza 2 R; "
                           "latch una volta per giro; 720 gradi = raddoppio del "
                           "rivestimento del collegamento (SO(3) -> SU(2))",
        "sangue": SANGUE, "ottica": OTTICA_SANGUE,
        "riassunto": rs, "confronto_contrappeso": cf,
        "cosa_non_dice": [
            "non produce spinta: F <= 2P/c",
            "il ciclo di 720 non aumenta la coppia",
            "il riscaldamento del sangue e' stimato a sorgente lineare; la "
            "perfusione lo riduce, la coagulazione lo peggiora: qui non e' "
            "modellata nessuna delle due",
        ]})

    print("\n-- celle valide: %d / %d" % (rs["valide"], rs["celle"]))
    for lam in LAM:
        b = rs["migliore_%d_nm" % lam]
        if not b:
            print("   %4d nm: nessuna cella" % lam)
        elif b.get("nessuna_valida"):
            print("   %4d nm: NESSUNA valida. La meno peggiore: R = %.1f um, "
                  "z = %.0f um, dT = %.0f K (di cui %.0f K nel sangue), T_720 = %.1f s"
                  % (lam, b["R_um"], b["z_um"], b["dT_K"], b["dT_percorso_K"],
                     b["T_720_s"]))
        else:
            print("   %4d nm: R = %.1f um, z = %.0f um, P_target = %.0f uW, "
                  "P_inc = %.1f mW, T_720 = %.2f s, dT = %.1f K "
                  "(sangue %.1f), modulazione %.0f%%"
                  % (lam, b["R_um"], b["z_um"], b["P_target_uW"], b["P_inc_mW"],
                     b["T_720_s"], b["dT_K"], b["dT_percorso_K"],
                     b["modulazione"] * 100))

    print("\n-- costo dell'idea (800 nm, z = 10 um, Pe_r = 10)")
    print("   %-7s %13s %13s %13s %11s" % ("R (um)", "T_720 senza", "T_720 con",
                                           "rallentamento", "modulazione"))
    for r in cf:
        print("   %-7.1f %13.2f %13.2f %13.2f %10.0f%%"
              % (r["R_um"], r["T_720_senza_s"], r["T_720_con_s"],
                 r["rallentamento"], r["modulazione"] * 100))

    c = ciclo(0.3e-6, "vaterite", 800, 2e-6)
    rif = {"R_um": 0.3, "z_um": 2.0, "lam_nm": 800,
           "P_target_uW": c["P_target"] * 1e6, "P_inc_mW": c["P_inc"] * 1e3,
           "tau": c["tau"], "omega": c["omega"],
           "T_360_s": c["T_360_s"], "T_720_s": c["T_720_s"],
           "E_720_mJ": c["E_720_J"] * 1e3, "dT_K": c["dT"],
           "modulazione": c["modulazione_omega_su_2"],
           "costo_inerzia": c["costo_inerzia"],
           "costo_attrito": c["costo_attrito"],
           "omega_latch": c["omega_latch"], "omega_rottura": c["omega_rottura"]}
    with open(os.path.join(DATI, "V71_sintesi.json"), encoding="utf-8") as fh:
        d = json.load(fh)
    d["cella_riferimento_3d"] = rif
    d["relazione"] = ("T_720 = 4 pi beta / (Pe_r kT): fissato Pe_r, il tempo del "
                      "ciclo dipende solo dall'attrito, cioe' da R^3 e da mu.  "
                      "Non dipende dalla lunghezza d'onda.")
    with open(os.path.join(DATI, "V71_sintesi.json"), "w", encoding="utf-8") as fh:
        json.dump(d, fh, indent=1, ensure_ascii=False, default=float)

    print("\n-- cella di riferimento per il 3D: R = 0.3 um, 800 nm, z = 2 um")
    print("   P_target = %.0f uW, P_inc = %.2f mW, omega = %.2f rad/s, "
          "T_360 = %.2f s, T_720 = %.2f s, dT = %.1f K"
          % (c["P_target"] * 1e6, c["P_inc"] * 1e3, c["omega"],
             c["T_360_s"], c["T_720_s"], c["dT"]))
    print("\nFatto.")


if __name__ == "__main__":
    main()
