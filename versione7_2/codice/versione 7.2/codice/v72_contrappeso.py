"""
MOTORE FOTONICO DEL VUOTO - V7.2
v72_contrappeso.py - il ciclo di 720 gradi con contrappeso, dentro il sangue.

La V7.2 non cambia l'idea della V7.1: la porta a numeri che reggono una
revisione.  Le differenze sono dichiarate qui e finiscono nel JSON.

CORREZIONI RISPETTO ALLA V7.1
-----------------------------
1. 1064 nm entra nella griglia.  La V7.1 aveva i dati ottici ma non li usava.
2. Il costo energetico distingue tre casi:
     - 720 gradi senza contrappeso;
     - 720 gradi con contrappeso innestato per tutto il ciclo (limite
       conservativo della V7.1);
     - 720 gradi con latch a meta' ciclo (il contrappeso entra dopo 360
       gradi, come dice l'idea).  Quest'ultimo e' il caso fisico.
   Il fattore 3.88 della V7.1 e' il caso conservativo, non il caso fisico.
3. La modulazione di omega e' 1 - beta_rot/beta_tot = 74.2%, non
   1 - I_rot/I_tot = 84.2%.  L'84.2% e' la frazione di inerzia: conta nel
   transitorio (nanosecondi), non nella velocita' di regime.
4. Il subarmonico a f_rot/2 e' una proprieta' dello STATO DEL LATCH, non
   dell'orientazione.  L'orientazione e' 360-periodica; solo un osservabile
   dispari sotto 360 gradi (il collegamento A/B) ha periodo 720.
5. La viscosita' del sangue e' un parametro, non una costante: 4 mPa s e'
   il valore ad alto shear; a basso shear e con R < globulo rosso il valore
   effettivo puo' essere 2.5 volte piu' alto.
6. La dispersione di taglia smorza gli zeri della forma di ritardo.  I picchi
   di dT della V7.1 sono zeri di 1-cos(delta), non risonanze Mie del modello.
7. La soglia di 20 K non e' una soglia di sicurezza biologica: a 310 K + 20 K
   il sangue e' a 57 C.  La V7.2 riporta anche la soglia prudente di 5 K
   (42 C) e la distingue dal criterio operativo della V7.1.

Uso:
    python v72_contrappeso.py
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

# Soglie termiche.  20 K e' il criterio operativo della V7.1 (57 C); 5 K e'
# la soglia prudente per esposizione prolungata (42 C).
D_T_OPERATIVA = 20.0
D_T_SICURA = 5.0
SIGMA_Y_BINGHAM = 0.01        # Pa, ordine di grandezza del sangue a basso shear
RBC_RAGGIO = 3.5e-6           # m, globulo rosso (6-8 um di diametro)


def profondita_balistica(lam_nm):
    d = OTTICA_SANGUE[lam_nm]
    return 1.0 / (d["mu_a"] + d["mu_s"])


def trasmissione_sangue(lam_nm, z):
    d = OTTICA_SANGUE[lam_nm]
    return math.exp(-(d["mu_a"] + d["mu_s"]) * z)


def lunghezza_perfusione(kappa=0.55, rho_b=1060.0, cp_b=3600.0,
                         omega_b=8.3e-4):
    """L = sqrt(kappa / (rho_b cp_b omega_b)); 5 mL/100 g/min -> ~1.4 cm."""
    return math.sqrt(kappa / (rho_b * cp_b * omega_b))


# ---------------------------------------------------------------------------
# 2. IL CONTRAPPESO
# ---------------------------------------------------------------------------
def contrappeso(R, materiale, forma="sfera", d_rel=2.0, frazione_volume=0.5,
                mu=None):
    """Inerzia e attrito del rotore CON il contrappeso innestato.

    I    = I_rot + m_cw (d^2 + 2/5 r_cw^2)
    beta = beta_rot + 8 pi mu r_cw^3 + 6 pi mu r_cw d^2
    """
    mu = SANGUE["mu"] if mu is None else mu
    mat = m.MATERIALI[materiale]
    g = m.geometria(forma, R)
    V_cw = frazione_volume * g["V"]
    r_cw = (3.0 * V_cw / (4.0 * math.pi)) ** (1.0 / 3.0)
    m_rot = mat["rho"] * g["V"]
    m_cw = mat["rho"] * V_cw
    d = d_rel * R

    I_rot = g["kI"] * m_rot * R ** 2
    I_cw = m_cw * (d ** 2 + 0.4 * r_cw ** 2)
    beta_rot = m.beta_acqua(R, forma) * (mu / m.MU_ACQUA)
    beta_cw = 8.0 * math.pi * mu * r_cw ** 3 + \
              6.0 * math.pi * mu * r_cw * d ** 2

    I_tot = I_rot + I_cw
    beta_tot = beta_rot + beta_cw
    return {"R": R, "forma": forma, "materiale": materiale, "mu": mu,
            "frazione_volume": frazione_volume, "d_rel": d_rel,
            "r_cw": r_cw, "d": d, "m_rot": m_rot, "m_cw": m_cw,
            "I_rot": I_rot, "I_cw": I_cw, "I_tot": I_tot,
            "beta_rot": beta_rot, "beta_cw": beta_cw, "beta_tot": beta_tot,
            "costo_inerzia": I_tot / I_rot,
            "costo_attrito": beta_tot / beta_rot,
            # modulazione di omega a regime: dipende da beta, non da I
            "modulazione_omega": 1.0 - beta_rot / beta_tot,
            "frazione_inerzia": 1.0 - I_rot / I_tot,
            "tau_spin_s": I_tot / beta_tot,
            "D_r": m.KT / beta_tot,
            "tau_r_s": beta_tot / (2.0 * m.KT)}


def contrappeso_assente(R, materiale, forma="sfera", mu=None):
    mu = SANGUE["mu"] if mu is None else mu
    mat = m.MATERIALI[materiale]
    g = m.geometria(forma, R)
    I = g["kI"] * mat["rho"] * g["V"] * R ** 2
    b = m.beta_acqua(R, forma) * (mu / m.MU_ACQUA)
    return {"R": R, "I_rot": I, "I_tot": I, "beta_rot": b, "beta_tot": b,
            "m_cw": 0.0, "d": 0.0, "r_cw": 0.0, "mu": mu,
            "costo_inerzia": 1.0, "costo_attrito": 1.0,
            "modulazione_omega": 0.0, "frazione_inerzia": 0.0,
            "tau_spin_s": I / b, "D_r": m.KT / b,
            "tau_r_s": b / (2.0 * m.KT)}


def omega_latch(R, materiale, forma="sfera", d_rel=2.0, contatto_rel=0.2):
    """Velocita' oltre cui il perno non tiene: m_cw omega^2 d <= sigma_y A."""
    c = contrappeso(R, materiale, forma, d_rel)
    A = (contatto_rel * R) ** 2
    return math.sqrt(m.MATERIALI[materiale]["sy"] * A / (c["m_cw"] * c["d"]))


def coppia_soglia_bingham(R, sigma_y=SIGMA_Y_BINGHAM):
    """Soglia di snervamento: tau_y ~ (8/3) pi R^3 sigma_y (ordine di grandezza)."""
    return (8.0 / 3.0) * math.pi * R ** 3 * sigma_y


# ---------------------------------------------------------------------------
# 3. IL CICLO
# ---------------------------------------------------------------------------
def media_ritardo(delta, dispersione=0.10, n=21):
    """Media di (1-cos d) e |sin d| su una dispersione uniforme di taglia.

    Con dispersione=0.10 riproduce la V7.1.  Con 0.20 si mostra quanto gli
    zeri della forma di ritardo siano sensibili alla polidispersita'.
    """
    if delta < 0.05:
        return 1.0 - math.cos(delta), math.sin(delta)
    c1 = c2 = 0.0
    for i in range(n):
        t = (i - (n - 1) / 2.0) / ((n - 1) / 2.0) * dispersione
        d = delta * (1.0 + t)
        c1 += 1.0 - math.cos(d)
        c2 += abs(math.sin(d))
    return c1 / n, c2 / n


def coppia_ritardo(R, materiale, forma, lam_nm, w0, dispersione=0.10):
    """Coppia per watt del canale a ritardo, con dispersione di taglia.

    tau/P = eta_int * c_circ / omega_lambda,  c_circ = <1 - cos(delta)>.
    """
    mat = m.MATERIALI[materiale]
    g = m.geometria(forma, R)
    eta = m.eta_int(R, forma, w0)
    dn = m.delta_n(materiale, forma)
    if mat["k"] > 0.05 or dn <= 0.0:
        return {"applicabile": False, "tau_per_W": 0.0, "delta": 0.0,
                "eta_int": eta, "c_circ": 0.0,
                "motivo": ("materiale opaco (k = %.2g)" % mat["k"])
                          if dn <= 0.0 else "nessuna anisotropia"}
    lam = lam_nm * 1e-9
    delta = 2.0 * math.pi * dn * g["percorso"] / lam
    c_circ, c_lin = media_ritardo(delta, dispersione)
    omega_lam = 2.0 * math.pi * m.C.c / lam
    tau_per_W = eta * c_circ / omega_lam
    return {"applicabile": tau_per_W > 0.0, "tau_per_W": tau_per_W,
            "delta": delta, "eta_int": eta, "c_circ": c_circ,
            "c_lin": c_lin, "dispersione": dispersione}


def potenza_assorbita(R, materiale, forma, lam_nm, w0):
    """Potenza assorbita dalla particella per watt incidente sul bersaglio."""
    g = m.geometria(forma, R)
    eta = m.eta_int(R, forma, w0)
    P_int = eta
    I0 = 2.0 / (math.pi * w0 ** 2)
    q = m.q_ottico(materiale, g["Req"], lam_nm * 1e-9)
    P_abs_mie = I0 * g["A_proj"] * max(q["Q_abs"], 0.0)
    P_ass = min(P_abs_mie, P_int * m.alpha_assorbimento(
        materiale, g["percorso"], lam_nm * 1e-9))
    return P_ass, eta, q


def ciclo(R, materiale, lam_nm, z, Pe_target=10.0, forma="sfera", d_rel=2.0,
          innestato=True, w0=None, mu=None, dispersione=0.10,
          latch_fraction=1.0):
    """Ciclo di 720 gradi alla profondita' z nel sangue, tarato su Pe_r.

    latch_fraction = frazione del ciclo di 720 gradi in cui il contrappeso e'
    innestato.  1.0 = limite conservativo V7.1 (sempre innestato);
    0.5 = il latch entra dopo 360 gradi (caso fisico dichiarato).
    """
    w0 = w0 or m.W0
    T_s = trasmissione_sangue(lam_nm, z)
    tr = coppia_ritardo(R, materiale, forma, lam_nm, w0, dispersione)
    if not tr["applicabile"] or tr["tau_per_W"] <= 0.0:
        return {"R": R, "materiale": materiale, "lam_nm": lam_nm, "z": z,
                "valido": False, "valido_sicuro": False,
                "motivi": [tr.get("motivo", "nessuna coppia")],
                "tau": 0.0, "omega": 0.0, "Pe_r": 0.0,
                "T_720_s": float("inf"), "E_720_J": float("inf"),
                "dT": float("inf"), "P_inc": float("inf"),
                "P_target": float("inf"), "T_sangue": T_s, "delta": 0.0,
                "modulazione_omega": 0.0, "frazione_inerzia": 0.0,
                "costo_inerzia": 1.0, "costo_attrito": 1.0,
                "omega_rottura": 0.0, "omega_latch": 0.0}

    P_target = Pe_target * m.KT / tr["tau_per_W"]
    P_inc = P_target / T_s
    tau = tr["tau_per_W"] * P_target

    c_no = contrappeso_assente(R, materiale, forma, mu=mu)
    c_cw = contrappeso(R, materiale, forma, d_rel, mu=mu)
    c = c_cw if innestato else c_no
    omega_no = tau / c_no["beta_tot"] if c_no["beta_tot"] > 0 else float("inf")
    omega_cw = tau / c_cw["beta_tot"] if c_cw["beta_tot"] > 0 else float("inf")
    omega = tau / c["beta_tot"] if c["beta_tot"] > 0 else float("inf")

    T360_no = 2.0 * math.pi / omega_no if omega_no > 0 else float("inf")
    T360_cw = 2.0 * math.pi / omega_cw if omega_cw > 0 else float("inf")
    T720_no = 2.0 * T360_no
    T720_cw = 2.0 * T360_cw
    f = latch_fraction if innestato else 0.0
    T720_latch = 2.0 * ((1.0 - f) * T360_no + f * T360_cw)

    P_abs_part, eta, q = potenza_assorbita(R, materiale, forma, lam_nm, w0)
    P_abs_part *= P_target
    g = m.geometria(forma, R)
    dT_part = P_abs_part / (4.0 * math.pi * SANGUE["kappa"] * g["Req"])
    # sangue lungo il cammino: sorgente LINEARE, non sferica
    P_path = max(P_inc - P_target, 0.0)
    dT_path = (P_path / (2.0 * math.pi * SANGUE["kappa"] * z)) * \
              math.log(max(z, w0) / w0) if z > 0 else 0.0
    dT = max(dT_part, dT_path)

    w_lim, _ = m.omega_rottura(R, materiale, forma)
    w_latch = omega_latch(R, materiale, forma, d_rel)
    tau_y = coppia_soglia_bingham(R)
    tau_ratio = tau / tau_y if tau_y > 0 else float("inf")
    # Re rotazionale: rho omega R^2 / mu (non beta)
    Re_rot = SANGUE["rho"] * omega * R ** 2 / (SANGUE["mu"] if mu is None else mu)

    # Diffusione rotazionale efficace del ciclo a due stati.  La media e'
    # quella di 1/beta pesata nel tempo, non 1/beta medio.
    D_r_eff = (1.0 - f) * (m.KT / c_no["beta_tot"]) + \
              f * (m.KT / c_cw["beta_tot"])
    tau_r_eff = 1.0 / (2.0 * D_r_eff) if D_r_eff > 0 else float("inf")

    motivi = []
    if dT > D_T_OPERATIVA:
        motivi.append("dT > 20 K")
    if omega > 0.5 * w_lim:
        motivi.append("omega > 0.5 omega_rottura")
    if omega > w_latch:
        motivi.append("omega > omega_latch")
    if T_s < math.exp(-1.0):
        motivi.append("z oltre la profondita' balistica")
    if P_inc > 1.0:
        motivi.append("P_inc > 1 W (sangue ustionato in superficie)")

    motivi_sicuri = list(motivi)
    if dT > D_T_SICURA and "dT > 5 K" not in motivi_sicuri:
        motivi_sicuri.append("dT > 5 K (soglia prudente)")

    E_720 = P_inc * T720_latch
    E_720_no = P_inc * T720_no
    E_360_no = P_inc * T360_no

    return {
        "R": R, "materiale": materiale, "lam_nm": lam_nm, "z": z,
        "innestato": innestato, "latch_fraction": f, "T_sangue": T_s,
        "P_target": P_target, "P_inc": P_inc, "P_ratio": P_inc / P_target,
        "delta": tr["delta"], "eta_int": tr["eta_int"], "c_circ": tr["c_circ"],
        "dispersione": dispersione, "tau": tau, "Pe_r": Pe_target,
        "I_tot": c["I_tot"], "beta_tot": c["beta_tot"],
        "beta_no": c_no["beta_tot"], "beta_cw": c_cw["beta_tot"],
        "omega": omega, "omega_no": omega_no, "omega_cw": omega_cw,
        "omega_rottura": w_lim, "omega_latch": w_latch,
        "T_360_s": 2.0 * math.pi / omega if omega > 0 else float("inf"),
        "T_360_no_s": T360_no, "T_360_cw_s": T360_cw,
        "T_720_no_s": T720_no, "T_720_cw_s": T720_cw,
        "T_720_s": T720_latch, "T_720_latch_s": T720_latch,
        "E_720_J": E_720, "E_720_no_J": E_720_no, "E_360_no_J": E_360_no,
        "costo_energia_720": E_720 / E_720_no if E_720_no > 0 else None,
        "costo_energia_360": E_720 / E_360_no if E_360_no > 0 else None,
        "dT": dT, "dT_particella": dT_part, "dT_percorso": dT_path,
        "modulazione_omega": c_cw["modulazione_omega"],
        "frazione_inerzia": c_cw["frazione_inerzia"],
        "costo_inerzia": c["costo_inerzia"],
        "costo_attrito": c["costo_attrito"],
        "tau_spin_s": c["tau_spin_s"], "D_r": c["D_r"],
        "tau_r_s": c["tau_r_s"],
        "D_r_eff": D_r_eff, "tau_r_eff_s": tau_r_eff,
        "tau_bingham": tau_y, "tau_su_bingham": tau_ratio,
        "Re_rot": Re_rot,
        "valido": len(motivi) == 0, "motivi": motivi,
        "valido_sicuro": len(motivi_sicuri) == 0,
        "motivi_sicuri": motivi_sicuri,
    }


# ---------------------------------------------------------------------------
# 4. GRIGLIA
# ---------------------------------------------------------------------------
R_GRID = tuple(float(x) for x in np.logspace(math.log10(0.15e-6),
                                             math.log10(50e-6), 32))
Z_GRID = (2e-6, 5e-6, 10e-6, 15e-6, 20e-6, 30e-6, 50e-6, 100e-6, 300e-6,
          1e-3)
LAM = (532, 800, 1064)
MAT = ("vaterite", "silice")
R_FINE = tuple(float(x) for x in np.logspace(math.log10(0.15e-6),
                                             math.log10(50e-6), 80))
Z_FINE = (2e-6, 5e-6, 10e-6, 15e-6, 20e-6)


def griglia():
    return [ciclo(R, mat, lam, z, latch_fraction=0.5)
            for mat in MAT for lam in LAM for z in Z_GRID for R in R_GRID]


# ---------------------------------------------------------------------------
# 5. SINTESI
# ---------------------------------------------------------------------------
def _migliore(righe, lam, chiave="valido"):
    cand = [r for r in righe if r["lam_nm"] == lam and
            r["materiale"] == "vaterite" and r[chiave]]
    if not cand:
        cand2 = [r for r in righe if r["lam_nm"] == lam and
                 r["materiale"] == "vaterite" and r["omega"] > 0]
        if not cand2:
            return None
        b = min(cand2, key=lambda r: max(r["dT"] / D_T_OPERATIVA,
                                         math.exp(-1.0) / r["T_sangue"]))
        return {"nessuna_valida": True, "R_um": b["R"] * 1e6,
                "z_um": b["z"] * 1e6, "dT_K": b["dT"],
                "dT_percorso_K": b["dT_percorso"],
                "T_720_s": b["T_720_s"], "motivi": b["motivi"]}
    b = min(cand, key=lambda r: r["T_720_s"])
    return {"R_um": b["R"] * 1e6, "materiale": b["materiale"],
            "z_um": b["z"] * 1e6, "P_target_uW": b["P_target"] * 1e6,
            "P_inc_mW": b["P_inc"] * 1e3, "T_720_s": b["T_720_s"],
            "T_720_no_s": b["T_720_no_s"],
            "E_720_mJ": b["E_720_J"] * 1e3,
            "dT_K": b["dT"], "dT_percorso_K": b["dT_percorso"],
            "tau": b["tau"], "delta": b["delta"],
            "modulazione": b["modulazione_omega"],
            "frazione_inerzia": b["frazione_inerzia"],
            "costo_inerzia": b["costo_inerzia"],
            "costo_attrito": b["costo_attrito"],
            "costo_energia_720": b["costo_energia_720"]}


def _finestre(righe, lam, z, chiave="valido", materiale="vaterite"):
    rr = sorted([r for r in righe if r["lam_nm"] == lam and
                 abs(r["z"] - z) < 1e-15 and r["materiale"] == materiale],
                key=lambda r: r["R"])
    out = []
    start = prev = None
    for r in rr:
        if r[chiave]:
            if start is None:
                start = r
            prev = r
        else:
            if start is not None:
                out.append((start, prev))
                start = prev = None
    if start is not None:
        out.append((start, prev))
    res = []
    for a, b in out:
        sub = [x for x in rr if a["R"] <= x["R"] <= b["R"] and x[chiave]]
        res.append({"R_min_um": a["R"] * 1e6, "R_max_um": b["R"] * 1e6,
                    "n_celle": len(sub),
                    "max_dT_K": max(x["dT"] for x in sub),
                    "min_T_720_s": min(x["T_720_s"] for x in sub)})
    return res


def riassunto(righe):
    out = {"celle": len(righe),
           "valide": sum(1 for r in righe if r["valido"]),
           "valide_sicure": sum(1 for r in righe if r["valido_sicuro"])}
    out["profondita"] = {
        str(lam): {"profondita_balistica_um": profondita_balistica(lam) * 1e6,
                   "mu_a_cm": OTTICA_SANGUE[lam]["mu_a"] / 100.0,
                   "mu_s_cm": OTTICA_SANGUE[lam]["mu_s"] / 100.0}
        for lam in LAM}
    for lam in LAM:
        out["migliore_%d_nm" % lam] = _migliore(righe, lam)
        out["migliore_sicuro_%d_nm" % lam] = _migliore(
            righe, lam, "valido_sicuro")
    out["finestre"] = {
        str(lam): {("%.0f_um" % (z * 1e6)): _finestre(righe, lam, z)
                   for z in Z_GRID} for lam in LAM}
    out["finestre_sicure"] = {
        str(lam): {("%.0f_um" % (z * 1e6)): _finestre(righe, lam, z,
                                                       "valido_sicuro")
                   for z in Z_GRID} for lam in LAM}
    return out


def finestre_fini():
    """Finestre termiche su una griglia fitta, per il testo e le figure."""
    out = {}
    for lam in LAM:
        out[str(lam)] = {}
        for z in Z_FINE:
            righe = [ciclo(R, "vaterite", lam, z, latch_fraction=0.5)
                     for R in R_FINE]
            out[str(lam)]["%.0f_um" % (z * 1e6)] = {
                "valide": _finestre(righe, lam, z, "valido"),
                "sicure": _finestre(righe, lam, z, "valido_sicuro"),
            }
    return out


def confronto_lambda(z=2e-6):
    righe = []
    for R in (0.3e-6, 0.5e-6, 1e-6, 2e-6, 5e-6):
        for lam in LAM:
            c = ciclo(R, "vaterite", lam, z, latch_fraction=0.5)
            righe.append({
                "R_um": R * 1e6, "lam_nm": lam, "z_um": z * 1e6,
                "P_target_uW": c["P_target"] * 1e6,
                "P_inc_mW": c["P_inc"] * 1e3,
                "delta": c["delta"], "T_720_s": c["T_720_s"],
                "E_720_mJ": c["E_720_J"] * 1e3,
                "dT_K": c["dT"], "dT_percorso_K": c["dT_percorso"],
                "T_sangue": c["T_sangue"],
                "valido": c["valido"], "valido_sicuro": c["valido_sicuro"],
            })
    return righe


def costo_energetico(z=2e-6, lam=800):
    righe = []
    for R in (0.3e-6, 0.5e-6, 1e-6, 2e-6):
        no = ciclo(R, "vaterite", lam, z, innestato=False)
        full = ciclo(R, "vaterite", lam, z, innestato=True, latch_fraction=1.0)
        half = ciclo(R, "vaterite", lam, z, innestato=True, latch_fraction=0.5)
        righe.append({
            "R_um": R * 1e6, "lam_nm": lam, "z_um": z * 1e6,
            "T_360_no_s": no["T_360_no_s"],
            "T_360_cw_s": full["T_360_cw_s"],
            "T_720_no_s": no["T_720_no_s"],
            "T_720_full_s": full["T_720_s"],
            "T_720_half_s": half["T_720_s"],
            "E_360_no_mJ": no["E_360_no_J"] * 1e3,
            "E_720_no_mJ": no["E_720_no_J"] * 1e3,
            "E_720_full_mJ": full["E_720_J"] * 1e3,
            "E_720_half_mJ": half["E_720_J"] * 1e3,
            "costo_720_full": full["E_720_J"] / no["E_720_no_J"],
            "costo_720_half": half["E_720_J"] / no["E_720_no_J"],
            "costo_360_half": half["E_720_J"] / no["E_360_no_J"],
            "costo_inerzia": full["costo_inerzia"],
            "costo_attrito": full["costo_attrito"],
            "modulazione_omega": full["modulazione_omega"],
            "frazione_inerzia": full["frazione_inerzia"],
            "tau_spin_s": full["tau_spin_s"],
        })
    return righe


def viscosita(lam=800, z=2e-6):
    righe = []
    for mu_mpas in (4.0, 7.0, 10.0):
        for R in (0.3e-6, 1e-6):
            c = ciclo(R, "vaterite", lam, z, mu=mu_mpas * 1e-3,
                      latch_fraction=0.5)
            righe.append({
                "R_um": R * 1e6, "mu_mPa_s": mu_mpas,
                "beta_tot": c["beta_tot"],
                "omega_cw": c["omega_cw"],
                "T_720_s": c["T_720_s"],
                "E_720_mJ": c["E_720_J"] * 1e3,
                "dT_K": c["dT"],
                "tau_r_s": c["tau_r_s"],
                "Re_rot": c["Re_rot"],
                "tau_bingham": c["tau_bingham"],
                "tau_su_bingham": c["tau_su_bingham"],
                "valido": c["valido"],
                "valido_sicuro": c["valido_sicuro"],
            })
    return righe


def rumore(lam=800, z=2e-6):
    righe = []
    for R in (0.3e-6, 0.5e-6, 1e-6):
        c = ciclo(R, "vaterite", lam, z, latch_fraction=0.5)
        T360 = c["T_720_s"] / 2.0
        D_r = c["D_r_eff"]
        A_omega = 0.5 * (c["omega_no"] - c["omega_cw"])
        row = {"R_um": R * 1e6, "lam_nm": lam, "z_um": z * 1e6,
               "D_r_rad2_s": D_r, "tau_r_s": c["tau_r_eff_s"],
               "T_360_s": T360, "T_720_s": c["T_720_s"],
               "sigma_theta_360_rad": math.sqrt(2.0 * D_r * T360),
               "sigma_theta_360_deg": math.degrees(
                   math.sqrt(2.0 * D_r * T360)),
               "A_omega_rad_s": A_omega}
        for binned in (20, 100):
            dt = T360 / binned
            sigma_omega = math.sqrt(2.0 * D_r / dt)
            snr_bin = A_omega / sigma_omega if sigma_omega > 0 \
                else float("inf")
            n10 = (10.0 / snr_bin) ** 2 if snr_bin > 0 else float("inf")
            row["dt_%d_s" % binned] = dt
            row["sigma_omega_%d_rad_s" % binned] = sigma_omega
            row["snr_bin_%d" % binned] = snr_bin
            row["N_cicli_snr10_dt%d" % binned] = n10
            row["tempo_snr10_dt%d_s" % binned] = \
                n10 * c["T_720_s"] if n10 != float("inf") else float("inf")
        righe.append(row)
    return righe


LATCH = {
    "cosa_non_e": "Non e' un cavo nel sangue.  I due cavi della scena Blender "
                  "sono un diagramma del collegamento topologico (nastro di "
                  "Dirac): servono a mostrare che lo stato A e lo stato B si "
                  "scambiano a 360 gradi e tornano a 720.",
    "braccio": "Un linker di lunghezza d = 2 R (DNA origami, PEG rigido, "
               "tether ottico).  A R = 0.3 um il braccio e' lungo 0.6 um: "
               "non e' una bacchetta macroscopica, e' un vincolo molecolare.",
    "contrappeso": "Una bead di volume V/2 (r_cw = 0.794 R, cioe' 0.238 um "
                   "per R = 0.3 um) alla distanza d = 2 R.  Non e' un peso "
                   "meccanico: e' un secondo corpo idrodinamicamente accoppiato.",
    "latch_baseline": "Latch ottico: un secondo spot (o una seconda modalita' "
                      "del fascio) aggancia la bead una volta per giro.  La "
                      "sincronizzazione e' data dal segnale di orientazione; "
                      "il costo energetico del secondo spot non e' incluso in "
                      "P_inc e va aggiunto al bilancio.",
    "meccanismi": [
        {"nome": "spot ottico secondario", "fisica": "trap ottica che si "
         "accende/spegne sulla bead satellite", "tempo": "< 1 us",
         "reversibile": True, "forza": "pN", "criticita": "richiede un "
         "secondo fascio e un feedback di fase"},
        {"nome": "linker DNA con toehold", "fisica": "ibridazione reversibile "
         "A/B", "tempo": "ms-s", "reversibile": True, "forza": "10-100 pN",
         "criticita": "affaticamento e deriva termica"},
        {"nome": "bead magnetica", "fisica": "gradiente di campo che aggancia "
         "la bead", "tempo": "< 1 ms", "reversibile": True,
         "forza": "pN-nN", "criticita": "riscaldamento e ingombro"},
        {"nome": "foto-switch azobenzene", "fisica": "cambio di conformazione "
         "indotto dalla luce", "tempo": "us-ms", "reversibile": True,
         "forza": "pN", "criticita": "fatica fotochimica"},
    ],
}


LIMITI = [
    "Sangue non newtoniano: 4 mPa s e' il valore ad alto shear.  A omega ~ 4 "
    "1/s e con R < globulo rosso la viscosita' effettiva puo' essere 10 mPa s; "
    "beta e T_720 raddoppiano, E_720 raddoppia.  La finestra in R non cambia "
    "perche' P_target non dipende da mu.",
    "Mezzo continuo: a R = 0.3 um il rotore e' piu' piccolo di un globulo "
    "rosso (6-8 um).  Stokes descrive una viscosita' media, non l'urto con "
    "le cellule; il valore efficace va misurato in vitro.",
    "Soglia di snervamento: a basso shear il sangue ha uno sforzo di "
    "snervamento ~0.005-0.01 Pa.  La coppia ottica lo supera di un fattore "
    "~20-40 a R = 0.3 um, ma il transitorio di avvio va verificato.",
    "Termica: la formula di dT_path e' una sorgente lineare in mezzo infinito. "
    "La perfusione (Pennes) ha lunghezza L ~ 1.4 cm per 5 mL/100 g/min: su "
    "2-100 um il sink volumetrico e' piccolo, ma la parete del capillare e il "
    "flusso contano.  La coagulazione locale aumenta assorbimento e viscosita'.",
    "20 K non e' una soglia di sicurezza: 310 K + 20 K = 57 C, coagulazione.  "
    "La soglia prudente per esposizione prolungata e' 5 K (42 C).",
    "I picchi di dT vs R sono zeri di <1-cos(delta)> (forma di ritardo), non "
    "risonanze Mie del modello.  Una vera trattazione Mie aggiunge risonanze "
    "morfologiche; la dispersione di taglia li smorza.",
    "Il latch non e' dimostrato: la V7.2 assume che esista e ne calcola il "
    "costo.  Senza un meccanismo fisico e un tempo di risposta misurato, il "
    "modello resta cinematico.",
]


CORREZIONI_V71 = [
    "La modulazione di omega e' 74.2% (1 - beta_rot/beta_tot), non 84.2%.  "
    "L'84.2% e' 1 - I_rot/I_tot, la frazione di inerzia; conta nel "
    "transitorio (tau_spin ~ 7 ns), non nella velocita' di regime.",
    "Il fattore 3.88 e' il caso conservativo con contrappeso innestato per "
    "tutto il 720.  Con il latch a meta' ciclo (il caso fisico) T_720 e' "
    "2.44 volte un 720 senza contrappeso e 4.88 volte un 360 senza.",
    "Il subarmonico a f_rot/2 e' nello stato del latch (osservabile dispari "
    "sotto 360 gradi), non nell'orientazione.  L'orientazione ha periodo 360.",
    "1064 nm non sposta il bordo a basso R, che e' fissato dalla forma di "
    "ritardo; elimina il buco a 2.34 um e riduce dT a R >= 1 um.",
    "La soglia di 20 K e' un criterio operativo, non una soglia biologica.",
]


def _scrivi(nome, payload):
    p = os.path.join(DATI, nome)
    with open(p, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1, ensure_ascii=False, default=float)
    print("  scritto %-28s %6.1f kB" % (nome, os.path.getsize(p) / 1024.0))


def main():
    print("=" * 74)
    print("V7.2 - ciclo di 720 gradi con contrappeso, dentro il sangue")
    print("=" * 74)

    print("\n-- profondita' balistica nel sangue intero (serve luce coerente)")
    for lam in LAM:
        d = OTTICA_SANGUE[lam]
        print("   %4d nm: mu_a = %6.0f 1/m, mu_s = %6.0f 1/m -> 1/e = %5.1f um"
              % (lam, d["mu_a"], d["mu_s"],
                 profondita_balistica(lam) * 1e6))

    print("\n-- il prezzo del contrappeso (vaterite, sfera, d = 1.5 R)")
    print("   %-8s %11s %14s %12s %12s"
          % ("R (um)", "I_tot/I_rot", "beta_tot/beta_rot",
             "omega_rott", "omega_latch"))
    for R in (1e-6, 2e-6, 5e-6, 10e-6):
        c = contrappeso(R, "vaterite")
        print("   %-8.1f %11.2f %14.2f %12.2e %12.2e"
              % (R * 1e6, c["costo_inerzia"], c["costo_attrito"],
                 m.omega_rottura(R, "vaterite")[0],
                 omega_latch(R, "vaterite")))

    righe = griglia()
    _scrivi("V72_griglia.json", {
        "sangue": SANGUE, "ottica": OTTICA_SANGUE,
        "R_um": [r * 1e6 for r in R_GRID],
        "z_um": [z * 1e6 for z in Z_GRID],
        "lambda_nm": list(LAM), "Pe_r_target": 10.0,
        "soglie": {"dT_operativa_K": D_T_OPERATIVA,
                   "dT_sicura_K": D_T_SICURA},
        "celle": righe})

    rs = riassunto(righe)
    cf = confronto_lambda()
    ce = costo_energetico()
    vis = viscosita()
    rum = rumore()
    ff = finestre_fini()
    _scrivi("V72_sintesi.json", {
        "descrizione": "V7.2: ciclo di 720 gradi con contrappeso di mezzo "
                       "volume dentro sangue intero a 310 K, tarato su "
                       "Pe_r = 10; griglia a 532/800/1064 nm",
        "interpretazione": "contrappeso = corpo di volume V/2 a distanza 2 R; "
                           "latch reversibile una volta per giro; 720 gradi = "
                           "periodo dello STATO del collegamento "
                           "(SO(3) -> SU(2)), non dell'orientazione",
        "sangue": SANGUE, "ottica": OTTICA_SANGUE,
        "soglie": {"dT_operativa_K": D_T_OPERATIVA,
                   "dT_sicura_K": D_T_SICURA,
                   "nota": "20 K = 57 C, coagulazione; 5 K = 42 C, soglia "
                           "prudente per esposizione prolungata"},
        "perfusione": {"lunghezza_perfusione_m": lunghezza_perfusione(),
                       "nota": "5 mL/100 g/min, sangue a 310 K"},
        "riassunto": rs,
        "finestre_fini": ff,
        "confronto_lambda": cf,
        "costo_energetico": ce,
        "viscosita": vis,
        "rumore": rum,
        "latch": LATCH,
        "limiti": LIMITI,
        "correzioni_v71": CORREZIONI_V71,
    })

    # cella di riferimento per la scena 3D e per il testo
    c = ciclo(0.3e-6, "vaterite", 800, 2e-6, latch_fraction=0.5)
    rif = {"R_um": 0.3, "z_um": 2.0, "lam_nm": 800,
           "P_target_uW": c["P_target"] * 1e6,
           "P_inc_mW": c["P_inc"] * 1e3, "tau": c["tau"],
           "omega": c["omega"], "omega_cw": c["omega_cw"],
           "omega_no": c["omega_no"],
           "T_360_s": c["T_360_s"], "T_360_no_s": c["T_360_no_s"],
           "T_360_cw_s": c["T_360_cw_s"],
           "T_720_s": c["T_720_s"], "T_720_no_s": c["T_720_no_s"],
           "E_720_mJ": c["E_720_J"] * 1e3,
           "E_720_no_mJ": c["E_720_no_J"] * 1e3,
           "dT_K": c["dT"], "modulazione": c["modulazione_omega"],
           "frazione_inerzia": c["frazione_inerzia"],
           "costo_inerzia": c["costo_inerzia"],
           "costo_attrito": c["costo_attrito"],
           "costo_energia_720": c["costo_energia_720"],
           "omega_latch": c["omega_latch"],
           "omega_rottura": c["omega_rottura"]}
    with open(os.path.join(DATI, "V72_sintesi.json"), encoding="utf-8") as fh:
        d = json.load(fh)
    d["cella_riferimento_3d"] = rif
    d["relazione"] = ("T_720 = T_360_no + T_360_cw = 2 pi (beta_no + beta_cw) "
                      "/(Pe_r kT) nel caso fisico a meta' ciclo.  Il caso "
                      "conservativo con contrappeso sempre innestato vale "
                      "4 pi beta_cw/(Pe_r kT).")
    with open(os.path.join(DATI, "V72_sintesi.json"), "w", encoding="utf-8") as fh:
        json.dump(d, fh, indent=1, ensure_ascii=False, default=float)

    print("\n-- celle valide: %d / %d (sicure: %d)"
          % (rs["valide"], rs["celle"], rs["valide_sicure"]))
    for lam in LAM:
        b = rs["migliore_%d_nm" % lam]
        if not b:
            print("   %4d nm: nessuna cella" % lam)
        elif b.get("nessuna_valida"):
            print("   %4d nm: NESSUNA valida. La meno peggiore: R = %.2f um, "
                  "z = %.0f um, dT = %.1f K, T_720 = %.1f s"
                  % (lam, b["R_um"], b["z_um"], b["dT_K"], b["T_720_s"]))
        else:
            print("   %4d nm: R = %.2f um, z = %.0f um, P_target = %.0f uW, "
                  "P_inc = %.2f mW, T_720 = %.2f s, dT = %.1f K "
                  "(sangue %.1f), modulazione %.1f%%, E_720 = %.1f mJ"
                  % (lam, b["R_um"], b["z_um"], b["P_target_uW"],
                     b["P_inc_mW"], b["T_720_s"], b["dT_K"],
                     b["dT_percorso_K"], b["modulazione"] * 100,
                     b["E_720_mJ"]))

    print("\n-- costo dell'idea (800 nm, z = 2 um, Pe_r = 10)")
    print("   %-7s %11s %11s %11s %11s %11s"
          % ("R (um)", "T360 no", "T720 no", "T720 full", "T720 half",
             "E half/no"))
    for r in ce:
        print("   %-7.1f %11.2f %11.2f %11.2f %11.2f %11.2f"
              % (r["R_um"], r["T_360_no_s"], r["T_720_no_s"],
                 r["T_720_full_s"], r["T_720_half_s"],
                 r["costo_720_half"]))

    print("\n-- confronto lunghezze d'onda (R = 0.3 um, z = 2 um, latch 1/2)")
    print("   %-6s %10s %11s %9s %10s %9s"
          % ("lam", "P_target", "P_inc", "T_720", "E_720", "dT"))
    for r in cf:
        if abs(r["R_um"] - 0.3) < 1e-9:
            print("   %-6d %8.0f uW %8.2f mW %7.2f s %7.1f mJ %7.1f K"
                  % (r["lam_nm"], r["P_target_uW"], r["P_inc_mW"],
                     r["T_720_s"], r["E_720_mJ"], r["dT_K"]))

    print("\n-- rumore rotazionale (800 nm, z = 2 um, 100 bin/ciclo)")
    print("   %-6s %10s %10s %12s %12s %12s"
          % ("R (um)", "D_r", "tau_r", "sigma_theta", "SNR/bin", "N(SNR10)"))
    for r in rum:
        print("   %-6.1f %10.3f %10.2f %9.0f deg %12.2f %12.1f"
              % (r["R_um"], r["D_r_rad2_s"], r["tau_r_s"],
                 r["sigma_theta_360_deg"], r["snr_bin_100"],
                 r["N_cicli_snr10_dt100"]))

    print("\nFatto.")


if __name__ == "__main__":
    main()
