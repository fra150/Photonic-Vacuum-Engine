"""
MOTORE FOTONICO DEL VUOTO - V7.3
v73_taglio.py - il cancello sul flusso: soglia di shear, non finestra.

Questo script non aggiunge un modello: decide se il rotore puo' vivere in un
vaso con flusso.  La risposta e' no, e la decisione non dipende dalle costanti
otiche che stiamo ancora cercando di calibrare (lambda, Delta_n, Q, eta_int):
dipende solo da tau_ott, beta_cw e mu.

QUATTRO USCITE
--------------
1. gamma_coppia = (tau_ott - M_y) / beta_cw.
   Corpo vincolato (rotore + bead su linker): il taglio impone un rate di
   deformazione, coefficiente beta_cw.  M_y = 4 pi R^3 sigma_y e' il momento
   di snervamento del Bingham, trattato come scostamento additivo.
   Plasma, sigma_y = 0: 2.700 s^-1.
   Sfera libera (verifica di convenzione): 2 tau_ott / beta_no = 21.0 s^-1.

2. gamma_forza = F_trap / (6 pi mu R_eff z_c), con F_trap = Q n P / c.
   F_trap non dipende da R; gamma_forza va come 1/R (R_eff proporzionale a R).
   Q = 0.1-0.3: 5.35-16.05 s^-1.

3. T_J = 2 pi (r_e + 1/r_e) / gamma, r_e = 1.897.
   Solo verifica di coerenza: T_J = T720 a gamma = 5.20 s^-1.  L'orbita di
   Jeffery non sta nel piano flusso-gradiente per un asse parallelo alla
   parete: non e' un risultato, e' un controllo.

4. U0 = F_trap w0 / 4; tau_esc = tau_D exp(U0/kT), tau_D = xi / k_rad.
   Q = 0.1: U0 = 8.6 kT, tau_esc ~ 2e2 s, al limite di un esperimento di 61 s.
   Q = 0.2-0.3: 17-26 kT, fuga non osservabile.

SOGLIA
------
gamma_crit = min(output 1, output 2) = 2.700 s^-1.
Il minimo va dichiarato: non e' una finestra.  Il valore del minimo diminuisce
con R (gamma_coppia ~ 1/(mu R^3)), mentre il ramo che lo realizza cambia con R;
non c'e' alcun incrocio a 0.422 um -- quello era un artefatto Rayleigh e resta
ritrattato.  Che i tre numeri cadano in 2.7-5.4 s^-1 (fattore 2.0) e' una
coincidenza; il valore robusto e' il minimo, 2.7 s^-1, indipendente da lambda,
Delta_n_eff, Q e eta_int.  Dipende da mu e R (gamma_coppia ~ 1/(mu R^3)):
sono gli input da dichiarare, non da nascondere.

CONVENZIONI
-----------
- corpo vincolato: guida il rate di deformazione, coefficiente beta_cw;
- sfera libera: guida la vorticity gamma/2, coefficiente beta_no/2;
- Bingham: 4 pi R^3 sigma_y e' il criterio di parete conservativo (area x
  braccio).  La soluzione esatta della sfera in Bingham da' (8/3) pi R^3
  sigma_y, piu' basso di 1.5x; qui si usa il valore conservativo e si
  dichiara il fattore.
- CFL = limite del dominio, non ambiente: esiste solo con flusso.  Plasma +
  stasi si realizza solo in vitro.

Uso:
    python codice/v73_taglio.py
"""

import json
import math
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(BASE, ".."))
if BASE not in sys.path:
    sys.path.insert(0, BASE)

import v73_contrappeso as c  # noqa: E402

DATI = os.path.join(ROOT, "dati")
os.makedirs(DATI, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. INPUT DICHIARATI
# ---------------------------------------------------------------------------
T0 = c.T0                       # 310 K
KT = c.KT                       # J
C_LUCE = c.C_LUCE
N_MED = 1.33                    # indice del mezzo (plasma/tampone)
PE_R = c.PE_R                   # 10
TAU_OTT = PE_R * KT             # 4.280e-20 N m

R_REF = 0.5e-6                  # m
LAM_REF = 800                   # nm
Z_REF = 2.0e-6                  # m
Z_C = 2.5e-6                    # m, centro della bead sotto la parete
W0 = 0.5e-6                     # m, waist di riferimento
MU_PLASMA = 1.3e-3              # Pa s, newtoniano
MU_RANGE = (1.3e-3, 1.5e-3)     # plasma/tampone
MU_BLOOD = (7.0e-3, 8.0e-3, 9.0e-3)
SIGMA_Y_BLOOD = (5.0e-3, 10.0e-3, 20.0e-3)
Q_BAND = (0.10, 0.20, 0.30)
SIGMA_Y_RANGE = (5.0e-3, 10.0e-3)

# Drag del corpo a due sfere, calibrato in V7.2: 17.6 pN a 800 um/s.
XI_REF = 17.6e-12 / 800e-6      # 2.2e-8 N s/m
R_EFF_FACTOR = XI_REF / (6.0 * math.pi * MU_PLASMA * R_REF)  # 1.7956

FISIO = {
    "venula 30 um": (100.0, 300.0),      # s^-1, ordini dichiarati
    "capillare 6 um": (100.0, 500.0),
    "arteriola": (1500.0, 2000.0),
}
FISIO_PUNTI = {"venula": 267.0, "capillare": 667.0, "arteriola": 2000.0}

# Riferimento V7.2 corretta (R=0.5 um, 800 nm, z=2 um), ricalcolato qui:
_rif = c.ciclo(R_REF, LAM_REF, Z_REF, mu=MU_PLASMA)
P_REF = _rif["P_inc"]           # 0.662775 mW
T720_REF = _rif["T720_half"]    # 2.926 s
E720_REF = _rif["E720_half"]    # 1.940 mJ

# ---------------------------------------------------------------------------
# 2. LE QUATTRO USCITE
# ---------------------------------------------------------------------------
def beta_no(R, mu):
    return 8.0 * math.pi * mu * R ** 3


def beta_cw(R, mu):
    return c.BETA_RATIO * beta_no(R, mu)


def correzione_idrodinamica():
    """Interazione a due sfere: (r_cw/d)^2 applicato al termine della bead.

    Il rapporto ideale 3.881 = 1 + r_cw^3 + 0.75 r_cw d^2 (r_cw = 0.5^(1/3),
    d = 2R) trascura l'interazione idrodinamica fra rotore e bead.  La
    correzione (r_cw/d)^2 = 0.1575 si applica al solo termine traslazionale
    della bead (limite inferiore), a tutta la bead (valore centrale) o
    all'intero attrito cw (limite superiore).  Il calcolo di rumore usa
    3.881, il valore conservativo.
    """
    r_cw = 0.5 ** (1.0 / 3.0)
    d_rel = 2.0
    eps = (r_cw / d_rel) ** 2
    term_rot = 1.0
    term_bead_rot = r_cw ** 3
    term_bead_trans = 0.75 * r_cw * d_rel ** 2
    ratio = c.BETA_RATIO
    return {
        "r_cw": r_cw, "d_rel": d_rel, "termine": eps,
        "term_rot": term_rot, "term_bead_rot": term_bead_rot,
        "term_bead_trans": term_bead_trans,
        "ratio_ideale": ratio,
        "ratio_bead": term_rot + term_bead_rot
                      + term_bead_trans * (1.0 + eps),
        "ratio_bead_totale": term_rot
                             + (term_bead_rot + term_bead_trans) * (1.0 + eps),
        "ratio_limite_superiore": ratio * (1.0 + eps),
        "nota": "correzione non validata; il rumore usa 3.881",
    }


def momento_snervamento(R, sigma_y):
    """M_y = 4 pi R^3 sigma_y (criterio di parete conservativo)."""
    return 4.0 * math.pi * R ** 3 * sigma_y


def gamma_coppia(R=R_REF, mu=MU_PLASMA, sigma_y=0.0, tau=TAU_OTT):
    """Soglia di coppia per il corpo vincolato, e la variante sfera libera."""
    My = momento_snervamento(R, sigma_y)
    drive = tau - My
    b_cw = beta_cw(R, mu)
    b_no = beta_no(R, mu)
    if drive <= 0.0:
        return {"gamma_vincolato": 0.0, "gamma_sfera_libera": 0.0,
                "M_y": My, "tau": tau, "drive": drive,
                "beta_cw": b_cw, "beta_no": b_no, "bloccato": True}
    return {"gamma_vincolato": drive / b_cw,
            "gamma_sfera_libera": 2.0 * tau / b_no,
            "M_y": My, "tau": tau, "drive": drive,
            "beta_cw": b_cw, "beta_no": b_no, "bloccato": False}


def forza_trappola(P_inc, Q, n=N_MED):
    """F = Q n P / c.  Indipendente da R (a Q fissato)."""
    return Q * n * P_inc / C_LUCE


def gamma_forza(R=R_REF, P_inc=P_REF, Q=0.10, mu=MU_PLASMA, z_c=Z_C):
    """gamma = F_trap / (6 pi mu R_eff z_c), con R_eff = 1.7956 R."""
    F = forza_trappola(P_inc, Q)
    R_eff = R_EFF_FACTOR * R
    xi = 6.0 * math.pi * mu * R_eff
    v = F / xi
    return {"F_N": F, "F_pN": F * 1e12, "v_m_s": v, "v_um_s": v * 1e6,
            "gamma": v / z_c, "R_eff": R_eff, "xi": xi, "Q": Q}


def jeffery(R=None, T720=T720_REF):
    """T_J = 2 pi (r_e + 1/r_e)/gamma; r_e = (3 + 2^-1/3)/2 = 1.897."""
    r_e = (3.0 + 2.0 ** (-1.0 / 3.0)) / 2.0
    coeff = 2.0 * math.pi * (r_e + 1.0 / r_e)
    return {"r_e": r_e, "coeff": coeff, "T_J_s": coeff,
            "gamma_a_T720": coeff / T720}


def profondita_buca(F, w0=W0):
    """U0 = F w0 / 4."""
    return F * w0 / 4.0


def tempo_fuga(F, mu=MU_PLASMA, R=R_REF, w0=W0, T=T0):
    """tau_esc = tau_D exp(U0/kT), tau_D = xi/k_rad, k_rad = F/w0."""
    R_eff = R_EFF_FACTOR * R
    xi = 6.0 * math.pi * mu * R_eff
    k_rad = F / w0
    tau_D = xi / k_rad
    U0 = profondita_buca(F, w0)
    x = U0 / (KT * T / T0)
    return {"U0_J": U0, "U0_kT": x, "tau_D_s": tau_D,
            "tau_esc_s": tau_D * math.exp(min(x, 700.0)),
            "tau_esc_log10": math.log10(max(tau_D, 1e-300)) + x / math.log(10.0)}


def sedimentazione(R=R_REF, F_trap=None, w0=W0,
                   rho_p=2710.0, rho_f=1025.0, g=9.81):
    """Spostamento assiale da galleggiamento, k_ax = 0.13 k_rad."""
    if F_trap is None:
        F_trap = forza_trappola(P_REF, 0.10)
    V = 4.0 / 3.0 * math.pi * R ** 3
    Fb = (rho_p - rho_f) * V * g
    k_rad = F_trap / w0
    k_ax = 0.13 * k_rad
    dz = Fb / k_ax
    return {"F_buoy_N": Fb, "F_buoy_fN": Fb * 1e15,
            "k_rad": k_rad, "k_ax": k_ax,
            "dz_m": dz, "dz_nm": dz * 1e9, "dz_su_w0": dz / w0}


def raggio_parete(sigma_y, tau=TAU_OTT):
    """R tale che 4 pi R^3 sigma_y = tau."""
    return (tau / (4.0 * math.pi * sigma_y)) ** (1.0 / 3.0)


def confronto_fisiologico(gamma_lim):
    return {k: {"gamma": v, "rapporto": v / gamma_lim}
            for k, v in FISIO_PUNTI.items()}


# ---------------------------------------------------------------------------
# 3. I DUE RAMI
# ---------------------------------------------------------------------------
def ramo_A():
    """Plasma/tampone in vitro: mu newtoniano, nessun snervamento, nessun CFL."""
    return {
        "nome": "A - plasma/tampone in vitro (stasi)",
        "mu_mPa_s": [m * 1e3 for m in MU_RANGE],
        "reologia": "newtoniana",
        "sigma_y_mPa": 0.0,
        "cfl": "assente (stasi)",
        "cellule": "assenti; nessun rumore di collisione",
        "parete": "vetro/PDMS (pozzo termico)",
        "T720_s": T720_REF,
        "E720_mJ": E720_REF * 1e3,
        "tre_cicli_s": 3.0 * T720_REF,
        "gamma_lim_s-1": gamma_coppia(mu=MU_PLASMA)["gamma_vincolato"],
        "nota": "protocollo appaiato 3+3; i 3 cicli servono ai sistematici, "
                "non alla statistica (la statistica non migliora con N); "
                "termica non vincolante",
    }


def ramo_B():
    """Sangue intero in stasi: mu 7-9 mPa s, Bingham 5-20 mPa."""
    out = []
    for mu, sy in zip(MU_BLOOD, SIGMA_Y_BLOOD):
        My = momento_snervamento(R_REF, sy)
        fattore = 1.0 + My / TAU_OTT
        T720 = T720_REF * mu / MU_PLASMA
        P = P_REF * fattore
        E = P * T720
        out.append({
            "mu_mPa_s": mu * 1e3,
            "sigma_y_mPa": sy * 1e3,
            "M_y_Nm": My,
            "fattore_potenza": fattore,
            "T720_s": T720,
            "P_target_mW": P * 1e3,
            "E720_mJ": E * 1e3,
            "tre_cicli_s": 3.0 * T720,
            "gamma_coppia_s-1": gamma_coppia(mu=mu, sigma_y=sy)["gamma_vincolato"],
            "E720_su_V72": E / E720_REF,
        })
    return {
        "nome": "B - sangue intero in vivo (stasi)",
        "mu_mPa_s": [m * 1e3 for m in MU_BLOOD],
        "reologia": "Bingham",
        "sigma_y_mPa": [s * 1e3 for s in SIGMA_Y_BLOOD],
        "cfl": "assente (stasi)",
        "cellule": "ematicrito pieno, rouleaux, aggregazione",
        "parete": "endotelio (pozzo termico)",
        "casi": out,
        "nota": "potenza aumentata per compensare lo snervamento; il tempo "
                "scala solo con mu, l'energia porta il fattore (1+M_y/tau); "
                "la distinzione A/B resta 3.881 anche con lo scostamento",
    }


# ---------------------------------------------------------------------------
# 4. SINTESI
# ---------------------------------------------------------------------------
def sintesi():
    g = gamma_coppia()
    gf = [gamma_forza(Q=q) for q in Q_BAND]
    j = jeffery()
    gamma_lim = min(g["gamma_vincolato"], min(x["gamma"] for x in gf))
    fuga = [{"Q": x["Q"], "F_pN": x["F_pN"],
             "U0_kT": profondita_buca(x["F_N"]) / KT,
             "tau_esc_s": tempo_fuga(x["F_N"])["tau_esc_s"]}
            for x in gf]
    sed = [{"Q": x["Q"], **sedimentazione(F_trap=x["F_N"])} for x in gf]
    parete = [{"sigma_y_mPa": s * 1e3,
               "R_um": raggio_parete(s) * 1e6}
              for s in (5e-3, 6.5e-3, 10e-3, 20e-3)]
    return {
        "descrizione": "V7.3: cancello sul flusso, quattro uscite, due rami",
        "input": {
            "R_um": R_REF * 1e6, "lambda_nm": LAM_REF, "z_um": Z_REF * 1e6,
            "z_c_um": Z_C * 1e6, "w0_um": W0 * 1e6,
            "P_inc_mW": P_REF * 1e3, "Pe_r": PE_R,
            "tau_ott_Nm": TAU_OTT, "mu_plasma_mPa_s": MU_PLASMA * 1e3,
            "mu_blood_mPa_s": [m * 1e3 for m in MU_BLOOD],
            "sigma_y_blood_mPa": [s * 1e3 for s in SIGMA_Y_BLOOD],
            "Q_band": list(Q_BAND), "xi_ref_Ns_m": XI_REF,
            "R_eff_factor": R_EFF_FACTOR,
            "T720_ref_s": T720_REF, "E720_ref_mJ": E720_REF * 1e3,
        },
        "output_1_coppia": {
            "gamma_vincolato_s-1": g["gamma_vincolato"],
            "gamma_sfera_libera_s-1": g["gamma_sfera_libera"],
            "M_y_Nm": g["M_y"], "drive_Nm": g["drive"],
            "beta_cw": g["beta_cw"], "beta_no": g["beta_no"],
            "convenzione": "vincolato: tau/beta_cw; libera: 2 tau/beta_no",
        },
        "output_2_forza": [
            {"Q": x["Q"], "F_pN": x["F_pN"], "v_um_s": x["v_um_s"],
             "gamma_s-1": x["gamma"], "R_eff_um": x["R_eff"] * 1e6,
             "xi_Ns_m": x["xi"]} for x in gf],
        "output_3_jeffery": {
            "r_e": j["r_e"], "coeff_s": j["coeff"],
            "T_J_a_T720_s": j["coeff"],
            "gamma_a_T720_s-1": j["gamma_a_T720"],
            "ruolo": "verifica di coerenza; orbita non nel piano flusso-gradiente",
        },
        "output_4_fuga": fuga,
        "gamma_limite_s-1": gamma_lim,
        "gamma_limite_regola": "min(output 1, output 2); il minimo va dichiarato",
        "coincidenza": {
            "tre_stime_s-1": [g["gamma_vincolato"], gf[0]["gamma"],
                              j["gamma_a_T720"]],
            "intervallo": [g["gamma_vincolato"], gf[0]["gamma"]],
            "fattore": gf[0]["gamma"] / g["gamma_vincolato"],
            "nota": "2.7-5.4 s^-1 e' una coincidenza, non una finestra; "
                    "Jeffery 5.20 cade fra 2.70 e 5.40",
        },
        "confronto_fisiologico": confronto_fisiologico(gamma_lim),
        "fisiologia": {
            k: {"gamma_range_s-1": list(FISIO[k]),
                "gamma_punto_s-1": FISIO_PUNTI[{"venula 30 um": "venula",
                                                "capillare 6 um": "capillare",
                                                "arteriola": "arteriola"}[k]]}
            for k in FISIO
        },
        "sensibilita_mu": {
            "plasma_1.3": gamma_coppia(mu=1.3e-3)["gamma_vincolato"],
            "plasma_1.5": gamma_coppia(mu=1.5e-3)["gamma_vincolato"],
            "sangue_7_sy5": gamma_coppia(mu=7e-3, sigma_y=5e-3)["gamma_vincolato"],
            "sangue_8_sy10": gamma_coppia(mu=8e-3, sigma_y=10e-3)["gamma_vincolato"],
            "sangue_9_sy20": gamma_coppia(mu=9e-3, sigma_y=20e-3)["gamma_vincolato"],
        },
        "sedimentazione": sed,
        "criterio_parete": parete,
        "correzione_idrodinamica": correzione_idrodinamica(),
        "ramo_A": ramo_A(),
        "ramo_B": ramo_B(),
        "convenzioni": {
            "bingham": "M_y = 4 pi R^3 sigma_y (parete conservativa); la "
                       "soluzione esatta della sfera in Bingham da' "
                       "(8/3) pi R^3 sigma_y, fattore 1.5 piu' basso",
            "cfl": "limite del dominio, non ambiente; esiste solo con flusso",
            "plasma_stasi": "plasma + stasi si realizza solo in vitro",
            "parete_termica": "il segno della sorgente immagine va dichiarato: "
                              "vetro/endotelio = pozzo; parete adiabatica = "
                              "specchio",
        },
        "verdetto": "flusso chiuso: gamma_limite 2.70 s^-1 sotto ogni shear "
                    "fisiologico (venula x99, capillare x247); il dispositivo "
                    "vive in stasi ex vivo o in stasi in vivo",
    }


def stampa(s):
    print("=" * 78)
    print("V7.3 - v73_taglio.py: soglia di shear, non finestra")
    print("=" * 78)
    i = s["input"]
    print("\n-- input")
    print("   R=%.2f um  lambda=%d nm  z=%.1f um  z_c=%.1f um  w0=%.2f um"
          % (i["R_um"], i["lambda_nm"], i["z_um"], i["z_c_um"], i["w0_um"]))
    print("   P=%.4f mW  Pe_r=%.1f  tau=%.4e N m  mu=%.2f mPa s"
          % (i["P_inc_mW"], i["Pe_r"], i["tau_ott_Nm"], i["mu_plasma_mPa_s"]))
    print("   xi=%.2e N s/m  R_eff/R=%.4f  T720=%.3f s  E720=%.3f mJ"
          % (i["xi_ref_Ns_m"], i["R_eff_factor"], i["T720_ref_s"],
             i["E720_ref_mJ"]))

    o1 = s["output_1_coppia"]
    print("\n-- output 1: coppia")
    print("   gamma_vincolato   = %.4f s^-1  (tau/beta_cw)" % o1["gamma_vincolato_s-1"])
    print("   gamma_sfera_libera= %.4f s^-1  (2 tau/beta_no, convenzione)"
          % o1["gamma_sfera_libera_s-1"])
    print("   M_y=%.3e N m  drive=%.3e N m  beta_cw=%.3e  beta_no=%.3e"
          % (o1["M_y_Nm"], o1["drive_Nm"], o1["beta_cw"], o1["beta_no"]))

    print("\n-- output 2: forza")
    for x in s["output_2_forza"]:
        print("   Q=%.2f  F=%.3f pN  v=%.2f um/s  gamma=%.2f s^-1  R_eff=%.3f um"
              % (x["Q"], x["F_pN"], x["v_um_s"], x["gamma_s-1"],
                 x["R_eff_um"]))

    o3 = s["output_3_jeffery"]
    print("\n-- output 3: Jeffery (solo verifica)")
    print("   r_e=%.4f  coeff=%.3f s  T_J=T720 a gamma=%.2f s^-1"
          % (o3["r_e"], o3["coeff_s"], o3["gamma_a_T720_s-1"]))

    print("\n-- output 4: profondita' e fuga")
    for x in s["output_4_fuga"]:
        print("   Q=%.2f  U0=%.2f kT  tau_esc=%.3e s"
              % (x["Q"], x["U0_kT"], x["tau_esc_s"]))

    print("\n-- soglia")
    print("   gamma_crit = min(%.3f, %.3f) = %.3f s^-1"
          % (s["coincidenza"]["tre_stime_s-1"][0],
             s["coincidenza"]["tre_stime_s-1"][1], s["gamma_limite_s-1"]))
    print("   le tre stime: %.2f / %.2f / %.2f s^-1 (fattore %.2f)"
          % (*s["coincidenza"]["tre_stime_s-1"],
             s["coincidenza"]["fattore"]))
    for k, v in s["confronto_fisiologico"].items():
        print("   %-14s gamma=%6.1f s^-1  rapporto x%.1f"
              % (k, v["gamma"], v["rapporto"]))

    print("\n-- sensibilita' a mu e sigma_y")
    for k, v in s["sensibilita_mu"].items():
        print("   %-18s gamma=%.4f s^-1" % (k, v))

    print("\n-- sedimentazione (Q=0.1-0.3)")
    for x in s["sedimentazione"]:
        print("   Q=%.2f  F_buoy=%.2f fN  dz=%.0f nm (%.1f%% w0)"
              % (x["Q"], x["F_buoy_fN"], x["dz_nm"], 100.0 * x["dz_su_w0"]))

    print("\n-- criterio di parete Bingham: 4 pi R^3 sigma_y = tau")
    for x in s["criterio_parete"]:
        print("   sigma_y=%.1f mPa -> R=%.2f um"
              % (x["sigma_y_mPa"], x["R_um"]))

    h = s["correzione_idrodinamica"]
    print("\n-- correzione idrodinamica a due sfere")
    print("   (r_cw/d)^2 = %.4f" % h["termine"])
    print("   ratio ideale %.4f -> bead %.4f (%+.1f%%), bead totale %.4f "
          "(%+.1f%%), limite %.4f (%+.1f%%)"
          % (h["ratio_ideale"], h["ratio_bead"],
             100.0 * (h["ratio_bead"] / h["ratio_ideale"] - 1.0),
             h["ratio_bead_totale"],
             100.0 * (h["ratio_bead_totale"] / h["ratio_ideale"] - 1.0),
             h["ratio_limite_superiore"],
             100.0 * (h["ratio_limite_superiore"] / h["ratio_ideale"] - 1.0)))
    print("   %s" % h["nota"])

    print("\n-- ramo A: plasma/tampone in vitro")
    a = s["ramo_A"]
    print("   mu=%s mPa s  T720=%.2f s  E720=%.2f mJ  3 cicli=%.1f s"
          % (a["mu_mPa_s"], a["T720_s"], a["E720_mJ"], a["tre_cicli_s"]))
    print("   %s" % a["nota"])

    print("\n-- ramo B: sangue intero in stasi")
    for x in s["ramo_B"]["casi"]:
        print("   mu=%.0f mPa s  sigma_y=%.0f mPa  T720=%.1f s  E720=%.1f mJ"
              "  3 cicli=%.0f s  E/E72=%.1f"
              % (x["mu_mPa_s"], x["sigma_y_mPa"], x["T720_s"], x["E720_mJ"],
                 x["tre_cicli_s"], x["E720_su_V72"]))

    print("\n-- verdetto")
    print("   %s" % s["verdetto"])
    print("   %s" % s["convenzioni"]["plasma_stasi"])


def main():
    s = sintesi()
    p = os.path.join(DATI, "V73_taglio.json")
    with open(p, "w", encoding="utf-8") as fh:
        json.dump(s, fh, indent=1, ensure_ascii=False, default=float)
    stampa(s)
    print("\n  scritto %s (%.1f kB)" % (p, os.path.getsize(p) / 1024.0))
    return s


if __name__ == "__main__":
    main()
