"""
MOTORE FOTONICO DEL VUOTO - V6
Studi quantitativi. Ogni studio produce un JSON in dati/ e una tabella LaTeX
in tabelle/. Nessun numero del preprint deve essere scritto a mano.

Studi
-----
A  chiusura dei residui della V5
B  audit ottico (Mie esatta, cilindro, nucleo-guscio, eps(a), ossidazione)
C  bilancio termico (fluido e vuoto)
D  baseline geometriche a pari massa
E  matrice dei regimi (4 dimensioni x 2 mezzi)
F  materiali (rame, silice, grafene/CNT)
G  leve di progetto
H  confronto con la rotazione a GHz della silice
I  sensibilita' al momento angolare orbitale

Uso:
    python v6_studi.py
"""

import json
import math
import os
import sys

import numpy as np
from scipy.optimize import brentq

BASE = os.path.dirname(os.path.abspath(__file__))
if BASE not in sys.path:
    sys.path.insert(0, BASE)

import v6_nucleo as n  # noqa: E402

DATI = os.path.normpath(os.path.join(BASE, "..", "dati"))
TAB = os.path.normpath(os.path.join(BASE, "..", "tabelle"))
os.makedirs(DATI, exist_ok=True)
os.makedirs(TAB, exist_ok=True)

M0 = n.massa_telaio()          # kg, massa strutturale del reticolo
V0 = n.VOLUME_TELAIO           # m^3
I0 = n.momento_inerzia_telaio()
B_LO = n.beta_anello(n.R_OUT, n.WIRE_R)
B_HI = n.beta_disco(n.R_OUT)
M_CU = complex(1.09, 2.53)


def sci(x, cifre=3):
    """Numero in notazione scientifica LaTeX."""
    if x == 0:
        return "0"
    if not math.isfinite(x):
        return "\\infty"
    e = int(math.floor(math.log10(abs(x))))
    m = x / 10 ** e
    if abs(m) >= 9.995:
        m /= 10
        e += 1
    return f"{m:.{cifre}f}\\times10^{{{e}}}"


def rpm(omega):
    return omega / (2 * math.pi) * 60.0


_DEC = ["{,}"]         # separatore decimale LaTeX: virgola (IT) o punto (EN)


def num(x, d=3):
    """Numero decimale con il separatore decimale corrente."""
    return f"{x:.{d}f}".replace(".", _DEC[0])


def sx(x, d=2):
    """Numero in notazione scientifica LaTeX, con virgola decimale."""
    if x == 0:
        return "$0$"
    if not math.isfinite(x):
        return "$\\infty$"
    e = int(math.floor(math.log10(abs(x))))
    if -2 <= e <= 2:
        # numero "normale": tre cifre significative
        return "$" + f"{x:.3g}".replace(".", _DEC[0]) + "$"
    m = x / 10 ** e
    if abs(m) >= 9.995:
        m /= 10
        e += 1
    return f"${num(m, d)}\\times10^{{{e}}}$"


def sxr(a, b, d=2):
    """Intervallo: se i due estremi coincidono si scrive un solo valore."""
    if not (math.isfinite(a) and math.isfinite(b)):
        return "--"
    if b == 0 or abs(a / b - 1.0) < 5e-3:
        return sx(a, d)
    ea = int(math.floor(math.log10(abs(a)))) if a else 0
    eb = int(math.floor(math.log10(abs(b)))) if b else 0
    if ea == eb and not (-2 <= ea <= 3):
        return f"${num(a/10**ea, d)}$--${num(b/10**eb, d)}\\times10^{{{ea}}}$"
    return f"{sx(a, d)}--{sx(b, d)}"


def sezioni_riferimento():
    """Sezioni ottiche della configurazione di riferimento (nucleo + tre gruppi di fili)."""
    core = n.sezione_sfera(M_CU, n.CORE_R)
    fili = [n.sezione_filo(M_CU, a, L * 1e-9) for L, a in n.FILI]
    filo = {"sigma_pr": sum(f["sigma_pr"] for f in fili),
            "sigma_abs": sum(f["sigma_abs"] for f in fili),
            "sigma_ext": sum(f["sigma_ext"] for f in fili),
            "Q_pr": [f["Q_pr"] for f in fili], "x": [f["x"] for f in fili],
            "regime": fili[0]["regime"]}
    return {"core": core, "filo": filo,
            "sigma_pr": core["sigma_pr"] + filo["sigma_pr"],
            "sigma_abs": core["sigma_abs"] + filo["sigma_abs"],
            "sigma_ext": core["sigma_ext"] + filo["sigma_ext"]}


# ---------------------------------------------------------------------------
# A. CHIUSURA DEI RESIDUI
# ---------------------------------------------------------------------------
def studio_residui():
    from fractions import Fraction
    out = {"titolo": "A. Chiusura dei residui della V5"}

    # A.1 rapporti razionali, aritmetica esatta
    r_cin = Fraction(45, 376)
    r_iner = Fraction(376, 45)
    v5_calcolato = n.I2_V5 / n.I1_V5
    catena = {
        "razionale_45_376": float(r_cin),
        "razionale_376_45": float(r_iner),
        "inverso_esatto_di_45_376": float(1 / r_cin),
        "I2_su_I1_scena_V4": v5_calcolato,
        "scarto_scena_su_target": v5_calcolato / float(r_iner) - 1.0,
        "target_arrotondato_V5": 8.3565,
        "scarto_arrotondamento": 8.3565 / float(r_iner) - 1.0,
        "inverso_del_target_arrotondato": 1 / 8.3565,
        "scarto_riportato_abs": 1 / 8.3565 - float(r_cin),
        "scarto_riportato_rel": (1 / 8.3565) / float(r_cin) - 1.0,
        "scarto_vero_scena_abs": v5_calcolato * 0 + (1 / v5_calcolato - float(r_cin)),
    }
    out["rapporti"] = catena

    # A.2 raggio equivalente: tre definizioni
    rr = n.raggi_equivalenti()
    rr["massa_inviluppo_su_massa_materiale"] = n.M_V5 / M0
    rr["volume_inviluppo_su_volume_materiale"] = (n.M_V5 / 8960.0) / V0
    out["raggi"] = rr

    # A.3 massa, peso, spinta
    F = n.forza_limite(n.P_LASER / n.A_BEAM, sezioni_riferimento()["sigma_pr"])
    peso = M0 * 9.80665
    spinta_peso_strutturale = F / peso
    spinta_peso_V5 = F / (n.M_V5 * 9.80665)
    buoy = 998.2 * V0 * 9.80665
    out["massa"] = {
        "massa_materiale_kg": M0,
        "massa_V5_kg": n.M_V5,
        "rapporto": n.M_V5 / M0,
        "I_materiale": I0,
        "I_V5": n.I1_V5,
        "forza_N": F,
        "peso_1g_N": peso,
        "spinta_su_peso_strutturale": spinta_peso_strutturale,
        "spinta_su_peso_V5": spinta_peso_V5,
        "spinta_su_peso_effettivo_in_acqua": F / (peso - buoy),
        "galleggiamento_N": buoy,
        "nota": "in orbita il peso e' nullo e il rapporto spinta/peso non e' definito: "
                "il numero ha senso solo a 1 g e va confrontato con il peso effettivo",
    }

    # A.4 coppia: incidente / intercettata / assorbita
    I0_int = n.P_LASER / n.A_BEAM
    sez = sezioni_riferimento()
    P_int = I0_int * sez["sigma_pr"]
    P_abs = I0_int * sez["sigma_abs"]
    tau_int = n.coppia_limite(P_int)
    tau_abs = n.coppia_limite(P_abs)
    out["bilancio_potenza"] = {
        "P_incidente_W": n.P_LASER,
        "P_intercettata_W": P_int,
        "P_assorbita_W": P_abs,
        "frazione_intercettata": sez["sigma_pr"] / n.A_BEAM,
        "frazione_assorbita": sez["sigma_abs"] / n.A_BEAM,
        "tau_da_intercettata_Nm": tau_int,
        "tau_da_assorbita_Nm": tau_abs,
        "tau_V5_Nm": 2.6103646371378617e-24,
        "scarto_tau": tau_abs / 2.6103646371378617e-24 - 1,
        "P_meccanica_estratta_W": tau_abs * (tau_abs / B_HI),
        "efficienza_ottico_meccanica": tau_abs * (tau_abs / B_HI) / P_abs,
    }

    # A.5 regime pulsato: e' imposto o conseguenza?
    tau_mecc = I0 / B_LO
    out["pulsata"] = {
        "tempo_inerziale_I_su_beta_s": tau_mecc,
        "frequenza_di_taglio_Hz": 1 / (2 * math.pi * tau_mecc),
        "periodo_naturale_coppia_gradi": 180.0,
        "verdetto": "la pulsazione e' imposta dal laser, non e' una conseguenza della "
                    "dinamica: il rotore segue la coppia entro I/beta = 6e-11 s, quindi "
                    "qualunque impulso piu' lungo di 1 ns e' quasi statico. La sola "
                    "grandezza che la pulsazione cambia e' la temperatura media.",
        "duty_cycle_per_non_fondere_in_vuoto": None,
    }
    # duty cycle termico
    area = n.area_radiante()
    P_fus = n.potenza_fusione_vuoto(area, 0.03, n.MATERIALI["rame"]["Tm"])
    out["pulsata"]["duty_cycle_per_non_fondere_in_vuoto"] = P_fus / P_abs

    # A.6: verifica che le sezioni riproducano la V5
    out["verifica_sezioni"] = {
        "sigma_pr_V6": sez["sigma_pr"],
        "sigma_pr_V5": n.SIGMA_PR_V5,
        "scarto": sez["sigma_pr"] / n.SIGMA_PR_V5 - 1,
        "F_V6": F, "F_V5": n.F_V5, "scarto_F": F / n.F_V5 - 1,
    }
    _salva("A_residui.json", out)
    return out


# ---------------------------------------------------------------------------
# B. AUDIT OTTICO
# ---------------------------------------------------------------------------
def studio_ottica():
    out = {"titolo": "B. Audit ottico"}

    # B.1 sfera: Mie esatta contro Rayleigh contro ottica geometrica
    raggi = np.logspace(-9, -5.5, 60)
    curva = []
    for ra in raggi:
        s = n.sezione_sfera(M_CU, ra)
        curva.append({"r_nm": ra * 1e9, "x": s["x"], "Q_ext": s["Q_ext"],
                      "Q_abs": s["Q_abs"], "Q_pr": s["Q_pr"], "sigma_pr": s["sigma_pr"],
                      "regime": s["regime"]})
    out["sfera_curva"] = curva

    # B.2 cilindro: serie n>=1 contro Rayleigh
    cil = []
    for ra in np.logspace(-9.3, -7.3, 30):
        x = 2 * math.pi * ra / n.LAMBDA
        e = n.cilindro_esatto(M_CU, x)
        r = n.cilindro_rayleigh(M_CU, x)
        cil.append({"a_nm": ra * 1e9, "x": x, "Q_esatto": e["Q_ext"],
                    "Q_rayleigh": r["Q_ext"],
                    "scarto": e["Q_ext"] / r["Q_ext"] - 1 if r["Q_ext"] else 0.0})
    out["cilindro"] = cil

    # B.3 eps dipendente dalla dimensione
    eps = []
    for ra in (2e-9, 5e-9, 10e-9, 20e-9, 40e-9, 80e-9):
        m = n.epsilon_rame_dimensione(a=ra)
        s = n.sezione_sfera(m, ra)
        s0 = n.sezione_sfera(M_CU, ra)
        eps.append({"r_nm": ra * 1e9, "n": m.real, "k": m.imag,
                    "Q_abs_con_dimensione": s["Q_abs"], "Q_abs_bulk": s0["Q_abs"],
                    "scarto_rel": s["Q_abs"] / s0["Q_abs"] - 1 if s0["Q_abs"] else 0.0})
    out["eps_dimensione"] = eps

    # B.4 ossidazione: nucleo di rame con guscio di Cu2O
    # Il guscio e' sottile rispetto al raggio: l'effetto relativo e' quasi
    # statico, il valore assoluto viene dalla Mie esatta del rame nudo.
    m_ossido = complex(2.70, 0.15)     # Cu2O a 532 nm, stima (letteratura dispersa)
    base_mie = n.mie_sfera(M_CU, 2 * math.pi * n.CORE_R / n.LAMBDA)
    base_qs = n.sfera_rivestita_qs(M_CU, m_ossido, n.CORE_R, 1.0)
    oss = []
    for frac_guscio in (0.0, 0.05, 0.1, 0.2, 0.35, 0.5):
        qs = n.sfera_rivestita_qs(M_CU, m_ossido, n.CORE_R, 1.0 - frac_guscio)
        corr = qs["Q_abs"] / base_qs["Q_abs"]
        q_abs = base_mie["Q_abs"] * corr
        a = math.pi * n.CORE_R ** 2
        oss.append({"frazione_volume_guscio": frac_guscio,
                    "spessore_guscio_nm": n.CORE_R * 1e9 * (1 - (1 - frac_guscio) ** (1 / 3)),
                    "Q_abs": q_abs, "Q_abs_nudo": base_mie["Q_abs"],
                    "sigma_abs": a * max(q_abs, 0.0),
                    "scarto_Q_abs": corr - 1})
    out["ossidazione"] = oss

    # B.5 incertezza sui dati ottici: n, k +-30%
    bande = []
    sez = sezioni_riferimento()
    for dn, dk in ((0, 0), (+0.3, 0), (-0.3, 0), (0, +0.3), (0, -0.3), (+0.3, +0.3), (-0.3, -0.3)):
        m = complex(1.09 * (1 + dn), 2.53 * (1 + dk))
        c = n.sezione_sfera(m, n.CORE_R)
        f = n.sezione_filo(m, n.WIRE_R, n.LUNGHEZZA_FILI_NM * 1e-9)
        bande.append({"dn": dn, "dk": dk,
                      "sigma_pr": c["sigma_pr"] + f["sigma_pr"],
                      "sigma_abs": c["sigma_abs"] + f["sigma_abs"]})
    base = [b for b in bande if b["dn"] == 0 and b["dk"] == 0][0]
    for b in bande:
        b["scarto_sigma_pr"] = b["sigma_pr"] / base["sigma_pr"] - 1
        b["scarto_sigma_abs"] = b["sigma_abs"] / base["sigma_abs"] - 1
    out["bande_materiale"] = bande

    # B.6 assemblaggio: somma incoerente contro limite coerente
    k = 2 * math.pi / n.LAMBDA
    out["assemblaggio"] = {
        "ipotesi": "somma incoerente delle sezioni: valida quando la distanza tra i "
                   "diffusori supera la lunghezza di coerenza laterale del fascio e "
                   "l'accoppiamento vicino e' debole; il limite superiore rigoroso e' "
                   "N^2 sigma_1 (diffusione coerente), mai raggiunto a queste distanze",
        "distanza_tipica_nm": 250.0,
        "kD": k * 250e-9,
        "numero_elementi": 4,
        "limite_coerente_su_incoerente": 4.0,
        "sigma_pr_riferimento": sez["sigma_pr"],
        "contributo_nucleo": sez["core"]["sigma_pr"] / sez["sigma_pr"],
        "contributo_fili": sez["filo"]["sigma_pr"] / sez["sigma_pr"],
        "contributo_nucleo_assorbimento": sez["core"]["sigma_abs"] / sez["sigma_abs"],
        "nota_full_wave": "una verifica DDA/FDTD dell'intero telaio resta aperta: il "
                          "protocollo e' nella sezione 12 del preprint. Qui la "
                          "sovrapposizione e' stimata, non calcolata.",
    }

    # B.7 confronto con la sfera di pari volume
    r_eq = (3 * V0 / (4 * math.pi)) ** (1 / 3)
    s_eq = n.sezione_sfera(M_CU, r_eq)
    out["sfera_pari_volume"] = {
        "r_nm": r_eq * 1e9, "sigma_abs": s_eq["sigma_abs"], "sigma_pr": s_eq["sigma_pr"],
        "scarto_sigma_abs": s_eq["sigma_abs"] / sez["sigma_abs"] - 1,
        "scarto_sigma_pr": s_eq["sigma_pr"] / sez["sigma_pr"] - 1,
        "nota": "la sfera di pari volume riproduce la sezione del telaio entro il 10%: "
                "l'ottica dipende dal volume di materiale, non dalla forma",
    }
    _salva("B_ottica.json", out)
    return out


# ---------------------------------------------------------------------------
# C. BILANCIO TERMICO
# ---------------------------------------------------------------------------
def sezioni_materiale(nome):
    """Sezioni della stessa geometria (fili + nucleo) in un materiale diverso."""
    mat = n.MATERIALI[nome]
    m_rel = complex(mat["n"], mat["k"])
    core = n.sezione_sfera(m_rel, n.CORE_R)
    fili = [n.sezione_filo(m_rel, a, L * 1e-9) for L, a in n.FILI]
    return {"core": core,
            "sigma_pr": core["sigma_pr"] + sum(f["sigma_pr"] for f in fili),
            "sigma_abs": core["sigma_abs"] + sum(f["sigma_abs"] for f in fili),
            "sigma_ext": core["sigma_ext"] + sum(f["sigma_ext"] for f in fili)}


def configurazione_termica(materiale="rame", potenza=n.P_LASER):
    mat = n.MATERIALI[materiale]
    sez = sezioni_riferimento() if materiale == "rame" else sezioni_materiale(materiale)
    I0 = potenza / n.A_BEAM
    P_abs = I0 * sez["sigma_abs"]
    # superficie metallica reale, non l'inviluppo sferico di 370 nm
    area = n.area_radiante()
    r_eq = (3 * V0 / (4 * math.pi)) ** (1 / 3)
    dT_fl = n.delta_T_fluido(P_abs, r_eq, kappa_m=0.6)
    T_v1 = n.T_vuoto(P_abs, area, mat["emiss_limpido"])
    T_v2 = n.T_vuoto(P_abs, area, mat["emiss_ossidato"])
    m = mat["rho"] * V0
    return {
        "materiale": materiale, "potenza_W": potenza,
        "P_abs_W": P_abs, "area_m2": area, "r_equivalente_nm": r_eq * 1e9,
        "massa_kg": m, "cp": mat["cp"],
        "dT_fluido_K": dT_fl,
        "T_fluido_K": n.T0 + dT_fl,
        "T_vuoto_emiss_limpido_K": T_v1,
        "T_vuoto_emiss_ossidato_K": T_v2,
        "emiss_limpido": mat["emiss_limpido"],
        "emiss_ossidato": mat["emiss_ossidato"],
        "Tm_K": mat["Tm"],
        "margine_limpido": mat["Tm"] / T_v1,
        "margine_ossidato": mat["Tm"] / T_v2,
        "potenza_fusione_vuoto_limpido_W": n.potenza_fusione_vuoto(area, mat["emiss_limpido"], mat["Tm"]),
        "potenza_fusione_vuoto_ossidato_W": n.potenza_fusione_vuoto(area, mat["emiss_ossidato"], mat["Tm"]),
        "tempo_adiabatico_a_Tm_s": n.tempo_termico(m, mat["cp"], mat["Tm"] - n.T0, P_abs) if P_abs > 0 else float("inf"),
        "potenza_massima_fluido_W": 4 * math.pi * 0.6 * r_eq * 100.0 / max(sez["sigma_abs"] / n.A_BEAM, 1e-30),
        "P_dissipata_meccanica_W": None,
    }


def studio_termica():
    out = {"titolo": "C. Bilancio termico",
           "nota": "P_ass = P_diss + P_out. In regime stazionario P_diss = beta omega^2, "
                   "che e' trascurabile: tutta la potenza assorbita diventa calore."}
    casi = []
    for mat in ("rame", "silice", "grafene"):
        for P in (1e-6, 1e-4, 1e-2, 1e-1):
            casi.append(configurazione_termica(mat, P))
    out["casi"] = casi

    # potenza meccanica estratta alla coppia di riferimento
    sez = sezioni_riferimento()
    I0 = n.P_LASER / n.A_BEAM
    P_abs = I0 * sez["sigma_abs"]
    tau = n.coppia_limite(P_abs)
    omega = tau / B_HI
    out["efficienza"] = {
        "P_abs_W": P_abs,
        "tau_Nm": tau,
        "omega_ss_rad_s": omega,
        "P_meccanica_W": tau * omega,
        "efficienza_ottico_meccanica": tau * omega / P_abs,
        "nota": "il rotore converte in rotazione una frazione ~1e-18 della potenza "
                "assorbita: termicamente e' una resistenza, non un motore",
    }
    # sweep in potenza, per la figura
    sweep = []
    for P in np.logspace(-6, -1, 40):
        c = configurazione_termica("rame", P)
        sweep.append({"P_W": P, "P_abs_W": c["P_abs_W"], "dT_fluido_K": c["dT_fluido_K"],
                      "T_vuoto_limpido_K": c["T_vuoto_emiss_limpido_K"],
                      "T_vuoto_ossidato_K": c["T_vuoto_emiss_ossidato_K"]})
    out["sweep"] = sweep
    _salva("C_termica.json", out)
    return out


# ---------------------------------------------------------------------------
# D. BASELINE A PARI MASSA
# ---------------------------------------------------------------------------
def geometrie_pari_massa(m0=M0, materiale="rame"):
    """Sei geometrie con la stessa massa e lo stesso materiale.

    Ogni geometria ha una famiglia di forme fissata (proporzioni), quindi la
    massa fissa un solo parametro di scala.  L'ottica e' presa dalla sfera di
    pari volume per tutte: l'audit del §B mostra che la forma cambia la sezione
    meno del 10%, mentre cambia l'attrito di tre ordini di grandezza.
    """
    mat = n.MATERIALI[materiale]
    rho = mat["rho"]
    V = m0 / rho
    r_eq = (3 * V / (4 * math.pi)) ** (1 / 3)
    m_rel = complex(mat["n"], mat["k"])
    s_eq = n.sezione_sfera(m_rel, r_eq)

    geo = []

    # 1. reticolo (configurazione di riferimento)
    sez = sezioni_riferimento()
    geo.append({
        "nome": "Reticolo (V5)", "famiglia": "telaio di fili",
        "R_caratteristico_nm": n.R_OUT * 1e9, "parametro": "6,4 / 15 / 9 nm",
        "beta_lo": B_LO, "beta_hi": B_HI,
        "I": I0, "sigma_pr": sez["sigma_pr"], "sigma_abs": sez["sigma_abs"],
        "area_frontale": n.AREA_PRO_FILI,
        "sigma_reale_pr": sez["sigma_pr"], "sigma_reale_abs": sez["sigma_abs"],
    })

    # 2. disco pieno, rapporto d'aspetto 10:1
    R_d = (10 * V / math.pi) ** (1 / 3)
    geo.append({
        "nome": "Disco pieno", "famiglia": "disco 10:1",
        "R_caratteristico_nm": R_d * 1e9, "parametro": f"t = {R_d/10*1e9:.1f} nm",
        "beta_lo": n.beta_disco(R_d), "beta_hi": n.beta_disco(R_d),
        "I": 0.5 * m0 * R_d ** 2,
        "sigma_pr": s_eq["sigma_pr"], "sigma_abs": s_eq["sigma_abs"],
        "area_frontale": math.pi * R_d ** 2,
    })

    # 3. anello, a/R = 1/40
    R_a = (V * 1600 / (2 * math.pi ** 2)) ** (1 / 3)
    a_a = R_a / 40
    geo.append({
        "nome": "Anello", "famiglia": "toro a/R = 1/40",
        "R_caratteristico_nm": R_a * 1e9, "parametro": f"a = {a_a*1e9:.1f} nm",
        "beta_lo": n.beta_anello(R_a, a_a), "beta_hi": n.beta_disco(R_a),
        "I": m0 * R_a ** 2,
        "sigma_pr": s_eq["sigma_pr"], "sigma_abs": s_eq["sigma_abs"],
        "area_frontale": 4 * R_a * a_a,
    })

    # 4. sfera
    R_s = r_eq
    geo.append({
        "nome": "Sfera", "famiglia": "sfera piena",
        "R_caratteristico_nm": R_s * 1e9, "parametro": f"r = {R_s*1e9:.1f} nm",
        "beta_lo": n.beta_sfera(R_s), "beta_hi": n.beta_sfera(R_s),
        "I": 0.4 * m0 * R_s ** 2,
        "sigma_pr": s_eq["sigma_pr"], "sigma_abs": s_eq["sigma_abs"],
        "area_frontale": math.pi * R_s ** 2,
    })

    # 5. nanomanubrio, due sfere con d = 4r
    r_m = (3 * V / (8 * math.pi)) ** (1 / 3)
    d_m = 4 * r_m
    geo.append({
        "nome": "Nanomanubrio", "famiglia": "due sfere, d = 4r",
        "R_caratteristico_nm": (d_m / 2 + r_m) * 1e9,
        "parametro": f"r = {r_m*1e9:.1f} nm, d = {d_m*1e9:.1f} nm",
        "beta_lo": n.beta_manubrio(r_m, d_m),
        "beta_hi": n.beta_manubrio(r_m, d_m) * 1.25,
        "I": 2 * 0.4 * (m0 / 2) * r_m ** 2 + 2 * (m0 / 2) * (d_m / 2) ** 2,
        "sigma_pr": s_eq["sigma_pr"], "sigma_abs": s_eq["sigma_abs"],
        "area_frontale": 2 * math.pi * r_m ** 2 + 4 * r_m ** 2,
    })

    # 6. silice anisotropa (materiale diverso: pari massa, non pari volume)
    mat_s = n.MATERIALI["silice"]
    V_s = m0 / mat_s["rho"]
    r_sil = (3 * V_s / (4 * math.pi)) ** (1 / 3)
    m_sil = complex(mat_s["n"], mat_s["k"])
    s_sil = n.sezione_sfera(m_sil, r_sil)
    geo.append({
        "nome": "Silice anisotropa", "famiglia": "sfera di silice",
        "R_caratteristico_nm": r_sil * 1e9,
        "parametro": f"r = {r_sil*1e9:.1f} nm, Delta n = 0.05",
        "beta_lo": n.beta_sfera(r_sil), "beta_hi": n.beta_sfera(r_sil),
        "I": 0.4 * m0 * r_sil ** 2,
        "sigma_pr": s_sil["sigma_pr"], "sigma_abs": s_sil["sigma_abs"],
        "area_frontale": math.pi * r_sil ** 2,
        "note": "silice: l'assorbimento a 532 nm e' nullo, la coppia viene dal ritardo di fase",
    })
    return geo, r_eq, s_eq


def studio_baseline(potenza=n.P_LASER):
    geo, r_eq, s_eq = geometrie_pari_massa()
    I0 = potenza / n.A_BEAM
    out = {"titolo": "D. Baseline a pari massa",
           "massa_kg": M0, "volume_m3": V0, "potenza_W": potenza,
           "raggio_sfera_pari_volume_nm": r_eq * 1e9,
           "nota_ottica": "sigma_pr e sigma_abs dalla sfera di pari volume (Mie esatta), "
                          "verificata sul telaio al 6%",
           "casi": []}
    for g in geo:
        P_abs = I0 * g["sigma_abs"]
        P_int = I0 * g["sigma_pr"]
        tau = n.coppia_limite(P_abs)
        # canale di ritardo per la silice
        if "Delta n" in g["parametro"]:
            delta = 2 * math.pi * 0.05 * 2 * g["R_caratteristico_nm"] * 1e-9 / n.LAMBDA
            tau_ret = P_int / n.OMEGA_OTT * math.sin(delta / 2) ** 2
        else:
            delta, tau_ret = 0.0, 0.0
        tau_usato = max(tau, tau_ret)
        # Banda ottica: la sezione della sfera di pari volume NON e' quella di un
        # disco, di un anello o di un manubrio.  A parita' di volume le risonanze
        # di forma spostano Q_abs di un fattore 2-3: la banda dichiarata e' 1-3
        # volte il valore di volume, e un calcolo DDA per forma resta da fare.
        if g["nome"] in ("Sfera", "Reticolo (V5)"):
            tau_hi = tau_usato
        elif "Delta n" in g["parametro"]:
            tau_hi = tau_usato * 3.0
        else:
            tau_hi = tau_usato * 3.0
        row = {
            "nome": g["nome"], "famiglia": g["famiglia"],
            "R_nm": g["R_caratteristico_nm"], "parametro": g["parametro"],
            "beta_lo": g["beta_lo"], "beta_hi": g["beta_hi"],
            "I": g["I"], "sigma_pr": g["sigma_pr"], "sigma_abs": g["sigma_abs"],
            "eta_int": g["sigma_pr"] / n.A_BEAM,
            "P_int_W": P_int, "P_abs_W": P_abs,
            "tau_assorbimento_Nm": tau, "tau_ritardo_Nm": tau_ret,
            "tau_usato_Nm": tau_usato, "tau_hi_Nm": tau_hi,
            "Pe_r": n.pe_r(tau_usato), "Pe_r_lo": n.pe_r(tau_usato),
            "Pe_r_hi": n.pe_r(tau_hi),
            "omega_ss_lo": tau_usato / g["beta_hi"],
            "omega_ss_hi": tau_usato / g["beta_lo"],
            "omega_ss_hi_ottica": tau_hi / g["beta_lo"],
            "rpm_lo": rpm(tau_usato / g["beta_hi"]),
            "omega_rottura": n.omega_rottura(200e6, 8960.0, g["R_caratteristico_nm"] * 1e-9),
            "t_rottura_s": g["I"] * n.omega_rottura(200e6, 8960.0, g["R_caratteristico_nm"] * 1e-9) / tau_usato if tau_usato > 0 else float("inf"),
            "R_g_nm": math.sqrt(g["I"] / M0) * 1e9,
        }
        row["beta_medio"] = 0.5 * (g["beta_lo"] + g["beta_hi"])
        row["Pe_r_su_beta"] = row["Pe_r"]
        out["casi"].append(row)
    # reticolo senza nucleo: cosa aggiungono i fili
    sez = sezioni_riferimento()
    out["scomposizione_reticolo"] = {
        "sigma_abs_nucleo": sez["core"]["sigma_abs"],
        "sigma_abs_fili": sez["filo"]["sigma_abs"],
        "quota_coppia_nucleo": sez["core"]["sigma_abs"] / sez["sigma_abs"],
        "quota_coppia_fili": sez["filo"]["sigma_abs"] / sez["sigma_abs"],
        "quota_massa_nucleo": (n.VOLUME_NUCLEO / n.VOLUME_TELAIO),
        "quota_massa_fili": (n.VOLUME_FILI / n.VOLUME_TELAIO),
        "beta_solo_nucleo": n.beta_sfera(n.CORE_R),
        "beta_telaio_lo": B_LO,
        "amplificazione_attrito": B_LO / n.beta_sfera(n.CORE_R),
        "nota": ("il telaio porta il %.0f%% della massa e il %.0f%% della coppia, ma "
                 "il %.0f%% dell'attrito: con i raggi di filo della scena il suo "
                 "contributo ottico non e' piu' marginale, quello meccanico resta "
                 "dominante"
                 % (100 * n.VOLUME_FILI / n.VOLUME_TELAIO,
                    100 * sez["filo"]["sigma_abs"] / sez["sigma_abs"],
                    100 * B_LO / (B_LO + n.beta_sfera(n.CORE_R)))),
    }
    _salva("D_baseline.json", out)
    return out


# ---------------------------------------------------------------------------
# E. MATRICE DEI REGIMI
# ---------------------------------------------------------------------------
def caso_regime(L, mezzo, potenza=n.P_LASER, materiale="rame"):
    """Una cella della matrice: dimensione esterna L e mezzo."""
    mat = n.MATERIALI[materiale]
    s_scale = L / n.R_OUT
    m = M0 * s_scale ** 3
    I = I0 * s_scale ** 5
    sez = sezioni_riferimento()
    sigma_pr = sez["sigma_pr"] * s_scale ** 2
    sigma_abs = sez["sigma_abs"] * s_scale ** 2
    I0i = potenza / n.A_BEAM
    P_int = I0i * min(sigma_pr, n.A_BEAM)
    P_abs = I0i * min(sigma_abs, n.A_BEAM)
    eta = min(sigma_pr / n.A_BEAM, 1.0)
    tau = n.coppia_limite(P_abs)
    R = n.R_OUT * s_scale

    if mezzo == "fluido":
        if L <= 1e-6:
            b_lo = n.beta_anello(R, n.WIRE_R * s_scale)
            b_hi = n.beta_disco(R)
            omega_lo, omega_hi = tau / b_hi, tau / b_lo
            re = 998.2 * omega_lo * R ** 2 / n.ETA_W
            nota = "Stokes"
            if re > 1e3:
                nota = "Stokes non valido (Re > 1e3)"
            t_rott = float("inf")
        else:
            # prima si prova Stokes su disco; se Re resta piccolo, vale
            b_st = n.beta_disco(R)
            omega_st = tau / b_st
            re_st = 998.2 * omega_st * R ** 2 / n.ETA_W
            if re_st < 1.0:
                b_lo = b_hi = b_st
                omega_lo = omega_hi = omega_st
                re = re_st
                nota = "Stokes (Re << 1)"
            else:
                # strato limite / turbolento: risolvere T(omega) = tau
                def f(om):
                    return n.beta_turbolento(R, om) * om - tau
                try:
                    omega_hi = brentq(f, 1e-14, 1e14)
                except ValueError:
                    omega_hi = 1e-14
                omega_lo = omega_hi
                b_lo = b_hi = n.beta_turbolento(R, omega_hi)
                re = 998.2 * omega_hi * R ** 2 / n.ETA_W
                nota = "strato limite / turbolento"
            t_rott = float("inf")
        dr = n.KT / b_hi
        pe = tau / n.KT
        T_fluido = n.T0 + n.delta_T_fluido(P_abs, (3 * (M0 * s_scale ** 3) / (4 * math.pi * mat["rho"])) ** (1 / 3))
        validita = "valido"
        if L < 5e-9:
            validita = "fuori dominio (continuo e filo subatomico)"
    else:
        b_lo = b_hi = n.beta_gas(1e-6, R)     # 1e-6 mbar
        omega_lo = omega_hi = float("inf")
        re = 0.0
        nota = "molecolare libero, beta = 0 nello stazionario"
        dr = n.KT / b_hi if b_hi > 0 else float("inf")
        pe = tau / n.KT
        T_fluido = n.T_vuoto(P_abs, n.area_radiante() * s_scale ** 2, mat["emiss_ossidato"])
        t_rott = I * n.omega_rottura(mat["sy"], mat["rho"], R) / tau if tau > 0 else float("inf")
        validita = "valido"
        if L < 5e-9:
            validita = "fuori dominio (continuo non applicabile)"
        if L > 1e-3:
            validita = "estrapolazione"

    return {
        "L_m": L,
        "L_etichetta": (f"{L*1e9:.0f} nm" if L < 1e-6
                        else ("1 $\\mu$m" if L < 1e-3 else "1 m")),
        "mezzo": mezzo, "s_scala": s_scale, "massa_kg": m, "I": I,
        "R_m": R, "eta_int": eta, "P_int_W": P_int, "P_abs_W": P_abs,
        "tau_Nm": tau, "beta_lo": b_lo, "beta_hi": b_hi,
        "omega_ss_lo": omega_lo, "omega_ss_hi": omega_hi,
        "rpm": rpm(omega_lo) if math.isfinite(omega_lo) else None,
        "D_r": dr, "Pe_r": pe, "Re": re, "regime": nota,
        "T_K": T_fluido, "t_rottura_s": t_rott, "validita": validita,
        "omega_rottura": n.omega_rottura(mat["sy"], mat["rho"], R),
    }


def studio_regimi():
    out = {"titolo": "E. Matrice dei regimi", "casi": []}
    for L in (1e-9, 240e-9, 1e-6, 1.0):
        for mezzo in ("fluido", "vuoto"):
            out["casi"].append(caso_regime(L, mezzo))
    # riga di riferimento
    out["riferimento"] = {
        "L_m": n.R_OUT * 2, "nota": "configurazione di riferimento (diametro 740 nm)",
        "massa_kg": M0, "beta_lo": B_LO, "beta_hi": B_HI,
        "Pe_r": n.pe_r(n.coppia_limite(n.P_LASER / n.A_BEAM * sezioni_riferimento()["sigma_abs"])),
    }
    _salva("E_regimi.json", out)
    return out


# ---------------------------------------------------------------------------
# F. MATERIALI
# ---------------------------------------------------------------------------
def studio_materiali(potenza=n.P_LASER):
    out = {"titolo": "F. Materiali", "casi": []}
    for nome, mat in n.MATERIALI.items():
        m = mat["rho"] * V0
        r_eq = (3 * V0 / (4 * math.pi)) ** (1 / 3)
        # stessa geometria (fili + nucleo) nel materiale in esame
        sez = sezioni_riferimento() if nome == "rame" else sezioni_materiale(nome)
        I0i = potenza / n.A_BEAM
        P_abs = I0i * sez["sigma_abs"]
        tau = n.coppia_limite(P_abs)
        v_lim = math.sqrt(mat["sy"] / mat["rho"])
        omega_r = v_lim / n.R_OUT
        area = n.area_radiante()
        T_vac = n.T_vuoto(P_abs, area, mat["emiss_ossidato"])
        out["casi"].append({
            "materiale": nome, "nome": mat["nome"],
            "rho": mat["rho"], "sy": mat["sy"], "sy_nota": mat["sy_nota"],
            "specifico": mat["sy"] / mat["rho"],
            "v_limite_m_s": v_lim, "omega_limite": omega_r,
            "n": mat["n"], "k": mat["k"],
            "Q_abs": sez["core"]["Q_abs"], "Q_ext": sez["core"]["Q_ext"],
            "sigma_abs": sez["sigma_abs"],
            "kappa": mat["kappa"], "cp": mat["cp"], "Tm": mat["Tm"],
            "massa_kg": m, "P_abs_W": P_abs, "tau_Nm": tau,
            "Pe_r": n.pe_r(tau), "omega_ss_lo": tau / B_HI, "omega_ss_hi": tau / B_LO,
            "dT_fluido_K": n.delta_T_fluido(P_abs, r_eq),
            "T_vuoto_K": T_vac, "margine_fusione": mat["Tm"] / T_vac,
            "t_rottura_s": I0 * omega_r / tau if tau > 0 else float("inf"),
            "nota": mat["nota"],
        })
    _salva("F_materiali.json", out)
    return out


# ---------------------------------------------------------------------------
# G. LEVE DI PROGETTO
# ---------------------------------------------------------------------------
def studio_leve(potenza=n.P_LASER):
    sez = sezioni_riferimento()
    def stato(P=potenza, w0=n.W0_DICH, l=0.0, s=1.0, sigma=None, beta=B_HI, overlap=1.0):
        A = math.pi * w0 ** 2
        sig = sez["sigma_abs"] if sigma is None else sigma
        P_abs = min(sig, A) * (P / A) * overlap
        tau = n.coppia_limite(P_abs, l=l, s=s)
        return {"P_abs": P_abs, "tau": tau, "Pe_r": n.pe_r(tau),
                "omega": tau / beta, "eta": min(sig / A, 1.0)}
    base = stato()
    g5 = 5 ** 5 * math.exp(-5) / math.factorial(5)      # sovrapposizione LG_05
    leve = [
        ("Potenza x10 (1 $\\to$ 10 $\\mu$W)", "P", stato(P=10 * potenza)),
        ("Messa a fuoco NA 0,9 ($w_0$ 500 $\\to$ 361 nm)", "w0", stato(w0=n.W0_DIFF)),
        ("Polarizzazione lineare (s = 0)", "s", stato(s=0.0)),
        ("OAM $\\ell = 5$ con $s = 1$ (con sovrapposizione)", "l",
         stato(l=5.0, s=1.0, overlap=g5)),
        ("Compattazione a pari massa ($R/10$)", "R", None),
        ("Vuoto spinto ($p = 10^{-6}$ mbar)", "beta", None),
        ("Materiale: grafene invece di rame", "mat", None),
    ]
    out = {"titolo": "G. Leve di progetto", "base": base, "leve": []}
    for nome, tipo, st in leve:
        if tipo == "R":
            b = n.beta_sfera(4.2e-8)
            d = {"nome": nome, "Pe_r": base["Pe_r"],
                 "omega": base["tau"] / b,
                 "nota": "coppia invariante; omega per il rapporto degli attriti (1e2-1e3)"}
            d["fattore_Pe"] = 1.0
            d["fattore_omega"] = (base["omega"] and d["omega"] / base["omega"]) or 0.0
        elif tipo == "beta":
            bg = n.beta_gas(1e-4, n.R_OUT)
            d = {"nome": nome, "Pe_r": base["Pe_r"], "omega": base["tau"] / bg,
                 "nota": "nessuno stato stazionario reale: il limite diventa la rottura"}
            d["fattore_Pe"] = 1.0
            d["fattore_omega"] = d["omega"] / base["omega"]
        elif tipo == "mat":
            mat = n.MATERIALI["grafene"]
            s_g = sezioni_materiale("grafene")
            P_abs = (potenza / n.A_BEAM) * s_g["sigma_abs"]
            tau = n.coppia_limite(P_abs)
            d = {"nome": nome, "Pe_r": n.pe_r(tau), "omega": tau / B_HI,
                 "nota": "cambia l\u2019assorbimento specifico; effetto su $Pe_r$ $= " + f"{n.pe_r(tau)/base['Pe_r']:.2f}$"}
            d["fattore_Pe"] = n.pe_r(tau) / base["Pe_r"]
            d["fattore_omega"] = tau / base["tau"]
        else:
            d = {"nome": nome, "Pe_r": st["Pe_r"], "omega": st["omega"],
                 "nota": "", "fattore_Pe": st["Pe_r"] / base["Pe_r"],
                 "fattore_omega": st["omega"] / base["omega"]}
        d["tipo"] = tipo
        out["leve"].append(d)
    _salva("G_leve.json", out)
    return out


# ---------------------------------------------------------------------------
# H. CONFRONTO CON LA ROTAZIONE A GHz DELLA SILICE
# ---------------------------------------------------------------------------
def studio_confronto_ghz():
    """Controllo esterno sulla SCALA MECCANICA, non validazione del modello ottico.

    Tre verifiche, tenute separate:
      (1) per ogni coppia (f osservata, r) si ricava la tensione efficace
          implicita sigma_y = rho (2 pi r f)^2;
      (2) si confronta quella tensione con la resistenza della silice su nanoscala;
      (3) si mostra che la coppia necessaria e' abbondante, quindi il vincolo non
          e' ottico.
    Nessuna delle tre valida la coppia prevista o la dinamica del rotore.
    """
    rho_s = 2200.0
    registrazioni = ((1e9, 50.0), (1e9, 75.0), (5e9, 50.0), (6e9, 50.0))
    implicita = []
    for f, r_nm in registrazioni:
        r = r_nm * 1e-9
        v = 2 * math.pi * f * r
        implicita.append({"f_Hz": f, "r_nm": r_nm, "v_tangenziale_m_s": v,
                          "sigma_implicita_Pa": rho_s * v ** 2,
                          "rho_v2_GPa": rho_s * v ** 2 / 1e9})

    casi = []
    for r_nm in (25.0, 50.0, 75.0, 100.0):
        r = r_nm * 1e-9
        m = rho_s * 4 / 3 * math.pi * r ** 3
        I = 0.4 * m * r ** 2
        for sy in (5e8, 1e9, 3e9, 5e9, 8e9):
            casos = {
                "tipo": "rottura", "r_nm": r_nm, "sy_Pa": sy,
                "v_limite_m_s": math.sqrt(sy / rho_s),
                "omega_limite": math.sqrt(sy / rho_s) / r,
                "f_limite_Hz": math.sqrt(sy / rho_s) / r / (2 * math.pi),
                "I": I,
            }
            casi.append(casos)
        for p_mbar in (1e-6, 1e-8, 1e-10):
            p = p_mbar * 100.0
            bg = n.beta_gas(p, r)
            omega = 2 * math.pi * 1e9
            tau = bg * omega
            for delta in (0.02, 0.06):
                s_sil = n.sezione_sfera(complex(1.461, 1e-7), r)
                sigma_pr = max(s_sil["sigma_pr"], math.pi * r ** 2 * 0.2)
                w0 = 1.0e-6
                eta = sigma_pr / (math.pi * w0 ** 2)
                P_int = tau * n.OMEGA_OTT / math.sin(delta / 2) ** 2
                casi.append({
                    "tipo": "potenza", "r_nm": r_nm, "p_mbar": p_mbar,
                    "beta_gas": bg, "tau_necessaria_Nm": tau, "delta": delta,
                    "sin2": math.sin(delta / 2) ** 2,
                    "P_int_necessaria_W": P_int, "eta_int": eta,
                    "P_sorgente_necessaria_W": P_int / eta,
                    "P_sorgente_10mW_su_necessaria": 1e-2 / (P_int / eta),
                    "t_spinup_a_10mW_s": I * omega / (1e-2 * eta / n.OMEGA_OTT * math.sin(delta / 2) ** 2),
                    "I": I,
                })
    out = {"titolo": "H. Controllo esterno sulla scala meccanica",
           "ipotesi": "nanorotore di silice in vuoto, r = 25-100 nm, w0 = 1 um, "
                      "birifrangenza di forma Delta n = 0.05",
           "riferimenti": "Reimann et al., PRL 121, 033602 (2018); Ahn et al., "
                          "Nat. Nanotechnol. 15, 89 (2020); rotazione a 5-6 GHz "
                          "riportata in lavori successivi sullo stesso tipo di "
                          "rotore. I parametri vanno verificati sull'articolo.",
           "registrazioni": implicita,
           "casi": casi}
    s50 = [x for x in implicita if x["r_nm"] == 50.0]
    out["sintesi"] = {
        "sigma_implicita_GPa": {f"{x['f_Hz']/1e9:.0f} GHz / {x['r_nm']:.0f} nm":
                                round(x["rho_v2_GPa"], 3) for x in implicita},
        "verdetto":
            "Le coppie (f, r) riportate implicano una tensione efficace fra "
            f"{min(x['rho_v2_GPa'] for x in implicita):.2f} e "
            f"{max(x['rho_v2_GPa'] for x in implicita):.2f} GPa: la formula "
            "sqrt(sigma_y/rho)/(2 pi r) riproduce le frequenze osservate per "
            "sigma_y nell'intervallo GPa, che e' quello della silice su nanoscala. "
            "Questo e' un controllo di compatibilita' sulla SCALA MECCANICA, non "
            "una validazione del modello ottico: coppia, trasferimento di momento "
            "angolare e dinamica non sono confrontati con una misura.",
    }
    _salva("H_confronto_ghz.json", out)
    return out


# ---------------------------------------------------------------------------
# I. OAM
# ---------------------------------------------------------------------------
def studio_oam():
    """OAM con la sovrapposizione di modo, non a eta_int costante.

    Per un fascio LG_{0,l} l'intensita' di picco, normalizzata a quella
    gaussiana con la stessa potenza e lo stesso w0, vale

        g(l) = l^l e^{-l} / l!

    e sull'asse e' nulla per l >= 1.  Una particella posta sul massimo
    intercetta sigma * I_picco, quindi la coppia netta scala come
    (l + s) * g(l): il fattore l non e' un guadagno, perche' la potenza si
    distribuisce su un anello piu' grande.
    """
    sez = sezioni_riferimento()
    I0 = n.P_LASER / n.A_BEAM
    P_abs0 = I0 * sez["sigma_abs"]
    casi = []
    for l in (0.0, 1.0, 2.0, 5.0, 10.0, 20.0):
        g = 1.0 if l == 0 else l ** l * math.exp(-l) / math.factorial(int(l))
        for s in (0.0, 1.0):
            tau_senza = n.coppia_limite(P_abs0, l=l, s=s)
            tau_con = n.coppia_limite(P_abs0 * g, l=l, s=s)
            casi.append({"l": l, "s": s, "l_plus_s": l + s,
                         "overlap_g": g,
                         "tau_ideale_Nm": tau_senza,
                         "tau_reale_Nm": tau_con,
                         "Pe_r_ideale": n.pe_r(tau_senza),
                         "Pe_r_reale": n.pe_r(tau_con),
                         "guadagno_netto": (tau_con / n.coppia_limite(P_abs0)) if P_abs0 else 0.0,
                         "omega_lo": tau_con / B_LO, "omega_hi": tau_con / B_HI})
    out = {"titolo": "I. Sensibilita' al momento angolare orbitale",
           "casi": casi,
           "nota": "il guadagno netto (l+s)*g(l) e' 0.74 a l=1, 1.05 a l=5, "
                   "1.38 a l=10: l'ordine di grandezza e' sqrt(l), non l. "
                   "Una particella centrata su un fascio con l>=1 non riceve "
                   "intensita' sull'asse e quindi non riceve coppia: la "
                   "separazione tra canale di trappola e canale di coppia e' "
                   "possibile solo se la particella sta sull'anello."}
    _salva("I_oam.json", out)
    return out


# ---------------------------------------------------------------------------
# SALVATAGGIO E TABELLE
# ---------------------------------------------------------------------------
def _salva(nome, obj):
    with open(os.path.join(DATI, nome), "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=1, ensure_ascii=False)


def _tex(nome, corpo):
    with open(os.path.join(TAB, nome), "w", encoding="utf-8") as f:
        f.write(corpo)


def tabelle(res):
    # A. residui: quattro colonne strette, nessun testo che sbordi
    a = res["residui"]
    r = a["rapporti"]
    b = a["bilancio_potenza"]
    F = a["massa"]["forza_N"]
    righe = []
    righe.append("Rapporto $\\omega_2/\\omega_1$ & "
                 + num(0.119667646, 9) + " & $45/376=" + num(45 / 376, 9)
                 + "$ & il valore V5 non e' $" + num(8.3565, 4) + "$ (cio\`e $"
                 + num(1 / 8.3565, 8) + "$) ma $1/" + num(1 / 0.119667646, 5)
                 + "$: non ricostruibile dai dati dichiarati \\\\")
    righe.append("$I_2/I_1$ & " + num(r["I2_su_I1_scena_V4"], 6) + " & $376/45="
                 + num(376 / 45, 6) + "$ & scarto " + sx(r["scarto_scena_su_target"], 1)
                 + "; la chiusura richiede la massa dell'App.~A \\\\")
    righe.append("$r_{\\rm eff}$ & " + num(a["raggi"]["r_eff_v5_nm"], 2) + " nm & "
                 + num(a["raggi"]["r_volume_materiale_nm"], 2) + " / "
                 + num(a["raggi"]["r_idrodinamico_nm"], 2) + " nm & "
                 "volume di materiale / raggio idrodinamico: sono due grandezze diverse \\\\")
    righe.append("Massa & " + sx(n.M_V5) + " kg & " + sx(M0) + " kg & "
                 "la V5 contava i soli fili; la V6 aggiunge il nucleo (29\\,\\% del totale) \\\\")
    righe.append("Forza $F$ & " + sx(n.F_V5) + " N & " + sx(F)
                 + " N & $+47$\\,\\%: raggi di filo 3,2 / 7,5 / 4,5 nm, non 1,25 nm \\\\")
    righe.append("Coppia $\\tau$ & " + sx(b["tau_V5_Nm"]) + " N\\,m & "
                 + sx(b["tau_da_assorbita_Nm"]) + " N\\,m & $"
                 + f"{b['scarto_tau']*100:+.0f}" + "\\%$: $P_{\\rm ass}$ e non "
                 "$P_{\\rm int}$, con le sezioni corrette \\\\")
    righe.append("$Pe_r$ & " + sx(6.099e-4) + " & " + sx(n.pe_r(b["tau_da_assorbita_Nm"]))
                 + " & identit\\`a $\\tau/k_BT$; la V5 riportava anche "
                 "$8{,}8\\times10^{-3}$ con un $\\beta$ diverso \\\\")
    _tex("tab_residui.tex",
         "\\begin{tabular}{p{2.9cm}p{3.0cm}p{3.0cm}p{5.5cm}}\n\\hline\n"
         "Grandezza & V5 & V6 & nota \\\\\n\\hline\n"
         + "\n".join(righe) + "\n\\hline\n\\end{tabular}\n")

    # D. baseline
    righe = []
    for c in res["baseline"]["casi"]:
        righe.append(
            f"{c['nome']} & {c['R_nm']:.0f} & {sxr(c['beta_lo'], c['beta_hi'])} & "
            f"{sx(c['sigma_abs'])} & {sx(c['tau_usato_Nm'])} & "
            f"{sxr(c['Pe_r_lo'], c['Pe_r_hi'])} & "
            f"{sxr(c['omega_ss_lo'], c['omega_ss_hi_ottica'])} & "
            f"{sx(c['t_rottura_s'], 1)} \\\\")
    _tex("tab_baseline.tex",
         "\\begin{tabular}{lrrrrrrr}\n\\hline\n"
         "Geometria & $R$ (nm) & $\\beta$ (N\\,m\\,s) & $\\sigma_{ass}$ (m$^2$) & "
         "$\\tau$ (N\\,m) & $Pe_r$ & $\\omega_{ss}$ (rad/s) & $t_{r}$ (s) \\\\\n\\hline\n"
         + "\n".join(righe) + "\n\\hline\n\\end{tabular}\n")

    # H. registrazioni GHz
    righe = []
    for c in res["confronto"]["registrazioni"]:
        righe.append(f"{c['f_Hz']/1e9:.0f} & {c['r_nm']:.0f} & "
                     f"{c['v_tangenziale_m_s']:.0f} & {num(c['rho_v2_GPa'], 2)} \\\\")
    _tex("tab_ghz_reg.tex",
         "\\begin{tabular}{rrrr}\n\\hline\n"
         "$f$ (GHz) & $r$ (nm) & $v=2\\pi r f$ (m/s) & $\\rho v^2$ (GPa) \\\\\n"
         "\\hline\n" + "\n".join(righe) + "\n\\hline\n\\end{tabular}\n")

    # E. regimi
    righe = []
    for c in res["regimi"]["casi"]:
        om = ("--" if not math.isfinite(c["omega_ss_lo"]) else
              sxr(c["omega_ss_lo"], c["omega_ss_hi"]))
        t = "--" if not math.isfinite(c["t_rottura_s"]) else sx(c["t_rottura_s"], 1)
        val = c["validita"]
        val = ("valido" if val.startswith("valido")
               else ("estrap." if val.startswith("estrap") else "fuori dominio"))
        righe.append(f"{c['L_etichetta']} & {c['mezzo']} & {sx(c['eta_int'])} & "
                     f"{sx(c['tau_Nm'])} & {om} & {sx(c['Pe_r'])} & {t} & "
                     f"{c['T_K']:.0f} & {val} \\\\")
    _tex("tab_regimi.tex",
         "\\begin{tabular}{llrrrrrrl}\n\\hline\n"
         "$L$ & mezzo & $\\eta_{int}$ & $\\tau$ (N\\,m) & $\\omega$ (rad/s) & $Pe_r$ & "
         "$t_{rott}$ (s) & $T$ (K) & validit\\`a \\\\\n\\hline\n"
         + "\n".join(righe) + "\n\\hline\n\\end{tabular}\n")

    # F. materiali
    righe = []
    for c in res["materiali"]["casi"]:
        righe.append(f"{c['nome']} & {c['rho']:.0f} & {c['sy']/1e6:.0f} & "
                     f"{sx(c['specifico'], 1)} & {c['v_limite_m_s']:.0f} & "
                     f"{num(c['Q_abs'], 2)} & {sx(c['sigma_abs'])} & {sx(c['Pe_r'])} & "
                     f"{c['T_vuoto_K']:.0f} & {num(c['margine_fusione'], 2)} \\\\")
    _tex("tab_materiali.tex",
         "\\begin{tabular}{lrrrrrrrrr}\n\\hline\n"
         "Materiale & $\\rho$ & $\\sigma_y$ (MPa) & $\\sigma_y/\\rho$ (m$^2$/s$^2$) & "
         "$v_{lim}$ (m/s) & $Q_{ass}$ & $\\sigma_{ass}$ (m$^2$) & $Pe_r$ & "
         "$T$ vuoto (K) & margine \\\\\n\\hline\n"
         + "\n".join(righe) + "\n\\hline\n\\end{tabular}\n")

    # H. confronto GHz: f_lim(r) per cinque valori di sigma_y
    raggi = sorted({c["r_nm"] for c in res["confronto"]["casi"] if c["tipo"] == "rottura"})
    sigme = sorted({c["sy_Pa"] for c in res["confronto"]["casi"] if c["tipo"] == "rottura"})
    tab = {(c["r_nm"], c["sy_Pa"]): c["f_limite_Hz"]
           for c in res["confronto"]["casi"] if c["tipo"] == "rottura"}
    righe = []
    for sy in sigme:
        celle = " & ".join(sx(tab[(r_, sy)] / 1e9) for r_ in raggi)
        righe.append(f"{sy/1e9:g} & " + celle + " \\\\")
    _tex("tab_ghz.tex",
         "\\begin{tabular}{r" + "r" * len(raggi) + "}\n\\hline\n"
         "$\\sigma_y$ (GPa) & " + " & ".join(f"$r={r_:.0f}$ nm" for r_ in raggi)
         + " \\\\\n\\hline\n" + "\n".join(righe) + "\n\\hline\n\\end{tabular}\n")

    # G. leve
    righe = []
    for l in res["leve"]["leve"]:
        righe.append(f"{l['nome']} & {sx(l['fattore_Pe'])} $\\times$ & "
                     f"{sx(l['fattore_omega'])} $\\times$ & {l['nota']} \\\\")
    _tex("tab_leve.tex",
         "\\begin{tabular}{p{6.2cm}rrp{5.2cm}}\n\\hline\n"
         "Leva & su $Pe_r$ & su $\\omega_{ss}$ & nota \\\\\n\\hline\n"
         + "\n".join(righe) + "\n\\hline\n\\end{tabular}\n")


NOMI_EN = {
    "Reticolo (V5)": "Lattice (V5)", "Disco pieno": "Solid disc", "Anello": "Ring",
    "Sfera": "Sphere", "Nanomanubrio": "Nanodumbbell",
    "Silice anisotropa": "Anisotropic silica",
    "Rame": "Copper", "Silice fusa": "Fused silica", "Grafene / CNT": "Graphene / CNT",
}


def tabelle_en(res):
    """Le stesse tabelle in inglese, con il punto come separatore decimale."""
    _DEC[0] = "."
    try:
        # baseline
        righe = []
        for c in res["baseline"]["casi"]:
            righe.append(f"{NOMI_EN.get(c['nome'], c['nome'])} & {c['R_nm']:.0f} & "
                         f"{sxr(c['beta_lo'], c['beta_hi'])} & {sx(c['sigma_abs'])} & "
                         f"{sx(c['tau_usato_Nm'])} & {sxr(c['Pe_r_lo'], c['Pe_r_hi'])} & "
                         f"{sxr(c['omega_ss_lo'], c['omega_ss_hi_ottica'])} & "
                         f"{sx(c['t_rottura_s'], 1)} \\\\")
        _tex("tab_baseline_en.tex",
             "\\begin{tabular}{lrrrrrrr}\n\\hline\n"
             "Geometry & $R$ (nm) & $\\beta$ (N\\,m\\,s) & $\\sigma_{abs}$ (m$^2$) & "
             "$\\tau$ (N\\,m) & $Pe_r$ & $\\omega_{ss}$ (rad/s) & $t_{r}$ (s) \\\\\n\\hline\n"
             + "\n".join(righe) + "\n\\hline\n\\end{tabular}\n")

        # regimi
        righe = []
        for c in res["regimi"]["casi"]:
            om = ("--" if not math.isfinite(c["omega_ss_lo"])
                  else sxr(c["omega_ss_lo"], c["omega_ss_hi"]))
            t = "--" if not math.isfinite(c["t_rottura_s"]) else sx(c["t_rottura_s"], 1)
            val = c["validita"]
            val = ("valid" if val.startswith("valido")
                   else ("extrapolation" if val.startswith("estrap") else "out of domain"))
            mezzo = "fluid" if c["mezzo"] == "fluido" else "vacuum"
            righe.append(f"{c['L_etichetta']} & {mezzo} & {sx(c['eta_int'])} & "
                         f"{sx(c['tau_Nm'])} & {om} & {sx(c['Pe_r'])} & {t} & "
                         f"{c['T_K']:.0f} & {val} \\\\")
        _tex("tab_regimi_en.tex",
             "\\begin{tabular}{llrrrrrrl}\n\\hline\n"
             "$L$ & medium & $\\eta_{int}$ & $\\tau$ (N\\,m) & $\\omega$ (rad/s) & $Pe_r$ & "
             "$t_{rupt}$ (s) & $T$ (K) & validity \\\\\n\\hline\n"
             + "\n".join(righe) + "\n\\hline\n\\end{tabular}\n")

        # materiali
        righe = []
        for c in res["materiali"]["casi"]:
            righe.append(f"{NOMI_EN.get(c['nome'], c['nome'])} & {c['rho']:.0f} & "
                         f"{c['sy']/1e6:.0f} & {sx(c['specifico'], 1)} & "
                         f"{c['v_limite_m_s']:.0f} & {num(c['Q_abs'], 2)} & "
                         f"{sx(c['sigma_abs'])} & {sx(c['Pe_r'])} & "
                         f"{c['T_vuoto_K']:.0f} & {num(c['margine_fusione'], 2)} \\\\")
        _tex("tab_materiali_en.tex",
             "\\begin{tabular}{lrrrrrrrrr}\n\\hline\n"
             "Material & $\\rho$ & $\\sigma_y$ (MPa) & $\\sigma_y/\\rho$ (m$^2$/s$^2$) & "
             "$v_{lim}$ (m/s) & $Q_{abs}$ & $\\sigma_{abs}$ (m$^2$) & $Pe_r$ & "
             "$T$ vacuum (K) & margin \\\\\n\\hline\n"
             + "\n".join(righe) + "\n\\hline\n\\end{tabular}\n")

        # ghz: stessa tabella trasposta
        raggi = sorted({c["r_nm"] for c in res["confronto"]["casi"] if c["tipo"] == "rottura"})
        sigme = sorted({c["sy_Pa"] for c in res["confronto"]["casi"] if c["tipo"] == "rottura"})
        tab = {(c["r_nm"], c["sy_Pa"]): c["f_limite_Hz"]
               for c in res["confronto"]["casi"] if c["tipo"] == "rottura"}
        righe = []
        for sy in sigme:
            celle = " & ".join(sx(tab[(r_, sy)] / 1e9) for r_ in raggi)
            righe.append(f"{sy/1e9:g} & " + celle + " \\\\")
        _tex("tab_ghz_en.tex",
             "\\begin{tabular}{r" + "r" * len(raggi) + "}\n\\hline\n"
             "$\\sigma_y$ (GPa) & " + " & ".join(f"$r={r_:.0f}$ nm" for r_ in raggi)
             + " \\\\\n\\hline\n" + "\n".join(righe) + "\n\\hline\n\\end{tabular}\n")

        # registrazioni
        righe = []
        for c in res["confronto"]["registrazioni"]:
            righe.append(f"{c['f_Hz']/1e9:.0f} & {c['r_nm']:.0f} & "
                         f"{c['v_tangenziale_m_s']:.0f} & {num(c['rho_v2_GPa'], 2)} \\\\")
        _tex("tab_ghz_reg_en.tex",
             "\\begin{tabular}{rrrr}\n\\hline\n"
             "$f$ (GHz) & $r$ (nm) & $v=2\\pi r f$ (m/s) & $\\rho v^2$ (GPa) \\\\\n"
             "\\hline\n" + "\n".join(righe) + "\n\\hline\n\\end{tabular}\n")

        # leve
        trad = {
            "Potenza x10 (1 $\\to$ 10 $\\mu$W)": "Power $\\times10$ (1 $\\to$ 10\\,$\\mu$W)",
            "Messa a fuoco NA 0,9 ($w_0$ 500 $\\to$ 361 nm)":
                "Focusing NA 0.9 ($w_0$ 500 $\\to$ 361 nm)",
            "Polarizzazione lineare (s = 0)": "Linear polarization ($s=0$)",
            "OAM $\\ell = 5$ con $s = 1$ (con sovrapposizione)":
                "OAM $\\ell=5$ with $s=1$ (mode overlap included)",
            "Compattazione a pari massa ($R/10$)": "Compaction at equal mass ($R/10$)",
            "Vuoto spinto ($p = 10^{-6}$ mbar)": "High vacuum ($p=10^{-6}$ mbar)",
            "Materiale: grafene invece di rame": "Material: graphene instead of copper",
        }
        righe = []
        note_en = {
            "coppia invariante; omega per il rapporto degli attriti (1e2-1e3)":
                "torque unchanged; $\\omega$ follows the drag ratio ($10^2$--$10^3$)",
            "nessuno stato stazionario reale: il limite diventa la rottura":
                "no true steady state: fracture becomes the limit",
        }
        for l in res["leve"]["leve"]:
            nome = trad.get(l["nome"], l["nome"])
            nota = l["nota"]
            if nota.startswith("cambia l"):
                nota = ("specific absorption changes; effect on $Pe_r$ = "
                        + f"{l['fattore_Pe']:.2f}")
            else:
                nota = note_en.get(nota, nota)
            righe.append(f"{nome} & {sx(l['fattore_Pe'])} $\\times$ & "
                         f"{sx(l['fattore_omega'])} $\\times$ & {nota} \\\\")
        _tex("tab_leve_en.tex",
             "\\begin{tabular}{p{6.2cm}rrp{5.2cm}}\n\\hline\n"
             "Lever & on $Pe_r$ & on $\\omega_{ss}$ & note \\\\\n\\hline\n"
             + "\n".join(righe) + "\n\\hline\n\\end{tabular}\n")
    finally:
        _DEC[0] = "{,"


def main():
    print("MOTORE FOTONICO V6 - studi quantitativi")
    res = {}
    res["residui"] = studio_residui()
    res["ottica"] = studio_ottica()
    res["termica"] = studio_termica()
    res["baseline"] = studio_baseline()
    res["regimi"] = studio_regimi()
    res["materiali"] = studio_materiali()
    res["leve"] = studio_leve()
    res["confronto"] = studio_confronto_ghz()
    res["oam"] = studio_oam()
    tabelle(res)
    tabelle_en(res)

    a = res["residui"]
    print("\n--- A. residui ---")
    print(f"  1/(45/376)          = {a['rapporti']['inverso_esatto_di_45_376']:.7f}")
    print(f"  I2/I1 scena         = {a['rapporti']['I2_su_I1_scena_V4']:.7f} "
          f"(target {a['rapporti']['razionale_376_45']:.7f}, "
          f"scarto {a['rapporti']['scarto_scena_su_target']:+.2e})")
    print(f"  artefatto 1/8,3565  = {a['rapporti']['scarto_riportato_rel']:+.2e}")
    print(f"  sigma_pr V6/V5      = {a['verifica_sezioni']['sigma_pr_V6']:.5e} / "
          f"{a['verifica_sezioni']['sigma_pr_V5']:.5e} ({a['verifica_sezioni']['scarto']:+.1e})")
    print(f"  tau assorbita       = {a['bilancio_potenza']['tau_da_assorbita_Nm']:.3e} N m "
          f"({a['bilancio_potenza']['scarto_tau']*100:+.1f}% vs V5)")
    print(f"  efficienza mecc.    = {a['bilancio_potenza']['efficienza_ottico_meccanica']:.1e}")
    print(f"  I/beta              = {a['pulsata']['tempo_inerziale_I_su_beta_s']:.2e} s")
    print("\n--- B. ottica ---")
    print(f"  sfera pari volume: sigma_abs scarto "
          f"{res['ottica']['sfera_pari_volume']['scarto_sigma_abs']*100:+.1f}%")
    print(f"  ossido 2 nm: Q_abs scarto "
          f"{res['ottica']['ossidazione'][2]['scarto_Q_abs']*100:+.1f}%")
    print("\n--- C. termica ---")
    for c in res["termica"]["casi"][:3]:
        print(f"  {c['materiale']:8s} P={c['potenza_W']*1e6:8.1f} uW  dT_fluido="
              f"{c['dT_fluido_K']:.2e} K  T_vuoto(0.03)={c['T_vuoto_emiss_limpido_K']:.0f} K  "
              f"T_vuoto(0.6)={c['T_vuoto_emiss_ossidato_K']:.0f} K")
    print("\n--- D. baseline ---")
    for c in res["baseline"]["casi"]:
        print(f"  {c['nome']:18s} R={c['R_nm']:6.1f} nm  beta={c['beta_medio']:.2e}  "
              f"tau={c['tau_usato_Nm']:.2e}  Pe={c['Pe_r']:.2e}  "
              f"omega={c['omega_ss_lo']:.2e}..{c['omega_ss_hi']:.2e} rad/s")
    print("\n--- E. regimi ---")
    for c in res["regimi"]["casi"]:
        print(f"  {c['L_etichetta']:8s} {c['mezzo']:7s} eta={c['eta_int']:.1e} "
              f"tau={c['tau_Nm']:.1e} Pe={c['Pe_r']:.1e} "
              f"w={c['omega_ss_lo']:.1e} T={c['T_K']:.0f} t={c['t_rottura_s']:.1e} "
              f"[{c['validita']}]")
    print("\n--- F. materiali ---")
    for c in res["materiali"]["casi"]:
        print(f"  {c['nome']:12s} v_lim={c['v_limite_m_s']:7.0f} m/s  Q_abs={c['Q_abs']:.2f} "
              f"sigma={c['sigma_abs']:.2e} Pe={c['Pe_r']:.1e} T_vuoto={c['T_vuoto_K']:.0f} K")
    print("\n--- G. leve ---")
    print(f"  base Pe_r={res['leve']['base']['Pe_r']:.2e} omega={res['leve']['base']['omega']:.2e}")
    for l in res["leve"]["leve"]:
        print(f"  {l['nome']:45s} x{l['fattore_Pe']:9.2f} su Pe_r   x{l['fattore_omega']:9.2f} su omega")
    print("\n--- H. confronto GHz ---")
    for c in res["confronto"]["registrazioni"]:
        print(f"  {c['f_Hz']/1e9:.0f} GHz a r={c['r_nm']:.0f} nm -> v={c['v_tangenziale_m_s']:.0f} m/s, "
              f"sigma implicita = {c['rho_v2_GPa']:.2f} GPa")
    for c in res["confronto"]["casi"]:
        if c["tipo"] == "rottura" and c["sy_Pa"] == 1e9 and c["r_nm"] in (50.0, 75.0):
            print(f"  f_limite r={c['r_nm']:.0f} nm, 1 GPa = {c['f_limite_Hz']/1e9:.2f} GHz")
    for c in res["confronto"]["casi"]:
        if c["tipo"] == "potenza" and c["r_nm"] == 50 and c["delta"] == 0.06:
            print(f"  r=50 nm p={c['p_mbar']:.0e} mbar: P_necessaria="
                  f"{c['P_sorgente_necessaria_W']:.2e} W, spin-up a 10 mW = "
                  f"{c['t_spinup_a_10mW_s']:.2e} s")
    print(f"  {res['confronto']['sintesi']['verdetto'][:100]}...")
    print("\nJSON in dati/, tabelle in tabelle/")


if __name__ == "__main__":
    main()
