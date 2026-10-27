"""
MOTORE FOTONICO V5.1 - VUOTO PROFONDO
Autore: Dr. Francesco Bulla

Studio del comportamento del motore Spin/Contro-spin nello spazio profondo:

    g = 0        nessun peso, nessuna sedimentazione, nessuna convezione
    eta = 0      nessun fluido, nessun attrito viscoso
    beta = 0     nessuna coppia dissipativa: tau_diss = -beta*omega = 0

Conseguenze dirette, e sono il punto di tutto il modulo:

  * non esiste una velocita' di regime omega_ss = tau/beta: il rotore non
    puo' raggiungere un equilibrio, accelera finche' non si rompe o finche'
    il laser non viene spento;
  * il tempo di vita a piena potenza e' quindi il tempo critico di rottura
    strutturale t_crit = omega_rottura / alpha;
  * il materiale decide tutto: il rame ha sigma_max = 200 MPa e
    rho = 8960 kg/m^3, il composito Grafene/CNT ha sigma_max = 50 GPa e
    rho = 2200 kg/m^3. Il rapporto omega_rottura(Grafene)/omega_rottura(Cu)
    vale sqrt((50e9/2200)/(200e6/8960)) = 31.9.

Due scale, come richiesto:
    nano  : ingombro 240 nm  (R = 120 nm)
    macro : ingombro 1 m     (R = 0.5 m)

Uso:
    python v51_vuoto.py
"""

import json
import math
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
PADRE = os.path.dirname(BASE)
if PADRE not in sys.path:
    sys.path.insert(0, PADRE)

import v5_nucleo as n  # noqa: E402

DATI = os.path.join(BASE, "dati")
os.makedirs(DATI, exist_ok=True)

# ---------------------------------------------------------------------------
# 0. COSTANTI
# ---------------------------------------------------------------------------
C_LIGHT = 299792458.0
LAMBDA = 532e-9
OMEGA_OTTICA = 2.0 * math.pi * C_LIGHT / LAMBDA      # rad/s, 3.5416e15

G_TERRA = 9.80665

# Geometria di riferimento V4 (rame)
R_DISC_V4 = 3.275 * 100e-9          # m, raggio esterno del rotore V4
M_V4 = 1.2463630573655231e-18 + 4.487710753596799e-18   # kg, massa totale
I1_V4 = 3.092e-32                   # kg m^2
I2_V4 = 2.583e-31                   # kg m^2

# Materiali (valori dichiarati nella richiesta)
MATERIALI = {
    "Rame": {
        "rho_kg_m3": 8960.0,
        "sigma_max_Pa": 200e6,          # snervamento
        "E_Pa": 110e9,
        "nota": "riferimento V4: duttile, si snerva a 200 MPa",
    },
    "Grafene/CNT": {
        "rho_kg_m3": 2200.0,
        "sigma_max_Pa": 50e9,           # 50 GPa, valore conservativo
        "E_Pa": 1e12,                   # 1 TPa
        "nota": "composito Grafene/Nanotubi di Carbonio, valore conservativo",
    },
}

SCALE = {
    "nano_240nm": {"ingombro_m": 240e-9, "R_m": 120e-9},
    "macro_1m": {"ingombro_m": 1.0, "R_m": 0.5},
}

# Potenze ottiche di prova
POTENZE = {"1_uW": 1.0e-6, "1_mW": 1.0e-3}

# Modalita' pulsata richiesta
DT_ON = 1.0e-3        # s, laser acceso
DT_OFF = 1.0         # s, laser spento
T_MISSIONE = 3600.0   # s, un'ora


# ---------------------------------------------------------------------------
# 1. SCALATURA GEOMETRICA
# ---------------------------------------------------------------------------
def fattore_scala(R):
    """Fattore k rispetto alla geometria V4 (R_DISC = 327.5 nm)."""
    return R / R_DISC_V4


def massa_totale(R, materiale):
    """Massa scalata come k^3 e corretta per la densita' del materiale."""
    k = fattore_scala(R)
    m = MATERIALI[materiale]["rho_kg_m3"]
    return M_V4 * k ** 3 * (m / MATERIALI["Rame"]["rho_kg_m3"])


def inerzie(R, materiale):
    """Momenti d'inerzia dei due rotori: I scala come rho * k^5."""
    k = fattore_scala(R)
    f = k ** 5 * (MATERIALI[materiale]["rho_kg_m3"]
                  / MATERIALI["Rame"]["rho_kg_m3"])
    return I1_V4 * f, I2_V4 * f


# ---------------------------------------------------------------------------
# 2. ROTTURA STRUTTURALE
# ---------------------------------------------------------------------------
def omega_rottura(R, materiale):
    """Tensione centrifuga sigma = rho omega^2 R^2 posta uguale a sigma_max."""
    mat = MATERIALI[materiale]
    return math.sqrt(mat["sigma_max_Pa"] / (mat["rho_kg_m3"] * R ** 2))


def coppia_massima(P, momento_fotonico=1.0, frazione=1.0):
    """tau_max = frazione * P / omega_ottica.

    Ogni fotone porta un momento angolare hbar: il flusso di momento
    angolare di un fascio di potenza P vale P/omega. Con momento_fotonico
    = 2 si intende il caso a due lobi (coppia di ordine 2 della V5).
    La frazione e' la parte di potenza effettivamente intercettata dalla
    struttura (v5_nucleo.scattering_globale): senza di essa la coppia e'
    quella del fascio intero, che nessun rotore di questa scala riceve.
    """
    return momento_fotonico * frazione * P / OMEGA_OTTICA


def tempo_critico(I_tot, tau_max, omega_lim):
    return omega_lim * I_tot / tau_max


# ---------------------------------------------------------------------------
# 3. PROPULSIONE
# ---------------------------------------------------------------------------
def forze(P, R):
    """Tre stime della spinta lungo z, come richiesto."""
    k = fattore_scala(R)
    sc = n.scattering_globale(k)
    F_scat = float(sc["F_z_N"]) * (P / n.P_LASER)
    return {
        "scattering_N": F_scat,
        "P_su_c_N": P / C_LIGHT,
        "2P_su_c_N": 2.0 * P / C_LIGHT,
        "scattering_dettaglio": {
            "regime_filo": str(sc["regime_filo"]),
            "regime_nucleo": str(sc["regime_nucleo"]),
            "Q_pr_filo": float(sc["Q_pr_filo"]),
            "Q_pr_nucleo": float(sc["Q_pr_nucleo"]),
            "sigma_pr_m2": float(sc["sigma_pr_m2"]),
            "frazione_fascio": float(sc["frazione_fascio"]),
        },
    }


# ---------------------------------------------------------------------------
# 4. CASO COMPLETO: una scala, un materiale, una potenza
# ---------------------------------------------------------------------------
def caso(scala, materiale, potenza, P):
    R = SCALE[scala]["R_m"]
    mat = MATERIALI[materiale]
    I1, I2 = inerzie(R, materiale)
    I_tot = I1 + I2
    m = massa_totale(R, materiale)
    omega_lim = omega_rottura(R, materiale)
    sc = n.scattering_globale(fattore_scala(R))
    f_int = float(sc["frazione_fascio"])
    tau = coppia_massima(P, frazione=f_int)
    tau_pieno = coppia_massima(P)
    alpha = tau / I_tot
    t_crit = omega_lim / alpha

    f = forze(P, R)
    a_scat = f["scattering_N"] / m
    a_pc = f["P_su_c_N"] / m
    a_2pc = f["2P_su_c_N"] / m

    # Un'ora a piena potenza: il rotore si rompe se t_crit < 3600 s.
    omega_1h = alpha * T_MISSIONE
    rottura_in_un_ora = omega_1h > omega_lim

    # Energia rotazionale immagazzinata al punto di rottura
    E_rot = 0.5 * I_tot * omega_lim ** 2
    # Velocita' periferica al bordo: e' il vero limite, non omega
    v_bordo = omega_lim * R

    return {
        "scala": scala,
        "materiale": materiale,
        "ingombro_m": SCALE[scala]["ingombro_m"],
        "R_m": R,
        "k": fattore_scala(R),
        "rho_kg_m3": mat["rho_kg_m3"],
        "sigma_max_Pa": mat["sigma_max_Pa"],
        "E_Pa": mat["E_Pa"],
        "massa_totale_kg": m,
        "volume_totale_m3": m / mat["rho_kg_m3"],
        "I1_kgm2": I1,
        "I2_kgm2": I2,
        "I_tot_kgm2": I_tot,
        "Potenza_W": P,
        "omega_ottica_rad_s": OMEGA_OTTICA,
        "frazione_fascio": f_int,
        "tau_max_Nm": tau,
        "tau_fascio_intero_Nm": tau_pieno,
        "alpha_rad_s2": alpha,
        "omega_rottura_rad_s": omega_lim,
        "giri_al_secondo_rottura": omega_lim / (2.0 * math.pi),
        "v_periferica_rottura_m_s": v_bordo,
        "v_periferica_frazione_c": v_bordo / C_LIGHT,
        "t_critico_s": t_crit,
        "t_critico_ore": t_crit / 3600.0,
        "omega_dopo_1h_rad_s": omega_1h,
        "rottura_entro_1h": rottura_in_un_ora,
        "energia_rotazionale_rottura_J": E_rot,
        "spinta_N": f["scattering_N"],
        "spinta_P_su_c_N": f["P_su_c_N"],
        "spinta_2P_su_c_N": f["2P_su_c_N"],
        "accelerazione_scattering_m_s2": a_scat,
        "accelerazione_scattering_g": a_scat / G_TERRA,
        "accelerazione_P_su_c_m_s2": a_pc,
        "accelerazione_P_su_c_g": a_pc / G_TERRA,
        "accelerazione_2P_su_c_m_s2": a_2pc,
        "accelerazione_2P_su_c_g": a_2pc / G_TERRA,
        "scattering": f["scattering_dettaglio"],
    }


# ---------------------------------------------------------------------------
# 5. SPIN-UP NEL VUOTO (beta = 0)
# ---------------------------------------------------------------------------
def spin_up(caso_d, punti=240):
    """omega(t) = alpha t fino alla rottura: nessun regime stazionario."""
    alpha = caso_d["alpha_rad_s2"]
    t_end = caso_d["t_critico_s"]
    tempi = [t_end * i / (punti - 1) for i in range(punti)]
    return {
        "materiale": caso_d["materiale"],
        "scala": caso_d["scala"],
        "alpha_rad_s2": alpha,
        "omega_rottura_rad_s": caso_d["omega_rottura_rad_s"],
        "t_critico_s": t_end,
        "curva": [[t, alpha * t] for t in tempi],
        "nota": "in assenza di attrito non esiste omega_ss = tau/beta: "
                "la curva e' una retta che termina sul muro di rottura",
    }


# ---------------------------------------------------------------------------
# 6. MODALITA' PULSATA
# ---------------------------------------------------------------------------
def pulsata(caso_d, dt_on=DT_ON, dt_off=DT_OFF, t_missione=T_MISSIONE):
    """Accensione/spegnimento: omega sale di alpha*dt_on a ogni impulso e
    resta costante nella pausa (inerzia pura, beta = 0)."""
    alpha = caso_d["alpha_rad_s2"]
    omega_lim = caso_d["omega_rottura_rad_s"]
    a = caso_d["accelerazione_scattering_m_s2"]

    periodo = dt_on + dt_off
    cicli = int(t_missione // periodo)
    omega, v, t = 0.0, 0.0, 0.0
    curva = [[0.0, 0.0, 0.0]]
    t_superamento = None
    cicli_superamento = None

    for i in range(cicli):
        omega += alpha * dt_on
        v += a * dt_on
        t += dt_on
        curva.append([t, omega, v])
        if t_superamento is None and omega > omega_lim:
            t_superamento = t
            cicli_superamento = i + 1
        t += dt_off
        curva.append([t, omega, v])

    # Potenza che porta omega esattamente al muro a fine missione
    P_sicura = (omega_lim * caso_d["I_tot_kgm2"] * OMEGA_OTTICA
                / (cicli * dt_on * caso_d["frazione_fascio"]))

    return {
        "materiale": caso_d["materiale"],
        "scala": caso_d["scala"],
        "dt_on_s": dt_on,
        "dt_off_s": dt_off,
        "cicli": cicli,
        "tempo_on_totale_s": cicli * dt_on,
        "alpha_rad_s2": alpha,
        "omega_finale_rad_s": omega,
        "omega_rottura_rad_s": omega_lim,
        "supera_rottura": t_superamento is not None,
        "t_superamento_s": t_superamento,
        "cicli_superamento": cicli_superamento,
        "frazione_muro_raggiunta": omega / omega_lim,
        "delta_v_m_s": v,
        "delta_v_km_s": v / 1000.0,
        "delta_v_frazione_c": v / C_LIGHT,
        "potenza_sicura_W": P_sicura,
        "curva": curva,
        "nota": "nella pausa omega resta costante: nel vuoto non c'e' "
                "nulla che possa frenare il rotore",
    }


# ---------------------------------------------------------------------------
# 7. MODALITA' PULSATA SICURA (potenza ridotta: si arriva al muro solo alla fine)
# ---------------------------------------------------------------------------
def pulsata_sicura(caso_d, dt_on=DT_ON, dt_off=DT_OFF, t_missione=T_MISSIONE):
    """Stessa sequenza di impulsi, con la potenza scelta in modo che dopo
    un'ora intera omega resti a meta' del muro di rottura (margine 2)."""
    alpha = caso_d["alpha_rad_s2"]
    omega_lim = caso_d["omega_rottura_rad_s"]
    a = caso_d["accelerazione_scattering_m_s2"]

    periodo = dt_on + dt_off
    cicli = int(t_missione // periodo)
    fattore = 0.5 * omega_lim / (alpha * dt_on * cicli)   # P_sicura / P
    alpha_s = alpha * fattore
    a_s = a * fattore

    omega, v, t = 0.0, 0.0, 0.0
    curva = [[0.0, 0.0, 0.0]]
    for _ in range(cicli):
        omega += alpha_s * dt_on
        v += a_s * dt_on
        t += dt_on
        curva.append([t, omega, v])
        t += dt_off
        curva.append([t, omega, v])

    return {
        "materiale": caso_d["materiale"],
        "scala": caso_d["scala"],
        "dt_on_s": dt_on,
        "dt_off_s": dt_off,
        "cicli": cicli,
        "potenza_W": caso_d["Potenza_W"] * fattore,
        "alpha_rad_s2": alpha_s,
        "omega_finale_rad_s": omega,
        "omega_rottura_rad_s": omega_lim,
        "frazione_muro_raggiunta": omega / omega_lim,
        "delta_v_m_s": v,
        "delta_v_km_s": v / 1000.0,
        "curva": curva,
        "nota": "potenza pari a %.3f volte 1 uW: la scala a gradini si ferma "
                "a meta' del muro dopo i %d impulsi dell'ora"
                % (fattore, cicli),
    }


# ---------------------------------------------------------------------------
# 8. CONFRONTO RAME vs GRAFENE
# ---------------------------------------------------------------------------
def confronto(scala, P):
    tabella = {}
    for materiale in ("Rame", "Grafene/CNT"):
        d = caso(scala, materiale, "prova", P)
        tabella[materiale] = {
            "massa_totale_kg": d["massa_totale_kg"],
            "I_tot_kgm2": d["I_tot_kgm2"],
            "omega_rottura_rad_s": d["omega_rottura_rad_s"],
            "v_periferica_rottura_m_s": d["v_periferica_rottura_m_s"],
            "alpha_rad_s2": d["alpha_rad_s2"],
            "t_critico_s": d["t_critico_s"],
            "t_critico_ore": d["t_critico_ore"],
            "energia_rotazionale_J": d["energia_rotazionale_rottura_J"],
            "accelerazione_scattering_g": d["accelerazione_scattering_g"],
        }
    r = tabella["Grafene/CNT"]["omega_rottura_rad_s"] / tabella["Rame"]["omega_rottura_rad_s"]
    tabella["rapporto_omega_rottura"] = r
    tabella["nota"] = ("sqrt((sigma_graf/rho_graf)/(sigma_cu/rho_cu)) = %.4f; "
                       "il Grafene vince di oltre un ordine di grandezza" % r)
    return tabella


# ---------------------------------------------------------------------------
# 8. USCITA
# ---------------------------------------------------------------------------
def scrivi(nome, contenuto):
    percorso = os.path.join(DATI, nome)
    with open(percorso, "w", encoding="utf-8") as fh:
        json.dump(contenuto, fh, indent=2, ensure_ascii=False, allow_nan=False)
    print("  %-34s %7.1f kB" % (nome, os.path.getsize(percorso) / 1024.0))


def main():
    print("[V5.1] vuoto profondo: g = 0, eta = 0, beta = 0")

    # --- V1: parametri dei materiali ---------------------------------------
    materiali = {
        "riferimento": "Rame V4 vs composito Grafene/CNT",
        "Rame": MATERIALI["Rame"],
        "Grafene/CNT": MATERIALI["Grafene/CNT"],
        "rapporto_densita": (MATERIALI["Grafene/CNT"]["rho_kg_m3"]
                             / MATERIALI["Rame"]["rho_kg_m3"]),
        "rapporto_resistenza": (MATERIALI["Grafene/CNT"]["sigma_max_Pa"]
                                / MATERIALI["Rame"]["sigma_max_Pa"]),
        "guadagno_omega": math.sqrt(
            (MATERIALI["Grafene/CNT"]["sigma_max_Pa"]
             / MATERIALI["Grafene/CNT"]["rho_kg_m3"])
            / (MATERIALI["Rame"]["sigma_max_Pa"] / MATERIALI["Rame"]["rho_kg_m3"])),
        "ambiente": {"g_m_s2": 0.0, "eta_Pa_s": 0.0, "beta_Nms": 0.0,
                     "conseguenza": "omega_ss = tau/beta non esiste: "
                                    "nessuna velocita' di regime"},
    }
    scrivi("V1_parametri_materiali.json", materiali)

    # --- V2: spin-up nel vuoto, tutte le combinazioni -----------------------
    spin = {"descrizione": "omega(t) = alpha t fino alla rottura, beta = 0",
            "casi": {}}
    for scala in SCALE:
        for materiale in ("Rame", "Grafene/CNT"):
            for nome_P, P in POTENZE.items():
                chiave = "%s|%s|%s" % (scala, materiale, nome_P)
                d = caso(scala, materiale, nome_P, P)
                s = spin_up(d)
                s["dettaglio"] = {kk: vv for kk, vv in d.items()
                                  if kk not in ("scattering",)}
                spin["casi"][chiave] = s
    scrivi("V2_spin_up.json", spin)

    # --- V3: modalita' pulsata ---------------------------------------------
    puls = {"descrizione": "1 ms acceso, 1 s spento, un'ora di missione; per "
                           "ogni caso anche la variante a potenza scelta per "
                           "fermarsi a meta' del muro di rottura (margine 2)",
            "casi": {}}
    for scala in SCALE:
        for materiale in ("Rame", "Grafene/CNT"):
            chiave = "%s|%s" % (scala, materiale)
            d = caso(scala, materiale, "1_uW", POTENZE["1_uW"])
            p = pulsata(d)
            p["sicura"] = pulsata_sicura(d)
            puls["casi"][chiave] = p
    scrivi("V3_pulsata.json", puls)

    # --- V4: propulsione ----------------------------------------------------
    prop = {"descrizione": "spinta e accelerazione lineare nel vuoto",
            "casi": {}}
    for scala in SCALE:
        for materiale in ("Rame", "Grafene/CNT"):
            for nome_P, P in POTENZE.items():
                chiave = "%s|%s|%s" % (scala, materiale, nome_P)
                d = caso(scala, materiale, nome_P, P)
                prop["casi"][chiave] = {
                    "spinta_scattering_N": d["spinta_N"],
                    "spinta_P_su_c_N": d["spinta_P_su_c_N"],
                    "spinta_2P_su_c_N": d["spinta_2P_su_c_N"],
                    "massa_kg": d["massa_totale_kg"],
                    "a_scattering_m_s2": d["accelerazione_scattering_m_s2"],
                    "a_scattering_g": d["accelerazione_scattering_g"],
                    "a_P_su_c_m_s2": d["accelerazione_P_su_c_m_s2"],
                    "a_P_su_c_g": d["accelerazione_P_su_c_g"],
                    "a_2P_su_c_m_s2": d["accelerazione_2P_su_c_m_s2"],
                    "a_2P_su_c_g": d["accelerazione_2P_su_c_g"],
                    "sostenibile": bool(d["accelerazione_scattering_g"] < 10.0),
                    "delta_v_1h_scattering_m_s":
                        d["accelerazione_scattering_m_s2"] * T_MISSIONE,
                    "delta_v_1h_scattering_km_s":
                        d["accelerazione_scattering_m_s2"] * T_MISSIONE / 1000.0,
                    "delta_v_1h_scattering_frazione_c":
                        d["accelerazione_scattering_m_s2"] * T_MISSIONE / C_LIGHT,
                }
    scrivi("V4_propulsione.json", prop)

    # --- V5: confronto e sintesi -------------------------------------------
    sintesi = {"confronto": {}, "sintesi": {}}
    for scala in SCALE:
        sintesi["confronto"][scala] = confronto(scala, POTENZE["1_uW"])
    for scala in SCALE:
        for materiale in ("Rame", "Grafene/CNT"):
            d = caso(scala, materiale, "1_uW", POTENZE["1_uW"])
            p = pulsata(d)
            sintesi["sintesi"]["%s|%s" % (scala, materiale)] = {
                "omega_rottura_rad_s": d["omega_rottura_rad_s"],
                "t_critico_s": d["t_critico_s"],
                "t_critico_ore": d["t_critico_ore"],
                "rottura_entro_1h": d["rottura_entro_1h"],
                "pulsata_supera": p["supera_rottura"],
                "pulsata_t_superamento_s": p["t_superamento_s"],
                "pulsata_delta_v_m_s": p["delta_v_m_s"],
                "potenza_sicura_W": p["potenza_sicura_W"],
                "accelerazione_g": d["accelerazione_scattering_g"],
            }
    scrivi("V5_sintesi.json", sintesi)

    # --- stdout leggibile ---------------------------------------------------
    print("\n  %-26s %10s %12s %14s %12s" %
          ("caso", "omega_lim", "t_crit [s]", "a [g]", "dv/1h [m/s]"))
    for scala in SCALE:
        for materiale in ("Rame", "Grafene/CNT"):
            d = caso(scala, materiale, "1_uW", POTENZE["1_uW"])
            p = pulsata(d)
            print("  %-26s %10.3e %12.3e %14.3e %12.3e"
                  % ("%s|%s" % (scala, materiale), d["omega_rottura_rad_s"],
                     d["t_critico_s"], d["accelerazione_scattering_g"],
                     p["delta_v_m_s"]))


if __name__ == "__main__":
    main()
