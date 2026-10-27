"""
MOTORE FOTONICO V5 - MODULO SPAZIALE (MICROGRAVITA', g = 0)
Autore: Dr. Francesco Bulla / Sistema GNRX

Estende v5_nucleo.py al caso di assenza di gravita' e confronta, per lo Spin e
per il Contro-spin, tre configurazioni:

  S1   scala 10^-9 m          (k = 10^-3 rispetto alla geometria V4)
  S2   dispositivo V4         massa 1.246e-18 kg (Spin) e 4.488e-18 kg
                              (Contro-spin): e' il "10^-18" del preprint V4
  S3   scala 1 m              (k = 10^6), con tre sistemi di materiali
                              1. grafene / nanotubi di carbonio
                              2. carburo di silicio (variante: diamante
                                 nanocristallino)
                              3. metamateriale dielettrico silicio-su-ossido

Per ogni caso e per ogni rotore si calcolano:

  - attrito rotazionale beta sulla geometria reale (telaio di fili);
  - coppia ottica tau = P_ass/omega_ott = F*lambda/(2 pi), l'unico modello di
    momento angolare ceduto usato nel documento, sia con la forza dichiarata
    nella V4 sia con quella ricavata dallo scattering (v5_nucleo);
  - regime libero omega = tau / beta e regime accoppiato che rispetta
    omega2/omega1 = -I1/I2 (45/376);
  - diffusione rotazionale D_r e numero di Peclet rotazionale Pe_r;
  - confronto terra / spazio: peso, spinta/peso, velocita' di sedimentazione,
    numero di Peclet gravitazionale, numero di Grashof;
  - margine meccanico: tensione centrifuga rho*omega^2*R^2 contro snervamento.

Nessun numero e' scritto a mano: tutto deriva da v5_nucleo.py piu' le
proprieta' dei materiali dichiarate qui sotto.

Uso:
    python v5_spazio.py
"""

import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v5_nucleo as n  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "dati")
os.makedirs(OUT, exist_ok=True)

G_TERRA = 9.80665
RHO_FLUIDO = n.RHO

# ---------------------------------------------------------------------------
# 1. GEOMETRIA DEI DUE ROTORI (vertici della scena V4, convertiti in SI)
#    1 unita' di scena = 100 nm = 1e-7 m
# ---------------------------------------------------------------------------
UA = 1e-7

ROTORI = {
    "Spin": {
        "volume_m3": 0.1391 * UA ** 3,     # 1.391e-22 m^3
        "r2_m2": 2.4807 * UA ** 2,         # 2.4807e-14 m^2, <x^2 + y^2>
        "raggio_m": 3.70 * UA,             # 3.70e-7 m, semiasse maggiore
        "lunghezza_filo_m": 43.240 * UA,
        "massa_v4_kg": 1.2463630573655231e-18,
    },
    "Contro-spin": {
        "volume_m3": 0.5009 * UA ** 3,
        "r2_m2": 5.7567 * UA ** 2,
        "raggio_m": 3.275 * UA,
        "lunghezza_filo_m": 42.986 * UA,
        "massa_v4_kg": 4.487710753596799e-18,
    },
}

RHO_CU = 8960.0
SIGMA_CU = 2.0e8          # Pa, snervamento del rame ricotto

# ---------------------------------------------------------------------------
# 2. MATERIALI PER IL CASO A 1 METRO
#    rho  densita' (kg/m^3)
#    sig  resistenza a trazione (Pa)
#    E    modulo elastico (Pa)
#    n,k  indice di rifrazione complesso a 532 nm
#    cp   calore specifico (J/kg/K)
# ---------------------------------------------------------------------------
MATERIALI = {
    "Grafene/CNT": {
        "rho": 1350.0, "sig": 5.0e9, "E": 200.0e9,
        "n": 2.70, "k": 1.40, "cp": 700.0,
        "nota": "reticolo di nanotubi in fibra; densita' efficace 1.35 g/cm3",
    },
    "SiC": {
        "rho": 3210.0, "sig": 3.0e9, "E": 410.0e9,
        "n": 2.65, "k": 0.02, "cp": 750.0,
        "nota": "SiC CVD; variante diamante nanocristallino: rho 3510, "
                "sig 1.5e9, E 1.05e12, n 2.42",
    },
    "Metamateriale Si/SiO2": {
        "rho": 2250.0, "sig": 0.5e9, "E": 120.0e9,
        "n": 2.40, "k": 0.03, "cp": 700.0,
        "nota": "multistrato silicio su ossido, mezzo efficace Maxwell-Garnett",
    },
}

# ---------------------------------------------------------------------------
# 3. MATTONI FISICI
# ---------------------------------------------------------------------------
def attrito_rotazionale(raggio_m, eta=n.ETA):
    """beta = (32/3) eta R^3 sul telaio di fili (frazione bagnata FILL_FRAME).

    Identico a v5_nucleo.beta_telaio: un'unica definizione di attrito per
    tutto il progetto.
    """
    return n.beta_telaio(raggio_m, eta)


def raggio_equivalente(volume_m3):
    return (3.0 * volume_m3 / (4.0 * math.pi)) ** (1.0 / 3.0)


def forza_dichiarata(k):
    """Forza di progetto della V4, saturata quando il dispositivo copre il fascio."""
    return (n.P_LASER / n.C_LIGHT) * (n.ETA_S_DECL + n.ETA_C_DECL * math.cos(n.ALPHA_DECL)) \
        * min(1.0, k ** 2)


def forza_scattering(k):
    """Forza da calcolo di scattering (Mie + Rayleigh cilindrico), v5_nucleo."""
    return n.scattering_globale(k)["F_z_N"]


def peso_netto(massa, volume, g, rho_f=RHO_FLUIDO):
    """Peso apparente in fluido: (m - rho_f V) g.  In assenza di gravita' e' 0."""
    return (massa - rho_f * volume) * g


# ---------------------------------------------------------------------------
# 4. CASO DI UN ROTORE A UNA DATA SCALA E CON UN DATO MATERIALE
# ---------------------------------------------------------------------------
def caso_rotore(nome, k, materiale, eta=n.ETA):
    geo = ROTORI[nome]
    rho = materiale["rho"]
    sig = materiale["sig"]

    volume = geo["volume_m3"] * k ** 3
    r2 = geo["r2_m2"] * k ** 2
    raggio = geo["raggio_m"] * k
    r_h = math.sqrt(r2)                      # raggio idrodinamico <r^2>^(1/2), convenzione V4
    massa = rho * volume
    inerzia = massa * r2
    r_eq = raggio_equivalente(volume)        # sfera compatta di pari volume (riferimento)

    beta = attrito_rotazionale(raggio, eta)
    # Einstein: la stessa beta di tau_diss = -beta*omega governa la diffusione.
    # Usare kT/(8 pi eta R^3) insieme a una beta ridotta da FILL_FRAME gonfia
    # Pe_r di 1/FILL_FRAME = 15.5: qui i due lati sono coerenti per costruzione.
    D_r = n.KT / beta

    F_dich = forza_dichiarata(k)
    F_mie = forza_scattering(k)
    # Unico modello di momento angolare ceduto: tau = F lambda/(2 pi).
    # La coppia non contiene il raggio del rotore: dipende solo dalla potenza
    # effettivamente intercettata, cioe' dalla forza di scattering.
    tau_dich = n.coppia_oam(F_dich)
    tau_mie = n.coppia_oam(F_mie)

    omega_dich = tau_dich / beta
    omega_mie = tau_mie / beta

    # margine meccanico alla velocita' di regime
    def tensione(w):
        return rho * (w * raggio) ** 2

    def omega_limite():
        return math.sqrt(sig / (rho * raggio ** 2))

    w_lim = omega_limite()

    # trasporto traslazionale: Stokes sul raggio idrodinamico
    v_sed = abs(peso_netto(massa, volume, G_TERRA)) / (6 * math.pi * eta * r_h)
    Pe_grav = massa * G_TERRA * r_h / n.KT

    # tempo di risposta inerziale
    t_inerzia = inerzia / beta

    # potenza ottica per portare la coppia a kT (Pe_r = 1)
    potenza_pe1_dich = n.P_LASER * n.KT / tau_dich
    potenza_pe1_mie = n.P_LASER * n.KT / tau_mie

    return {
        "rotore": nome,
        "materiale": materiale.get("nome", "Cu"),
        "k": k,
        "geometria": {
            "raggio_m": raggio,
            "raggio_idrodinamico_m": r_h,
            "raggio_sfera_equivalente_m": r_eq,
            "volume_m3": volume,
            "r2_m2": r2,
            "lunghezza_filo_m": geo["lunghezza_filo_m"] * k,
        },
        "inerzia": {
            "massa_kg": massa,
            "momento_kgm2": inerzia,
            "densita_kgm3": rho,
        },
        "attrito": {
            "beta_Nms": beta,
            "D_rot_einstein_1_s": D_r,
            "tau_rot_s": 1.0 / (2.0 * D_r),
            "tempo_inerziale_s": t_inerzia,
            "nota": "D_r = kT/beta: stessa beta del termine dissipativo",
        },
        "coppia": {
            "forza_dichiarata_N": F_dich,
            "coppia_dichiarata_Nm": tau_dich,
            "forza_scattering_N": F_mie,
            "coppia_scattering_Nm": tau_mie,
            "rapporto_scattering_dichiarata": F_mie / F_dich if F_dich else 0.0,
            "coppia_su_kT_dichiarata": tau_dich / n.KT,
            "coppia_su_kT_scattering": tau_mie / n.KT,
        },
        "regime_libero": {
            "omega_dichiarata_rad_s": omega_dich,
            "omega_scattering_rad_s": omega_mie,
            "Pe_r_dichiarata": omega_dich / D_r,
            "Pe_r_scattering": omega_mie / D_r,
            "giri_anno_dichiarata": omega_dich / (2 * math.pi) * 3.15576e7,
            "giri_anno_scattering": omega_mie / (2 * math.pi) * 3.15576e7,
            "energia_su_kT_dichiarata": 0.5 * inerzia * omega_dich ** 2 / n.KT,
            "potenza_per_Pe_r_1_dichiarata_W": potenza_pe1_dich,
            "potenza_per_Pe_r_1_scattering_W": potenza_pe1_mie,
            "verdetto_dichiarata": "deterministico" if omega_dich > D_r else "termico",
            "verdetto_scattering": "deterministico" if omega_mie > D_r else "termico",
            "nota": "Pe_r = omega/D_r = tau/kT: criterio indipendente da beta.  "
                    "Con tau = F lambda/(2 pi) cresce con l'area intercettata e "
                    "satura quando la struttura copre il fascio.",
        },
        "terra": {
            "g_m_s2": G_TERRA,
            "peso_N": massa * G_TERRA,
            "peso_apparente_N": peso_netto(massa, volume, G_TERRA),
            "spinta_peso_dichiarata": F_dich / (massa * G_TERRA),
            "spinta_peso_scattering": F_mie / (massa * G_TERRA),
            "v_sedimentazione_m_s": v_sed,
            "Pe_gravitazionale": Pe_grav,
            "galleggia": rho < RHO_FLUIDO,
        },
        "spazio": {
            "g_m_s2": 0.0,
            "peso_apparente_N": 0.0,
            "spinta_peso": None,
            "spinta_peso_nota": "indefinito: in g = 0 il peso e' nullo e il "
                                "rapporto spinta/peso diverge",
            "v_sedimentazione_m_s": 0.0,
            "Pe_gravitazionale": 0.0,
            "numero_di_grashof": 0.0,
            "nota": "beta e tau non dipendono da g: omega di regime e' identico "
                    "a terra.  In microgravita' spariscono peso, sedimentazione "
                    "e convezione, non l'attrito viscoso.",
        },
        "meccanica": {
            "resistenza_materiale_Pa": sig,
            "omega_limite_rad_s": w_lim,
            "velocita_limite_m_s": w_lim * raggio,
            "tensione_a_regime_Pa": tensione(omega_dich),
            "margine": sig / tensione(omega_dich) if tensione(omega_dich) > 0 else None,
            "margine_scattering": sig / tensione(omega_mie) if tensione(omega_mie) > 0 else None,
        },
    }


# ---------------------------------------------------------------------------
# 5. REGIME ACCOPPIATO: omega2 / omega1 = -I1 / I2 = -45/376
# ---------------------------------------------------------------------------
def regime_accoppiato(spin, contro):
    """Velocita' angolari stazionarie del sistema accoppiato.

    Vincolo di bilancio del momento angolare:  I1*w1 + I2*w2 = 0, cioe'
    w2 = -(I1/I2) w1.  Sostituendo nei due bilanci con la coppia di scambio
    tau_c si ottiene

        w1 = (tau1 + tau2) / (beta1 - beta2 * I1/I2)

    che e' positivo e ben definito perche' beta1 > beta2 * I1/I2.
    """
    I1 = spin["inerzia"]["momento_kgm2"]
    I2 = contro["inerzia"]["momento_kgm2"]
    b1 = spin["attrito"]["beta_Nms"]
    b2 = contro["attrito"]["beta_Nms"]
    rapporto = I1 / I2

    def risolvi(tau1, tau2):
        den = b1 - b2 * rapporto
        w1 = (tau1 + tau2) / den
        w2 = -rapporto * w1
        tau_c = b1 * w1 - tau1
        return dict(
            omega_spin_rad_s=w1,
            omega_contro_rad_s=w2,
            omega_contro_su_spin=w2 / w1,
            coppia_di_scambio_Nm=tau_c,
            potenza_dissipata_W=0.5 * (b1 * w1 ** 2 + b2 * w2 ** 2),
            denominatore_Nms=den,
        )

    return {
        "I1_su_I2": rapporto,
        "target_45_376": n.RATIO,
        "scarto_rapporto_inerzia": rapporto - n.RATIO,
        "dichiarata": risolvi(spin["coppia"]["coppia_dichiarata_Nm"],
                              contro["coppia"]["coppia_dichiarata_Nm"]),
        "scattering": risolvi(spin["coppia"]["coppia_scattering_Nm"],
                              contro["coppia"]["coppia_scattering_Nm"]),
        "nota": "il rapporto I1/I2 dipende solo dalla geometria (la densita' si "
                "semplifica): 45/376 e' una proprieta' della forma, non del "
                "materiale ne' della scala.",
    }


# ---------------------------------------------------------------------------
# 6. I TRE CASI DI SCALA
# ---------------------------------------------------------------------------
MATERIALE_CU = {"rho": RHO_CU, "sig": SIGMA_CU, "E": 130e9, "n": 0.0, "k": 0.0,
                "cp": 385.0, "nome": "Cu"}


def caso_scala(nome, k, etichetta):
    spin = caso_rotore("Spin", k, MATERIALE_CU)
    contro = caso_rotore("Contro-spin", k, MATERIALE_CU)
    rif = n.scala_geometrica(nome, k)
    return {
        "nome": nome,
        "etichetta": etichetta,
        "k": k,
        "ingombro_m": 1.9e-6 * k,
        "ambiente": "microgravita' con fluido (g = 0)",
        "spin": spin,
        "contro": contro,
        "accoppiato": regime_accoppiato(spin, contro),
        "riferimento_nucleo": {
            "R_disco_m": rif["R_disco_m"],
            "beta_telaio_Nms": rif["beta_telaio_Nms"],
            "coppia_ottica_Nm": rif["tau_ottica_Nm"],
            "omega_ss_rad_s": rif["omega_ss_rad_s"],
            "nota": "v5_nucleo applica R_DISC a entrambi i rotori; qui ogni "
                    "rotore usa il proprio raggio esterno (3.70 u.a. Spin, "
                    "3.275 u.a. Contro-spin).  Lo scarto fra i due omega e' il "
                    "solo fattore geometrico (R_rotore/R_DISC)^3.",
        },
    }


def casi_scala():
    return {
        "10^-9_m": caso_scala("10^-9_m", 1e-3, "dispositivo da 1 nm"),
        "massa_10^-18_kg": caso_scala("massa_10^-18_kg", 1.0,
                                      "dispositivo V4: Spin 1.246e-18 kg, "
                                      "Contro-spin 4.488e-18 kg"),
        "1_m": caso_scala("1_m", 1e6, "pseudo-motore da 1 m di ingombro"),
    }


# ---------------------------------------------------------------------------
# 7. I TRE MATERIALI A 1 METRO
# ---------------------------------------------------------------------------
def materiali_a_1m():
    k = 1e6
    risultato = {"k": k, "ingombro_m": 1.9e-6 * k, "materiali": {}}
    for nome_mat, dati in MATERIALI.items():
        mat = dict(dati)
        mat["nome"] = nome_mat
        spin = caso_rotore("Spin", k, mat)
        contro = caso_rotore("Contro-spin", k, mat)
        risultato["materiali"][nome_mat] = {
            "proprieta": dati,
            "spin": spin,
            "contro": contro,
            "accoppiato": regime_accoppiato(spin, contro),
            "massa_totale_kg": spin["inerzia"]["massa_kg"] + contro["inerzia"]["massa_kg"],
            "indice_specifico_kJ_kg": dati["sig"] / dati["rho"] / 1000.0,
        }
    return risultato


# ---------------------------------------------------------------------------
# 8. CONFRONTO TERRA / SPAZIO - SINTESI
# ---------------------------------------------------------------------------
def confronto_terra_spazio():
    righe = []
    for nome, caso in casi_scala().items():
        for chiave in ("spin", "contro"):
            r = caso[chiave]
            righe.append({
                "caso": nome,
                "rotore": r["rotore"],
                "raggio_m": r["geometria"]["raggio_m"],
                "massa_kg": r["inerzia"]["massa_kg"],
                "beta_Nms": r["attrito"]["beta_Nms"],
                "omega_libero_dich_rad_s": r["regime_libero"]["omega_dichiarata_rad_s"],
                "Pe_r_dich": r["regime_libero"]["Pe_r_dichiarata"],
                "Pe_r_scattering": r["regime_libero"]["Pe_r_scattering"],
                "spinta_peso_terra_dich": r["terra"]["spinta_peso_dichiarata"],
                "spinta_peso_terra_scattering": r["terra"]["spinta_peso_scattering"],
                "spinta_peso_spazio": "non definito (peso nullo)",
                "v_sedimentazione_m_s": r["terra"]["v_sedimentazione_m_s"],
                "Pe_gravitazionale": r["terra"]["Pe_gravitazionale"],
            })
    return righe


def riepilogo_accoppiato():
    fuori = {}
    for nome, caso in casi_scala().items():
        fuori[nome] = {
            "dichiarata": caso["accoppiato"]["dichiarata"],
            "scattering": caso["accoppiato"]["scattering"],
            "I1_su_I2": caso["accoppiato"]["I1_su_I2"],
            "target_45_376": caso["accoppiato"]["target_45_376"],
        }
    fuori["1_m_materiali"] = {
        mat: dati["accoppiato"] for mat, dati in materiali_a_1m()["materiali"].items()
    }
    return fuori


# ---------------------------------------------------------------------------
# 9. SCRITTURA
# ---------------------------------------------------------------------------
def scrivi(nome, payload):
    percorso = os.path.join(OUT, nome)
    with open(percorso, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False, default=float)
    print("  scritto %-30s %6.1f kB" % (nome, os.path.getsize(percorso) / 1024.0))


def main():
    print("=" * 68)
    print("MOTORE FOTONICO V5 - MODULO SPAZIALE (g = 0)")
    print("=" * 68)

    casi = casi_scala()

    verifica_masse = {}
    for chiave in ("Spin", "Contro-spin"):
        r = casi["massa_10^-18_kg"]["spin" if chiave == "Spin" else "contro"]
        m_calc = r["inerzia"]["massa_kg"]
        m_v4 = ROTORI[chiave]["massa_v4_kg"]
        verifica_masse[chiave] = {
            "massa_calcolata_kg": m_calc,
            "massa_preprint_v4_kg": m_v4,
            "scarto_relativo": m_calc / m_v4 - 1,
        }

    scrivi("SP1_spazio_scale.json", {
        "descrizione": "Spin e Contro-spin in microgravita' alle tre scale richieste",
        "ambiente": "g = 0, fluido presente (eta = %.1e Pa s)" % n.ETA,
        "criterio_determinismo": {
            "definizione": "Pe_r = omega_ss / D_r = tau / kT",
            "D_r": "kT / beta, con la stessa beta del termine dissipativo "
                   "tau_diss = -beta*omega (relazione di Einstein)",
            "controprova": "Pe_r = omega_ss/D_r = (tau/beta)/(kT/beta) = tau/kT: "
                           "il criterio non dipende dall'attrito.",
            "attenzione": "accoppiare beta ridotta da FILL_FRAME con "
                          "D_r = kT/(8 pi eta R^3) gonfia Pe_r: il fattore di "
                          "gonfiamento e' quantificato in SP2 (riconciliazione).",
        },
        "verifica_masse_v4": verifica_masse,
        "casi": casi,
    })

    confronto = confronto_terra_spazio()
    soglia = n.soglia_rotazionale()
    pe_cont = casi["massa_10^-18_kg"]["contro"]["regime_libero"]["Pe_r_scattering"]
    scrivi("SP2_terra_vs_spazio.json", {
        "descrizione": "cio' che cambia fra terra e microgravita'",
        "forza_dichiarata_N": n.F_DECL,
        "nota": "beta e tau non dipendono da g: la rotazione e' identica. "
                "In assenza di gravita' spariscono peso, sedimentazione e "
                "convezione; l'attrito viscoso resta.",
        "riconciliazione_Pe_r": {
            "Pe_r_su_tau": pe_cont,
            "Pe_r_Stokes": soglia["Pe_rot"],
            "D_rot_Stokes_1_s": soglia["D_rot_1_s"],
            "raggio_equivalente_m": n.R_EFF_V4,
            "fattore": soglia["Pe_rot"] / pe_cont,
            "causa": "la stessa coppia con due definizioni di D_r: "
                     "tau/kT usa la beta ridotta da FILL_FRAME = %.4f, mentre "
                     "D_r = kT/(8 pi eta r_eff^3) usa il raggio equivalente "
                     "nudo.  Il fattore e' 1/FILL_FRAME; il criterio del testo "
                     "resta Pe_r = tau/kT." % n.FILL_FRAME,
        },
        "righe": confronto,
    })

    mat_1m = materiali_a_1m()
    scrivi("SP3_materiali_1m.json", mat_1m)

    accoppiato = riepilogo_accoppiato()
    scrivi("SP4_rapporto_accoppiato.json", {
        "descrizione": "regime stazionario con il vincolo I1 w1 + I2 w2 = 0",
        "formula": "w1 = (tau1 + tau2) / (beta1 - beta2 * I1/I2),  w2 = -(I1/I2) w1",
        "target_rapporto": n.RATIO,
        "risultati": accoppiato,
    })

    sintesi = {
        "descrizione": "sintesi del modulo spaziale",
        "scale": list(casi),
        "materiali_1m": list(MATERIALI),
        "valori_chiave": {
            "spinta_peso_terra_v4_scattering": casi["massa_10^-18_kg"]["spin"]["terra"]["spinta_peso_scattering"],
            "spinta_peso_terra_v4_dichiarata": casi["massa_10^-18_kg"]["spin"]["terra"]["spinta_peso_dichiarata"],
            "omega_accoppiato_v4_dich": casi["massa_10^-18_kg"]["accoppiato"]["dichiarata"]["omega_spin_rad_s"],
            "omega_accoppiato_v4_scattering": casi["massa_10^-18_kg"]["accoppiato"]["scattering"]["omega_spin_rad_s"],
            "Pe_r_1nm_scattering": casi["10^-9_m"]["spin"]["regime_libero"]["Pe_r_scattering"],
            "Pe_r_v4_scattering": casi["massa_10^-18_kg"]["spin"]["regime_libero"]["Pe_r_scattering"],
            "Pe_r_1m_scattering": casi["1_m"]["spin"]["regime_libero"]["Pe_r_scattering"],
        },
    }
    scrivi("SP5_sintesi_spazio.json", sintesi)

    # ---- stampa di controllo ------------------------------------------------
    print("\n%-26s %-11s %12s %12s %10s" %
          ("caso", "rotore", "omega_dich", "omega_scatt", "Pe_r_scatt"))
    for nome, caso in casi.items():
        for chiave in ("spin", "contro"):
            r = caso[chiave]
            print("%-26s %-11s %12.3e %12.3e %10.3e" % (
                nome, r["rotore"],
                r["regime_libero"]["omega_dichiarata_rad_s"],
                r["regime_libero"]["omega_scattering_rad_s"],
                r["regime_libero"]["Pe_r_scattering"]))

    print("\nRegime accoppiato (vincolo I1 w1 + I2 w2 = 0):")
    for nome, val in accoppiato.items():
        if nome == "1_m_materiali":
            continue
        d = val["dichiarata"]
        print("  %-18s w1 = %10.3e rad/s   w2 = %10.3e rad/s   w2/w1 = %+.5f"
              % (nome, d["omega_spin_rad_s"], d["omega_contro_rad_s"],
                 d["omega_contro_su_spin"]))

    print("\nMateriali a 1 m:")
    for nome_mat, dati in mat_1m["materiali"].items():
        s = dati["spin"]
        print("  %-24s massa %8.3f kg  I1 %9.3e kg m2  omega_lim %9.3e rad/s"
              % (nome_mat, s["inerzia"]["massa_kg"], s["inerzia"]["momento_kgm2"],
                 s["meccanica"]["omega_limite_rad_s"]))

    print("\nTerra contro spazio (dispositivo V4):")
    v4 = casi["massa_10^-18_kg"]
    for chiave in ("spin", "contro"):
        r = v4[chiave]
        print("  %-12s peso %.3e N | peso apparente %.3e N | spinta/peso dich %.3f"
              " | scatt %.4f | v_sed %.3e m/s"
              % (r["rotore"], r["terra"]["peso_N"], r["terra"]["peso_apparente_N"],
                 r["terra"]["spinta_peso_dichiarata"],
                 r["terra"]["spinta_peso_scattering"],
                 r["terra"]["v_sedimentazione_m_s"]))

    print("\nFatto: 5 file JSON in", OUT)


if __name__ == "__main__":
    main()
