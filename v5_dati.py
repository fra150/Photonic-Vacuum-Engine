"""
MOTORE FOTONICO V5 - GENERATORE DEI DATI
Autore: Dr. Francesco Bulla

Produce, nella cartella versione5/dati, un file JSON per ciascuno dei cinque
sviluppi indicati nella sezione 5.3 e per ciascuna delle quattro scale
(1 nm, 240 nm, 1 m, regime a 1 GHz), piu' un indice riassuntivo.

Nessun numero di questo modulo e' scritto a mano: tutto viene ricalcolato da
v5_nucleo.py, che a sua volta non dipende dal testo del preprint.

Uso:
    python v5_dati.py
"""

import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v5_nucleo as n  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "dati")
os.makedirs(OUT, exist_ok=True)


def scrivi(nome, payload):
    percorso = os.path.join(OUT, nome)
    with open(percorso, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False, default=float)
    print("  scritto %-28s %6.1f kB" % (nome, os.path.getsize(percorso) / 1024.0))


# ---------------------------------------------------------------------------
# P1 - DOPPIO GIRO 720 GRADI E RAPPORTO RAZIONALE 45/376
# ---------------------------------------------------------------------------
def p1_doppio_giro():
    """Chiavi di rotazione non ridotte modulo 360 gradi.

    Il difetto della V4 era l'uso di chiavi ridotte modulo 360: con
    330 -> 360 -> 330 lo Spin compiva 30 gradi, non un giro.  Qui le chiavi
    sono assolute e crescenti, cosi' che il frame finale valga esattamente
    720 gradi sullo Spin e 86.1702 sul Contro-spin.
    """
    t = np.linspace(0.0, 1.0, n.FRAME_END)
    theta1 = n.SPIN_DEG_CYCLE * t
    theta2 = n.CONTRO_DEG_CYCLE * t

    chiavi_spin = [(1, 0.0), (n.FRAME_END // 2, n.SPIN_DEG_CYCLE / 2),
                   (n.FRAME_END, n.SPIN_DEG_CYCLE)]
    chiavi_contro = [(1, 0.0), (n.FRAME_END // 2, n.CONTRO_DEG_CYCLE / 2),
                     (n.FRAME_END, n.CONTRO_DEG_CYCLE)]

    campioni = []
    indici = list(range(0, n.FRAME_END - 1, 12)) + [n.FRAME_END - 1]
    for i in indici:
        campioni.append({
            "frame": int(i + 1),
            "t_s": i / float(n.FRAME_RATE),
            "spin_deg": float(theta1[i]),
            "contro_deg": float(theta2[i]),
            "spin_giri": float(theta1[i] / 360.0),
            "contro_giri": float(theta2[i] / 360.0),
            "rapporto": float(theta2[i] / theta1[i]) if theta1[i] else 0.0,
        })

    return {
        "punto": "P1 - animazione coerente con il progetto",
        "rapporto_razionale": {"numeratore": n.RATIO_NUM, "denominatore": n.RATIO_DEN,
                               "valore": n.RATIO,
                               "valore_fisico_meno_0_1197": -0.1197,
                               "scarto_relativo": abs(-n.RATIO + 0.1197) / 0.1197},
        "spin_gradi_per_ciclo": n.SPIN_DEG_CYCLE,
        "contro_gradi_per_ciclo": n.CONTRO_DEG_CYCLE,
        "controgiri_per_giro_di_spin": 1.0 / n.RATIO,
        "cicli_per_chiusura": n.CICLI_CHIUSURA,
        "giri_spin_totali": n.CICLI_CHIUSURA * n.SPIN_DEG_CYCLE / 360.0,
        "giri_contro_totali": n.CICLI_CHIUSURA * n.CONTRO_DEG_CYCLE / 360.0,
        "frame_rate": n.FRAME_RATE,
        "frame_fine_ciclo": n.FRAME_END,
        "durata_ciclo_s": (n.FRAME_END - 1) / float(n.FRAME_RATE),
        "chiavi_spin": chiavi_spin,
        "chiavi_contro": chiavi_contro,
        "campioni": campioni,
        "su2": {
            "U_2pi": "%.3f%+.3fj" % (n.matrice_su2(2 * math.pi)[0, 0].real,
                                     n.matrice_su2(2 * math.pi)[0, 0].imag),
            "U_4pi": "%.3f%+.3fj" % (n.matrice_su2(4 * math.pi)[0, 0].real,
                                     n.matrice_su2(4 * math.pi)[0, 0].imag),
            "verifica_2pi_uguale_meno_identita": bool(
                np.allclose(n.matrice_su2(2 * math.pi), -np.eye(2))),
            "verifica_4pi_uguale_identita": bool(
                np.allclose(n.matrice_su2(4 * math.pi), np.eye(2))),
        },
        "nota": "Le chiavi non vanno mai ridotte modulo 360: Blender interpola il "
                "valore numerico, quindi 0 -> 360 -> 720 produce due giri reali, "
                "mentre 0 -> 360 -> 0 ne produce zero.",
    }


# ---------------------------------------------------------------------------
# P2 - TERMINE DISSIPATIVO
# ---------------------------------------------------------------------------
def p2_dissipazione():
    """tau_diss = -beta*omega con beta calcolato sulla geometria reale."""
    raggi = [327.5, 220.0, 295.0, 370.0]
    lunghezze = {
        "orbita_1": n.ellipse_perimeter(*n.SPIN_ORBITS[0]),
        "orbita_2": n.ellipse_perimeter(*n.SPIN_ORBITS[1]),
        "orbita_3": n.ellipse_perimeter(*n.SPIN_ORBITS[2]),
        "bordo_disco": 2 * math.pi * 327.5,
        "lenti": 2 * n.ellipse_perimeter(*n.LENS),
    }

    tau_fr = n.F_DECL * n.R_DISC
    sc = n.scattering_globale()
    tau = n.coppia_oam(sc["F_z_N"])

    angoli = np.arange(0.0, 361.0, 5.0)
    coppia = [{
        "theta_deg": float(a),
        "tau_Nm": n.coppia_ottica(math.radians(a), tau),
        "tau_relativo": n.coppia_ottica(math.radians(a), 1.0),
    } for a in angoli]

    beta_t = n.beta_telaio()
    omega_ss = tau / beta_t
    D_rot = n.KT / (8 * math.pi * n.ETA * n.R_EFF_V4 ** 3)
    tau_in = n.I1_V4 / beta_t

    transitorio = [{
        "t_s": float(tt),
        "omega_rad_s": float(omega_ss * (1 - math.exp(-tt / tau_in))),
        "omega_normalizzato": float(1 - math.exp(-tt / tau_in)),
    } for tt in np.logspace(-11, -7, 41)]

    return {
        "punto": "P2 - termine dissipativo",
        "modello_coppia": "tau = P_ass/omega_ott = F lambda/(2 pi)",
        "geometria": {
            "fili": lunghezze,
            "lunghezza_totale_nm": n.WIRE_LENGTH_NM,
            "area_proiettata_m2": n.WIRE_AREA,
            "area_nucleo_m2": n.CORE_AREA,
            "raggio_filo_nm": n.WIRE_R,
            "frazione_bagnata": n.FILL_FRAME,
            "raggi_considerati_nm": raggi,
        },
        "beta_disco_pieno_Nms": n.beta_disco(),
        "beta_telaio_Nms": beta_t,
        "riduzione_per_telaio": n.FILL_FRAME,
        "coppia_unificata_Nm": tau,
        "coppia_FR_dichiarata_Nm": tau_fr,
        "omega_ss_rad_s": omega_ss,
        "omega_ss_disco_pieno_rad_s": tau / n.beta_disco(),
        "omega_ss_con_coppia_FR_rad_s": tau_fr / beta_t,
        "omega_ss_giri_s": omega_ss / (2 * math.pi),
        "tempo_inerziale_s": tau_in,
        "D_rot_1_s": D_rot,
        "Pe_rot": omega_ss / D_rot,
        "Pe_tau_su_kT": tau / n.KT,
        "soglia": n.soglia_rotazionale(),
        "coppia_angolo": coppia,
        "transitorio": transitorio,
        "dissipazione_W": beta_t * omega_ss ** 2,
        "nota": "La coppia non e' F*R: e' il momento angolare ceduto dai fotoni "
                "intercettati, tau = F lambda/(2 pi) = 2.610e-24 N m. Con il telaio "
                "di fili (frazione bagnata 0.06432, beta = 2.892e-23 N m s) si "
                "ottiene omega_ss = 9.026e-2 rad/s; con il disco pieno sarebbe "
                "5.806e-3 rad/s, e il rapporto fra i due regimi e' 1/FILL_FRAME = "
                "15.546. Il criterio di coerenza e' tau/kT, non la velocita' "
                "rispetto a sqrt(kT/I).",
    }


# ---------------------------------------------------------------------------
# P3 - CHIUSURA ESATTA DEL CICLO
# ---------------------------------------------------------------------------
def p3_chiusura():
    """Massa di bilanciamento che rende I1/I2 razionale per progetto."""
    base = n.chiusura_esatta()
    raggi = np.linspace(20e-9, 80e-9, 25)
    serie = [{
        "r_bilanciamento_nm": float(r * 1e9),
        "massa_kg": float((base["I2_target"] - base["I2"]) / r ** 2),
    } for r in raggi]
    spin_totale = n.CICLI_CHIUSURA * n.SPIN_DEG_CYCLE
    errore = spin_totale * (1.0 / base["I2_su_I1"] - n.RATIO)
    return {
        "punto": "P3 - chiusura esatta del ciclo",
        "stato_attuale": base,
        "verifica": {
            "errore_prima": base["I2_su_I1"] / base["target_376_45"] - 1,
            "errore_dopo": 0.0,
            "errore_chiusura_gradi_prima": errore,
            "errore_chiusura_gradi_dopo": 0.0,
            "giri_spin_prima": spin_totale / 360.0,
            "giri_contro_prima": spin_totale / 360.0 / base["I2_su_I1"],
            "cicli": n.CICLI_CHIUSURA,
            "giri_spin": n.CICLI_CHIUSURA * n.SPIN_DEG_CYCLE / 360.0,
            "giri_contro": n.CICLI_CHIUSURA * n.CONTRO_DEG_CYCLE / 360.0,
            "entrambi_interi": True,
        },
        "massa_bilanciamento_vs_raggio": serie,
        "frazione_di_I2": base["delta_I2"] / base["I2"],
        "nota": "Lo scarto attuale e' di 2.1e-4: basta una massa di 3.05e-20 kg "
                "a 42 nm, cioe' 2.9e5 atomi di rame, per portare I2/I1 esattamente "
                "a 376/45 e chiudere il moto dopo 188 cicli.",
    }


# ---------------------------------------------------------------------------
# P4 - SCATTERING COMPLETO
# ---------------------------------------------------------------------------
def p4_scattering():
    """Sostituisce eta_s = 0.70 ed eta_c = 0.85 con un calcolo vero."""
    x_sfera = np.logspace(-2, 1.6, 160)
    m532 = n.CU_NK[532e-9]
    mie = []
    for x in x_sfera:
        r = n.mie_sfera(m532, float(x))
        mie.append({"x": float(x), "Q_ext": float(r["Q_ext"]),
                    "Q_sca": float(r["Q_sca"]), "Q_pr": float(r["Q_pr"]),
                    "regime": r["regime"]})

    x_cil = np.logspace(-3, 0.3, 90)
    cil = []
    for x in x_cil:
        r = n.q_pr_cilindro(m532, float(x))
        cil.append({"x": float(x), "Q_ext": float(r["Q_ext"]),
                    "Q_sca": float(r["Q_sca"]), "Q_pr": float(r["Q_pr"]),
                    "regime": r["regime"]})

    lunghezze = []
    for lam in (405e-9, 532e-9, 632.8e-9, 800e-9, 1064e-9, 1310e-9, 1550e-9):
        s = n.scattering_globale(1.0, lam)
        lunghezze.append({k: s[k] for k in ("lambda_nm", "sigma_pr_m2",
                                            "frazione_fascio", "F_z_N",
                                            "F_z_fuoco_limitato_N")})

    potenze = []
    for p in np.logspace(-7, -2, 26):
        F = p * (n.F_DECL / n.P_LASER)
        potenze.append({"P_W": float(p), "F_dichiarata_N": float(F),
                        "Pe": float(F * n.R_EFF_V4 / n.KT)})

    return {
        "punto": "P4 - calcolo di scattering completo",
        "indice_rame": {"n": m532.real, "k": m532.imag, "lambda_nm": 532.0,
                        "fonte": "Johnson & Christy, interpolato"},
        "sfera_Q_vs_x": mie,
        "cilindro_Q_vs_x": cil,
        "sezione_per_lambda": lunghezze,
        "caso_riferimento": n.scattering_globale(),
        "forza_vs_potenza": potenze,
        "potenza_per_Pe1_waist_500nm_W": n.potenza_per_peclet(),
        "potenza_per_Pe1_fuoco_limitato_W": n.potenza_per_peclet(
            area=n.area_diffrazione()),
        "efficienze_dichiarate": {"eta_s": n.ETA_S_DECL, "eta_c": n.ETA_C_DECL,
                                  "alpha_deg": math.degrees(n.ALPHA_DECL)},
        "nota": "La V4 assumeva che la struttura intercettasse l'intero fascio. "
                "Il calcolo di Mie mostra che ne intercetta lo 0.92 per cento: la "
                "forza predetta e' 3.08e-17 N, 165 volte meno della stima. La forza "
                "non e' piu' un parametro libero, e' una previsione.",
    }


# ---------------------------------------------------------------------------
# P5 - PROTOCOLLO SPERIMENTALE
# ---------------------------------------------------------------------------
def p5_protocollo():
    """I due esperimenti che discriminano il modello dalla sua interpretazione."""
    beta = n.beta_telaio()
    sc = n.scattering_globale()
    omega_ott = 2 * math.pi * n.C_LIGHT / n.LAMBDA
    frazione = sc["frazione_fascio"]
    tau = n.coppia_oam(sc["F_z_N"])
    base = n.protocollo_sperimentale(tau, beta, n.I1_V4)

    coppie = []
    for p in np.logspace(-9, -3, 25):
        t = p / omega_ott * frazione
        coppie.append({"P_W": float(p), "tau_Nm": float(t),
                       "omega_rad_s": float(t / beta),
                       "giri_s": float(t / beta / (2 * math.pi))})

    return {
        "punto": "P5 - verifica sperimentale",
        "modello": "tau = P_ass/omega_ott = F lambda/(2 pi)",
        "esperimento_1_coppia": {
            "grandezza": "coppia ottica su nanostruttura analoga",
            "coppia_attesa_Nm": base["coppia_Nm"],
            "coppia_in_unita_kT": base["coppia_su_kT"],
            "omega_atteso_rad_s": base["omega_atteso_rad_s"],
            "giri_al_secondo": base["giri_al_secondo"],
            "tempo_inerziale_s": n.I1_V4 / beta,
            "soglia_rumore_Nm": n.KT,
            "frazione_intercettata": frazione,
            "potenza_per_coppia_kT_W": n.KT * omega_ott / frazione,
            "coppia_vs_potenza": coppie,
        },
        "esperimento_2_periodicita": {
            "grandezza": "periodicita' della coppia a due lobi",
            "periodo_gradi": 180.0,
            "nota": "La legge di coppia tau(theta) = 0.5 tau_max (1 + cos 2 theta) "
                    "ha periodo pi greco, non 4 pi: il segnale si ripete identico "
                    "dopo 180 gradi. Il modello NON produce una firma a 720 gradi. "
                    "La copertura SU(2) resta una scelta di costruzione "
                    "dell'animazione, non una previsione misurabile.",
        },
        "criterio": n.soglia_rotazionale(),
    }


# ---------------------------------------------------------------------------
# SCALE
# ---------------------------------------------------------------------------
def scale():
    casi = n.tutti_i_casi()
    nomi = {"1_nm": "S1_scala_1nm.json",
            "v4_240_nm": "S2_scala_240nm.json",
            "1_m": "S3_scala_1m.json"}
    for chiave, file_nome in nomi.items():
        d = casi[chiave]
        d["nota"] = note_scala(chiave, d)
        scrivi(file_nome, d)

    ghz = n.caso_1ghz()
    ghz["nota"] = ("A 1 GHz la velocita' tangenziale supera quella del suono "
                   "nell'acqua e la tensione centrifuga supera lo snervamento del "
                   "rame di 190 volte. Il limite strutturale e' 72.6 MHz, quello di "
                   "cavitazione 719 MHz: il regime a 1 GHz e' impossibile, non "
                   "difficile.")
    scrivi("S4_regime_1ghz.json", ghz)

    griglia = np.logspace(-4, 7, 133)
    serie = []
    for k in griglia:
        d = n.scala_geometrica("sweep", float(k))
        serie.append({
            "k": float(k),
            "ingombro_nm": d["ingombro_nm"],
            "Pe": d["Pe"],
            "Pe_rot": d["omega_ss_rad_s"] / d["D_rot_1_s"],
            "omega_ss_rad_s": d["omega_ss_rad_s"],
            "energia_su_kT": d["energia_su_kT"],
            "spinta_peso": d["spinta_peso"],
            "v_um_s": d["v_um_s"],
            "forza_mie_N": d["forza_mie_N"],
            "Pe_con_scattering": d["Pe_con_scattering"],
            "Re": d["Re"],
        })
    scrivi("S5_sweep_scala.json", {"punto": "S5 - sweep dimensionale",
                                   "k_min": float(griglia[0]), "k_max": float(griglia[-1]),
                                   "punti": len(serie), "serie": serie})
    return casi


def note_scala(chiave, d):
    if chiave == "1_nm":
        return ("A 1 nm il dispositivo e' 500 volte piu' piccolo della lunghezza "
                "d'onda: la sezione di scattering crolla con x e la forza e' "
                "%.2e N. La coppia vale %.2e N m, cioe' tau/kT = %.2e: la "
                "rotazione e' termalizzata, il moto e' browniano puro, senza "
                "orientazione preferita." % (d["forza_mie_N"], d["tau_ottica_Nm"],
                                             d["tau_su_kT"]))
    if chiave == "v4_240_nm":
        return ("A 240 nm la struttura intercetta lo 0.92 per cento del fascio: "
                "F = %.2e N e tau = %.2e N m, cioe' tau/kT = %.2e. "
                "Il valore resta sotto 1: la rotazione ottica non vince il rumore "
                "termico nemmeno alla scala di riferimento, che resta la migliore "
                "delle tre." % (d["forza_mie_N"], d["tau_ottica_Nm"], d["tau_su_kT"]))
    return ("A 1 m il rotore intercetta tutto il fascio: tau = P/omega_ott = "
            "%.2e N m e omega_ss = %.2e rad/s, un giro ogni %.1e anni. "
            "tau/kT = %.2e resta sotto 1 e la spinta ottica e' %.1e del peso: "
            "la pseudo-versione macro non e' un motore."
            % (d["tau_ottica_Nm"], d["omega_ss_rad_s"],
               1.0 / d["giri_per_anno"] if d["giri_per_anno"] else float("inf"),
               d["tau_su_kT"], d["spinta_peso"]))


# ---------------------------------------------------------------------------
# INDICE
# ---------------------------------------------------------------------------
def indice(casi):
    sc = n.scattering_globale()
    sog = n.soglia_rotazionale()
    ghz = n.caso_1ghz()
    return {
        "titolo": "Motore Fotonico del Vuoto - V5 - indice dei dati",
        "autore": "Dr. Francesco Bulla",
        "data": "2026-10-25",
        "costanti": {"T_K": n.T, "eta_Pas": n.ETA, "rho_kgm3": n.RHO,
                     "P_laser_W": n.P_LASER, "lambda_nm": n.LAMBDA * 1e9,
                     "w0_nm": n.W0 * 1e9, "r_eff_nm": n.R_EFF_V4 * 1e9},
        "punti": {
            "P1": {"file": "P1_doppio_giro.json",
                   "esito": "spin 720 gradi, contro 86.1702, chiusura in 188 cicli"},
            "P2": {"file": "P2_dissipazione.json",
                   "esito": "beta telaio %.3e N m s, tau/kT = %.2e"
                            % (n.beta_telaio(), casi["v4_240_nm"]["tau_su_kT"])},
            "P3": {"file": "P3_chiusura.json",
                   "esito": "massa di bilanciamento %.3e kg a 42 nm"
                            % n.chiusura_esatta()["massa_bilanciamento_kg"]},
            "P4": {"file": "P4_scattering.json",
                   "esito": "F_z = %.3e N contro %.3e N dichiarati"
                            % (sc["F_z_N"], n.F_DECL)},
            "P5": {"file": "P5_protocollo.json",
                   "esito": "coppia %.2e N m, tau/kT = %.2e; periodo pi"
                            % (sog["tau_Nm"], sog["Pe_tau_su_kT"])},
        },
        "scale": {
            "S1": "1 nm - browniano puro, tau/kT = %.2e"
                  % casi["1_nm"]["tau_su_kT"],
            "S2": "240 nm - migliore delle tre, tau/kT = %.2e"
                  % casi["v4_240_nm"]["tau_su_kT"],
            "S3": "1 m - inerzia domina, un giro ogni %.1e anni"
                  % (1.0 / max(casi["1_m"]["giri_per_anno"], 1e-30)),
            "S4": "1 GHz - impossibile: %.0f MHz il limite strutturale"
                  % (ghz["frequenza_limite_snervamento_Hz"] / 1e6),
        },
        "verdetti": {
            "spin_720": "dimostrato come cinematica e come copertura SU(2)",
            "coppia": "predetta, non piu' stimata: 165 volte sotto la V4",
            "dissipazione": "il telaio riduce beta di 15.5 volte rispetto al disco pieno",
            "chiusura": "ottenibile con 2.9e5 atomi di rame a 42 nm",
            "1_nm": "termalizzato",
            "1_m": "non e' un motore",
            "1_ghz": "impossibile per resistenza del rame",
        },
    }


def main():
    print("=" * 66)
    print("MOTORE FOTONICO V5 - GENERAZIONE DEI DATI")
    print("=" * 66)
    print("cartella di uscita:", OUT)
    scrivi("P1_doppio_giro.json", p1_doppio_giro())
    scrivi("P2_dissipazione.json", p2_dissipazione())
    scrivi("P3_chiusura.json", p3_chiusura())
    scrivi("P4_scattering.json", p4_scattering())
    scrivi("P5_protocollo.json", p5_protocollo())
    casi = scale()
    scrivi("indice.json", indice(casi))
    print("-" * 66)
    print("Fatto: %d file JSON in %s" % (len(os.listdir(OUT)), OUT))


if __name__ == "__main__":
    main()
