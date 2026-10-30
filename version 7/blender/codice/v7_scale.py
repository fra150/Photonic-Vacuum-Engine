"""
MOTORE FOTONICO DEL VUOTO - V7.1
v7_scale.py - i due estremi: il motore da 10^-9 m e il motore da 1 m.

Il resto della V7 vive fra 0.1 e 30 um.  Qui si chiudono le due code, con lo
stesso formalismo e gli stessi criteri: nessuna formula nuova, solo lo stesso
codice di v7_mappa.py applicato a R = 0.5 nm (diametro 1 nm) e R = 0.5 m
(diametro 1 m).  E' il test che dice se la finestra e' una proprieta' della
fisica o della scala.

Uso:
    python v7_scale.py
"""

import json
import math
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
if BASE not in sys.path:
    sys.path.insert(0, BASE)

import v7_mappa as m  # noqa: E402

DATI = os.path.normpath(os.path.join(BASE, "..", "dati"))
os.makedirs(DATI, exist_ok=True)

SCALE = (("1_nm", 0.5e-9, "diametro 1 nm"),
         ("1_m", 0.5, "diametro 1 m"))

MATERIALI = ("cu_lucido", "cu_nero", "silice", "vaterite")
FORME = ("sfera", "disco")
CANALI = ("assorbimento", "ritardo")
MEZZI = (("acqua", None), ("vuoto", 1e-6), ("vuoto", 1e-8))

P_LAVORO = 10e-3          # W, il tetto dei criteri


def potenza_pe10(canale, R, materiale, forma, **kw):
    """Potenza incidente per Pe_r = 10 (la coppia e' lineare in P)."""
    tau1 = m.coppia(canale, 1.0, R, materiale, forma, **kw)["tau"]
    return 10.0 * m.KT / tau1 if tau1 > 0 else float("inf")


def caso(R, materiale, forma, mezzo, p_mbar, canale):
    c = m.cella(canale, P_LAVORO, R, materiale, forma, mezzo, p_mbar)
    c["P_pe10"] = potenza_pe10(canale, R, materiale, forma)
    c["P_pe10_mW"] = c["P_pe10"] * 1e3
    return c


def blocco(etichetta, R, descrizione):
    righe = []
    for materiale in MATERIALI:
        for forma in FORME:
            for canale in CANALI:
                for mezzo, p in MEZZI:
                    righe.append(caso(R, materiale, forma, mezzo, p, canale))
    return {"nome": etichetta, "R_m": R, "descrizione": descrizione,
            "P_lavoro_W": P_LAVORO, "righe": righe}


def scala_1nm(R=0.5e-9):
    """Perche' a 1 nm nessun canale arriva a Pe_r = 10.

    Due fattori indipendenti, entrambi geometrici:
      eta_int ~ 2R^2/w0^2        la sezione d'urto e' 1e-6 del fascio
      delta   ~ 2 pi dn 2R/lambda  il ritardo e' 1e-3 rad, e la coppia va
                                   come delta^2 per la polarizzazione circolare
    Il prodotto e' ~ 1e-12 rispetto alla scala micrometrica.
    """
    b = blocco("1_nm", R, "diametro 1 nm")
    righe = []
    for materiale in MATERIALI:
        for forma in FORME:
            c = caso(R, materiale, forma, "acqua", None, "assorbimento")
            eta = m.eta_int(R, forma)
            righe.append({
                "materiale": materiale, "forma": forma,
                "eta_int": eta,
                "eta_int_su_w0_quadro": eta / ((R / m.W0) ** 2),
                "Pe_r_a_10mW": c["Pe_r"], "omega_a_10mW": c["omega"],
                "dT_a_10mW": c["dT"], "P_pe10_W": c["P_pe10"],
                "nota": "Pe_r < 1: il rotore non esiste, e' rumore browniano",
            })
    b["analisi"] = righe
    b["conclusione"] = (
        "a 1 nm nessuno dei tre canali raggiunge Pe_r = 10 nemmeno a 10 mW. "
        "Il limite non e' la coppia ma la sezione d'urto relativa al fascio "
        "(eta_int ~ (R/w0)^2 ~ 1e-6) e, per il ritardo di fase, il ritardo "
        "stesso (delta ~ 1e-3 rad).  Il 'motore da 1 nm' della V5 era un "
        "rapporto fra numeri, non un rotore: a quella scala il momento "
        "angolare ottico per fotone e' intatto, ma il fotone non incontra "
        "l'oggetto.")
    return b


def scala_1m(R=0.5):
    """Cosa cambia a 1 m: Pe_r enorme, rotazione inesistente.

    La coppia cresce come P (non come R): a 1 m la potenza si distribuisce su
    un'area enorme e eta_int = 1 da sempre.  L'attrito cresce come R^3, quindi
    omega = tau/beta crolla come R^-3.  Il rapporto omega/omega_ott e' la
    misura di quanto il rotore sia 'trascinato' e non guidato.
    """
    b = blocco("1_m", R, "diametro 1 m")
    righe = []
    for materiale in MATERIALI:
        for forma in FORME:
            c = caso(R, materiale, forma, "acqua", None, "ritardo")
            if not c["tau"]:
                continue
            righe.append({
                "materiale": materiale, "forma": forma,
                "massa_kg": m.MATERIALI[materiale]["rho"] * m.geometria(forma, R)["V"],
                "tau_10mW": c["tau"], "Pe_r": c["Pe_r"],
                "omega": c["omega"], "omega_rottura": c["omega_rottura"],
                "omega_frac": c["omega_frac"], "eta_mecc": c["eta_mecc"],
                "t_spin_a_meta_rottura_s": c["t_spin"],
                "giri_anno": c["omega"] / (2 * math.pi) * 3.15576e7,
                "P_pe10_W": c["P_pe10"],
            })
    b["analisi"] = righe
    b["conclusione"] = (
        "a 1 m Pe_r e' enorme (la coppia e' deterministica) ma la rotazione e' "
        "inesistente: omega ~ 1e-13 rad/s, un giro ogni milioni di anni, e "
        "omega/omega_ott ~ 1e-29.  Il corpo e' un oggetto ORIENTATO, non un "
        "rotore.  La coppia ottica non scala con R: e' la potenza che conta, "
        "e a 1 m la potenza per unita' di area e' nulla.")
    return b


def sweep(n=31, materiale="vaterite", forma="sfera", canale="ritardo",
          mezzo="acqua", p_mbar=None, P=P_LAVORO):
    """Da 1 nm a 1 m in 12 decadi, a potenza fissata."""
    raggi = np.logspace(-9, 0, n)
    righe = []
    for R in raggi:
        c = m.cella(canale, P, float(R), materiale, forma, mezzo, p_mbar)
        righe.append({
            "R_m": float(R),
            "eta_int": c["eta_int"], "tau": c["tau"], "Pe_r": c["Pe_r"],
            "omega": c["omega"], "omega_frac": c["omega_frac"],
            "eta_mecc": c["eta_mecc"], "dT": c["dT"], "T": c["T"],
            "valido": c["valido"], "motivi": c["motivi"],
        })
    # massimo di Pe_r e massimo di omega: dove sta la finestra
    i_pe = max(range(n), key=lambda i: righe[i]["Pe_r"])
    i_w = max(range(n), key=lambda i: righe[i]["omega"])
    return {"canale": canale, "materiale": materiale, "forma": forma,
            "mezzo": mezzo, "P_W": P, "righe": righe,
            "max_Pe_r": {"R_m": righe[i_pe]["R_m"], "Pe_r": righe[i_pe]["Pe_r"]},
            "max_omega": {"R_m": righe[i_w]["R_m"], "omega": righe[i_w]["omega"]}}


def cards():
    """I numeri che finiscono nelle schede renderizzate (dati/C1_cards.json).

    Il modello 3D non deve contenere numeri scritti a mano: questa funzione li
    calcola e v7_blender.py li legge.
    """
    out = {}

    # scheda del rotore: vaterite R = 2 um, acqua, 30 mW
    c = m.cella("ritardo", 30e-3, 2e-6, "vaterite", "sfera", "acqua", s=1.0)
    co = m.coppia("ritardo", 30e-3, 2e-6, "vaterite", "sfera")
    out["rotore"] = {
        "titolo": "V7 - rotore a ritardo di fase",
        "righe": [
            "vaterite, R = 2 um, acqua 310 K",
            "P = 30 mW, polarizzazione circolare",
            "ritardo delta = %.2f rad" % co["delta"],
            "tau = %.2e N m" % c["tau"],
            "f = %.1f Hz" % (c["omega"] / (2 * math.pi)),
            "Pe_r = %.0f" % c["Pe_r"],
            "dT = %.1f K" % c["dT"],
        ],
    }

    # scheda 1 nm
    c1 = caso(0.5e-9, "cu_nero", "sfera", "acqua", None, "assorbimento")
    out["1_nm"] = {
        "titolo": "V7 - il motore da 1 nm",
        "righe": [
            "R = 0.5 nm, diametro 1 nm",
            "lambda = 532 nm: R/lambda = %.1e" % (0.5e-9 / m.LAMBDA),
            "eta_int = %.1e" % m.eta_int(0.5e-9, "sfera"),
            "Pe_r a 10 mW = %.1e" % c1["Pe_r"],
            "per Pe_r = 10 servono %.1f kW" % (c1["P_pe10"] / 1e3),
            "nessun rotore: e' rumore browniano",
        ],
    }

    # scheda 1 m
    c2 = caso(0.5, "vaterite", "sfera", "acqua", None, "ritardo")
    out["1_m"] = {
        "titolo": "V7 - il motore da 1 m",
        "righe": [
            "disco/sfera R = 0.5 m, diametro 1 m",
            "P = 10 mW, acqua 310 K",
            "tau = %.2e N m" % c2["tau"],
            "Pe_r = %.0f (deterministico)" % c2["Pe_r"],
            "omega = %.1e rad/s" % c2["omega"],
            "eta_mecc = %.1e" % c2["eta_mecc"],
            "oggetto orientato, non rotore",
        ],
    }

    # scheda manubrio (test F).  Lo scattering e' lecito qui perche' il
    # manubrio non e' assiale; la coppia e' un limite superiore.
    c3 = caso(0.1e-6, "silice", "manubrio", "vuoto", 1e-8, "scattering")
    w_lim, _ = m.omega_rottura(0.1e-6, "silice", "manubrio")
    out["manubrio"] = {
        "titolo": "V7.1 - manubrio di silice, 200 nm",
        "righe": [
            "silice, R = 100 nm, vuoto 1e-8 mbar",
            "omega_rottura = %.2f GHz" % (w_lim / 2 / math.pi / 1e9),
            "I = %.1e kg m^2" % c3["I"],
            "tau a 10 mW = %.1e N m" % c3["tau"],
            "t_spin a meta' limite = %.2e s" % c3["t_spin"],
            "limite assoluto 2P/omega: %.1e s" % c3["t_spin_limite"],
            "il limite e' il materiale, non la coppia",
        ],
    }
    return out


def _scrivi(nome, payload):
    p = os.path.join(DATI, nome)
    with open(p, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1, ensure_ascii=False, default=float)
    print("  scritto %-26s %6.1f kB" % (nome, os.path.getsize(p) / 1024.0))


def main():
    print("=" * 70)
    print("V7 - SCALA: 10^-9 m e 1 m")
    print("=" * 70)
    n1 = scala_1nm()
    _scrivi("S1_scala_1nm.json", n1)
    m1 = scala_1m()
    _scrivi("S3_scala_1m.json", m1)

    print("\n-- 1 nm (diametro), 10 mW, acqua, canale di assorbimento")
    print("   %-10s %-6s %10s %10s %10s" % ("materiale", "forma", "eta_int", "Pe_r", "P(Pe=10)"))
    for r in n1["analisi"]:
        print("   %-10s %-6s %10.2e %10.2e %10.2e W"
              % (r["materiale"], r["forma"], r["eta_int"], r["Pe_r_a_10mW"], r["P_pe10_W"]))

    print("\n-- 1 m (diametro), 10 mW, acqua, canale a ritardo")
    print("   %-10s %-6s %12s %10s %12s %10s" % ("materiale", "forma", "tau (Nm)", "Pe_r", "omega (rad/s)", "eta_mecc"))
    for r in m1["analisi"]:
        print("   %-10s %-6s %12.3e %10.2e %12.3e %10.2e"
              % (r["materiale"], r["forma"], r["tau_10mW"], r["Pe_r"], r["omega"], r["eta_mecc"]))

    sw = sweep()
    _scrivi("S5_sweep_scala.json", sw)
    print("\n-- sweep 1 nm -> 1 m (vaterite sfera, acqua, 10 mW, ritardo)")
    print("   max Pe_r  = %.3e a R = %.2e m" % (sw["max_Pe_r"]["Pe_r"], sw["max_Pe_r"]["R_m"]))
    print("   max omega = %.3e rad/s a R = %.2e m" % (sw["max_omega"]["omega"], sw["max_omega"]["R_m"]))

    sunto = {
        "versione": "V7.1",
        "descrizione": "i due estremi della V7.1, stesso formalismo e stessi criteri",
        "criteri": m.CRITERI,
        "scale": {"1_nm": n1, "1_m": m1},
        "sweep": sw,
        "risposta": {
            "1_nm": "Pe_r ~ 1e-7 a 10 mW: nessun rotore.  eta_int ~ (R/w0)^2 "
                    "domina, non la coppia per fotone.",
            "1_m": "Pe_r ~ 1e5 ma omega ~ 1e-13 rad/s e eta_mecc ~ 1e-29: "
                   "oggetto orientato, non rotore.",
            "finestra": "la rotazione massima sta fra 1 e 100 um; fuori da li' "
                        "o non c'e' coppia (sotto) o non c'e' velocita' (sopra).",
        },
    }
    _scrivi("S6_sintesi_scala.json", sunto)
    cd = cards()
    _scrivi("C1_cards.json", cd)
    print("\n-- schede per il modello 3D")
    for k, v in cd.items():
        print("   %-10s %s" % (k, " | ".join(v["righe"][:2])))
    print("\nFatto.")


if __name__ == "__main__":
    main()
