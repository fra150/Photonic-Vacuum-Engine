"""
MOTORE FOTONICO DEL VUOTO - V8.0
v8_mappa.py - nucleo fisico e mappa della finestra valida.

Autore: Dr. Francesco Bulla

La V6 ha chiuso il motore.  La V7 apre il rotore: un oggetto che gira per
RITARDO DI FASE (come la vaterite di Friese 1998), non per assorbimento
(come il rame della V5).  Il contributo non e' "il motore funziona": e' "il
rotore funziona, e la mappa dice dove".

V8.0 technical correction of V7.1: vaterite Δn=0.06 effective, T_m=723 K,
beam-validity R<=2 µm, sphere Epstein drag (16/3), Rayleigh crossover 50 nm,
Pe_r defined with T_bath=310 K, physical efficiency η_lim, trap flag.

Correzioni V8.0 (rispetto alla V7.1)
------------------------------------
1. VATERITE.  dn = 0.06 e' la birifrangenza EFFICACE di un aggregato
   policristallino (range 0.06-0.10), non 0.17: 0.17 e' la calcite, non la
   vaterite.  k = 1e-6 (trasparente, conservativo).  T_m = 723 K e' l'inizio
   della transizione vaterite -> calcite, non 1000 K.
2. VALIDITA' DEL FASCIO.  R <= 2 um: oltre, la particella e' piu' grande
   della regione in cui il fascio gaussiano e' ben descritto (z_R ~ 2 um).
   Ogni cella con R > 2 um e' invalida per geometria, non per la fisica del
   rotore.
3. ATTRITO MOLECOLARE.  Sfera: beta = (16/3) p R^4 / v_bar (Epstein, drag
   rotazionale, sigma_t = 1).  Per le forme non sferiche si usa lo stesso
   valore e si marca beta_shape_incerto: non si inventano fattori di forma.
   Le affermazioni quantitative in vuoto valgono solo per sfere.
4. EMISSIVITA' A TEMPERATURA.  eps_eff(materiale, R, T): il peso di Planck
   usa T, e in temperatura l'emissivita' viene ricalcolata a ogni passo.
   Il Cu nero resta il parametro dichiarato 0.30, indipendente da T.
5. Pe_r.  Pe = tau / (k_B T_bath) con T_bath = 310 K (stesso valore di T0,
   ma la definizione e' esplicita).  Non cambia i numeri; cambia il testo.
6. EFFICIENZA FISICA.  eta_lim = tau omega_rottura / P e' il rendimento
   limite del rotore, non eta_mecc = omega/omega_ottico.  Nessuna cella
   VALIDA supera 1e-7: valido_rotore (V7.1) resta per confronto, ma
   valido_rotore_fisico e valido_fisico sono falsi ovunque.  Solo alcune
   celle invalide di Cu nero (assorbimento, P > 10 mW o fuori criterio
   termico) arrivano a ~1.1e-7, e non contano come rotori.
7. TRAPPOLA.  F_scat = P/c, peso = rho V g, trap_ratio = F_scat/peso.  In
   vuoto la trappola e' necessaria; il numero dice quanto e' marginale.
8. R_RAYLEIGH = 50 nm: e' un raccordo di taglia, non un regime fisico.  La
   mappa parte da R = 0.1 um, quindi la correzione e' essenzialmente
   inattiva; il bordo inferiore della finestra lo decide la soglia di coppia.

Correzioni V7.1 (storiche, mantenute)
-------------------------------------
1. SCATTERING ASSIALE RIMOSSO.  Per una particella assialmente simmetrica
   illuminata lungo l'asse, J_z si conserva: il campo diffuso non porta via
   momento angolare netto e la coppia di scattering e' nulla (Marston e
   Crichton, PRA 30, 2508, 1984).  Il canale resta solo per il manubrio, che
   non e' assiale; per sfere, dischi, anelli e reticoli vale "non applicabile".
2. INERZIA DEL MANUBRIO CORRETTA.  kI era in m^2 e veniva moltiplicato ancora
   per R^2: I risultava piu' piccolo di ~R^2 e t_spin di oltre 10 ordini.
3. RITARDO DI FASE REALISTICO.  Nella V7.1, sotto R = 0.15 um il ritardo
   efficace scala come R^3 (polarizzabilita' di Rayleigh), non come R; la
   formula a singolo passaggio sovrastima (1-cos delta) fino a ~5x a 100 nm.
   La V8 porta il raccordo a 50 nm e lo chiama per quello che e': un
   raccordo di taglia, non un regime fisico.  Le oscillazioni
   di delta con la taglia non sono un effetto fisico da riportare in mappa:
   si media su una dispersione larga (sigma_delta = 0.5 delta).  Il caso
   monodisperso delta = pi resta come limite ideale, non come previsione.
4. CU NERO.  eps_IR = 0.30 e' un parametro DICHIARATO, non un risultato di
   Mie: un film da 50 nm su uno specchio emette poco nell'IR e servono strati
   di spessore micron.  La mappa usa 0.30 (conservativo); il bordo termico
   scala come eps^-1/2 e finisce intorno a 30 um.
5. PRESSIONE RESIDUA.  P_MBAR era 0.1 invece di 100: la pressione effettiva
   era 1000 volte piu' bassa dell'etichetta "1e-8 mbar".  Corretto.  Non
   cambia le finestre di orientamento (che non usano beta), ma cambia
   omega, eta_mecc e il conteggio dei rotori completi in vuoto.

Le tre funzioni pubbliche di questo modulo
------------------------------------------
  eta_int(R, forma, w0)                  frazione di fascio intercettata
  coppia(canale, P, R, ...)             coppia per i tre canali
  temperatura(P_abs, R, ...)            delta T in acqua / T in vuoto

Tutto il resto e' contorno: geometrie, attrito, meccanica, criteri, griglia.

I tre canali di coppia, in un solo formalismo
---------------------------------------------
1. ASSORBIMENTO   tau = (l + s) P_abs / omega_ott
   E' il canale della V5 e della V6.  Riscalda e non usa il fascio.

2. RITARDO DI FASE   oggetto trasparente e anisotropo, ritardo delta
       tau = (P_irr / omega_ott) * [ |S3| (1-cos d) + |sin d| sqrt(S1^2+S2^2) ]
   con gli Stokes del fascio in ingresso.  La somma e' lineare, non in
   quadratura: il termine di spin non dipende dall'orientamento e quello
   lineare si somma ad esso.
     - circolare  (S3 = +-1):  tau = (P/omega)(1 - cos d)   [limite di Beth 1936:
       lamina lambda/2 che inverte l'elica -> 2 hbar per fotone]
     - lineare    (S1^2+S2^2 = 1): tau = (P/omega) sin d     [lamina lambda/2 a
       45 gradi -> 0: nessun momento angolare trasferito, e infatti sin(pi)=0;
       lamina lambda/4 a 45 gradi -> P/omega: la lineare diventa circolare]
   La formula sin^2(d/2) usata in letteratura per il caso lineare NON
   riproduce il caso lambda/2 a 45 gradi (darebbe P/omega invece di 0) e non
   e' usata qui.  Per il caso circolare le due coincidono: 1-cos d = 2 sin^2(d/2).

3. SCATTERING ANISOTROPO   tau = (P_sca/omega_ott) (1 - g) (l + s)
   (1-g) e' la frazione di luce non diretta in avanti; il fattore 2 sullo
   spin viene dal backscattering, che inverte l'elica.  E' il caso limite e
   vale SOLO per una forma non assiale (manubrio).  Per una sfera o un disco
   con l'asse lungo il fascio la simmetria conserva J_z e la coppia e' nulla.

Cosa la V7 NON fa
-----------------
Non produce spinta: F <= 2P/c per qualunque fascio.  La propulsione e' fuori.
Non include l'attrito interno del materiale: in vuoto il rotore e' un volano
che accelera fino alla rottura, e la mappa riporta il tempo di spin-up.

Uso:
    python v8_mappa.py            # self-test + mappa completa
    python v8_mappa.py test       # solo self-test
"""

import csv
import cmath
import json
import math
import os
import sys
from functools import lru_cache

import numpy as np
from scipy import constants as C
from scipy.special import spherical_jn, spherical_yn

BASE = os.path.dirname(os.path.abspath(__file__))
DATI = os.path.normpath(os.path.join(BASE, "..", "dati_v8"))
os.makedirs(DATI, exist_ok=True)

# ---------------------------------------------------------------------------
# 0. COSTANTI
# ---------------------------------------------------------------------------
LAMBDA = 532e-9                        # m
OMEGA_OTT = 2.0 * math.pi * C.c / LAMBDA     # 3.5449e15 rad/s
NU_OTT = C.c / LAMBDA                  # 5.635e14 Hz
W0 = 500e-9                            # m, waist dichiarato (fascio gaussiano)
NA = 0.9
W0_DIFF = 0.61 * LAMBDA / NA           # 3.6e-7 m, limite di diffrazione
T0 = 310.0                             # K, acqua e ambiente
T_BATH = 310.0                         # K, temperatura del bagno per Pe_r
KT = C.k * T0
SIGMA_SB = C.sigma

# Raccordo di taglia per il ritardo di fase.  Sotto R_RAYLEIGH si usa
# delta ~ R^3 come limite di polarizzabilita'; NON e' un confine di regime
# fisico e non e' il bordo inferiore della mappa (che parte da R = 0.1 um).
# Il crossover a 50 nm e' quindi essenzialmente inattivo sulla griglia: il
# bordo inferiore della finestra lo decide la soglia di coppia (Pe_r >= 10).
R_RAYLEIGH = 50e-9                     # m
SIGMA_RITARDO = 0.50
N_DISP = 31

# Dominio del fascio gaussiano: oltre 2 um la particella e' piu' grande della
# regione in cui il fascio e' ben descritto (z_R ~ 2 um per w0 = 500 nm).
R_MAX_VALID = 2e-6                     # m

MU_ACQUA = 6.95e-4                     # Pa s, acqua a 310 K (non 1.2e-3 a 293 K)
K_ACQUA = 0.628                        # W/m/K, acqua a 310 K
M_GAS = 4.65e-26                       # kg, N2
P_MBAR = 100.0                         # 1 mbar = 100 Pa (V7.1: corretto da 0.1)

DDL = 1e-9                             # m, un nanometro

# ---------------------------------------------------------------------------
# 1. MATERIALI
#    n, k a 532 nm; alpha532 = frazione della potenza intercettata assorbita
#    (1 - R per i metalli, dal solo indice; il nero e' un rivestimento
#    dichiarato, non una proprieta' del rame).
# ---------------------------------------------------------------------------
MATERIALI = {
    "cu_lucido": {
        "nome": "Cu lucido",
        "n": 1.09, "k": 2.53, "alpha532": None,     # calcolato da n, k
        "rho": 8960.0, "sy": 200e6, "kappa": 401.0, "cp": 385.0,
        "Tm": 1358.0, "dn": 0.0,
        "nota": "assorbe il 40% perche' e' opaco; il resto lo riflette. "
                "Trasparente: no.  Canale a ritardo non applicabile.",
    },
    "cu_nero": {
        "nome": "Cu nero (strato spesso, eps_IR dichiarato)",
        "n": 1.09, "k": 2.53, "alpha532": 0.95,
        "rho": 8960.0, "sy": 200e6, "kappa": 401.0, "cp": 385.0,
        "Tm": 1358.0, "dn": 0.0,
        "nota": "Il rame nero non e' una proprieta' di Mie: eps_IR = 0.30 e' "
                "un parametro dichiarato, conservativo per un rivestimento "
                "di spessore micron.  Un film da 50 nm su uno specchio "
                "assorbe poco nell'IR (eps ~ 0.05-0.1) e non basta: servono "
                "strati spessi.  Il bordo termico scala come eps^-1/2, "
                "quindi con 0.30 sta intorno a 30 um, non a 21 um.",
        "emiss_ir": 0.30,
        "emiss_ir_incertezza": (0.20, 0.60),
    },
    "silice": {
        "nome": "Silice fusa",
        "n": 1.461, "k": 1.0e-7, "alpha532": None,
        "rho": 2200.0, "sy": 1.0e9, "kappa": 1.4, "cp": 740.0,
        "Tm": 1983.0, "dn": 0.0,
        "nota": "isotropa: il ritardo di fase e' solo di FORMA (manubrio, "
                "disco, reticolo), non di materiale.",
    },
    "vaterite": {
        "nome": "Vaterite (CaCO3)",
        "n": 1.65, "k": 1.0e-6, "alpha532": None,
        "rho": 2710.0, "sy": 100e6, "kappa": 2.2, "cp": 800.0,
        "Tm": 723.0, "dn": 0.06,
        "nota": "dn = 0.06 e' la birifrangenza EFFICACE di un aggregato "
                "policristallino (range 0.06-0.10).  Il valore 0.17 e' quello "
                "della calcite, non della vaterite: usarlo qui sovrastimava la "
                "coppia di un fattore ~3.  k = 1e-6 (trasparente, "
                "conservativo).  T_m = 723 K e' l'inizio della transizione "
                "vaterite -> calcite, non 1000 K.  sy = 100 MPa resta un "
                "valore conservativo per un cristallo policristallino fragile.",
    },
    "diamante": {
        "nome": "Diamante + assorbitore",
        "n": 2.42, "k": 0.02, "alpha532": 0.10,
        "rho": 3510.0, "sy": 1.5e9, "kappa": 2000.0, "cp": 520.0,
        "Tm": 3800.0, "dn": 0.0,
        "nota": "isotropo (cubico): nessun ritardo intrinseco.  Entra solo "
                "nel canale di assorbimento, e solo se gli altri non bastano.",
    },
}

# Proprieta' ottiche nell'infrarosso termico (2-25 um), per l'emissivita'.
# Silice, vaterite e diamante: valori tabulati approssimati (n, k).
TAB_IR = {
    "silice": [
        (2.0, 1.44, 0.0), (2.5, 1.43, 0.0), (3.0, 1.42, 0.0), (3.5, 1.40, 0.0),
        (4.0, 1.38, 1e-3), (5.0, 1.34, 5e-3), (6.0, 1.30, 2e-2), (7.0, 1.25, 5e-2),
        (8.0, 1.18, 0.12), (9.0, 1.05, 0.35), (9.6, 0.80, 0.70), (10.6, 0.95, 0.45),
        (12.0, 1.25, 0.20), (14.0, 1.45, 0.28), (17.0, 1.60, 0.35), (20.0, 1.70, 0.30),
        (25.0, 1.85, 0.25),
    ],
    "vaterite": [
        (2.0, 1.68, 0.0), (2.5, 1.67, 0.0), (3.0, 1.66, 0.0), (3.5, 1.65, 1e-3),
        (4.0, 1.63, 3e-3), (5.0, 1.60, 8e-2), (6.0, 1.55, 0.25), (7.0, 1.45, 0.55),
        (8.0, 1.50, 0.30), (9.0, 1.55, 0.25), (10.0, 1.60, 0.22), (12.0, 1.70, 0.30),
        (14.0, 1.75, 0.40), (17.0, 1.80, 0.45), (20.0, 1.85, 0.40), (25.0, 1.90, 0.35),
    ],
    "diamante": [
        (2.0, 2.38, 0.0), (2.5, 2.38, 0.0), (3.0, 2.38, 0.0), (3.5, 2.38, 5e-4),
        (4.0, 2.38, 2e-3), (5.0, 2.38, 0.02), (6.0, 2.39, 0.10), (7.0, 2.40, 0.30),
        (8.0, 2.40, 0.20), (9.0, 2.41, 0.08), (10.0, 2.41, 0.05), (12.0, 2.41, 0.03),
        (14.0, 2.41, 0.02), (17.0, 2.41, 0.02), (20.0, 2.41, 0.03), (25.0, 2.41, 0.04),
    ],
}
LAM_IR = (2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0,
          14.0, 17.0, 20.0, 25.0)

# Ritardo di fase per forma (birifrangenza di forma, adimensionale).
# Per la silice l'anisotropia e' geometrica: una sfera e' isotropa, un
# manubrio no.  Valori dichiarati, incertezza +-50%.
DN_FORMA = {
    ("silice", "sfera"): 0.0,
    ("silice", "disco"): 0.020,
    ("silice", "anello"): 0.010,
    ("silice", "manubrio"): 0.030,
    ("silice", "reticolo"): 0.010,
}
DN_FORMA_DEFAULT = 0.0
DN_INCERTEZZA = 0.5

# ---------------------------------------------------------------------------
# 2. GEOMETRIE
#    R e' il raggio esterno.  Ogni forma dichiara:
#      V        volume di materiale
#      A_proj   area proiettata sul fascio
#      A_rad    superficie radiante REALE (non 4 pi R^2, correzione V6.1)
#      kI       coefficiente ADIMENSIONALE: I = kI m R^2.  Per il manubrio
#               kI = 2(0.4 r^2 + d^2)/(2 m_lobo) / R^2 = 0.4(r/R)^2+(d/R)^2;
#               la V7.0 lo lasciava in m^2 e lo moltiplicava ancora per R^2.
#      c_mecc   coefficiente di tensione: sigma_max = c rho v^2
#      percorso cammino ottico lungo il fascio (per il ritardo di fase)
#      elementi lista di bersagli per l'integrale di sovrapposizione
# ---------------------------------------------------------------------------
FORME = ("sfera", "disco", "anello", "manubrio", "reticolo")


def geometria(forma, R):
    if forma == "sfera":
        V = 4.0 / 3.0 * math.pi * R ** 3
        elementi = [("disco", R, 0.0)]
        A_proj = math.pi * R ** 2
        A_rad = 4.0 * math.pi * R ** 2
        kI, c_mecc, percorso = 0.4, 0.40, 2.0 * R
    elif forma == "disco":
        t = 0.2 * R
        V = math.pi * R ** 2 * t
        elementi = [("disco", R, 0.0)]
        A_proj = math.pi * R ** 2
        A_rad = 2.0 * math.pi * R ** 2 + 2.0 * math.pi * R * t
        kI, c_mecc, percorso = 0.5, 0.40, t
    elif forma == "anello":
        a = 0.1 * R
        V = 2.0 * math.pi ** 2 * R * a ** 2
        elementi = [("anello", R, a)]
        A_proj = 4.0 * math.pi * R * a
        A_rad = 4.0 * math.pi ** 2 * R * a
        kI, c_mecc, percorso = 1.0, 1.0, 2.0 * a
    elif forma == "manubrio":
        r, d = 0.25 * R, 0.75 * R
        V = 2.0 * (4.0 / 3.0 * math.pi * r ** 3)
        elementi = [("disco_offset", r, d), ("disco_offset", r, -d)]
        A_proj = 2.0 * math.pi * r ** 2
        A_rad = 2.0 * 4.0 * math.pi * r ** 2
        # I = 2 m_lobo (0.4 r^2 + d^2), quindi kI = 0.4(r/R)^2+(d/R)^2
        kI = 0.4 * (r / R) ** 2 + (d / R) ** 2
        c_mecc, percorso = 0.50, 2.0 * r
    elif forma == "reticolo":
        a = 0.03 * R
        raggi = (0.45 * R, 0.75 * R, R)
        elementi = [("anello", r, a) for r in raggi]
        elementi += [("razze", R, a, 6)]
        V = sum(2.0 * math.pi ** 2 * r * a ** 2 for r in raggi)
        V += 6.0 * (R * (0.06 * R) * (0.06 * R))
        A_proj = sum(4.0 * math.pi * r * a for r in raggi)
        A_proj += 6.0 * (R * 0.06 * R)
        A_rad = sum(4.0 * math.pi ** 2 * r * a for r in raggi)
        A_rad += 6.0 * 2.0 * (R * 0.06 * R)
        kI, c_mecc, percorso = 0.80, 1.00, 2.0 * a
    else:
        raise ValueError("forma sconosciuta: %s" % forma)

    Req = (3.0 * V / (4.0 * math.pi)) ** (1.0 / 3.0)
    return {"R": R, "V": V, "A_proj": A_proj, "A_rad": A_rad, "kI": kI,
            "c_mecc": c_mecc, "percorso": percorso, "Req": Req,
            "elementi": elementi}


# ---------------------------------------------------------------------------
# 3. OTTICA: MIE ESATTA (sezione d'urto, emissivita' spettrale)
# ---------------------------------------------------------------------------
X_GEOM = 300.0       # oltre: ottica geometrica con Q = 2, 1, 1, g = 0.5


def _fresnel_medio(m, n=180):
    """Riflettanza media su una sfera in ottica geometrica (media su 2 sin cos)."""
    num = den = 0.0
    for i in range(n):
        th = (i + 0.5) / n * math.pi / 2.0
        ct, st = math.cos(th), math.sin(th)
        st2 = st / m
        ct2 = cmath.sqrt(1.0 - st2 * st2)
        rs = (ct - m * ct2) / (ct + m * ct2)
        rp = (m * ct - ct2) / (m * ct + ct2)
        R = 0.5 * (abs(rs) ** 2 + abs(rp) ** 2)
        w = 2.0 * st * ct
        num += R * w
        den += w
    return num / den if den > 0 else 0.0


def mie_sfera(m, x):
    """Coefficienti di Mie per una sfera omogenea (Bohren & Huffman, cap. 4).

    Restituisce Q_ext, Q_sca, Q_abs e il parametro di asimmetria g.
    Oltre x = 300 si usa il limite geometrico con la riflettanza di Fresnel
    mediata sulla sfera: Q_abs = 1 - R, Q_sca = 1 + R, Q_ext = 2.  Non e'
    Q_abs = 1: un corpo metallico grande resta uno specchio.
    """
    if x <= 0.0:
        return {"Q_ext": 0.0, "Q_sca": 0.0, "Q_abs": 0.0, "g": 0.0,
                "regime": "nullo"}
    if x > X_GEOM:
        Rb = _fresnel_medio(m)
        return {"Q_ext": 2.0, "Q_sca": 1.0 + Rb, "Q_abs": 1.0 - Rb, "g": 0.5,
                "regime": "ottica geometrica (Fresnel)"}
    nmax = max(int(2.0 + x + 4.0 * x ** (1.0 / 3.0)) + 10, 12)
    nx = m * x
    D = np.zeros(nmax + 2, dtype=complex)
    for nn in range(nmax, 0, -1):
        D[nn - 1] = nn / nx - 1.0 / (D[nn] + nn / nx)

    def psi(nn, z):
        return z * spherical_jn(nn, z)

    def xi(nn, z):
        return z * (spherical_jn(nn, z) + 1j * spherical_yn(nn, z))

    a = np.zeros(nmax + 2, dtype=complex)
    b = np.zeros(nmax + 2, dtype=complex)
    for nn in range(1, nmax + 1):
        a[nn] = ((D[nn] / m + nn / x) * psi(nn, x) - psi(nn - 1, x)) / \
               ((D[nn] / m + nn / x) * xi(nn, x) - xi(nn - 1, x))
        b[nn] = ((D[nn] * m + nn / x) * psi(nn, x) - psi(nn - 1, x)) / \
               ((D[nn] * m + nn / x) * xi(nn, x) - xi(nn - 1, x))
    idx = np.arange(1, nmax + 1)
    q_ext = (2.0 / x ** 2) * float(np.sum((2 * idx + 1) * (a[1:nmax + 1] + b[1:nmax + 1]).real))
    q_sca = (2.0 / x ** 2) * float(np.sum((2 * idx + 1) *
                                          (np.abs(a[1:nmax + 1]) ** 2 + np.abs(b[1:nmax + 1]) ** 2)))
    g = 0.0
    if q_sca > 1e-30:
        t1 = 0.0
        for nn in range(1, nmax):
            t1 += (nn * (nn + 2.0) / (nn + 1.0)) * \
                  float((a[nn] * np.conj(a[nn + 1]) + b[nn] * np.conj(b[nn + 1])).real)
        t2 = 0.0
        for nn in range(1, nmax + 1):
            t2 += ((2.0 * nn + 1.0) / (nn * (nn + 1.0))) * \
                  float((a[nn] * np.conj(b[nn])).real)
        g = (4.0 / (x ** 2 * q_sca)) * (t1 + t2)
    return {"Q_ext": q_ext, "Q_sca": q_sca, "Q_abs": q_ext - q_sca,
            "g": g, "regime": "Mie esatta"}


def nk_rame_ir(lam):
    """Drude: il rame nell'IR e' un metallo quasi libero."""
    wp, gamma = 1.36e16, 1.0 / 27e-15
    w = 2.0 * math.pi * C.c / lam
    eps = 1.0 - wp ** 2 / (w ** 2 + 1j * gamma * w)
    m = cmath.sqrt(complex(eps))
    return complex(abs(m.real), abs(m.imag))


def nk_ir(nome, lam):
    if nome.startswith("cu"):
        return nk_rame_ir(lam)
    tab = TAB_IR.get(nome)
    if tab is None:
        return complex(1.5, 0.0)
    lam_um = lam * 1e6
    if lam_um <= tab[0][0]:
        return complex(tab[0][1], tab[0][2])
    if lam_um >= tab[-1][0]:
        return complex(tab[-1][1], tab[-1][2])
    for i in range(len(tab) - 1):
        l0, n0, k0 = tab[i]
        l1, n1, k1 = tab[i + 1]
        if l0 <= lam_um <= l1:
            t = (lam_um - l0) / (l1 - l0)
            return complex(n0 + t * (n1 - n0), k0 + t * (k1 - k0))
    return complex(1.5, 0.0)


def planck(lam, T):
    return (2.0 * C.h * C.c ** 2 / lam ** 5) / (math.exp(C.h * C.c / (lam * C.k * T)) - 1.0)


@lru_cache(maxsize=16384)
def _q_abs_ir(materiale, R, lam):
    """Q_abs IR della sfera, memorizzato: il peso di Planck cambia con T,
    i coefficienti di Mie no.  Separare i due evita di ricalcolare Mie a ogni
    iterazione termica."""
    return max(mie_sfera(nk_ir(materiale, lam), 2.0 * math.pi * R / lam)["Q_abs"], 0.0)


@lru_cache(maxsize=8192)
def eps_eff(materiale, R, T=T0):
    """Emissivita' efficace mediata sul peso di Planck a temperatura T.

    Per i materiali tabulati e' calcolata da Mie:
    eps_eff = somma_lambda Q_abs(lambda) * B(lambda, T) * dlambda / somma ...
    Il peso spettrale usa planck(lam, T), non planck(lam, T0): in vuoto la
    particella puo' stare a centinaia di kelvin e lo spettro termico si
    sposta verso il vicino IR, dove n e k cambiano.
    Per il rame la Drude da' Q_abs ~ 3e-3 a R = 1 um: una sfera metallica
    sub-lunghezza d'onda e' uno specchio anche nell'IR, e questo e' il punto.
    Per il Cu nero NON e' Mie: e' il parametro dichiarato emiss_ir = 0.30,
    indipendente da T, perche' l'emissivita' dipende dallo spessore e dalla
    formulazione del rivestimento.
    """
    mat = MATERIALI[materiale]
    if materiale == "cu_nero":
        return mat["emiss_ir"]
    num = den = 0.0
    for i, lam_um in enumerate(LAM_IR):
        lam = lam_um * 1e-6
        lo = LAM_IR[i - 1] if i > 0 else LAM_IR[0]
        hi = LAM_IR[i + 1] if i + 1 < len(LAM_IR) else LAM_IR[-1]
        dl = 0.5 * (hi - lo) * 1e-6
        w = planck(lam, T) * dl
        num += _q_abs_ir(materiale, R, lam) * w
        den += w
    return min(num / den, 1.0) if den > 0 else 0.02


@lru_cache(maxsize=4096)
def q_ottico(materiale, R, lam=LAMBDA):
    """Sezione d'urto della sfera di pari volume: Q_abs, Q_sca, g."""
    mat = MATERIALI[materiale]
    m = complex(mat["n"], mat["k"])
    return mie_sfera(m, 2.0 * math.pi * R / lam)


def alpha_assorbimento(materiale, percorso, lam=LAMBDA):
    """Frazione della potenza intercettata che viene assorbita a 532 nm.

    Due regimi, perche' sono due fisiche diverse:
      opaco (k > 0.05): tutta la potenza che entra viene assorbita nello
                        spessore ottico, quindi alpha = 1 - R di Fresnel;
      trasparente:      assorbe solo per assorbimento volumetrico,
                        alpha = 1 - exp(-4 pi k L / lambda).
    Usare 1 - R anche per un dielettrico trasparente (errore della prima
    stesura) dava alpha = 0.94 per la vaterite e un riscaldamento di 400 K:
    e' il test E che l'ha trovato.
    """
    mat = MATERIALI[materiale]
    if mat["alpha532"] is not None:
        return mat["alpha532"]
    n, k = mat["n"], mat["k"]
    if k > 0.05:
        R = ((n - 1.0) ** 2 + k ** 2) / ((n + 1.0) ** 2 + k ** 2)
        return 1.0 - R
    return 1.0 - math.exp(-4.0 * math.pi * k * percorso / lam)


# ---------------------------------------------------------------------------
# 4. FUNZIONE PUBBLICA 1: eta_int
# ---------------------------------------------------------------------------
@lru_cache(maxsize=8192)
def _eta_disco_offset(r, d, w0, N=241):
    """Sovrapposizione gaussiana-disco decentrato (integrazione 2D)."""
    if abs(d) > 3.0 * w0 + r:
        return 0.0
    L = max(3.0 * r, 2.5 * w0)
    ax = np.linspace(d - L, d + L, N)
    ay = np.linspace(-L, L, N)
    X, Y = np.meshgrid(ax, ay)
    mask = ((X - d) ** 2 + Y ** 2) <= r ** 2
    gauss = (2.0 / (math.pi * w0 ** 2)) * np.exp(-2.0 * (X ** 2 + Y ** 2) / w0 ** 2)
    return float(np.sum(mask * gauss) * (ax[1] - ax[0]) * (ay[1] - ay[0]))


def _eta_razze(R, a, n_razze, w0):
    """Sovrapposizione con n razze radiali di semilarghezza a, lunghezza R.

    Vicino all'asse le razze si sovrappongono: la copertura angolare e'
    f(r) = min(1, 2 n a / (2 pi r)) e l'integrale si fa in forma chiusa.
    Sommare n integrali indipendenti gonfierebbe eta di un fattore n.
    """
    r_m = min(R, n_razze * a / math.pi)
    e = math.exp(-2.0 * r_m ** 2 / w0 ** 2)
    eR = math.exp(-2.0 * R ** 2 / w0 ** 2)
    if R <= r_m:
        return 1.0 - eR
    return (1.0 - e) + (n_razze * a / math.pi) * (e - eR)


@lru_cache(maxsize=8192)
def eta_int(R, forma="sfera", w0=W0):
    """Frazione di potenza del fascio che investe il bersaglio.

    E' l'integrale di sovrapposizione fra il profilo gaussiano e la sagoma
    proiettata: NON si assume uguale a 1.  Per una sfera o un disco centrato
    e' esatto, 1 - exp(-2R^2/w0^2).  Per un anello sottile il fascio passa
    nel foro e eta_int -> 0 per R >> w0: e' il motivo per cui il reticolo
    non vince.
    """
    g = geometria(forma, R)
    eta = 0.0
    for el in g["elementi"]:
        if el[0] == "disco":
            eta += 1.0 - math.exp(-2.0 * el[1] ** 2 / w0 ** 2)
        elif el[0] == "anello":
            r, a = el[1], el[2]
            eta += (8.0 * a * r / w0 ** 2) * math.exp(-2.0 * r ** 2 / w0 ** 2)
        elif el[0] == "disco_offset":
            eta += _eta_disco_offset(el[1], el[2], w0)
        elif el[0] == "razze":
            eta += _eta_razze(el[1], el[2], el[3], w0)
        else:
            raise ValueError(el[0])
    return min(eta, 1.0)


# ---------------------------------------------------------------------------
# 5. FUNZIONE PUBBLICA 2: coppia
# ---------------------------------------------------------------------------
def delta_n(materiale, forma):
    """Birifrangenza efficace: intrinseca (cristallo) + di forma."""
    mat = MATERIALI[materiale]
    return mat["dn"] + DN_FORMA.get((materiale, forma), DN_FORMA_DEFAULT)


def delta_ritardo(R, materiale, forma, lam=LAMBDA):
    """Ritardo efficace (mediana) per una particella di raggio R.

    Sopra R_RAYLEIGH e' il ritardo di una lamina di spessore percorso:
        delta = 2 pi dn percorso / lambda.
    Sotto R_RAYLEIGH si raccorda con delta_bulk * (R/R_RAYLEIGH)^2: e' un
    raccordo di taglia a basso R (polarizzabilita' ~ volume), non un confine
    di regime fisico.  Con R_RAYLEIGH = 50 nm e la mappa che parte da 0.1 um,
    il fattore non si applica mai sulla griglia; il bordo inferiore della
    finestra e' deciso dalla soglia di coppia.
    """
    dn = delta_n(materiale, forma)
    g = geometria(forma, R)
    delta = 2.0 * math.pi * dn * g["percorso"] / lam
    if R < R_RAYLEIGH:
        delta *= (R / R_RAYLEIGH) ** 2
    return delta


def _media_ritardo(delta, sigma=SIGMA_RITARDO, n=N_DISP):
    """Media di (1-cos delta) e |sin delta| su una dispersione larga.

    Il caso monodisperso produce oscillazioni di coppia con la taglia che non
    sono una previsione: sono il posto dove cade il campione.  Si media su una
    distribuzione gaussiana del ritardo con sigma_rel = 0.5 (dispersione
    efficace di taglia/orientamento).  Per delta >> 1 la media tende a 1 per
    (1-cos) e a 2/pi per |sin|; per delta << 1 riproduce il valore puntuale.
    La media di (1-cos) ha forma chiusa:
        <1-cos(delta(1+sigma z))> = 1 - cos(delta) exp(-sigma^2 delta^2 / 2).
    Il caso ideale delta = pi (lamina lambda/2) e' riportato a parte: non e'
    quello che la mappa usa.
    """
    if delta <= 0.0:
        return 0.0, 0.0
    c1 = 1.0 - math.cos(delta) * math.exp(-0.5 * (sigma * delta) ** 2)
    if sigma * delta > 8.0:
        return c1, 2.0 / math.pi
    x, w = np.polynomial.hermite_e.hermegauss(n)
    w = w / math.sqrt(2.0 * math.pi)
    c2 = 0.0
    for xi, wi in zip(x, w):
        c2 += wi * abs(math.sin(delta * math.exp(sigma * xi)))
    return c1, c2


def coppia(canale, P, R, materiale, forma="sfera", *, l=0.0, s=1.0,
           S1=0.0, S2=0.0, S3=1.0, w0=W0, lam=LAMBDA):
    """Coppia sul bersaglio per uno dei tre canali.

    Parametri
    ---------
    canale : "assorbimento" | "ritardo" | "scattering"
    P      : potenza incidente sul bersaglio (W), prima dell'intercettazione
    l      : momento angolare orbitale per fotone (fascio LG)
    s      : +-1, momento angolare di spin per fotone circolare
    S1,S2,S3 : parametri di Stokes normalizzati del fascio in ingresso
               (circolare: S3=+-1; lineare: S1^2+S2^2=1)

    Ritorna un dizionario con tau, le potenze in gioco e il ritardo delta.

    Polarizzazione mista: la coppia e' la somma dei contributi di spin e
    lineare, non la somma in quadratura.  Per un orientamento ottimale
        tau/(P/omega) = |S3| (1-cos delta) + |sin delta| sqrt(S1^2+S2^2),
    perche' il termine di spin non dipende dall'orientamento e il termine
    lineare si somma ad esso.  La mappa usa luce circolare (S3=+-1), dove
    conta solo il primo termine.
    """
    mat = MATERIALI[materiale]
    g = geometria(forma, R)
    eta = eta_int(R, forma, w0)
    P_int = P * eta
    I0 = 2.0 * P / (math.pi * w0 ** 2)

    q = q_ottico(materiale, g["Req"])
    # potenza assorbita: il minore fra il limite di Mie (illuminazione
    # uniforme, valido per R << w0) e il tetto geometrico (tutta la potenza
    # intercettata, valido per R >> w0)
    P_abs_mie = I0 * g["A_proj"] * max(q["Q_abs"], 0.0)
    P_ass = min(P_abs_mie, P_int * alpha_assorbimento(materiale, g["percorso"]))

    out = {"canale": canale, "eta_int": eta, "P_int": P_int, "P_ass": P_ass,
           "tau": 0.0, "P_lavoro": 0.0, "delta": 0.0, "Q_abs": q["Q_abs"],
           "Q_sca": q["Q_sca"], "g_aniso": q["g"], "applicabile": True,
           "canale_utile": True}

    # 1. Lo scattering su una forma assialmente simmetrica non trasferisce
    #    J_z: la coppia e' nulla per simmetria, non piccola.  Il canale e'
    #    lecito solo per il manubrio (forma non assiale).
    if canale == "scattering" and forma != "manubrio":
        out["applicabile"] = False
        out["tau"] = 0.0
        out["motivo"] = ("forma assialmente simmetrica: J_z si conserva, "
                         "nessuna coppia di scattering (Marston & Crichton 1984)")
        return out

    if canale == "assorbimento":
        out["tau"] = (l + s) * P_ass / OMEGA_OTT
        out["P_lavoro"] = P_ass
        # Un dielettrico trasparente assorbe, ma non abbastanza da essere un
        # canale di coppia: entra solo come perdita termica.
        if mat["k"] <= 0.05:
            out["canale_utile"] = False
            out["nota_utile"] = "assorbimento volumetrico trascurabile"

    elif canale == "ritardo":
        dn = delta_n(materiale, forma)
        if mat["k"] > 0.05 or dn <= 0.0:
            # opaco o isotropo: il canale non esiste
            out["applicabile"] = False
            out["motivo"] = ("materiale opaco (k = %.2g)" % mat["k"]) if dn <= 0.0 \
                else "nessuna anisotropia"
            return out
        delta_bulk = 2.0 * math.pi * dn * g["percorso"] / lam
        delta = delta_ritardo(R, materiale, forma, lam)
        c_circ, c_lin = _media_ritardo(delta)
        tau = (P_int / OMEGA_OTT) * (c_circ * abs(S3) +
                                     c_lin * math.hypot(S1, S2))
        out["tau"] = tau
        out["P_lavoro"] = P_int
        out["delta"] = delta
        out["delta_bulk"] = delta_bulk
        out["c_circ"] = c_circ
        out["c_lin"] = c_lin
        out["delta_n"] = dn

    elif canale == "scattering":
        # Solo manubrio (vedi sopra): la formula e' un limite superiore, con
        # anisotropia di forma esplicita.  Non si applica a sfere o dischi.
        P_sca_mie = I0 * g["A_proj"] * max(q["Q_sca"], 0.0)
        P_sca = min(P_sca_mie, P_int)
        out["tau"] = (P_sca / OMEGA_OTT) * max(1.0 - q["g"], 0.0) * (l + s)
        out["P_lavoro"] = P_sca
        out["P_sca"] = P_sca
        out["nota"] = ("forma non assiale (manubrio): limite superiore, "
                       "non una soluzione di Mie del manubrio")

    else:
        raise ValueError("canale sconosciuto: %s" % canale)

    return out


# ---------------------------------------------------------------------------
# 6. FUNZIONE PUBBLICA 3: temperatura
# ---------------------------------------------------------------------------
def temperatura(P_ass, R, materiale, forma="sfera", mezzo="acqua", T0=T0):
    """Delta T nel fluido, oppure T di equilibrio radiativo in vuoto.

    Acqua:  dT = P_ass / (4 pi kappa R_eq)      (sfera isoterma, Nusselt = 1)
    Vuoto:  P_ass = eps_eff(T) sigma A_rad (T^4 - T0^4).
            eps_eff e' da Mie per i materiali tabulati e viene ricalcolata a
            ogni passo con la T corrente; per il Cu nero e' il parametro
            dichiarato 0.30 (il rivestimento non e' una sfera di Mie).
    """
    g = geometria(forma, R)
    if mezzo == "acqua":
        dT = P_ass / (4.0 * math.pi * K_ACQUA * g["Req"])
        return {"dT": dT, "T": T0 + dT, "eps_eff": None, "eps_fonte": None,
                "modello": "conduzione sferica, Nusselt = 1"}
    A = g["A_rad"]
    fonte = "dichiarato (non Mie)" if materiale == "cu_nero" else "Mie"
    T = T0
    e = eps_eff(materiale, g["Req"], T)
    for _ in range(6):
        e = eps_eff(materiale, g["Req"], T)
        T4 = T0 ** 4 + P_ass / max(e * SIGMA_SB * A, 1e-300)
        T = T4 ** 0.25
    return {"dT": T - T0, "T": T, "eps_eff": e, "eps_fonte": fonte,
            "modello": "equilibrio radiativo, emissivita' da Mie" if fonte == "Mie"
            else "equilibrio radiativo, eps_IR dichiarato"}


# ---------------------------------------------------------------------------
# 7. ATTrito ROTAZIONALE
# ---------------------------------------------------------------------------
def beta_acqua(R, forma="sfera"):
    """Coefficiente di attrito rotazionale in acqua a 310 K."""
    if forma == "sfera":
        return 8.0 * math.pi * MU_ACQUA * R ** 3
    if forma == "disco":
        return (32.0 / 3.0) * MU_ACQUA * R ** 3
    if forma == "anello":
        a = 0.1 * R
        return 8.0 * math.pi ** 2 * MU_ACQUA * R ** 3 / (2.0 * math.log(2.0 * R / a) + 1.0)
    if forma == "manubrio":
        r, d = 0.25 * R, 0.75 * R
        return 3.0 * math.pi * MU_ACQUA * r * (2.0 * d) ** 2 + 16.0 * math.pi * MU_ACQUA * r ** 3
    # reticolo: limite inferiore da corpo snello (anello piu' esterno)
    a = 0.03 * R
    return 8.0 * math.pi ** 2 * MU_ACQUA * R ** 3 / (2.0 * math.log(2.0 * R / a) + 1.0)


def beta_gas(R, forma="sfera", p_mbar=1e-6, T=T0):
    """Attrito molecolare libero (Epstein) su una sfera.

    beta = (16/3) p R^4 / v_bar, con v_bar = sqrt(8 k T / (pi m)).
    Il prefattore 16/3 vale per una sfera con accomodamento tangenziale
    unitario (sigma_t = 1).  Per forme non sferiche si usa lo stesso valore:
    non si inventano fattori di forma.  La cella marca il risultato con
    beta_shape_incerto e le affermazioni quantitative in vuoto valgono solo
    per le sfere.
    """
    p = p_mbar * P_MBAR
    v_bar = math.sqrt(8.0 * C.k * T / (math.pi * M_GAS))
    return (16.0 / 3.0) * p * R ** 4 / v_bar


def beta(R, forma, mezzo, p_mbar=None):
    if mezzo == "acqua":
        return beta_acqua(R, forma)
    return beta_gas(R, forma, p_mbar if p_mbar is not None else 1e-6)


# ---------------------------------------------------------------------------
# 8. MECCANICA
# ---------------------------------------------------------------------------
def omega_rottura(R, materiale, forma="sfera"):
    """sigma_max = c rho v^2  ->  v_lim = sqrt(sy/(c rho)),  omega = v_lim/R.

    Disco pieno: c = 0.40 (3+nu)/8, non 1 come l'anello.  Il tetto meccanico
    del disco e' 1.58 volte quello dell'anello, ma il momento d'inerzia e'
    diverso: omega_rottura cambia per entrambi i motivi, non solo per v_lim.
    """
    mat = MATERIALI[materiale]
    g = geometria(forma, R)
    v_lim = math.sqrt(mat["sy"] / (g["c_mecc"] * mat["rho"]))
    return v_lim / R, v_lim


# ---------------------------------------------------------------------------
# 9. CRITERI (scritti prima di calcolare)
# ---------------------------------------------------------------------------
CRITERI = {
    "Pe_r_min": 10.0,        # il rotore e' orientato, non rumoroso
    "dT_max_acqua": 20.0,    # K
    "T_frac_vuoto": 0.7,     # T <= 0.7 T_lim
    "omega_frac_rottura": 0.5,
    "P_max": 10e-3,          # W sul bersaglio
    "eta_mecc_min": 1e-7,    # criterio di rotore, non di fisica
    "t_spin_min_vuoto": 1.0, # s, soglia operativa dichiarata (volano)
}


def cella(canale, P, R, materiale, forma, mezzo, p_mbar=None, **kw):
    """Una cella della mappa: coppia, rotazione, termica, criteri.

    La validita' geometrica (R <= 2 um) e' un requisito, non un criterio
    termico o meccanico: fuori dominio del fascio la cella non e' valida
    nemmeno se tutti gli altri numeri sono buoni.
    """
    mat = MATERIALI[materiale]
    g = geometria(forma, R)
    co = coppia(canale, P, R, materiale, forma, **kw)
    th = temperatura(co["P_ass"], R, materiale, forma, mezzo)
    b = beta(R, forma, mezzo, p_mbar)
    I = g["kI"] * mat["rho"] * g["V"] * R ** 2
    w_lim, v_lim = omega_rottura(R, materiale, forma)

    valido_geometria = (R <= R_MAX_VALID)
    beta_shape_incerto = (mezzo == "vuoto" and forma != "sfera")

    tau = co["tau"]
    Pe = tau / (C.k * T_BATH)
    w = tau / b
    # Tempo per raggiungere meta' del limite meccanico a coppia costante.
    # Il limite inferiore assoluto e' quello con la coppia massima 2P/omega:
    # serve a distinguere "il materiale limita" da "la coppia non basta".
    t_spin = (I * 0.5 * w_lim / tau) if tau > 0 else float("inf")
    tau_max = 2.0 * P / OMEGA_OTT
    t_spin_limite = (I * 0.5 * w_lim / tau_max) if tau_max > 0 else float("inf")
    eta_mecc = w / OMEGA_OTT

    # Rendimento FISICO del rotore: coppia disponibile per la massima
    # velocita' meccanica, diviso la potenza incidente.  Non e' eta_mecc,
    # che confronta omega con omega_ottico ed e' solo una comodita' numerica.
    eta_lim = (tau * w_lim / P) if P > 0 else 0.0
    omega_max = min(w, w_lim)
    P_mech = tau * omega_max
    eta_phys = (P_mech / P) if P > 0 else 0.0

    # Bilancio trappola: forza di scattering (limite 2P/c, qui P/c) contro
    # il peso.  In vuoto la trappola e' necessaria; il rapporto dice quanto.
    F_scat = P / C.c
    peso = mat["rho"] * g["V"] * 9.81
    trap_ratio = (F_scat / peso) if peso > 0 else float("inf")
    trappola_necessaria = (mezzo == "vuoto")

    motivi = []
    if not valido_geometria:
        motivi.append("R > 2 um: fuori dominio del fascio (z_R ~ 2 um)")
    if not co["applicabile"]:
        motivi.append(co.get("motivo", "canale non applicabile"))
    if not co.get("canale_utile", True):
        motivi.append(co.get("nota_utile", "canale non utile"))
    if Pe < CRITERI["Pe_r_min"]:
        motivi.append("Pe_r < 10")
    if mezzo == "acqua":
        if th["dT"] > CRITERI["dT_max_acqua"]:
            motivi.append("dT > 20 K")
        if w > CRITERI["omega_frac_rottura"] * w_lim:
            motivi.append("omega > 0.5 omega_rottura")
    else:
        if th["T"] > CRITERI["T_frac_vuoto"] * mat["Tm"]:
            motivi.append("T > 0.7 T_lim")
        if t_spin < CRITERI["t_spin_min_vuoto"]:
            motivi.append("t_spin < 1 s (rottura prima dell'uso)")
    if P > CRITERI["P_max"]:
        motivi.append("P > 10 mW")

    valido = (len(motivi) == 0)
    return {
        "canale": canale, "materiale": materiale, "forma": forma, "mezzo": mezzo,
        "p_mbar": p_mbar, "R": R, "P": P,
        "eta_int": co["eta_int"], "P_ass": co["P_ass"], "P_lavoro": co["P_lavoro"],
        "tau": tau, "Pe_r": Pe, "beta": b, "omega": w, "omega_rottura": w_lim,
        "omega_frac": w / w_lim if w_lim > 0 else float("inf"),
        "t_spin": t_spin, "t_spin_limite": t_spin_limite, "tau_max": tau_max,
        "eta_mecc": eta_mecc, "canale_utile": co.get("canale_utile", True),
        "dT": th["dT"], "T": th["T"], "eps_eff": th["eps_eff"],
        "eps_fonte": th.get("eps_fonte"),
        "I": I, "valido": valido, "valido_rotore": valido and eta_mecc >= CRITERI["eta_mecc_min"],
        "valido_rotore_fisico": valido and eta_lim >= 1e-7,
        "valido_fisico": valido and eta_lim >= 1e-7,
        "valido_geometria": valido_geometria,
        "beta_shape_incerto": beta_shape_incerto,
        "eta_lim": eta_lim, "omega_max": omega_max, "P_mech": P_mech,
        "eta_phys": eta_phys, "F_scat": F_scat, "peso": peso,
        "trap_ratio": trap_ratio, "trappola_necessaria": trappola_necessaria,
        "motivi": motivi,
    }


# ---------------------------------------------------------------------------
# 10. GRIGLIA
# ---------------------------------------------------------------------------
R_GRID = tuple(float(x) for x in np.logspace(math.log10(0.1e-6), math.log10(30e-6), 18))
P_GRID = tuple(float(x) for x in np.logspace(math.log10(1e-6), math.log10(100e-3), 15))
MAT_GRID = ("cu_lucido", "cu_nero", "silice", "vaterite")
CANALI = ("assorbimento", "ritardo", "scattering")
MEZZI = (("acqua", None), ("vuoto", 1e-6), ("vuoto", 1e-8))


def mappa(scrittura=True):
    righe = []
    for materiale in MAT_GRID:
        for forma in FORME:
            for mezzo, p in MEZZI:
                for canale in CANALI:
                    for R in R_GRID:
                        for P in P_GRID:
                            c = cella(canale, P, R, materiale, forma, mezzo, p)
                            c["p"] = p
                            righe.append(c)
    if scrittura:
        scrivi_griglia(righe)
    return righe


def _p_per_pe(tau_per_w, pe_target=10.0):
    return pe_target * KT / tau_per_w if tau_per_w > 0 else float("inf")


def frontiere(righe=None):
    """Frontiera di ogni combinazione: P_min(Pe=10), P_max(termico/meccanico).

    Piu' utile della griglia grezza: dice subito dove la finestra esiste e
    quanto e' larga.
    """
    out = []
    for materiale in MAT_GRID:
        for forma in FORME:
            for mezzo, p in MEZZI:
                for canale in CANALI:
                    base = cella(canale, 1.0, 1e-6, materiale, forma, mezzo, p)
                    if not base["tau"] or not base.get("canale_utile", True):
                        nota = base["motivi"][0] if base["motivi"] else "canale non utile"
                        out.append({"materiale": materiale, "forma": forma,
                                    "mezzo": mezzo, "p_mbar": p, "canale": canale,
                                    "esiste": False, "nota": nota})
                        continue
                    tao_w = base["tau"] / base["P"]        # coppia per watt
                    p10 = _p_per_pe(tao_w)
                    righe_R = []
                    for R in R_GRID:
                        c1 = cella(canale, 1e-3, R, materiale, forma, mezzo, p)
                        if not c1["tau"] or not c1.get("canale_utile", True):
                            continue
                        Rcell = 1e-3
                        tao_R = c1["tau"] / Rcell
                        p10_R = _p_per_pe(tao_R)
                        # tetto termico
                        if mezzo == "acqua":
                            dT1 = c1["dT"]
                            p_term = (CRITERI["dT_max_acqua"] / dT1) * Rcell if dT1 > 0 else float("inf")
                            # tetto meccanico: omega = tau/beta lineare in P
                            b = c1["beta"]
                            w_lim = c1["omega_rottura"]
                            p_mecc = (CRITERI["omega_frac_rottura"] * w_lim * b / tao_R) if tao_R > 0 else float("inf")
                        else:
                            p_term = _p_termico_vuoto(c1, R, materiale, forma)
                            w_lim = c1["omega_rottura"]
                            I = c1["I"]
                            p_mecc = (I * 0.5 * w_lim / CRITERI["t_spin_min_vuoto"]) / tao_R \
                                if tao_R > 0 else float("inf")
                        p_max_fisico = min(p_term, p_mecc)
                        righe_R.append({"R": R, "p10": p10_R, "p_term": p_term,
                                        "p_mecc": p_mecc,
                                        "p_max_fisico": p_max_fisico,
                                        "p_max": min(p_max_fisico, CRITERI["P_max"]),
                                        "eta_mecc_a_p10": (p10_R * tao_R / c1["beta"]) / OMEGA_OTT
                                        if c1["beta"] > 0 else 0.0,
                                        "eta_lim_a_p10": tao_R * w_lim,
                                        "valido_geometria": c1["valido_geometria"]})
                    finestre = []
                    for r in righe_R:
                        if not r["valido_geometria"]:
                            continue
                        lo, hi = r["p10"], r["p_max"]
                        if lo <= hi and lo <= CRITERI["P_max"]:
                            finestre.append({"R": r["R"], "P_lo": lo, "P_hi": hi,
                                             "larghezza_decadi": math.log10(hi / lo) if lo > 0 else 0.0})
                    out.append({
                        "materiale": materiale, "forma": forma, "mezzo": mezzo,
                        "p_mbar": p, "canale": canale, "esiste": True,
                        "P_pe10_a_R1um": p10, "celle": righe_R,
                        "finestre": finestre,
                        "ha_finestra": len(finestre) > 0,
                        "R_min_finestra": min([f["R"] for f in finestre]) if finestre else None,
                        "R_max_finestra": max([f["R"] for f in finestre]) if finestre else None,
                    })
    return out


def _p_termico_vuoto(c, R, materiale, forma):
    """Potenza assorbita che porta a 0.7 T_lim, poi convertita in P incidente.

    L'emissivita' e' valutata alla temperatura del criterio (0.7 T_lim), non a
    quella del bagno: e' lo stesso principio della V8 per temperatura().
    """
    mat = MATERIALI[materiale]
    g = geometria(forma, R)
    T = CRITERI["T_frac_vuoto"] * mat["Tm"]
    e = eps_eff(materiale, g["Req"], T)
    P_ass = e * SIGMA_SB * g["A_rad"] * (T ** 4 - T0 ** 4)
    # P_ass e' lineare in P: dalla cella di riferimento si ricava il fattore
    if c["P_ass"] <= 0:
        return float("inf")
    return P_ass * (c["P"] / c["P_ass"])


def scrivi_griglia(righe):
    path = os.path.join(DATI, "M1_griglia.csv")
    campi = ("canale", "materiale", "forma", "mezzo", "p", "R", "P", "eta_int",
             "P_ass", "tau", "Pe_r", "beta", "omega", "omega_rottura",
             "omega_frac", "t_spin", "t_spin_limite", "eta_mecc", "eta_lim",
             "omega_max", "P_mech", "eta_phys", "dT", "T",
             "eps_eff", "eps_fonte", "canale_utile", "valido", "valido_rotore",
             "valido_rotore_fisico", "valido_fisico", "valido_geometria",
             "beta_shape_incerto", "F_scat", "peso", "trap_ratio",
             "trappola_necessaria")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(campi)
        for r in righe:
            w.writerow(["%.6g" % r[k] if isinstance(r[k], float) else r[k] for k in campi])
    return path


# ---------------------------------------------------------------------------
# 11. SELF-TEST
# ---------------------------------------------------------------------------
def _test():
    print("=== SELF-TEST V8.0 ===")
    print("-- Mie: sfera di rame a 532 nm")
    for R in (10e-9, 100e-9, 1e-6, 10e-6):
        q = q_ottico("cu_lucido", R)
        print("   R=%6.0f nm  x=%7.2f  Q_abs=%.4f  Q_sca=%.4f  g=%.3f  [%s]"
              % (R * 1e9, 2 * math.pi * R / LAMBDA, q["Q_abs"], q["Q_sca"], q["g"], q["regime"]))
    print("-- eta_int (integrale di sovrapposizione, w0 = %.0f nm)" % (W0 * 1e9))
    for forma in FORME:
        vals = [eta_int(R, forma) for R in (0.1e-6, 1e-6, 10e-6, 30e-6)]
        print("   %-9s " % forma + "  ".join("%.3e" % v for v in vals))
    print("-- emissivita' efficace (Mie, tranne Cu nero dichiarato)")
    for mat in ("cu_lucido", "cu_nero", "silice", "vaterite", "diamante"):
        print("   %-10s " % mat + "  ".join("%.4f" % eps_eff(mat, R)
                                           for R in (0.5e-6, 3e-6, 15e-6, 30e-6))
              + ("  [dichiarato]" if mat == "cu_nero" else ""))
    print("-- coppia a ritardo di fase: casi limite")
    for (S1, S2, S3, nome) in ((0, 0, 1, "circolare"), (1, 0, 0, "lineare")):
        for dn_frac, etichetta in ((0.5, "lambda/4"), (1.0, "lambda/2")):
            # costruisco un delta noto con un finto materiale
            print("   %-10s %-8s tau/(P/omega) = %.4f"
                  % (nome, etichetta, _tau_norm(S1, S2, S3, dn_frac * math.pi)))
    print("-- meccanica: omega_rottura")
    for mat in ("cu_lucido", "silice", "vaterite"):
        for forma in ("sfera", "disco", "anello"):
            w, v = omega_rottura(1e-6, mat, forma)
            print("   %-9s %-7s v_lim=%8.1f m/s  omega_lim=%10.3e rad/s  (%.2f GHz)"
                  % (mat, forma, v, w, w / 2 / math.pi / 1e9))
    print("-- scattering: solo forme non assiali")
    for forma in ("sfera", "disco", "manubrio"):
        c = coppia("scattering", 10e-3, 1e-6, "silice", forma)
        print("   %-9s applicabile=%s  tau=%.2e N m  %s"
              % (forma, c["applicabile"], c["tau"], c.get("motivo", "")))
    print("-- inerzia del manubrio (R = 100 nm)")
    g = geometria("manubrio", 0.1e-6)
    I = g["kI"] * MATERIALI["silice"]["rho"] * g["V"] * (0.1e-6) ** 2
    print("   kI = %.4f (adimensionale)  I = %.3e kg m^2  [atteso ~1.7e-33]"
          % (g["kI"], I))
    print("-- dispersione larga: P(Pe=10), vaterite sfera, acqua")
    for R in (0.1e-6, 0.5e-6, 1e-6, 2e-6, 10e-6):
        c = cella("ritardo", 1e-3, R, "vaterite", "sfera", "acqua")
        p10 = 10.0 * KT / (c["tau"] / 1e-3) if c["tau"] > 0 else float("inf")
        print("   R=%6.2f um  delta=%6.3f  c_circ=%.3f  P10=%.1f uW"
              % (R * 1e6, coppia("ritardo", 1e-3, R, "vaterite", "sfera")["delta"],
                 coppia("ritardo", 1e-3, R, "vaterite", "sfera")["c_circ"],
                 p10 * 1e6))


def _tau_norm(S1, S2, S3, delta):
    """Coppia normalizzata a P/omega per un ritardo noto (senza dispersione).

    Massimo sull'orientamento: il termine di spin non dipende dall'angolo, il
    termine lineare si somma.  Per gli stati puri (S3=1 o S1=1) coincide con
    la formula usata nei test A e B.
    """
    c1 = 1.0 - math.cos(delta)
    c2 = abs(math.sin(delta))
    return c1 * abs(S3) + c2 * math.hypot(S1, S2)


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        _test()
        return
    _test()
    print("\n=== MAPPA ===")
    righe = mappa()
    n_tot = len(righe)
    n_val = sum(1 for r in righe if r["valido"])
    n_val_acqua = sum(1 for r in righe if r["valido"] and r["mezzo"] == "acqua")
    n_val_vuoto = sum(1 for r in righe if r["valido"] and r["mezzo"] == "vuoto")
    n_rot = sum(1 for r in righe if r["valido_rotore"])
    n_rot_fis = sum(1 for r in righe if r["valido_rotore_fisico"])
    n_eta_lim = sum(1 for r in righe if r["eta_lim"] >= 1e-7)
    n_le2 = sum(1 for r in righe if r["R"] <= R_MAX_VALID)
    n_gt2 = n_tot - n_le2
    print("   celle: %d" % n_tot)
    print("   valide: %d  (acqua %d, vuoto %d)" % (n_val, n_val_acqua, n_val_vuoto))
    print("   rotori formali (valido_rotore): %d" % n_rot)
    print("   rotori fisici (eta_lim >= 1e-7): %d  [celle con eta_lim>=1e-7: %d]"
          % (n_rot_fis, n_eta_lim))
    print("   celle con R <= 2 um: %d ; R > 2 um: %d" % (n_le2, n_gt2))
    fr = frontiere()
    valide = [f for f in fr if f.get("ha_finestra")]
    n_fin = sum(len(f["finestre"]) for f in valide)
    print("   combinazioni con finestra: %d / %d  (finestre totali: %d)"
          % (len(valide), len(fr), n_fin))

    riassunto = []
    for f in valide:
        riassunto.append({
            "canale": f["canale"], "materiale": f["materiale"], "forma": f["forma"],
            "mezzo": f["mezzo"], "p_mbar": f["p_mbar"],
            "R_min_um": f["R_min_finestra"] * 1e6 if f["R_min_finestra"] else None,
            "R_max_um": f["R_max_finestra"] * 1e6 if f["R_max_finestra"] else None,
            "P_pe10_a_R1um_W": f["P_pe10_a_R1um"],
            "finestre": f["finestre"],
        })
    riassunto.sort(key=lambda x: -len(x["finestre"]))
    for r in riassunto[:20]:
        print("   %-12s %-10s %-9s %-5s  R = %5.1f-%5.1f um  finestre=%2d  P(Pe=10,1um)=%.2e W"
              % (r["canale"], r["materiale"], r["forma"], r["mezzo"],
                 r["R_min_um"] or 0, r["R_max_um"] or 0, len(r["finestre"]),
                 r["P_pe10_a_R1um_W"]))

    payload = {
        "versione": "V8.0",
        "nota": "mappa della finestra valida: tre canali, cinque forme, quattro "
                "materiali, tre mezzi.  Ogni numero viene da v8_mappa.py.  "
                "V8.0 e' la correzione tecnica della V7.1: vaterite con dn "
                "efficace 0.06 e T_m 723 K, validita' del fascio R <= 2 um, "
                "attrito di Epstein (16/3) per sfere, emissivita' dipendente "
                "da T, Pe_r con T_bath esplicita, rendimento fisico eta_lim, "
                "flag di trappola.",
        "correzioni": [
            "vaterite: Δn = 0.06 efficace (range 0.06-0.10), non 0.17 della "
            "calcite; k = 1e-6; T_m = 723 K (vaterite -> calcite)",
            "validita' del fascio: R <= 2 µm (z_R ~ 2 µm); oltre la cella e' "
            "invalida per geometria",
            "attrito molecolare di Epstein per sfera: beta = (16/3) p R^4 / "
            "v_bar; forme non sferiche marcate beta_shape_incerto",
            "emissivita' eps_eff(materiale, R, T) con peso di Planck a T e "
            "iterazione in temperatura; Cu nero 0.30 dichiarato e costante",
            "Pe_r = tau / (k_B T_bath), T_bath = 310 K (definizione esplicita)",
            "rendimento fisico eta_lim = tau omega_rottura / P; "
            "valido_rotore_fisico distinto da valido_rotore",
            "bilancio trappola: F_scat = P/c, peso = rho V g, trap_ratio",
            "R_RAYLEIGH = 50 nm come raccordo di taglia (inattivo sulla "
            "griglia, che parte da 0.1 um)",
            "tabelle/conteggi/unita': celle valide distinte per mezzo, rotori "
            "formali e fisici, R <= 2 um vs R > 2 um, P_MBAR = 100 Pa/mbar",
        ],
        "costanti": {"lambda_nm": LAMBDA * 1e9, "omega_ott": OMEGA_OTT,
                     "nu_ott_Hz": NU_OTT, "w0_nm": W0 * 1e9, "T0_K": T0,
                     "T_bath_K": T_BATH, "mu_acqua_Pas": MU_ACQUA,
                     "R_rayleigh_nm": R_RAYLEIGH * 1e9,
                     "R_max_valid_um": R_MAX_VALID * 1e6,
                     "sigma_ritardo": SIGMA_RITARDO},
        "criteri": CRITERI,
        "griglia": {"R_um": [r * 1e6 for r in R_GRID],
                    "P_W": list(P_GRID),
                    "materiali": list(MAT_GRID), "forme": list(FORME),
                    "canali": list(CANALI),
                    "mezzi": [{"mezzo": m, "p_mbar": p} for m, p in MEZZI]},
        "conteggi": {"celle_totali": n_tot, "celle_valide": n_val,
                     "valide_acqua": n_val_acqua, "valide_vuoto": n_val_vuoto,
                     "combinazioni_con_finestra": len(valide),
                     "combinazioni_totali": len(fr),
                     "finestre_totali": n_fin,
                     "rotori_formali": n_rot,
                     "rotori_fisici_eta_lim": n_rot_fis,
                     "celle_eta_lim_ge_1e-7": n_eta_lim,
                     "celle_R_le_2um": n_le2, "celle_R_gt_2um": n_gt2},
        "frontiere": fr,
        "riassunto": riassunto,
    }
    with open(os.path.join(DATI, "M1_mappa.json"), "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1, ensure_ascii=False, default=float)
    print("   scritto dati_v8/M1_mappa.json e dati_v8/M1_griglia.csv")


if __name__ == "__main__":
    main()
