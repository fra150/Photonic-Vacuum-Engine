"""
MOTORE FOTONICO DEL VUOTO - V7.3
v73_contrappeso.py - il ciclo di 720 gradi corretto: un contatore Z2 nel sangue.

Questa versione non aggiunge un'idea: toglie gli errori della V7.2 e rifà il
conto.  Le correzioni del ciclo sono elencate in CORREZIONI e finiscono
nel JSON; le sette correzioni V7.3 (cancello, termica, rumore esatto) sono
negli script v73_taglio.py, v73_termica.py e v73_rumore.py.

CORREZIONI RISPETTO ALLA V7.2
-----------------------------
1. CALORE DI CAMMINO.  La sorgente e' l'assorbimento, non la potenza tolta al
   fascio.  La V7.2 usava P_inc - P_target, cioe' contava come calore anche la
   luce diffusa.  Con g ~ 0.98 la luce diffusa e' quasi tutta in avanti e non
   viene assorbita nei primi micron.  Qui:
       dT_path = [mu_a P_inc / (2 pi kappa)] *
                 [(1 - exp(-mu_t z))/(mu_t z)] * ln(1 + z/w0)
   con mu_t = mu_a + mu_s.  Il fattore (1-exp(-mu_t z))/(mu_t z) e' la media
   del deposito lungo il cammino; ln(1+z/w0) e' il taglio cilindrico.
   mu_a(532) = 2.35e4 1/m (HbO2, 150 g/L, OMLC), non 7500.

2. LAVORO VISCOSO.  A Pe_r fissato la coppia e' tau = Pe_r kT.  Il lavoro per
   720 gradi e' tau * 4 pi = 5.38e-19 J = 126 kT, indipendente da beta.  Il
   fattore 4.88 della V7.2 e' un rapporto di TEMPI (e di energia ottica
   fornita), non di lavoro meccanico.  L'efficienza vera e' ~1e-16.

3. RUMORE E PROTOCOLLO.  La varianza per giro e' 4 pi / Pe_r in entrambi gli
   stati, quindi sigma_theta(360) = sqrt(4 pi / Pe_r) = 1.12 rad = 64 gradi,
   non 80.  Le durate dei giri alternati differiscono di beta_cw/beta_no =
   3.88, con CV = 1/sqrt(pi Pe_r) = 18%.  Il test sulla durata e' appaiato
   (3 giri veloci contro 3 lenti); la statistica non migliora con N, perche'
   la prima passa esatta e' una inverse Gaussian e la sua shape scala con N.
   I 119-595 cicli della V7.2 sono un artefatto della FFT cieca.

4. INPUT.  Viscosita' del plasma 1.3 mPa s nel cell-free layer (CFL), sangue
   intero 4 mPa s fuori; interpolazione volumetrica con delta_CFL = 3 um.
   Delta_n della vaterite 0.08 (0.06-0.10), non 0.172 della calcite.
   R_m scala come 1/Delta_n; P_target a piccolo R come 1/Delta_n^2.

5. INCOERENZE.  Il predicato di validita' include tau > tau_B (Bingham).
   z e' la profondita' della superficie prossimale del rotore sotto
   l'interfaccia sangue/vetro.  Lo spot secondario da 1 mW non raddoppia la
   potenza: con 0.2 mW di fascio primario il totale e' 1.2 mW, cioe' ~6x.
   Il transitorio non e' "sei ordini" sotto T_360: il rapporto va calcolato.
   Il titolo di Nature 394, 348 (1998) e' "Optical alignment and spinning of
   laser-trapped microscopic particles".  Dirac §11 non contiene il belt
   trick: si citano Nakahara e Misner-Thorne-Wheeler, e il belt trick resta
   un'analogia.

6. TOPOLOGIA.  Il periodo di 720 gradi viene da un contatore Z2 (il latch
   A/B), non dalla copertura doppia SO(3) -> SU(2).  La copertura doppia vale
   per un corpo legato all'ambiente; qui il linker collega due parti dello
   stesso rotore.  La simmetria SO(3) -> SU(2) e' un'analogia, non il
   meccanismo.

Uso:
    python codice/v73_contrappeso.py
"""

import json
import math
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(BASE, ".."))
if BASE not in sys.path:
    sys.path.insert(0, BASE)

import v7_mappa as m  # noqa: E402

DATI = os.path.join(ROOT, "dati")
PAPER = os.path.join(ROOT, "paper")
TAB = os.path.join(PAPER, "it", "tabelle")
TAB_EN = os.path.join(PAPER, "en", "tabelle")
FIG = os.path.join(PAPER, "it", "figures")
FIG_EN = os.path.join(PAPER, "en", "figures")
FIG_ROOT = os.path.join(ROOT, "figure")
for d in (DATI, TAB, TAB_EN, FIG, FIG_EN, FIG_ROOT,
          os.path.join(FIG_ROOT, "it"), os.path.join(FIG_ROOT, "en")):
    os.makedirs(d, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. COSTANTI E INPUT
# ---------------------------------------------------------------------------
T0 = 310.0
KT = 1.380649e-23 * T0
C_LUCE = 299792458.0
W0 = 0.5e-6

# Sangue: plasma nel CFL, sangue intero fuori.  kappa e rho sono valori medi.
SANGUE = {
    "mu_p": 1.3e-3,          # Pa s, plasma a 310 K
    "mu_b": 4.0e-3,          # Pa s, sangue intero ad alto shear
    "rho": 1060.0,           # kg/m3
    "kappa": 0.55,           # W/m/K, sangue intero
    "cp": 3600.0,            # J/kg/K
    "T": T0,
    "delta_cfl": 3.0e-6,     # m, spessore del cell-free layer
    "nota_mu": "plasma 1.3 mPa s nel CFL; sangue intero 4 mPa s fuori; "
               "interpolazione volumetrica; 7-10 mPa s restano come "
               "sensibilita' a basso shear.",
}

# Attenuazione del sangue intero a 150 g/L di Hb.
# mu_a: OMLC/Prahl, mu_a = 2.303 * eps * 150/64500 cm^-1 = 0.005356 eps.
# mu_s: valori efficaci con g ~ 0.98 (Bosschaart 2014); mu_s' = mu_s(1-g).
OTTICA = {
    532: {"mu_a": 2.35e4, "mu_s": 8.0e4, "g": 0.98,
          "eps": 43876.0,
          "fonte": "HbO2, eps(532)=43,876 cm^-1/M, 150 g/L (OMLC)"},
    800: {"mu_a": 4.37e2, "mu_s": 4.5e4, "g": 0.98,
          "eps": 816.0,
          "fonte": "HbO2, eps(800)=816 cm^-1/M, 150 g/L (OMLC)"},
    1064: {"mu_a": 2.00e2, "mu_s": 3.5e4, "g": 0.98,
           "eps": 374.0,
           "fonte": "coda HbO2 ~1000 nm + acqua; valore dichiarato"},
}

DN_VATERITE = 0.08
DN_RANGE = (0.06, 0.10)
PE_R = 10.0
SIGMA_Y = 0.01              # Pa, snervamento del sangue intero a basso shear
SIGMA_Y_RANGE = (0.005, 0.01)
D_T_OPERATIVA = 20.0        # K, 57 C: criterio operativo, non biologico
D_T_SICURA = 5.0            # K, 42 C: soglia prudente
SIGMA_DELTA = 0.10          # dispersione relativa di delta (10%)
LATCH_RATIO = None          # calcolato da beta_ratio()
P_MAX = 1.0                 # W, limite operativo sul fascio primario

R_GRID = tuple(float(x) for x in np.logspace(math.log10(0.15e-6),
                                             math.log10(20e-6), 36))
R_FINE = tuple(float(x) for x in np.logspace(math.log10(0.15e-6),
                                             math.log10(20e-6), 90))
Z_GRID = (2e-6, 5e-6, 10e-6, 15e-6, 20e-6)
LAM = (532, 800, 1064)

# ---------------------------------------------------------------------------
# 1b. LE SEI CORREZIONI, IN FORMA LEGGIBILE E MACCHINA
# ---------------------------------------------------------------------------
CORREZIONI = [
    {"voce": "calore di cammino",
     "v72": "P_path = P_inc - P_target, cioe' tutta la luce tolta al fascio "
            "diventa calore; mu_a(532) = 7500 1/m",
     "v73": "P_path = mu_a P_inc z con il fattore di media "
            "(1-e^(-mu_t z))/(mu_t z); mu_a(532) = 2.35e4 1/m",
     "perche": "con g ~ 0.98 la luce diffusa resta in avanti e non viene "
               "assorbita nei primi micron; contarla come calore gonfiava "
               "il riscaldamento e falsava il confronto fra lambda"},
    {"voce": "lavoro viscoso",
     "v72": "costo 4.88 = beta_tot/beta_rot presentato come costo "
            "energetico del contrappeso",
     "v73": "W_720 = tau 4 pi = 126 k_BT, indipendente da beta; il 4.88 e' "
            "un rapporto di durate (e di energia ottica), non di lavoro",
     "perche": "a Pe_r fissato la coppia e' tau = Pe_r k_BT; il lavoro "
               "meccanico di un giro e' tau per 2 pi.  L'efficienza vera "
               "e' ~1e-16"},
    {"voce": "rumore e protocollo",
     "v72": "sigma_theta(360) = 80 deg; N = 119-595 giri da una FFT cieca; "
            "contrasto 1.00 letto come prestazione",
     "v73": "sigma_theta(360) = sqrt(4 pi/Pe_r) = 1.12 rad = 64 deg; "
            "CV della durata = 1/sqrt(pi Pe_r) = 18%; il test appaiato 3+3 "
            "e' obbligatorio e la statistica non migliora con N",
     "perche": "la varianza per giro e' 4 pi/Pe_r in entrambi gli stati; "
               "il contrasto 1.00 e' la definizione del latch, non una "
               "misura; il direttore e' pi-periodico e serve un marcatore "
               "a 2 pi (la bead)"},
    {"voce": "input del sangue e della vaterite",
     "v72": "mu = 4 mPa s ovunque; Delta_n = 0.172 (calcite); nessun CFL",
     "v73": "mu_eff(R,z): plasma 1.3 mPa s nel CFL, sangue intero 4 mPa s "
            "fuori, frazione volumetrica con delta_CFL = 3 um; "
            "Delta_n = 0.08 (0.06-0.10)",
     "perche": "il rotore sta nel CFL; la vaterite non ha la birifrangenza "
               "della calcite.  R_m va come 1/Delta_n e P_target a piccolo R "
               "come 1/Delta_n^2"},
    {"voce": "incoerenze interne",
     "v72": "nessun tau > tau_B nella validita'; z non definito; spot "
            "secondario 'circa doppio'; transitorio 'sei ordini'; titolo "
            "Nature sbagliato; Dirac §11",
     "v73": "tau > tau_B nel predicato; z = profondita' della superficie "
            "prossimale; 1 mW + 0.2 mW = 1.2 mW (~6x); rapporto calcolato; "
            "titolo corretto; Nakahara/MTW",
     "perche": "un modello valido solo dove il sangue scorre come un "
               "liquido semplice non e' un modello del sangue; i numeri "
               "vanno confrontati con la stessa convenzione"},
    {"voce": "topologia",
     "v72": "il periodo di 720 gradi attribuito alla copertura doppia "
            "SO(3) -> SU(2)",
     "v73": "contatore Z2 classico (latch A/B); la copertura doppia resta "
            "un'analogia",
     "perche": "il linker collega due parti dello stesso rotore, non il "
               "rotore all'ambiente; il segnale del direttore e' "
               "pi-periodico e la bead e' il marcatore a 2 pi"},
]


# ---------------------------------------------------------------------------
# 2. IL SANGUE: VISCOSITA' LOCALE E GEOMETRIA DEL CFL
# ---------------------------------------------------------------------------
def frazione_cfl(R, z, delta=None):
    """Frazione del volume della sfera dentro il CFL (0 <= depth <= delta).

    z e' la profondita' della superficie PROSSIMALE del rotore sotto
    l'interfaccia sangue/vetro.  Il centro e' a z + R.  La fetta di sfera
    oltre il piano a distanza h dal centro ha volume
        V_cap = pi (R+h)^2 (2R-h)/3,
    quindi f = (2 + 3x - x^3)/4 con x = h/R.
    """
    delta = SANGUE["delta_cfl"] if delta is None else delta
    zc = z + R
    x = (delta - zc) / R
    if x >= 1.0:
        return 1.0
    if x <= -1.0:
        return 0.0
    return 0.25 * (2.0 + 3.0 * x - x ** 3)


def mu_eff(R, z, delta=None):
    """Viscosita' volumetrica media: plasma nel CFL, sangue intero fuori."""
    f = frazione_cfl(R, z, delta)
    return f * SANGUE["mu_p"] + (1.0 - f) * SANGUE["mu_b"]


def profondita_balistica(lam):
    d = OTTICA[lam]
    return 1.0 / (d["mu_a"] + d["mu_s"])


def trasmissione(lam, z):
    d = OTTICA[lam]
    return math.exp(-(d["mu_a"] + d["mu_s"]) * z)


# ---------------------------------------------------------------------------
# 3. COPPIA: RITARDO DI FASE
# ---------------------------------------------------------------------------
def delta_ritardo(R, lam_nm, dn=DN_VATERITE, path_factor=2.0):
    """Ritardo a singolo passaggio: delta = 2 pi dn (path_factor R) / lambda."""
    return 2.0 * math.pi * dn * path_factor * R / (lam_nm * 1e-9)


def c_circ(delta, sigma=SIGMA_DELTA):
    """Media gaussiana di 1 - cos(delta) su una dispersione relativa sigma.

    <1 - cos(delta(1+sigma z))> = 1 - cos(delta) exp(-sigma^2 delta^2 / 2).
    Per delta -> 0 riproduce delta^2/2; per sigma -> 0 riproduce 1-cos.
    """
    return 1.0 - math.cos(delta) * math.exp(-0.5 * (sigma * delta) ** 2)


def coppia_per_watt(R, lam_nm, dn=DN_VATERITE, sigma=SIGMA_DELTA,
                    path_factor=2.0, w0=W0):
    """Coppia per unita' di potenza.

    Nel limite di piccolo ritardo, tau/P = eta (1-cos delta)/omega_lambda
    con 1-cos delta ~ delta^2/2 ~ (4 pi dn R / lambda)^2 / 2 e
    omega_lambda = 2 pi c / lambda.  Quindi
        tau/P ~ eta * 8 pi dn^2 R^2 / (c lambda),
    cioe' tau/P va come 1/lambda e P_target, a Pe fissato, va come lambda
    (NON lambda^2: il fattore lambda^2 del ritardo e' diviso da 1/omega).
    Il controllo numerico e' in verifica_scaling().
    """
    eta = m.eta_int(R, "sfera", w0)
    delta = delta_ritardo(R, lam_nm, dn, path_factor)
    cc = c_circ(delta, sigma)
    omega_lam = 2.0 * math.pi * C_LUCE / (lam_nm * 1e-9)
    return {"eta": eta, "delta": delta, "c_circ": cc,
            "tau_per_W": eta * cc / omega_lam,
            "omega_lam": omega_lam}


def verifica_scaling(R=0.5e-6, tol=0.02):
    """Autotest delle due dipendenze che il testo usa.

    (a) tau/P va come 1/lambda, non 1/lambda^2.  A parita' di R, dn, eta_int
        e senza dispersione, il rapporto fra 1064 e 800 nm deve essere
        1064/800 = 1.33.
    (b) P_target a piccolo R va come 1/dn^2.  Con delta piccolo,
        (tau/P)(dn=0.10)/(tau/P)(dn=0.06) deve valere (0.10/0.06)^2.
    """
    a = coppia_per_watt(R, 800, sigma=0.0)["tau_per_W"]
    b = coppia_per_watt(R, 1064, sigma=0.0)["tau_per_W"]
    r_lam = b / a
    atteso_lam = 800.0 / 1064.0
    ok_lam = abs(r_lam / atteso_lam - 1.0) < tol

    R_dn = 0.2e-6
    c = coppia_per_watt(R_dn, 800, dn=0.06, sigma=0.0)["tau_per_W"]
    d = coppia_per_watt(R_dn, 800, dn=0.10, sigma=0.0)["tau_per_W"]
    r_dn = d / c
    atteso_dn = (0.10 / 0.06) ** 2
    ok_dn = abs(r_dn / atteso_dn - 1.0) < tol
    return {"lambda_1064_su_800": r_lam, "lambda_atteso": atteso_lam,
            "lambda_ok": ok_lam,
            "dn_010_su_006": r_dn, "dn_atteso": atteso_dn, "dn_ok": ok_dn}


# ---------------------------------------------------------------------------
# 4. CONTRAPPESO: INERZIA E ATTRITO
# ---------------------------------------------------------------------------
def beta_ratio(R=None, d_rel=2.0, frazione_volume=0.5):
    """beta_tot/beta_rot per una bead di volume V/2 a distanza d = 2 R.

    beta_cw = 8 pi mu r_cw^3 + 6 pi mu r_cw d^2;  beta_rot = 8 pi mu R^3.
    Il rapporto e' indipendente da R e da mu.
    """
    r_cw = frazione_volume ** (1.0 / 3.0)
    # divido per 8 pi mu R^3
    return 1.0 + r_cw ** 3 + 0.75 * r_cw * d_rel ** 2


def inerzia_ratio(R=None, d_rel=2.0, frazione_volume=0.5, rho=2710.0):
    """I_tot/I_rot per la stessa geometria.  Indipendente da R.

    I_rot = (2/5) m R^2.  La bead ha massa (m/2) e raggio r_cw R, quindi
        I_cw = (m/2)[d^2 + (2/5) r_cw^2 R^2].
    Dividendo per I_rot:
        I_tot/I_rot = 1 + frazione_volume (d_rel^2 + 0.4 r_cw^2) / 0.4.
    """
    r_cw = frazione_volume ** (1.0 / 3.0)
    return 1.0 + frazione_volume * (d_rel ** 2 + 0.4 * r_cw ** 2) / 0.4


BETA_RATIO = beta_ratio()
INERZIA_RATIO = inerzia_ratio()


def beta_no(R, mu):
    return 8.0 * math.pi * mu * R ** 3


def beta_cw(R, mu):
    return BETA_RATIO * beta_no(R, mu)


# ---------------------------------------------------------------------------
# 5. TERMICA
# ---------------------------------------------------------------------------
def potenza_assorbita_particella(P_inc, R, lam_nm, w0=W0):
    """Potenza assorbita dalla sfera di vaterite (Mie + tetto geometrico).

    Q_abs viene da Mie con n=1.65, k=1e-4 (valore dichiarato: la vaterite
    reale puo' avere k minore, e il termine scala linearmente con k).
    Il tetto geometrico usa alpha = 1 - exp(-4 pi k L / lambda), L = 2R.
    """
    d = OTTICA[lam_nm]
    q = m.q_ottico("vaterite", R, lam_nm * 1e-9)
    I0 = 2.0 * P_inc / (math.pi * w0 ** 2)
    P_mie = I0 * math.pi * R ** 2 * max(q["Q_abs"], 0.0)
    k = m.MATERIALI["vaterite"]["k"]
    alpha = 1.0 - math.exp(-4.0 * math.pi * k * (2.0 * R) /
                           (lam_nm * 1e-9))
    P_cap = m.eta_int(R, "sfera", w0) * P_inc * alpha
    P_abs = min(P_mie, P_cap)
    dT = P_abs / (4.0 * math.pi * SANGUE["kappa"] * R)
    return {"P_abs": P_abs, "P_abs_su_P": P_abs / P_inc if P_inc > 0 else 0.0,
            "Q_abs": q["Q_abs"], "Q_sca": q["Q_sca"], "dT": dT}


def dT_cammino(P_inc, lam_nm, z, w0=W0):
    """Sorgente lineare di assorbimento nel sangue.

    dT = [mu_a P_inc/(2 pi kappa)] * [(1-e^(-mu_t z))/(mu_t z)] *
         ln(1 + z/w0)
    Solo mu_a scalda; mu_s toglie luce al fascio ma non deposita calore nei
    primi micron (g ~ 0.98).  Il fattore (1-e^(-mu_t z))/(mu_t z) e' la media
    del deposito; per mu_t z << 1 tende a 1.  w0 e' il waist locale: se il
    fascio viene focalizzato con w0 proporzionale a lambda, il taglio
    cilindrico cambia e va passato esplicitamente.
    """
    if z <= 0.0:
        return 0.0
    d = OTTICA[lam_nm]
    mu_a = d["mu_a"]
    mu_t = mu_a + d["mu_s"]
    x = mu_t * z
    f = (1.0 - math.exp(-x)) / x if x > 1e-12 else 1.0
    return (mu_a * P_inc / (2.0 * math.pi * SANGUE["kappa"])) * f * \
        math.log1p(z / w0)


# ---------------------------------------------------------------------------
# 6. IL CICLO
# ---------------------------------------------------------------------------
def ciclo(R, lam_nm, z, dn=DN_VATERITE, mu=None, Pe=PE_R,
          sigma=SIGMA_DELTA, path_factor=2.0, sigma_y=SIGMA_Y,
          delta_cfl=None, w0=W0):
    """Ciclo di 720 gradi a profondita' z, tarato su Pe_r.

    z = profondita' della superficie prossimale sotto l'interfaccia.
    delta_cfl permette di valutare la stessa cella con uno spessore di CFL
    diverso da quello nominale (usato da continuita()).
    w0 e' il waist locale; il default W0 = 0.5 um riproduce la V7.2.
    Il latch e' a meta' ciclo: primo giro senza contrappeso, secondo con.
    """
    mu = mu_eff(R, z, delta_cfl) if mu is None else mu
    tr = coppia_per_watt(R, lam_nm, dn, sigma, path_factor, w0)
    if tr["tau_per_W"] <= 0.0:
        return None

    tau = Pe * KT
    P_target = tau / tr["tau_per_W"]
    T_s = trasmissione(lam_nm, z)
    P_inc = P_target / T_s

    b_no = beta_no(R, mu)
    b_cw = beta_cw(R, mu)
    b_half = 0.5 * (b_no + b_cw)
    omega_no = tau / b_no
    omega_cw = tau / b_cw
    T360_no = 2.0 * math.pi * b_no / tau
    T360_cw = 2.0 * math.pi * b_cw / tau
    T720_no = 2.0 * T360_no
    T720_cw = 2.0 * T360_cw
    T720_half = T360_no + T360_cw
    E720_no = P_inc * T720_no
    E720_cw = P_inc * T720_cw
    E720_half = P_inc * T720_half

    part = potenza_assorbita_particella(P_inc, R, lam_nm, w0)
    dT_part = part["dT"]
    dT_path = dT_cammino(P_inc, lam_nm, z, w0)
    dT = dT_part + dT_path

    tau_B = (8.0 / 3.0) * math.pi * R ** 3 * sigma_y
    ratio_B = tau / tau_B if tau_B > 0 else float("inf")
    f_cfl = frazione_cfl(R, z, delta_cfl)
    tau_B_eff = (8.0 / 3.0) * math.pi * R ** 3 * sigma_y * (1.0 - f_cfl)
    ratio_B_eff = tau / tau_B_eff if tau_B_eff > 0 else float("inf")

    # Inerzia e transitorio (per la scena e per il testo)
    rho_v = m.MATERIALI["vaterite"]["rho"]
    I_rot = 0.4 * rho_v * (4.0 / 3.0 * math.pi * R ** 3) * R ** 2
    I_tot = INERZIA_RATIO * I_rot
    tau_spin_no = I_rot / b_no
    tau_spin_cw = I_tot / b_cw

    # Rumore: varianza per giro = 4 pi / Pe in entrambi gli stati
    sigma_theta_360 = math.sqrt(4.0 * math.pi / Pe)
    sigma_theta_720 = math.sqrt(8.0 * math.pi / Pe)
    cv_durata = 1.0 / math.sqrt(math.pi * Pe)
    separazione = b_cw / b_no
    # soglia geometrica fra le due durate
    soglia = math.sqrt(T360_no * T360_cw)
    # errore gaussiano di classificazione per singolo giro: approssimazione
    # della V7.2, ritrattata in V7.3.  Serve solo a riprodurre i numeri
    # storici; il calcolo valido e' la prima passa esatta di v73_rumore.py.
    z_no = (soglia - T360_no) / (cv_durata * T360_no)
    z_cw = (T360_cw - soglia) / (cv_durata * T360_cw)
    err_no = 0.5 * math.erfc(z_no / math.sqrt(2.0))
    err_cw = 0.5 * math.erfc(z_cw / math.sqrt(2.0))

    # Efficienza meccanica: lavoro tau*4pi diviso energia ottica fornita
    W720 = tau * 4.0 * math.pi
    eta_mech = W720 / E720_half if E720_half > 0 else 0.0

    # Reynolds rotazionale (massimo, nel primo giro senza contrappeso)
    Re_rot = SANGUE["rho"] * omega_no * R ** 2 / mu

    motivi = []
    if dT > D_T_OPERATIVA:
        motivi.append("dT > 20 K")
    if ratio_B < 1.0:
        motivi.append("tau < tau_B (Bingham)")
    if omega_no > 0.5 * m.omega_rottura(R, "vaterite")[0]:
        motivi.append("omega > 0.5 omega_rottura")
    if T_s < math.exp(-1.0):
        motivi.append("z oltre la profondita' balistica")
    if P_inc > P_MAX:
        motivi.append("P_inc > 1 W")

    motivi_sicuri = list(motivi)
    if dT > D_T_SICURA and "dT > 5 K" not in motivi_sicuri:
        motivi_sicuri.append("dT > 5 K (soglia prudente)")

    motivi_cfl = [x for x in motivi if "Bingham" not in x]
    if ratio_B_eff < 1.0:
        motivi_cfl.append("tau < tau_B,eff (Bingham nel CFL)")

    return {
        "R": R, "R_um": R * 1e6, "lam_nm": lam_nm, "z": z, "z_um": z * 1e6,
        "mu": mu, "mu_mPa_s": mu * 1e3, "f_cfl": f_cfl,
        "T_s": T_s, "delta": tr["delta"], "c_circ": tr["c_circ"],
        "eta_int": tr["eta"], "tau_per_W": tr["tau_per_W"],
        "P_target": P_target, "P_target_mW": P_target * 1e3,
        "P_inc": P_inc, "P_inc_mW": P_inc * 1e3,
        "tau": tau, "beta_no": b_no, "beta_cw": b_cw, "beta_half": b_half,
        "omega_no": omega_no, "omega_cw": omega_cw,
        "T360_no": T360_no, "T360_cw": T360_cw,
        "T720_no": T720_no, "T720_cw": T720_cw, "T720_half": T720_half,
        "E720_no": E720_no, "E720_cw": E720_cw, "E720_half": E720_half,
        "costo_energia_720": E720_half / E720_no if E720_no > 0 else None,
        "costo_energia_360": E720_half / (0.5 * E720_no) if E720_no > 0 else None,
        "dT_part": dT_part, "dT_path": dT_path, "dT": dT,
        "P_abs_part": part["P_abs"], "P_abs_su_P": part["P_abs_su_P"],
        "Q_abs": part["Q_abs"],
        "tau_B": tau_B, "tau_su_tau_B": ratio_B,
        "tau_B_eff": tau_B_eff, "tau_su_tau_B_eff": ratio_B_eff,
        "I_rot": I_rot, "I_tot": I_tot,
        "tau_spin_no": tau_spin_no, "tau_spin_cw": tau_spin_cw,
        "sigma_theta_360": sigma_theta_360,
        "sigma_theta_360_deg": math.degrees(sigma_theta_360),
        "sigma_theta_720": sigma_theta_720,
        "sigma_theta_720_deg": math.degrees(sigma_theta_720),
        "cv_durata": cv_durata, "separazione_durate": separazione,
        "soglia_durate": soglia,
        "err_no": err_no, "err_cw": err_cw,
        "z_no": z_no, "z_cw": z_cw,
        "T_misura_3cicli": 3.0 * T720_half,
        "W720": W720, "W720_kT": W720 / KT, "eta_mech": eta_mech,
        "Re_rot": Re_rot,
        "valido": len(motivi) == 0, "motivi": motivi,
        "valido_sicuro": len(motivi_sicuri) == 0,
        "motivi_sicuri": motivi_sicuri,
        "valido_cfl": len(motivi_cfl) == 0, "motivi_cfl": motivi_cfl,
        "omega_rottura": m.omega_rottura(R, "vaterite")[0],
    }


# ---------------------------------------------------------------------------
# 7. GRIGLIA E SINTESI
# ---------------------------------------------------------------------------
def griglia():
    return [ciclo(R, lam, z) for lam in LAM for z in Z_GRID for R in R_GRID]


def _migliore(righe, lam, chiave="valido"):
    cand = [r for r in righe if r["lam_nm"] == lam and r[chiave]]
    if not cand:
        return None
    b = min(cand, key=lambda r: r["E720_half"])
    return {k: b[k] for k in (
        "R_um", "z_um", "P_target_mW", "P_inc_mW", "T720_half",
        "E720_half", "dT", "dT_part", "dT_path", "tau_su_tau_B",
        "tau_su_tau_B_eff", "mu_mPa_s", "eta_mech", "c_circ", "delta")}


def _finestre(righe, lam, z, chiave="valido"):
    rr = sorted([r for r in righe if r["lam_nm"] == lam and
                 abs(r["z"] - z) < 1e-15], key=lambda r: r["R"])
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
        res.append({"R_min_um": a["R_um"], "R_max_um": b["R_um"],
                    "n": len(sub),
                    "max_dT_K": max(x["dT"] for x in sub),
                    "min_T720_s": min(x["T720_half"] for x in sub)})
    return res


def riassunto(righe):
    out = {"celle": len(righe),
           "valide": sum(1 for r in righe if r["valido"]),
           "valide_sicure": sum(1 for r in righe if r["valido_sicuro"]),
           "valide_cfl": sum(1 for r in righe if r["valido_cfl"])}
    out["profondita"] = {
        str(lam): {"mu_a": OTTICA[lam]["mu_a"],
                   "mu_s": OTTICA[lam]["mu_s"],
                   "mu_t": OTTICA[lam]["mu_a"] + OTTICA[lam]["mu_s"],
                   "ell_um": profondita_balistica(lam) * 1e6,
                   "T_a_2um": trasmissione(lam, 2e-6),
                   "T_a_10um": trasmissione(lam, 10e-6)}
        for lam in LAM}
    for lam in LAM:
        out["migliore_%d" % lam] = _migliore(righe, lam, "valido")
        out["migliore_sicuro_%d" % lam] = _migliore(righe, lam, "valido_sicuro")
        out["migliore_cfl_%d" % lam] = _migliore(righe, lam, "valido_cfl")
    out["finestre"] = {
        str(lam): {("%.0f_um" % (z * 1e6)): _finestre(righe, lam, z,
                                                      "valido")
                   for z in Z_GRID} for lam in LAM}
    out["finestre_sicure"] = {
        str(lam): {("%.0f_um" % (z * 1e6)): _finestre(righe, lam, z,
                                                      "valido_sicuro")
                   for z in Z_GRID} for lam in LAM}
    out["finestre_cfl"] = {
        str(lam): {("%.0f_um" % (z * 1e6)): _finestre(righe, lam, z,
                                                      "valido_cfl")
                   for z in Z_GRID} for lam in LAM}
    return out


def confronto_lambda(z=2e-6):
    righe = []
    for R in (0.3e-6, 0.4e-6, 0.5e-6, 0.6e-6, 0.8e-6):
        for lam in LAM:
            c = ciclo(R, lam, z)
            righe.append(c)
    return righe


def costo_energetico(z=2e-6, lam=800):
    righe = []
    for R in (0.3e-6, 0.4e-6, 0.5e-6, 0.6e-6, 0.8e-6, 1.0e-6):
        c = ciclo(R, lam, z)
        righe.append(c)
    return righe


def viscosita(lam=800, z=2e-6):
    righe = []
    for mu_mpas in (1.3, 2.0, 4.0, 7.0, 10.0):
        for R in (0.3e-6, 0.5e-6, 0.6e-6):
            c = ciclo(R, lam, z, mu=mu_mpas * 1e-3)
            righe.append(c)
    return righe


def rumore(lam=800, z=2e-6):
    righe = []
    for R in (0.3e-6, 0.5e-6, 0.6e-6, 0.8e-6):
        c = ciclo(R, lam, z)
        righe.append(c)
    return righe


def continuita(lam=800, z=2e-6):
    """Dove il rotore sta davvero nel CFL, e dove il continuo regge.

    Per ogni spessore del CFL e per ogni raggio: frazione di volume nel CFL,
    viscosita' effettiva, tau/tau_B,eff, e se la sfera e' INTERAMENTE nel CFL
    (z + 2R <= delta).  La soglia di continuo e' R >= 1 um: sotto, il rotore
    e' confrontabile con la membrana e con il gap fra globuli rossi.
    """
    out = []
    for delta in (1.5e-6, 2.0e-6, 3.0e-6):
        for R in (0.3e-6, 0.5e-6, 0.6e-6, 0.8e-6, 1.0e-6, 1.2e-6):
            c = ciclo(R, lam, z, delta_cfl=delta)
            f = frazione_cfl(R, z, delta)
            dentro = (z + 2.0 * R) <= delta
            out.append({
                "delta_cfl_um": delta * 1e6, "R_um": R * 1e6,
                "f_cfl": f, "mu_mPa_s": c["mu_mPa_s"],
                "tau_su_tau_B": c["tau_su_tau_B"],
                "tau_su_tau_B_eff": c["tau_su_tau_B_eff"],
                "interamente_nel_CFL": bool(dentro),
                "continuo_valido": bool(R >= 1.0e-6),
                "valido_cfl": c["valido_cfl"],
            })
    return out


def fuga_bingham(lam=800, z=2e-6,
                 R_list=(1.0e-6, 1.5e-6, 2.0e-6, 3.0e-6, 5.0e-6),
                 margine=1.2):
    """Che cosa costa superare lo snervamento alzando Pe invece di R.

    A Pe=10 la coppia e' 10 kT; se tau_B > 10 kT il rotore non parte.
    Ma la coppia non e' fissata: si puo' salire a Pe > tau_B/kT.  L'energia
    per ciclo NON cambia (E = W/eta, W = 4 pi Pe kT e P proporzionale a Pe),
    cambiano potenza, tempo e calore.  Qui si calcola il punto minimo
    Pe = margine * tau_B/kT e si vede dove il calore lo esclude.
    """
    out = []
    for R in R_list:
        c10 = ciclo(R, lam, z)
        Pe_B = c10["tau_B"] / KT
        Pe = margine * Pe_B
        c = ciclo(R, lam, z, Pe=Pe)
        out.append({
            "R_um": R * 1e6, "Pe_B": Pe_B, "Pe": Pe,
            "tau_B": c10["tau_B"],
            "P_target_mW": c["P_target_mW"], "P_inc_mW": c["P_inc_mW"],
            "T720_half": c["T720_half"], "E720_half": c["E720_half"],
            "dT": c["dT"], "dT_path": c["dT_path"],
            "tau_su_tau_B": c["tau_su_tau_B"],
            "valido": c["valido"], "motivi": c["motivi"],
        })
    return out


def spot_secondario(lam, R=0.5e-6, z=2e-6, P_add=1.0e-3):
    """Costo del secondo spot: potenza e calore aggiunti.

    Il secondo fascio non cambia P_target (la coppia viene dal primo), ma
    si somma alla potenza incidente e al riscaldamento di cammino.  Il
    rapporto P_add/P_inc e' quello che il testo deve citare, non un
    'circa doppio' scritto a mano.
    """
    c = ciclo(R, lam, z)
    dT_add = dT_cammino(P_add, lam, z)
    return {
        "lam_nm": lam, "R_um": R * 1e6, "z_um": z * 1e6,
        "P_inc_mW": c["P_inc_mW"], "P_add_mW": P_add * 1e3,
        "P_tot_mW": c["P_inc_mW"] + P_add * 1e3,
        "rapporto_tot_primario": 1.0 + P_add / c["P_inc"],
        "dT_path_primario": c["dT_path"],
        "dT_path_aggiunto": dT_add,
        "dT_path_totale": c["dT_path"] + dT_add,
        "dT_totale": c["dT"] + dT_add,
    }


def zeri():
    righe = []
    for dn in (0.06, 0.08, 0.10):
        for lam in LAM:
            Rm = lam * 1e-9 / (2.0 * dn)
            tau_B = (8.0 / 3.0) * math.pi * Rm ** 3 * SIGMA_Y
            # Coppia ottica massima a 1 W sul bersaglio: eta=1, c_circ=1.
            tau_opt_1W = 1.0 / (2.0 * math.pi * C_LUCE / (lam * 1e-9))
            righe.append({"dn": dn, "lam_nm": lam, "R_m_um": Rm * 1e6,
                          "tau_B": tau_B,
                          "tau_opt_1W": tau_opt_1W,
                          "tau_opt_su_tau_B": tau_opt_1W / tau_B,
                          "tau_su_tau_B": PE_R * KT / tau_B,
                          "delta": 2 * math.pi})
    return righe


# ---------------------------------------------------------------------------
# 8. SCRITTURA JSON
# ---------------------------------------------------------------------------
def _scrivi(nome, payload):
    p = os.path.join(DATI, nome)
    with open(p, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1, ensure_ascii=False, default=float)
    print("  scritto %-28s %6.1f kB" % (nome, os.path.getsize(p) / 1024.0))


def main():
    print("=" * 78)
    print("V7.3 - ciclo di 720 gradi: contatore Z2, sangue, CFL, Bingham")
    print("=" * 78)
    print("\n-- parametri geometrici")
    print("   beta_tot/beta_rot = %.4f   I_tot/I_rot = %.4f"
          % (BETA_RATIO, INERZIA_RATIO))
    vs = verifica_scaling()
    print("   scaling: tau/P(1064)/tau/P(800) = %.4f (atteso %.4f, %s); "
          "P(0.10)/P(0.06) = %.4f (atteso %.4f, %s)"
          % (vs["lambda_1064_su_800"], vs["lambda_atteso"],
             "ok" if vs["lambda_ok"] else "FALLITO",
             vs["dn_010_su_006"], vs["dn_atteso"],
             "ok" if vs["dn_ok"] else "FALLITO"))
    print("   Delta_n = %.2f (range %.2f-%.2f); Pe_r = %.1f"
          % (DN_VATERITE, DN_RANGE[0], DN_RANGE[1], PE_R))
    print("\n-- ottica del sangue")
    for lam in LAM:
        d = OTTICA[lam]
        print("   %4d nm: mu_a=%7.0f 1/m  mu_s=%6.0f 1/m  "
              "1/mu_t=%5.1f um  T(2um)=%.3f  T(10um)=%.3f"
              % (lam, d["mu_a"], d["mu_s"], profondita_balistica(lam) * 1e6,
                 trasmissione(lam, 2e-6), trasmissione(lam, 10e-6)))

    righe = griglia()
    rs = riassunto(righe)
    cf = confronto_lambda()
    ce = costo_energetico()
    vis = viscosita()
    rum = rumore()
    cont = continuita()
    fuga = fuga_bingham()
    spot = [spot_secondario(lam) for lam in LAM]
    zr = zeri()
    rif = ciclo(0.5e-6, 800, 2e-6)

    _scrivi("V73_griglia.json", {
        "descrizione": "V7.3: griglia R-z, 532/800/1064 nm, vaterite, "
                       "Pe_r=10, latch a meta' ciclo, CFL e Bingham",
        "sangue": SANGUE, "ottica": OTTICA,
        "dn": DN_VATERITE, "dn_range": DN_RANGE,
        "Pe_r": PE_R, "sigma_y": SIGMA_Y,
        "soglie": {"dT_operativa": D_T_OPERATIVA, "dT_sicura": D_T_SICURA},
        "R_um": [r * 1e6 for r in R_GRID],
        "z_um": [z * 1e6 for z in Z_GRID],
        "lambda_nm": list(LAM),
        "celle": righe})

    _scrivi("V73_sintesi.json", {
        "descrizione": "V7.3: correzioni V7.2, ciclo 720, costo, rumore, "
                       "finestre, zeridi forma, confronto lambda",
        "correzioni": CORREZIONI,
        "sangue": SANGUE, "ottica": OTTICA,
        "dn": DN_VATERITE, "dn_range": DN_RANGE,
        "Pe_r": PE_R, "sigma_y": SIGMA_Y,
        "soglie": {"dT_operativa": D_T_OPERATIVA, "dT_sicura": D_T_SICURA},
        "beta_ratio": BETA_RATIO, "inerzia_ratio": INERZIA_RATIO,
        "verifica_scaling": verifica_scaling(),
        "riassunto": rs,
        "confronto_lambda": cf,
        "costo_energetico": ce,
        "viscosita": vis,
        "rumore": rum,
        "rumore_gaussiano_ritrattato": {
            "campi": ["confronto_lambda[].err_no", "confronto_lambda[].err_cw",
                      "confronto_lambda[].z_no", "confronto_lambda[].z_cw",
                      "confronto_lambda[].T_misura_3cicli",
                      "confronto_lambda[].soglia_durate"],
            "nota": "valori della V7.2 tenuti solo per riprodurre la "
                    "ritrattazione; il rumore valido e' la prima passa "
                    "esatta in V73_rumore.json",
        },
        "continuita": cont,
        "fuga_bingham": fuga,
        "spot_secondario": spot,
        "zeri": zr,
        "cella_riferimento_3d": {
            "R_um": rif["R_um"], "z_um": rif["z_um"],
            "lam_nm": rif["lam_nm"],
            "P_target_uW": rif["P_target"] * 1e6,
            "P_inc_mW": rif["P_inc_mW"],
            "T_360_no_s": rif["T360_no"],
            "T_360_cw_s": rif["T360_cw"],
            "T_720_s": rif["T720_half"],
            "E_720_mJ": rif["E720_half"] * 1e3,
            "dT_K": rif["dT"], "dT_part_K": rif["dT_part"],
            "dT_path_K": rif["dT_path"],
            "tau_su_tau_B": rif["tau_su_tau_B"],
            "beta_cw_su_beta_no": rif["beta_cw"] / rif["beta_no"],
            "inerzia_su_rot": INERZIA_RATIO,
            "modulazione": 1.0 - rif["omega_cw"] / rif["omega_no"],
            "mu_mPa_s": rif["mu_mPa_s"],
        },
        "relazione": "tau = Pe_r kT; W720 = 4 pi tau; T720_half = "
                     "2 pi (beta_no + beta_cw)/tau; E720 = P_inc T720; "
                     "eta_mech = W720/E720.",
    })

    print("\n-- celle valide: %d / %d (sicure %d, CFL %d)"
          % (rs["valide"], rs["celle"], rs["valide_sicure"],
             rs["valide_cfl"]))
    for lam in LAM:
        b = rs["migliore_%d" % lam]
        if b:
            print("   %4d nm: R=%.2f um z=%.0f um  Ptar=%.3f mW  "
                  "Pinc=%.3f mW  T720=%.2f s  E=%.2f mJ  dT=%.2f K "
                  "(part %.2f, path %.2f)  tau/tauB=%.2f  mu=%.2f mPa s"
                  % (lam, b["R_um"], b["z_um"], b["P_target_mW"],
                     b["P_inc_mW"], b["T720_half"], b["E720_half"] * 1e3,
                     b["dT"], b["dT_part"], b["dT_path"],
                     b["tau_su_tau_B"], b["mu_mPa_s"]))
        else:
            print("   %4d nm: nessuna cella valida" % lam)

    print("\n-- riferimento R=0.5 um, z=2 um, 800 nm")
    c = ciclo(0.5e-6, 800, 2e-6)
    for k in ("P_target_mW", "P_inc_mW", "T360_no", "T720_half",
              "E720_half", "dT_part", "dT_path", "dT",
              "tau_su_tau_B", "tau_su_tau_B_eff", "eta_mech",
              "sigma_theta_360_deg", "cv_durata", "W720_kT",
              "mu_mPa_s"):
        print("   %-24s %.6g" % (k, c[k]))

    print("\n-- confronto lambda a R=0.5 um, z=2 um")
    for r in cf:
        if abs(r["R_um"] - 0.5) < 1e-9:
            print("   %4d nm: Ptar=%.3f mW Pinc=%.3f mW T720=%.2f s "
                  "E=%.2f mJ dT=%.2f K (part %.2f path %.2f) tau/tauB=%.2f"
                  % (r["lam_nm"], r["P_target_mW"], r["P_inc_mW"],
                     r["T720_half"], r["E720_half"], r["dT"],
                     r["dT_part"], r["dT_path"], r["tau_su_tau_B"]))

    print("\n-- costo del contrappeso (800 nm, z=2 um)")
    for r in ce:
        print("   R=%.1f um: T360no=%.2f T720no=%.2f T720cw=%.2f "
              "T720half=%.2f  Ehalf=%.2f mJ  Ehalf/Eno=%.2f"
              % (r["R_um"], r["T360_no"], r["T720_no"], r["T720_cw"],
                 r["T720_half"], r["E720_half"] * 1e3,
                 r["costo_energia_720"]))

    print("\n-- rumore e test di durata")
    for r in rum:
        print("   R=%.1f um: sigma360=%.2f deg  sigma720=%.2f deg  "
              "CV=%.3f  beta_cw/beta_no=%.2f  err_no=%.2e err_cw=%.2e"
              % (r["R_um"], r["sigma_theta_360_deg"],
                 r["sigma_theta_720_deg"], r["cv_durata"],
                 r["separazione_durate"], r["err_no"], r["err_cw"]))

    print("\n-- zeri di forma e Bingham")
    for r in zr:
        print("   dn=%.2f lam=%4d nm: R_m=%5.2f um  tau_opt(1W)/tauB=%.3f  "
              "Pe kT/tauB=%.2f"
              % (r["dn"], r["lam_nm"], r["R_m_um"],
                 r["tau_opt_su_tau_B"], r["tau_su_tau_B"]))

    print("\n-- CFL, continuo e Bingham (800 nm, z=2 um)")
    for r in cont:
        if abs(r["delta_cfl_um"] - 3.0) < 1e-9:
            print("   delta=%.1f R=%.1f um: f=%.2f mu=%.2f  tau/tauB=%.2f "
                  "tau/tauB,eff=%.2f  tutto CFL=%s  cont=%s  val=%s"
                  % (r["delta_cfl_um"], r["R_um"], r["f_cfl"],
                     r["mu_mPa_s"], r["tau_su_tau_B"],
                     r["tau_su_tau_B_eff"], r["interamente_nel_CFL"],
                     r["continuo_valido"], r["valido_cfl"]))

    print("\n-- secondo spot da 1 mW (z=2 um, R=0.5 um)")
    for s in spot:
        print("   %4d nm: Pinc=%.3f mW  Ptot=%.3f mW  tot/prim=%.2f  "
              "dT_path +%.2f K  dT tot=%.2f K"
              % (s["lam_nm"], s["P_inc_mW"], s["P_tot_mW"],
                 s["rapporto_tot_primario"], s["dT_path_aggiunto"],
                 s["dT_totale"]))

    print("\n-- fuga dallo snervamento alzando Pe (800 nm, z=2 um)")
    for f in fuga:
        print("   R=%.1f um: Pe_B=%.1f  Pe=%.1f  Pinc=%.2f mW  T720=%.3f s  "
              "E=%.2f mJ  dT=%.1f K  val=%s %s"
              % (f["R_um"], f["Pe_B"], f["Pe"], f["P_inc_mW"],
                 f["T720_half"], f["E720_half"] * 1e3, f["dT"],
                 f["valido"], f["motivi"]))

    print("\nFatto.")


if __name__ == "__main__":
    main()
