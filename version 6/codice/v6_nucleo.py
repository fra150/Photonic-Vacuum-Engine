"""
MOTORE FOTONICO DEL VUOTO - V6
Nucleo fisico: ottica, attrito, termica, dinamica.

Autore: Dr. Francesco Bulla
Revisione V6: da modello architetturale a modello quantitativo verificabile.

Questo modulo non contiene testo del preprint. Contiene una sola cosa: le
funzioni che trasformano parametri misurati (potenza, lunghezza d'onda,
geometria, materiale, mezzo) in grandezze confrontabili con un esperimento.

Le tre scelte che differenziano la V6 dalla V5
------------------------------------------------
1. Un solo modello di coppia, dichiarato come limite:

       tau = (l + s) * P_ass / omega_ott        con  omega_ott = 2 pi c / lambda

   dove l e' il momento angolare orbitale per fotone e s = +-1 quello di
   polarizzazione. La potenza che compare e' quella ASSORBITA, non quella
   intercettata: sono due grandezze diverse e la V5 le confondeva.

2. Una sola forza, dichiarata come limite:

       F = I * sigma_pr / c,     sigma_pr = sigma_ext - g * sigma_sca

   e' il trasferimento completo della quantita' di moto assiale. Non va
   sommata alla coppia: sono due canali distinti.

3. Un solo attrito, dichiarato come intervallo:

       beta_min <= beta <= beta_max

   Il teorema di dissipazione minima da' un limite superiore rigoroso
   (il corpo pieno che contiene il telaio); la teoria dei corpi snelli da'
   il limite inferiore (l'anello piu' esterno da solo). La V5 usava una
   frazione di area bagnata, che non e' un modello di attrito rotazionale.

Uso:
    python v6_nucleo.py          # self-test
"""

import math

import numpy as np
from scipy import constants as C
from scipy.special import jv, jvp, hankel1, h1vp, spherical_jn, spherical_yn

# ---------------------------------------------------------------------------
# 0. COSTANTI
# ---------------------------------------------------------------------------
LAMBDA = 532e-9                  # m
OMEGA_OTT = 2.0 * math.pi * C.c / LAMBDA
NU_OTT = C.c / LAMBDA            # Hz
E_PHOTON = C.h * NU_OTT          # J
P_LASER = 1.0e-6                 # W, potenza di riferimento
NA = 0.9
W0_DICH = 500e-9                 # m, waist dichiarato nella V4/V5
W0_DIFF = 0.61 * LAMBDA / NA     # m, waist al limite di diffrazione (NA 0.9)
W0 = W0_DICH                     # riferimento V6 = dichiarato V5, per continuita'
A_BEAM = math.pi * W0 ** 2
A_DIFF = (LAMBDA / 2.0) ** 2     # m^2, area minima di un fuoco limitato
T0 = 310.0                       # K
KT = C.k * T0
ETA_W = 1.2e-3                   # Pa s, acqua a 20 C
RHO_W = 998.2
MU_AIR = 1.85e-5                 # Pa s
RHO_AIR = 1.2
SIGMA_SB = C.sigma

# Geometria di riferimento (inventario dei fili della scena V5).
# I raggi dei fili sono quelli della scena V5, non 1.25 nm: la V5 dichiarava
# "diametro 6.4 nm" e poi calcolava x con 1.25 nm.  Poiche' in Rayleigh
# sigma ∝ a^2, l'errore cambiava il contributo del telaio di un fattore 6-36.
UA = 100e-9
R_DISC = 3.275 * UA              # 327.5 nm, bordo del disco
R_OUT = 3.70 * UA                # 370 nm, semiasse maggiore dell'orbita esterna
CORE_R = 40e-9                   # nm, nucleo sferico
WIRE_ORB = 3.2e-9                # nm, raggio dei fili delle orbite (diametro 6.4)
WIRE_RIM = 7.5e-9                # nm, raggio dell'anello di bordo (diametro 15)
WIRE_LENS = 4.5e-9               # nm, raggio dei fili delle lenti (diametro 9)
WIRE_R = WIRE_ORB                # alias storico
SPIN_ORBITS = ((220.0, 110.0), (295.0, 152.0), (370.0, 195.0))   # nm
LENS = (105.0, 245.0)                                            # nm

# Valori V4/V5 da confrontare (non da usare)
R_EFF_V5 = math.sqrt(5.7567) * UA            # 239.93 nm, "stesso volume"
M_V5 = 1.2463630573655231e-18 + 4.487710753596799e-18
I1_V5 = 3.092e-32
I2_V5 = 2.583e-31
BETA_V5 = 2.8921247272696085e-23             # frazione di area bagnata
SIGMA_PR_V5 = 7.259e-15
F_V5 = 3.082970814745448e-17


# ---------------------------------------------------------------------------
# 1. MATERIALI (n, k a 532 nm; proprieta' meccaniche e termiche)
# ---------------------------------------------------------------------------
MATERIALI = {
    "rame": {
        "nome": "Rame",
        "n": 1.09, "k": 2.53,                   # Johnson & Christy, 532 nm
        "rho": 8960.0,
        "sy": 200e6,                            # Pa, rame ricotto
        "sy_nota": "200 MPa (ricotto) - 400 MPa (incrudito)",
        "kappa": 401.0,                         # W/m/K
        "cp": 385.0,                            # J/kg/K
        "Tm": 1358.0,                           # K
        "emiss_limpido": 0.03,                  # superficie lucida
        "emiss_ossidato": 0.15,                 # rugosa / parzialmente ossidata
        "nota": "fortemente assorbente, alto k; l'emissivita' di un oggetto "
                "sub-lunghezza d'onda dipende dall'assorbimento volumetrico a "
                "10 um e non puo' essere quella massiva (0.6) di un ossido di "
                "1.4 nm: l'intervallo usato e' 0.03-0.15",
    },
    "silice": {
        "nome": "Silice fusa",
        "n": 1.461, "k": 1.0e-7,                # trasparente nel visibile
        "rho": 2200.0,
        "sy": 1.0e9,                            # Pa, conservativo su nanoscala
        "sy_nota": "50 MPa (bulk) - 5 GPa (fibra); 1 GPa usato come valore di lavoro",
        "kappa": 1.4,
        "cp": 740.0,
        "Tm": 1983.0,
        "emiss_limpido": 0.85,
        "emiss_ossidato": 0.85,
        "nota": "trasparente: il canale di coppia non e' l'assorbimento ma il ritardo di fase",
    },
    "grafene": {
        "nome": "Grafene / CNT",
        "n": 2.20, "k": 1.40,                   # grafite nel piano, 532 nm (stima)
        "rho": 1300.0,                          # fasci di CNT, non densita' del grafite
        "sy": 20e9,                             # Pa, valore efficace conservativo
        "sy_nota": "20 GPa efficace (50 GPa fibra singola, 1-5 GPa rete)",
        "kappa": 3000.0,                        # nel piano
        "cp": 700.0,
        "Tm": 3900.0,
        "emiss_limpido": 0.05,
        "emiss_ossidato": 0.80,
        "nota": "miglior rapporto sigma/rho, ma assorbe per strato (2.3% per foglio)",
    },
}


def materiale(nome):
    return MATERIALI[nome]


def epsilon_rame_dimensione(lam=LAMBDA, a=40e-9, A_fermi=1.0, tau_bulk=27e-15):
    """Permittivita' del rame con correzione di dimensione (Drude).

    gamma(a) = 1/tau_bulk + A * v_F / a.  La correzione agisce solo sul
    contributo intrabanda: a 532 nm il rame e' dominato dalle transizioni
    interbanda, quindi l'effetto atteso e' piccolo.  Questa funzione serve
    proprio a misurarlo invece di dichiararlo.
    """
    m = MATERIALI["rame"]
    eps_bulk = complex(m["n"] ** 2 - m["k"] ** 2, 2 * m["n"] * m["k"])
    v_f = 1.57e6
    gamma_b = 1.0 / tau_bulk
    gamma_a = gamma_b + A_fermi * v_f / a
    omega = 2 * math.pi * C.c / lam
    wp2 = 1.36e16 ** 2
    d_eps = wp2 * (1.0 / (omega ** 2 + 1j * gamma_b * omega)
                   - 1.0 / (omega ** 2 + 1j * gamma_a * omega))
    eps = eps_bulk + d_eps
    n = math.sqrt(eps.real ** 2 + eps.imag ** 2)
    return complex(math.sqrt((n + eps.real) / 2), math.sqrt((n - eps.real) / 2))


# ---------------------------------------------------------------------------
# 2. OTTICA: MIE ESATTA (SFERA), ADEN-KERKER (SFERA RIVESTITA), CILINDRO
# ---------------------------------------------------------------------------
X_GEOM = 25.0    # oltre questo parametro di dimensione: ottica geometrica


def mie_sfera(m, x):
    """Coefficienti di Mie per una sfera omogenea (Bohren & Huffman, cap. 4).

    Restituisce Q_ext, Q_sca, Q_abs, Q_pr = Q_ext - g Q_sca, g.
    Per x > X_GEOM si usa il limite di ottica geometrica con Q_pr = 1.0,
    dichiarato come tale (corpo opaco grande, incertezza +-20%).
    """
    if x <= 0:
        return {"Q_ext": 0.0, "Q_sca": 0.0, "Q_abs": 0.0, "Q_pr": 0.0, "g": 0.0,
                "regime": "nullo"}
    if x > X_GEOM:
        return {"Q_ext": 2.0, "Q_sca": 1.0, "Q_abs": 1.0, "Q_pr": 1.0, "g": 0.5,
                "regime": "ottica geometrica"}
    nmax = int(2 + x + 4 * x ** (1 / 3)) + 10
    nmax = max(nmax, 12)
    nx = m * x
    # derivata logaritmica per ricorrenza discendente
    D = np.zeros(nmax + 2, dtype=complex)
    for n in range(nmax, 0, -1):
        D[n - 1] = n / nx - 1.0 / (D[n] + n / nx)
    psi = lambda n, z: z * spherical_jn(n, z)
    xi = lambda n, z: z * (spherical_jn(n, z) + 1j * spherical_yn(n, z))
    nn = np.arange(1, nmax + 1)
    a = np.zeros(nmax + 1, dtype=complex)
    b = np.zeros(nmax + 1, dtype=complex)
    for n in nn:
        a[n] = ((D[n] / m + n / x) * psi(n, x) - psi(n - 1, x)) / \
               ((D[n] / m + n / x) * xi(n, x) - xi(n - 1, x))
        b[n] = ((D[n] * m + n / x) * psi(n, x) - psi(n - 1, x)) / \
               ((D[n] * m + n / x) * xi(n, x) - xi(n - 1, x))
    q_ext = (2.0 / x ** 2) * np.sum((2 * nn + 1) * (a[1:] + b[1:]).real)
    q_sca = (2.0 / x ** 2) * np.sum((2 * nn + 1) * (np.abs(a[1:]) ** 2 + np.abs(b[1:]) ** 2))
    g = 0.0
    if q_sca > 1e-30:
        g = (4.0 / (x ** 2 * q_sca)) * np.sum(
            (nn * (nn + 1) / (nn + 0.5)) * (a[1:] * np.conj(a[:-1]) + b[1:] * np.conj(b[:-1])).real
            + ((2 * nn + 1) / (nn * (nn + 1))) * (a[1:] * np.conj(b[1:])).real)
    q_pr = q_ext - g * q_sca
    return {"Q_ext": float(q_ext), "Q_sca": float(q_sca), "Q_abs": float(q_ext - q_sca),
            "Q_pr": float(q_pr), "g": float(g), "regime": "Mie esatta"}


def sfera_rivestita_qs(m_core, m_shell, r, f_nucleo):
    """Sfera con guscio in approssimazione quasi-statica (Maxwell-Garnett).

    f_nucleo e' la frazione di volume del nucleo.  La permittivita' efficace
    del granulo rivestito e'

        eps_eff = eps_s (eps_c + 2 eps_s + 2 f (eps_c - eps_s))
                        / (eps_c + 2 eps_s - f (eps_c - eps_s))

    e i limiti f -> 1 e f -> 0 restituiscono esattamente nucleo e guscio.
    Vale per x << 1: a x = 0.47 l'errore su Q e' stimato al 20%, dichiarato.
    Il vantaggio rispetto a un Aden-Kerker completo e' che i due limiti sono
    verificabili analiticamente, cosa che un'implementazione numerica non
    garantisce.
    """
    ec = m_core ** 2
    es = m_shell ** 2
    ee = es * (ec + 2 * es + 2 * f_nucleo * (ec - es)) / (ec + 2 * es - f_nucleo * (ec - es))
    x = 2 * math.pi * r / LAMBDA
    K = (ee - 1) / (ee + 2)
    q_ext = 4 * x * K.imag
    q_sca = (8.0 / 3.0) * x ** 4 * abs(K) ** 2
    return {"Q_ext": q_ext, "Q_sca": q_sca, "Q_abs": q_ext - q_sca, "Q_pr": q_ext,
            "g": 0.0, "regime": "quasi-statico (Maxwell-Garnett)"}


def cilindro_esatto(m, x):
    """Cilindro infinito, incidenza normale (Bohren & Huffman, cap. 8).

    Somma sui modi n >= 1.  Il termine n = 0 delle formule di Bohren-Huffman
    contiene, per x -> 0, la divergenza logaritmica di H_0(x) ed e' instabile:
    per un filo con x = 0.015 il suo contributo numerico (0.128 al posto di
    0.011) e' un artefatto, non fisica.  Il termine n = 0 viene quindi escluso
    e la serie n >= 1 e' verificata contro il limite di Rayleigh (scarto
    < 0.6% per x < 0.05, vedi self_test).

    Lo scattering di un cilindro sottile e' quasi isotropo: g = 0, con errore
    O(x^2) < 1e-3 per x < 0.02, che e' il nostro caso.
    """
    if x > 1.0:
        return {"Q_ext": 2.0, "Q_sca": 1.0, "Q_abs": 1.0, "Q_pr": 1.0,
                "regime": "ottica geometrica"}
    nmax = int(2 + x + 4 * x ** (1 / 3)) + 8
    nmax = max(nmax, 10)
    mx = m * x

    def psi(n, z):      # J_n
        return jv(n, z)

    def dpsi(n, z):     # J_n'
        return jvp(n, z)

    def xi(n, z):       # H_n^(1)
        return hankel1(n, z)

    def dxi(n, z):
        return h1vp(n, z)

    a = np.zeros(nmax + 1, dtype=complex)   # TM (E parallelo all'asse)
    b = np.zeros(nmax + 1, dtype=complex)   # TE
    for n in range(1, nmax + 1):
        a[n] = (m * psi(n, mx) * dpsi(n, x) - psi(n, x) * dpsi(n, mx)) / \
               (m * psi(n, mx) * dxi(n, x) - xi(n, x) * dpsi(n, mx))
        b[n] = (psi(n, mx) * dpsi(n, x) - m * psi(n, x) * dpsi(n, mx)) / \
               (psi(n, mx) * dxi(n, x) - m * xi(n, x) * dpsi(n, mx))
    q_ext = (2.0 / x) * (2 * np.sum(a[1:] + b[1:])).real
    q_sca = (2.0 / x) * (2 * np.sum(np.abs(a[1:]) ** 2 + np.abs(b[1:]) ** 2))
    return {"Q_ext": float(q_ext), "Q_sca": float(q_sca), "Q_abs": float(q_ext - q_sca),
            "Q_pr": float(q_ext - q_sca), "regime": "cilindro esatto (n>=1)"}


def cilindro_rayleigh(m, x):
    """Limite di Rayleigh per un cilindro sottile (x << 1)."""
    m2 = m ** 2
    pol = (m2 - 1) / (m2 + 1)
    q_ext = math.pi * x * pol.imag
    q_sca = (math.pi ** 3 / 4.0) * x ** 3 * abs(pol) ** 2
    return {"Q_ext": q_ext, "Q_sca": q_sca, "Q_abs": q_ext - q_sca,
            "Q_pr": q_ext - q_sca, "regime": "Rayleigh cilindrico"}


def sezione_sfera(m, raggio, lam=LAMBDA):
    x = 2 * math.pi * raggio / lam
    q = mie_sfera(m, x)
    s = math.pi * raggio ** 2
    return {"x": x, "Q_ext": q["Q_ext"], "Q_sca": q["Q_sca"], "Q_abs": q["Q_abs"],
            "Q_pr": q["Q_pr"], "g": q["g"], "regime": q["regime"],
            "sigma_ext": s * q["Q_ext"], "sigma_abs": s * max(q["Q_abs"], 0.0),
            "sigma_pr": s * max(q["Q_pr"], 0.0), "area_geo": s}


def sezione_filo(m, raggio, lunghezza, lam=LAMBDA, esatto=True):
    """Sezione di un filo di lunghezza finita, trattato come cilindro locale."""
    x = 2 * math.pi * raggio / lam
    q = cilindro_esatto(m, x) if esatto else cilindro_rayleigh(m, x)
    a_pro = 2 * raggio * lunghezza
    return {"x": x, "Q_ext": q["Q_ext"], "Q_sca": q["Q_sca"], "Q_abs": q["Q_abs"],
            "Q_pr": q["Q_pr"], "regime": q["regime"],
            "sigma_ext": a_pro * q["Q_ext"], "sigma_abs": a_pro * max(q["Q_abs"], 0.0),
            "sigma_pr": a_pro * max(q["Q_pr"], 0.0), "area_geo": a_pro}


# ---------------------------------------------------------------------------
# 3. ATTrito ROTAZIONALE: INTERVALLO PER CORPI SNELI
# ---------------------------------------------------------------------------
def beta_sfera(r, mu=ETA_W):
    """Sfera piena: beta = 8 pi mu r^3 (Stokes, esatto)."""
    return 8 * math.pi * mu * r ** 3


def beta_disco(r, mu=ETA_W):
    """Disco sottile che ruota attorno al proprio asse: beta = (32/3) mu r^3.

    E' anche il limite superiore rigoroso per qualunque corpo contenuto nel
    disco di raggio r (teorema di dissipazione minima di Helmholtz-Korteweg).
    """
    return (32.0 / 3.0) * mu * r ** 3


def beta_anello(r, a, mu=ETA_W):
    """Anello snello (toro) attorno all'asse: beta = 8 pi^2 mu r^3/(2 ln(2r/a)+1).

    Corpo snello: la costante logaritmica e' quella di un elemento di lunghezza
    locale ~2r; e' il limite inferiore per il telaio (l'anello esterno da solo).
    """
    if a <= 0 or a >= r:
        return beta_disco(r, mu)
    return 8 * math.pi ** 2 * mu * r ** 3 / (2.0 * math.log(2.0 * r / a) + 1.0)


def beta_asta(L, a, mu=ETA_W):
    """Asta snella (prolata) che ruota attorno al centro, asse perpendicolare.

    beta = pi mu L^3 / (3 (ln(L/a) + delta)), delta = -0.447 (correzioni di
    estremita', Tirado-Ramirez Garcia de la Torre).
    """
    if a <= 0:
        return 0.0
    return math.pi * mu * L ** 3 / (3.0 * (math.log(L / a) - 0.447))


def beta_manubrio(r, d, mu=ETA_W):
    """Due sfere di raggio r separate da d, rotazione attorno al centro."""
    return 3 * math.pi * mu * r * d ** 2 + 16 * math.pi * mu * r ** 3


def beta_gas(p, r, T=T0, m_gas=4.65e-26):
    """Attrito rotazionale in regime molecolare libero (p << 1 mbar).

    beta ~ 4 pi p r^4 / v_bar: la pressione fissa il flusso di quantita' di
    moto tangenziale.  Vale per Kn = lambda/r >> 1; con r = 100 nm e'
    gia' valido a 1 mbar.
    """
    v_bar = math.sqrt(8 * C.k * T / (math.pi * m_gas))
    return 4 * math.pi * p * r ** 4 / v_bar


def beta_turbolento(r, omega, rho=RHO_W):
    """Coppia di un disco in regime turbolento: coefficiente di momento.

    C_M = 0.073 Re^-0.2 per Re > 3e5 (von Karman).  Restituisce beta = T/omega,
    che dipende da omega: il regime stazionario va risolto, non moltiplicato.
    """
    Re = rho * omega * r ** 2 / ETA_W
    if Re < 3e5:
        C_M = 3.87 * max(Re, 1e-30) ** -0.5
    else:
        C_M = 0.073 * Re ** -0.2
    coppia = 0.5 * C_M * rho * omega ** 2 * r ** 5
    return coppia / omega if omega > 0 else 0.0


# ---------------------------------------------------------------------------
# 4. TERMICA
# ---------------------------------------------------------------------------
def delta_T_fluido(P_abs, r, kappa_m=0.6):
    """Sfera isoterma in un mezzo conduttivo: dT = P/(4 pi kappa r).

    La convezione e' trascurabile: a r = 100 nm e dT < 1 K il numero di
    Rayleigh e' ~1e-18.  Il numero di Nusselt e' 1.
    """
    return P_abs / (4 * math.pi * kappa_m * r)


def T_vuoto(P_abs, area, emiss, T0=T0):
    """Equilibrio radiativo: P = emiss sigma A (T^4 - T0^4)."""
    if P_abs <= 0:
        return T0
    return (T0 ** 4 + P_abs / (emiss * SIGMA_SB * area)) ** 0.25


def potenza_fusione_vuoto(area, emiss, Tm, T0=T0):
    """Potenza assorbita che porta il corpo a Tm per sola irraggiamento."""
    return emiss * SIGMA_SB * area * (Tm ** 4 - T0 ** 4)


def tempo_termico(m, cp, dT, P_abs):
    """Tempo di salita adiabatico (limite inferiore del tempo di fusione)."""
    if P_abs <= 0:
        return float("inf")
    return m * cp * dT / P_abs


# ---------------------------------------------------------------------------
# 5. DINAMICA: IL MODELLO UNICO
# ---------------------------------------------------------------------------
def coppia_limite(P_ass, l=0.0, s=1.0, lam=LAMBDA):
    """tau = (l+s) P_ass/omega_ott.  Unica formula di coppia della V6."""
    return (l + s) * P_ass * lam / (2 * math.pi * C.c)


def forza_limite(I, sigma_pr):
    """F = I sigma_pr/c.  Unica formula di forza della V6."""
    return I * sigma_pr / C.c


def pe_r(tau, kT=KT):
    """Numero di Peclet rotazionale.  Identita': Pe_r = tau/kT, indipendente da beta.

    Pe_r = omega_ss/D_r con D_r = kT/beta, quindi Pe_r = tau/kT.  La V5
    riportava due valori diversi (0.0088 e 0.00061) perche' usava due beta
    diversi nei due punti: l'identita' e' anche un test di coerenza.
    """
    return tau / kT


def omega_stazionario(tau, beta):
    return tau / beta


def omega_rottura(sy, rho, r):
    """Velocita' angolare a cui la tensione centrifuga rho omega^2 r^2 = sy."""
    return math.sqrt(sy / rho) / r


def tempo_rottura(I, tau, sy, rho, r):
    """Tempo per raggiungere omega_rottura in assenza di attrito (vuoto)."""
    if tau <= 0:
        return float("inf")
    return I * omega_rottura(sy, rho, r) / tau


def tempi(I, beta):
    """Tempo inerziale I/beta e tempo di smorzamento."""
    return {"I_su_beta": I / beta if beta > 0 else float("inf")}


def numero_peclet_rot(tau):
    return pe_r(tau)


def regime_fluidodinamico(omega, r, mu=ETA_W, rho=RHO_W):
    re = rho * abs(omega) * r ** 2 / mu
    if re < 1e-3:
        return re, "Stokes"
    if re < 1.0:
        return re, "Stokes (inerzia trascurabile)"
    if re < 3e5:
        return re, "laminare con inerzia"
    return re, "turbolento"


# ---------------------------------------------------------------------------
# 6. INVENTARIO DELLA GEOMETRIA V5 (per il confronto a pari massa)
# ---------------------------------------------------------------------------
def perimetro_ellisse(a, b):
    """Formula di Ramanujan (errore relativo < 1e-9)."""
    return math.pi * (3 * (a + b) - math.sqrt((3 * a + b) * (a + 3 * b)))


LUNGHEZZA_FILI_NM = (sum(perimetro_ellisse(a, b) for a, b in SPIN_ORBITS)
                     + 2 * math.pi * 327.5
                     + 2 * perimetro_ellisse(*LENS))

# inventario per gruppi: (lunghezza in nm, raggio in m)
FILI = (
    (sum(perimetro_ellisse(a, b) for a, b in SPIN_ORBITS), WIRE_ORB),
    (2 * math.pi * 327.5, WIRE_RIM),
    (2 * perimetro_ellisse(*LENS), WIRE_LENS),
)

AREA_PRO_FILI = sum(L * 1e-9 * 2 * a for L, a in FILI)          # m^2, proiettata
AREA_SUP_FILI = sum(L * 1e-9 * 2 * math.pi * a for L, a in FILI)  # m^2, laterale
VOLUME_FILI = sum(L * 1e-9 * math.pi * a ** 2 for L, a in FILI)   # m^3
VOLUME_NUCLEO = 4.0 / 3.0 * math.pi * CORE_R ** 3
VOLUME_TELAIO = VOLUME_FILI + VOLUME_NUCLEO
AREA_NUCLEO = 4 * math.pi * CORE_R ** 2


def area_radiante():
    """Superficie metallica effettivamente esposta, non 4 pi R_out^2.

    E' la grandezza che entra nel bilancio radiativo: usare l'inviluppo
    sferico di 370 nm sopravvaluta l'area di un fattore 6 e sottostima la
    temperatura di un fattore 1.6.
    """
    return AREA_SUP_FILI + AREA_NUCLEO


def massa_telaio(rho=8960.0):
    return rho * VOLUME_TELAIO


def momento_inerzia_telaio(rho=8960.0):
    """Momento d'inerzia rispetto all'asse del disco.

    Anelli: m R^2, con la massa per unita' di lunghezza del gruppo corretto.
    Fili delle orbite ellittiche assimilati ad anelli di raggio pari alla media
    dei semiassi (errore < 3%).  Nucleo: (2/5) m r^2.
    """
    I = 0.0
    for a, b in SPIN_ORBITS:
        L = perimetro_ellisse(a, b) * 1e-9
        Rm = 0.5 * (a + b) * 1e-9
        I += rho * math.pi * WIRE_ORB ** 2 * L * Rm ** 2
    I += rho * math.pi * WIRE_RIM ** 2 * (2 * math.pi * 327.5e-9) * (327.5e-9) ** 2
    for _ in range(2):
        L = perimetro_ellisse(*LENS) * 1e-9
        I += rho * math.pi * WIRE_LENS ** 2 * L * (0.5 * (LENS[0] + LENS[1]) * 1e-9) ** 2
    I += 0.4 * rho * VOLUME_NUCLEO * CORE_R ** 2
    return I


def raggi_equivalenti():
    """Le tre definizioni di raggio che la V5 confondeva.

    r_volume_materiale : sfera con lo stesso volume di materiale
    r_inviluppo        : sfera con lo stesso volume dell'inviluppo cilindrico
    r_idrodinamico     : sfera con lo stesso beta (quella giusta per D e D_r)
    """
    v_mat = VOLUME_TELAIO
    r_mat = (3 * v_mat / (4 * math.pi)) ** (1 / 3)
    v_env = math.pi * R_OUT ** 2 * (2 * R_OUT)      # inviluppo cilindrico
    r_env = (3 * v_env / (4 * math.pi)) ** (1 / 3)
    b_lo = beta_anello(R_OUT, WIRE_ORB)
    b_hi = beta_disco(R_OUT)
    r_idro = ((b_lo + b_hi) / 2 / (8 * math.pi * ETA_W)) ** (1 / 3)
    return {"r_volume_materiale_nm": r_mat * 1e9,
            "r_inviluppo_nm": r_env * 1e9,
            "r_idrodinamico_nm": r_idro * 1e9,
            "r_eff_v5_nm": R_EFF_V5 * 1e9}


# ---------------------------------------------------------------------------
# 7. SELF-TEST
# ---------------------------------------------------------------------------
def _test():
    print("=== SELF-TEST V6 ===")
    m_cu = complex(1.09, 2.53)
    print("-- cilindro: serie n>=1 contro Rayleigh")
    for ra in (0.5e-9, 1.25e-9, 5e-9, 20e-9):
        x = 2 * math.pi * ra / LAMBDA
        e = cilindro_esatto(m_cu, x)
        r = cilindro_rayleigh(m_cu, x)
        print(f"  a={ra*1e9:5.2f} nm x={x:.5f}  Q esatto={e['Q_ext']:.5e}  "
              f"Rayleigh={r['Q_ext']:.5e}  scarto={e['Q_ext']/r['Q_ext']-1:+.2e}")
    print("-- sfera: Mie esatta")
    for ra in (10e-9, 20e-9, 40e-9, 100e-9, 327.5e-9):
        s = sezione_sfera(m_cu, ra)
        print(f"  r={ra*1e9:6.1f} nm  x={s['x']:.3f}  Q_ext={s['Q_ext']:.3f}  "
              f"Q_abs={s['Q_abs']:.3f}  Q_pr={s['Q_pr']:.3f}  [{s['regime']}]")
    print(f"-- lunghezza fili = {LUNGHEZZA_FILI_NM:.1f} nm")
    print(f"-- volume telaio  = {VOLUME_TELAIO*1e27:.3f} nm^3 "
          f"(fili {VOLUME_FILI*1e27:.2f}, nucleo {VOLUME_NUCLEO*1e27:.2f})")
    print(f"-- massa telaio   = {massa_telaio():.4e} kg   (V5 usava {M_V5:.4e} kg)")
    print(f"-- I telaio       = {momento_inerzia_telaio():.4e} kg m^2  (V5: {I1_V5:.4e})")
    for k, v in raggi_equivalenti().items():
        print(f"   {k:24s} = {v:.2f} nm")
    b_lo, b_hi = beta_anello(R_OUT, WIRE_ORB), beta_disco(R_OUT)
    print(f"-- beta intervallo = [{b_lo:.3e}, {b_hi:.3e}]  V5 = {BETA_V5:.3e}  "
          f"rapporto V5/basso = {BETA_V5/b_lo:.3f}")
    # sezioni della configurazione di riferimento
    core = sezione_sfera(m_cu, CORE_R)
    s_pr_f = sum(sezione_filo(m_cu, a, L * 1e-9)["sigma_pr"] for L, a in FILI)
    s_abs_f = sum(sezione_filo(m_cu, a, L * 1e-9)["sigma_abs"] for L, a in FILI)
    filo = {"sigma_pr": s_pr_f, "sigma_abs": s_abs_f}
    s_pr = core["sigma_pr"] + filo["sigma_pr"]
    s_abs = core["sigma_abs"] + filo["sigma_abs"]
    I0 = P_LASER / A_BEAM
    print(f"-- sigma_pr = {s_pr:.4e} m^2 (V5 {SIGMA_PR_V5:.4e}), "
          f"sigma_abs = {s_abs:.4e}, eta_int = {s_pr/A_BEAM*100:.4f}%")
    print(f"-- F = {forza_limite(I0, s_pr):.4e} N  (V5 {F_V5:.4e})")
    P_abs = I0 * s_abs
    tau = coppia_limite(P_abs)
    print(f"-- P_inc={P_LASER*1e6:.2f} uW  P_int={I0*s_pr*1e9:.3f} nW  "
          f"P_abs={P_abs*1e9:.3f} nW")
    print(f"-- tau = {tau:.4e} N m  (V5 2.610e-24)   Pe_r = {pe_r(tau):.3e}")
    print(f"-- omega_ss intervallo = [{tau/b_hi:.3e}, {tau/b_lo:.3e}] rad/s "
          f"(V5 {tau/BETA_V5:.3e})")
    print(f"-- identita' Pe_r: {pe_r(tau):.6e} = {(tau/b_lo)/(KT/b_lo):.6e}")


if __name__ == "__main__":
    _test()
