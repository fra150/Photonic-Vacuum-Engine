"""
MOTORE FOTONICO V5 - NUCLEO FISICO MULTI-SCALA
Autore: Dr. Francesco Bulla / Sistema GNRX

Questo modulo implementa, in un unico posto e senza dipendenze dal testo del
preprint, le grandezze che servono a:

  1. dimostrare il doppio giro (720 gradi) sullo Spin e il rapporto razionale
     45/376 con il Contro-spin;
  2. introdurre il termine dissipativo tau_diss = -beta*omega con l'attrito
     rotazionale calcolato sulla geometria reale (telaio di fili, non disco pieno);
  3. chiudere esattamente il ciclo per progetto, agendo sulla massa di
     bilanciamento del Contro-spin;
  4. sostituire le efficienze dichiarate eta_s = 0.70 e eta_c = 0.85 con un
     calcolo di scattering (Mie per le sfere, Rayleigh cilindrico per i fili);
  5. definire il protocollo di verifica sperimentale (coppia ottica e
     periodicita' a pi greco, senza firma a 4 pi greco).

Le stesse funzioni vengono valutate su quattro scale:
  - 10^-9 m  (dispositivo da 1 nm, k = 10^-3 rispetto alla V4)
  - 240 nm   (geometria V4 di riferimento, r_eff = 239.93 nm)
  - 1 m      (pseudo-motore macroscopico, k = 10^6 rispetto alla V4)
  - 1 GHz    (geometria V4 con rotazione a 10^9 Hz)

Uso:
    python v5_nucleo.py
"""

import json
import math
import os

import numpy as np
from scipy import constants as const
from scipy.special import spherical_jn, spherical_yn

OUT = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# 0. COSTANTI E PARAMETRI DI RIFERIMENTO (V4)
# ---------------------------------------------------------------------------
T = 310.0                      # K, temperatura corporea
ETA = 1.2e-3                   # Pa s, viscosita' dinamica del fluido biologico
RHO = 998.2                    # kg/m^3, acqua a 20 C
NU = ETA / RHO                 # m^2/s, viscosita' cinematica
P_LASER = 1.0e-6               # W, potenza ottica incidente
LAMBDA = 532e-9                # m, luce verde
C_LIGHT = const.c              # m/s
KT = const.k * T               # J

UA = 100e-9                    # m, 1 unita' di scena = 100 nm
R_EFF_V4 = math.sqrt(5.7567) * UA          # 239.93 nm, raggio idrodinamico V4
M_V4 = 1.2463630573655231e-18 + 4.487710753596799e-18   # kg
R_DISC = 3.275 * UA            # 327.5 nm, raggio esterno del disco
W0 = 500e-9                    # m, raggio del fascio focalizzato (NA 0.9 a 532 nm)
A_BEAM = math.pi * W0 ** 2     # m^2, area del fascio

# Efficienze dichiarate nella V4 (stima, non calcolo)
ETA_S_DECL = 0.70
ETA_C_DECL = 0.85
ALPHA_DECL = math.radians(15.1)   # angolo di incidenza che riproduce F_z = 5.074 fN
F_DECL = (P_LASER / C_LIGHT) * (ETA_S_DECL + ETA_C_DECL * math.cos(ALPHA_DECL))

# Momenti d'inerzia della scena V4 (kg m^2)
I1_V4 = 3.092e-32
I2_V4 = 2.583e-31

# Geometria dei fili, letta dalla scena: tre orbite ellittiche dello Spin,
# bordo del disco del Contro-spin, due lenti verticali, nucleo sferico.
SPIN_ORBITS = ((220.0, 110.0), (295.0, 152.0), (370.0, 195.0))   # nm, semiassi
LENS = (105.0, 245.0)                                            # nm, semiassi
CORE_R = 40.0                                                    # nm
WIRE_R = 1.25                                                    # nm


def ellipse_perimeter(a, b):
    """Perimetro di un'ellisse con la formula di Ramanujan (errore < 1e-4)."""
    return math.pi * (3 * (a + b) - math.sqrt((3 * a + b) * (a + 3 * b)))


WIRE_LENGTH_NM = (
    sum(ellipse_perimeter(a, b) for a, b in SPIN_ORBITS)
    + 2 * math.pi * 327.5
    + 2 * ellipse_perimeter(*LENS)
)
WIRE_AREA = WIRE_LENGTH_NM * 1e-9 * 2 * WIRE_R * 1e-9    # m^2, area proiettata
CORE_AREA = math.pi * (CORE_R * 1e-9) ** 2               # m^2, sezione del nucleo


# ---------------------------------------------------------------------------
# 1. TRASPORTO A UNA SCALA QUALSIASI
# ---------------------------------------------------------------------------
def trasporto(r_eff, m_tot, F, rho=RHO, eta=ETA, kT=KT):
    """Grandezze di trasporto per una sfera equivalente di raggio r_eff."""
    drag = 6 * math.pi * eta * r_eff
    v = F / drag
    D = kT / drag
    D_rot = kT / (8 * math.pi * eta * r_eff ** 3)
    tau_rot = 1.0 / (2 * D_rot)
    return {
        "r_eff_nm": r_eff * 1e9,
        "massa_kg": m_tot,
        "drag_Ns_m": drag,
        "v_um_s": v * 1e6,
        "D_m2_s": D,
        "Re": rho * F / (3 * math.pi * eta ** 2),
        "Pe": F * r_eff / kT,
        "v_termica_m_s": math.sqrt(kT / m_tot),
        "D_rot_1_s": D_rot,
        "tau_rot_s": tau_rot,
        "l_pers_nm": v * tau_rot * 1e9,
        "t_cross_s": 2 * D / v ** 2,
        "L_crit_m": D / v,
        "Deff_over_D": 1 + v ** 2 / (6 * D_rot * D),
        "forza_N": F,
    }


def forza_ottica(k=1.0, eta_s=ETA_S_DECL, eta_c=ETA_C_DECL, alpha=ALPHA_DECL):
    """Spinta lungo z a potenza e waist fissati.

    Per k < 1 il dispositivo e' piu' piccolo del fascio e intercetta la
    frazione k^2 della potenza; per k >= 1 intercetta tutto il fascio e la
    spinta satura al valore (P/c)(eta_s + eta_c cos alpha).
    """
    return (P_LASER / C_LIGHT) * (eta_s + eta_c * math.cos(alpha)) * min(1.0, k ** 2)


def coppia_oam(F, lam=LAMBDA):
    """Coppia ceduta dai fotoni: tau = P_ass/omega_ott = F lambda/(2 pi).

    E' l'unico modello di momento angolare usato nel documento. Non dipende
    dal raggio del rotore: dipende solo dalla potenza effettivamente
    intercettata, cioe' dalla forza di scattering.
    """
    return F * lam / (2.0 * math.pi)


def scala_geometrica(nome, k):
    """Costruisce il caso fisico ottenuto scalando la V4 del fattore k.

    La coppia ottica non e' piu' F_dichiarata * R: e' il momento angolare
    ceduto dai fotoni intercettati, tau = F lambda/(2 pi), con F calcolata
    dal calcolo di scattering (Mie per il nucleo, Rayleigh per i fili).
    La forza dichiarata resta nel dizionario come termine di confronto.
    """
    r_eff = R_EFF_V4 * k
    m = M_V4 * k ** 3
    F_dich = forza_ottica(k)
    sc = scattering_globale(k)
    F = sc["F_z_N"]
    d = trasporto(r_eff, m, F)
    F_mie_adattato = (P_LASER / area_diffrazione()) * sc["sigma_effettiva_m2"] / C_LIGHT
    beta_t = beta_telaio(R_DISC * k)
    tau = coppia_oam(F)
    tau_adattato = coppia_oam(F_mie_adattato)
    omega = tau / beta_t
    d.update({
        "nome": nome,
        "k": k,
        "ingombro_nm": 1900.0 * k,
        "R_disco_m": R_DISC * k,
        "I1_kgm2": I1_V4 * k ** 5,
        "I2_kgm2": I2_V4 * k ** 5,
        "peso_N": m * 9.80665,
        "forza_dichiarata_N": F_dich,
        "spinta_peso": F / (m * 9.80665),
        "tau_ottica_Nm": tau,
        "tau_su_kT": tau / KT,
        "beta_disco_Nms": beta_disco(R_DISC * k),
        "beta_telaio_Nms": beta_t,
        "omega_ss_rad_s": omega,
        "giri_al_secondo": omega / (2 * math.pi),
        "giri_per_anno": omega / (2 * math.pi) * 3.15576e7,
        "energia_rotazionale_J": 0.5 * I1_V4 * k ** 5 * omega ** 2,
        "energia_su_kT": 0.5 * I1_V4 * k ** 5 * omega ** 2 / KT,
        "forza_mie_N": F,
        "forza_mie_fuoco_adattato_N": F_mie_adattato,
        "Pe_con_scattering": tau / KT,
        "Pe_fuoco_adattato": tau_adattato / KT,
        "potenza_per_spinta_peso_W": P_LASER / (F / (m * 9.80665)),
    })
    return d


# ---------------------------------------------------------------------------
# 2. ATTRITO ROTAZIONALE SULLA GEOMETRIA REALE
# ---------------------------------------------------------------------------
def fill_factor():
    """Frazione di area bagnata del telaio rispetto al disco pieno.

    Il disco della V4 e' un telaio di fili di raggio 1.25 nm, non un disco
    pieno: l'attrito rotazionale scala con l'area effettivamente bagnata.
    """
    a_disc = math.pi * (R_DISC * 1e9) ** 2          # nm^2
    return WIRE_AREA * 1e18 / a_disc


FILL_FRAME = fill_factor()


def beta_disco(R=R_DISC, eta=ETA):
    """Attrito rotazionale di un disco pieno di raggio R (Kogan)."""
    return (32.0 / 3.0) * eta * R ** 3


def beta_telaio(R=R_DISC, eta=ETA):
    """Attrito rotazionale del telaio di fili (disco pieno x frazione bagnata)."""
    return beta_disco(R, eta) * FILL_FRAME


# ---------------------------------------------------------------------------
# 3. COPPIA OTTICA DIPENDENTE DALL'ANGOLO E DINAMICA
# ---------------------------------------------------------------------------
def coppia_ottica(theta, tau_max, eta_s=ETA_S_DECL, eta_c=ETA_C_DECL):
    """Coppia a due lobi: tau(theta) = tau_max * (eta_s + eta_c cos 2 theta)/(eta_s+eta_c)."""
    return tau_max * (eta_s + eta_c * math.cos(2 * theta)) / (eta_s + eta_c)


def dinamica_rotore(I, beta, tau_max, omega0=0.0, theta0=0.0, t_end=60.0, dt=1e-3,
                    rumore=False, D_rot=0.0):
    """Integra I*theta'' = tau(theta) - beta*theta' [+ rumore browniano]."""
    rng = np.random.default_rng(20261025)
    n = int(t_end / dt)
    th = np.empty(n)
    om = np.empty(n)
    th[0], om[0] = theta0, omega0
    for i in range(1, n):
        acc = (coppia_ottica(th[i - 1], tau_max) - beta * om[i - 1]) / I
        om[i] = om[i - 1] + acc * dt
        th[i] = th[i - 1] + om[i] * dt
        if rumore and D_rot > 0:
            om[i] += math.sqrt(2 * D_rot * dt) * rng.standard_normal()
    return np.linspace(0, t_end, n), th, om


# ---------------------------------------------------------------------------
# 4. DOPPIO GIRO (720 GRADI), SU(2) E CHIUSURA RAZIONALE
# ---------------------------------------------------------------------------
RATIO_NUM, RATIO_DEN = 45, 376          # omega2/omega1 = -45/376 = -0.119680...
RATIO = RATIO_NUM / RATIO_DEN

# Un ciclo di animazione: Spin 720 gradi (due giri), Contro-spin 720*45/376
SPIN_DEG_CYCLE = 720.0
CONTRO_DEG_CYCLE = SPIN_DEG_CYCLE * RATIO     # 86.1702 gradi
CICLI_CHIUSURA = RATIO_DEN // 2               # 188 cicli -> 376 giri / 45 giri
FRAME_RATE = 24
FRAME_END = 240


def matrice_su2(theta, asse="z"):
    """U(theta) = exp(-i theta sigma/2): doppia copertura di SO(3)."""
    c, s = math.cos(theta / 2), math.sin(theta / 2)
    if asse == "z":
        return np.array([[c - 1j * s, 0], [0, c + 1j * s]])
    if asse == "x":
        return np.array([[c, -1j * s], [-1j * s, c]])
    return np.array([[c, -s], [s, c]])


def chiusura_esatta(I1=I1_V4, I2=I2_V4, r_bal=42.0e-9):
    """Massa di bilanciamento che rende I2/I1 = 376/45 esattamente."""
    target = RATIO_DEN / RATIO_NUM
    I2_target = target * I1
    dm = I2_target - I2
    m_bal = dm / r_bal ** 2
    m_cu = 63.546 * const.atomic_mass
    return {
        "I1": I1,
        "I2": I2,
        "I2_su_I1": I2 / I1,
        "target_376_45": target,
        "scarto_relativo": I2 / I1 / target - 1,
        "I2_target": I2_target,
        "delta_I2": dm,
        "r_bilanciamento_m": r_bal,
        "massa_bilanciamento_kg": m_bal,
        "atomi_Cu_equivalenti": m_bal / m_cu,
    }


# ---------------------------------------------------------------------------
# 5. SCATTERING: MIE (SFERE) E RAYLEIGH (FILI)
# ---------------------------------------------------------------------------
# Indici di rifrazione del rame (Johnson & Christy, interpolati)
CU_NK = {
    532e-9: 1.09 + 2.53j,
    1064e-9: 0.36 + 6.60j,
    1550e-9: 0.55 + 8.90j,
}

X_GEOMETRICA = 40.0   # oltre questo parametro di dimensione si usa l'ottica geometrica


def _log_derivata(nmax, z):
    """Derivata logaritmica D_n(z) per ricorrenza discendente (Bohren-Huffman)."""
    D = np.zeros(nmax + 2, dtype=complex)
    for n in range(nmax, 0, -1):
        D[n - 1] = n / z - 1.0 / (D[n] + n / z)
    return D


def mie_sfera(m, x, nmax=None):
    """Coefficienti di Mie per una sfera: m = n + ik relativo, x = 2 pi a / lambda.

    Per x > X_GEOMETRICA il calcolo esatto perde senso (servirebbero decine di
    migliaia di termini): si usa il limite di ottica geometrica, Q_pr -> 2.
    """
    if x > X_GEOMETRICA:
        return {"Q_ext": 2.0, "Q_sca": 1.0, "Q_abs": 1.0, "g": 0.85, "Q_pr": 1.15,
                "regime": "ottica geometrica"}
    if nmax is None:
        nmax = int(2 + x + 4 * x ** (1.0 / 3.0)) + 8
    nmax = max(min(nmax, 200), 12)
    mx = m * x
    D = _log_derivata(nmax, mx)
    psi = lambda n, z: z * spherical_jn(n, z)
    xi = lambda n, z: z * (spherical_jn(n, z) + 1j * spherical_yn(n, z))
    a = np.zeros(nmax + 1, dtype=complex)
    b = np.zeros(nmax + 1, dtype=complex)
    for n in range(1, nmax + 1):
        a[n] = ((D[n] / m + n / x) * psi(n, x) - psi(n - 1, x)) / \
               ((D[n] / m + n / x) * xi(n, x) - xi(n - 1, x))
        b[n] = ((D[n] * m + n / x) * psi(n, x) - psi(n - 1, x)) / \
               ((D[n] * m + n / x) * xi(n, x) - xi(n - 1, x))
    n = np.arange(1, nmax + 1)
    q_ext = (2.0 / x ** 2) * np.sum((2 * n + 1) * (a[1:] + b[1:]).real)
    q_sca = (2.0 / x ** 2) * np.sum((2 * n + 1) * (np.abs(a[1:]) ** 2 + np.abs(b[1:]) ** 2))
    g = (4.0 / (x ** 2 * q_sca)) * np.sum(
        (n * (n + 1) / (n + 0.5)) * (a[1:] * np.conj(a[:-1]) + b[1:] * np.conj(b[:-1])).real
        + ((2 * n + 1) / (n * (n + 1))) * (a[1:] * np.conj(b[1:])).real)
    q_pr = q_ext - g * q_sca
    return {"Q_ext": q_ext, "Q_sca": q_sca, "Q_abs": q_ext - q_sca, "g": g, "Q_pr": q_pr,
            "regime": "Mie esatta"}


def q_pr_cilindro(m, x):
    """Estinzione di un cilindro sottile in regime di Rayleigh (x << 1).

    Per un cilindro di raggio a investito perpendicolarmente:
        Q_ext = pi*x*Im[(m^2-1)/(m^2+1)],  Q_sca ~ O(x^3), Q_pr ~ Q_ext.
    """
    m2 = m ** 2
    pol = (m2 - 1) / (m2 + 1)
    if x > 1.0:
        return {"Q_ext": 2.0, "Q_sca": 0.9, "Q_abs": 1.1,
                "Q_pr": 1.1, "pol": pol, "regime": "ottica geometrica"}
    q_ext = math.pi * x * pol.imag
    q_sca = (math.pi ** 3 / 4.0) * x ** 3 * abs(pol) ** 2
    return {"Q_ext": q_ext, "Q_sca": q_sca, "Q_abs": q_ext - q_sca,
            "Q_pr": q_ext - q_sca, "pol": pol, "regime": "Rayleigh cilindrico"}


def scattering_globale(k=1.0, lam=LAMBDA):
    """Sezione di pressione di radiazione della struttura e forza risultante.

    I fili contribuiscono con la sezione cilindrica, il nucleo con la sua
    sezione geometrica per Q_pr di sfera.  La forza e' I*sigma_pr/c, con
    I = P / A_beam; il fattore geometrico k^2 scala le aree.
    """
    m_rel = CU_NK[min(CU_NK, key=lambda w: abs(w - lam))]
    x_wire = 2 * math.pi * (WIRE_R * 1e-9 * k) / lam
    x_core = 2 * math.pi * (CORE_R * 1e-9 * k) / lam
    cyl = q_pr_cilindro(m_rel, x_wire)
    sph = mie_sfera(m_rel, x_core) if x_core > 1e-3 else {"Q_pr": 0.0, "regime": "sotto-Raleigh"}
    sigma_wire = WIRE_AREA * k ** 2 * max(cyl["Q_pr"], 0.0)
    sigma_core = CORE_AREA * k ** 2 * max(sph["Q_pr"], 0.0)
    sigma_grezza = sigma_wire + sigma_core
    # La struttura non puo' intercettare piu' luce di quella che il fascio porta:
    # se e' piu' grande del fascio la sezione efficace satura all'area del fascio.
    sigma = min(sigma_grezza, A_BEAM)
    intensita = P_LASER / A_BEAM
    F = intensita * sigma / C_LIGHT
    a_min = area_diffrazione(lam)
    F_min_area = (P_LASER / a_min) * min(sigma_grezza, a_min) / C_LIGHT
    return {
        "lambda_nm": lam * 1e9,
        "n_Cu": m_rel.real, "k_Cu": m_rel.imag,
        "x_filo": x_wire, "x_nucleo": x_core,
        "regime_filo": cyl["regime"], "regime_nucleo": sph["regime"],
        "Q_pr_filo": cyl["Q_pr"], "Q_pr_nucleo": sph["Q_pr"],
        "sigma_filo_m2": sigma_wire, "sigma_nucleo_m2": sigma_core,
        "sigma_pr_m2": sigma_grezza,
        "sigma_effettiva_m2": sigma,
        "intensita_W_m2": intensita,
        "frazione_fascio": sigma / A_BEAM,
        "area_fascio_m2": A_BEAM,
        "area_min_diffrazione_m2": a_min,
        "F_z_N": F,
        "F_z_fuoco_limitato_N": F_min_area,
        "F_z_dichiarata_N": F_DECL,
        "efficienza_rispetto_dichiarata": F / F_DECL,
        "scarto_fuoco_limitato": F_min_area / F_DECL,
    }


def area_diffrazione(lam=LAMBDA):
    """Area minima di un fuoco limitato dalla diffrazione: (lambda/2)^2."""
    return (lam / 2.0) ** 2


def potenza_per_peclet(pe=1.0, r_eff=R_EFF_V4, area=None):
    """Potenza ottica necessaria per Pe = 1 a parita' di sezione di scattering."""
    area = A_BEAM if area is None else area
    sigma = scattering_globale()["sigma_effettiva_m2"]
    F_necessaria = pe * KT / r_eff
    return F_necessaria * area * C_LIGHT / sigma


# ---------------------------------------------------------------------------
# 6. PROTOCOLLO DI VERIFICA SPERIMENTALE
# ---------------------------------------------------------------------------
def soglia_rotazionale(I=I1_V4, beta=None, tau=None, R=R_DISC, r_eff=R_EFF_V4):
    """Confronto tra rotazione ottica e diffusione rotazionale.

    Il criterio corretto non e' la velocita' angolare termica istantanea
    omega_th = sqrt(kT/I), che vive sulla scala inerziale I/beta, ma il numero
    di Peclet rotazionale Pe_r = omega_ss / D_r: se Pe_r > 1 la rotazione
    ottica sopravvive alla diffusione browniana dell'orientazione.
    """
    beta = beta_telaio() if beta is None else beta
    sc = scattering_globale()
    F_reale = sc["F_z_N"]
    frazione = sc["frazione_fascio"]
    omega_ott = 2.0 * math.pi * C_LIGHT / LAMBDA
    tau = coppia_oam(F_reale) if tau is None else tau
    D_rot = KT / (8 * math.pi * ETA * r_eff ** 3)
    omega_ss = tau / beta
    tau_need = beta * D_rot
    F_need = tau_need / (LAMBDA / (2 * math.pi))
    return {
        "I_kgm2": I,
        "beta_Nms": beta,
        "tau_Nm": tau,
        "omega_ss_rad_s": omega_ss,
        "D_rot_1_s": D_rot,
        "Pe_rot": omega_ss / D_rot,
        "Pe_tau_su_kT": tau / KT,
        "tempo_inerziale_s": I / beta,
        "angolo_ottico_per_tau_rot_rad": omega_ss / (2 * D_rot),
        "spread_termico_per_tau_rot_rad": 1.0,
        "omega_termica_rad_s": math.sqrt(KT / I),
        "tau_necessaria_Nm": tau_need,
        "forza_necessaria_N": F_need,
        "potenza_pe_tau_1_W": KT * omega_ott / frazione,
        "potenza_pe_rot_1_W": P_LASER * F_need / F_reale,
        "coerenza": "deterministica" if omega_ss > D_rot else "termica",
    }


def protocollo_sperimentale(tau, beta, I, kT=KT):
    """Numeri attesi nei due esperimenti discriminanti."""
    omega = tau / beta
    tempo_equipartizione = I * omega ** 2 / kT if omega else float("inf")
    return {
        "coppia_Nm": tau,
        "beta_Nms": beta,
        "omega_atteso_rad_s": omega,
        "giri_al_secondo": omega / (2 * math.pi),
        "energia_rotazionale_J": 0.5 * I * omega ** 2,
        "energia_su_kT": 0.5 * I * omega ** 2 / kT,
        "soglia_coppia_kT": kT,
        "coppia_su_kT": tau / kT,
        "periodo_coppia_gradi": 180.0,
        "nota_periodicita": "tau(theta) = 0.5 tau_max (1 + cos 2 theta): il "
                            "segnale si ripete ogni pi greco, non ogni 4 pi.",
    }


# ---------------------------------------------------------------------------
# 7. SCALE
# ---------------------------------------------------------------------------
SCALE = {
    "1_nm": 1e-3,
    "v4_240_nm": 1.0,
    "1_m": 1e6,
}


def tutti_i_casi():
    casi = {}
    for nome, k in SCALE.items():
        casi[nome] = scala_geometrica(nome, k)
    return casi


def caso_1ghz():
    """Geometria V4 con rotazione a 1 GHz: conseguenze meccaniche e di regime."""
    omega = 2 * math.pi * 1e9
    R = R_DISC
    v_tang = omega * R
    acc = omega ** 2 * R
    rho_cu = 8960.0
    sigma = rho_cu * v_tang ** 2             # Pa, tensione centrifuga rho v^2
    snervamento = 2.0e8                      # Pa, rame ricotto
    c_suono = 1480.0                         # m/s, acqua
    v_y = math.sqrt(snervamento / rho_cu)    # velocita' tangenziale ammessa
    f_y = v_y / (2 * math.pi * R)
    f_c = c_suono / (2 * math.pi * R)
    beta = beta_telaio()
    vol = 4.0 / 3.0 * math.pi * R_EFF_V4 ** 3
    cp = 4186.0
    riscaldamento = beta * omega ** 2 / (RHO * cp * vol)
    return {
        "frequenza_Hz": 1e9,
        "omega_rad_s": omega,
        "v_tangenziale_m_s": v_tang,
        "Mach_in_acqua": v_tang / c_suono,
        "accelerazione_m_s2": acc,
        "accelerazione_in_g": acc / 9.80665,
        "tensione_centrifuga_Pa": sigma,
        "snervamento_Cu_Pa": snervamento,
        "margine_strutturale": snervamento / sigma,
        "velocita_limite_snervamento_m_s": v_y,
        "frequenza_limite_snervamento_Hz": f_y,
        "frequenza_limite_cavitazione_Hz": f_c,
        "Re_rotazione": v_tang * 2 * R / NU,
        "coppia_per_attrito_Nm": beta * omega,
        "potenza_dissipata_W": beta * omega ** 2,
        "riscaldamento_K_al_secondo": riscaldamento,
        "verdetto": "impossibile" if (sigma > snervamento or v_tang > c_suono) else "plausibile",
    }


if __name__ == "__main__":
    print("=" * 68)
    print("MOTORE FOTONICO V5 - CONTROLLO DEL NUCLEO FISICO")
    print("=" * 68)
    print("Forza dichiarata (eta_s=0.70, eta_c=0.85): %.3e N" % F_DECL)
    print("Frazione bagnata del telaio: %.4f" % FILL_FRAME)
    print("beta disco pieno: %.3e N m s | beta telaio: %.3e N m s"
          % (beta_disco(), beta_telaio()))
    for nome, d in tutti_i_casi().items():
        print("\n[%s]  k = %.3g" % (nome, d["k"]))
        for chiave in ("ingombro_nm", "r_eff_nm", "massa_kg", "drag_Ns_m", "v_um_s",
                       "D_m2_s", "Re", "Pe", "D_rot_1_s", "tau_rot_s", "l_pers_nm",
                       "t_cross_s", "L_crit_m", "Deff_over_D", "spinta_peso",
                       "tau_su_kT", "omega_ss_rad_s", "giri_per_anno",
                       "energia_su_kT", "forza_mie_N", "Pe_con_scattering",
                       "Pe_fuoco_adattato"):
            print("   %-22s % .4e" % (chiave, d[chiave]))
    print("\nCoppia/inerzia: SU(2) U(2pi) = -I ?", np.allclose(matrice_su2(2 * math.pi), -np.eye(2)))
    print("                 U(4pi) = +I ?", np.allclose(matrice_su2(4 * math.pi), np.eye(2)))
    print("Contro-spin per ciclo: %.4f gradi; chiusura in %d cicli"
          % (CONTRO_DEG_CYCLE, CICLI_CHIUSURA))
    print("\nChiusura esatta:", json.dumps(chiusura_esatta(), indent=2, default=float))
    print("\nScattering 532 nm:", json.dumps(scattering_globale(), indent=2, default=float))
    print("Potenza per Pe=1 (waist 500 nm): %.3e W" % potenza_per_peclet())
    print("Potenza per Pe=1 (fuoco limitato): %.3e W"
          % potenza_per_peclet(area=area_diffrazione()))
    print("\n1 GHz:", json.dumps(caso_1ghz(), indent=2, default=float))
