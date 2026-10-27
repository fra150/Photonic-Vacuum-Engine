"""
MOTORE FOTONICO V4 - MODULO DI SIMULAZIONE FISICA & FORMALIZZAZIONE
Autore: Dr. Francesco Bulla / Sistema GNRX
Funzionalità:
 1. Definizione rigorosa delle scale dimensionali nanometriche.
 2. Formalizzazione algebrica Spin 720° nello spazio di Hilbert SU(2).
 3. Bilancio di conservazione del momento angolare (Torque Cancellation).
 4. Calcolo della spinta fotonica tramite tensore di Maxwell e vettore di Poynting.
 5. Simulazione dinamica delle traiettorie dei flussi di fotoni.
 6. Attrito viscoso di Stokes e moto browniano: confronto vuoto / fluido biologico.

Uso:
    python motore_fotonico_fisica.py
Gli output (grafici PNG) vengono salvati nella stessa cartella di questo file.
"""

import os

import numpy as np
import scipy.constants as const
from scipy.integrate import solve_ivp
import sympy as sp
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT_DIR = os.path.dirname(os.path.abspath(__file__))


def out_path(name):
    return os.path.join(OUT_DIR, name)


# ==============================================================================
# PUNTO 3: DEFINIZIONE DELLE SCALE DIMENSIONALI (Costanti e Parametri Nanometrici)
# ==============================================================================
print("=" * 60)
print("1. DEFINIZIONE DELLE SCALE DIMENSIONALI")
print("=" * 60)

SCALE_PARAMS = {
    # Scale fisiche del DNA origami e rivestimento metallico
    "r_dna_core": 5.0e-9,           # Raggio del nucleo centrale (5 nm)
    "r_ring_primary": 35.0e-9,      # Raggio degli anelli di Spin (35 nm)
    "r_ring_counter": 42.0e-9,      # Raggio degli anelli di Contro-spin (42 nm)
    "wire_thickness": 2.5e-9,       # Spessore fascio di DNA + rame (2.5 nm)
    "copper_density": 8960,         # Densità del rame (kg/m^3)

    # Parametri operativi ottici
    "wavelength_laser": 532e-9,     # Luce verde incidente (532 nm)
    "laser_power_input": 1.0e-6,    # Potenza laser incidente sul nanosistema (1 microWatt)

    # Ambiente viscoso (es. saliva / soluzione acquosa)
    "fluid_viscosity": 1.2e-3,      # Viscosità dinamica media (Pa*s)
    "fluid_density": 998.2,         # Densità del fluido (kg/m^3), acqua a 20 C
    "temperature": 310.0,           # Temperatura corporea (K)
}

# Calcolo del volume e massa degli anelli di rame
vol_ring_1 = 2 * (np.pi ** 2) * (SCALE_PARAMS["r_ring_primary"]) * (SCALE_PARAMS["wire_thickness"] / 2) ** 2
vol_ring_2 = 2 * (np.pi ** 2) * (SCALE_PARAMS["r_ring_counter"]) * (SCALE_PARAMS["wire_thickness"] / 2) ** 2
mass_ring_1 = 3 * vol_ring_1 * SCALE_PARAMS["copper_density"]  # 3 anelli
mass_ring_2 = 3 * vol_ring_2 * SCALE_PARAMS["copper_density"]  # 3 anelli

# Momenti d'inerzia (I = m * r^2 per anelli sottili rispetto all'asse polare)
I_spin = mass_ring_1 * (SCALE_PARAMS["r_ring_primary"] ** 2)
I_contro = mass_ring_2 * (SCALE_PARAMS["r_ring_counter"] ** 2)

print(f"[*] Raggio Rotore Spin (DNA+Cu):       {SCALE_PARAMS['r_ring_primary'] * 1e9:.2f} nm")
print(f"[*] Raggio Rotore Contro-spin:         {SCALE_PARAMS['r_ring_counter'] * 1e9:.2f} nm")
print(f"[*] Massa Rotore Primario:             {mass_ring_1:.3e} kg")
print(f"[*] Massa Rotore Contro-spin:          {mass_ring_2:.3e} kg")
print(f"[*] Rapporto Momenti Inerzia (I1 / I2): {I_spin / I_contro:.4f}\n")

# ==============================================================================
# PUNTO 4: FORMALIZZAZIONE MATEMATICA DELLO SPIN 720° ED EQUAZIONI DI BILANCIO
# ==============================================================================
print("=" * 60)
print("2. FORMALIZZAZIONE MATEMATICA (SU(2) E BILANCIO DI COPPIA)")
print("=" * 60)

# Definizione simbolica con SymPy
theta = sp.Symbol('theta', real=True)
hbar = sp.Symbol('hbar', positive=True)

# Matrici di Pauli
sigma_x = sp.Matrix([[0, 1], [1, 0]])
sigma_y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
sigma_z = sp.Matrix([[1, 0], [0, -1]])
identity_2 = sp.Matrix([[1, 0], [0, 1]])

# Operatore di rotazione per spinori SU(2) attorno all'asse Z: U(theta) = exp(-i * theta/2 * sigma_z)
U_rot = sp.Matrix([
    [sp.exp(-sp.I * theta / 2), 0],
    [0, sp.exp(sp.I * theta / 2)]
])

# Valutazione a 360° (2*pi) e 720° (4*pi)
U_360 = U_rot.subs(theta, 2 * sp.pi)
U_720 = U_rot.subs(theta, 4 * sp.pi)

print("Operatore SU(2) a 360° (2*pi):")
sp.pprint(U_360)
print("-> Nota: U(360°) = -I (Inversione di segno dello spinore)\n")

print("Operatore SU(2) a 720° (4*pi):")
sp.pprint(U_720)
print("-> Nota: U(720°) = +I (Ritorno all'identità: ciclo spinoriale completo confermato matematicamente)\n")

# Equazione di conservazione della coppia (Torque Cancellation)
# tau_spin + tau_contro + tau_luce = 0
# omega_2 = - (I_spin / I_contro) * omega_1 per evitare auto-rotazione del nanorobot
k_rot = float(I_spin / I_contro)
print(f"[*] Condizione di stabilità angolare: omega_contro(t) = - {k_rot:.4f} * omega_spin(t)")

# ==============================================================================
# PUNTO 2 & 4: CALCOLO DELLA SPINTA FOTONICA (Poynting & Tensore di Maxwell)
# ==============================================================================
print("\n" + "=" * 60)
print("3. CALCOLO DELLA SPINTA FOTONICA DA FLUSSI ASIMMETRICI")
print("=" * 60)

# La forza fotonica F_net è data dal tasso di trasferimento di impulso:
# F_z = (P / c) * (1 + R_contro * cos(alpha)) * eta_capture
c = const.c  # Velocità della luce nel vuoto
P_in = SCALE_PARAMS["laser_power_input"]
R_contro = 0.85      # Efficienza riflettente/respingente del Contro-spin plasmonico
eta_spin = 0.70      # Efficienza di cattura/focalizzazione del sistema Spin
alpha_collim = np.radians(15.0)  # Angolo medio di deflessione assiale dei fotoni

# Spinta propulsiva calcolata
F_photonic = (P_in / c) * (eta_spin + R_contro * np.cos(alpha_collim))

print(f"[*] Potenza ottica incidente:          {P_in * 1e6:.1f} uW")
print(f"[*] Forza di propulsione fotonica F_z: {F_photonic * 1e15:.3f} fN ({F_photonic * 1e12:.6f} pN)")
accel = F_photonic / (mass_ring_1 + mass_ring_2)
print(f"[*] Accelerazione teorica nel vuoto:   {accel:.2e} m/s^2")

# ==============================================================================
# PUNTO 1 & 2: SIMULAZIONE DINAMICA DEI FLUSSI FOTONICI (Traiettorie Runge-Kutta)
# ==============================================================================
print("\n" + "=" * 60)
print("4. SIMULAZIONE DEI FLUSSI DINAMICI DEI FOTONI (Tracciamento)")
print("=" * 60)


def photon_field(t, state):
    """
    Equazioni differenziali del moto per i pacchetti fotonici interagenti
    con il campo combinato di Spin (attrazione/guida) e Contro-spin (repulsione/spinta).
    state = [x, y, z]
    """
    x, y, z = state
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2) + 1e-12

    # Campo di rotazione vorticoso (Spin) attorno all'asse Z
    v_vortex_x = -y / (r + 0.1)
    v_vortex_y = x / (r + 0.1)

    # Campo repulsivo assiale (Contro-spin indirizzato verso +Z)
    repulsion_z = 1.8 / (1.0 + np.exp(-z * 5.0))
    radial_confinement = -0.4 * np.array([x, y, 0])

    dx_dt = v_vortex_x + radial_confinement[0]
    dy_dt = v_vortex_y + radial_confinement[1]
    dz_dt = repulsion_z
    return [dx_dt, dy_dt, dz_dt]


# Generazione di 16 linee di flusso (fotoni incidenti)
trajectories = []
t_span = (0, 3.0)
t_eval = np.linspace(0, 3.0, 150)
for angle in np.linspace(0, 2 * np.pi, 16, endpoint=False):
    init_pos = [0.8 * np.cos(angle), 0.8 * np.sin(angle), -1.2]
    sol = solve_ivp(photon_field, t_span, init_pos, t_eval=t_eval, method='RK45')
    trajectories.append(sol.y)

# Grafico scientifico dei risultati
fig = plt.figure(figsize=(9, 7))
ax = fig.add_subplot(111, projection='3d')
for traj in trajectories:
    ax.plot(traj[0], traj[1], traj[2], color='#00d4ff', alpha=0.7, lw=1.2)

# Disegno indicativo del nanosistema al centro
u = np.linspace(0, 2 * np.pi, 50)
ax.plot(0.4 * np.cos(u), 0.4 * np.sin(u), 0, color='#e67e22', lw=2.5, label='Rotore Spin (Rame)')
ax.plot(0.55 * np.cos(u), 0.55 * np.sin(u), 0.1, color='#f1c40f', lw=2.5, label='Rotore Contro-spin')
ax.scatter([0], [0], [0], color='#ffffff', s=80, edgecolors='black', label='Core Centrale')

# Vettore spinta risultante
ax.quiver(0, 0, 0.2, 0, 0, 1.2, color='#2ecc71', lw=3.0, arrow_length_ratio=0.2, label='Vettore Spinta F_z')
ax.set_title("Motore Fotonico v4 - Studio dei Flussi e Spinta Risultante", fontsize=12)
ax.set_xlabel("X (nm normalizzato)")
ax.set_ylabel("Y (nm normalizzato)")
ax.set_zlabel("Z - Asse di Propulsione")
ax.legend(loc='upper right', fontsize=8)
plt.tight_layout()
plt.savefig(out_path("simulazione_motore_fotonico_v3.png"), dpi=200)
plt.close(fig)
print("[*] Grafico di simulazione dinamica salvato: 'simulazione_motore_fotonico_v3.png'")

# ==============================================================================
# SEZIONE INTEGRATIVA: ATTRITO VISCOSO (STOKES) VS VUOTO SPAZIALE
# ==============================================================================
print("\n" + "=" * 60)
print("5. ANALISI COMPARATIVA: VUOTO VS MEZZO VISCOSO (FLUIDO BIOLOGICO)")
print("=" * 60)

eta = SCALE_PARAMS["fluid_viscosity"]
T = SCALE_PARAMS["temperature"]


def trasporto(r_eff, m_tot, F=F_photonic, rho=SCALE_PARAMS["fluid_density"]):
    """Grandezze di trasporto in fluido per una sfera equivalente di raggio r_eff."""
    kT = const.k * T
    drag = 6 * np.pi * eta * r_eff
    v = F / drag
    D = kT / drag
    D_rot = kT / (8 * np.pi * eta * r_eff ** 3)
    return dict(drag=drag, v=v, D=D,
                Re=rho * F / (3 * np.pi * eta ** 2), Pe=F * r_eff / kT,
                v_th=np.sqrt(kT / m_tot),
                D_rot=D_rot, tau_rot=1 / (2 * D_rot), l_pers=v / (2 * D_rot),
                t_cross=2 * D / v ** 2, L_crit=D / v,
                Deff_over_D=1 + v ** 2 / (6 * D_rot * D))


# Raggio efficace della geometria V4: <r^2>^(1/2) = sqrt(5.7567) u.a. = 239.93 nm.
V4_RADIUS = np.sqrt(5.7567) * 100e-9
V4_MASS = 1.2463630573655231e-18 + 4.487710753596799e-18

gen1 = trasporto(SCALE_PARAMS["r_ring_counter"], mass_ring_1 + mass_ring_2)
v4 = trasporto(V4_RADIUS, V4_MASS)

for nome, d in (("prima generazione", gen1), ("V4", v4)):
    print(nome, {k: f"{x:.4g}" for k, x in d.items()})

# Geometria di riferimento del modello: la V4 della scena (cfr. Tabella 9 del preprint).
r_eff = V4_RADIUS
m_total = V4_MASS
drag_coeff = v4["drag"]
v_terminal = v4["v"]
D_brown = v4["D"]
v_thermal = v4["v_th"]
Re = v4["Re"]
Pe = v4["Pe"]

print(f"[*] Viscosità dinamica del mezzo:        {eta:.2e} Pa*s")
print(f"[*] Raggio idrodinamico efficace:        {r_eff * 1e9:.1f} nm")
print(f"[*] Coefficiente di drag 6*pi*eta*r:     {drag_coeff:.3e} N*s/m")
print(f"[*] Velocità terminale di deriva (Stokes): {v_terminal * 1e6:.3f} um/s")
print(f"[*] Velocità termica browniana sqrt(kT/m): {v_thermal:.3f} m/s")
print(f"[*] Coefficiente di diffusione D:        {D_brown:.3e} m^2/s")
print(f"[*] Numero di Reynolds:                  {Re:.3e} (regime di Stokes valido)")
print(f"[*] Numero di Péclet (drift/diffusione): {Pe:.4f}")
print(f"[*] Spostamento browniano in 1 s:        {np.sqrt(2 * D_brown * 1.0) * 1e6:.3f} um")
print(f"[*] Spostamento direzionale in 1 s:      {v_terminal * 1e6:.3f} um")
print(f"[*] Diffusione rotazionale D_r:          {v4['D_rot']:.3f} 1/s")
print(f"[*] Tempo di persistenza tau_r:          {v4['tau_rot']:.4f} s")
print(f"[*] Lunghezza di persistenza l = v*tau_r: {v4['l_pers'] * 1e9:.1f} nm")
print(f"[*] Incrocio deriva/diffusione t*:       {v4['t_cross']:.3f} s")
print(f"[*] Lunghezza critica L* = D/v:          {v4['L_crit'] * 1e6:.3f} um")
print(f"[*] Diffusione efficace D_eff/D:         {v4['Deff_over_D']:.4f}")

# Confronto dei regimi: distanza percorsa nel tempo (deriva vs diffusione)
t_axis = np.linspace(0.0, 10.0, 400)
fig2, ax2 = plt.subplots(figsize=(9, 5))
ax2.plot(t_axis, v_terminal * t_axis * 1e6, color='#2ecc71', lw=2.0,
         label='Deriva fotonica a direzione costante')
ax2.plot(t_axis, np.sqrt(2 * D_brown * t_axis) * 1e6, color='#3498db', lw=2.0,
         label=r'Diffusione browniana $\sqrt{2Dt}$')
ax2.plot(t_axis, np.sqrt(2 * D_brown * v4["Deff_over_D"] * t_axis) * 1e6,
         color='#8e44ad', lw=1.6, ls='--',
         label=r'Diffusione efficace $\sqrt{2D_{\mathrm{eff}}t}$')
ax2.axvline(v4["t_cross"], color='0.5', lw=0.9, ls=':')
ax2.annotate('t* = %.2f s' % v4["t_cross"], (v4["t_cross"], 0.0),
             xytext=(v4["t_cross"] * 1.05, 0.25), fontsize=9, color='0.35')
ax2.set_title("Motore Fotonico V4 \u2014 Deriva direzionale e diffusione termica "
              "(geometria V4, r = %.0f nm)" % (r_eff * 1e9), fontsize=11)
ax2.set_xlabel("Tempo (s)")
ax2.set_ylabel(r"Spostamento ($\mu$m)")
ax2.grid(alpha=0.3)
ax2.legend(loc='upper left', fontsize=9)
plt.tight_layout()
plt.savefig(out_path("confronto_regimi_stokes_browniano.png"), dpi=200)
plt.savefig(os.path.join(OUT_DIR, "versione5", "Figura_13_regimi_stokes_browniano.png"),
            dpi=200)
plt.close(fig2)
print("[*] Grafico dei regimi salvato: 'confronto_regimi_stokes_browniano.png' "
      "e 'versione5/Figura_13_regimi_stokes_browniano.png'")

print("\n[!] LIMITI DEL MODELLO")
print("    - Nel vuoto la propulsione e' libera (F_z / m), ma senza mezzo di reazione")
print("      il bilancio complessivo va valutato sul sistema chiuso (fotoni inclusi).")
print("    - In fluido biologico, per T = 310 K, Pe = 0.28: a scala r domina la diffusione,")
print("      e la diffusione rotazionale (tau_r ~ 0.05 s) rimescola la direzione della spinta:")
print("      il moto risultante e' una diffusione efficace con deriva trascurabile.")
print("    - Le efficienze eta_spin e R_contro sono parametri di input, non derivati")
print("      da un calcolo di scattering completo: i valori di forza vanno letti come")
print("      stime d'ordine di grandezza, non come previsioni quantitative.")
print("\n[+] Calcoli teorici e fisici completati con successo!\n")
