"""
MOTORE FOTONICO V5 - COPPIA UNIFICATA (momento angolare ceduto)
Autore: Dr. Francesco Bulla

Questo modulo sostituisce le tre definizioni di coppia ottica che convivevano
nella V5 con una sola:

    tau = P_ass / omega_ott        (momento angolare ceduto dal fascio)
        = F * lambda / (2 pi)      (equivalente, con P_ass = F c)

dove F e' la forza di scattering calcolata da v5_nucleo.scattering_globale(k)
e omega_ott = 2 pi c / lambda.  La versione "a fascio intero" della stessa
coppia e' P / omega_ott: e' il limite superiore che si ottiene quando il
rotore intercetta tutta la potenza.

Il modulo produce:
  * dati/U1_coppia_unificata.json   - la tabella completa alle tre scale
  * le figure di §4.2, §4.5, §5.1, §5.3 e §5.5 ricalcolate su questo modello

Uso:
    python v5_coppia_unificata.py
"""

import json
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
if BASE not in sys.path:
    sys.path.insert(0, BASE)

import v5_nucleo as n  # noqa: E402

DATI = os.path.join(BASE, "dati")
FIG = os.path.join(BASE, "figure")
os.makedirs(DATI, exist_ok=True)

OMEGA_OTT = 2.0 * math.pi * n.C_LIGHT / n.LAMBDA
LAMBDA = n.LAMBDA

BLU = "#1f4e79"
ARANCIO = "#c55a11"
ROSSO = "#9c2f22"
VERDE = "#3c6e47"
GRIGIO = "#7a7a7a"


# ---------------------------------------------------------------------------
# 1. LA COPPIA UNIFICATA
# ---------------------------------------------------------------------------
def coppia_oam(F):
    """Coppia da momento angolare ceduto: tau = F * lambda / (2 pi)."""
    return F * LAMBDA / (2.0 * math.pi)


def caso(k, nome):
    """Tutte le grandezze di una scala sotto il modello unificato."""
    sc = n.scattering_globale(k)
    F = float(sc["F_z_N"])
    F_decl = float(sc["F_z_dichiarata_N"])
    f_int = float(sc["frazione_fascio"])
    r_eff = n.R_EFF_V4 * k
    m = n.M_V4 * k ** 3
    beta = n.beta_telaio(n.R_DISC * k)
    D_r = n.KT / (8.0 * math.pi * n.ETA * r_eff ** 3)

    tau = coppia_oam(F)
    tau_pieno = n.P_LASER / OMEGA_OTT          # limite a fascio intero
    omega = tau / beta
    omega_pieno = tau_pieno / beta

    return {
        "nome": nome,
        "k": k,
        "r_eff_nm": r_eff * 1e9,
        "massa_kg": m,
        "forza_scattering_N": F,
        "frazione_intercettata": f_int,
        "potenza_assorbita_W": n.P_LASER * f_int,
        "beta_telaio_Nms": beta,
        "D_rot_1_s": D_r,
        "tau_oam_Nm": tau,
        "tau_fascio_intero_Nm": tau_pieno,
        "tau_dichiarata_Nm": F_decl * n.R_DISC * k,
        "tau_su_kT": tau / n.KT,
        "omega_ss_rad_s": omega,
        "omega_ss_fascio_intero_rad_s": omega_pieno,
        "Pe_tau_su_kT": tau / n.KT,
        "Pe_omega_su_Dr": omega / D_r,
        "verdetto": "termico" if tau / n.KT < 1.0 else "deterministico",
    }


# ---------------------------------------------------------------------------
# 2. LE TRE SCALE
# ---------------------------------------------------------------------------
def tabella():
    casi = [caso(1e-3, "1 nm"), caso(1.0, "240 nm"), caso(1e6, "1 m")]
    return casi


# ---------------------------------------------------------------------------
# 3. VUOTO PROFONDO (v5.1): coppia ridotta dall'intercettazione
# ---------------------------------------------------------------------------
def vuoto(R, materiale="Rame"):
    """t_crit con la coppia a fascio intero e con quella ridotta."""
    k = R / n.R_DISC
    sc = n.scattering_globale(k)
    f_int = float(sc["frazione_fascio"])
    sigma = 200e6 if materiale == "Rame" else 50e9
    rho = 8960.0 if materiale == "Rame" else 2200.0
    kk = k ** 3 * rho / 8960.0
    I_tot = (n.I1_V4 + n.I2_V4) * k ** 5 * rho / 8960.0
    m = n.M_V4 * kk
    omega_lim = math.sqrt(sigma / rho) / R
    tau_pieno = n.P_LASER / OMEGA_OTT
    out = {}
    for etichetta, tau in (("fascio_intero", tau_pieno),
                           ("intercettata", tau_pieno * f_int)):
        out[etichetta] = {
            "tau_Nm": tau,
            "alpha_rad_s2": tau / I_tot,
            "t_crit_s": omega_lim * I_tot / tau,
        }
    out.update({
        "R_m": R,
        "materiale": materiale,
        "massa_kg": m,
        "I_tot_kgm2": I_tot,
        "omega_rottura_rad_s": omega_lim,
        "frazione_intercettata": f_int,
    })
    return out


# ---------------------------------------------------------------------------
# 4. FIGURE
# ---------------------------------------------------------------------------
def _pulisci(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _salva(fig, nome):
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, nome), bbox_inches="tight")
    plt.close(fig)
    print("  %-42s %6.1f kB"
          % (nome, os.path.getsize(os.path.join(FIG, nome)) / 1024.0))


def fig_p2(casi):
    """Coppia a due lobi e transitorio, con la coppia unificata."""
    c = casi[1]
    theta = np.linspace(0, 360, 361)
    due_lobi = 0.5 * (1 + np.cos(2 * np.pi * theta / 180.0))

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 3.9))
    a1.plot(theta, due_lobi, color=BLU, lw=2.2)
    a1.fill_between(theta, 0, due_lobi, color=BLU, alpha=0.10)
    a1.axhline(0.0, color="#888888", lw=0.8)
    a1.set_xlabel("angolo del rotore $\\theta$ (gradi)")
    a1.set_ylabel("$\\tau(\\theta)/\\tau_{max}$")
    a1.set_title("(a) Due lobi, periodo $\\pi$")
    a1.set_xlim(0, 360)
    a1.set_ylim(-0.05, 1.08)
    a1.annotate("massimi a 0 e $\\pi$\n$\\tau=0$ a $\\pi/2$ e $3\\pi/2$",
                xy=(180, 1.0), xytext=(196, 0.62), fontsize=8.4, color=ROSSO,
                arrowprops=dict(arrowstyle="->", color=ROSSO, lw=1.0))

    beta = c["beta_telaio_Nms"]
    I1 = n.I1_V4
    t_in = I1 / beta
    tt = np.logspace(math.log10(t_in) - 2, math.log10(t_in) + 2, 60)
    a2.semilogx(tt, 1 - np.exp(-tt / t_in), color=ARANCIO, lw=2.2)
    a2.axhline(1.0, color=GRIGIO, ls=":", lw=0.9)
    a2.set_xlabel("tempo (s)")
    a2.set_ylabel("$\\omega/\\omega_{ss}$")
    a2.set_title("(b) Transitorio, $\\tau_{in}$ = %.2e s" % t_in)
    a2.set_ylim(0, 1.12)
    txt = ("$\\tau$ = %.2e N m ($P_{ass}/\\omega_{ott}$)\n"
           "$\\beta$ telaio = %.2e N m s\n"
           "$\\omega_{ss}$ = %.3e rad/s\n"
           "$Pe_r=\\tau/kT$ = %.2e"
           % (c["tau_oam_Nm"], beta, c["omega_ss_rad_s"], c["Pe_tau_su_kT"]))
    a2.text(0.28, 0.06, txt, transform=a2.transAxes, fontsize=8.2, color=GRIGIO,
            bbox=dict(boxstyle="round,pad=0.4", fc="#f7f7f7", ec="#dddddd"))
    fig.suptitle("P2 - Termine dissipativo e coppia da momento angolare",
                 fontsize=11.5, color=BLU, y=1.02)
    _salva(fig, "fig_P2_dissipazione.png")


def fig_p5(casi):
    """Coppia contro potenza e periodicita' reale del modello."""
    c = casi[1]
    P = np.logspace(-9, -1, 200)
    tau = (P / OMEGA_OTT) * c["frazione_intercettata"]
    rumore = n.KT

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 3.9))
    a1.loglog(P, tau, color=BLU, lw=2.2, label="$\\tau=P_{ass}/\\omega_{ott}$")
    a1.axhline(rumore, color=ROSSO, ls="--", lw=1.4, label="rumore $kT$")
    P_inc = rumore * OMEGA_OTT / c["frazione_intercettata"]
    a1.axvline(P_inc, color=GRIGIO, ls=":", lw=1.0)
    a1.plot([P_inc], [rumore], marker="o", ms=6, color=ROSSO)
    a1.annotate("$P$ = %.2e W" % P_inc, xy=(P_inc, rumore),
                xytext=(P_inc * 3, rumore * 30), fontsize=8.4, color=ROSSO,
                arrowprops=dict(arrowstyle="->", color=ROSSO, lw=1.0))
    a1.set_xlabel("potenza del laser (W)")
    a1.set_ylabel("coppia (N m)")
    a1.set_title("(a) Coppia contro rumore termico")
    a1.legend(loc="lower right", fontsize=8.4)

    th = np.linspace(0, 720, 721)
    seg = 0.5 * (1 + np.cos(2 * np.pi * th / 180.0))
    a2.plot(th, seg, color=ARANCIO, lw=2.0)
    a2.axvline(360, color=GRIGIO, ls=":", lw=0.9)
    a2.set_xlabel("angolo del rotore (gradi)")
    a2.set_ylabel("$\\tau(\\theta)/\\tau_{max}$")
    a2.set_title("(b) Periodo $\\pi$: nessuna firma a 720$^{\\circ}$")
    a2.set_xlim(0, 720)
    a2.set_ylim(-0.05, 1.08)
    a2.annotate("lo stesso valore a 0 e a 360$^{\\circ}$",
                xy=(360, 1.0), xytext=(384, 0.55), fontsize=8.4, color=ROSSO,
                arrowprops=dict(arrowstyle="->", color=ROSSO, lw=1.0))
    fig.suptitle("P5 - Soglia di rumore e periodicita' della coppia",
                 fontsize=11.5, color=BLU, y=1.02)
    _salva(fig, "fig_P5_protocollo.png")


def fig_sweep():
    """Peclet e velocita' di regime su nove ordini di grandezza."""
    k = np.logspace(-3, 6, 200)
    w, pe_tau, pe_st = [], [], []
    for kk in k:
        sc = n.scattering_globale(kk)
        F = float(sc["F_z_N"])
        beta = n.beta_telaio(n.R_DISC * kk)
        r = n.R_EFF_V4 * kk
        Dr = n.KT / (8 * math.pi * n.ETA * r ** 3)
        tau = coppia_oam(F)
        w.append(tau / beta)
        pe_tau.append(tau / n.KT)
        pe_st.append((tau / beta) / Dr)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.8, 3.9))
    a1.loglog(k, pe_tau, color=BLU, lw=2.2, label="$Pe_r=\\tau/kT$")
    a1.loglog(k, pe_st, color=ARANCIO, lw=1.8, ls="--",
              label="$Pe_r=\\omega_{ss}/D_r$ (Stokes)")
    a1.axhline(1.0, color=ROSSO, ls=":", lw=1.2)
    a1.set_xlabel("fattore di scala $k$ (1 = 240 nm)")
    a1.set_ylabel("$Pe_r$")
    a1.set_title("(a) P\\'eclet rotazionale: sempre sotto 1")
    a1.legend(loc="upper left", fontsize=8.2)

    a2.loglog(k, w, color=VERDE, lw=2.2)
    a2.set_xlabel("fattore di scala $k$")
    a2.set_ylabel("$\\omega_{ss}$ (rad/s)")
    a2.set_title("(b) $\\omega_{ss}=\\tau/\\beta$: nessuna invarianza")
    a2.annotate("$k^{-1}$", xy=(3e3, w[150]), fontsize=10, color=GRIGIO)
    fig.suptitle("S1-S3 - Legge di scala con la coppia da momento angolare",
                 fontsize=11.5, color=BLU, y=1.02)
    _salva(fig, "fig_S1_sweep_scala.png")


def fig_sp1(casi):
    """Determinismo alle tre configurazioni."""
    nomi = [c["nome"] for c in casi]
    om = [c["omega_ss_rad_s"] for c in casi]
    pe = [c["Pe_tau_su_kT"] for c in casi]
    x = np.arange(3)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.8, 3.9))
    a1.bar(x - 0.2, om, 0.4, color=BLU, label="$\\omega_{ss}=\\tau/\\beta$")
    a1.bar(x + 0.2, [c["omega_ss_fascio_intero_rad_s"] for c in casi], 0.4,
           color=ARANCIO, alpha=0.75, label="fascio intero")
    a1.set_yscale("log")
    a1.set_xticks(x)
    a1.set_xticklabels(nomi)
    a1.set_ylabel("rad/s")
    a1.set_title("(a) Velocit\\`a di regime")
    a1.legend(fontsize=8.2)

    a2.bar(x, pe, 0.55, color=VERDE)
    a2.axhline(1.0, color=ROSSO, ls=":", lw=1.2)
    a2.set_yscale("log")
    a2.set_xticks(x)
    a2.set_xticklabels(nomi)
    a2.set_ylabel("$Pe_r=\\tau/kT$")
    a2.set_title("(b) Nessuna scala supera la soglia")
    for i, v in enumerate(pe):
        a2.annotate("%.1e" % v, xy=(i, v), xytext=(i - 0.28, v * 8),
                    fontsize=8.0, color=GRIGIO)
    fig.suptitle("SP1 - Determinismo con la coppia unificata",
                 fontsize=11.5, color=BLU, y=1.02)
    _salva(fig, "fig_SP1_scala_determinismo.png")


def fig_sp3():
    """Materiali a 1 m: limite di rottura contro regime."""
    materiali = ["Grafene/CNT", "SiC", "Si/SiO$_2$"]
    rho = [1350.0, 3210.0, 2250.0]
    sig = [5e9, 3e9, 0.5e9]
    R = 0.5
    om_lim = [math.sqrt(s / (r * R ** 2)) for s, r in zip(sig, rho)]
    c = caso(1e6, "1 m")
    regime = c["tau_oam_Nm"] / n.beta_telaio(n.R_DISC * 1e6)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.6, 3.9))
    masse = [n.M_V4 * (R / n.R_DISC) ** 3 * r / 8960.0 for r in rho]
    a1.bar(materiali, masse, color=[VERDE, BLU, ARANCIO])
    a1.set_ylabel("massa totale (kg)")
    a1.set_title("(a) Massa dei due rotori")
    a2.bar(materiali, om_lim, color=[VERDE, BLU, ARANCIO])
    a2.axhline(regime, color=ROSSO, ls="--", lw=1.2)
    a2.set_yscale("log")
    a2.set_ylabel("$\\omega$ limite (rad/s)")
    a2.set_title("(b) Limite di rottura contro $\\omega_{ss}$")
    a2.annotate("regime %.1e" % regime, xy=(1.5, regime), xytext=(0.6, regime * 40),
                fontsize=8.2, color=ROSSO)
    fig.suptitle("SP3 - Tre materiali a 1 m", fontsize=11.5, color=BLU, y=1.02)
    _salva(fig, "fig_SP3_materiali_1m.png")


# ---------------------------------------------------------------------------
# 5. MAIN
# ---------------------------------------------------------------------------
def main():
    print("[V5] coppia unificata: tau = P_ass / omega_ott = F lambda / 2 pi")
    casi = tabella()

    vuoto_nano = vuoto(120e-9, "Rame")
    vuoto_nano_g = vuoto(120e-9, "Grafene/CNT")
    vuoto_macro = vuoto(0.5, "Rame")

    out = {
        "modello": "tau = P_ass/omega_ott = F lambda/(2 pi)",
        "omega_ottica_rad_s": OMEGA_OTT,
        "lambda_m": LAMBDA,
        "kT_J": n.KT,
        "casi": casi,
        "soglia_potenza_Pe1_W": n.KT * OMEGA_OTT,
        "soglia_potenza_Pe1_intercettata_W":
            n.KT * OMEGA_OTT / casi[1]["frazione_intercettata"],
        "vuoto_profondo": {
            "nano_rame": vuoto_nano,
            "nano_grafene": vuoto_nano_g,
            "macro_rame": vuoto_macro,
        },
    }
    percorso = os.path.join(DATI, "U1_coppia_unificata.json")
    with open(percorso, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, ensure_ascii=False)
    print("  %-42s %6.1f kB" % ("U1_coppia_unificata.json",
                                os.path.getsize(percorso) / 1024.0))

    print("\n  %-8s %12s %12s %12s %12s %10s"
          % ("scala", "tau [Nm]", "omega [rad/s]", "Pe(tau/kT)", "Pe(Stokes)", "esito"))
    for c in casi:
        print("  %-8s %12.3e %12.3e %12.3e %12.3e %10s"
              % (c["nome"], c["tau_oam_Nm"], c["omega_ss_rad_s"],
                 c["Pe_tau_su_kT"], c["Pe_omega_su_Dr"], c["verdetto"]))
    print("\n  soglia Pe_r = 1: P = %.3e W (fascio intero), %.3e W (intercettata)"
          % (out["soglia_potenza_Pe1_W"],
             out["soglia_potenza_Pe1_intercettata_W"]))
    for nome, d in (("nano rame", vuoto_nano), ("nano grafene", vuoto_nano_g),
                    ("macro rame", vuoto_macro)):
        print("  %-13s f=%.4e  t_crit(pieno)=%.3e s  t_crit(intercettata)=%.3e s"
              % (nome, d["frazione_intercettata"],
                 d["fascio_intero"]["t_crit_s"], d["intercettata"]["t_crit_s"]))

    fig_p2(casi)
    fig_p5(casi)
    fig_sweep()
    fig_sp1(casi)
    fig_sp3()


if __name__ == "__main__":
    main()
