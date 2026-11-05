"""
MOTORE FOTONICO DEL VUOTO - V7.2
v72_rumore.py - moto browniano rotazionale e tempo di integrazione.

La V7.1 disegnava la modulazione di omega come una curva deterministica.  Qui
si integra l'equazione di Langevin overdamped

    d(theta) = omega(theta) dt + sqrt(2 D_r dt) xi(t)

con omega(theta) modulato una volta per giro (latch a meta' ciclo) e uno stato
A/B che e' dispari sotto 360 gradi (periodo 720).  Il codice produce:

  - la traccia rumorosa;
  - la media ripiegata su theta mod 2 pi (omega) e mod 4 pi (stato A/B);
  - lo spettro dello stato, con il subarmonico a f_rot/2;
  - il numero di cicli necessari per SNR = 10.

Uso:
    python v72_rumore.py
"""

import json
import math
import os
import sys

import numpy as np
from scipy import signal
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.abspath(__file__))
V7 = os.path.normpath(os.path.join(BASE, "..", "..", "versione 7", "codice"))
for p in (V7, BASE):
    if p not in sys.path:
        sys.path.insert(0, p)

import v72_contrappeso as v72  # noqa: E402

DATI = os.path.normpath(os.path.join(BASE, "..", "dati"))
FIG = os.path.normpath(os.path.join(BASE, "..", "figure"))
os.makedirs(DATI, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 130, "savefig.dpi": 190, "font.size": 9.5,
    "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.5,
    "axes.edgecolor": "#444444", "legend.frameon": False,
    "figure.facecolor": "white",
})
C = {532: "#c1121f", 800: "#0077b6", 1064: "#588157"}


def omega_di(theta, omega_no, omega_cw):
    """Latch a meta' ciclo: omega_cw per 0..pi, omega_no per pi..2pi."""
    fase = np.mod(theta, 2.0 * math.pi)
    return np.where(fase < math.pi, omega_cw, omega_no)


def stato_di(theta):
    """Stato A/B del collegamento: periodo 4 pi (720 gradi)."""
    return np.where(np.mod(theta, 4.0 * math.pi) < 2.0 * math.pi, 1.0, -1.0)


def simula(R=0.3e-6, lam=800, z=2e-6, Pe=10.0, n_cicli=200,
           passi_per_ciclo=2000, seed=7, mu=None):
    c = v72.ciclo(R, "vaterite", lam, z, Pe_target=Pe, mu=mu,
                  latch_fraction=0.5)
    tau = c["tau"]
    omega_no = tau / c["beta_no"]
    omega_cw = tau / c["beta_cw"]
    T360_no = 2.0 * math.pi / omega_no
    T360_cw = 2.0 * math.pi / omega_cw
    T360 = 0.5 * (T360_no + T360_cw)       # media del latch a meta' ciclo
    D_r = c["D_r_eff"]

    dt = T360 / passi_per_ciclo
    n = int(n_cicli * passi_per_ciclo)
    rng = np.random.default_rng(seed)
    theta = np.empty(n + 1)
    theta[0] = 0.0
    xi = rng.standard_normal(n)
    for i in range(n):
        th = theta[i]
        theta[i + 1] = th + omega_di(th, omega_no, omega_cw) * dt + \
            math.sqrt(2.0 * D_r * dt) * xi[i]

    omega_meas = np.diff(theta) / dt
    stato = stato_di(theta[:-1])
    omega_det = omega_di(theta[:-1], omega_no, omega_cw)

    # media ripiegata: omega su theta mod 2 pi, stato su theta mod 4 pi
    def ripiega(valori, angolo, periodo, nbin):
        fase = np.mod(angolo, periodo) / periodo
        idx = np.minimum((fase * nbin).astype(int), nbin - 1)
        somma = np.bincount(idx, weights=valori, minlength=nbin)
        conta = np.bincount(idx, minlength=nbin)
        return somma / np.maximum(conta, 1), conta

    nbin_w = 40
    nbin_s = 80
    med_w, _ = ripiega(omega_meas, theta[:-1], 2.0 * math.pi, nbin_w)
    med_s, _ = ripiega(stato, theta[:-1], 4.0 * math.pi, nbin_s)
    fase_w = (np.arange(nbin_w) + 0.5) / nbin_w * 2.0 * math.pi
    fase_s = (np.arange(nbin_s) + 0.5) / nbin_s * 4.0 * math.pi
    det_w = omega_di(fase_w, omega_no, omega_cw)
    det_s = stato_di(fase_s)
    rms_w = math.sqrt(float(np.mean((med_w - det_w) ** 2)))
    contrasto_s = 0.5 * float(np.mean(med_s[fase_s < 2.0 * math.pi]) -
                              np.mean(med_s[fase_s >= 2.0 * math.pi]))
    rms_s = math.sqrt(float(np.mean((med_s - det_s) ** 2)))

    # medie ripiegate per N crescente: il protocollo giusto e' ripiegare su
    # theta, non sperare che l'FFT cieco separi una riga allargata da D_r.
    folded = {}
    for nc in (10, 30, 100):
        m = min(int(nc * passi_per_ciclo), len(stato))
        mw, _ = ripiega(omega_meas[:m], theta[:m], 2.0 * math.pi, nbin_w)
        ms, _ = ripiega(stato[:m], theta[:m], 4.0 * math.pi, nbin_s)
        c_s = 0.5 * float(np.mean(ms[fase_s < 2.0 * math.pi]) -
                          np.mean(ms[fase_s >= 2.0 * math.pi]))
        folded["omega_%d" % nc] = mw.tolist()
        folded["stato_%d" % nc] = ms.tolist()
        folded["contrasto_%d" % nc] = c_s

    # spettro dello stato e di omega per N cicli crescenti
    fs = 1.0 / dt
    spettri = {}
    for nc in (10, 30, 100):
        m = min(int(nc * passi_per_ciclo), len(stato))
        f, P = signal.welch(stato[:m], fs=fs, nperseg=min(2 ** 14, m),
                            noverlap=None, detrend="constant")
        spettri["stato_%d" % nc] = {"f": f.tolist(), "P": P.tolist()}
        f2, P2 = signal.welch(omega_meas[:m] - omega_meas[:m].mean(),
                              fs=fs, nperseg=min(2 ** 14, m),
                              noverlap=None, detrend="constant")
        spettri["omega_%d" % nc] = {"f": f2.tolist(), "P": P2.tolist()}

    f_rot = 1.0 / T360
    f_sub = 0.5 / T360

    def prominenza(f, P, f0, banda=0.15):
        mask = (f > f0 - banda) & (f < f0 + banda)
        if not np.any(mask):
            return 0.0
        i = np.argmin(np.abs(f - f0))
        fondo = np.median(P[mask])
        return float(P[i] / max(fondo, 1e-30))

    prom = {k: prominenza(np.array(v["f"]), np.array(v["P"]), f_sub)
            for k, v in spettri.items() if k.startswith("stato_")}

    return {
        "R_um": R * 1e6, "lam_nm": lam, "z_um": z * 1e6,
        "Pe_r": Pe, "n_cicli": n_cicli, "passi_per_ciclo": passi_per_ciclo,
        "seed": seed,
        "mu_Pa_s": (v72.SANGUE["mu"] if mu is None else mu),
        "beta_tot": c["beta_tot"],
        "omega_no": omega_no, "omega_cw": omega_cw,
        "T360_no_s": T360_no, "T360_cw_s": T360_cw, "T360_s": T360,
        "T720_s": 2.0 * T360, "D_r_eff": D_r,
        "tau_r_eff_s": 1.0 / (2.0 * D_r),
        "f_rot_Hz": f_rot, "f_sub_Hz": f_sub,
        "rms_omega_folded": rms_w,
        "contrasto_stato": contrasto_s,
        "rms_stato_folded": rms_s,
        "prominenza_sub_10": prom.get("stato_10", 0.0),
        "prominenza_sub_30": prom.get("stato_30", 0.0),
        "prominenza_sub_100": prom.get("stato_100", 0.0),
        "theta": theta.tolist(),
        "omega_meas": omega_meas.tolist(),
        "stato": stato.tolist(),
        "fase_w": fase_w.tolist(), "med_w": med_w.tolist(),
        "det_w": det_w.tolist(),
        "fase_s": fase_s.tolist(), "med_s": med_s.tolist(),
        "det_s": det_s.tolist(),
        "folded": folded,
        "spettri": spettri,
    }


def figura(d, nome="fig_V72_rumore.png"):
    th = np.array(d["theta"])
    om = np.array(d["omega_meas"])
    giri = th / (2.0 * math.pi)
    n_show = int(4.0 * d["passi_per_ciclo"])
    fig, axes = plt.subplots(2, 2, figsize=(12.6, 8.4))
    fig.subplots_adjust(hspace=0.42, wspace=0.24)

    ax = axes[0, 0]
    ax.plot(giri[:n_show], om[:n_show] / d["omega_no"], color="#0077b6",
            lw=0.6, alpha=0.65, label="$\\omega$ misurato")
    ax.plot(giri[:n_show], omega_di(th[:n_show], d["omega_no"], d["omega_cw"])
            / d["omega_no"], color="#c1121f", lw=1.0, ls="--",
            label="deterministico")
    ax.set_xlabel("giri")
    ax.set_ylabel("$\\omega/\\omega_{\\rm no}$")
    ax.set_title("(a) traccia rumorosa (4 giri)")
    ax.legend(fontsize=7.5)

    ax = axes[0, 1]
    for nc, col, ls in ((10, "#c1121f", ":"), (30, "#e07a1f", "--"),
                        (100, "#0077b6", "-")):
        ax.plot(np.array(d["fase_w"]) / (2 * math.pi),
                np.array(d["folded"]["omega_%d" % nc]) / d["omega_no"],
                color=col, lw=1.2, ls=ls, label="%d cicli" % nc)
    ax.plot(np.array(d["fase_w"]) / (2 * math.pi),
            np.array(d["det_w"]) / d["omega_no"],
            color="#333333", lw=0.9, ls="--", label="deterministico")
    ax.set_xlabel("fase del giro (360°)")
    ax.set_ylabel("$\\omega/\\omega_{\\rm no}$")
    ax.set_title("(b) media ripiegata: $\\omega$ ha periodo 360°")
    ax.legend(fontsize=7.5)

    ax = axes[1, 0]
    for nc, col, ls in ((10, "#c1121f", ":"), (30, "#e07a1f", "--"),
                        (100, "#0077b6", "-")):
        ax.plot(np.array(d["fase_s"]) / (4 * math.pi),
                np.array(d["folded"]["stato_%d" % nc]),
                color=col, lw=1.4, ls=ls, label="%d cicli" % nc)
    ax.plot(np.array(d["fase_s"]) / (4 * math.pi),
            np.array(d["det_s"]), color="#333333", lw=0.9, ls="--",
            label="deterministico")
    ax.set_xlabel("fase del ciclo di 720° (0 = A, 0.5 = B, 1 = A)")
    ax.set_ylabel("stato A/B")
    ax.set_ylim(-1.25, 1.25)
    ax.set_title("(c) lo stato 720° ripiegato su $\\theta$: 10 cicli bastano")
    ax.legend(fontsize=7.5, loc="lower right")

    ax = axes[1, 1]
    for k, col, ls in (("stato_10", "#c1121f", ":"),
                       ("stato_30", "#e07a1f", "--"),
                       ("stato_100", "#0077b6", "-")):
        f = np.array(d["spettri"][k]["f"])
        P = np.array(d["spettri"][k]["P"])
        m = f < 3.0
        ax.semilogy(f[m], P[m] + 1e-30, color=col, lw=1.2, ls=ls,
                    label="%s cicli" % k.split("_")[1])
    ax.axvline(d["f_rot_Hz"], color="#555555", lw=0.8, ls=":")
    ax.axvline(d["f_sub_Hz"], color="#555555", lw=0.8, ls=":")
    y0 = ax.get_ylim()[0]
    ax.text(d["f_rot_Hz"] * 1.02, y0 * 3, "$f_{\\rm rot}$",
            fontsize=7.5, color="#444444")
    ax.text(d["f_sub_Hz"] * 0.55, y0 * 3, "$f_{\\rm rot}/2$",
            fontsize=7.5, color="#444444")
    ax.set_xlim(0, 3.0)
    ax.set_xlabel("frequenza (Hz)")
    ax.set_ylabel("PSD dello stato A/B")
    ax.set_title("(d) FFT cieca: riga larga, serve integrazione lunga")
    ax.legend(fontsize=7.5)

    fig.suptitle("V7.2 - moto browniano: con il segnale di orientazione lo "
                 "stato 720° si ripiega in 10 cicli;\nsenza orientazione "
                 "l'FFT cieca è allargata da $D_r$", y=1.03)
    p = os.path.join(FIG, nome)
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    print("  %-30s %6.1f kB" % (nome, os.path.getsize(p) / 1024.0))


def main():
    print("V7.2 - rumore rotazionale")
    d = simula()
    # Il JSON tiene la sintesi, le medie ripiegate e gli spettri; la traccia
    # completa (400k punti) si rigenera dal seed, non serve salvarla.
    slim = {k: v for k, v in d.items()
            if k not in ("theta", "omega_meas", "stato")}
    with open(os.path.join(DATI, "V72_rumore.json"), "w", encoding="utf-8") as fh:
        json.dump(slim, fh, indent=1, ensure_ascii=False, default=float)
    print("  scritto V72_rumore.json")
    figura(d)
    print("  D_r = %.3f rad^2/s, f_rot = %.3f Hz, f_sub = %.3f Hz"
          % (d["D_r_eff"], d["f_rot_Hz"], d["f_sub_Hz"]))
    print("  contrasto stato = %.3f, prominenza sub (10/30/100) = "
          "%.2f / %.2f / %.2f"
          % (d["contrasto_stato"], d["prominenza_sub_10"],
             d["prominenza_sub_30"], d["prominenza_sub_100"]))
    print("Fatto.")


if __name__ == "__main__":
    main()
