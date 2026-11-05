"""
MOTORE FOTONICO DEL VUOTO - V7.2
v72_tabelle.py - tabelle LaTeX generate da dati/V72_sintesi.json.

Nessun numero e' scritto a mano: se cambia il modello, cambiano le tabelle.

Uso:
    python v72_tabelle.py
"""

import json
import math
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
if BASE not in sys.path:
    sys.path.insert(0, BASE)
import v72_contrappeso as v72  # noqa: E402
DATI = os.path.normpath(os.path.join(BASE, "..", "dati"))
TAB = os.path.normpath(os.path.join(BASE, "..", "tabelle"))
os.makedirs(TAB, exist_ok=True)


def f(x, nd=2):
    if x is None:
        return "--"
    try:
        if not math.isfinite(float(x)):
            return "--"
    except Exception:
        return "--"
    return ("%." + str(nd) + "f") % float(x)


def e(x, nd=1):
    if x is None:
        return "--"
    try:
        if not math.isfinite(float(x)):
            return "--"
    except Exception:
        return "--"
    return ("%." + str(nd) + "e") % float(x)


def intervalli(ws):
    if not ws:
        return "--"
    return "; ".join("%.2f--%.2f" % (w["R_min_um"], w["R_max_um"]) for w in ws)


def scrivi(nome, testo):
    p = os.path.join(TAB, nome)
    with open(p, "w", encoding="utf-8") as fh:
        fh.write(testo)
    print("  scritto %-30s %5.1f kB" % (nome, os.path.getsize(p) / 1024.0))


def main():
    print("V7.2 - tabelle")
    d = json.load(open(os.path.join(DATI, "V72_sintesi.json"),
                       encoding="utf-8"))

    # 1. costo energetico
    righe = []
    for r in d["costo_energetico"]:
        righe.append(" & ".join([
            f(r["R_um"], 1), f(r["T_360_no_s"]), f(r["T_720_no_s"]),
            f(r["T_720_half_s"]), f(r["T_720_full_s"]),
            f(r["E_360_no_mJ"]), f(r["E_720_half_mJ"]),
            f(r["E_720_full_mJ"]), f(r["costo_720_half"]),
            f(r["costo_720_full"]),
        ]) + r" \\")
    scrivi("tab_V72_costo.tex", r"""\begin{tabular}{rrrrrrrrrr}
\hline
$R$ & $T_{360}^{\rm no}$ & $T_{720}^{\rm no}$ & $T_{720}^{1/2}$ &
$T_{720}^{\rm full}$ & $E_{360}^{\rm no}$ & $E_{720}^{1/2}$ &
$E_{720}^{\rm full}$ & $\frac{E_{720}^{1/2}}{E_{720}^{\rm no}}$ &
$\frac{E_{720}^{\rm full}}{E_{720}^{\rm no}}$ \\
\micro m & s & s & s & s & mJ & mJ & mJ & -- & -- \\
\hline
""" + "\n".join(righe) + r"""
\hline
\end{tabular}
""")

    # 2. confronto lunghezze d'onda
    righe = []
    for r in d["confronto_lambda"]:
        righe.append(" & ".join([
            f(r["R_um"], 1), "%d" % r["lam_nm"], f(r["P_target_uW"], 0),
            f(r["P_inc_mW"], 2), f(r["delta"], 2), f(r["dT_K"], 1),
            f(r["T_720_s"], 2), f(r["E_720_mJ"], 1),
            "sì" if r["valido_sicuro"] else ("no" if not r["valido"] else "solo op."),
        ]) + r" \\")
    scrivi("tab_V72_lambda.tex", r"""\begin{tabular}{rrrrrrrrl}
\hline
$R$ & $\lambda$ & $P_{\rm target}$ & $P_{\rm inc}$ & $\delta$ &
$\Delta T$ & $T_{720}$ & $E_{720}$ & sicuro \\
\micro m & nm & \micro W & mW & rad & K & s & mJ & \\
\hline
""" + "\n".join(righe) + r"""
\hline
\end{tabular}
""")

    # 3. finestre termiche
    righe = []
    for lam in ("532", "800", "1064"):
        for z in ("2_um", "5_um", "10_um", "15_um", "20_um"):
            w = d["finestre_fini"][lam][z]
            v = w["valide"]
            s = w["sicure"]
            dmax = max([x["max_dT_K"] for x in v]) if v else None
            righe.append(" & ".join([
                lam, z.replace("_um", ""), intervalli(v), intervalli(s),
                f(dmax, 1),
            ]) + r" \\")
    scrivi("tab_V72_finestre.tex", r"""\begin{tabular}{rrp{4.2cm}p{4.2cm}r}
\hline
$\lambda$ & $z$ & valido ($\Delta T<20$ K) & sicuro ($\Delta T<5$ K) &
$\max \Delta T$ \\
nm & \micro m & \micro m & \micro m & K \\
\hline
""" + "\n".join(righe) + r"""
\hline
\end{tabular}
""")

    # 4. rumore
    righe = []
    for r in d["rumore"]:
        righe.append(" & ".join([
            f(r["R_um"], 1), f(r["D_r_rad2_s"], 3), f(r["tau_r_s"], 2),
            f(r["sigma_theta_360_deg"], 0), f(r["N_cicli_snr10_dt20"], 0),
            f(r["N_cicli_snr10_dt100"], 0), f(r["T_720_s"], 2),
        ]) + r" \\")
    scrivi("tab_V72_rumore.tex", r"""\begin{tabular}{rrrrrrr}
\hline
$R$ & $D_r$ & $\tau_r$ & $\sigma_\theta(360^\circ)$ &
$N_{10}$ (20 bin) & $N_{10}$ (100 bin) & $T_{720}$ \\
\micro m & rad$^2$/s & s & deg & cicli & cicli & s \\
\hline
""" + "\n".join(righe) + r"""
\hline
\end{tabular}
""")

    # 5. viscosita'
    righe = []
    for r in d["viscosita"]:
        righe.append(" & ".join([
            f(r["R_um"], 1), f(r["mu_mPa_s"], 0), f(r["T_720_s"], 1),
            f(r["E_720_mJ"], 1), f(r["dT_K"], 1), e(r["Re_rot"], 1),
            f(r["tau_su_bingham"], 1),
            "sì" if r["valido_sicuro"] else ("no" if not r["valido"] else "solo op."),
        ]) + r" \\")
    scrivi("tab_V72_viscosita.tex", r"""\begin{tabular}{rrrrrrrl}
\hline
$R$ & $\mu$ & $T_{720}$ & $E_{720}$ & $\Delta T$ & $Re$ &
$\tau/\tau_B$ & sicuro \\
\micro m & mPa\,s & s & mJ & K & -- & -- & \\
\hline
""" + "\n".join(righe) + r"""
\hline
\end{tabular}
""")

    # 6. dispersione di taglia: gli zeri della forma si smorzano
    righe = []
    for lam, R in ((800, 2.34e-6), (1064, 3.11e-6), (532, 1.55e-6)):
        for disp in (0.10, 0.15, 0.20, 0.30):
            c = v72.ciclo(R, "vaterite", lam, 2e-6, dispersione=disp,
                          latch_fraction=0.5)
            righe.append(" & ".join([
                "%d" % lam, f(R * 1e6, 2), "%.0f" % (disp * 100),
                f(c["P_target"] * 1e6, 0), f(c["dT"], 1),
            ]) + r" \\")
    scrivi("tab_V72_dispersione.tex", r"""\begin{tabular}{rrrrr}
\hline
$\lambda$ & $R$ & dispersione & $P_{\rm target}$ & $\Delta T$ \\
nm & \micro m & \% & \micro W & K \\
\hline
""" + "\n".join(righe) + r"""
\hline
\end{tabular}
""")

    print("Fatto: tabelle in", TAB)


if __name__ == "__main__":
    main()
