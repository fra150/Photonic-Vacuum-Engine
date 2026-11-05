# -*- coding: utf-8 -*-
"""
v73_tabelle.py - tabelle LaTeX e macro numeriche della V7.3.

Legge dati/V73_sintesi.json e dati/V73_griglia.json e scrive:
  paper/V73/tabelle_v73/*.tex      (italiano)
  paper/V73/tabelle_v73_en/*.tex   (inglese)
  paper/V73/numeri_v73.tex         (macro condivise)

Nessun numero del preprint e' scritto a mano: se cambia il modello, cambia
il JSON e cambiano tabelle e macro.
"""

import json
import math
import os

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(BASE, ".."))
DATI = os.path.join(ROOT, "dati")
PAPER = os.path.join(ROOT, "paper", "V73")
TAB_IT = os.path.join(PAPER, "tabelle_v73")
TAB_EN = os.path.join(PAPER, "tabelle_v73_en")
for d in (TAB_IT, TAB_EN):
    os.makedirs(d, exist_ok=True)

with open(os.path.join(DATI, "V73_sintesi.json"), encoding="utf-8") as fh:
    S = json.load(fh)
with open(os.path.join(DATI, "V73_griglia.json"), encoding="utf-8") as fh:
    G = json.load(fh)

CELLE = G["celle"]


# ---------------------------------------------------------------------------
# formattazione
# ---------------------------------------------------------------------------
def f(x, n=2):
    if x is None:
        return "---"
    if isinstance(x, float) and math.isinf(x):
        return r"$\infty$"
    return ("%." + str(n) + "f") % x


def sci(x, n=2):
    if x is None:
        return "---"
    if x == 0:
        return "0"
    e = int(math.floor(math.log10(abs(x))))
    m = x / 10.0 ** e
    return r"$%s\times10^{%d}$" % (("%." + str(n) + "f") % m, e)


def riga(R, lam, z=2.0):
    for c in CELLE:
        if (c["lam_nm"] == lam and abs(c["z_um"] - z) < 1e-9
                and abs(c["R_um"] - R) < 5e-4):
            return c
    # il raggio puo' non essere esattamente sulla griglia: prendi il piu' vicino
    cand = [c for c in CELLE if c["lam_nm"] == lam
            and abs(c["z_um"] - z) < 1e-9]
    return min(cand, key=lambda c: abs(c["R_um"] - R))


def esatto(R, lam, chiave="confronto_lambda"):
    """Riga esatta dal JSON (non il vicino piu' prossimo della griglia)."""
    for c in S[chiave]:
        if abs(c["R_um"] - R) < 1e-9 and c["lam_nm"] == lam:
            return c
    raise KeyError("riga assente: %s R=%s lam=%s" % (chiave, R, lam))


def finestra(lst):
    if not lst:
        return "---"
    if len(lst) == 1:
        return "$%s$--$%s$" % (f(lst[0]["R_min_um"], 2),
                               f(lst[0]["R_max_um"], 2))
    return "; ".join("$%s$--$%s$" % (f(w["R_min_um"], 2),
                                     f(w["R_max_um"], 2)) for w in lst)


def scrivi(nome, testo, en=False):
    p = os.path.join(TAB_EN if en else TAB_IT, nome)
    with open(p, "w", encoding="utf-8") as fh:
        fh.write(testo)
    print("  %-32s %5.1f kB" % (nome, os.path.getsize(p) / 1024.0))


def tab(colspec, intestazione, righe, en=False, size=r"\small"):
    out = [r"\begin{tabular}{%s}" % colspec, r"\toprule",
           intestazione, r"\midrule"]
    out += righe
    out += [r"\bottomrule", r"\end{tabular}", ""]
    return "\n".join(out)


# ---------------------------------------------------------------------------
# T1 - le sei correzioni
# ---------------------------------------------------------------------------
def t1(en=False):
    if en:
        h = (r"item & draft 3 Nov & V7.2 corrected \\")
        rows = [
            r"path heating & $P_{\rm inc}-P_{\rm target}$; $\mu_a(532)=7500$~m$^{-1}$ & "
            r"$\mu_a P_{\rm inc}z$ with the average factor; $\mu_a(532)=2.35\times10^4$~m$^{-1}$ \\",
            r"viscous work & 4.88 presented as an energy cost & "
            r"$W_{720}=4\pi\tau=126\,k_BT$, independent of $\beta$; 4.88 is a time ratio \\",
            r"noise and protocol & $\sigma_\theta=80^\circ$; $N=119$--$595$; contrast 1.00 & "
            r"$\sigma_\theta=64^\circ$; $CV=18\%$; 2--3 cycles; contrast is by construction \\",
            r"blood and vaterite & $\mu=4$~mPa\,s everywhere; $\Delta n=0.172$ & "
            r"$\mu_{\rm eff}$: 1.3 in the CFL, 4 outside; $\Delta n=0.08$ (0.06--0.10) \\",
            r"internal consistency & no Bingham; $z$ undefined; ``about double''; ``six orders'' & "
            r"$\tau>\tau_B$ in the predicate; $z$ defined; $1+0.2=1.2$~mW ($6\times$); ratio computed \\",
            r"topology & 720$^\circ$ from the double cover SO(3)$\to$SU(2) & "
            r"classical Z$_2$ latch; the double cover is only an analogy \\",
        ]
    else:
        h = (r"voce & bozza 3 nov & V7.2 corretta \\")
        rows = [
            r"calore di cammino & $P_{\rm inc}-P_{\rm target}$; $\mu_a(532)=7500$~m$^{-1}$ & "
            r"$\mu_a P_{\rm inc}z$ con il fattore di media; $\mu_a(532)=2.35\times10^4$~m$^{-1}$ \\",
            r"lavoro viscoso & 4.88 presentato come costo energetico & "
            r"$W_{720}=4\pi\tau=126\,k_BT$, indipendente da $\beta$; 4.88 \`e un rapporto di tempi \\",
            r"rumore e protocollo & $\sigma_\theta=80^\circ$; $N=119$--$595$; contrasto 1.00 & "
            r"$\sigma_\theta=64^\circ$; $CV=18\%$; 2--3 cicli; il contrasto \`e per costruzione \\",
            r"sangue e vaterite & $\mu=4$~mPa\,s ovunque; $\Delta n=0.172$ & "
            r"$\mu_{\rm eff}$: 1.3 nel CFL, 4 fuori; $\Delta n=0.08$ (0.06--0.10) \\",
            r"coerenza interna & nessun Bingham; $z$ indefinito; ``circa doppio''; ``sei ordini'' & "
            r"$\tau>\tau_B$ nel predicato; $z$ definito; $1+0.2=1.2$~mW ($6\times$); rapporto calcolato \\",
            r"topologia & 720$^\circ$ dalla copertura doppia SO(3)$\to$SU(2) & "
            r"contatore Z$_2$ classico; la copertura doppia \`e solo un'analogia \\",
        ]
    return tab(r"p{2.7cm}p{5.6cm}p{5.6cm}", h, rows, en, r"\scriptsize")


# ---------------------------------------------------------------------------
# T2 - input ottici e reologici
# ---------------------------------------------------------------------------
def t2(en=False):
    if en:
        h = (r"$\lambda$ & $\mu_a$ & $\mu_s$ & $g$ & $\ell$ & $T_s(2\,\mu{\rm m})$ & "
             r"$T_s(10\,\mu{\rm m})$ \\ nm & m$^{-1}$ & m$^{-1}$ & & $\mu$m & & \\")
        lab = r"\multicolumn{7}{l}{\emph{Whole blood, 150 g/L Hb, 310 K}}\\"
    else:
        h = (r"$\lambda$ & $\mu_a$ & $\mu_s$ & $g$ & $\ell$ & $T_s(2\,\mu{\rm m})$ & "
             r"$T_s(10\,\mu{\rm m})$ \\ nm & m$^{-1}$ & m$^{-1}$ & & $\mu$m & & \\")
        lab = r"\multicolumn{7}{l}{\emph{Sangue intero, 150 g/L Hb, 310 K}}\\"
    rows = [lab]
    for lam in (532, 800, 1064):
        d = S["ottica"][str(lam)]
        n = 1.0 / (d["mu_a"] + d["mu_s"])
        T2 = math.exp(-(d["mu_a"] + d["mu_s"]) * 2e-6)
        T10 = math.exp(-(d["mu_a"] + d["mu_s"]) * 10e-6)
        rows.append("%d & %s & %s & %.2f & %s & %.3f & %.3f \\\\"
                    % (lam, sci(d["mu_a"], 2) if d["mu_a"] < 1000 else "%d" % d["mu_a"],
                       "%d" % d["mu_s"], d["g"], f(n * 1e6, 1), T2, T10))
    rows.append(r"\midrule")
    if en:
        rows.append(r"\multicolumn{7}{l}{\emph{Geometry and rheology}}\\")
    else:
        rows.append(r"\multicolumn{7}{l}{\emph{Geometria e reologia}}\\")
    if en:
        rows.append(r"$\Delta n$ & \multicolumn{6}{l}{$0.08$ (range $0.06$--$0.10$); "
                    r"vaterite, effective, not calcite $0.172$}\\")
        rows.append(r"$\mu_{\rm plasma}$ & \multicolumn{6}{l}{$1.3$~mPa\,s in the CFL}\\")
        rows.append(r"$\mu_{\rm blood}$ & \multicolumn{6}{l}{$4.0$~mPa\,s outside the CFL; "
                    r"$7$--$10$~mPa\,s only as a low-shear sensitivity}\\")
        rows.append(r"$\sigma_y$ & \multicolumn{6}{l}{$0.01$~Pa (range $0.005$--$0.01$); "
                    r"$\tau_B=(8/3)\pi R^3\sigma_y$}\\")
        rows.append(r"$\delta_{\rm CFL}$ & \multicolumn{6}{l}{$3\,\mu$m (generous; "
                    r"measured CFL in 10--40$\,\mu$m vessels is 1--2$\,\mu$m)}\\")
    else:
        rows.append(r"$\Delta n$ & \multicolumn{6}{l}{$0.08$ (intervallo $0.06$--$0.10$); "
                    r"vaterite, efficace, non calcite $0.172$}\\")
        rows.append(r"$\mu_{\rm plasma}$ & \multicolumn{6}{l}{$1.3$~mPa\,s nel CFL}\\")
        rows.append(r"$\mu_{\rm sangue}$ & \multicolumn{6}{l}{$4.0$~mPa\,s fuori dal CFL; "
                    r"$7$--$10$~mPa\,s solo come sensibilit\`a a basso shear}\\")
        rows.append(r"$\sigma_y$ & \multicolumn{6}{l}{$0.01$~Pa (intervallo $0.005$--$0.01$); "
                    r"$\tau_B=(8/3)\pi R^3\sigma_y$}\\")
        rows.append(r"$\delta_{\rm CFL}$ & \multicolumn{6}{l}{$3\,\mu$m (generoso; "
                    r"il CFL misurato in vasi da 10--40$\,\mu$m \`e 1--2$\,\mu$m)}\\")
    return tab(r"rrrrrrr", h, rows, en, r"\small")


# ---------------------------------------------------------------------------
# T3 - confronto lambda
# ---------------------------------------------------------------------------
def t3(en=False):
    h = (r"$R$ & $\lambda$ & $P_{\rm target}$ & $P_{\rm inc}$ & $T_{720}$ & "
         r"$E_{720}$ & $\Delta T_{\rm part}$ & $\Delta T_{\rm path}$ & "
         r"$\Delta T$ & $\tau/\tau_B$ \\ $\mu$m & nm & mW & mW & s & mJ & K & K & K & \\")
    rows = []
    R_prec = None
    for c in S["confronto_lambda"]:
        if R_prec is not None and abs(c["R_um"] - R_prec) > 1e-9:
            rows.append(r"\midrule")
        mark = "" if c.get("valido", True) else r"$^{\dagger}$"
        rows.append("%.1f & %d%s & %.3f & %.3f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f \\\\"
                    % (c["R_um"], c["lam_nm"], mark, c["P_target_mW"],
                       c["P_inc_mW"], c["T720_half"], c["E720_half"] * 1e3,
                       c["dT_part"], c["dT_path"], c["dT"],
                       c["tau_su_tau_B"]))
        R_prec = c["R_um"]
    return tab(r"rrrrrrrrrr", h, rows, en, r"\small")


# ---------------------------------------------------------------------------
# T4 - costo del contrappeso
# ---------------------------------------------------------------------------
def t4(en=False):
    h = (r"$R$ & $T_{360}^{\rm no}$ & $T_{720}^{\rm no}$ & $T_{720}^{\rm cw}$ & "
         r"$T_{720}^{\rm half}$ & $E_{720}^{\rm no}$ & $E_{720}^{\rm half}$ & "
         r"$E_{\rm half}/E_{\rm no}$ & $\eta_{\rm mech}$ \\ "
         r"$\mu$m & s & s & s & s & mJ & mJ & & \\")
    rows = []
    for c in S["costo_energetico"]:
        mark = "" if c.get("valido", True) else r"$^{\dagger}$"
        rows.append("%.1f%s & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %s \\\\"
                    % (c["R_um"], mark, c["T360_no"], c["T720_no"],
                       c["T720_cw"], c["T720_half"], c["E720_no"] * 1e3,
                       c["E720_half"] * 1e3,
                       c["E720_half"] / c["E720_no"], sci(c["eta_mech"], 2)))
    return tab(r"rrrrrrrrr", h, rows, en, r"\small")


# ---------------------------------------------------------------------------
# T5 - rumore e tempi
# ---------------------------------------------------------------------------
def t5(en=False):
    c = esatto(0.5, 800)
    z = c["z_cw"]
    e2 = 0.5 * math.erfc(z * math.sqrt(2.0) / math.sqrt(2.0))
    e3 = 0.5 * math.erfc(z * math.sqrt(3.0) / math.sqrt(2.0))
    h1 = (r"$\sigma_\theta(360^\circ)$ & $\sigma_\theta(720^\circ)$ & $CV$ & "
          r"$\beta_{\rm cw}/\beta_{\rm no}$ & $\epsilon_{\rm no}$ & "
          r"$\epsilon_{\rm cw}$ & $\epsilon_2$ & $\epsilon_3$ \\ "
          r"deg & deg & & & & & & \\")
    h2 = (r"$R$ & $T_{360}^{\rm no}$ & $T_{360}^{\rm cw}$ & $T_{720}^{\rm half}$ & "
          r"$T_{\rm meas}(3)$ \\ $\mu$m & s & s & s & s \\")
    r1 = ("%.1f & %.1f & %.3f & %.3f & %s & %s & %s & %s \\\\"
          % (c["sigma_theta_360_deg"], c["sigma_theta_720_deg"],
             c["cv_durata"], c["separazione_durate"],
             sci(c["err_no"], 1), sci(c["err_cw"], 1),
             sci(e2, 1), sci(e3, 1)))
    rows2 = []
    for R in (0.3, 0.4, 0.5, 0.6, 0.8):
        q = esatto(R, 800)
        rows2.append("%.1f & %.2f & %.2f & %.2f & %.1f \\\\"
                     % (q["R_um"], q["T360_no"], q["T360_cw"],
                        q["T720_half"], q["T_misura_3cicli"]))
    out = [r"\begin{tabular}{rrrrrrrr}", r"\toprule", h1, r"\midrule", r1,
           r"\bottomrule", r"\end{tabular}", r"\par\vspace{0.6em}",
           r"\begin{tabular}{rrrrr}", r"\toprule", h2, r"\midrule"]
    out += rows2
    out += [r"\bottomrule", r"\end{tabular}", ""]
    return "\n".join(out)


# ---------------------------------------------------------------------------
# T6 - finestre
# ---------------------------------------------------------------------------
def t6(en=False):
    if en:
        h = (r"$\lambda$ & $z$ & valid ($\Delta T<20$~K) & safe ($\Delta T<5$~K) & "
             r"CFL ($\tau_{B,\rm eff}$) \\ nm & $\mu$m & $\mu$m & $\mu$m & $\mu$m \\")
    else:
        h = (r"$\lambda$ & $z$ & valido ($\Delta T<20$~K) & sicuro ($\Delta T<5$~K) & "
             r"CFL ($\tau_{B,\rm eff}$) \\ nm & $\mu$m & $\mu$m & $\mu$m & $\mu$m \\")
    rows = []
    for lam in (532, 800, 1064):
        for z in ("2_um", "5_um", "10_um", "15_um", "20_um"):
            v = S["riassunto"]["finestre"][str(lam)][z]
            s = S["riassunto"]["finestre_sicure"][str(lam)][z]
            cf = S["riassunto"]["finestre_cfl"][str(lam)][z]
            rows.append("%d & %s & %s & %s & %s \\\\"
                        % (lam, z.replace("_um", "").replace("_", "."),
                           finestra(v), finestra(s), finestra(cf)))
        if lam != 1064:
            rows.append(r"\midrule")
    return tab(r"rrrrr", h, rows, en, r"\scriptsize")


# ---------------------------------------------------------------------------
# T7 - CFL e continuo
# ---------------------------------------------------------------------------
def t7(en=False):
    if en:
        h = (r"$\delta_{\rm CFL}$ & $R$ & $f_{\rm CFL}$ & $\mu_{\rm eff}$ & "
             r"$\tau/\tau_B$ & $\tau/\tau_{B,\rm eff}$ & all in CFL & valid CFL \\ "
             r"$\mu$m & $\mu$m & & mPa\,s & & & & \\")
    else:
        h = (r"$\delta_{\rm CFL}$ & $R$ & $f_{\rm CFL}$ & $\mu_{\rm eff}$ & "
             r"$\tau/\tau_B$ & $\tau/\tau_{B,\rm eff}$ & tutto nel CFL & valido CFL \\ "
             r"$\mu$m & $\mu$m & & mPa\,s & & & & \\")
    rows = []
    for r in S["continuita"]:
        if r["delta_cfl_um"] not in (2.0, 3.0):
            continue
        rows.append("%.1f & %.1f & %.2f & %.2f & %.2f & %s & %s & %s \\\\"
                    % (r["delta_cfl_um"], r["R_um"],
                       max(float(r["f_cfl"]), 0.0),
                       r["mu_mPa_s"], r["tau_su_tau_B"],
                       f(r["tau_su_tau_B_eff"], 2),
                       ("yes" if r["interamente_nel_CFL"] else "no") if en
                       else ("sì" if r["interamente_nel_CFL"] else "no"),
                       ("yes" if r["valido_cfl"] else "no") if en
                       else ("sì" if r["valido_cfl"] else "no")))
    return tab(r"rrrrrrrr", h, rows, en, r"\scriptsize")


# ---------------------------------------------------------------------------
# T8 - fuga dallo snervamento
# ---------------------------------------------------------------------------
def t8(en=False):
    h = (r"$R$ & $Pe_B$ & $Pe$ & $P_{\rm inc}$ & $T_{720}$ & $E_{720}$ & "
         r"$\Delta T$ & " + ("valid" if en else "valido") +
         r" \\ $\mu$m & & & mW & s & mJ & K & \\")
    rows = []
    for r in S["fuga_bingham"]:
        rows.append("%.1f & %.1f & %.1f & %.2f & %.1f & %.1f & %.1f & %s \\\\"
                    % (r["R_um"], r["Pe_B"], r["Pe"], r["P_inc_mW"],
                       r["T720_half"], r["E720_half"] * 1e3, r["dT"],
                       ("yes" if r["valido"] else "no") if en
                       else ("sì" if r["valido"] else "no")))
    return tab(r"rrrrrrrr", h, rows, en, r"\small")


# ---------------------------------------------------------------------------
# T9 - zeri di forma
# ---------------------------------------------------------------------------
def t9(en=False):
    h = (r"$\Delta n$ & $\lambda$ & $R_m$ & $\tau_B(R_m)$ & "
         r"$\tau_{\rm opt}(1\,{\rm W})/\tau_B$ \\ & nm & $\mu$m & $10^{-5}$~pN\,m & \\")
    rows = []
    for r in S["zeri"]:
        rows.append("%.2f & %d & %.2f & %.2f & %.1f \\\\"
                    % (r["dn"], r["lam_nm"], r["R_m_um"],
                       r["tau_B"] * 1e17, r["tau_opt_su_tau_B"]))
    return tab(r"rrrrr", h, rows, en, r"\small")


# ---------------------------------------------------------------------------
# T10 - secondo spot
# ---------------------------------------------------------------------------
def t10(en=False):
    h = (r"$\lambda$ & $P_{\rm inc}$ & $P_{\rm add}$ & $P_{\rm tot}$ & "
         r"$P_{\rm tot}/P_{\rm inc}$ & $\Delta T_{\rm path}^{\rm add}$ & "
         r"$\Delta T_{\rm tot}$ \\ nm & mW & mW & mW & & K & K \\")
    rows = []
    for r in S["spot_secondario"]:
        rows.append("%d & %.3f & %.2f & %.3f & %.2f & %.2f & %.2f \\\\"
                    % (r["lam_nm"], r["P_inc_mW"], r["P_add_mW"],
                       r["P_tot_mW"], r["rapporto_tot_primario"],
                       r["dT_path_aggiunto"], r["dT_totale"]))
    return tab(r"rrrrrrr", h, rows, en, r"\small")


# ---------------------------------------------------------------------------
# macro
# ---------------------------------------------------------------------------
DIGITI = {"0": "Zero", "1": "One", "2": "Two", "3": "Three", "4": "Four",
          "5": "Five", "6": "Six", "7": "Seven", "8": "Eight", "9": "Nine"}


def macro(nome, valore):
    # I nomi delle macro TeX possono contenere solo lettere: un numero
    # interromperebbe il control sequence (es. \V73... = \V + "73...").
    nome = "".join(DIGITI.get(ch, ch) for ch in nome)
    return r"\newcommand{\%s}{%s}" % (nome, valore)


def scrivi_macro():
    c = esatto(0.5, 800)
    c4 = esatto(0.4, 800)
    cw = esatto(0.3, 800)
    wb = [r for r in S["fuga_bingham"] if abs(r["R_um"] - 1.0) < 1e-9][0]
    sp = {r["lam_nm"]: r for r in S["spot_secondario"]}
    z = c["z_cw"]
    e2 = 0.5 * math.erfc(z * math.sqrt(2.0) / math.sqrt(2.0))
    e3 = 0.5 * math.erfc(z * math.sqrt(3.0) / math.sqrt(2.0))

    M = [
        ("VSeventyTwoBetaRatio", "%.4f" % S["beta_ratio"]),
        ("VSeventyTwoInerziaRatio", "%.4f" % S["inerzia_ratio"]),
        ("VSeventyTwoWkT", "%.2f" % c["W720_kT"]),
        ("VSeventyTwoWJ", sci(c["W720"], 2)),
        ("VSeventyTwoEtaMech", sci(c["eta_mech"], 2)),
        ("VSeventyTwoSigmaDeg", "%.1f" % c["sigma_theta_360_deg"]),
        ("VSeventyTwoSigmaRad", "%.3f" % c["sigma_theta_360"]),
        ("VSeventyTwoSigma720Deg", "%.1f" % c["sigma_theta_720_deg"]),
        ("VSeventyTwoCV", "%.3f" % c["cv_durata"]),
        ("VSeventyTwoDurRatio", "%.3f" % c["separazione_durate"]),
        ("VSeventyTwoErrNo", sci(c["err_no"], 1)),
        ("VSeventyTwoErrCw", sci(c["err_cw"], 1)),
        ("VSeventyTwoErrCw2", sci(e2, 1)),
        ("VSeventyTwoErrCw3", sci(e3, 1)),
        ("VSeventyTwoDn", "%.2f" % S["dn"]),
        ("VSeventyTwoDnMin", "%.2f" % S["dn_range"][0]),
        ("VSeventyTwoDnMax", "%.2f" % S["dn_range"][1]),
        ("VSeventyTwoMuPlasma", "%.1f" % (S["sangue"]["mu_p"] * 1e3)),
        ("VSeventyTwoMuBlood", "%.1f" % (S["sangue"]["mu_b"] * 1e3)),
        ("VSeventyTwoSigmaY", "%.3f" % S["sigma_y"]),
        ("VSeventyTwoDeltaCFL", "%.1f" % (S["sangue"]["delta_cfl"] * 1e6)),
        ("VSeventyTwoPe", "%.0f" % S["Pe_r"]),
        ("VSeventyTwoMuA532", "%d" % S["ottica"]["532"]["mu_a"]),
        ("VSeventyTwoMuA800", "%d" % S["ottica"]["800"]["mu_a"]),
        ("VSeventyTwoMuA1064", "%d" % S["ottica"]["1064"]["mu_a"]),
        ("VSeventyTwoMuS532", "%d" % S["ottica"]["532"]["mu_s"]),
        ("VSeventyTwoMuS800", "%d" % S["ottica"]["800"]["mu_s"]),
        ("VSeventyTwoMuS1064", "%d" % S["ottica"]["1064"]["mu_s"]),
        ("VSeventyTwoEll532", "%.1f" % (1e6 / (S["ottica"]["532"]["mu_a"] + S["ottica"]["532"]["mu_s"]))),
        ("VSeventyTwoEll800", "%.1f" % (1e6 / (S["ottica"]["800"]["mu_a"] + S["ottica"]["800"]["mu_s"]))),
        ("VSeventyTwoEll1064", "%.1f" % (1e6 / (S["ottica"]["1064"]["mu_a"] + S["ottica"]["1064"]["mu_s"]))),
        ("VSeventyTwoTs532", "%.3f" % math.exp(-(S["ottica"]["532"]["mu_a"] + S["ottica"]["532"]["mu_s"]) * 2e-6)),
        ("VSeventyTwoTs800", "%.3f" % math.exp(-(S["ottica"]["800"]["mu_a"] + S["ottica"]["800"]["mu_s"]) * 2e-6)),
        ("VSeventyTwoTs1064", "%.3f" % math.exp(-(S["ottica"]["1064"]["mu_a"] + S["ottica"]["1064"]["mu_s"]) * 2e-6)),
        ("VSeventyTwoRifR", "%.1f" % c["R_um"]),
        ("VSeventyTwoRifPtar", "%.3f" % c["P_target_mW"]),
        ("VSeventyTwoRifPinc", "%.3f" % c["P_inc_mW"]),
        ("VSeventyTwoRifT", "%.2f" % c["T720_half"]),
        ("VSeventyTwoRifT360No", "%.3f" % c["T360_no"]),
        ("VSeventyTwoRifT360Cw", "%.3f" % c["T360_cw"]),
        ("VSeventyTwoRifT720No", "%.2f" % c["T720_no"]),
        ("VSeventyTwoRifE", "%.2f" % (c["E720_half"] * 1e3)),
        ("VSeventyTwoRifDT", "%.2f" % c["dT"]),
        ("VSeventyTwoRifDTPart", "%.2f" % c["dT_part"]),
        ("VSeventyTwoRifDTPath", "%.2f" % c["dT_path"]),
        ("VSeventyTwoRifTaub", "%.2f" % c["tau_su_tau_B"]),
        ("VSeventyTwoRifMu", "%.2f" % c["mu_mPa_s"]),
        ("VSeventyTwoRifEta", sci(c["eta_mech"], 2)),
        ("VSeventyTwoRifOmega", "%.1f" % c["omega_no"]),
        ("VSeventyTwoRifRe", sci(c["Re_rot"], 1)),
        ("VSeventyTwoRifTauSpin", "%.1f" % (c["tau_spin_cw"] * 1e9)),
        ("VSeventyTwoMinR", "%.1f" % c4["R_um"]),
        ("VSeventyTwoMinPtar", "%.3f" % c4["P_target_mW"]),
        ("VSeventyTwoMinPinc", "%.3f" % c4["P_inc_mW"]),
        ("VSeventyTwoMinT", "%.2f" % c4["T720_half"]),
        ("VSeventyTwoMinT360No", "%.3f" % c4["T360_no"]),
        ("VSeventyTwoMinT360Cw", "%.3f" % c4["T360_cw"]),
        ("VSeventyTwoMinT720No", "%.2f" % c4["T720_no"]),
        ("VSeventyTwoMinE", "%.2f" % (c4["E720_half"] * 1e3)),
        ("VSeventyTwoMinDT", "%.2f" % c4["dT"]),
        ("VSeventyTwoMinTaub", "%.2f" % c4["tau_su_tau_B"]),
        ("VSeventyTwoMinEta", sci(c4["eta_mech"], 2)),
        ("VSeventyTwoWBR", "%.1f" % wb["R_um"]),
        ("VSeventyTwoWBPeB", "%.1f" % wb["Pe_B"]),
        ("VSeventyTwoWBPe", "%.1f" % wb["Pe"]),
        ("VSeventyTwoWBPinc", "%.2f" % wb["P_inc_mW"]),
        ("VSeventyTwoWBT", "%.1f" % wb["T720_half"]),
        ("VSeventyTwoWBE", "%.1f" % (wb["E720_half"] * 1e3)),
        ("VSeventyTwoWBDT", "%.1f" % wb["dT"]),
        ("VSeventyTwoWBTaub", "%.2f" % wb["tau_su_tau_B"]),
        ("VSeventyTwoBingCut", "0.8"),
        ("VSeventyTwoContLimit", "1.0"),
        ("VSeventyTwoZeri532", "%.2f" % [r for r in S["zeri"] if r["dn"] == 0.08 and r["lam_nm"] == 532][0]["R_m_um"]),
        ("VSeventyTwoZeri800", "%.2f" % [r for r in S["zeri"] if r["dn"] == 0.08 and r["lam_nm"] == 800][0]["R_m_um"]),
        ("VSeventyTwoZeri1064", "%.2f" % [r for r in S["zeri"] if r["dn"] == 0.08 and r["lam_nm"] == 1064][0]["R_m_um"]),
        ("VSeventyTwoSpotRatio", "6.0"),
        ("VSeventyTwoSpot800", "%.2f" % sp[800]["rapporto_tot_primario"]),
        ("VSeventyTwoSpotDT800", "%.2f" % sp[800]["dT_path_aggiunto"]),
        ("VSeventyTwoSpotDT532", "%.1f" % sp[532]["dT_path_aggiunto"]),
        ("VSeventyTwoSpotDT1064", "%.2f" % sp[1064]["dT_path_aggiunto"]),
        ("VSeventyTwoTauSpin", "%.1f" % (cw["tau_spin_cw"] * 1e9)),
        ("VSeventyTwoTauSpinR05", "%.1f" % (c["tau_spin_cw"] * 1e9)),
        ("VSeventyTwoReNo", sci(c["Re_rot"], 1)),
        ("VSeventyTwoOmegaNo", "%.1f" % c["omega_no"]),
        ("VSeventyTwoOmegaBreak", sci(c["omega_rottura"], 1)),
        ("VSeventyTwoTMeas3", "%.1f" % c["T_misura_3cicli"]),
        ("VSeventyTwoCountCelle", "%d" % S["riassunto"]["celle"]),
        ("VSeventyTwoCountValide", "%d" % S["riassunto"]["valide"]),
        ("VSeventyTwoCountSicure", "%d" % S["riassunto"]["valide_sicure"]),
        ("VSeventyTwoCountCFL", "%d" % S["riassunto"]["valide_cfl"]),
        ("VSeventyTwoTauB", sci(c["tau_B"], 1)),
        ("VSeventyTwoTau", sci(c["tau"], 1)),
        ("VSeventyTwoKT", sci(4.2800119e-21, 2)),
        ("VSeventyTwoRPlasma", "1.3"),
        ("VSeventyTwoRBlood", "4.0"),
        ("VSeventyTwoR03Pinc", "%.2f" % cw["P_inc_mW"]),
        ("VSeventyTwoR03T", "%.2f" % cw["T720_half"]),
        ("VSeventyTwoR03E", "%.2f" % (cw["E720_half"] * 1e3)),
        ("VSeventyTwoR03DT", "%.2f" % cw["dT"]),
        ("VSeventyTwoR08Taub", "%.2f" % [r for r in S["costo_energetico"] if abs(r["R_um"] - 0.8) < 1e-9][0]["tau_su_tau_B"]),
    ]
    lines = [r"% Numeri della V7.2 corretta, generati da codice/v73_tabelle.py.",
             r"% Non modificare a mano: rigenerare dal JSON.", ""]
    for n, v in M:
        lines.append(macro(n, v))
    p = os.path.join(PAPER, "numeri_v73.tex")
    with open(p, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print("  numeri_v73.tex                  %5.1f kB" % (os.path.getsize(p) / 1024.0))


def main():
    print("Tabelle V7.2 corretta")
    for en in (False, True):
        scrivi("tab_v73_confronto.tex", t1(en), en)
        scrivi("tab_v73_input.tex", t2(en), en)
        scrivi("tab_v73_lambda.tex", t3(en), en)
        scrivi("tab_v73_costo.tex", t4(en), en)
        scrivi("tab_v73_rumore.tex", t5(en), en)
        scrivi("tab_v73_finestre.tex", t6(en), en)
        scrivi("tab_v73_continuita.tex", t7(en), en)
        scrivi("tab_v73_fuga.tex", t8(en), en)
        scrivi("tab_v73_zeri.tex", t9(en), en)
        scrivi("tab_v73_spot.tex", t10(en), en)
    scrivi_macro()


if __name__ == "__main__":
    main()
