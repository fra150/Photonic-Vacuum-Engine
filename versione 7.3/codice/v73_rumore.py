"""
MOTORE FOTONICO DEL VUOTO - V7.3
v73_rumore.py - il rumore del latch, calcolato con la prima passa esatta.

La V7.2 trattava la durata di un giro come una gaussiana con CV = 17.8% e
otteneva, per la classificazione A/B, errori dell'ordine di 10^-7 con due o
tre cicli.  Quel calcolo e' sbagliato nelle code: la durata di un giro e' una
prima passa (inverse Gaussian), non una gaussiana.  Le code della prima passa
sono molto piu' pesanti, e -- questa e' la parte sorprendente -- l'errore di
classificazione NON diminuisce aumentando il numero di cicli, perche' anche
la forma (shape) della somma scala con N.

MODELLO ESATTO
--------------
Un giro di 360 gradi sotto coppia costante e rumore browniano ha tempo di
prima passa T ~ IG(m, lambda), con
    m = beta Theta / tau,   lambda = Theta^2 / (2 D_r),   D_r = kT/beta,
cioe' CV^2 = 2 kT/(tau Theta) = 2/(Pe_r Theta).  Per Theta = 2 pi e
Pe_r = 10: CV = 1/sqrt(pi Pe_r) = 0.178.

Somma di N giri nello stesso stato: ancora IG, con media N m e shape N lambda.
Il test A/B corretto e' appaiato: si confronta la somma dei giri veloci
(prima meta' di ogni ciclo) con la somma dei giri lenti (seconda meta').
La soglia geometrica e' sqrt(N m_fast * N m_slow); l'errore totale e'
    eps = P(S_fast > soglia) + P(S_slow < soglia).
Con la IG esatta questo errore e' indipendente da N: la soglia e la shape
scalano insieme.  A CV = 17.8% e r = 3.881 si ottiene eps = 1.03e-4 con la
soglia geometrica e 9.73e-5 con la soglia ottimizzata.  Il budget 10^-4 e'
esattamente saturo: non c'e' margine.

Il test "durata totale del ciclo" (latch contro nessun latch) e' peggiore:
la separazione media e' (1+r)/2 = 2.44, non 3.881, e l'errore a CV = 17.8%
e' 1.12e-2.  Il test appaiato e' obbligatorio.

RITRATTAZIONI
-------------
- "CV fino a 30.2% per 10^-4": approssimazione log-normale a CV uguali.
  La IG esatta da' 17.83%; a CV = 30.2% l'errore e' 2.07e-2, non 10^-4.
- "8.8e-7 con tre cicli": approssimazione gaussiana.  La IG esatta da'
  9.7e-5, e non dipende dal numero di cicli.
- "2.9e-3 per un giro lento": la coda inferiore esatta e' 6.8e-5.
- "due o tre cicli bastano": vero solo per il controllo dei sistematici,
  non per la statistica; la statistica non migliora con N.

CORREZIONE IDRODINAMICA
-----------------------
Il rapporto 3.881 = 1 + r_cw^3 + 0.75 r_cw d^2 (r_cw = 0.5^(1/3), d = 2R)
e' ideale.  L'interazione a due sfere aggiunge (r_cw/d)^2 = 0.157 al termine
della bead: 3.881 -> 4.26 (+9.7%).  Se la correzione si applica all'intero
attrito cw: 4.49 (+15.7%), che e' il limite superiore citato.  Il calcolo di
rumore usa 3.881, il valore conservativo; con 4.26 il limite di CV sale a
19.3%, con 4.49 a 19.9%.  Nessuno dei due arriva a 30%.

Uso:
    python codice/v73_rumore.py
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

import v73_contrappeso as c  # noqa: E402
import v73_taglio as t  # noqa: E402

DATI = os.path.join(ROOT, "dati")
TAB_IT = os.path.join(ROOT, "paper", "it", "tabelle")
TAB_EN = os.path.join(ROOT, "paper", "en", "tabelle")
FIG_IT = os.path.join(ROOT, "paper", "it", "figures")
FIG_EN = os.path.join(ROOT, "paper", "en", "figures")
FIG_ROOT_IT = os.path.join(ROOT, "figure", "it")
FIG_ROOT_EN = os.path.join(ROOT, "figure", "en")
for d in (DATI, TAB_IT, TAB_EN, FIG_IT, FIG_EN, FIG_ROOT_IT, FIG_ROOT_EN):
    os.makedirs(d, exist_ok=True)

LANG = os.environ.get("V72_FIG_LANG", "it")
PE_R = c.PE_R
CV_IDEALE = 1.0 / math.sqrt(math.pi * PE_R)       # 0.1784
RATIO = c.BETA_RATIO                              # 3.8811
R_CW = 0.5 ** (1.0 / 3.0)
D_REL = 2.0
CORR_DUE_SFERE = (R_CW / D_REL) ** 2              # 0.1575
TERM_ROT = 1.0
TERM_BEAD_ROT = R_CW ** 3
TERM_BEAD_TRANS = 0.75 * R_CW * D_REL ** 2
RATIO_BEAD = TERM_ROT + TERM_BEAD_ROT + TERM_BEAD_TRANS * (1.0 + CORR_DUE_SFERE)
RATIO_BEAD_TOT = TERM_ROT + (TERM_BEAD_ROT + TERM_BEAD_TRANS) * (1.0 + CORR_DUE_SFERE)
RATIO_UP = RATIO * (1.0 + CORR_DUE_SFERE)

R_TEMPI = (0.3, 0.4, 0.5, 0.6, 0.8)
D_T_CRITERIO = 1.0e-4


# ---------------------------------------------------------------------------
# 1. PRIMA PASSA ESATTA (inverse Gaussian), stabile in log-spazio
# ---------------------------------------------------------------------------
def _norm_cdf(z):
    return 0.5 * math.erfc(-z / math.sqrt(2.0))


def _norm_logcdf(z):
    from scipy.special import log_ndtr
    return float(log_ndtr(-z))


def ig_cdf(x, mean, cv):
    """CDF di una IG con media `mean` e CV `cv`, stabile anche nelle code.

    F(x) = Phi(z1) + exp(2 lambda/m) Phi(-z2),
    z1 = sqrt(lambda/x)(x/m - 1), z2 = sqrt(lambda/x)(x/m + 1),
    lambda = m / cv^2.
    Il secondo termine e' calcolato in log-spazio: exp(2 lambda/m) puo'
    essere enorme (10^60) ma il prodotto con Phi(-z2) e' finito.
    """
    if x <= 0.0:
        return 0.0
    m = float(mean)
    lam = m / (cv * cv)
    s = math.sqrt(lam / x)
    z1 = s * (x / m - 1.0)
    z2 = s * (x / m + 1.0)
    lt = 2.0 * lam / m + _norm_logcdf(z2)
    if lt > 700.0:
        return 1.0
    return min(1.0, max(0.0, _norm_cdf(z1) + math.exp(lt)))


def errore_appaiato(cv, ratio=RATIO, N=3, ottimizza=True):
    """Errore del test appaiato: N giri veloci contro N giri lenti.

    eps(soglia) = P(S_fast > soglia) + P(S_slow < soglia).
    Con la IG esatta eps non dipende da N; N=3 e' il protocollo.
    """
    m_fast = float(N)
    m_slow = float(N) * ratio
    if not ottimizza:
        th = N * math.sqrt(ratio)
        return {"eps": (1.0 - ig_cdf(th, m_fast, cv)) + ig_cdf(th, m_slow, cv),
                "soglia": th, "soglia_geo": N * math.sqrt(ratio)}

    def f(x):
        return (1.0 - ig_cdf(x, m_fast, cv)) + ig_cdf(x, m_slow, cv)

    from scipy.optimize import minimize_scalar
    lo = m_fast * (1.0 + 1e-9)
    hi = m_slow * (1.0 - 1e-9)
    res = minimize_scalar(f, bounds=(lo, hi), method="bounded",
                          options={"xatol": 1e-12})
    th_geo = N * math.sqrt(ratio)
    return {"eps": float(res.fun), "soglia": float(res.x),
            "eps_geo": float(f(th_geo)), "soglia_geo": th_geo,
            "eps_fast": float(1.0 - ig_cdf(th_geo, m_fast, cv)),
            "eps_slow": float(ig_cdf(th_geo, m_slow, cv))}


def errore_totale(cv, ratio=RATIO, N=3):
    """Test 'durata totale': N cicli con latch contro 2N giri veloci.

    Media del latch: N(1+ratio); media senza latch: 2N.  La separazione
    e' (1+ratio)/2 = 2.44, non ratio.  Anche questo errore e' indipendente
    da N per la IG esatta.
    """
    m_latch = float(N) * (1.0 + ratio)
    m_no = 2.0 * float(N)

    def f(x):
        return (1.0 - ig_cdf(x, m_no, cv)) + ig_cdf(x, m_latch, cv)

    from scipy.optimize import minimize_scalar
    res = minimize_scalar(f, bounds=(m_no * (1.0 + 1e-9),
                                     m_latch * (1.0 - 1e-9)),
                          method="bounded", options={"xatol": 1e-12})
    return {"eps": float(res.fun), "soglia": float(res.x)}


def limite_cv(ratio=RATIO, N=3):
    """CV massimo per eps_tot <= 10^-4 (test appaiato e test totale)."""
    from scipy.optimize import brentq

    def g(cv):
        return math.log10(errore_appaiato(cv, ratio, N)["eps"]) + 4.0

    def h(cv):
        return math.log10(errore_totale(cv, ratio, N)["eps"]) + 4.0

    out = {}
    for nome, fun, hi in (("appaiato", g, 0.6), ("totale", h, 0.9)):
        try:
            out[nome] = float(brentq(fun, 0.05, hi, xtol=1e-10))
        except ValueError:
            out[nome] = None
    return out


# ---------------------------------------------------------------------------
# 2. TEMPI E GEOMETRIA
# ---------------------------------------------------------------------------
def tempi(R_um):
    """T360 no/cw, T720 (somma), T_meas(3), dal motore a 800 nm, z=2 um."""
    rif = c.ciclo(R_um * 1e-6, 800, 2.0e-6, mu=c.SANGUE["mu_p"])
    return {
        "R_um": R_um,
        "T360_no_s": rif["T720_no"] / 2.0,
        "T360_cw_s": rif["T720_cw"] / 2.0,
        "T720_s": rif["T720_half"],
        "T_meas3_s": 3.0 * rif["T720_half"],
        "omega_no": rif["omega_no"],
        "omega_cw": rif["omega_cw"],
        "ratio": rif["omega_no"] / rif["omega_cw"],
    }


def specifica_latch():
    """Requisiti del latch, con stato dichiarato."""
    t_ref = tempi(0.5)
    lavoro_giro_kT = 2.0 * math.pi * PE_R
    tau_bead = (4.0 / 3.0 * math.pi * (R_CW * 0.5e-6) ** 3 * 2000.0) / (
        6.0 * math.pi * c.SANGUE["mu_p"] * R_CW * 0.5e-6)
    lim = limite_cv()
    return {
        "ratio": {
            "valore": RATIO,
            "target": 3.9,
            "corretto_bead": RATIO_BEAD,
            "corretto_bead_tot": RATIO_BEAD_TOT,
            "limite_superiore": RATIO_UP,
            "stato": "geometria ideale verificata; correzione idrodinamica "
                     "non validata",
        },
        "cv": {
            "ideale": CV_IDEALE,
            "limite_10e4_appaiato": lim["appaiato"],
            "limite_10e4_totale": lim["totale"],
            "target": 0.30,
            "stato": "il target 30% NON e' supportato: a CV=30.2% "
                     "l'errore esatto e' 2.07e-2; il limite e' 17.8%",
        },
        "tempo_commutazione": {
            "requisito_s": 0.1 * t_ref["T360_cw_s"],
            "rilassamento_viscoso_bead_s": tau_bead,
            "stato": "il vincolo viscoso e' ~10^-7 s; l'attuazione non e' "
                     "modellata",
        },
        "energia": {
            "lavoro_ottico_per_giro_kT": lavoro_giro_kT,
            "barriera_latch_kT": [10.0, 60.0],
            "stato": "finestra 10-60 kT dentro i 62.8 kT per giro; non "
                     "dimostrata",
        },
        "fatica": {
            "cicli_richiesti": 1.0e6,
            "durata_giorni_a_T720": 1.0e6 * t_ref["T720_s"] / 86400.0,
            "stato": "non misurata",
        },
        "protocollo": {
            "tipo": "appaiato: 3 giri veloci contro 3 giri lenti",
            "soglia": "geometrica; l'ottimizzata guadagna il 6%",
            "N": 3,
            "nota": "la statistica non migliora con N; i 3 cicli servono "
                    "per i sistematici",
        },
    }


def sintesi():
    lim = limite_cv()
    casi = []
    for cv in (0.15, CV_IDEALE, 0.20, 0.25, 0.302):
        p = errore_appaiato(cv)
        q = errore_totale(cv)
        casi.append({
            "cv": cv,
            "sigma_theta_360_deg": math.degrees(cv * 2.0 * math.pi),
            "eps_appaiato": p["eps"],
            "eps_appaiato_geo": p["eps_geo"],
            "soglia": p["soglia"],
            "eps_fast": p["eps_fast"],
            "eps_slow": p["eps_slow"],
            "eps_totale": q["eps"],
        })
    # errore a 30.2% esplicito per la ritrattazione
    e302 = errore_appaiato(0.302)
    return {
        "descrizione": "V7.3: rumore del latch con prima passa esatta; "
                       "test appaiato; ritrattazione del 30.2%",
        "input": {
            "Pe_r": PE_R, "CV_ideale": CV_IDEALE, "ratio": RATIO,
            "R_cw": R_CW, "d_rel": D_REL,
            "correzione_due_sfere": CORR_DUE_SFERE,
            "ratio_bead": RATIO_BEAD,
            "ratio_bead_tot": RATIO_BEAD_TOT,
            "ratio_limite_superiore": RATIO_UP,
            "N_cicli": 3,
            "criterio_eps": D_T_CRITERIO,
        },
        "modello": {
            "distribuzione": "inverse Gaussian (prima passa di un giro)",
            "somma": "IG con media N m e shape N m/cv^2",
            "test": "appaiato: N giri veloci vs N giri lenti",
            "invarianza_N": "l'errore e' indipendente da N; la shape e la "
                            "soglia scalano insieme",
        },
        "casi": casi,
        "limite_cv": lim,
        "ritrattazione_30_2": {
            "cv": 0.302,
            "eps_appaiato": e302["eps"],
            "eps_totale": errore_totale(0.302)["eps"],
            "nota": "30.2% era l'approssimazione log-normale a CV uguali; "
                    "la IG esatta da' 2.07e-2, non 10^-4",
        },
        "tempi": [tempi(R) for R in R_TEMPI],
        "specifica_latch": specifica_latch(),
        "ritrattazioni": [
            "CV fino a 30.2% per 10^-4: ritrattato; il limite esatto e' "
            "17.83% (appaiato), 11.58% (durata totale)",
            "8.8e-7 con tre cicli: ritrattato; l'errore esatto e' 9.7e-5 e "
            "non dipende da N",
            "2.9e-3 per il giro lento: ritrattato; la coda esatta e' 6.8e-5",
            "due o tre cicli bastano: vero per i sistematici, non per la "
            "statistica; la statistica non migliora con N",
            "durata totale del ciclo: non basta; l'errore e' 1.1e-2 a CV "
            "ideale, il test appaiato e' obbligatorio",
        ],
        "verdetto": "il latch ideale e' al limite: CV 17.8% satura il budget "
                    "10^-4; ogni rumore aggiuntivo lo porta sopra.  Il latch "
                    "non e' dimostrato e non ha margine.",
    }


# ---------------------------------------------------------------------------
# 3. TABELLE
# ---------------------------------------------------------------------------
def tabella_rumore(lingua="it"):
    s = sintesi()
    it = lingua == "it"
    if it:
        intest = (r"$CV$ & $\sigma_\theta(360^\circ)$ & "
                  r"$\epsilon_{\rm pair}$ & $\epsilon_{\rm pair}^{\rm geo}$ & "
                  r"$\epsilon_{\rm tot}$ \\"
                  "\n" r" & deg & & & \\")
        nota = (r"Nota: prima passa esatta (inverse Gaussian).  "
                r"$\epsilon_{\rm pair}$ confronta 3 giri veloci con 3 giri "
                r"lenti; $\epsilon_{\rm tot}$ confronta la durata totale.  "
                r"L'errore e' indipendente dal numero di cicli: la shape "
                r"scala con $N$.  A $CV=17.84\%$ (ideale) il test appaiato "
                r"da' $1.01\times10^{-4}$, appena sopra il budget; il limite "
                r"per $10^{-4}$ e' $CV=17.83\%$ (appaiato) e $11.58\%$ "
                r"(totale).  Il valore $30.2\%$ della V7.2 era "
                r"l'approssimazione log-normale e va ritrattato: a "
                r"$CV=30.2\%$ l'errore esatto e' $2.1\times10^{-2}$.")
    else:
        intest = (r"$CV$ & $\sigma_\theta(360^\circ)$ & "
                  r"$\epsilon_{\rm pair}$ & $\epsilon_{\rm pair}^{\rm geo}$ & "
                  r"$\epsilon_{\rm tot}$ \\"
                  "\n" r" & deg & & & \\")
        nota = (r"Note: exact first passage (inverse Gaussian).  "
                r"$\epsilon_{\rm pair}$ compares 3 fast turns with 3 slow "
                r"turns; $\epsilon_{\rm tot}$ compares the total duration.  "
                r"The error is independent of the number of cycles: the "
                r"shape scales with $N$.  At $CV=17.84\%$ (ideal) the paired "
                r"test gives $1.01\times10^{-4}$, just above budget; the "
                r"$10^{-4}$ limit is $CV=17.83\%$ (paired) and $11.58\%$ "
                r"(total).  The V7.2 value $30.2\%$ was the log-normal "
                r"approximation and is retracted: at $CV=30.2\%$ the exact "
                r"error is $2.1\times10^{-2}$.")
    lines = [r"\begin{tabular}{@{}rrrrr@{}}", r"\toprule", intest, r"\midrule"]
    for x in s["casi"]:
        lines.append(r"%.3f & %.1f & $%.2f\times10^{%d}$ & $%.2f\times10^{%d}$ & "
                     r"$%.2f\times10^{%d}$ \\"
                     % (x["cv"], x["sigma_theta_360_deg"],
                        _mant(x["eps_appaiato"])[0], _exp(x["eps_appaiato"]),
                        _mant(x["eps_appaiato_geo"])[0], _exp(x["eps_appaiato_geo"]),
                        _mant(x["eps_totale"])[0], _exp(x["eps_totale"])))
    lines += [
        r"\midrule",
        r"\multicolumn{5}{p{12.6cm}}{\scriptsize " + nota + r"} \\",
        r"\bottomrule", r"\end{tabular}",
    ]
    return "\n".join(lines) + "\n"


def tabella_tempi(lingua="it"):
    s = sintesi()
    it = lingua == "it"
    intest = (r"$R$ & $T_{360}^{\rm no}$ & $T_{360}^{\rm cw}$ & "
              r"$T_{720}$ & $T_{\rm meas}(3)$ \\"
              "\n" r"$\mu$m & s & s & s & s \\")
    if not it:
        intest = (r"$R$ & $T_{360}^{\rm no}$ & $T_{360}^{\rm cw}$ & "
                  r"$T_{720}$ & $T_{\rm meas}(3)$ \\"
                  "\n" r"$\mu$m & s & s & s & s \\")
    lines = [r"\begin{tabular}{@{}rrrrr@{}}", r"\toprule", intest, r"\midrule"]
    for x in s["tempi"]:
        lines.append(r"%.1f & %.2f & %.2f & %.2f & %.1f \\"
                     % (x["R_um"], x["T360_no_s"], x["T360_cw_s"],
                        x["T720_s"], x["T_meas3_s"]))
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines) + "\n"


def tabella_latch(lingua="it"):
    s = sintesi()["specifica_latch"]
    it = lingua == "it"
    if it:
        intest = (r"Requisito & Valore & Stato \\")
        righe = [
            ("rapporto", r"$%.3f$ (ideale); $%.2f$--$%.2f$ corretto"
             % (s["ratio"]["valore"], s["ratio"]["corretto_bead"],
                s["ratio"]["limite_superiore"]), "geometria; correzione non validata"),
            ("$CV$", r"$%.1f\%%$ ideale; limite $%.2f\%%$"
             % (100 * s["cv"]["ideale"], 100 * s["cv"]["limite_10e4_appaiato"]),
             "nessun margine: il target 30\\% non e' supportato"),
            ("tempo di commutazione", r"$<%.2f$ s" % s["tempo_commutazione"]["requisito_s"],
             "viscoso $\\sim10^{-7}$ s; attuazione non modellata"),
            ("energia", r"$10$--$60$ kT" , "dentro i 62.8 kT per giro; non dimostrata"),
            ("fatica", r"$10^6$ cicli ($%.0f$ giorni)" % s["fatica"]["durata_giorni_a_T720"],
             "non misurata"),
            ("protocollo", "appaiato, 3+3 giri", "obbligatorio; la statistica non migliora con $N$"),
        ]
    else:
        intest = (r"Requirement & Value & Status \\")
        righe = [
            ("ratio", r"$%.3f$ (ideal); $%.2f$--$%.2f$ corrected"
             % (s["ratio"]["valore"], s["ratio"]["corretto_bead"],
                s["ratio"]["limite_superiore"]), "geometry; correction not validated"),
            ("$CV$", r"$%.1f\%%$ ideal; limit $%.2f\%%$"
             % (100 * s["cv"]["ideale"], 100 * s["cv"]["limite_10e4_appaiato"]),
             "no margin: the 30\\% target is not supported"),
            ("switching time", r"$<%.2f$ s" % s["tempo_commutazione"]["requisito_s"],
             "viscous $\\sim10^{-7}$ s; actuation not modelled"),
            ("energy", r"$10$--$60$ kT", "within 62.8 kT per turn; not demonstrated"),
            ("fatigue", r"$10^6$ cycles ($%.0f$ days)" % s["fatica"]["durata_giorni_a_T720"],
             "not measured"),
            ("protocol", "paired, 3+3 turns", "mandatory; statistics do not improve with $N$"),
        ]
    lines = [r"\begin{tabular}{@{}p{2.3cm}p{5.6cm}p{6.2cm}@{}}",
             r"\toprule", intest, r"\midrule"]
    for a, b, c3 in righe:
        lines.append(r"%s & %s & %s \\" % (a, b, c3))
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines) + "\n"


def _mant(x):
    if x <= 0:
        return 0.0, 0
    e = int(math.floor(math.log10(x)))
    return x / 10.0 ** e, e


def _exp(x):
    return _mant(x)[1]


def scrivi(nome, testo):
    with open(nome, "w", encoding="utf-8") as fh:
        fh.write(testo)
    print("  scritto %-54s %5.1f kB" % (nome, os.path.getsize(nome) / 1024.0))


# ---------------------------------------------------------------------------
# 4. FIGURA
# ---------------------------------------------------------------------------
def figura(lingua="it"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    it = lingua == "it"
    s = sintesi()
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.9))

    cvs = np.linspace(0.08, 0.40, 400)
    ep = np.array([errore_appaiato(float(x), RATIO)["eps"] for x in cvs])
    et = np.array([errore_totale(float(x), RATIO)["eps"] for x in cvs])
    ax = axes[0]
    ax.semilogy(cvs, ep, color="#2166ac", lw=1.9,
                label=("appaiato 3+3" if it else "paired 3+3"))
    ax.semilogy(cvs, et, color="#b2182b", lw=1.7, ls="--",
                label=("durata totale" if it else "total duration"))
    ax.axhline(1e-4, color="#444444", ls=":", lw=1.0)
    ax.axvline(CV_IDEALE, color="#1b7837", ls="-.", lw=1.0)
    ax.axvline(0.302, color="#999999", ls=":", lw=1.0)
    ax.text(CV_IDEALE + 0.004, 3e-8, ("CV ideale\n17.8%" if it else
                                     "ideal CV\n17.8%"), fontsize=7.5,
            color="#1b7837")
    ax.text(0.304, 2e-6, ("30.2%\n(ritrattato)" if it else
                          "30.2%\n(retracted)"), fontsize=7.5, color="#666666")
    ax.text(0.088, 1.4e-4, r"$10^{-4}$", fontsize=8, color="#444444")
    ax.set_xlabel("$CV$")
    ax.set_ylabel(("errore totale di classificazione" if it else
                   "total classification error"))
    ax.set_xlim(0.08, 0.40)
    ax.set_ylim(1e-9, 0.5)
    ax.grid(True, alpha=0.22, lw=0.5)
    ax.legend(fontsize=8, loc="lower right")
    ax.set_title(("(a) errore contro CV" if it else "(a) error versus CV"),
                 fontsize=9.5)

    ax = axes[1]
    ratios = np.linspace(3.0, 5.2, 400)
    e = np.array([errore_appaiato(CV_IDEALE, float(x))["eps"] for x in ratios])
    ax.semilogy(ratios, e, color="#2166ac", lw=1.9)
    ax.axhline(1e-4, color="#444444", ls=":", lw=1.0)
    for r, col, lab in ((RATIO, "#b2182b", "3.881"),
                        (RATIO_BEAD, "#1b7837", "4.26"),
                        (RATIO_UP, "#666666", "4.49")):
        ax.axvline(r, color=col, ls="--", lw=1.0)
        ax.text(r + 0.02, 2e-7, lab, color=col, fontsize=7.5, rotation=90)
    ax.set_xlabel(("rapporto dei tempi $r$" if it else "time ratio $r$"))
    ax.set_ylabel(("errore appaiato a CV ideale" if it else
                   "paired error at ideal CV"))
    ax.set_xlim(3.0, 5.2)
    ax.set_ylim(1e-9, 1.0)
    ax.grid(True, alpha=0.22, lw=0.5)
    ax.set_title(("(b) effetto del rapporto" if it else "(b) effect of the ratio"),
                 fontsize=9.5)

    fig.suptitle(("Latch: la prima passa esatta non perdona; il 30.2% e' "
                  "un artefatto gaussiano" if it else
                  "Latch: the exact first passage is unforgiving; 30.2% is a "
                  "Gaussian artefact"), fontsize=10)
    fig.tight_layout()
    for d in ((FIG_IT, FIG_ROOT_IT) if it else (FIG_EN, FIG_ROOT_EN)):
        p = os.path.join(d, "v73_rumore_%s.pdf" % lingua)
        fig.savefig(p, bbox_inches="tight")
        fig.savefig(os.path.join(d, "v73_rumore_%s.png" % lingua),
                    bbox_inches="tight", dpi=220)
        print("  scritto %-54s %5.1f kB" % (p, os.path.getsize(p) / 1024.0))
    plt.close(fig)


# ---------------------------------------------------------------------------
def main():
    s = sintesi()
    p = os.path.join(DATI, "V73_rumore.json")
    with open(p, "w", encoding="utf-8") as fh:
        json.dump(s, fh, indent=1, ensure_ascii=False, default=float)
    print("  scritto %-54s %5.1f kB" % (p, os.path.getsize(p) / 1024.0))
    scrivi(os.path.join(TAB_IT, "tab_v73_rumore.tex"), tabella_rumore("it"))
    scrivi(os.path.join(TAB_EN, "tab_v73_rumore.tex"), tabella_rumore("en"))
    scrivi(os.path.join(TAB_IT, "tab_v73_tempi.tex"), tabella_tempi("it"))
    scrivi(os.path.join(TAB_EN, "tab_v73_tempi.tex"), tabella_tempi("en"))
    scrivi(os.path.join(TAB_IT, "tab_v73_latch.tex"), tabella_latch("it"))
    scrivi(os.path.join(TAB_EN, "tab_v73_latch.tex"), tabella_latch("en"))
    figura("it")
    figura("en")

    print("\n-- correzione idrodinamica")
    print("   (r_cw/d)^2 = %.4f" % CORR_DUE_SFERE)
    print("   ratio ideale        = %.4f" % RATIO)
    print("   ratio bead          = %.4f  (%+.1f%%)"
          % (RATIO_BEAD, 100.0 * (RATIO_BEAD / RATIO - 1.0)))
    print("   ratio bead totale   = %.4f  (%+.1f%%)"
          % (RATIO_BEAD_TOT, 100.0 * (RATIO_BEAD_TOT / RATIO - 1.0)))
    print("   limite superiore    = %.4f  (%+.1f%%)"
          % (RATIO_UP, 100.0 * (RATIO_UP / RATIO - 1.0)))
    print("\n-- errore appaiato (3+3 giri)")
    for x in s["casi"]:
        print("   CV=%.3f  sigma=%.1f deg  eps=%.3e  geo=%.3e  tot=%.3e"
              % (x["cv"], x["sigma_theta_360_deg"], x["eps_appaiato"],
                 x["eps_appaiato_geo"], x["eps_totale"]))
    print("\n-- limiti CV per 10^-4")
    print("   appaiato: %.4f  totale: %.4f"
          % (s["limite_cv"]["appaiato"], s["limite_cv"]["totale"]))
    print("   a CV=30.2%%: appaiato %.3e  totale %.3e"
          % (s["ritrattazione_30_2"]["eps_appaiato"],
             s["ritrattazione_30_2"]["eps_totale"]))
    print("\n-- specifica del latch")
    sp = s["specifica_latch"]
    print("   ratio %.3f  CV ideale %.3f  CV limite %.3f"
          % (sp["ratio"]["valore"], sp["cv"]["ideale"],
             sp["cv"]["limite_10e4_appaiato"]))
    print("   commutazione < %.2f s (viscoso %.2e s)"
          % (sp["tempo_commutazione"]["requisito_s"],
             sp["tempo_commutazione"]["rilassamento_viscoso_bead_s"]))
    print("   energia 10-60 kT dentro %.1f kT per giro"
          % sp["energia"]["lavoro_ottico_per_giro_kT"])
    print("\n   %s" % s["verdetto"])
    return s


if __name__ == "__main__":
    main()
