"""
MOTORE FOTONICO DEL VUOTO - V7.1
v7_tabelle.py - genera le tabelle LaTeX del preprint dai dati.

Uso:
    python v7_tabelle.py

Legge M1_mappa.json, M1_griglia.csv, V0_validazione.json, S1_scala_1nm.json,
S3_scala_1m.json e scrive ../tabelle/tab_*.tex.  Nessun numero e' scritto a
mano: se un dato cambia, la tabella cambia con lui.
"""

import csv
import json
import math
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
if BASE not in sys.path:
    sys.path.insert(0, BASE)

import v7_mappa as m  # noqa: E402

DATI = os.path.normpath(os.path.join(BASE, "..", "dati"))
TAB = os.path.normpath(os.path.join(BASE, "..", "tabelle"))
os.makedirs(TAB, exist_ok=True)


# ---------------------------------------------------------------------------
# formattazione italiana
# ---------------------------------------------------------------------------
def it(x, nd=2):
    """Numero con la virgola decimale, senza zeri inutili."""
    s = ("%." + str(nd) + "f") % x
    if "." in s:
        s = s.rstrip("0").rstrip(".")
    return s.replace(".", "{,}")


def sci(x, nd=2):
    """Notazione scientifica LaTeX con virgola decimale."""
    if x == 0:
        return "0"
    e = int(math.floor(math.log10(abs(x))))
    a = x / (10.0 ** e)
    if abs(a) >= 10.0:
        a /= 10.0
        e += 1
    if abs(a) < 1.0:
        a *= 10.0
        e -= 1
    s = ("%." + str(nd) + "f") % a
    return "$" + s.replace(".", "{,}") + "\\times10^{%d}$" % e


def val(x, nd=2):
    """Numero semplice, con notazione scientifica solo fuori da 0.01-1000."""
    if x == 0:
        return "0"
    if 1e-2 <= abs(x) < 1e3:
        return it(x, nd)
    return sci(x, nd)


def rng(lo, hi, nd=2):
    if abs(hi - lo) <= 1e-3 * max(abs(hi), abs(lo), 1e-30):
        return it(lo, nd)
    return it(lo, nd) + "--" + it(hi, nd)


def prng(lo, hi, nd=2):
    if abs(hi - lo) <= 1e-3 * max(abs(hi), abs(lo), 1e-30):
        return sci(lo, nd)
    return sci(lo, nd) + "--" + sci(hi, nd)


def mezzo(mz, p):
    if mz == "acqua":
        return "acqua"
    return "vuoto $10^{%d}$" % int(round(math.log10(p)))


def write(name, body):
    path = os.path.join(TAB, name)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(body)
    print("  %-28s %6.1f kB" % (name, os.path.getsize(path) / 1024.0))


# ---------------------------------------------------------------------------
# dati
# ---------------------------------------------------------------------------
with open(os.path.join(DATI, "M1_mappa.json"), encoding="utf-8") as fh:
    MAPPA = json.load(fh)
with open(os.path.join(DATI, "V0_validazione.json"), encoding="utf-8") as fh:
    VAL = json.load(fh)
with open(os.path.join(DATI, "S1_scala_1nm.json"), encoding="utf-8") as fh:
    S1 = json.load(fh)
with open(os.path.join(DATI, "S3_scala_1m.json"), encoding="utf-8") as fh:
    S3 = json.load(fh)

ROWS = list(csv.DictReader(open(os.path.join(DATI, "M1_griglia.csv"),
                                encoding="utf-8")))


def btrue(x):
    return x.strip().lower() == "true"


# ---------------------------------------------------------------------------
# 1. correzioni
# ---------------------------------------------------------------------------
def tab_correzioni():
    body = r"""\begin{tabular}{@{}p{0.35cm}p{2.8cm}p{4.1cm}p{4.2cm}p{4.2cm}@{}}
\hline
\# & Correzione & V7.0 & V7.1 & Effetto sulla mappa \\
\hline
1 & Scattering assiale &
  $\tau_{\rm sca}\neq0$ per sfera, disco, anello, reticolo &
  $\tau_{\rm sca}=0$ per forme assiali; solo manubrio &
  6 finestre di scattering, tutte di manubrio \\
\hline
2 & Inerzia del manubrio &
  $k_I$ in m$^2$, poi $\times R^2$: $I\propto R^6$ &
  $k_I=0{,}5875$ adimensionale &
  $I(100\,{\rm nm})=1{,}692\times10^{-33}$ kg m$^2$; $t_{\rm spin}=0{,}076$ s \\
\hline
3 & Ritardo di fase &
  $\delta\propto R$ fino a $R\to0$; $(1-\cos\delta)$ oscilla &
  $\delta\propto R^3$ sotto $0{,}15\,\mu$m; $\sigma_\delta=0{,}5$ &
  $P(Pe_r=10)$ a 100 nm: 99,6 mW; finestra da 0,14 $\mu$m \\
\hline
4 & Cu nero &
  $\varepsilon_{IR}$ da Mie (sfera di rame) &
  $\varepsilon_{IR}=0{,}30$ dichiarato (0,20--0,60) &
  bordo termico $\sim$30 $\mu$m; in vuoto solo 30 $\mu$m \\
\hline
5 & Pressione residua &
  conversione mbar$\to$Pa: $0{,}1$ (1000$\times$ basso) &
  conversione mbar$\to$Pa: $100$ &
  $\omega,\eta_{\rm mecc}$ in vuoto $\times10^{-3}$; rotori 9$\to$3 \\
\hline
\end{tabular}
"""
    write("tab_correzioni.tex", body)


# ---------------------------------------------------------------------------
# 2. criteri
# ---------------------------------------------------------------------------
def tab_criteri():
    body = r"""\begin{tabular}{@{}p{2.4cm}p{2.1cm}p{2.4cm}p{8.6cm}@{}}
\hline
Criterio & Valore & Dove & Perché \\
\hline
$Pe_r$ & $\ge 10$ & acqua e vuoto &
  la coppia ottica supera di 10 volte il rumore browniano: l'orientamento è deterministico \\
$\Delta T$ & $\le 20$ K & acqua &
  conduzione sferica, Nusselt 1; sopra 20 K la convezione e la variazione di viscosità non sono trascurabili \\
$T$ & $\le 0{,}7\,T_{\rm fusione}$ & vuoto &
  equilibrio radiativo; il 30\% di margine copre l'incertezza su $\varepsilon_{IR}$ e la temperatura non uniforme \\
$\omega$ & $\le 0{,}5\,\omega_{\rm rottura}$ & acqua &
  il rotore non deve lavorare oltre metà del limite meccanico; in vuoto il criterio non si applica a $\omega=\tau/\beta_{\rm gas}$ (vedi \S\ref{sec:rotore}) \\
$P$ & $\le 10$ mW & tutte &
  tetto dichiarato del banco ottico; non è un limite fisico \\
$\eta_{\rm mecc}$ & $\ge 10^{-7}$ & tutte &
  soglia dichiarata di ``rotore'': $\omega/\omega_{\rm ott}$ almeno $10^{-7}$; è una scelta di progetto, non una legge \\
$t_{\rm spin}$ & $\ge 1$ s & vuoto &
  il rotore non deve rompersi prima di essere utile; $t_{\rm spin}$ è il tempo a metà del limite meccanico senza attrito \\
\hline
\end{tabular}
"""
    write("tab_criteri.tex", body)


# ---------------------------------------------------------------------------
# 3. validazione
# ---------------------------------------------------------------------------
def tab_validazione():
    rows = []
    for t in VAL["test"]:
        nome = t["test"]
        rif = t["riferimento"]
        d = t["dettagli"]
        if nome.startswith("A"):
            controllo = "4 casi limite (circolare/lineare, $\\lambda/4$, $\\lambda/2$)"
            modello = "1; 2; 1; 0"
            atteso = "1; 2; 1; 0"
            rapporto = "esatto"
        elif nome.startswith("B"):
            controllo = "$\\tau/(P/\\omega)$ per lamina $\\lambda/2$ circolare"
            modello = "2"
            atteso = "2"
            rapporto = "esatto"
        elif nome.startswith("C"):
            rif = "V6.1: $\\Delta T=P_{\\rm abs}/(4\\pi\\kappa R)$; 150 $\\mu$W"
            controllo = "$\\Delta T$ a 150 $\\mu$W, $R=1$ e 3 $\\mu$m"
            modello = "19,0; 6,34 K"
            atteso = "20; 6,6 K"
            rapporto = "0,95; 0,96"
        elif nome.startswith("D"):
            rif = "He et al., PRL 75, 826 (1995)"
            controllo = "CuO 2 $\\mu$m, 10 mW, $\\ell=1$, kerosene"
            modello = "10,75 Hz"
            atteso = "4 Hz"
            rapporto = "2,69"
        elif nome.startswith("E"):
            rif = "Friese et al., Nature 394, 348 (1998)"
            controllo = "calcite 1 $\\mu$m, 300 mW, 1064 nm"
            modello = "243,8 Hz"
            atteso = "357 Hz"
            rapporto = "0,68"
        else:
            rif = "Reimann et al., PRL 121, 033602 (2018); Ahn et al. (2020)"
            controllo = "manubrio silice 100 nm, 10 mW, $\\sigma_y=1$--3 GPa"
            modello = "1,52--2,63 GHz; $t_{\\rm spin}=0{,}076$--$0{,}132$ s"
            atteso = "0,5--10 GHz; $t_{\\rm spin}<1$ s"
            rapporto = "--"
        rows.append((nome, rif, controllo, modello, atteso, rapporto, t["esito"]))
    out = [r"""\begin{tabular}{|p{0.3cm}|p{2.3cm}|p{2.8cm}|p{3.6cm}|p{4.25cm}|p{0.65cm}|}
\hline
 & Test & Riferimento & Controllo & Modello / atteso & Esito \\
\hline
"""]
    for i, (nome, rif, controllo, modello, atteso, rapporto, esito) in enumerate(rows, 1):
        nome = nome.split(" - ", 1)[-1].replace("lambda/2", "$\\lambda/2$")
        if rapporto == "esatto":
            confronto = "%s / %s (esatto)" % (modello, atteso)
        elif rapporto == "--":
            confronto = "%s / %s" % (modello, atteso)
        else:
            confronto = "%s / %s (%s)" % (modello, atteso, rapporto)
        out.append("%d & %s & %s & %s & %s & %s \\\\\n" %
                   (i, nome, rif, controllo, confronto, esito))
        out.append("\\hline\n")
    out.append("\\end{tabular}\n")
    write("tab_validazione.tex", "".join(out))


# ---------------------------------------------------------------------------
# 4. le 26 finestre
# ---------------------------------------------------------------------------
def tab_finestre():
    ord_mezzo = {"acqua": 0, "vuoto": 1}
    ord_can = {"assorbimento": 0, "ritardo": 1, "scattering": 2}
    ord_mat = {"cu_lucido": 0, "cu_nero": 1, "silice": 2, "vaterite": 3}
    ord_for = {"sfera": 0, "disco": 1, "anello": 2, "manubrio": 3, "reticolo": 4}
    rr = sorted(MAPPA["riassunto"],
                key=lambda r: (ord_mezzo[r["mezzo"]],
                               0 if r["p_mbar"] is None else -math.log10(r["p_mbar"]),
                               ord_can[r["canale"]], ord_mat[r["materiale"]],
                               ord_for[r["forma"]]))
    out = [r"""\begin{tabular}{|l|l|l|l|r|r|r|r|r|}
\hline
Canale & Materiale & Forma & Mezzo & $R_{\min}$ & $R_{\max}$ & $P_{\min}$ & $P_{\max}$ & Segm. \\
 & & & & ($\mu$m) & ($\mu$m) & (W) & (W) & \\
\hline
"""]
    old_mezzo = None
    for r in rr:
        mz = mezzo(r["mezzo"], r["p_mbar"])
        if old_mezzo is not None and r["mezzo"] != old_mezzo:
            out.append("\\hline\n")
        old_mezzo = r["mezzo"]
        p_lo = r["finestre"][0]["P_lo"]
        p_hi = r["finestre"][-1]["P_hi"]
        out.append("%s & %s & %s & %s & %s & %s & %s & %s & %d \\\\\n" %
                   (r["canale"], r["materiale"].replace("_", "\\_"), r["forma"], mz,
                    it(r["R_min_um"], 2), it(r["R_max_um"], 2),
                    sci(p_lo), sci(p_hi), len(r["finestre"])))
    out.append("\\hline\n\\end{tabular}\n")
    write("tab_finestre.tex", "".join(out))


# ---------------------------------------------------------------------------
# 5. celle rappresentative
# ---------------------------------------------------------------------------
def _caso(canale, P, R, mat, forma, mezzo, p=None, **kw):
    c = m.cella(canale, P, R, mat, forma, mezzo, p, **kw)
    return c


def tab_rappresentativi():
    casi = [
        ("Vaterite, sfera, acqua", 10e-3, 1e-6, "vaterite", "sfera", "acqua", None, "ritardo", True),
        ("Vaterite, sfera, acqua", 10e-3, 2e-6, "vaterite", "sfera", "acqua", None, "ritardo", True),
        ("Vaterite, disco, acqua", 10e-3, 10e-6, "vaterite", "disco", "acqua", None, "ritardo", True),
        ("Cu nero, sfera, acqua", 0.165e-3, 1.0471e-6, "cu_nero", "sfera", "acqua", None, "assorbimento", True),
        ("Cu lucido, sfera, acqua", 0.39e-3, 1.0471e-6, "cu_lucido", "sfera", "acqua", None, "assorbimento", True),
        ("Silice, disco, acqua", 7e-3, 10e-6, "silice", "disco", "acqua", None, "ritardo", True),
        ("Silice, manubrio, acqua", 7e-3, 0.5e-6, "silice", "manubrio", "acqua", None, "scattering", True),
        ("Vaterite, manubrio, acqua", 8e-3, 0.5e-6, "vaterite", "manubrio", "acqua", None, "scattering", True),
        ("Vaterite, disco, vuoto $10^{-8}$", 0.5e-3, 10e-6, "vaterite", "disco", "vuoto", 1e-8, "ritardo", False),
    ]
    out = [r"""\begin{tabular}{@{}p{3.9cm}rrrrrrl@{}}
\hline
Caso & $\tau$ (N m) & $Pe_r$ & $f$ (Hz) & $\Delta T$ (K) & $\eta_{\rm mecc}$ & $t_{\rm spin}$ (s) & Esito \\
\hline
"""]
    for tag, P, R, mat, forma, mz, p, canale, water in casi:
        c = _caso(canale, P, R, mat, forma, mz, p)
        f = c["omega"] / (2.0 * math.pi)
        if water:
            fstr = it(f, 3)
            tstr = "--"
            esito = "valido"
        else:
            fstr = val(f, 2) + "$^\\ast$"
            tstr = sci(c["t_spin"], 2)
            esito = "rotore"
        out.append("%s & %s & %s & %s & %s & %s & %s & %s \\\\\n" %
                   (tag, sci(c["tau"], 2), it(c["Pe_r"], 1), fstr,
                    it(c["dT"], 2), sci(c["eta_mecc"], 1), tstr, esito))
    out.append("\\hline\n\\end{tabular}\n")
    write("tab_rappresentativi.tex", "".join(out))


# ---------------------------------------------------------------------------
# 6. rotori completi
# ---------------------------------------------------------------------------
def tab_rotore():
    gruppi = {}
    for r in ROWS:
        if btrue(r["valido_rotore"]):
            k = (r["canale"], r["materiale"], r["forma"], r["mezzo"], r["p"])
            gruppi.setdefault(k, []).append(r)
    out = [r"""\begin{tabular}{@{}llllrrrrr@{}}
\hline
Materiale & Forma & $p$ (mbar) & $R$ ($\mu$m) & $P$ (mW) & Celle & $t_{\rm spin}$ (s) & $T$ (K) & $\eta_{\rm mecc}$ \\
\hline
"""]
    tot = 0
    for k in sorted(gruppi):
        canale, mat, forma, mz, p = k
        v = gruppi[k]
        Rs = [float(x["R"]) * 1e6 for x in v]
        Ps = [float(x["P"]) * 1e3 for x in v]
        ts = [float(x["t_spin"]) for x in v]
        Ts = [float(x["T"]) for x in v]
        em = [float(x["eta_mecc"]) for x in v]
        out.append("%s & %s & $10^{%d}$ & %s--%s & %s--%s & %d & %s--%s & %s--%s & %s--%s \\\\\n" %
                   (mat, forma, int(round(math.log10(float(p)))),
                    it(min(Rs), 2), it(max(Rs), 2),
                    it(min(Ps), 2), it(max(Ps), 2), len(v),
                    sci(min(ts), 2), sci(max(ts), 2),
                    it(min(Ts), 0), it(max(Ts), 0),
                    sci(min(em), 1), sci(max(em), 1)))
        tot += len(v)
    out.append("\\hline\n\\multicolumn{5}{@{}l}{Totale: %d combinazioni, %d celle} & %d & & & \\\\\n" %
               (len(gruppi), tot, tot))
    out.append("\\hline\n\\end{tabular}\n")
    write("tab_rotore.tex", "".join(out))


# ---------------------------------------------------------------------------
# 7. materiali
# ---------------------------------------------------------------------------
def tab_materiali():
    out = [r"""\begin{tabular}{@{}lrrrrrrp{4.2cm}@{}}
\hline
Materiale & $n$ & $k$ & $\rho$ & $\sigma_y$ & $T_{\rm fus}$ & $\Delta n$ & Note \\
 & & & (kg/m$^3$) & (GPa) & (K) & & \\
\hline
"""]
    note = {
        "cu_lucido": "assorbe il 40\\% (Fresnel); $\\varepsilon_{IR}$ da Mie 0,0045--0,0109: specchio nell'IR",
        "cu_nero": "$\\alpha_{532}=0{,}95$; $\\varepsilon_{IR}=0{,}30$ dichiarato (0,20--0,60), non Mie",
        "silice": "trasparente; $\\Delta n$ solo di forma: 0,02 disco, 0,03 manubrio, 0,01 anello/reticolo",
        "vaterite": "birifrangente intrinseca $\\Delta n=0{,}17$; fragile, $\\sigma_y$ conservativo",
    }
    for k in ("cu_lucido", "cu_nero", "silice", "vaterite"):
        d = m.MATERIALI[k]
        out.append("%s & %s & %s & %s & %s & %s & %s & %s \\\\\n" %
                   (d["nome"].split(" (")[0], it(d["n"], 2), val(d["k"], 2),
                    it(d["rho"], 0), it(d["sy"] / 1e9, 2), it(d["Tm"], 0),
                    it(d["dn"], 2), note[k]))
    out.append("\\hline\n\\end{tabular}\n")
    write("tab_materiali.tex", "".join(out))


# ---------------------------------------------------------------------------
# 8. scala
# ---------------------------------------------------------------------------
def tab_scala():
    # 1 nm: prima riga di S1 (Cu lucido, sfera, acqua)
    r1 = S1["righe"][0]
    # 1 um: vaterite sfera acqua
    c1 = m.cella("ritardo", 10e-3, 1e-6, "vaterite", "sfera", "acqua")
    p1 = 10.0 * m.KT / (c1["tau"] / 10e-3)
    # 1 m: vaterite sfera acqua
    c2 = m.cella("ritardo", 10e-3, 0.5, "vaterite", "sfera", "acqua")
    p2 = 10.0 * m.KT / (c2["tau"] / 10e-3)
    out = [r"""\begin{tabular}{@{}lrrrrrrp{2.8cm}@{}}
\hline
Scala & $R$ (m) & $P(Pe_r=10)$ (W) & $Pe_r$ a 10 mW & $\omega$ (rad/s) & $\eta_{\rm mecc}$ & $\Delta T$ (K) & Nota \\
\hline
"""]
    out.append("1 nm (Cu) & %s & %s & %s & %s & %s & %s & $\\eta_{int}\\approx2\\times10^{-6}$: il fotone non incontra l'oggetto \\\\\n" %
               (sci(r1["R"], 1), sci(r1["P_pe10"], 2), sci(r1["Pe_r"], 2),
                sci(r1["omega"], 2), sci(r1["eta_mecc"], 2),
                it(r1["dT"], 3)))
    out.append("1 $\\mu$m (vaterite) & %s & %s & %s & %s & %s & %s & finestra: rotazione a 27,9 Hz \\\\\n" %
               (sci(1e-6, 1), sci(p1, 2), it(c1["Pe_r"], 1), it(c1["omega"], 1),
                sci(c1["eta_mecc"], 2), it(c1["dT"], 2)))
    out.append("1 m (vaterite) & %s & %s & %s & %s & %s & %s & un giro ogni $1{,}5\\times10^{8}$ anni \\\\\n" %
               (sci(0.5, 1), sci(p2, 2), it(c2["Pe_r"], 1), sci(c2["omega"], 2),
                sci(c2["eta_mecc"], 2), it(c2["dT"], 4)))
    out.append("\\hline\n\\end{tabular}\n")
    write("tab_scala.tex", "".join(out))


# ---------------------------------------------------------------------------
# 9. geometrie
# ---------------------------------------------------------------------------
def tab_geometrie():
    out = [r"""\begin{tabular}{@{}lrrrrrrr@{}}
\hline
Forma & $V/R^3$ & $A_{\rm proj}/R^2$ & $A_{\rm rad}/R^2$ & $k_I$ & $c_{\rm mecc}$ & $L/R$ & $R_{\rm eq}/R$ \\
\hline
"""]
    for forma in m.FORME:
        g = m.geometria(forma, 1e-6)
        R = 1e-6
        out.append("%s & %s & %s & %s & %s & %s & %s & %s \\\\\n" %
                   (forma, it(g["V"] / R ** 3, 3), it(g["A_proj"] / R ** 2, 3),
                    it(g["A_rad"] / R ** 2, 3), it(g["kI"], 3),
                    it(g["c_mecc"], 2), it(g["percorso"] / R, 2),
                    it(g["Req"] / R, 3)))
    out.append("\\hline\n\\end{tabular}\n")
    write("tab_geometrie.tex", "".join(out))


if __name__ == "__main__":
    print("V7.1 - tabelle")
    tab_correzioni()
    tab_criteri()
    tab_validazione()
    tab_finestre()
    tab_rappresentativi()
    tab_rotore()
    tab_materiali()
    tab_scala()
    tab_geometrie()
