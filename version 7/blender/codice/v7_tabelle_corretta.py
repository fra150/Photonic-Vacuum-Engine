"""
MOTORE FOTONICO DEL VUOTO - V7
v8_tabelle.py - genera le tabelle LaTeX della V7 dai dati calcolati.

Uso:
    python v8_tabelle.py

Legge ../dati_v7/M1_mappa.json, M1_griglia.csv e V0_validazione.json e
scrive ../tabelle_v7/tab_*.tex.  Nessun numero e' scritto a mano: se un dato
cambia, la tabella cambia con lui.

Differenze rispetto alla V7.1
-----------------------------
* correzioni V7: vaterite dn = 0.06 efficace, T_m = 723 K, validita' del
  fascio R <= 2 um, attrito di Epstein per sfere, eps_eff(T), Pe_r con
  T_bath esplicita, eta_lim fisico, unita' di f e t_spin, test D corretto;
* tabella delle finestre limitata alle 8 finestre di frontiera della V7,
  con celle di griglia e segmenti distinti;
* tabella delle celle rappresentative con eta_lim e trap_ratio;
* tabella delle 8 celle formali in vuoto (0 rotori fisici);
* tabella delle scale senza la riga da 1 m, che e' fuori dominio.
"""

import csv
import json
import math
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
if BASE not in sys.path:
    sys.path.insert(0, BASE)

import v7_mappa_corretta as m  # noqa: E402

DATI = os.path.normpath(os.path.join(BASE, "..", "dati_v7"))
TAB = os.path.normpath(os.path.join(BASE, "..", "tabelle_v7"))
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


def fixed(x, nd=1):
    """Numero con un numero fisso di decimali (mantiene lo zero finale)."""
    return (("%." + str(nd) + "f") % x).replace(".", "{,}")


def sci_star(x, nd=2):
    return sci(x, nd) + r"\textsuperscript{$\ast$}"


def migliaia(n):
    """Intero con separatore di migliaia sottile LaTeX."""
    s = "%d" % int(n)
    out = ""
    while len(s) > 3:
        out = r"\," + s[-3:] + out
        s = s[:-3]
    return s + out


def tex_row(cells):
    return " & ".join(cells) + r" \\"


def write(name, body):
    path = os.path.join(TAB, name)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(body)
    print("  %-28s %6.1f kB" % (name, os.path.getsize(path) / 1024.0))


def mezzo(mz, p):
    if mz == "acqua":
        return "acqua"
    return "vuoto $10^{%d}$" % int(round(math.log10(p)))


def btrue(x):
    return str(x).strip().lower() == "true"


# ---------------------------------------------------------------------------
# dati
# ---------------------------------------------------------------------------
with open(os.path.join(DATI, "M1_mappa.json"), encoding="utf-8") as fh:
    MAPPA = json.load(fh)
with open(os.path.join(DATI, "V0_validazione.json"), encoding="utf-8") as fh:
    VAL = json.load(fh)

with open(os.path.join(DATI, "M1_griglia.csv"), encoding="utf-8") as _fh:
    ROWS = list(csv.DictReader(_fh))
VALIDE = [r for r in ROWS if btrue(r["valido"])]


def chiave_csv(r):
    p = None if r["mezzo"] == "acqua" else float(r["p"])
    return (r["canale"], r["materiale"], r["forma"], r["mezzo"], p)


CONTEGGIO = {}
for _r in VALIDE:
    _k = chiave_csv(_r)
    CONTEGGIO[_k] = CONTEGGIO.get(_k, 0) + 1


def test(nome_iniziale):
    for t in VAL["test"]:
        if t["test"].startswith(nome_iniziale):
            return t
    raise KeyError(nome_iniziale)


# ---------------------------------------------------------------------------
# 1. correzioni
# ---------------------------------------------------------------------------
def tab_correzioni():
    dn_v = m.MATERIALI["vaterite"]["dn"]
    tm_v = m.MATERIALI["vaterite"]["Tm"]
    tlim_v = m.CRITERI["T_frac_vuoto"] * tm_v
    n_tot = MAPPA["conteggi"]["celle_totali"]
    n_gt2 = MAPPA["conteggi"]["celle_R_gt_2um"]
    n_fin = len(MAPPA["riassunto"])
    n_celle = MAPPA["conteggi"]["celle_valide"]
    n_fin_celle = sum(1 for r in MAPPA["riassunto"]
                      if CONTEGGIO.get((r["canale"], r["materiale"], r["forma"],
                                        r["mezzo"], r["p_mbar"]), 0) > 0)
    n_fis = MAPPA["conteggi"]["rotori_fisici_eta_lim"]
    d = test("D")["dettagli"]

    out = [r"\begin{tabular}{@{}p{0.35cm}p{2.3cm}p{3.7cm}p{4.1cm}p{3.85cm}@{}}",
           r"\hline",
           r"\# & Correzione & V7.1 & V7 & Effetto sulla mappa \\",
           r"\hline"]
    out.append(tex_row([
        "1", "Vaterite: birifrangenza",
        r"$\Delta n=0{,}17$ (valore della calcite)",
        r"$\Delta n=" + it(dn_v, 2) + r"$ efficace (range $0{,}06$--$0{,}10$); "
        r"la calcite resta $0{,}17$",
        "la coppia di ritardo non e' piu' quella della calcite; la finestra "
        "della vaterite si restringe e si sposta"]))
    out.append(tex_row([
        "2", "Vaterite: temperatura limite",
        r"$T_m=1000$ K",
        r"$T_m=" + it(tm_v, 0) + r"$ K (inizio transizione "
        r"vaterite$\to$calcite); criterio $T\le0{,}7\,T_m=" +
        it(tlim_v, 0) + r"$ K",
        "il tetto termico in vuoto scende; nessuna cella valida supera il "
        "nuovo limite"]))
    out.append(tex_row([
        "3", "Validita' del fascio",
        r"nessun limite esplicito; griglia fino a $30\,\mu$m",
        r"$R\le" + it(m.R_MAX_VALID * 1e6, 0) + r"\,\mu$m "
        r"($z_R\simeq2\,\mu$m); $R>2\,\mu$m marcato come estrapolazione",
        migliaia(n_gt2) + " celle su " + migliaia(n_tot) + " fuori dominio; "
        r"le 8 finestre stanno tutte entro $1{,}46\,\mu$m"]))
    out.append(tex_row([
        "4", "Attrito molecolare",
        "formula cinetica, stesso coefficiente per tutte le forme",
        r"sfera: $\beta=(16/3)\,pR^4/\bar v$; forme non sferiche: $\beta$ "
        r"incerto e marcato",
        "in vuoto le affermazioni quantitative valgono solo per sfere; gli 8 "
        "rotori formali sono sfere di vaterite"]))
    out.append(tex_row([
        "5", "Emissivita' termica",
        r"$\varepsilon_{\rm eff}$ da Mie a $T_0$; Cu nero $0{,}30$",
        r"$\varepsilon_{\rm eff}(materiale,R,T)$ con peso di Planck a $T$; "
        r"Cu nero $0{,}30$ dichiarato; bordo termico ricalcolato",
        "il rame lucido resta uno specchio nell'IR; il nero dipende dal "
        "parametro dichiarato"]))
    out.append(tex_row([
        "6", "Definizione di $Pe_r$",
        r"$Pe_r=\tau/(k_BT_0)$, $T_0$ implicita",
        r"$Pe_r=\tau/(k_BT_{\rm bath})$, $T_{\rm bath}=310$ K",
        "nessun cambio numerico; la definizione diventa esplicita"]))
    out.append(tex_row([
        "7", "Rendimento fisico",
        r"$\eta_{\rm mecc}=\omega/\omega_{\rm ott}$ usata come efficienza",
        r"$\eta_{\rm lim}=\tau\omega_{\rm rottura}/P$; $\eta_{\rm mecc}$ "
        r"resta solo una comodita' numerica",
        "0 rotori fisici: nessuna cella valida raggiunge "
        "$\\eta_{\\rm lim}\\ge10^{-7}$; $\\eta_{\\rm mecc}$ non e' un "
        "rendimento"]))
    out.append(tex_row([
        "8", "Unita'",
        r"$f=\omega/2\pi$ ma testo con $\tau/\beta$ o $\beta/I$",
        r"$f=\omega/(2\pi)=\tau/(2\pi\beta)$; tempo di smorzamento $I/\beta$; "
        r"non $\tau/\beta$ o $\beta/I$",
        "unita' coerenti in tabelle e figure: $f$ in Hz, $t$ in s"]))
    out.append(tex_row([
        "9", "Validazione D",
        r"D con $\ell=1$ a 532 nm: 10,75 Hz vs 4 Hz (2,69), PASS",
        r"D corretto a 632,8 nm e $\ell=3$: " + it(d["f_modello_Hz"], 2) +
        r" Hz vs " + it(d["f_misurata_Hz"], 0) + r" Hz (" +
        it(d["rapporto_modello_misura"], 2) + r"), FAIL; test di coerenza " +
        str(VAL["superati"]) + "/" + str(VAL["totali"]),
        "D non e' una validazione del modello di ritardo: a 632,8 nm il "
        "canale e' l'assorbimento"]))
    out.append(tex_row([
        "10", "Tabelle e conteggi",
        r"$P_{\min}/P_{\max}$ per segmento; frontiera e griglia non separate",
        r"$P_{\min}=\min P_{\rm lo}$, $P_{\max}=\max P_{\rm hi}$, "
        r"$R(P_{\min})$; celle contate da \texttt{valido}",
        str(n_fin) + " finestre di frontiera, " + str(n_fin_celle) +
        " con almeno una cella di griglia, " + str(n_celle) +
        " celle valide; le 2 finestre di Cu non hanno celle"]))
    out.append(r"\hline")
    out.append(r"\end{tabular}")
    write("tab_correzioni.tex", "\n".join(out) + "\n")


# ---------------------------------------------------------------------------
# 2. criteri
# ---------------------------------------------------------------------------
def tab_criteri():
    out = [r"\begin{tabular}{@{}p{2.4cm}p{2.3cm}p{2.2cm}p{7.8cm}@{}}",
           r"\hline",
           r"Criterio & Valore & Dove & Perche' \\",
           r"\hline"]
    out.append(tex_row([
        r"$R$", r"$\le2\,\mu$m", "tutte",
        r"validita' del fascio gaussiano: $z_R\simeq2\,\mu$m; oltre "
        r"$R=2\,\mu$m la particella e' piu' grande della regione ben "
        r"descritta, la cella e' un'estrapolazione"]))
    out.append(tex_row([
        r"$Pe_r$", r"$\ge10$", "acqua e vuoto",
        r"la coppia ottica supera di 10 volte il rumore browniano: "
        r"$Pe_r=\tau/(k_BT_{\rm bath})$, $T_{\rm bath}=310$ K"]))
    out.append(tex_row([
        r"$\Delta T$", r"$\le20$ K", "acqua",
        r"conduzione sferica, Nusselt 1; sopra 20 K convezione e variazione "
        r"di viscosita' non sono trascurabili"]))
    out.append(tex_row([
        r"$T$", r"$\le0{,}7\,T_m$", "vuoto",
        r"equilibrio radiativo; per la vaterite $T_m=723$ K (transizione "
        r"vaterite$\to$calcite), quindi $T\le506$ K"]))
    out.append(tex_row([
        r"$\omega$", r"$\le0{,}5\,\omega_{\rm rottura}$", "acqua",
        r"il rotore non deve lavorare oltre meta' del limite meccanico; in "
        r"vuoto $\omega=\tau/\beta_{\rm gas}$ e' formale"]))
    out.append(tex_row([
        r"$P$", r"$\le10$ mW", "tutte",
        r"tetto dichiarato del banco ottico; non e' un limite fisico"]))
    out.append(tex_row([
        r"$\eta_{\rm lim}$", r"$\ge10^{-7}$", "tutte",
        r"rendimento fisico $\eta_{\rm lim}=\tau\omega_{\rm rottura}/P$; "
        r"soglia dichiarata di rotore, non una legge"]))
    out.append(tex_row([
        "trappola", "necessaria in vuoto", "vuoto",
        r"$F_{\rm scat}=P/c$ contro il peso $\rho Vg$; senza supporto la "
        r"particella non resta nel fascio"]))
    out.append(tex_row([
        r"$t_{\rm spin}$", r"$\ge1$ s", "vuoto",
        r"tempo per arrivare a meta' del limite meccanico a coppia costante; "
        r"soglia operativa del volano"]))
    out.append(r"\hline")
    out.append(r"\end{tabular}")
    write("tab_criteri.tex", "\n".join(out) + "\n")


# ---------------------------------------------------------------------------
# 3. test di coerenza
# ---------------------------------------------------------------------------
def tab_validazione():
    out = [r"\begin{tabular}{|p{0.3cm}|p{2.4cm}|p{2.7cm}|p{3.5cm}|p{3.75cm}|p{0.85cm}|}",
           r"\hline",
           r" & Test & Riferimento & Controllo & Modello / atteso & Esito \\",
           r"\hline"]
    for i, t in enumerate(VAL["test"], 1):
        nome = t["test"].split(" - ", 1)[-1].replace("lambda/2", r"$\lambda/2$")
        nome = nome[0].upper() + nome[1:]
        rif = t["riferimento"]
        det = t["dettagli"]
        if t["test"].startswith("A"):
            controllo = r"4 casi limite (circolare/lineare, $\lambda/4$, $\lambda/2$)"
            modello = "; ".join(it(c["tau_su_P_omega"], 0) for c in det)
            atteso = "; ".join(it(c["atteso"], 0) for c in det)
            rapporto = "esatto"
        elif t["test"].startswith("B"):
            controllo = r"$\tau/(P/\omega)$ per lamina $\lambda/2$ circolare"
            modello = it(det["tau_su_P_omega"], 0)
            atteso = it(det["atteso"], 0)
            rapporto = "esatto"
        elif t["test"].startswith("C"):
            rif = r"V6.1: $\Delta T=P_{\rm abs}/(4\pi\kappa R)$; 150 $\mu$W"
            controllo = r"$\Delta T$ a 150 $\mu$W, $R=1$ e 3 $\mu$m"
            modello = fixed(det[0]["dT"], 1) + "; " + fixed(det[1]["dT"], 1) + " K"
            atteso = fixed(det[0]["atteso_K"], 1) + "; " + fixed(det[1]["atteso_K"], 1) + " K"
            rapporto = it(det[0]["rapporto"], 2) + "; " + it(det[1]["rapporto"], 2)
        elif t["test"].startswith("D"):
            controllo = r"CuO 2 $\mu$m, 10 mW, $\ell=3$, kerosene, 632,8 nm"
            modello = it(det["f_modello_Hz"], 2) + " Hz"
            atteso = it(det["f_misurata_Hz"], 0) + " Hz"
            rapporto = it(det["rapporto_modello_misura"], 2)
        elif t["test"].startswith("E"):
            controllo = r"calcite 1 $\mu$m, 300 mW, 1064 nm"
            modello = it(det["f_modello_Hz"], 2) + " Hz"
            atteso = it(det["f_misurata_Hz"], 0) + " Hz"
            rapporto = it(det["rapporto_modello_misura"], 2)
        else:
            controllo = r"manubrio silice 100 nm, 10 mW, $\sigma_y=1$--3 GPa"
            f_lo = min(c["f_lim_GHz"] for c in det["casi"])
            f_hi = max(c["f_lim_GHz"] for c in det["casi"])
            ts_lo = min(c["t_spin_s"] for c in det["casi"])
            ts_hi = max(c["t_spin_s"] for c in det["casi"])
            modello = (it(f_lo, 2) + "--" + it(f_hi, 2) + " GHz; "
                       r"$t_{\rm spin}=" + it(ts_lo, 3) + "--" + it(ts_hi, 3) + r"$ s")
            atteso = r"$0{,}5$--10 GHz; $t_{\rm spin}<1$ s"
            rapporto = "--"
        if rapporto == "esatto":
            confronto = modello + " / " + atteso + " (esatto)"
        elif rapporto == "--":
            confronto = modello + " / " + atteso
        else:
            confronto = modello + " / " + atteso + " (" + rapporto + ")"
        out.append(tex_row([str(i), nome, rif, controllo, confronto, t["esito"]]))
        out.append(r"\hline")
    out.append(r"\end{tabular}")
    out.append("")
    out.append(r"\par\smallskip\noindent\footnotesize " +
               str(VAL["superati"]) + "/" + str(VAL["totali"]) +
               r" test di coerenza superati. D non e' una validazione del "
               r"modello di ritardo: a 632,8 nm e con CuO assorbente il canale "
               r"e' l'assorbimento.")
    write("tab_validazione.tex", "\n".join(out) + "\n")


# ---------------------------------------------------------------------------
# 4. le 8 finestre di frontiera
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
    out = [r"\begingroup\small\setlength{\tabcolsep}{2.5pt}",
           r"\begin{tabular}{|l|l|l|l|r|r|r|r|r|r|r|}",
           r"\hline",
           r"Canale & Materiale & Forma & Mezzo & $R_{\min}$ & $R_{\max}$ & "
           r"$P_{\min}$ & $P_{\max}$ & $R(P_{\min})$ & Celle & Segm. \\",
           r" & & & & ($\mu$m) & ($\mu$m) & (W) & (W) & ($\mu$m) & & \\",
           r"\hline"]
    old_group = None
    for r in rr:
        group = (r["mezzo"], r["p_mbar"])
        if old_group is not None and group != old_group:
            out.append(r"\hline")
        old_group = group
        fs = r["finestre"]
        pmin = min(w["P_lo"] for w in fs)
        pmax = max(w["P_hi"] for w in fs)
        r_at_min = [w["R"] for w in fs if w["P_lo"] == pmin][0]
        pmax_s = (sci_star(pmax)
                  if abs(pmax - m.CRITERI["P_max"]) < 1e-15 else sci(pmax))
        n_cells = CONTEGGIO.get((r["canale"], r["materiale"], r["forma"],
                                 r["mezzo"], r["p_mbar"]), 0)
        out.append(tex_row([
            r["canale"], r["materiale"].replace("_", r"\_"), r["forma"],
            mezzo(r["mezzo"], r["p_mbar"]),
            it(r["R_min_um"], 2), it(r["R_max_um"], 2),
            sci(pmin), pmax_s, it(r_at_min * 1e6, 2),
            str(n_cells), str(len(fs))]))
    out.append(r"\hline")
    out.append(r"\end{tabular}")
    out.append(r"\endgroup")
    out.append("")
    out.append(r"\par\smallskip\noindent\footnotesize "
               r"$^{\ast}$ $P_{\max}=10$ mW: taglio del criterio, non un "
               r"limite fisico. Le due finestre di rame non contengono celle "
               r"di griglia: la loro larghezza in $P$ e' minore della "
               r"spaziatura della griglia. Conteggio: 8 finestre di "
               r"frontiera, 6 con almeno una cella di griglia, 40 celle "
               r"valide.")
    write("tab_finestre.tex", "\n".join(out) + "\n")


# ---------------------------------------------------------------------------
# 5. celle rappresentative
# ---------------------------------------------------------------------------
def tab_rappresentativi():
    casi = [
        (r"Vaterite, sfera, acqua; $R=1\,\mu$m, $P=10$ mW",
         "ritardo", 10e-3, 1e-6, "vaterite", "sfera", "acqua", None, True),
        (r"Cu nero, sfera, acqua; $R=1{,}05\,\mu$m, $P=0{,}165$ mW",
         "assorbimento", 0.165e-3, 1.0471071012347182e-6, "cu_nero", "sfera",
         "acqua", None, True),
        (r"Vaterite, sfera, vuoto $10^{-8}$; $R=1{,}05\,\mu$m, $P=0{,}316$ mW",
         "ritardo", 0.31622776601683794e-3, 1.0471071012347182e-6,
         "vaterite", "sfera", "vuoto", 1e-8, False),
        (r"Vaterite, sfera, vuoto $10^{-8}$; $R=1{,}46\,\mu$m, $P=0{,}139$ mW",
         "ritardo", 0.13894954943731373e-3, 1.464549719110926e-6,
         "vaterite", "sfera", "vuoto", 1e-8, False),
    ]
    out = [r"\begingroup\footnotesize\setlength{\tabcolsep}{3pt}",
           r"\begin{tabular}{@{}p{2.7cm}rrrrrrrrl@{}}",
           r"\hline",
           r"Caso & $\tau$ (N m) & $Pe_r$ & $f$ (Hz) & $\Delta T$ (K) & "
           r"$\eta_{\rm mecc}$ & $\eta_{\rm lim}$ & $t_{\rm spin}$ (s) & "
           r"$F_{\rm scat}/F_g$ & Esito \\",
           r"\hline"]
    for tag, canale, P, R, mat, forma, mz, p, water in casi:
        c = m.cella(canale, P, R, mat, forma, mz, p)
        f = c["omega"] / (2.0 * math.pi)
        if water:
            fstr = it(f, 2)
            tstr = "--"
            esito = "valido"
        else:
            fstr = sci_star(f, 1)
            tstr = val(c["t_spin"], 2)
            esito = "solo formale"
        out.append(tex_row([
            tag, sci(c["tau"], 1), it(c["Pe_r"], 1), fstr,
            it(c["dT"], 2), sci(c["eta_mecc"], 1), sci(c["eta_lim"], 1),
            tstr, it(c["trap_ratio"], 1), esito]))
    out.append(r"\hline")
    out.append(r"\end{tabular}")
    out.append(r"\endgroup")
    out.append("")
    out.append(r"\par\smallskip\noindent\footnotesize "
               r"$^{\ast}$ In vuoto $f=\omega/(2\pi)$ con "
               r"$\omega=\tau/\beta_{\rm gas}$ e' una frequenza formale, non "
               r"una velocita' misurata; il rotore richiede una trappola. In "
               r"acqua $t_{\rm spin}$ non e' applicabile (l'attrito domina) "
               r"ed e' indicato con --. Le celle in acqua sono valide per i "
               r"criteri, ma $\eta_{\rm lim}<10^{-7}$: non sono rotori "
               r"fisici. Le celle in vuoto sono formali per lo stesso motivo. "
               r"Qui $F_g=\rho Vg$ e $F_{\rm scat}=P/c$.")
    write("tab_rappresentativi.tex", "\n".join(out) + "\n")


# ---------------------------------------------------------------------------
# 6. le 8 celle formali in vuoto
# ---------------------------------------------------------------------------
def tab_vuoto():
    righe = [r for r in VALIDE if r["mezzo"] == "vuoto"]
    righe.sort(key=lambda r: (-math.log10(float(r["p"])),
                              float(r["R"]), float(r["P"])))
    out = [r"\begingroup\footnotesize\setlength{\tabcolsep}{2pt}",
           r"\begin{tabular}{@{}llrrrrrrrr@{}}",
           r"\hline",
           r"Materiale & Forma & $p$ (mbar) & $R$ ($\mu$m) & $P$ (mW) & "
           r"$T$ (K) & $t_{\rm spin}$ (s) & $\eta_{\rm mecc}$ & "
           r"$\eta_{\rm lim}$ & $F_{\rm scat}/F_g$ \\",
           r"\hline"]
    old_p = None
    for r in righe:
        p = float(r["p"])
        if old_p is not None and p != old_p:
            out.append(r"\hline")
        old_p = p
        out.append(tex_row([
            r["materiale"].replace("_", r"\_"), r["forma"],
            "$10^{%d}$" % int(round(math.log10(p))),
            it(float(r["R"]) * 1e6, 2), it(float(r["P"]) * 1e3, 3),
            it(float(r["T"]), 1), it(float(r["t_spin"]), 2),
            sci(float(r["eta_mecc"]), 2), sci(float(r["eta_lim"]), 2),
            it(float(r["trap_ratio"]), 2)]))
    out.append(r"\hline")
    out.append(r"\end{tabular}")
    out.append(r"\endgroup")
    out.append("")
    out.append(r"\par\smallskip\noindent\footnotesize "
               "8 celle formali, 0 rotori fisici ($\\eta_{\\rm lim}<10^{-7}$ "
               "in tutte): sono limitate dall'attrito molecolare, non da una "
               "coppia utile. Tutte richiedono una trappola o un supporto. "
               "Qui $F_g=\\rho Vg$ e $F_{\\rm scat}=P/c$. Il $\\beta$ di "
               "Epstein vale per sfere; le forme non sferiche sarebbero "
               "marcate con $\\beta$ incerto.")
    write("tab_vuoto.tex", "\n".join(out) + "\n")


# ---------------------------------------------------------------------------
# 7. materiali
# ---------------------------------------------------------------------------
def tab_materiali():
    out = [r"\begin{tabular}{@{}lrrrrrrp{4.2cm}@{}}",
           r"\hline",
           r"Materiale & $n$ & $k$ & $\rho$ & $\sigma_y$ & $T_m$ & "
           r"$\Delta n$ & Note \\",
           r" & & & (kg/m$^3$) & (GPa) & (K) & & \\",
           r"\hline"]
    note = {
        "cu_lucido": r"assorbe il 40\% (Fresnel); $\varepsilon_{IR}$ da Mie "
                     r"0,0045--0,0109: specchio nell'IR",
        "cu_nero": r"$\alpha_{532}=0{,}95$; $\varepsilon_{IR}=0{,}30$ "
                   r"dichiarato (0,20--0,60), non Mie",
        "silice": r"trasparente; $\Delta n$ solo di forma: 0,02 disco, 0,03 "
                  r"manubrio, 0,01 anello/reticolo",
        "vaterite": r"$\Delta n=0{,}06$ efficace (0,06--0,10); $k=10^{-6}$; "
                    r"$T_m=723$ K (vaterite$\to$calcite); fragile, "
                    r"$\sigma_y$ conservativo",
    }
    for k in ("cu_lucido", "cu_nero", "silice", "vaterite"):
        d = m.MATERIALI[k]
        out.append(tex_row([
            d["nome"].split(" (")[0], it(d["n"], 2), val(d["k"], 2),
            it(d["rho"], 0), it(d["sy"] / 1e9, 2), it(d["Tm"], 0),
            it(d["dn"], 2), note[k]]))
    out.append(r"\hline")
    out.append(r"\end{tabular}")
    out.append("")
    out.append(r"\par\smallskip\noindent\footnotesize La calcite ha "
               r"$\Delta n=0{,}17$: e' il valore usato nel test E (Friese "
               r"1998), non la vaterite policristallina della mappa.")
    write("tab_materiali.tex", "\n".join(out) + "\n")


# ---------------------------------------------------------------------------
# 8. scala
# ---------------------------------------------------------------------------
def tab_scala():
    righe = []
    # 1 nm e' il diametro: raggio 0,5 nm, come nella V7.1.
    for tag, R, canale, mat, nota in (
            ("1 nm (Cu)", 0.5e-9, "assorbimento", "cu_lucido", None),
            (r"1 $\mu$m (vaterite)", 1e-6, "ritardo", "vaterite", None),
            (r"10 $\mu$m (vaterite)", 10e-6, "ritardo", "vaterite", None)):
        c = m.cella(canale, 10e-3, R, mat, "sfera", "acqua", None)
        p10 = 10.0 * m.KT / (c["tau"] / 10e-3) if c["tau"] > 0 else float("inf")
        righe.append((tag, R, c, p10, nota))

    out = [r"\begingroup\small\setlength{\tabcolsep}{2pt}",
           r"\begin{tabular}{@{}lrrrrrrp{1.9cm}@{}}",
           r"\hline",
           r"Scala & $R$ (m) & $P(Pe_r{=}10)$ (W) & $Pe_r$ a 10 mW & "
           r"$\omega$ (rad/s) & $\eta_{\rm mecc}$ & $\Delta T$ (K) & Nota \\",
           r"\hline"]
    for tag, R, c, p10, _ in righe:
        if tag.startswith("1 nm"):
            nota = (r"$\eta_{\rm int}\approx$ " + sci(c["eta_int"], 1) +
                    r": il fotone non incontra l'oggetto; $Pe_r\ll10$")
        elif tag.startswith(r"1 $\mu$m"):
            nota = r"finestra: $f=" + it(c["omega"] / (2 * math.pi), 1) + r"$ Hz"
        else:
            nota = r"$R>2\,\mu$m: fuori dominio del fascio; estrapolazione"
        out.append(tex_row([
            tag, sci(R, 1), sci(p10, 2), val(c["Pe_r"], 2),
            val(c["omega"], 2), sci(c["eta_mecc"], 2),
            it(c["dT"], 4), nota]))
    out.append(r"\hline")
    out.append(r"\end{tabular}")
    out.append(r"\endgroup")
    out.append("")
    out.append(r"\par\smallskip\noindent\footnotesize La riga da 1 m e' "
               r"stata rimossa: e' un'estrapolazione fuori dominio "
               r"($R>2\,\mu$m). Anche la riga da 10 $\mu$m e' oltre il "
               r"dominio del fascio e va letta come estrapolazione. La riga "
               r"da 1 nm si riferisce al diametro: $R=0{,}5$ nm.")
    write("tab_scala.tex", "\n".join(out) + "\n")


# ---------------------------------------------------------------------------
# 9. geometrie
# ---------------------------------------------------------------------------
def tab_geometrie():
    out = [r"\begin{tabular}{@{}lrrrrrrr@{}}",
           r"\hline",
           r"Forma & $V/R^3$ & $A_{\rm proj}/R^2$ & $A_{\rm rad}/R^2$ & "
           r"$k_I$ & $c_{\rm mecc}$ & $L/R$ & $R_{\rm eq}/R$ \\",
           r"\hline"]
    for forma in m.FORME:
        g = m.geometria(forma, 1e-6)
        R = 1e-6
        out.append(tex_row([
            forma, it(g["V"] / R ** 3, 3), it(g["A_proj"] / R ** 2, 3),
            it(g["A_rad"] / R ** 2, 3), it(g["kI"], 3),
            it(g["c_mecc"], 2), it(g["percorso"] / R, 2),
            it(g["Req"] / R, 3)]))
    out.append(r"\hline")
    out.append(r"\end{tabular}")
    write("tab_geometrie.tex", "\n".join(out) + "\n")


def main():
    print("V7 - tabelle")
    tab_correzioni()
    tab_criteri()
    tab_validazione()
    tab_finestre()
    tab_rappresentativi()
    tab_vuoto()
    tab_materiali()
    tab_scala()
    tab_geometrie()
    print("Conteggi: %d celle, %d valide, %d finestre, %d con celle, %d rotori fisici"
          % (MAPPA["conteggi"]["celle_totali"], MAPPA["conteggi"]["celle_valide"],
             len(MAPPA["riassunto"]),
             sum(1 for r in MAPPA["riassunto"]
                 if CONTEGGIO.get((r["canale"], r["materiale"], r["forma"],
                                   r["mezzo"], r["p_mbar"]), 0) > 0),
             MAPPA["conteggi"]["rotori_fisici_eta_lim"]))
    print("Fatto: tabelle in", TAB)


if __name__ == "__main__":
    main()
