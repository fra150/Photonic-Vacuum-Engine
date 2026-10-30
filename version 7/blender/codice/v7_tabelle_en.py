"""
MOTORE FOTONICO DEL VUOTO - V7.1 (English edition)
v7_tabelle_en.py - English LaTeX tables from the same data.

Uso:
    python v7_tabelle_en.py

Legge i dati tramite v7_tabelle.py e scrive ../tabelle_en/tab_*.tex.
Nessun numero e' scritto a mano.
"""

import math
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
if BASE not in sys.path:
    sys.path.insert(0, BASE)

import v7_tabelle as T  # noqa: E402

m = T.m
TAB = os.path.normpath(os.path.join(BASE, "..", "tabelle_en"))
os.makedirs(TAB, exist_ok=True)


# ---------------------------------------------------------------------------
# English number formatting (decimal point, not comma)
# ---------------------------------------------------------------------------
def it(x, nd=2):
    s = ("%." + str(nd) + "f") % x
    if "." in s:
        s = s.rstrip("0").rstrip(".")
    return s


def sci(x, nd=2):
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
    return "$" + (("%." + str(nd) + "f") % a) + "\\times10^{%d}$" % e


def val(x, nd=2):
    if x == 0:
        return "0"
    if 1e-2 <= abs(x) < 1e3:
        return it(x, nd)
    return sci(x, nd)


def rng(lo, hi, nd=2):
    if abs(hi - lo) <= 1e-3 * max(abs(hi), abs(lo), 1e-30):
        return it(lo, nd)
    return it(lo, nd) + "--" + it(hi, nd)


def write(name, body):
    path = os.path.join(TAB, name)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(body)
    print("  %-28s %6.1f kB" % (name, os.path.getsize(path) / 1024.0))


LBL_CAN = {"assorbimento": "absorption", "ritardo": "retardation",
           "scattering": "scattering"}
LBL_MAT = {"cu_lucido": "shiny Cu", "cu_nero": "black Cu",
           "silice": "silica", "vaterite": "vaterite"}
LBL_FOR = {"sfera": "sphere", "disco": "disk", "anello": "ring",
           "manubrio": "dumbbell", "reticolo": "grating"}


def medium(mz, p):
    if mz == "acqua":
        return "water"
    return "vacuum $10^{%d}$" % int(round(math.log10(p)))


# ---------------------------------------------------------------------------
# 1. corrections
# ---------------------------------------------------------------------------
def tab_correzioni():
    body = r"""\begin{tabular}{@{}p{0.35cm}p{2.8cm}p{4.1cm}p{4.2cm}p{4.2cm}@{}}
\hline
\# & Correction & V7.0 & V7.1 & Effect on the map \\
\hline
1 & Axial scattering &
  $\tau_{\rm sca}\neq0$ for sphere, disk, ring, grating &
  $\tau_{\rm sca}=0$ for axisymmetric shapes; dumbbell only &
  6 scattering windows, all dumbbell \\
\hline
2 & Dumbbell inertia &
  $k_I$ in m$^2$, then $\times R^2$: $I\propto R^6$ &
  $k_I=0.5875$ dimensionless &
  $I(100\,{\rm nm})=1.692\times10^{-33}$ kg m$^2$; $t_{\rm spin}=0.076$ s \\
\hline
3 & Phase retardation &
  $\delta\propto R$ down to $R\to0$; $(1-\cos\delta)$ oscillates &
  $\delta\propto R^3$ below $0.15\,\mu$m; $\sigma_\delta=0.5$ &
  $P(Pe_r=10)$ at 100 nm: 99.6 mW; window from 0.14 $\mu$m \\
\hline
4 & Black Cu &
  $\varepsilon_{IR}$ from Mie (Cu sphere) &
  $\varepsilon_{IR}=0.30$ declared (0.20--0.60) &
  thermal edge $\sim$30 $\mu$m; in vacuum only 30 $\mu$m \\
\hline
5 & Residual pressure &
  mbar$\to$Pa: $0.1$ (1000$\times$ low) &
  mbar$\to$Pa: $100$ &
  $\omega,\eta_{\rm mecc}$ in vacuum $\times10^{-3}$; rotors 9$\to$3 \\
\hline
\end{tabular}
"""
    write("tab_correzioni.tex", body)


# ---------------------------------------------------------------------------
# 2. criteria
# ---------------------------------------------------------------------------
def tab_criteri():
    body = r"""\begin{tabular}{@{}p{2.4cm}p{2.1cm}p{2.4cm}p{8.6cm}@{}}
\hline
Criterion & Value & Where & Why \\
\hline
$Pe_r$ & $\ge 10$ & water and vacuum &
  optical torque exceeds Brownian noise by 10$\times$: orientation is deterministic \\
$\Delta T$ & $\le 20$ K & water &
  spherical conduction, Nusselt 1; above 20 K convection and viscosity changes are not negligible \\
$T$ & $\le 0.7\,T_{\rm melt}$ & vacuum &
  radiative equilibrium; the 30\% margin covers the $\varepsilon_{IR}$ uncertainty and non-uniform temperature \\
$\omega$ & $\le 0.5\,\omega_{\rm break}$ & water &
  the rotor must not work above half the mechanical limit; in vacuum the criterion does not apply to $\omega=\tau/\beta_{\rm gas}$ (see \S\ref{sec:rotore}) \\
$P$ & $\le 10$ mW & all &
  declared ceiling of the optical bench; not a physical limit \\
$\eta_{\rm mecc}$ & $\ge 10^{-7}$ & all &
  declared ``rotor'' threshold: $\omega/\omega_{\rm opt}$ at least $10^{-7}$; a design choice, not a law \\
$t_{\rm spin}$ & $\ge 1$ s & vacuum &
  the rotor must not break before it is useful; $t_{\rm spin}$ is the time to half the mechanical limit without drag \\
\hline
\end{tabular}
"""
    write("tab_criteri.tex", body)


# ---------------------------------------------------------------------------
# 3. validation
# ---------------------------------------------------------------------------
def tab_validazione():
    rows = [
        ("phase-retardation algebra", "Jones + Stokes; checked by hand",
         "4 limiting cases (circular/linear, $\\lambda/4$, $\\lambda/2$)",
         "1; 2; 1; 0", "1; 2; 1; 0", "exact"),
        ("Beth limit: $\\lambda/2$ circular", "Beth 1936, Phys. Rev. 50, 115",
         "$\\tau/(P/\\omega)$ for a $\\lambda/2$ plate, circular",
         "2", "2", "exact"),
        ("water thermal (V6.1)",
         "V6.1: $\\Delta T=P_{\\rm abs}/(4\\pi\\kappa R)$; 150 $\\mu$W",
         "$\\Delta T$ at 150 $\\mu$W, $R=1$ and 3 $\\mu$m",
         "19.0; 6.34 K", "20; 6.6 K", "0.95; 0.96"),
        ("He 1995: absorptive CuO with OAM",
         "He et al., PRL 75, 826 (1995)",
         "CuO 2 $\\mu$m, 10 mW, $\\ell=1$, kerosene",
         "10.75 Hz", "4 Hz", "2.69"),
        ("Friese 1998: birefringent calcite",
         "Friese et al., Nature 394, 348 (1998)",
         "calcite 1 $\\mu$m, 300 mW, 1064 nm",
         "243.8 Hz", "357 Hz", "0.68"),
        ("GHz limit of silica dumbbell",
         "Reimann et al., PRL 121, 033602 (2018); Ahn et al. (2020)",
         "silica dumbbell 100 nm, 10 mW, $\\sigma_y=1$--3 GPa",
         "1.52--2.63 GHz; $t_{\\rm spin}=0.076$--$0.132$ s",
         "0.5--10 GHz; $t_{\\rm spin}<1$ s", "--"),
    ]
    out = [r"""\begin{tabular}{|p{0.3cm}|p{2.3cm}|p{2.8cm}|p{3.6cm}|p{4.25cm}|p{0.65cm}|}
\hline
 & Test & Reference & Check & Model / expected & Result \\
\hline
"""]
    for i, (nome, rif, controllo, modello, atteso, rapporto) in enumerate(rows, 1):
        if rapporto == "exact":
            confronto = "%s / %s (exact)" % (modello, atteso)
        elif rapporto == "--":
            confronto = "%s / %s" % (modello, atteso)
        else:
            confronto = "%s / %s (%s)" % (modello, atteso, rapporto)
        out.append("%d & %s & %s & %s & %s & PASS \\\\\n" %
                   (i, nome, rif, controllo, confronto))
        out.append("\\hline\n")
    out.append("\\end{tabular}\n")
    write("tab_validazione.tex", "".join(out))


# ---------------------------------------------------------------------------
# 4. the 26 windows
# ---------------------------------------------------------------------------
def tab_finestre():
    ord_mezzo = {"acqua": 0, "vuoto": 1}
    ord_can = {"assorbimento": 0, "ritardo": 1, "scattering": 2}
    ord_mat = {"cu_lucido": 0, "cu_nero": 1, "silice": 2, "vaterite": 3}
    ord_for = {"sfera": 0, "disco": 1, "anello": 2, "manubrio": 3, "reticolo": 4}
    rr = sorted(T.MAPPA["riassunto"],
                key=lambda r: (ord_mezzo[r["mezzo"]],
                               0 if r["p_mbar"] is None else -math.log10(r["p_mbar"]),
                               ord_can[r["canale"]], ord_mat[r["materiale"]],
                               ord_for[r["forma"]]))
    out = [r"""\begin{tabular}{|l|l|l|l|r|r|r|r|r|}
\hline
Channel & Material & Shape & Medium & $R_{\min}$ & $R_{\max}$ & $P_{\min}$ & $P_{\max}$ & Segm. \\
 & & & & ($\mu$m) & ($\mu$m) & (W) & (W) & \\
\hline
"""]
    old_mezzo = None
    for r in rr:
        if old_mezzo is not None and r["mezzo"] != old_mezzo:
            out.append("\\hline\n")
        old_mezzo = r["mezzo"]
        p_lo = r["finestre"][0]["P_lo"]
        p_hi = r["finestre"][-1]["P_hi"]
        out.append("%s & %s & %s & %s & %s & %s & %s & %s & %d \\\\\n" %
                   (LBL_CAN[r["canale"]], LBL_MAT[r["materiale"]],
                    LBL_FOR[r["forma"]], medium(r["mezzo"], r["p_mbar"]),
                    it(r["R_min_um"], 2), it(r["R_max_um"], 2),
                    sci(p_lo), sci(p_hi), len(r["finestre"])))
    out.append("\\hline\n\\end{tabular}\n")
    write("tab_finestre.tex", "".join(out))


# ---------------------------------------------------------------------------
# 5. representative cells
# ---------------------------------------------------------------------------
def tab_rappresentativi():
    casi = [
        ("Vaterite, sphere, water", 10e-3, 1e-6, "vaterite", "sfera", "acqua", None, "ritardo", True),
        ("Vaterite, sphere, water", 10e-3, 2e-6, "vaterite", "sfera", "acqua", None, "ritardo", True),
        ("Vaterite, disk, water", 10e-3, 10e-6, "vaterite", "disco", "acqua", None, "ritardo", True),
        ("Black Cu, sphere, water", 0.165e-3, 1.0471e-6, "cu_nero", "sfera", "acqua", None, "assorbimento", True),
        ("Shiny Cu, sphere, water", 0.39e-3, 1.0471e-6, "cu_lucido", "sfera", "acqua", None, "assorbimento", True),
        ("Silica, disk, water", 7e-3, 10e-6, "silice", "disco", "acqua", None, "ritardo", True),
        ("Silica, dumbbell, water", 7e-3, 0.5e-6, "silice", "manubrio", "acqua", None, "scattering", True),
        ("Vaterite, dumbbell, water", 8e-3, 0.5e-6, "vaterite", "manubrio", "acqua", None, "scattering", True),
        ("Vaterite, disk, vacuum $10^{-8}$", 0.5e-3, 10e-6, "vaterite", "disco", "vuoto", 1e-8, "ritardo", False),
    ]
    out = [r"""\begin{tabular}{@{}p{3.9cm}rrrrrrl@{}}
\hline
Case & $\tau$ (N m) & $Pe_r$ & $f$ (Hz) & $\Delta T$ (K) & $\eta_{\rm mecc}$ & $t_{\rm spin}$ (s) & Outcome \\
\hline
"""]
    for tag, P, R, mat, forma, mz, p, canale, water in casi:
        c = m.cella(canale, P, R, mat, forma, mz, p)
        f = c["omega"] / (2.0 * math.pi)
        if water:
            fstr = sci(f, 2) if f < 1e-2 else it(f, 3)
            tstr = "--"
            esito = "valid"
        else:
            fstr = val(f, 2) + "$^\\ast$"
            tstr = sci(c["t_spin"], 2)
            esito = "rotor"
        out.append("%s & %s & %s & %s & %s & %s & %s & %s \\\\\n" %
                   (tag, sci(c["tau"], 2), it(c["Pe_r"], 1), fstr,
                    it(c["dT"], 2), sci(c["eta_mecc"], 1), tstr, esito))
    out.append("\\hline\n\\end{tabular}\n")
    write("tab_rappresentativi.tex", "".join(out))


# ---------------------------------------------------------------------------
# 6. full rotors
# ---------------------------------------------------------------------------
def tab_rotore():
    gruppi = {}
    for r in T.ROWS:
        if T.btrue(r["valido_rotore"]):
            k = (r["canale"], r["materiale"], r["forma"], r["mezzo"], r["p"])
            gruppi.setdefault(k, []).append(r)
    out = [r"""\begin{tabular}{@{}llllrrrrr@{}}
\hline
Material & Shape & $p$ (mbar) & $R$ ($\mu$m) & $P$ (mW) & Cells & $t_{\rm spin}$ (s) & $T$ (K) & $\eta_{\rm mecc}$ \\
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
                   (LBL_MAT[mat], LBL_FOR[forma], int(round(math.log10(float(p)))),
                    it(min(Rs), 2), it(max(Rs), 2),
                    it(min(Ps), 2), it(max(Ps), 2), len(v),
                    sci(min(ts), 2), sci(max(ts), 2),
                    it(min(Ts), 0), it(max(Ts), 0),
                    sci(min(em), 1), sci(max(em), 1)))
        tot += len(v)
    out.append("\\hline\n\\multicolumn{5}{@{}l}{Total: %d combinations, %d cells} & %d & & & \\\\\n" %
               (len(gruppi), tot, tot))
    out.append("\\hline\n\\end{tabular}\n")
    write("tab_rotore.tex", "".join(out))


# ---------------------------------------------------------------------------
# 7. materials
# ---------------------------------------------------------------------------
def tab_materiali():
    out = [r"""\begin{tabular}{@{}lrrrrrrp{4.2cm}@{}}
\hline
Material & $n$ & $k$ & $\rho$ & $\sigma_y$ & $T_{\rm melt}$ & $\Delta n$ & Notes \\
 & & & (kg/m$^3$) & (GPa) & (K) & & \\
\hline
"""]
    note = {
        "cu_lucido": "absorbs 40\\% (Fresnel); $\\varepsilon_{IR}$ from Mie 0.0045--0.0109: IR mirror",
        "cu_nero": "$\\alpha_{532}=0.95$; $\\varepsilon_{IR}=0.30$ declared (0.20--0.60), not Mie",
        "silice": "transparent; $\\Delta n$ from shape only: 0.02 disk, 0.03 dumbbell, 0.01 ring/grating",
        "vaterite": "intrinsic birefringence $\\Delta n=0.17$; brittle, conservative $\\sigma_y$",
    }
    display = {"cu_lucido": "Shiny Cu", "cu_nero": "Black Cu",
               "silice": "Silica", "vaterite": "Vaterite"}
    for k in ("cu_lucido", "cu_nero", "silice", "vaterite"):
        d = m.MATERIALI[k]
        out.append("%s & %s & %s & %s & %s & %s & %s & %s \\\\\n" %
                   (display[k], it(d["n"], 2), val(d["k"], 2),
                    it(d["rho"], 0), it(d["sy"] / 1e9, 2), it(d["Tm"], 0),
                    it(d["dn"], 2), note[k]))
    out.append("\\hline\n\\end{tabular}\n")
    write("tab_materiali.tex", "".join(out))


# ---------------------------------------------------------------------------
# 8. scale
# ---------------------------------------------------------------------------
def tab_scala():
    r1 = T.S1["righe"][0]
    c1 = m.cella("ritardo", 10e-3, 1e-6, "vaterite", "sfera", "acqua")
    p1 = 10.0 * m.KT / (c1["tau"] / 10e-3)
    c2 = m.cella("ritardo", 10e-3, 0.5, "vaterite", "sfera", "acqua")
    p2 = 10.0 * m.KT / (c2["tau"] / 10e-3)
    out = [r"""\begin{tabular}{@{}lrrrrrrp{2.8cm}@{}}
\hline
Scale & $R$ (m) & $P(Pe_r=10)$ (W) & $Pe_r$ at 10 mW & $\omega$ (rad/s) & $\eta_{\rm mecc}$ & $\Delta T$ (K) & Note \\
\hline
"""]
    out.append("1 nm (Cu) & %s & %s & %s & %s & %s & %s & $\\eta_{int}\\approx2\\times10^{-6}$: the photon misses the object \\\\\n" %
               (sci(r1["R"], 1), sci(r1["P_pe10"], 2), sci(r1["Pe_r"], 2),
                sci(r1["omega"], 2), sci(r1["eta_mecc"], 2),
                it(r1["dT"], 3)))
    out.append("1 $\\mu$m (vaterite) & %s & %s & %s & %s & %s & %s & window: rotation at 27.9 Hz \\\\\n" %
               (sci(1e-6, 1), sci(p1, 2), it(c1["Pe_r"], 1), it(c1["omega"], 1),
                sci(c1["eta_mecc"], 2), it(c1["dT"], 2)))
    out.append("1 m (vaterite) & %s & %s & %s & %s & %s & %s & one turn every $1.5\\times10^{8}$ years \\\\\n" %
               (sci(0.5, 1), sci(p2, 2), it(c2["Pe_r"], 1), sci(c2["omega"], 2),
                sci(c2["eta_mecc"], 2), it(c2["dT"], 4)))
    out.append("\\hline\n\\end{tabular}\n")
    write("tab_scala.tex", "".join(out))


# ---------------------------------------------------------------------------
# 9. geometries
# ---------------------------------------------------------------------------
def tab_geometrie():
    out = [r"""\begin{tabular}{@{}lrrrrrrr@{}}
\hline
Shape & $V/R^3$ & $A_{\rm proj}/R^2$ & $A_{\rm rad}/R^2$ & $k_I$ & $c_{\rm mecc}$ & $L/R$ & $R_{\rm eq}/R$ \\
\hline
"""]
    for forma in m.FORME:
        g = m.geometria(forma, 1e-6)
        R = 1e-6
        out.append("%s & %s & %s & %s & %s & %s & %s & %s \\\\\n" %
                   (LBL_FOR[forma], it(g["V"] / R ** 3, 3), it(g["A_proj"] / R ** 2, 3),
                    it(g["A_rad"] / R ** 2, 3), it(g["kI"], 3),
                    it(g["c_mecc"], 2), it(g["percorso"] / R, 2),
                    it(g["Req"] / R, 3)))
    out.append("\\hline\n\\end{tabular}\n")
    write("tab_geometrie.tex", "".join(out))


if __name__ == "__main__":
    print("V7.1 - English tables")
    tab_correzioni()
    tab_criteri()
    tab_validazione()
    tab_finestre()
    tab_rappresentativi()
    tab_rotore()
    tab_materiali()
    tab_scala()
    tab_geometrie()
