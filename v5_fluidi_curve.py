"""
MOTORE FOTONICO V5 - CURVE DEI FLUSSI DI ROTAZIONE
Autore: Dr. Francesco Bulla / Sistema GNRX

Estrae dalle immagini analitiche (vista ortografica dall'alto, solo tracciante)
le curve che descrivono l'interazione fra il flusso di Spin e quello di
Contro-spin nel vuoto:

    massa_ottica(t)      somma della luminanza del tracciante (proxy della massa)
    raggio_medio(t)      raggio pesato della distribuzione (spreading centrifugo)
    turbolenza(t)        coefficiente di variazione azimutale della luminanza
    penetrazione(t)      quanto ciascun tracciante supera la linea mediana

Uso:
    python v5_fluidi_curve.py
"""

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fluidi")
ANALISI = os.path.join(BASE, "analisi")
DATI = os.path.join(BASE, "dati")
FIGURE = os.path.join(BASE, "figure")
os.makedirs(DATI, exist_ok=True)
os.makedirs(FIGURE, exist_ok=True)

FPS = 24.0
PASSO = 4
ORTHO = 16.6          # larghezza in unita' della vista ortografica
LATO_PX = 640
SEPARAZIONE = 4.40
NM_PER_UNITA = 100.0
BINS = 36

ARANCIO = "#d2761f"
BLU = "#1f6fb5"
GRIGIO = "#4a4a4a"

plt.rcParams.update({"figure.dpi": 160, "savefig.dpi": 160, "font.size": 9.5,
                     "axes.titlesize": 10.5, "axes.labelsize": 9.5,
                     "axes.grid": True, "grid.color": "#dcdcdc",
                     "legend.frameon": False, "figure.facecolor": "white"})


def mediana(serie, finestra=3):
    """Filtro mediano a finestra dispari: toglie i picchi isolati di penetrazione."""
    serie = np.asarray(serie, dtype=float)
    if finestra < 3 or serie.size < finestra:
        return serie
    meta = finestra // 2
    imbottita = np.pad(serie, meta, mode="edge")
    return np.array([np.median(imbottita[i:i + finestra])
                     for i in range(serie.size)])


def carica(f):
    img = mpimg.imread(os.path.join(ANALISI, "tracciante_f%03d.png" % f))
    return img[:, :, :3].astype(float)


def campi(img):
    """Separa i due traccianti in base alla dominante di colore.

    Al fondo si sottrae un pavimento di 0.06: e' la somma della tinta del mondo
    (0.008, 0.010, 0.016) e dell'alone volumetrico a bassa densita' che copre
    l'intero dominio. Senza il pavimento il centroide del campo blu cade a
    x = 3.55 invece che sul centro del rotore (4.40).
    """
    r, g, b = img[:, :, 0], img[:, :, 1], img[:, :, 2]
    arancio = np.clip(r - b - 0.06, 0.0, None)
    blu = np.clip(b - r - 0.06, 0.0, None)
    return arancio, blu


def coordinate():
    h, w = 360, LATO_PX
    ppu = w / ORTHO
    x = (np.arange(w) - w / 2.0) / ppu
    y = (h / 2.0 - np.arange(h)) / ppu
    return np.meshgrid(x, y)


def misura(img, xx, yy):
    arancio, blu = campi(img)
    out = {"massa": {}, "raggio": {}, "dispersione": {}, "turbolenza": {},
           "penetrazione": {}}
    for nome, campo, centro, verso in (("spin", arancio, -SEPARAZIONE, +1.0),
                                       ("contro", blu, +SEPARAZIONE, -1.0)):
        peso = campo
        tot = float(peso.sum())
        if tot <= 0.0:
            for chiave in out:
                out[chiave][nome] = 0.0
            continue
        dx = xx - centro
        r = np.hypot(dx, yy)
        out["massa"][nome] = tot
        medio = float((peso * r).sum() / tot)
        out["raggio"][nome] = medio
        out["dispersione"][nome] = float(np.sqrt((peso * (r - medio) ** 2).sum() / tot))
        theta = np.arctan2(yy, dx)
        indici = ((theta + np.pi) / (2.0 * np.pi) * BINS).astype(int) % BINS
        profilo = np.bincount(indici.ravel(), weights=peso.ravel(), minlength=BINS)
        profilo = profilo / max(profilo.mean(), 1e-12)
        out["turbolenza"][nome] = float(profilo.std())
        soglia = 0.05 * peso.max()
        attivo = peso > soglia
        if attivo.any():
            out["penetrazione"][nome] = float(
                np.clip(verso * xx[attivo], 0.0, None).max())
        else:
            out["penetrazione"][nome] = 0.0
    return out


def main():
    frame = sorted(int(f[12:15])
                   for f in os.listdir(ANALISI) if f.startswith("tracciante_f"))
    xx, yy = coordinate()
    serie = {"t": [], "massa": {"spin": [], "contro": []},
             "raggio": {"spin": [], "contro": []},
             "dispersione": {"spin": [], "contro": []},
             "turbolenza": {"spin": [], "contro": []},
             "penetrazione": {"spin": [], "contro": []}}
    ultima = {}
    for f in frame:
        img = carica(f)
        m = misura(img, xx, yy)
        serie["t"].append((f - 1) / FPS)
        for chiave in ("massa", "raggio", "dispersione", "turbolenza",
                       "penetrazione"):
            for lato in ("spin", "contro"):
                serie[chiave][lato].append(m[chiave][lato])
        ultima[f] = img
        print("[curve] frame %3d  r_spin %.3f  r_contro %.3f  turb %.2f/%.2f"
              % (f, m["raggio"]["spin"], m["raggio"]["contro"],
                 m["turbolenza"]["spin"], m["turbolenza"]["contro"]))

    t = np.array(serie["t"])
    dati = {"nota": "curve estratte dalle immagini della simulazione Mantaflow "
                    "(vista ortografica dall'alto, solo tracciante). Unita' di "
                    "lunghezza: 1 unita' = 100 nm. Nessuna gravita', nessuna "
                    "dissipazione viscosa: il solo attrito e' la dissipazione "
                    "numerica alpha = 0.10. massa = luminanza integrata (unita' "
                    "arbitrarie), raggio = unita' di scena.",
            "fps": FPS, "passo_fotogrammi": PASSO,
            "t_s": list(t), "curve": serie}
    with open(os.path.join(DATI, "curve_flussi.json"), "w", encoding="utf-8") as fh:
        json.dump(dati, fh, indent=2, ensure_ascii=False)
    print("[curve] scritto curve_flussi.json")

    # normalizzazioni per il disegno
    for lato in ("spin", "contro"):
        serie["massa"][lato] = list(np.array(serie["massa"][lato])
                                    / max(serie["massa"][lato][-1], 1e-12))
        serie["raggio"][lato] = list(np.array(serie["raggio"][lato]) * NM_PER_UNITA)
        serie["dispersione"][lato] = list(
            np.array(serie["dispersione"][lato]) * NM_PER_UNITA)
        serie["penetrazione"][lato] = list(
            mediana(np.array(serie["penetrazione"][lato]) * NM_PER_UNITA, 3))

    # --- figura 1: dispersione, turbolenza, massa -------------------------
    fig, assi = plt.subplots(1, 3, figsize=(13.4, 3.6))
    a = assi[0]
    a.plot(t, serie["dispersione"]["spin"], color=ARANCIO, lw=1.8,
           label="Spin (+720°)")
    a.plot(t, serie["dispersione"]["contro"], color=BLU, lw=1.8,
           label="Contro-spin (−86.17°)")
    a.set_xlabel("tempo [s]")
    a.set_ylabel("dispersione radiale $\\sigma_r$ [nm]")
    a.set_title("Spessore dell'anello di tracciante")
    a.legend(loc="upper left")

    a = assi[1]
    a.plot(t, serie["turbolenza"]["spin"], color=ARANCIO, lw=1.8, label="Spin")
    a.plot(t, serie["turbolenza"]["contro"], color=BLU, lw=1.8, label="Contro-spin")
    a.set_xlabel("tempo [s]")
    a.set_ylabel("scarto azimutale normalizzato")
    a.set_title("Irregolarità azimutale del fronte")
    a.legend(loc="upper left")

    a = assi[2]
    a.plot(t, serie["massa"]["spin"], color=ARANCIO, lw=1.8, label="Spin")
    a.plot(t, serie["massa"]["contro"], color=BLU, lw=1.8, label="Contro-spin")
    a.set_xlabel("tempo [s]")
    a.set_ylabel("massa ottica (normalizzata a t$_{fine}$)")
    a.set_title("Inventario del tracciante")
    a.set_ylim(0.0, 1.15)
    a.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGURE, "fig_F1_curve_flussi.png"))
    plt.close(fig)
    print("[curve] scritto fig_F1_curve_flussi.png")

    # --- figura 2: penetrazione e profilo assiale -------------------------
    fig, assi = plt.subplots(1, 2, figsize=(11.0, 3.8),
                             gridspec_kw={"width_ratios": [1.25, 1.0]})
    a = assi[0]
    a.plot(t, serie["penetrazione"]["spin"],
           color=ARANCIO, lw=1.8, label="Spin → semipiano del Contro-spin")
    a.plot(t, serie["penetrazione"]["contro"],
           color=BLU, lw=1.8, label="Contro-spin → semipiano dello Spin")
    a.set_xlabel("tempo [s]")
    a.set_ylabel("penetrazione oltre la linea mediana [nm]")
    a.set_title("Interazione all'interfaccia")
    a.legend(loc="upper left")

    a = assi[1]
    for f, stile in ((frame[0], ":"), (frame[len(frame) // 2], "--"),
                     (frame[-1], "-")):
        img = ultima[f]
        arancio, blu = campi(img)
        profilo = (arancio.sum(axis=0) + blu.sum(axis=0))
        profilo = profilo / max(profilo.max(), 1e-12)
        asse_x = (np.arange(LATO_PX) - LATO_PX / 2.0) / (LATO_PX / ORTHO)
        a.plot(asse_x, profilo, stile, color=GRIGIO, lw=1.4,
               label="t = %.2f s" % ((f - 1) / FPS))
    a.axvline(0.0, color="#b00020", lw=0.9, ls="--")
    a.set_xlabel("x [unità, 1 = 100 nm]")
    a.set_ylabel("intensità integrata su y")
    a.set_title("Profilo assiale: due popolazioni distinte")
    a.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGURE, "fig_F2_interazione.png"))
    plt.close(fig)
    print("[curve] scritto fig_F2_interazione.png")


if __name__ == "__main__":
    main()
