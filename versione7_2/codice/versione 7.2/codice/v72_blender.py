"""
MOTORE FOTONICO DEL VUOTO - V7.2
v72_blender.py - il ciclo di 720 gradi: armatura, contrappeso, collegamento.

Tre fotogrammi: 0, 360 e 720 gradi.  L'armatura e il contrappeso tornano nello
stesso posto in tutti e tre: cio' che cambia e' lo STATO del collegamento.
A 360 gradi A e B sono scambiati, a 720 tornano.  E' l'unico contenuto fisico
del ciclo di 720 gradi: una proprieta' topologica del collegamento
(SO(3) -> SU(2)), non una sorgente di momento angolare.

ATTENZIONE: i due cavi sono un diagramma, non un meccanismo letterale.  Nel
sangue il "braccio" e' un linker molecolare e il "latch" e' un vincolo
reversibile (spot ottico secondario, linker DNA, bead magnetica).  Le schede
lo dicono; il modello 3D non deve suggerire che esistano cavi.

Unita' di scena: 1 unita' = 0.1 um.  Rotore R = 0.3 um = 3 unita';
contrappeso V/2 (r = 0.238 um = 2.38 unita') a d = 2 R = 6 unita'.

I numeri delle schede vengono da dati/V72_sintesi.json (cella_riferimento_3d),
cioe' da v72_contrappeso.py.  Nessun numero e' scritto a mano.

Uso:
    blender -b -P v72_blender.py
"""

import json
import math
import os
import sys

import bpy

BASE = os.path.dirname(os.path.abspath(__file__))
V7COD = os.path.normpath(os.path.join(BASE, "..", "..", "versione 7", "codice"))
DATI = os.path.normpath(os.path.join(BASE, "..", "dati"))
OUT = os.path.normpath(os.path.join(BASE, "..", "blender"))
for p in (V7COD, BASE):
    if p not in sys.path:
        sys.path.insert(0, p)
os.makedirs(OUT, exist_ok=True)

import v7_blender as v7b  # noqa: E402

v7b.OUT = OUT


def dati_riferimento():
    with open(os.path.join(DATI, "V72_sintesi.json"), encoding="utf-8") as fh:
        return json.load(fh)["cella_riferimento_3d"]


def scheda_fase(fase, rif):
    comune = "cavi = diagramma del collegamento, non meccanismo nel sangue"
    if fase == 0:
        righe = ["stato A: il collegamento e' A",
                 "contrappeso disinnestato",
                 "I = I_rot, beta = beta_rot",
                 "Pe_r = 10 nel sangue",
                 "800 nm, z = %.0f um" % rif["z_um"],
                 comune]
        titolo = "V7.2 - ciclo di 720 gradi: stato A (0°)"
    elif fase == 360:
        righe = ["stato B: A e B sono scambiati",
                 "il contrappeso si innesta",
                 "beta x %.2f, omega -%.1f%%"
                 % (rif["costo_attrito"], rif["modulazione"] * 100),
                 "T_360(cw) = %.2f s, T_720 = %.2f s"
                 % (rif["T_360_s"], rif["T_720_s"]),
                 "il latch e' un vincolo reversibile",
                 comune]
        titolo = "V7.2 - 360 gradi: stato B"
    else:
        righe = ["stato A: il collegamento torna",
                 "T_360(no) = %.2f s, T_360(cw) = %.2f s"
                 % (rif["T_360_no_s"], rif["T_360_cw_s"]),
                 "T_720 = %.2f s, E_720 = %.1f mJ"
                 % (rif["T_720_s"], rif["E_720_mJ"]),
                 "dT = %.1f K, P_target = %.0f uW"
                 % (rif["dT_K"], rif["P_target_uW"]),
                 "il 720 e' il periodo dello STATO",
                 comune]
        titolo = "V7.2 - 720 gradi: stato A, ciclo chiuso"
    return {"titolo": titolo, "righe": righe}


def corda(nome, punti, mat, bevel=0.045):
    return v7b.curva(nome, punti, bevel, mat)


def scena(fase_deg, rif):
    v7b.pulisci()
    v7b.mondo((0.008, 0.011, 0.018), 0.5)
    mat_rot = v7b.materiale("vaterite", (1.0, 0.80, 0.86), ruvidita=0.06,
                            alpha=0.34)
    mat_cw = v7b.materiale("contrappeso", (0.75, 0.78, 0.85), metallo=0.6,
                           ruvidita=0.35)
    mat_braccio = v7b.materiale("braccio", (0.35, 0.38, 0.45), metallo=0.8,
                                ruvidita=0.3)
    mat_cavo_a = v7b.materiale("cavo_a", (0.2, 0.85, 1.0),
                               emissione=(0.15, 0.8, 1.0), forza=3.0)
    mat_cavo_b = v7b.materiale("cavo_b", (1.0, 0.65, 0.2),
                               emissione=(1.0, 0.6, 0.15), forza=3.0)
    mat_anello = v7b.materiale("anello", (0.45, 0.48, 0.55), metallo=0.9,
                               ruvidita=0.25)
    mat_latch_a = v7b.materiale("latch_a", (0.15, 0.9, 1.0),
                                emissione=(0.1, 0.8, 1.0), forza=6.0)
    mat_latch_b = v7b.materiale("latch_b", (1.0, 0.55, 0.1),
                                emissione=(1.0, 0.45, 0.05), forza=6.0)

    R = 3.0
    d = 6.0                                # 2 R: il contrappeso non compenetra
    th = math.radians(fase_deg)
    tip = (d * math.cos(th), d * math.sin(th), 0.0)

    v7b.sfera("rotore", R, (0, 0, 0), mat_rot, segmenti=72)
    v7b.sfera("contrappeso", 2.38, (d * math.cos(th), d * math.sin(th), 0.0),
              mat_cw, segmenti=56)
    v7b.cilindro("braccio", 0.26, 0.9,
                 (3.35 * math.cos(th), 3.35 * math.sin(th), 0.0),
                 mat_braccio, rot=(0, math.pi / 2, th), vertici=24)

    # anello fisso (il telaio) e i due cavi che si scambiano ogni 360 gradi
    v7b.arco("telaio", 14.0, 0.0, 0, 2 * math.pi, 0.07, mat_anello, n=180)
    z_swap = 1.5 * math.cos(0.5 * th)      # 0 -> +1.5, 360 -> -1.5, 720 -> +1.5
    for segno, mat, nome in ((+1, mat_cavo_a, "cavo_a"),
                             (-1, mat_cavo_b, "cavo_b")):
        a_anchor = th + segno * math.radians(40)
        anchor = (14.0 * math.cos(a_anchor), 14.0 * math.sin(a_anchor), 0.0)
        p0 = (tip[0] + segno * 1.2 * -math.sin(th),
              tip[1] + segno * 1.2 * math.cos(th), 0.0)
        p1 = (0.5 * (p0[0] + anchor[0]), 0.5 * (p0[1] + anchor[1]),
              segno * z_swap)
        p2 = (0.82 * anchor[0] + 0.18 * p0[0],
              0.82 * anchor[1] + 0.18 * p0[1], segno * 0.45 * z_swap)
        corda(nome, [p0, p1, p2, anchor], mat)

    # indicatore A/B del latch: due semicirconferenze che si scambiano.
    # A 0 e 720 il colore A e' in alto; a 360 il colore B prende il sopravvento.
    swap = 1 if fase_deg == 360 else 0
    mat_alto = mat_latch_b if swap else mat_latch_a
    mat_basso = mat_latch_a if swap else mat_latch_b
    v7b.arco("latch_alto", 1.55, 0.22, 0.0, math.pi, 0.12, mat_alto, n=72)
    v7b.arco("latch_basso", 1.55, 0.22, math.pi, 2.0 * math.pi, 0.12,
             mat_basso, n=72)

    v7b.luce("chiave", (22, -30, 24), 24000, (1.0, 0.97, 0.93), 14,
             (0, 0, 1))
    v7b.luce("riempimento", (-26, -20, 10), 8000, (0.55, 0.7, 1.0), 16,
             (0, 0, 1))
    cam = v7b.camera("cam", (0.0, -46.0, 22.0), (0, 0, 0.5), 50)
    v7b.scheda(cam, scheda_fase(fase_deg, rif), dz=42.0, lente=50,
               aspetto=1100 / 850.0)


def scena_dettaglio(rif):
    """Primo piano: il latch A/B e la bead, per chiarire cosa non e' un cavo."""
    v7b.pulisci()
    v7b.mondo((0.006, 0.009, 0.016), 0.4)
    mat_rot = v7b.materiale("vaterite", (1.0, 0.80, 0.86), ruvidita=0.06,
                            alpha=0.30)
    mat_cw = v7b.materiale("contrappeso", (0.75, 0.78, 0.85), metallo=0.6,
                           ruvidita=0.35)
    mat_braccio = v7b.materiale("braccio", (0.35, 0.38, 0.45), metallo=0.8,
                                ruvidita=0.3)
    mat_a = v7b.materiale("latch_a", (0.15, 0.9, 1.0),
                          emissione=(0.1, 0.8, 1.0), forza=8.0)
    mat_b = v7b.materiale("latch_b", (1.0, 0.55, 0.1),
                          emissione=(1.0, 0.45, 0.05), forza=8.0)
    v7b.sfera("rotore", 3.0, (0, 0, 0), mat_rot, segmenti=72)
    v7b.sfera("bead", 2.38, (6.0, 0, 0), mat_cw, segmenti=56)
    v7b.cilindro("linker", 0.18, 0.9, (3.35, 0, 0), mat_braccio,
                 rot=(0, math.pi / 2, 0), vertici=24)
    v7b.arco("latch_a", 1.55, 0.22, 0.0, math.pi, 0.12, mat_a, n=72)
    v7b.arco("latch_b", 1.55, 0.22, math.pi, 2.0 * math.pi, 0.12, mat_b,
             n=72)
    v7b.luce("chiave", (10, -18, 16), 18000, (1.0, 0.97, 0.93), 12,
             (0, 0, 0))
    v7b.luce("riempimento", (-14, -10, 6), 6000, (0.55, 0.7, 1.0), 12,
             (0, 0, 0))
    cam = v7b.camera("cam", (0.0, -24.0, 10.0), (2.0, 0, 0.2), 65)
    v7b.scheda(cam, {
        "titolo": "V7.2 - cosa sono braccio e latch",
        "righe": ["braccio = linker molecolare di lunghezza 2 R",
                  "contrappeso = bead di volume V/2",
                  "latch = vincolo reversibile, non un cavo",
                  "A/B = stato del collegamento (periodo 720°)",
                  "realizzazione: spot ottico, DNA, bead magnetica"],
    }, dz=30.0, lente=65, aspetto=1100 / 850.0)


def main():
    print("V7.2 - Blender")
    rif = dati_riferimento()
    for fase, nome in ((0, "v72_ciclo_000.png"),
                       (360, "v72_ciclo_360.png"),
                       (720, "v72_ciclo_720.png")):
        scena(fase, rif)
        v7b.render(nome, ris=(1100, 850))
        if fase == 0:
            bpy.ops.wm.save_as_mainfile(
                filepath=os.path.join(OUT, "motore_fotonico_v72.blend"))
    scena_dettaglio(rif)
    v7b.render("v72_latch.png", ris=(1100, 850))
    print("Fatto: blender in", OUT)


if __name__ == "__main__":
    main()
