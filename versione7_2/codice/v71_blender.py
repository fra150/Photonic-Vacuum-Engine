"""
MOTORE FOTONICO DEL VUOTO - V7.1
v71_blender.py - il ciclo di 720 gradi: armatura, contrappeso, collegamento.

Tre fotogrammi: 0, 360 e 720 gradi.  L'armatura e il contrappeso tornano nello
stesso posto in tutti e tre: cio' che cambia e' il COLLEGAMENTO.  A 360 gradi i
due cavi sono scambiati, a 720 tornano.  E' l'unico contenuto fisico del ciclo
di 720 gradi: una proprieta' topologica del collegamento (SO(3) -> SU(2)),
non una sorgente di momento angolare.

Unita' di scena: 1 unita' = 0.1 um.  Rotore R = 0.3 um = 3 unita';
contrappeso V/2 (r = 0.238 um = 2.38 unita') a d = 1.5 R = 4.5 unita'.

I numeri delle schede vengono da dati/V71_sintesi.json (cella_riferimento_3d),
cioe' da v71_contrappeso.py.  Nessun numero e' scritto a mano.

Uso:
    blender -b -P v71_blender.py
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
    with open(os.path.join(DATI, "V71_sintesi.json"), encoding="utf-8") as fh:
        return json.load(fh)["cella_riferimento_3d"]


def scheda_fase(fase, rif):
    if fase == 0:
        righe = ["armatura a 0 gradi",
                 "contrappeso disinnestato",
                 "I = I_rot, beta = beta_rot",
                 "Pe_r = 10 nel sangue",
                 "800 nm, z = %.0f um" % rif["z_um"]]
        titolo = "V7.1 - ciclo di 720 gradi: 0"
    elif fase == 360:
        righe = ["armatura di nuovo a 0 gradi",
                 "il contrappeso si innesta",
                 "I x %.2f, beta x %.2f" % (rif["costo_inerzia"], rif["costo_attrito"]),
                 "omega cala del %.0f%%" % (rif["modulazione"] * 100),
                 "il collegamento e' invertito"]
        titolo = "V7.1 - 360 gradi: meta' ciclo"
    else:
        righe = ["armatura di nuovo a 0 gradi",
                 "il collegamento torna come all'inizio",
                 "T_360 = %.2f s, T_720 = %.2f s" % (rif["T_360_s"], rif["T_720_s"]),
                 "E_720 = %.1f mJ, dT = %.1f K" % (rif["E_720_mJ"], rif["dT_K"]),
                 "P_target = %.0f uW" % rif["P_target_uW"]]
        titolo = "V7.1 - 720 gradi: ciclo chiuso"
    return {"titolo": titolo, "righe": righe}


def corda(nome, punti, mat, bevel=0.045):
    return v7b.curva(nome, punti, bevel, mat)


def scena(fase_deg, rif):
    v7b.pulisci()
    v7b.mondo((0.008, 0.011, 0.018), 0.5)
    mat_rot = v7b.materiale("vaterite", (1.0, 0.80, 0.86), ruvidita=0.06, alpha=0.34)
    mat_cw = v7b.materiale("contrappeso", (0.75, 0.78, 0.85), metallo=0.6, ruvidita=0.35)
    mat_braccio = v7b.materiale("braccio", (0.35, 0.38, 0.45), metallo=0.8, ruvidita=0.3)
    mat_cavo_a = v7b.materiale("cavo_a", (0.2, 0.85, 1.0), emissione=(0.15, 0.8, 1.0),
                               forza=3.0)
    mat_cavo_b = v7b.materiale("cavo_b", (1.0, 0.65, 0.2), emissione=(1.0, 0.6, 0.15),
                               forza=3.0)
    mat_anello = v7b.materiale("anello", (0.45, 0.48, 0.55), metallo=0.9, ruvidita=0.25)

    R = 3.0
    d = 6.0                                # 2 R: il contrappeso non compenetra
    th = math.radians(fase_deg)
    tip = (d * math.cos(th), d * math.sin(th), 0.0)

    v7b.sfera("rotore", R, (0, 0, 0), mat_rot, segmenti=72)
    v7b.sfera("contrappeso", 2.38, (d * math.cos(th), d * math.sin(th), 0.0),
              mat_cw, segmenti=56)
    v7b.cilindro("braccio", 0.26, 0.9, (3.35 * math.cos(th), 3.35 * math.sin(th), 0.0),
                 mat_braccio, rot=(0, math.pi / 2, th), vertici=24)

    # anello fisso (il telaio) e i due cavi che si scambiano ogni 360 gradi
    v7b.arco("telaio", 14.0, 0.0, 0, 2 * math.pi, 0.07, mat_anello, n=180)
    z_swap = 1.5 * math.cos(0.5 * th)      # 0 -> +1.5, 360 -> -1.5, 720 -> +1.5
    for segno, mat, nome in ((+1, mat_cavo_a, "cavo_a"), (-1, mat_cavo_b, "cavo_b")):
        a_anchor = th + segno * math.radians(40)
        anchor = (14.0 * math.cos(a_anchor), 14.0 * math.sin(a_anchor), 0.0)
        p0 = (tip[0] + segno * 1.2 * -math.sin(th),
              tip[1] + segno * 1.2 * math.cos(th), 0.0)
        p1 = (0.5 * (p0[0] + anchor[0]), 0.5 * (p0[1] + anchor[1]), segno * z_swap)
        p2 = (0.82 * anchor[0] + 0.18 * p0[0], 0.82 * anchor[1] + 0.18 * p0[1],
              segno * 0.45 * z_swap)
        corda(nome, [p0, p1, p2, anchor], mat)

    v7b.luce("chiave", (22, -30, 24), 24000, (1.0, 0.97, 0.93), 14, (0, 0, 1))
    v7b.luce("riempimento", (-26, -20, 10), 8000, (0.55, 0.7, 1.0), 16, (0, 0, 1))
    cam = v7b.camera("cam", (0.0, -46.0, 22.0), (0, 0, 0.5), 50)
    v7b.scheda(cam, scheda_fase(fase_deg, rif), dz=42.0, lente=50,
               aspetto=1100 / 850.0)


def main():
    print("V7.1 - Blender")
    rif = dati_riferimento()
    for fase, nome in ((0, "v71_ciclo_000.png"), (360, "v71_ciclo_360.png"),
                       (720, "v71_ciclo_720.png")):
        scena(fase, rif)
        v7b.render(nome, ris=(1100, 850))
        if fase == 0:
            bpy.ops.wm.save_as_mainfile(
                filepath=os.path.join(OUT, "motore_fotonico_v71.blend"))
    print("Fatto: blender in", OUT)


if __name__ == "__main__":
    main()
