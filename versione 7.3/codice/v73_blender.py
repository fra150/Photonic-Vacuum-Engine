"""
MOTORE FOTONICO DEL VUOTO - V7.3
v73_blender.py - le scene del latch a 720 gradi e il cancello sul flusso.

Le proporzioni vengono da dati/V73_sintesi.json (cella_riferimento_3d) e da
dati/V73_taglio.json: R = 0.50 um, bead r = R/2^(1/3) = 0.397 um, d = 2R.
La scena non e' decorativa: il latch e' mostrato a 0, 360 e 720 gradi, con il
contrappeso che passa da un lato all'altro; il dettaglio quota R, r e d.

La scena del cancello mostra il profilo di Poiseuille in un capillare (frecce
proporzionali alla velocita') e il verdetto: flusso escluso, stasi ammessa.

Uso:
    blender -b -P v73_blender.py
"""

import json
import math
import os
import sys

import bpy

BASE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(BASE, ".."))
OUT = os.path.join(ROOT, "blender")
DATI = os.path.join(ROOT, "dati")
os.makedirs(OUT, exist_ok=True)

if BASE not in sys.path:
    sys.path.insert(0, BASE)

import v7_blender as b  # noqa: E402

with open(os.path.join(DATI, "V73_sintesi.json"), encoding="utf-8") as fh:
    S = json.load(fh)
with open(os.path.join(DATI, "V73_taglio.json"), encoding="utf-8") as fh:
    TG = json.load(fh)

C = S["cella_riferimento_3d"]
R = C["R_um"]                       # 0.50 um
BEAD = 0.5 ** (1.0 / 3.0) * R       # 0.397 um
D = 2.0 * R                         # 1.00 um
GAMMA_LIM = TG["gamma_limite_s-1"]  # 2.700 s^-1


# ---------------------------------------------------------------------------
# testo
# ---------------------------------------------------------------------------
def etichetta(nome, corpo, pos, size, colore, cam, emissione=2.2,
              align="CENTER"):
    cu = bpy.data.curves.new(nome, type="FONT")
    cu.body = corpo
    cu.size = size
    cu.align_x = align
    cu.align_y = "CENTER"
    ob = bpy.data.objects.new(nome, cu)
    bpy.context.collection.objects.link(ob)
    ob.location = pos
    ob.rotation_euler = cam.rotation_euler
    ob.data.materials.append(b.materiale(
        "mat_" + nome, colore, ruvidita=0.4, emissione=colore, forza=emissione))
    return ob


def scheda_v73(cam, titolo, righe, dz=6.0, lente=50.0, aspetto=1.3,
               y_frac=0.44):
    W = 36.0 * dz / lente
    H = W / aspetto
    corpo = "\n".join([titolo, ""] + ["  " + r for r in righe])
    cu = bpy.data.curves.new("scheda_v73", type="FONT")
    cu.body = corpo
    cu.size = 0.030 * H
    cu.align_x = "LEFT"
    cu.align_y = "TOP"
    cu.space_line = 0.95
    ob = bpy.data.objects.new("scheda_v73", cu)
    bpy.context.collection.objects.link(ob)
    ob.parent = cam
    ob.location = (-0.47 * W, y_frac * H, -dz)
    ob.data.materials.append(b.materiale(
        "testo_v73", (0.02, 0.02, 0.02), emissione=(0.92, 0.96, 1.0),
        forza=1.8))
    return ob


# ---------------------------------------------------------------------------
# geometria del latch
# ---------------------------------------------------------------------------
def materiali_latch():
    return {
        "cristallo": b.materiale("vaterite", (1.0, 0.80, 0.86),
                                 ruvidita=0.05, alpha=0.30),
        "bead": b.materiale("bead_silice", (0.35, 0.72, 1.0),
                            ruvidita=0.08, alpha=0.55,
                            emissione=(0.10, 0.35, 1.0), forza=0.8),
        "linker": b.materiale("linker", (0.85, 0.90, 1.0), ruvidita=0.2,
                              alpha=0.55),
        "marca": b.materiale("marca_oro", (1.0, 0.72, 0.15),
                             metallo=0.8, ruvidita=0.18,
                             emissione=(1.0, 0.55, 0.05), forza=2.0),
        "fascio": b.materiale("fascio", (0.25, 0.70, 1.0),
                              emissione=(0.20, 0.62, 1.0), forza=1.0,
                              alpha=0.10),
        "elica": b.materiale("elica", (0.15, 0.75, 1.0),
                             emissione=(0.10, 0.65, 1.0), forza=2.2),
        "giro": b.materiale("giro", (0.25, 0.95, 0.70),
                            emissione=(0.15, 0.95, 0.65), forza=3.2),
    }


def scena_latch(gradi, stato, file_nome, dettaglio=False):
    b.pulisci()
    b.mondo((0.008, 0.011, 0.018), 0.55)
    m = materiali_latch()

    th = math.radians(gradi)
    th_bead = math.radians(0.5 * gradi)   # A a 0, B a 360, A a 720
    cx, cy = math.cos(th), math.sin(th)
    bx, by = math.cos(th_bead), math.sin(th_bead)

    # rotore, bead, linker
    b.sfera("rotore", R, (0, 0, 0), m["cristallo"], segmenti=96)
    b.sfera("marca", 0.075, (0.88 * R * cx, 0.88 * R * cy, 0.0),
            m["marca"], segmenti=32)
    b.sfera("bead", BEAD, (D * bx, D * by, 0.0), m["bead"], segmenti=64)
    b.cilindro("linker", 0.028, D, (0.5 * D * bx, 0.5 * D * by, 0.0),
               m["linker"], vertici=24, rot=(0, math.pi / 2, th_bead))

    # fascio ottico lungo z e fronte elicoidale
    b.cilindro("fascio", 0.52, 6.4, (0, 0, 0), m["fascio"], vertici=64)
    b.cono("cono_fuoco", 0.08, 0.52, 2.4, (0, 0, 3.2), m["fascio"], vertici=64)
    elica = []
    for i in range(241):
        t = i / 240.0
        z = -3.0 + 6.0 * t
        a = 3.0 * 2.0 * math.pi * t
        elica.append((0.62 * math.cos(a), 0.62 * math.sin(a), z))
    b.curva("fronte_elicoide", elica, 0.010, m["elica"])

    if not dettaglio:
        # arco di rotazione: spirale, cosi' 720 gradi non si sovrappone
        punti = []
        n = max(64, int(abs(gradi) * 1.2))
        for i in range(n + 1):
            t = i / n
            a = th * t
            punti.append((1.55 * math.cos(a), 1.55 * math.sin(a),
                          0.06 + 0.22 * t * abs(gradi) / 360.0))
        b.curva("arco_giro", punti, 0.030, m["giro"])
        # indice radiale dello stato del rotore, sollevato sopra la bead
        b.curva("gambo_indice", [(0.82 * R * cx, 0.82 * R * cy, 0.24),
                                 (0.82 * R * cx, 0.82 * R * cy, 0.62)],
                0.016, m["giro"])
        b.curva("indice", [(0.78 * R * cx, 0.78 * R * cy, 0.62),
                           (1.46 * R * cx, 1.46 * R * cy, 0.62)],
                0.022, m["giro"])
        b.cono("punta_indice", 0.085, 0.0, 0.26,
               (1.50 * R * cx, 1.50 * R * cy, 0.62), m["giro"],
               rot=(0, math.pi / 2, th), vertici=20)

        # tacche a 0, 90, 180, 270, 360, 450, 540, 630, 720
        for k in range(0, int(abs(gradi)) // 90 + 1):
            a = math.radians(90.0 * k)
            r1, r2 = 1.38, 1.50
            b.curva("tacca_%d" % k,
                    [(r1 * math.cos(a), r1 * math.sin(a), 0.015 * k),
                     (r2 * math.cos(a), r2 * math.sin(a), 0.015 * k)],
                    0.012, m["marca"])

    b.luce("chiave", (5.5, -7.0, 6.5), 900, (1.0, 0.97, 0.93), 5.0, (0, 0, 0))
    b.luce("riempimento", (-6.0, -4.0, 2.5), 260, (0.55, 0.70, 1.0), 6.0, (0, 0, 0))
    b.luce("bordo", (0, 6.0, -3.0), 340, (0.85, 0.55, 1.0), 5.0, (0, 0, 0))

    if dettaglio:
        cam = b.camera("cam", (3.2, -5.0, 2.9), (0.55, 0, 0.0), 50)
    else:
        cam = b.camera("cam", (5.2, -8.2, 4.2), (0.35, 0, 0.1), 50)

    # quote nel dettaglio
    if dettaglio:
        b.curva("quota_R", [(0, 0, 0.62), (R * bx, R * by, 0.62)],
                0.012, m["marca"])
        etichetta("q_R", "R = %.2f um" % R,
                  (0.62 * bx, 0.05, 0.80), 0.085,
                  (1.0, 0.78, 0.25), cam)
        b.curva("quota_d", [(0, 0, -0.62), (D * bx, D * by, -0.62)],
                0.012, m["giro"])
        etichetta("q_d", "d = %.2f um" % D,
                  (0.5 * D * bx, 0.5 * D * by - 0.16, -0.82), 0.085,
                  (0.35, 1.0, 0.65), cam)
        etichetta("q_bead", "bead r = %.2f um" % BEAD,
                  (D * bx, D * by, 0.70), 0.085,
                  (0.45, 0.78, 1.0), cam)

    righe = [
        "R = %.2f um; bead r = %.2f um; d = %.2f um" % (R, BEAD, D),
        "T360 senza contrappeso = %.2f s" % C["T_360_no_s"],
        "T360 con contrappeso   = %.2f s" % C["T_360_cw_s"],
        "rapporto = %.3f; CV = %.3f" % (C["beta_cw_su_beta_no"], 0.1784),
        "rotore = %d gradi; bead = %.0f gradi" % (gradi, 0.5 * gradi),
        "stato %s" % stato,
    ]
    if dettaglio:
        righe = [
            "dettaglio del latch a 720 gradi",
            "rotore di vaterite R = %.2f um" % R,
            "bead di silice r = %.2f um a d = %.2f um" % (BEAD, D),
            "il contrappeso cambia il coefficiente di attrito",
            "beta_cw/beta_no = %.3f" % C["beta_cw_su_beta_no"],
        ]
    scheda_v73(cam, "Latch V7.3 - %s" % stato, righe,
               dz=2.4 if not dettaglio else 1.6, lente=50)

    b.imposta_render(os.path.join(OUT, file_nome),
                     ris=(1100, 850) if not dettaglio else (1000, 800))
    bpy.ops.render.render(write_still=True)
    print("  reso %-28s %6.1f kB"
          % (file_nome, os.path.getsize(os.path.join(OUT, file_nome)) / 1024.0))


# ---------------------------------------------------------------------------
# scena del cancello
# ---------------------------------------------------------------------------
def scatola(nome, dim, pos, mat, alpha=0.10):
    bpy.ops.mesh.primitive_cube_add(size=1.0, location=pos)
    o = bpy.context.object
    o.name = nome
    o.scale = (dim[0], dim[1], dim[2])
    o.data.materials.append(mat)
    return o


def freccia(nome, x0, y, z, lunghezza, mat):
    """Freccia lungo +x, con punta; lunghezza in unita' di scena."""
    if lunghezza < 0.08:
        return
    b.cilindro(nome + "_asta", 0.045, lunghezza, (x0 + 0.5 * lunghezza, y, z),
               mat, vertici=16, rot=(0, math.pi / 2, 0))
    b.cono(nome + "_punta", 0.11, 0.0, 0.22,
           (x0 + lunghezza + 0.11, y, z), mat,
           rot=(0, math.pi / 2, 0), vertici=16)


def scena_cancello():
    b.pulisci()
    b.mondo((0.006, 0.009, 0.015), 0.55)
    mat_vetro = b.materiale("vetro", (0.55, 0.75, 0.90), ruvidita=0.05,
                            alpha=0.10, trasmissione=0.6)
    mat_fluido = b.materiale("fluido", (0.20, 0.55, 0.85), ruvidita=0.1,
                             alpha=0.08)
    mat_freccia = b.materiale("freccia", (0.25, 0.85, 1.0),
                              emissione=(0.20, 0.80, 1.0), forza=2.0)
    mat_rotore = b.materiale("rotore", (1.0, 0.80, 0.86), ruvidita=0.08,
                             alpha=0.45)
    mat_bead = b.materiale("bead", (0.35, 0.72, 1.0), ruvidita=0.1,
                           alpha=0.65, emissione=(0.1, 0.35, 1.0), forza=0.8)
    mat_ok = b.materiale("ok", (0.20, 0.95, 0.45), emissione=(0.15, 0.95, 0.40),
                         forza=4.0)
    mat_no = b.materiale("no", (1.0, 0.20, 0.15), emissione=(1.0, 0.15, 0.10),
                         forza=4.0)
    mat_bianco = b.materiale("bianco", (0.92, 0.95, 1.0),
                             emissione=(0.85, 0.92, 1.0), forza=2.5)

    # ---- sinistra: capillare con flusso di Poiseuille
    x0c, Rc, Lc = -5.0, 3.0, 13.0
    b.cilindro("capillare", Rc, Lc, (x0c, 0, 0), mat_vetro, vertici=96,
               rot=(0, math.pi / 2, 0))
    b.cilindro("fluido", 0.97 * Rc, 0.99 * Lc, (x0c, 0, 0), mat_fluido,
               vertici=96, rot=(0, math.pi / 2, 0))

    # frecce: v(r) = 1 - (r/R)^2, lunghezza proporzionale
    for iy in range(-5, 6):
        for iz in range(-5, 6):
            y = 0.55 * iy
            z = 0.55 * iz
            rr = math.hypot(y, z)
            if rr > 0.90 * Rc:
                continue
            v = 1.0 - (rr / Rc) ** 2
            freccia("f_%d_%d" % (iy, iz), x0c - 0.44 * Lc, y, z,
                    0.35 + 4.4 * v, mat_freccia)

    # rotore bloccato alla parete
    b.sfera("rotore_flusso", R, (x0c + 3.4, 0.72 * Rc, 0.0),
            mat_rotore, segmenti=64)
    b.sfera("bead_flusso", BEAD, (x0c + 3.4 + D, 0.72 * Rc, 0.0),
            mat_bead, segmenti=48)
    b.curva("croce_a", [(x0c + 3.4 - 0.7, 0.72 * Rc, 0.9),
                        (x0c + 3.4 + 0.7, 0.72 * Rc, 1.9)], 0.045, mat_no)
    b.curva("croce_b", [(x0c + 3.4 - 0.7, 0.72 * Rc, 1.9),
                        (x0c + 3.4 + 0.7, 0.72 * Rc, 0.9)], 0.045, mat_no)

    # ---- destra: camera di stasi
    x0s = 7.0
    scatola("camera_stasi", (5.6, 5.6, 5.6), (x0s, 0, 0), mat_vetro, 0.10)
    b.cilindro("fluido_stasi", 2.55, 2.55, (x0s, 0, 0), mat_fluido,
               vertici=64, rot=(0, 0, 0))
    b.sfera("rotore_stasi", R, (x0s, 0, 0), mat_rotore, segmenti=64)
    b.sfera("bead_stasi", BEAD, (x0s + D, 0, 0), mat_bead, segmenti=48)
    b.cilindro("linker_stasi", 0.028, D, (x0s + 0.5 * D, 0, 0),
               mat_bead, vertici=16, rot=(0, math.pi / 2, 0))

    # ---- etichette
    cam = b.camera("cam", (1.0, -27.0, 12.0), (1.0, 0, 0.6), 42)
    etichetta("eti_flusso",
              "FLUSSO: ESCLUSO\ngamma = 100-500 s^-1", (x0c, 0, -4.4), 0.45,
              (1.0, 0.25, 0.18), cam, emissione=3.0)
    etichetta("eti_stasi",
              "STASI: AMMESSO\ngamma_lim = %.2f s^-1" % GAMMA_LIM,
              (x0s, 0, -4.4), 0.45, (0.25, 1.0, 0.45), cam, emissione=3.0)
    etichetta("eti_plasma", "plasma + stasi: solo in vitro", (1.0, 0, -6.6),
              0.34, (0.80, 0.85, 0.95), cam, emissione=2.0)

    b.luce("chiave", (6, -16, 14), 2600, (1.0, 0.97, 0.92), 12, (1, 0, 0))
    b.luce("riempimento", (-12, -10, 5), 900, (0.55, 0.70, 1.0), 12, (-4, 0, 0))
    b.luce("bordo", (12, 8, -4), 1200, (0.75, 0.55, 1.0), 10, (6, 0, 0))

    scheda_v73(cam, "Cancello V7.3 - quattro uscite",
               ["gamma_coppia = %.3f s^-1" % TG["output_1_coppia"]["gamma_vincolato_s-1"],
                "gamma_forza = %.2f-%.2f s^-1" % (
                    min(x["gamma_s-1"] for x in TG["output_2_forza"]),
                    max(x["gamma_s-1"] for x in TG["output_2_forza"])),
                "gamma_lim = %.3f s^-1" % GAMMA_LIM,
                "venula x%.1f; capillare x%.0f" % (
                    TG["confronto_fisiologico"]["venula"]["rapporto"],
                    TG["confronto_fisiologico"]["capillare"]["rapporto"]),
                "il latch vive in stasi, non nel flusso"],
               dz=10.5, lente=50, aspetto=1.35)

    b.imposta_render(os.path.join(OUT, "v73_cancello.png"), ris=(1300, 820))
    bpy.ops.render.render(write_still=True)
    print("  reso %-28s %6.1f kB"
          % ("v73_cancello.png",
             os.path.getsize(os.path.join(OUT, "v73_cancello.png")) / 1024.0))

    bpy.ops.wm.save_as_mainfile(
        filepath=os.path.join(OUT, "motore_fotonico_v73.blend"))


# ---------------------------------------------------------------------------
def main():
    print("V7.3 - Blender")
    scena_latch(0, "A (0 gradi)", "v73_latch_000.png")
    scena_latch(360, "B (360 gradi)", "v73_latch_360.png")
    scena_latch(720, "A (720 gradi)", "v73_latch_720.png")
    scena_latch(720, "dettaglio", "v73_latch_dettaglio.png", dettaglio=True)
    scena_cancello()
    print("Fatto: blender in", OUT)


if __name__ == "__main__":
    main()
