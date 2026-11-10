"""
MOTORE FOTONICO DEL VUOTO - V7
v7_blender.py - il rotore della V7: ritardo di fase, non assorbimento.

La scena non e' decorativa: le proporzioni sono quelle del modello e i numeri
delle schede vengono da dati/C1_cards.json, cioe' da v7_scale.py.  Nessun
numero e' scritto a mano in questo file.

  v7_ritardo_hero.png    la sfera di vaterite (R = 2 um) nel fuoco: asse
                         ottico, fascio elicoidale, freccia di rotazione
  v7_ritardo_seq_1/2     due fasi (60 e 120 gradi): il fascio non cambia,
                         cambia l'asse del cristallo.  E' la differenza fra
                         "trascinato dal campo" e "girato dal ritardo".
  v7_manubrio.png        il manubrio di silice da 200 nm del test F (GHz)
  v7_scala_1nm.png       la particella da 1 nm: 1/532 di lunghezza d'onda
  v7_scala_1m.png        il disco da 1 m, con la sagoma umana come scala

Unita' di scena
  1 unita' = 1 um   scene del ritardo e del manubrio
  1 unita' = 1 nm   scena da 1 nm
  1 unita' = 0.05 m scena da 1 m (il disco da 1 m ha raggio 10 unita')

Uso:
    blender -b -P v7_blender.py
"""

import json
import math
import os

import bpy
import mathutils

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(BASE, "..", "blender"))
DATI = os.path.normpath(os.path.join(BASE, "..", "dati"))
os.makedirs(OUT, exist_ok=True)


def schede():
    p = os.path.join(DATI, "C1_cards.json")
    if not os.path.exists(p):
        raise SystemExit("manca dati/C1_cards.json: eseguire prima v7_scale.py")
    with open(p, encoding="utf-8") as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------
# INFRASTRUTTURA
# ---------------------------------------------------------------------------
def pulisci():
    bpy.ops.wm.read_factory_settings(use_empty=True)


def motore_render():
    props = bpy.types.RenderSettings.bl_rna.properties["engine"]
    nomi = [e.identifier for e in props.enum_items]
    for cand in ("BLENDER_EEVEE_NEXT", "BLENDER_EEVEE", "BLENDER_WORKBENCH"):
        if cand in nomi:
            return cand
    return nomi[0]


def imposta_render(percorso, ris=(1100, 850), campioni=64):
    sc = bpy.context.scene
    sc.render.engine = motore_render()
    sc.render.resolution_x, sc.render.resolution_y = ris
    sc.render.resolution_percentage = 100
    sc.render.film_transparent = False
    sc.render.filepath = percorso
    sc.render.image_settings.file_format = "PNG"
    for attr in ("taa_render_samples", "taa_samples"):
        try:
            setattr(sc.eevee, attr, campioni)
        except Exception:
            pass
    try:
        sc.eevee.use_raytracing = True
    except Exception:
        pass


def _set(bsdf, nomi, val):
    for n in nomi:
        if n in bsdf.inputs:
            bsdf.inputs[n].default_value = val
            return True
    return False


def materiale(nome, colore, metallo=0.0, ruvidita=0.3, emissione=None,
              forza=0.0, alpha=1.0, trasmissione=0.0):
    mat = bpy.data.materials.new(nome)
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes.get("Principled BSDF")
    _set(bsdf, ("Base Color",), (*colore, 1.0))
    _set(bsdf, ("Metallic",), metallo)
    _set(bsdf, ("Roughness",), ruvidita)
    _set(bsdf, ("IOR",), 1.5)
    if alpha < 1.0:
        _set(bsdf, ("Alpha",), alpha)
    if trasmissione > 0:
        _set(bsdf, ("Transmission Weight", "Transmission"), trasmissione)
    if emissione is not None:
        _set(bsdf, ("Emission Color", "Emission"), (*emissione, 1.0))
        _set(bsdf, ("Emission Strength",), forza)
    for attr, val in (("blend_method", "BLEND"),
                      ("surface_render_method", "BLENDED")):
        try:
            setattr(mat, attr, val)
        except Exception:
            pass
    return mat


def mondo(colore=(0.010, 0.013, 0.020), forza=0.6):
    w = bpy.data.worlds.new("mondo")
    bpy.context.scene.world = w
    w.use_nodes = True
    w.node_tree.nodes["Background"].inputs[0].default_value = (*colore, 1.0)
    w.node_tree.nodes["Background"].inputs[1].default_value = forza


def luce(nome, pos, energia, colore=(1, 1, 1), size=6.0, target=(0, 0, 0)):
    d = bpy.data.lights.new(nome, type="AREA")
    d.energy = energia
    d.color = colore
    d.size = size
    o = bpy.data.objects.new(nome, d)
    o.location = pos
    bpy.context.collection.objects.link(o)
    verso = mathutils.Vector(target) - mathutils.Vector(pos)
    o.rotation_euler = verso.to_track_quat("-Z", "Y").to_euler()
    return o


def camera(nome, pos, target=(0, 0, 0), lente=55.0):
    d = bpy.data.cameras.new(nome)
    d.lens = lente
    o = bpy.data.objects.new(nome, d)
    o.location = pos
    bpy.context.collection.objects.link(o)
    verso = mathutils.Vector(target) - mathutils.Vector(pos)
    o.rotation_euler = verso.to_track_quat("-Z", "Y").to_euler()
    bpy.context.scene.camera = o
    return o


def scheda(cam, blocco, dz=5.0, lente=50.0, aspetto=1.3):
    """Testo ancorato alla camera, con i numeri presi dal calcolo.

    La posizione e la dimensione si ricavano dal campo visivo alla distanza dz.
    Con sensore 36 mm la larghezza visibile e' W = 36 dz / f e l'altezza
    H = W / aspetto (misurato con un render di prova: a dz = 5 e f = 50 il
    riquadro e' 3.6 x 2.4 per un render 600x400).  Prima stesura: i due fattori
    erano invertiti e il testo usciva 1.5 volte troppo grande e fuori campo.
    """
    W = 36.0 * dz / lente
    H = W / aspetto
    corpo = [blocco["titolo"], ""] + ["  " + r for r in blocco["righe"]]
    cu = bpy.data.curves.new("testo", type="FONT")
    cu.body = "\n".join(corpo)
    cu.size = 0.030 * H
    cu.align_x = "LEFT"
    cu.align_y = "TOP"
    cu.space_line = 0.95
    ob = bpy.data.objects.new("testo", cu)
    bpy.context.collection.objects.link(ob)
    ob.parent = cam
    ob.location = (-0.47 * W, 0.44 * H, -dz)
    ob.data.materials.append(materiale("testo_m", (0.02, 0.02, 0.02),
                                       emissione=(0.92, 0.96, 1.0), forza=1.6))
    return ob


def curva(nome, punti, bevel, mat, risoluzione=12):
    cu = bpy.data.curves.new(nome, type="CURVE")
    cu.dimensions = "3D"
    cu.bevel_depth = bevel
    cu.bevel_resolution = risoluzione
    sp = cu.splines.new("POLY")
    sp.points.add(len(punti) - 1)
    for i, p in enumerate(punti):
        sp.points[i].co = (p[0], p[1], p[2], 1.0)
    o = bpy.data.objects.new(nome, cu)
    bpy.context.collection.objects.link(o)
    o.data.materials.append(mat)
    return o


def arco(nome, raggio, z, a0, a1, bevel, mat, n=96):
    pts = []
    for i in range(n + 1):
        a = a0 + (a1 - a0) * i / n
        pts.append((raggio * math.cos(a), raggio * math.sin(a), z))
    return curva(nome, pts, bevel, mat)


def sfera(nome, raggio, pos, mat, segmenti=64):
    bpy.ops.mesh.primitive_uv_sphere_add(radius=raggio, location=pos,
                                         segments=segmenti, ring_count=segmenti // 2)
    o = bpy.context.object
    o.name = nome
    o.data.materials.append(mat)
    for p in o.data.polygons:
        p.use_smooth = True
    return o


def cilindro(nome, raggio, altezza, pos, mat, vertici=48, rot=(0, 0, 0)):
    bpy.ops.mesh.primitive_cylinder_add(radius=raggio, depth=altezza,
                                        location=pos, vertices=vertici)
    o = bpy.context.object
    o.name = nome
    o.rotation_euler = rot
    o.data.materials.append(mat)
    for p in o.data.polygons:
        p.use_smooth = True
    return o


def cono(nome, r1, r2, altezza, pos, mat, rot=(0, 0, 0), vertici=64):
    bpy.ops.mesh.primitive_cone_add(radius1=r1, radius2=r2, depth=altezza,
                                    location=pos, vertices=vertici)
    o = bpy.context.object
    o.name = nome
    o.rotation_euler = rot
    o.data.materials.append(mat)
    for p in o.data.polygons:
        p.use_smooth = True
    return o


# ---------------------------------------------------------------------------
# 1. IL ROTORE A RITARDO DI FASE (1 unita' = 1 um)
# ---------------------------------------------------------------------------
def scena_ritardo(fase_gradi=0.0):
    pulisci()
    mondo()
    mat_cristallo = materiale("vaterite", (1.0, 0.80, 0.86), ruvidita=0.05,
                              alpha=0.26)
    mat_asse = materiale("asse_ottico", (1.0, 0.6, 0.2), emissione=(1.0, 0.5, 0.12),
                         forza=7.0)
    mat_fascio = materiale("fascio", (0.25, 0.75, 1.0), emissione=(0.25, 0.7, 1.0),
                           forza=1.2, alpha=0.12)
    mat_elica = materiale("elica", (0.15, 0.72, 1.0), emissione=(0.10, 0.62, 1.0),
                          forza=2.0)
    mat_freccia = materiale("freccia", (0.3, 0.9, 0.75), emissione=(0.2, 0.9, 0.7),
                            forza=3.0)
    mat_lente = materiale("lente", (0.55, 0.62, 0.72), metallo=0.9, ruvidita=0.2)

    R = 2.0
    sfera("vaterite_R2um", R, (0, 0, 0), mat_cristallo, segmenti=96)

    fase = math.radians(fase_gradi)
    L = 3.2
    dx, dy = math.cos(fase), math.sin(fase)
    cilindro("asse_ottico", 0.04, 2 * L, (0, 0, 0), mat_asse,
             rot=(0, math.pi / 2, fase))
    for s in (+1, -1):
        cono("punta_asse", 0.13, 0.0, 0.55, (s * L * dx, s * L * dy, 0), mat_asse,
             rot=(0, s * math.pi / 2, fase), vertici=24)

    cilindro("lente", 1.30, 0.26, (0, 0, 7.4), mat_lente)
    cono("cono_luce", 0.14, 1.22, 7.2, (0, 0, 3.7), mat_fascio, vertici=80)
    elica = []
    for i in range(241):
        t = i / 240.0
        z = -4.6 + 9.2 * t
        a = 3.0 * 2 * math.pi * t
        elica.append((0.60 * math.cos(a), 0.60 * math.sin(a), z))
    curva("fronte_elicoidale", elica, 0.022, mat_elica)
    arco("freccia_arco", 3.05, -0.15, math.radians(-55), math.radians(195),
         0.045, mat_freccia)
    a_fine = math.radians(195)
    p = (3.05 * math.cos(a_fine), 3.05 * math.sin(a_fine), -0.15)
    cono("freccia_punta", 0.14, 0.0, 0.5, p, mat_freccia,
         rot=(0, math.pi / 2, a_fine + math.pi / 2), vertici=24)

    luce("chiave", (8, -11, 10), 1400, (1.0, 0.96, 0.92), 8, (0, 0, 1))
    luce("riempimento", (-10, -6, 4), 380, (0.55, 0.7, 1.0), 10, (0, 0, 1))
    luce("bordo", (0, 10, -5), 520, (0.8, 0.55, 1.0), 9, (0, 0, 1))
    cam = camera("cam", (9.0, -13.5, 8.2), (0, 0, 0.7), 50)
    scheda(cam, schede()["rotore"], dz=6.0, lente=50, aspetto=1100 / 850.0)


# ---------------------------------------------------------------------------
# 2. IL MANUBRIO DA 200 nm (1 unita' = 1 um)
# ---------------------------------------------------------------------------
def scena_manubrio():
    pulisci()
    mondo((0.008, 0.010, 0.016), 0.5)
    mat_sil = materiale("silice", (0.85, 0.93, 1.0), ruvidita=0.04, alpha=0.45)
    mat_fascio = materiale("fascio", (0.3, 0.75, 1.0), emissione=(0.3, 0.75, 1.0),
                           forza=1.4, alpha=0.12)
    mat_freccia = materiale("freccia", (0.3, 0.9, 0.7), emissione=(0.2, 0.9, 0.72),
                            forza=5.0)

    sfera("lobo_a", 0.025, (-0.075, 0, 0), mat_sil, segmenti=48)
    sfera("lobo_b", 0.025, (0.075, 0, 0), mat_sil, segmenti=48)
    cono("cono_luce", 0.008, 0.075, 0.55, (0, 0, 0.34), mat_fascio, vertici=64)
    arco("freccia", 0.100, 0.0, math.radians(-50), math.radians(200), 0.0016,
         mat_freccia)
    cono("punta", 0.007, 0.0, 0.025,
         (0.100 * math.cos(math.radians(200)), 0.100 * math.sin(math.radians(200)), 0),
         mat_freccia, rot=(0, math.pi / 2, math.radians(200) + math.pi / 2), vertici=20)

    luce("chiave", (0.35, -0.55, 0.45), 2.5, (1.0, 0.97, 0.93), 0.4)
    luce("riempimento", (-0.5, -0.4, 0.2), 0.8, (0.5, 0.7, 1.0), 0.5)
    cam = camera("cam", (0.32, -0.50, 0.24), (0, 0, 0.0), 50)
    scheda(cam, schede()["manubrio"], dz=0.62, lente=50, aspetto=1000 / 780.0)


# ---------------------------------------------------------------------------
# 3. LA SCALA: 1 nm (1 unita' = 1 nm)
# ---------------------------------------------------------------------------
def scena_1nm():
    pulisci()
    mondo((0.006, 0.008, 0.014), 0.5)
    mat_p = materiale("particella", (1.0, 0.85, 0.3), emissione=(1.0, 0.78, 0.2),
                      forza=6.0)
    mat_anello = materiale("anello", (0.35, 0.65, 1.0), emissione=(0.3, 0.6, 1.0),
                           forza=2.0)

    sfera("particella_1nm", 0.5, (0, 0, 0), mat_p, segmenti=64)
    arco("anello", 2.2, 0.0, 0, 2 * math.pi, 0.012, mat_anello, n=160)

    luce("chiave", (3.5, -5.0, 4.5), 320, (1, 1, 1), 2.5)
    cam = camera("cam", (0, -7.5, 3.0), (0, 0, 0), 50)
    scheda(cam, schede()["1_nm"], dz=8.0, lente=50, aspetto=1000 / 700.0)


# ---------------------------------------------------------------------------
# 4. LA SCALA: 1 m (1 unita' = 0.05 m)
# ---------------------------------------------------------------------------
def scena_1m():
    pulisci()
    mondo((0.006, 0.008, 0.014), 0.4)
    mat_disco = materiale("disco", (0.78, 0.80, 0.86), metallo=0.8, ruvidita=0.4)
    mat_ago = materiale("ago", (0.3, 0.8, 1.0), emissione=(0.3, 0.8, 1.0), forza=12.0)
    mat_umano = materiale("umano", (0.30, 0.32, 0.38), ruvidita=0.7)

    cilindro("disco_1m", 10.0, 0.9, (0, 0, 0), mat_disco, vertici=128,
             rot=(math.radians(90), 0, 0))
    cilindro("ago_fascio", 0.10, 26.0, (0, -13.0, 0), mat_ago, vertici=24,
             rot=(math.radians(90), 0, 0))
    # sagoma umana da 1.7 m = 34 unita': la scala del mondo macroscopico
    cilindro("umano_corpo", 1.1, 27.0, (15.5, 0, 10.5), mat_umano, vertici=32)
    sfera("umano_testa", 2.1, (15.5, 0, 26.2), mat_umano, segmenti=48)

    luce("chiave", (24, -40, 30), 18000, (1, 0.98, 0.95), 26, (0, 0, 8))
    luce("riempimento", (-30, -26, 10), 5000, (0.5, 0.7, 1.0), 30, (0, 0, 8))
    cam = camera("cam", (44, -56, 34), (0, 0, 11), 50)
    scheda(cam, schede()["1_m"], dz=70.0, lente=50, aspetto=1000 / 760.0)


# ---------------------------------------------------------------------------
# RENDER
# ---------------------------------------------------------------------------
def render(nome, ris=(1100, 850)):
    p = os.path.join(OUT, nome)
    imposta_render(p, ris)
    bpy.ops.render.render(write_still=True)
    print("  reso %-26s %6.1f kB" % (nome, os.path.getsize(p) / 1024.0))


def main():
    print("V7 - Blender")
    schede()                     # fallisce subito se mancano i dati

    scena_ritardo(0.0)
    render("v7_ritardo_hero.png")
    bpy.ops.wm.save_as_mainfile(filepath=os.path.join(OUT, "motore_fotonico_v7.blend"))

    for i, fase in enumerate((60.0, 120.0), start=1):
        scena_ritardo(fase)
        render("v7_ritardo_seq_%d.png" % i, ris=(950, 730))

    scena_manubrio()
    render("v7_manubrio.png", ris=(1000, 780))
    scena_1nm()
    render("v7_scala_1nm.png", ris=(1000, 700))
    scena_1m()
    render("v7_scala_1m.png", ris=(1000, 760))
    print("Fatto: blender in", OUT)


if __name__ == "__main__":
    main()
