"""
MOTORE FOTONICO DEL VUOTO - V6
Scena Blender: il rotore otticamente guidato e le sei geometrie a pari massa.

La scena non e' un'illustrazione decorativa: le proporzioni sono quelle del
modello (1 unita' = 100 nm), i fili hanno il raggio della geometria V5 e le
sei varianti hanno la stessa massa della configurazione di riferimento.

Render prodotti (in ../blender/):
  rotore_hero.png        il reticolo attraversato dal fascio elicoidale
  rotore_seq_0/1/2.png   tre orientazioni (0, 120, 240 gradi)
  base_1..6.png          le sei geometrie a pari massa, stesso ingombro in
                         immagine, dimensione reale nella didascalia

Uso:
    blender -b -P v6_blender.py
"""

import math
import os
import sys

import bpy
import mathutils

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(BASE, "..", "blender"))
os.makedirs(OUT, exist_ok=True)

UA = 1.0                 # 1 unita' di scena = 100 nm
R_DISC = 3.275
R_OUT = 3.70
CORE_R = 0.40
WIRE_ORB = 0.032         # 3.2 nm di raggio (diametro 6.4 nm)
WIRE_RIM = 0.075         # 7.5 nm (diametro 15 nm)
WIRE_LENS = 0.045        # 4.5 nm (diametro 9 nm)
WIRE_R = WIRE_ORB
SPIN_ORBITS = ((2.20, 1.10), (2.95, 1.52), (3.70, 1.95))
LENS = (1.05, 2.45)


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


def imposta_render(percorso, ris=(960, 760), campioni=48):
    sc = bpy.context.scene
    sc.render.engine = motore_render()
    sc.render.resolution_x, sc.render.resolution_y = ris
    sc.render.resolution_percentage = 100
    sc.render.film_transparent = False
    sc.render.filepath = percorso
    sc.render.image_settings.file_format = "PNG"
    try:
        sc.eevee.taa_render_samples = campioni
    except Exception:
        pass
    try:
        sc.eevee.use_raytracing = True
    except Exception:
        pass


def materiale(nome, colore, metallo=1.0, ruvidita=0.28, emissione=None, forza=0.0,
              alpha=1.0):
    m = bpy.data.materials.new(nome)
    m.use_nodes = True
    bsdf = m.node_tree.nodes.get("Principled BSDF")
    bsdf.inputs["Base Color"].default_value = (*colore, 1.0)
    bsdf.inputs["Metallic"].default_value = metallo
    bsdf.inputs["Roughness"].default_value = ruvidita
    if "Alpha" in bsdf.inputs:
        bsdf.inputs["Alpha"].default_value = alpha
    if emissione is not None:
        for chiave in ("Emission Color", "Emission"):
            if chiave in bsdf.inputs:
                bsdf.inputs[chiave].default_value = (*emissione, 1.0)
                break
        if "Emission Strength" in bsdf.inputs:
            bsdf.inputs["Emission Strength"].default_value = forza
    if alpha < 1.0:
        m.blend_method = "BLEND"
    return m


def mondo(colore=(0.035, 0.045, 0.060)):
    w = bpy.data.worlds.new("mondo")
    bpy.context.scene.world = w
    w.use_nodes = True
    w.node_tree.nodes["Background"].inputs[0].default_value = (*colore, 1.0)
    w.node_tree.nodes["Background"].inputs[1].default_value = 1.0


def luce(nome, pos, energia, colore=(1, 1, 1), size=6.0):
    d = bpy.data.lights.new(nome, type="AREA")
    d.energy = energia
    d.color = colore
    d.size = size
    o = bpy.data.objects.new(nome, d)
    o.location = pos
    bpy.context.collection.objects.link(o)
    verso = mathutils.Vector((0, 0, 0)) - mathutils.Vector(pos)
    o.rotation_euler = verso.to_track_quat("-Z", "Y").to_euler()
    return o


def camera(nome, pos, target=(0, 0, 0), lente=60.0):
    d = bpy.data.cameras.new(nome)
    d.lens = lente
    d.clip_start = 0.01
    d.clip_end = 1000
    o = bpy.data.objects.new(nome, d)
    o.location = pos
    bpy.context.collection.objects.link(o)
    o.rotation_euler = (mathutils.Vector(target) - mathutils.Vector(pos)).to_track_quat("-Z", "Y").to_euler()
    bpy.context.scene.camera = o
    return o


def toro(nome, maggiore, minore, materiale_, scala=(1, 1, 1), pos=(0, 0, 0),
         rot=(0, 0, 0), genitore=None, segmenti=96):
    bpy.ops.mesh.primitive_torus_add(major_radius=maggiore, minor_radius=minore,
                                     major_segments=segmenti, minor_segments=12,
                                     location=pos, rotation=rot)
    o = bpy.context.object
    o.name = nome
    o.scale = scala
    if materiale_:
        o.data.materials.append(materiale_)
    if genitore:
        o.parent = genitore
    return o


def sfera(nome, raggio, materiale_, pos=(0, 0, 0), genitore=None):
    bpy.ops.mesh.primitive_uv_sphere_add(radius=raggio, segments=64, ring_count=32,
                                         location=pos)
    o = bpy.context.object
    o.name = nome
    o.data.materials.append(materiale_)
    if genitore:
        o.parent = genitore
    bpy.ops.object.shade_smooth()
    return o


def cilindro(nome, raggio, altezza, materiale_, pos=(0, 0, 0), rot=(0, 0, 0)):
    bpy.ops.mesh.primitive_cylinder_add(radius=raggio, depth=altezza, vertices=96,
                                        location=pos, rotation=rot)
    o = bpy.context.object
    o.name = nome
    o.data.materials.append(materiale_)
    return o


def elica(nome, raggio_in, raggio_out, altezza, giri, materiale_, punti=400, spessore=0.02):
    """Fascio elicoidale: curva a spirale con bevel, per l'OAM."""
    curva = bpy.data.curves.new(nome, type="CURVE")
    curva.dimensions = "3D"
    curva.bevel_depth = spessore
    curva.bevel_resolution = 4
    sp = curva.splines.new("POLY")
    sp.points.add(punti - 1)
    for i in range(punti):
        t = i / (punti - 1)
        a = 2 * math.pi * giri * t
        r = raggio_in + (raggio_out - raggio_in) * t
        sp.points[i].co = (r * math.cos(a), r * math.sin(a), altezza * (t - 0.5), 1)
    o = bpy.data.objects.new(nome, curva)
    o.data.materials.append(materiale_)
    bpy.context.collection.objects.link(o)
    return o


# ---------------------------------------------------------------------------
# COSTRUZIONE
# ---------------------------------------------------------------------------
def costruisci_reticolo(m_rame, m_nucleo, m_beam, genitore=None, mostra_fascio=True):
    """Il reticolo della V5: tre orbite ellittiche, bordo, due lenti, nucleo."""
    objs = []
    for i, (a, b) in enumerate(SPIN_ORBITS):
        o = toro(f"orbita_{i}", 1.0, WIRE_ORB / max(a, 1e-6) * 1.0, m_rame,
                 scala=(a, b, 1.0), genitore=genitore)
        objs.append(o)
    objs.append(toro("bordo", R_DISC, WIRE_RIM / R_DISC, m_rame, genitore=genitore))
    for i in range(2):
        rot = (math.pi / 2 if i == 0 else 0.0, 0, 0)
        objs.append(toro(f"lente_{i}", 1.0, WIRE_LENS / LENS[0], m_rame,
                         scala=(LENS[0], LENS[1], 1.0),
                         rot=rot, genitore=genitore))
    objs.append(sfera("nucleo", CORE_R, m_nucleo, genitore=genitore))
    if mostra_fascio:
        # solo l'elica: legge come fascio con momento angolare orbitale
        objs.append(elica("elica_esterna", 0.90, 0.90, 16.0, 2.5, m_beam, spessore=0.045))
        objs.append(elica("elica_interna", 0.42, 0.42, 16.0, 3.5, m_beam, spessore=0.028))
    return objs


def scena_hero():
    pulisci()
    mondo()
    m_rame = materiale("rame", (0.80, 0.36, 0.15), 1.0, 0.30)
    m_nucleo = materiale("nucleo", (0.90, 0.45, 0.12), 1.0, 0.24,
                         emissione=(1.0, 0.42, 0.08), forza=0.55)
    m_beam = materiale("fascio", (0.15, 1.0, 0.55), 0.0, 0.3,
                       emissione=(0.10, 1.0, 0.45), forza=7.0)
    vuoto = bpy.data.objects.new("rotore", None)
    bpy.context.collection.objects.link(vuoto)
    costruisci_reticolo(m_rame, m_nucleo, m_beam, genitore=vuoto)

    # piano di fondo appena riflettente, per staccare l'oggetto
    bpy.ops.mesh.primitive_plane_add(size=90, location=(0, 0, -4.2))
    piano = bpy.context.object
    piano.data.materials.append(materiale("piano", (0.015, 0.018, 0.024), 0.35, 0.28))

    luce("key", (6, -9, 8), 1400, (1.0, 0.94, 0.86), 9)
    luce("fill", (-8, -6, 3), 500, (0.45, 0.70, 1.0), 10)
    luce("rim", (0, 9, 5), 900, (0.45, 1.0, 0.70), 8)
    cam = camera("cam", (9.5, -11.5, 7.2), (0, 0, -0.1), 68)
    cam.data.dof.use_dof = True
    cam.data.dof.focus_distance = 15.0
    cam.data.dof.aperture_fstop = 2.2

    imposta_render(os.path.join(OUT, "rotore_hero.png"), (1100, 820), 64)
    bpy.ops.render.render(write_still=True)
    print("  render: rotore_hero.png")

    # sequenza: tre orientazioni
    for k, ang in enumerate((0, 120, 240)):
        vuoto.rotation_euler = (math.radians(0), math.radians(0), math.radians(ang))
        cam2 = camera("cam_seq", (1.2, -12.6, 6.4), (0, 0, 0), 62)
        cam2.data.dof.use_dof = True
        cam2.data.dof.focus_distance = 12.8
        cam2.data.dof.aperture_fstop = 4.0
        imposta_render(os.path.join(OUT, f"rotore_seq_{k}.png"), (760, 760), 48)
        bpy.ops.render.render(write_still=True)
        print(f"  render: rotore_seq_{k}.png")


def testo(nome, corpo, loc, size=0.5, materiale_=None):
    """Etichetta su un piano verticale, rivolta alla camera."""
    cu = bpy.data.curves.new(nome, "FONT")
    cu.body = corpo
    cu.size = size
    cu.align_x = "CENTER"
    cu.align_y = "BOTTOM"
    o = bpy.data.objects.new(nome, cu)
    o.location = loc
    o.rotation_euler = (math.radians(90), 0, 0)
    if materiale_ is not None:
        o.data.materials.append(materiale_)
    bpy.context.collection.objects.link(o)
    return o


def scena_baseline():
    """Le sei geometrie, una per frame, con lo stesso ingombro in immagine."""
    pulisci()
    mondo()
    m_rame = materiale("rame", (0.80, 0.36, 0.15), 1.0, 0.28)
    m_silice = materiale("silice", (0.55, 0.80, 0.92), 0.0, 0.06)
    m_testo = materiale("testo", (0.80, 0.84, 0.90), 0.0, 0.5,
                        emissione=(0.78, 0.84, 0.95), forza=1.4)
    m_beam = materiale("fascio", (0.15, 1.0, 0.55), 0.0, 0.3,
                       emissione=(0.10, 1.0, 0.45), forza=6.0)

    # dimensioni in unita' di scena (1 u = 100 nm), dal calcolo a pari massa
    casi = [
        ("base_1", "reticolo", "LATTICE  /  RETICOLO\n740 nm  -  frame + lenses",
         lambda: costruisci_reticolo(m_rame, m_rame, m_beam, mostra_fascio=True)),
        ("base_2", "disco", "SOLID DISC  /  DISCO PIENO\n200 nm  -  aspect 10:1",
         lambda: cilindro("disco", 0.996, 0.100, m_rame)),
        ("base_3", "anello", "RING  /  ANELLO\n586 nm  -  a/R = 1/40",
         lambda: toro("anello", 2.931, 0.073, m_rame)),
        ("base_4", "sfera", "SPHERE  /  SFERA\n84 nm",
         lambda: sfera("sfera", 0.420, m_rame)),
        ("base_5", "manubrio", "NANODUMBBELL  /  NANOMANUBRIO\n200 nm  -  d = 4r",
         lambda: [sfera("s1", 0.333, m_rame, (-0.667, 0, 0)),
                  sfera("s2", 0.333, m_rame, (0.667, 0, 0))]),
        ("base_6", "silice", "SILICA  /  SILICE\n134 nm  -  anisotropic",
         lambda: sfera("silice", 0.671, m_silice)),
    ]
    luce("key", (5, -7, 7), 900, (1.0, 0.94, 0.86), 8)
    luce("fill", (-7, -5, 3), 320, (0.45, 0.70, 1.0), 9)
    luce("rim", (0, 7, 5), 600, (0.45, 1.0, 0.70), 7)

    for nome, tipo, etichetta, costruttore in casi:
        # pulizia degli oggetti della variante precedente (tranne luci e camera)
        for o in list(bpy.data.objects):
            if o.type in {"MESH", "CURVE", "EMPTY"}:
                bpy.data.objects.remove(o, do_unlink=True)
        costruttore()
        # fascio piu' corto nelle inquadrature singole, per non uscire dal frame
        for nm in ("elica_esterna", "elica_interna"):
            ob = bpy.data.objects.get(nm)
            if ob:
                ob.scale = (1, 1, 0.55)
        # piano di appoggio coerente con l'inquadratura
        dist = {"reticolo": 15.0, "disco": 5.4, "anello": 11.5, "sfera": 3.4,
                "manubrio": 6.0, "silice": 4.6}[tipo]
        z_piano = -{"reticolo": 2.6, "disco": 0.35, "anello": 3.2, "sfera": 0.75,
                    "manubrio": 1.1, "silice": 1.0}[tipo]
        bpy.ops.mesh.primitive_plane_add(size=dist * 6, location=(0, 0, z_piano))
        bpy.context.object.data.materials.append(
            materiale("piano", (0.015, 0.018, 0.024), 0.35, 0.30))
        # etichetta appoggiata al piano, davanti all'oggetto
        testo(f"lab_{tipo}", etichetta, (0.0, -0.105 * dist, z_piano + 0.02),
              0.026 * dist, m_testo)
        camera("cam", (0.575 * dist, -0.99 * dist, 0.46 * dist), (0, 0, 0.4 * z_piano), 62)
        imposta_render(os.path.join(OUT, f"{nome}.png"), (760, 760), 48)
        bpy.ops.render.render(write_still=True)
        print(f"  render: {nome}.png ({tipo})")


if __name__ == "__main__":
    print("BLENDER V6 - scena del rotore otticamente guidato")
    scena_hero()
    # salva la scena principale, cosi' si apre direttamente in Blender
    bpy.ops.wm.save_as_mainfile(filepath=os.path.join(OUT, "motore_fotonico_v6.blend"))
    print("  scena salvata: motore_fotonico_v6.blend")
    scena_baseline()
    print("fatto: render in", OUT)
