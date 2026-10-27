"""
MOTORE FOTONICO V5 - COSTRUZIONE DELLA SCENA BLENDER (da zero)
Autore: Dr. Francesco Bulla / Sistema GNRX

Eseguito in live sull'istanza di Blender in ascolto sulla porta 9877.
Ricostruisce il progetto in una scena nuova, pulita, senza riusare nulla
del file V4.

Geometria (1 unita' di scena = 100 nm), fedele al disegno di riferimento:
  - ogni rotore e' composto da tre anelli complanari (il piu' esterno in
    rame brunito, gli altri in oro pallido), un anello minore attorno al
    nucleo, la sfera di rame al centro e i marcatori blu (sferetta +
    freccia) distribuiti sugli anelli, orientati lungo la tangente, cioe'
    lungo il verso di rotazione.
  - Spin       : ruota di +720 gradi (due giri esatti) in 240 fotogrammi.
  - Contro-spin: ruota di -86.1702 gradi nello stesso ciclo; il rapporto
                 45/376 e' rispettato e la chiusura esatta si ottiene dopo
                 188 cicli.
  - I due rotori sono affiancati lungo x, ciascuno sul proprio asse:
    Spin a x = -SEPARAZIONE_X, Contro-spin a x = +SEPARAZIONE_X.  Ognuno
    ruota attorno al proprio centro, non attorno a un centro comune.

Nessuna chiave di rotazione viene ridotta modulo 360: le curve restano
assolute, cosi' l'animazione compie davvero il doppio giro.
"""

import math
import os

import bpy
from mathutils import Euler, Vector

# ---------------------------------------------------------------------------
# 0. Parametri condivisi con v5_nucleo.py
# ---------------------------------------------------------------------------
OUT_DIR = r"C:\Users\acese\Desktop\motore fotonico del vuoto\versione5"
BLEND_PATH = os.path.join(OUT_DIR, "motore_fotonico_v5.blend")

SCENE_NAME = "Motore_Fotonico_V5"
FPS = 24
FRAME_START = 1
FRAME_END = 240

SPIN_DEG_CYCLE = 720.0         # due giri esatti
CONTRO_DEG_CYCLE = -86.1702    # 2 giri * 45/376 * 360, segno opposto
RATIO_TXT = "45/376"

# Anelli del rotore, tutti complanari (piano XY), come nel disegno:
# (semiasse x, semiasse y, rotazione in pianta in gradi, e' l'anello scuro?)
ANELLI = (
    (2.54, 2.30, 0.0, True),
    (2.34, 1.97, 62.0, False),
    (2.08, 1.73, 121.0, False),
    (0.64, 0.54, 22.0, False),
)
CORE_R = 0.34
WIRE_R = 0.030
WIRE_SOTTILE = 0.024

# Marcatori blu: per ogni anello, gli angoli (gradi) lungo l'ellisse.
MARK_ANGOLI = (
    (25.0, 150.0, 265.0),
    (80.0, 205.0, 320.0),
    (40.0, 168.0),
    (200.0,),
)
MARK_SFERA_R = 0.105
MARK_CONO_R = 0.085
MARK_CONO_H = 0.360

SEPARAZIONE_X = 4.40           # distanza di ogni rotore dall'origine (u.a.)
X_SPIN = -SEPARAZIONE_X
X_CONTRO = +SEPARAZIONE_X


# ---------------------------------------------------------------------------
# 1. Utilita' di scena
# ---------------------------------------------------------------------------
def nuova_scena():
    """Crea una scena vuota e la rende attiva, eliminando tutto il resto."""
    scene = bpy.data.scenes.new(SCENE_NAME)
    if bpy.context.window is not None:
        bpy.context.window.scene = scene
    for altra in list(bpy.data.scenes):
        if altra is not scene:
            bpy.data.scenes.remove(altra)
    scene.name = SCENE_NAME

    for ob in list(bpy.data.objects):
        bpy.data.objects.remove(ob, do_unlink=True)
    for blocco in (bpy.data.meshes, bpy.data.curves, bpy.data.materials,
                   bpy.data.lights, bpy.data.cameras):
        for item in list(blocco):
            if item.users == 0:
                blocco.remove(item)

    scene.frame_start = FRAME_START
    scene.frame_end = FRAME_END
    scene.render.fps = FPS
    scene.render.fps_base = 1.0
    return scene


def materiale_metallo(nome, base, ruvidita, emissione, forza):
    """Metallo PBR con una debole emissione: senza, con un mondo scuro,
    i fili sottili riflettono solo il nero e diventano invisibili."""
    mat = bpy.data.materials.new(nome)
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes.get("Principled BSDF")
    bsdf.inputs["Base Color"].default_value = base
    bsdf.inputs["Metallic"].default_value = 1.0
    bsdf.inputs["Roughness"].default_value = ruvidita
    bsdf.inputs["Emission Color"].default_value = emissione
    bsdf.inputs["Emission Strength"].default_value = forza
    return mat


def materiale_rame():
    """Sfera centrale: rame lucido, come nel disegno."""
    return materiale_metallo("V5_Rame_Plasmonico",
                             (0.88, 0.66, 0.46, 1.0), 0.22,
                             (1.0, 0.62, 0.32, 1.0), 0.35)


def materiale_anello_scuro():
    """Anello orizzontale: rame brunito, piu' scuro e opaco."""
    return materiale_metallo("V5_Rame_Scuro",
                             (0.22, 0.095, 0.035, 1.0), 0.34,
                             (1.0, 0.42, 0.14, 1.0), 0.06)


def materiale_anello_chiaro():
    """Anelli verticali: oro pallido, quasi satinato."""
    return materiale_metallo("V5_Oro_Pallido",
                             (0.88, 0.76, 0.56, 1.0), 0.26,
                             (1.0, 0.84, 0.56, 1.0), 0.20)


def materiale_freccia():
    """Marcatori: blu acciaio opaco, come nel disegno."""
    mat = bpy.data.materials.new("V5_Blu_Marcatore")
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes.get("Principled BSDF")
    bsdf.inputs["Base Color"].default_value = (0.36, 0.47, 0.68, 1.0)
    bsdf.inputs["Metallic"].default_value = 0.30
    bsdf.inputs["Roughness"].default_value = 0.42
    bsdf.inputs["Emission Color"].default_value = (0.30, 0.45, 0.80, 1.0)
    bsdf.inputs["Emission Strength"].default_value = 0.18
    return mat


def materiale_emissivo(nome, colore, forza):
    mat = bpy.data.materials.new(nome)
    mat.use_nodes = True
    nodi = mat.node_tree.nodes
    nodi.clear()
    em = nodi.new(type="ShaderNodeEmission")
    em.inputs["Color"].default_value = colore
    em.inputs["Strength"].default_value = forza
    out = nodi.new(type="ShaderNodeOutputMaterial")
    mat.node_tree.links.new(em.outputs["Emission"], out.inputs["Surface"])
    return mat


def nuovo_empty(nome):
    ob = bpy.data.objects.new(nome, None)
    ob.empty_display_type = "ARROWS"
    ob.empty_display_size = 1.2
    bpy.context.scene.collection.objects.link(ob)
    return ob


def toro(nome, maggiore, minore, materiale, genitore, scala=(1, 1, 1),
         posizione=(0, 0, 0), rotazione=(0, 0, 0), segmenti=128):
    bpy.ops.mesh.primitive_torus_add(
        major_radius=maggiore, minor_radius=minore,
        major_segments=segmenti, minor_segments=12,
        location=posizione, rotation=rotazione,
    )
    ob = bpy.context.active_object
    ob.name = nome
    ob.scale = scala
    ob.data.materials.append(materiale)
    bpy.ops.object.shade_smooth()
    ob.parent = genitore
    return ob


def marcatore(nome, punto, tangente, materiale, genitore):
    """Sferetta + freccia orientata lungo la tangente dell'anello: indica
    il verso di rotazione locale, come nel disegno di riferimento."""
    bpy.ops.mesh.primitive_uv_sphere_add(radius=MARK_SFERA_R, segments=20,
                                         ring_count=12, location=punto)
    sfera = bpy.context.active_object
    sfera.name = nome + "_Sfera"
    bpy.ops.object.shade_smooth()
    sfera.data.materials.append(materiale)
    sfera.parent = genitore

    bpy.ops.mesh.primitive_cone_add(radius1=MARK_CONO_R, radius2=0.0,
                                    depth=MARK_CONO_H, location=punto)
    cono = bpy.context.active_object
    cono.name = nome + "_Freccia"
    cono.rotation_mode = "QUATERNION"
    cono.rotation_quaternion = tangente.to_track_quat("Z", "Y")
    cono.data.materials.append(materiale)
    bpy.ops.object.shade_smooth()
    cono.parent = genitore


def punto_ellisse(semiassi, angolo, matrice=None):
    """Punto e tangente di un'ellisse percorsa in senso antiorario, poi
    trasformati dalla matrice di rotazione in pianta dell'anello."""
    a, b = semiassi
    c, s = math.cos(angolo), math.sin(angolo)
    punto = Vector((a * c, b * s, 0.0))
    tangente = Vector((-a * s, b * c, 0.0)).normalized()
    if matrice is not None:
        punto = matrice @ punto
        tangente = (matrice @ tangente).normalized()
    return punto, tangente


def costruisci_rotore(prefisso, genitore, mat_scuro, mat_chiaro, mat_rame,
                      mat_freccia):
    """Un rotore completo, come nel disegno: tre anelli complanari (il piu'
    esterno scuro, gli altri chiari), un anello minore attorno al nucleo,
    la sfera di rame al centro e i marcatori blu lungo gli anelli."""
    matrici = []
    for i, (a, b, azimut, scuro) in enumerate(ANELLI):
        rotazione = (0.0, 0.0, math.radians(azimut))
        toro("%s_Anello_%d" % (prefisso, i + 1),
             maggiore=1.0, minore=WIRE_R if i == 0 else WIRE_SOTTILE,
             materiale=mat_scuro if scuro else mat_chiaro,
             genitore=genitore,
             scala=(a, b, 1.0), rotazione=rotazione)
        matrici.append(Euler(rotazione, "XYZ").to_matrix())

    # --- sfera di rame al centro ------------------------------------------
    bpy.ops.mesh.primitive_uv_sphere_add(radius=CORE_R, segments=48,
                                         ring_count=32, location=(0, 0, 0))
    nucleo = bpy.context.active_object
    nucleo.name = "%s_Nucleo_Cu" % prefisso
    nucleo.data.materials.append(mat_rame)
    bpy.ops.object.shade_smooth()
    nucleo.parent = genitore

    # --- marcatori blu -----------------------------------------------------
    for i, angoli in enumerate(MARK_ANGOLI):
        for k, gradi in enumerate(angoli):
            p, t = punto_ellisse(ANELLI[i][:2], math.radians(gradi), matrici[i])
            marcatore("%s_Marca_%d_%02d" % (prefisso, i + 1, k), p, t,
                      mat_freccia, genitore)


def testo_sovraimpresso(nome, corpo, camera):
    """Didascalia 3D solidale alla camera (dimensione calibrata sul campo
    inquadrato a 1.60 unita' di distanza con obiettivo 58 mm)."""
    testo = bpy.data.curves.new(type="FONT", name=nome)
    testo.body = corpo
    testo.size = 0.022
    testo.align_x = "CENTER"
    testo.align_y = "CENTER"
    ob = bpy.data.objects.new(nome, testo)
    bpy.context.scene.collection.objects.link(ob)
    ob.data.materials.append(bpy.data.materials["V5_Ologramma"])
    ob.parent = camera
    ob.location = (0.0, -0.360, -1.60)
    return ob


def etichetta_3d(nome, corpo, posizione, camera, materiale, dimensione=0.34):
    """Etichetta piatta parallela al piano immagine, sopra un rotore."""
    testo = bpy.data.curves.new(type="FONT", name=nome)
    testo.body = corpo
    testo.size = dimensione
    testo.align_x = "CENTER"
    testo.align_y = "CENTER"
    ob = bpy.data.objects.new(nome, testo)
    bpy.context.scene.collection.objects.link(ob)
    ob.data.materials.append(materiale)
    ob.location = posizione
    ob.rotation_euler = camera.rotation_euler.copy()
    return ob


# ---------------------------------------------------------------------------
# 2. Costruzione
# ---------------------------------------------------------------------------
def costruisci():
    scene = nuova_scena()

    mat_rame = materiale_rame()
    mat_scuro = materiale_anello_scuro()
    mat_chiaro = materiale_anello_chiaro()
    mat_freccia = materiale_freccia()
    materiale_emissivo("V5_Ologramma", (0.35, 0.95, 1.0, 1.0), 6.0)

    asse_spin = nuovo_empty("V5_Spin_AX1")
    asse_cont = nuovo_empty("V5_Contro_AX2")
    asse_spin.location = (X_SPIN, 0.0, 0.0)
    asse_cont.location = (X_CONTRO, 0.0, 0.0)

    # --- I due motori, separati, ognuno sul proprio asse --------------------
    costruisci_rotore("V5_Spin", asse_spin, mat_scuro, mat_chiaro,
                      mat_rame, mat_freccia)
    costruisci_rotore("V5_Contro", asse_cont, mat_scuro, mat_chiaro,
                      mat_rame, mat_freccia)

    # --- Animazione: chiavi assolute, mai ridotte modulo 360 ----------------
    bpy.context.preferences.edit.keyframe_new_interpolation_type = "LINEAR"
    for ob in (asse_spin, asse_cont):
        ob.rotation_mode = "XYZ"

    fotogrammi = (1, 60, 120, 180, 240)
    for k, f in enumerate(fotogrammi):
        frazione = k / 4.0
        asse_spin.rotation_euler = (0.0, 0.0, math.radians(SPIN_DEG_CYCLE * frazione))
        asse_spin.keyframe_insert(data_path="rotation_euler", index=2, frame=f)
        asse_cont.rotation_euler = (0.0, 0.0, math.radians(CONTRO_DEG_CYCLE * frazione))
        asse_cont.keyframe_insert(data_path="rotation_euler", index=2, frame=f)

    # --- Camera, luci, mondo -------------------------------------------------
    bpy.ops.object.camera_add(location=(1.6, -15.5, 6.0))
    camera = bpy.context.active_object
    camera.name = "V5_Camera"
    camera.data.lens = 38.0
    direzione = Vector((0.0, 0.0, 0.10)) - camera.location
    camera.rotation_euler = direzione.to_track_quat("-Z", "Y").to_euler()
    scene.camera = camera

    for nome, pos, energia, colore in (
        ("V5_Key", (5.5, -6.5, 8.0), 2800.0, (1.0, 0.96, 0.90)),
        ("V5_Key2", (-6.5, -6.0, 8.0), 1800.0, (1.0, 0.90, 0.82)),
        ("V5_Fill", (-7.0, -4.5, 3.5), 1000.0, (0.55, 0.78, 1.0)),
        ("V5_Rim", (0.0, 8.5, 5.0), 1800.0, (0.35, 0.95, 1.0)),
    ):
        bpy.ops.object.light_add(type="AREA", location=pos)
        luce = bpy.context.active_object
        luce.name = nome
        luce.data.energy = energia
        luce.data.color = colore
        luce.data.size = 6.0
        luce.rotation_euler = (Vector((0, 0, 0)) - luce.location).to_track_quat("-Z", "Y").to_euler()

    mondo = bpy.data.worlds.new("V5_Mondo")
    mondo.use_nodes = True
    sfondo = mondo.node_tree.nodes["Background"]
    sfondo.inputs["Color"].default_value = (0.020, 0.032, 0.055, 1.0)
    sfondo.inputs["Strength"].default_value = 1.4
    scene.world = mondo

    # --- Resa ----------------------------------------------------------------
    scene.render.engine = "BLENDER_EEVEE"
    scene.eevee.taa_render_samples = 64
    if hasattr(scene.eevee, "use_raytracing"):
        scene.eevee.use_raytracing = True
    scene.render.resolution_x = 1280
    scene.render.resolution_y = 720
    scene.render.resolution_percentage = 100
    scene.render.image_settings.file_format = "PNG"
    scene.view_settings.view_transform = "AgX"
    try:
        scene.view_settings.look = "AgX - Punchy"
    except TypeError:
        pass

    testo_sovraimpresso(
        "V5_Didascalia",
        "MOTORE FOTONICO V5   ·   Spin 720° (due giri)   |   "
        "Contro-spin \u221286.1702° per ciclo   ·   %s   ·   rotori separati"
        % RATIO_TXT,
        camera)

    etichetta_3d("V5_Etichetta_Spin", "SPIN  ·  720°",
                 (X_SPIN, 0.0, 2.75), camera,
                 materiale_emissivo("V5_Txt_Spin", (1.0, 0.62, 0.28, 1.0), 8.0))
    etichetta_3d("V5_Etichetta_Contro", "CONTRO-SPIN  ·  \u221286.1702°",
                 (X_CONTRO, 0.0, 2.75), camera,
                 materiale_emissivo("V5_Txt_Contro", (0.35, 0.95, 1.0, 1.0), 8.0))

    # --- Metadati di scena ---------------------------------------------------
    scene["v5_versione"] = "V5"
    scene["v5_autore"] = "Dr. Francesco Bulla"
    scene["v5_spin_gradi_ciclo"] = SPIN_DEG_CYCLE
    scene["v5_contro_gradi_ciclo"] = CONTRO_DEG_CYCLE
    scene["v5_rapporto"] = RATIO_TXT
    scene["v5_cicli_chiusura"] = 188
    scene["v5_fotogrammi_chiusura"] = 188 * (FRAME_END - FRAME_START + 1)
    scene["v5_scala_unita"] = "1 unita' di scena = 100 nm"
    scene["v5_separazione_rotori"] = SEPARAZIONE_X
    scene["v5_spin_x"] = X_SPIN
    scene["v5_contro_x"] = X_CONTRO

    print("[V5] oggetti:", len(bpy.data.objects))
    print("[V5] frame:", scene.frame_start, scene.frame_end, "@", scene.render.fps, "fps")
    print("[V5] Spin  : 0 -> %.1f gradi" % SPIN_DEG_CYCLE)
    print("[V5] Contro: 0 -> %.4f gradi" % CONTRO_DEG_CYCLE)
    return scene


if __name__ == "__main__":
    costruisci()
    bpy.ops.wm.save_as_mainfile(filepath=BLEND_PATH)
    print("[V5] salvato:", BLEND_PATH)
