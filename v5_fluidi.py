"""
MOTORE FOTONICO V5 - STUDIO DEI FLUSSI DI ROTAZIONE (Mantaflow)
Autore: Dr. Francesco Bulla / Sistema GNRX

Dominio gassoso con g = 0 (spazio profondo, nessuna galleggiamento) e due
rotori affiancati animati con le stesse velocita' angolari della scena V5:
    Spin        +720.0000 gradi in 120 fotogrammi  ->  +2.5342 rad/s
    Contro-spin  -86.1702 gradi in 120 fotogrammi  ->  -0.3033 rad/s

Ogni rotore ha tre anelli concentrici complanari che ruotano come ostacoli
(Collision) e otto iniettori disposti sull'orbita esterna, ciascuno con la
velocita' tangenziale omega*r del proprio rotore: il tracciante entra gia'
co-rotante e il campo si sviluppa per pura convezione, senza gravita'.
Gli iniettori del Contro-spin hanno segno opposto, quindi nella zona centrale
i due flussi si incontrano in controcorrente.

Il bake scrive la cache in versione5/fluidi/cache e salva fluidi_v5.blend.
"""

import math
import os
import shutil

import bpy

OUT_DIR = r"C:\Users\acese\Desktop\motore fotonico del vuoto\versione5"
FLUID_DIR = os.path.join(OUT_DIR, "fluidi")
CACHE_DIR = os.path.join(FLUID_DIR, "cache")
BLEND_PATH = os.path.join(FLUID_DIR, "fluidi_v5.blend")
os.makedirs(FLUID_DIR, exist_ok=True)

FPS = 24
FRAME_START = 1
FRAME_END = 120
RES_MAX = 128
DURATA_S = (FRAME_END - FRAME_START) / float(FPS)

DOM_HALF = (8.0, 4.0, 1.6)          # semiestensioni del dominio (1 unita' = 100 nm)
SEPARAZIONE = 4.40
ANELLI = ((2.54, 2.30, 0.0), (2.34, 1.97, 62.0), (2.08, 1.73, 121.0))
SPIN_DEG = 720.0
CONTRO_DEG = -86.1702
R_INIETTORI = 2.90
N_INIETTORI = 8

RAME = (0.72, 0.45, 0.20, 1.0)
SCURO = (0.22, 0.10, 0.04, 1.0)
FUMO_SPIN = (0.95, 0.48, 0.12)
FUMO_CONTRO = (0.20, 0.55, 0.95)


# ---------------------------------------------------------------------------
# 1. Scena pulita
# ---------------------------------------------------------------------------
def nuova_scena():
    scene = bpy.data.scenes.new("Fluidi_V5")
    try:
        bpy.context.window.scene = scene
    except Exception:                      # noqa: BLE001
        pass
    for altra in list(bpy.data.scenes):
        if altra is not scene:
            bpy.data.scenes.remove(altra)
    for ob in list(bpy.data.objects):
        bpy.data.objects.remove(ob, do_unlink=True)
    for blocco in (bpy.data.meshes, bpy.data.materials, bpy.data.lights,
                   bpy.data.cameras):
        for item in list(blocco):
            if item.users == 0:
                blocco.remove(item)
    scene.frame_start = FRAME_START
    scene.frame_end = FRAME_END
    scene.render.fps = FPS
    scene.render.fps_base = 1.0
    scene.render.engine = "BLENDER_EEVEE"
    scene.render.resolution_x = 1280
    scene.render.resolution_y = 720
    world = bpy.data.worlds.new("Mondo_Fluidi")
    world.use_nodes = True
    world.node_tree.nodes["Background"].inputs[0].default_value = (
        0.008, 0.010, 0.016, 1.0)
    scene.world = world
    return scene


def materiale_semplice(nome, colore, metallic=0.85, rough=0.32):
    mat = bpy.data.materials.new(nome)
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes["Principled BSDF"]
    bsdf.inputs["Base Color"].default_value = colore
    bsdf.inputs["Metallic"].default_value = metallic
    bsdf.inputs["Roughness"].default_value = rough
    return mat


def materiale_fumo():
    """Volume: densita' presa dal grid 'density' del dominio Mantaflow."""
    mat = bpy.data.materials.new("Fumo_Dominio")
    mat.use_nodes = True
    nt = mat.node_tree
    nt.nodes.clear()
    out = nt.nodes.new("ShaderNodeOutputMaterial")
    out.location = (420, 0)
    vol = nt.nodes.new("ShaderNodeVolumePrincipled")
    vol.location = (120, 0)
    vol.inputs["Anisotropy"].default_value = 0.35
    vol.inputs["Color"].default_value = (0.92, 0.94, 1.0, 1.0)
    vol.inputs["Emission Strength"].default_value = 1.4
    attr = nt.nodes.new("ShaderNodeAttribute")
    attr.location = (-460, 0)
    attr.attribute_name = "density"
    attr.attribute_type = "GEOMETRY"
    scala = nt.nodes.new("ShaderNodeMath")
    scala.location = (-180, 0)
    scala.operation = "MULTIPLY"
    scala.inputs[1].default_value = 90.0   # il grid vale 0..1: senza questo il fumo e' evanescente
    color = nt.nodes.new("ShaderNodeAttribute")
    color.location = (-460, -220)
    color.attribute_name = "color"
    color.attribute_type = "GEOMETRY"
    nt.links.new(attr.outputs["Fac"], scala.inputs[0])
    nt.links.new(scala.outputs["Value"], vol.inputs["Density"])
    nt.links.new(color.outputs["Color"], vol.inputs["Color"])
    nt.links.new(color.outputs["Color"], vol.inputs["Emission Color"])
    nt.links.new(vol.outputs["Volume"], out.inputs["Volume"])
    return mat


# ---------------------------------------------------------------------------
# 2. Dominio
# ---------------------------------------------------------------------------
def costruisci_dominio(scene):
    bpy.ops.mesh.primitive_cube_add(size=2.0, location=(0.0, 0.0, 0.0))
    dom = bpy.context.active_object
    dom.name = "Fluido_Dominio"
    dom.scale = DOM_HALF
    bpy.ops.object.transform_apply(location=False, rotation=False, scale=True)

    mod = dom.modifiers.new("Fluid", "FLUID")
    mod.fluid_type = "DOMAIN"
    ds = mod.domain_settings
    ds.domain_type = "GAS"
    ds.resolution_max = RES_MAX
    ds.use_adaptive_timesteps = True
    ds.cfl_condition = 4.0
    ds.gravity = (0.0, 0.0, 0.0)     # g = 0, spazio profondo
    ds.beta = 0.0                    # nessuna galleggiamento termico
    ds.alpha = 0.10                  # dissipazione del tracciante (archi leggibili)
    ds.vorticity = 0.12              # confinamento della vorticita'
    ds.use_dissolve_smoke = False
    ds.time_scale = 1.0
    ds.cache_type = "ALL"
    ds.cache_directory = CACHE_DIR
    ds.cache_frame_start = FRAME_START
    ds.cache_frame_end = FRAME_END
    ds.use_collision_border_front = True
    ds.use_collision_border_back = True
    ds.use_collision_border_left = True
    ds.use_collision_border_right = True
    ds.use_collision_border_top = True
    ds.use_collision_border_bottom = True

    dom.data.materials.append(materiale_fumo())
    dom.display_type = "WIRE"
    print("[fluidi] dominio: %s  griglia %dx%dx%d"
          % (DOM_HALF, RES_MAX, int(RES_MAX * DOM_HALF[1] / DOM_HALF[0]),
             int(RES_MAX * DOM_HALF[2] / DOM_HALF[0])))
    return dom


# ---------------------------------------------------------------------------
# 3. Rotori (ostacoli in movimento)
# ---------------------------------------------------------------------------
def chiavi_rotazione(oggetto, angolo_totale, azimut_deg=0.0):
    """Rotazione assoluta, senza riduzione modulo 360, interpolazione lineare."""
    for frame, frazione in ((FRAME_START, 0.0), (FRAME_END, 1.0)):
        oggetto.rotation_euler.z = math.radians(azimut_deg
                                                + angolo_totale * frazione)
        oggetto.keyframe_insert("rotation_euler", index=2, frame=frame)


def costruisci_rotore(nome, x_centro, angolo_totale, colori):
    for i, (a, b, azimut) in enumerate(ANELLI):
        bpy.ops.mesh.primitive_torus_add(
            major_radius=1.0, minor_radius=0.14, major_segments=72,
            minor_segments=12, location=(x_centro, 0.0, 0.0),
            rotation=(0.0, 0.0, math.radians(azimut)))
        toro = bpy.context.active_object
        toro.name = "%s_Anello_%d" % (nome, i + 1)
        toro.scale = (a, b, 1.0)
        bpy.ops.object.transform_apply(location=False, rotation=False,
                                       scale=True)
        toro.data.materials.append(colori[i])
        collisione = toro.modifiers.new("Collision", "COLLISION")
        collisione.settings.damping = 1.0
        collisione.settings.thickness_outer = 0.06
        chiavi_rotazione(toro, angolo_totale, azimut)

    bpy.ops.mesh.primitive_uv_sphere_add(radius=0.34, segments=24,
                                         ring_count=16,
                                         location=(x_centro, 0.0, 0.0))
    nucleo = bpy.context.active_object
    nucleo.name = "%s_Nucleo" % nome
    nucleo.data.materials.append(colori[0])
    nucleo.modifiers.new("Collision", "COLLISION")


def costruisci_iniettori(nome, x_centro, omega, colore_fumo):
    """Otto iniettori sull'orbita esterna, solidali al rotore.

    Ogni iniettore orbita con la stessa omega del rotore (chiavi di posizione
    assolute) e immette con la velocita' tangenziale omega*r: il tracciante
    resta cosi' depositato lungo un'elica, che e' la firma visibile dello spin.
    """
    velocita = abs(omega) * R_INIETTORI
    for k in range(N_INIETTORI):
        phi0 = 2.0 * math.pi * k / N_INIETTORI
        segno = 1.0 if omega > 0.0 else -1.0

        x = x_centro + R_INIETTORI * math.cos(phi0)
        y = R_INIETTORI * math.sin(phi0)
        vx = -math.sin(phi0) * velocita * segno
        vy = math.cos(phi0) * velocita * segno

        bpy.ops.mesh.primitive_uv_sphere_add(radius=0.16, segments=16,
                                             ring_count=10,
                                             location=(x, y, 0.0))
        sorg = bpy.context.active_object
        sorg.name = "%s_Iniettore_%d" % (nome, k + 1)
        mod = sorg.modifiers.new("Fluid", "FLUID")
        mod.fluid_type = "FLOW"
        fs = mod.flow_settings
        fs.flow_type = "SMOKE"
        fs.flow_behavior = "INFLOW"
        fs.flow_source = "MESH"
        fs.density = 1.0
        fs.temperature = 0.0
        fs.smoke_color = colore_fumo
        fs.surface_distance = 0.20      # MAI 0.0: azzera l'emissione
        fs.use_initial_velocity = True
        fs.velocity_coord = (vx, vy, 0.0)
        fs.velocity_factor = 1.0
        sorg.display_type = "WIRE"

        for frame in range(FRAME_START, FRAME_END + 1, 6):
            t = (frame - FRAME_START) / float(FPS)
            phi = phi0 + omega * t
            sorg.location = (x_centro + R_INIETTORI * math.cos(phi),
                             R_INIETTORI * math.sin(phi), 0.0)
            sorg.keyframe_insert("location", frame=frame)
        sorg.location = (x_centro + R_INIETTORI * math.cos(phi0 + omega * DURATA_S),
                         R_INIETTORI * math.sin(phi0 + omega * DURATA_S), 0.0)
        sorg.keyframe_insert("location", frame=FRAME_END)
    print("[fluidi] %s: omega = %+.4f rad/s  velocita' tangenziale %.3f unita'/s"
          % (nome, omega, velocita))


# ---------------------------------------------------------------------------
# 4. Camera e luci
# ---------------------------------------------------------------------------
def costruisci_camera(scene):
    bpy.ops.object.camera_add(location=(0.0, -13.5, 13.5))
    cam = bpy.context.active_object
    cam.name = "Fluidi_Camera"
    cam.data.lens = 38.0
    cam.rotation_euler = (math.radians(45.0), 0.0, 0.0)
    scene.camera = cam

    bpy.ops.object.light_add(type="AREA", location=(-2.0, -9.0, 9.0))
    key = bpy.context.active_object
    key.data.energy = 3000.0
    key.data.size = 9.0
    key.rotation_euler = (math.radians(35.0), 0.0, math.radians(-12.0))

    bpy.ops.object.light_add(type="AREA", location=(7.0, 2.0, 6.0))
    rim = bpy.context.active_object
    rim.data.energy = 1800.0
    rim.data.size = 7.0
    rim.rotation_euler = (math.radians(60.0), 0.0, math.radians(115.0))
    return cam


def main():
    if os.path.isdir(CACHE_DIR):
        shutil.rmtree(CACHE_DIR)
    os.makedirs(CACHE_DIR, exist_ok=True)

    scene = nuova_scena()
    bpy.context.preferences.edit.keyframe_new_interpolation_type = "LINEAR"
    dominio = costruisci_dominio(scene)

    colori = (materiale_semplice("Anello_Scuro", SCURO),
              materiale_semplice("Anello_Chiaro", RAME, 0.9, 0.22),
              materiale_semplice("Anello_Chiaro2", RAME, 0.9, 0.22))

    omega_spin = math.radians(SPIN_DEG) / DURATA_S
    omega_contro = math.radians(CONTRO_DEG) / DURATA_S

    costruisci_rotore("Spin", -SEPARAZIONE, SPIN_DEG, colori)
    costruisci_rotore("Contro", +SEPARAZIONE, CONTRO_DEG, colori)
    costruisci_iniettori("Spin", -SEPARAZIONE, omega_spin, FUMO_SPIN)
    costruisci_iniettori("Contro", +SEPARAZIONE, omega_contro, FUMO_CONTRO)
    costruisci_camera(scene)

    bpy.ops.wm.save_as_mainfile(filepath=BLEND_PATH)
    print("[fluidi] scena salvata:", BLEND_PATH)
    print("[fluidi] oggetti:", len(bpy.data.objects))
    print("[fluidi] avvio bake", FRAME_START, "-", FRAME_END,
          "| risoluzione", RES_MAX)
    for ob in bpy.context.view_layer.objects:
        ob.select_set(False)
    dominio.select_set(True)
    bpy.context.view_layer.objects.active = dominio
    bpy.ops.fluid.bake_all()
    print("[fluidi] bake completato | dati:", dominio.modifiers["Fluid"]
          .domain_settings.has_cache_baked_data)
    bpy.ops.wm.save_as_mainfile(filepath=BLEND_PATH)


if __name__ == "__main__":
    main()
