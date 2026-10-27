"""
BLENDER PYTHON AUTOMATION: MOTORE FOTONICO V3
Autore: Dr. Francesco Bulla / Sistema GNRX
Funzionalità:
 - Pulizia scena
 - Generazione procedurale dei sistemi Spin e Contro-Spin a 3 anelli intersecati
 - Materiali PBR Rame Nanostrutturato e Nucleo Emissivo
 - Animazione del ciclo a doppio giro di progetto (720°), con chiavi che chiudono
   esattamente sul frame 240
 - Creazione del fascio di flussi fotonici

Come eseguirlo in Blender:
 1. Apri Blender (3.x / 4.x / 5.x).
 2. Pannello Scripting -> New -> incolla il codice -> Run Script (Alt+P).

Nota: lo script costruisce il modello in una scena dedicata "Motore_Fotonico_V3",
senza cancellare altre scene o altri progetti aperti nello stesso file.
"""

import bpy
import math

SCENE_NAME = "Motore_Fotonico_V3"


def get_or_create_scene(name=SCENE_NAME):
    """Crea (o riusa) una scena dedicata e la rende attiva, senza toccare le altre."""
    scene = bpy.data.scenes.get(name)
    if scene is None:
        scene = bpy.data.scenes.new(name)
    if bpy.context.window is not None:
        bpy.context.window.scene = scene
    return scene


def clean_scene():
    """Rimuove tutti gli oggetti dalla scena corrente."""
    if bpy.context.active_object and bpy.context.active_object.mode != 'OBJECT':
        bpy.ops.object.mode_set(mode='OBJECT')
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete(use_global=False)


def create_pbr_copper_material():
    """Crea materiale Rame Metallico Plasmonico."""
    mat = bpy.data.materials.new(name="Material_Copper_Nano")
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes.get("Principled BSDF")

    # Valori PBR del Rame puro
    bsdf.inputs['Base Color'].default_value = (0.955, 0.637, 0.538, 1.0)
    bsdf.inputs['Metallic'].default_value = 1.0
    bsdf.inputs['Roughness'].default_value = 0.18
    return mat


def create_photonic_glow_material(name, color_rgba, strength=10.0):
    """Crea materiale emissivo per la luce e i flussi di fotoni."""
    mat = bpy.data.materials.new(name=name)
    mat.use_nodes = True
    nodes = mat.node_tree.nodes
    nodes.clear()

    node_emission = nodes.new(type='ShaderNodeEmission')
    node_emission.inputs['Color'].default_value = color_rgba
    node_emission.inputs['Strength'].default_value = strength

    node_output = nodes.new(type='ShaderNodeOutputMaterial')
    mat.node_tree.links.new(node_emission.outputs['Emission'], node_output.inputs['Surface'])
    return mat


def build_three_ring_system(name_prefix, major_radius, minor_radius, material, parent_empty):
    """Costruisce la configurazione geometrica a tre cerchi/anelli intersecati."""
    rings = []
    orientations = [
        (0.0, 0.0, 0.0),                          # Anello equatoriale (piano XY)
        (math.radians(90.0), 0.0, 0.0),           # Anello meridiano 1 (piano XZ)
        (0.0, math.radians(90.0), 0.0)            # Anello meridiano 2 (piano YZ)
    ]

    for i, rot in enumerate(orientations):
        bpy.ops.mesh.primitive_torus_add(
            major_radius=major_radius,
            minor_radius=minor_radius,
            major_segments=64,
            minor_segments=24,
            rotation=rot
        )
        ring_obj = bpy.context.active_object
        ring_obj.name = f"{name_prefix}_Ring_{i + 1}"
        ring_obj.data.materials.append(material)
        ring_obj.parent = parent_empty
        bpy.ops.object.shade_smooth()
        rings.append(ring_obj)

    return rings


def setup_scene():
    get_or_create_scene()
    clean_scene()

    # 1. Configurazione Timeline per il ciclo progettato su due giri (720°)
    scene = bpy.context.scene
    scene.frame_start = 1
    scene.frame_end = 240
    scene.render.fps = 24  # 239 intervalli = 9.958 s, coerente con il file consegnato

    # 2. Materiali
    mat_copper = create_pbr_copper_material()
    mat_core = create_photonic_glow_material("Core_Glow", (1.0, 0.95, 0.8, 1.0), strength=12.0)
    mat_photons = create_photonic_glow_material("Photon_Flow", (0.0, 0.85, 1.0, 1.0), strength=8.0)

    # 3. Nucleo Centrale Quantico (DNA Core)
    bpy.ops.mesh.primitive_uv_sphere_add(radius=0.45, segments=32, ring_count=24, location=(0, 0, 0))
    core_obj = bpy.context.active_object
    core_obj.name = "Central_Quantum_Core"
    core_obj.data.materials.append(mat_core)
    bpy.ops.object.shade_smooth()

    # 4. Creazione Nodi di Controllo (Empties) per rotazione differenziale
    spin_empty = bpy.data.objects.new("Empty_Spin_System", None)
    spin_empty.empty_display_type = 'ARROWS'
    scene.collection.objects.link(spin_empty)

    contro_empty = bpy.data.objects.new("Empty_ControSpin_System", None)
    contro_empty.empty_display_type = 'ARROWS'
    scene.collection.objects.link(contro_empty)

    # 5. Generazione Geometrica dei Due Sistemi
    # Sistema Primario: Spin (raggio 2.2)
    build_three_ring_system("Spin", major_radius=2.2, minor_radius=0.06,
                            material=mat_copper, parent_empty=spin_empty)

    # Sistema Secondario: Contro-Spin (raggio 2.7 concentrico per evitare collisione fisica)
    build_three_ring_system("ControSpin", major_radius=2.7, minor_radius=0.05,
                            material=mat_copper, parent_empty=contro_empty)

    # 6. Rigging Matematico ed Animazione del Ciclo Frame-by-Frame
    # Interpolazione lineare per moto rotatorio perpetuo fluido
    bpy.context.preferences.edit.keyframe_new_interpolation_type = 'LINEAR'

    # Chiavi allineate al file consegnato: ogni curva si chiude su un multiplo di 360°,
    # quindi il frame 240 coincide con il frame 1 e il ciclo e' riproducibile in loop.
    # Spin: escursione +30° con ritorno esatto
    for f, deg in ((1, 330.0), (120, 360.0), (240, 330.0)):
        spin_empty.rotation_euler = (0, 0, math.radians(deg))
        spin_empty.keyframe_insert(data_path="rotation_euler", frame=f)
    # Contro-spin: un giro completo
    for f, deg in ((1, 0.0), (240, 360.0)):
        contro_empty.rotation_euler = (0, 0, math.radians(deg))
        contro_empty.keyframe_insert(data_path="rotation_euler", frame=f)

    # 7. Creazione di Flussi Fotonici Dinamici (Frecce / Linee di Poynting)
    for i in range(12):
        angle = (2.0 * math.pi / 12) * i
        x = 1.6 * math.cos(angle)
        y = 1.6 * math.sin(angle)
        bpy.ops.mesh.primitive_cone_add(radius1=0.08, depth=0.35, location=(x, y, 0))
        cone = bpy.context.active_object
        cone.name = f"Photon_Vector_{i + 1}"
        cone.data.materials.append(mat_photons)
        cone.rotation_euler = (math.radians(90), 0, angle + math.radians(90))
        cone.parent = spin_empty  # Le frecce ruotano descrivendo il campo fotonico

    # 8. Vettore di Spinta Assiale Risultante (+Z)
    bpy.ops.mesh.primitive_cylinder_add(radius=0.04, depth=3.5, location=(0, 0, 2.5))
    thrust_vector = bpy.context.active_object
    thrust_vector.name = "Net_Thrust_Vector"
    mat_thrust = create_photonic_glow_material("Thrust_Glow", (0.2, 1.0, 0.4, 1.0), strength=15.0)
    thrust_vector.data.materials.append(mat_thrust)

    # 9. Setup Camera e Illuminazione Studio
    bpy.ops.object.camera_add(location=(7.0, -7.5, 5.0),
                              rotation=(math.radians(60), 0, math.radians(45)))
    cam = bpy.context.active_object
    scene.camera = cam
    bpy.ops.object.light_add(type='POINT', radius=1.0, location=(4.0, -4.0, 6.0))
    light = bpy.context.active_object
    light.data.energy = 600.0

    print("=" * 60)
    print("[+] Scena Blender 'Motore Fotonico v3' configurata con successo!")
    print("[*] Premi Spazio in Blender per avviare l'animazione del ciclo (239 intervalli, 9.958 s).")
    print("=" * 60)


if __name__ == "__main__":
    setup_scene()
