"""Render analitico dei fluidi V5: vista ortografica dall'alto del solo tracciante.

Produce una sequenza PNG (una ogni 4 fotogrammi) in fluidi/analisi, pensata per
l'estrazione delle curve: i rotori e gli iniettori sono nascosti, resta visibile
solo il dominio volumetrico, cosi' che il colore dei pixel sia solo tracciante
(arancio = Spin, blu = Contro-spin).

Uso:
    blender -b fluidi_v5.blend -P v5_fluidi_analisi.py
"""
import math
import os

import bpy

BASE = r"C:\Users\acese\Desktop\motore fotonico del vuoto\versione5\fluidi"
OUT = os.path.join(BASE, "analisi")
os.makedirs(OUT, exist_ok=True)

PASSO = 4

scene = bpy.context.scene
scene.render.engine = "BLENDER_EEVEE"
scene.render.resolution_x = 640
scene.render.resolution_y = 360
scene.render.resolution_percentage = 100
scene.render.image_settings.media_type = "IMAGE"
scene.render.image_settings.file_format = "PNG"
scene.render.film_transparent = False
scene.eevee.taa_render_samples = 24
if hasattr(scene.eevee, "volumetric_samples"):
    scene.eevee.volumetric_samples = 64
    scene.eevee.volumetric_start = 0.1
    scene.eevee.volumetric_end = 40.0

dominio = bpy.data.objects["Fluido_Dominio"]
for ob in bpy.context.view_layer.objects:
    if ob is not dominio:
        ob.hide_render = True

# volume piu' brillante: la misura e' fotometrica, serve margine dinamico
mat = dominio.data.materials[0]
vol = [n for n in mat.node_tree.nodes if n.type == "PRINCIPLED_VOLUME"][0]
vol.inputs["Emission Strength"].default_value = 4.0

bpy.ops.object.camera_add(location=(0.0, 0.0, 14.0))
cam = bpy.context.active_object
cam.name = "Analisi_Camera"
cam.data.type = "ORTHO"
cam.data.ortho_scale = 16.6
cam.rotation_euler = (0.0, 0.0, 0.0)
scene.camera = cam

print("[analisi] vista ortografica, scala %.2f, pixel per unita' %.3f"
      % (cam.data.ortho_scale, 640.0 / cam.data.ortho_scale))
for f in range(scene.frame_start, scene.frame_end + 1, PASSO):
    scene.frame_set(f)
    scene.render.filepath = os.path.join(OUT, "tracciante_f%03d.png" % f)
    bpy.ops.render.render(write_still=True)
    print("[analisi] scritto", os.path.basename(scene.render.filepath))
print("[analisi] fine")
