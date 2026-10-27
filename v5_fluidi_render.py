"""Render di controllo della simulazione fluidi V5 (cache gia' cotta).

Uso:
    blender -b fluidi_v5.blend -P v5_fluidi_render.py -- [--video]
"""
import os
import sys

import bpy

BASE = r"C:\Users\acese\Desktop\motore fotonico del vuoto\versione5\fluidi"
RENDER_DIR = os.path.join(BASE, "render")
os.makedirs(RENDER_DIR, exist_ok=True)

scene = bpy.context.scene
scene.render.engine = "BLENDER_EEVEE"
scene.render.resolution_x = 1280
scene.render.resolution_y = 720
scene.render.resolution_percentage = 100
scene.render.film_transparent = False
scene.eevee.taa_render_samples = 48
if hasattr(scene.eevee, "volumetric_samples"):
    scene.eevee.volumetric_samples = 96
    scene.eevee.volumetric_start = 0.1
    scene.eevee.volumetric_end = 60.0

dominio = bpy.data.objects["Fluido_Dominio"]
ds = dominio.modifiers["Fluid"].domain_settings
print("[render] cache:", ds.cache_type, ds.cache_directory)
print("[render] dati cotti:", ds.has_cache_baked_data,
      "| mesh:", ds.has_cache_baked_mesh)

video = "--video" in sys.argv
if video:
    scene.frame_start = 1
    scene.frame_end = 120
    scene.render.image_settings.media_type = "VIDEO"
    scene.render.image_settings.file_format = "FFMPEG"
    scene.render.ffmpeg.format = "MPEG4"
    scene.render.ffmpeg.codec = "H264"
    scene.render.ffmpeg.constant_rate_factor = "HIGH"
    scene.render.filepath = os.path.join(RENDER_DIR, "fluidi_v5")
    bpy.ops.render.render(animation=True)
    print("[render] video completato")
else:
    scene.render.image_settings.media_type = "IMAGE"
    scene.render.image_settings.file_format = "PNG"
    for f in (12, 36, 60, 84, 108, 120):
        scene.frame_set(f)
        scene.render.filepath = os.path.join(RENDER_DIR, "flusso_f%03d.png" % f)
        bpy.ops.render.render(write_still=True)
        print("[render] scritto", scene.render.filepath)
print("[render] fine")
