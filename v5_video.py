"""Render dell'animazione completa del sistema combinato V5 (240 fotogrammi).

Uso:
    blender -b motore_fotonico_v5.blend -P v5_video.py

Produce versione5/render/motore_fotonico_v5.mp4 (H.264, 24 fps, 1280x720).
Le chiavi di rotazione sono assolute (0 -> +720 gradi sullo Spin,
0 -> -86.1702 gradi sul Contro-spin), quindi il video mostra davvero
il doppio giro.
"""
import os
import time

import bpy

OUT_DIR = r"C:\Users\acese\Desktop\motore fotonico del vuoto\versione5"
RENDER_DIR = os.path.join(OUT_DIR, "render")
os.makedirs(RENDER_DIR, exist_ok=True)

scene = bpy.context.scene
scene.render.engine = "BLENDER_EEVEE"
scene.render.resolution_x = 1280
scene.render.resolution_y = 720
scene.render.resolution_percentage = 100
scene.render.fps = 24
scene.frame_start = 1
scene.frame_end = 240
scene.render.use_file_extension = True

scene.render.image_settings.media_type = "VIDEO"
scene.render.image_settings.file_format = "FFMPEG"
scene.render.ffmpeg.format = "MPEG4"
scene.render.ffmpeg.codec = "H264"
scene.render.ffmpeg.constant_rate_factor = "HIGH"
scene.render.ffmpeg.ffmpeg_preset = "GOOD"
scene.render.filepath = os.path.join(RENDER_DIR, "motore_fotonico_v5")

print("[video] scena      :", scene.name)
print("[video] oggetti    :", len(bpy.data.objects))
print("[video] fotogrammi :", scene.frame_start, "-", scene.frame_end)
print("[video] uscita     :", scene.render.filepath)

t0 = time.time()
bpy.ops.render.render(animation=True)
print("[video] completato in %.1f s" % (time.time() - t0))
for nome in ("motore_fotonico_v5.mp4", "motore_fotonico_v50001-0240.mp4"):
    percorso = os.path.join(RENDER_DIR, nome)
    print("[video] %-32s %s" % (nome, os.path.exists(percorso)))
