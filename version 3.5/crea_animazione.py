"""Punto 1 degli sviluppi futuri del paper
"Modellazione geometrica 3D del sistema Spin / Contro-spin per il Motore Fotonico v3"

ANIMAZIONE DINAMICA COMPLETA DEL SISTEMA COMBINATO
==================================================
Assembla la sequenza renderizzata da Blender (240 frame, 24 fps, 1200x700) in:
  - un foglio di contatto a 6 riquadri  -> Figura_8_sequenza_animazione.png
  - un GIF animato a 24 fps             -> animazione_sistema_combinato.gif

I frame PNG provengono da versione5/frames/ (render EEVEE della scena "Scene",
camera Cam_Rame_DNA).  Il ciclo e' chiuso: frame 1 e frame 240 coincidono,
quindi il GIF puo' essere riprodotto in loop senza salti.
"""

import glob
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

OUT = os.path.dirname(os.path.abspath(__file__))
FRAMES = os.path.join(OUT, "frames")

files = sorted(glob.glob(os.path.join(FRAMES, "frame_*.png")))
if not files:
    raise SystemExit("Nessun frame trovato in %s" % FRAMES)
print("Frame trovati:", len(files))

# ----------------------------------------------------------------------------
# 1. Foglio di contatto (6 riquadri)
# ----------------------------------------------------------------------------
SCELTA = [
    (1, "frame 1 - sistemi distanti\nd = 950 nm"),
    (45, "frame 45 - avvicinamento\nd = 734 nm"),
    (75, "frame 75 - fase finale\nd = 543 nm"),
    (90, "frame 90 - contatto\nd = 510 nm"),
    (160, "frame 160 - contatto mantenuto\nrotazione dei due sistemi"),
    (240, "frame 240 - separazione\nd = 950 nm (ciclo chiuso)"),
]

fig, axs = plt.subplots(3, 2, figsize=(13.0, 12.5))
fig.suptitle(
    "Motore Fotonico v4 - Animazione dinamica del sistema combinato Spin / Contro-spin\n"
    "240 frame, 24 fps, ciclo chiuso (frame 1 = frame 240)",
    fontsize=13, fontweight="bold",
)
for ax, (f, titolo) in zip(axs.ravel(), SCELTA):
    img = Image.open(os.path.join(FRAMES, "frame_%04d.png" % f))
    ax.imshow(img)
    ax.axis("off")
    ax.set_title(titolo, fontsize=10)

plt.tight_layout(rect=(0, 0, 1, 0.95))
p_story = os.path.join(OUT, "Figura_8_sequenza_animazione.png")
plt.savefig(p_story, dpi=150)
plt.close(fig)
print("Figura 8:", p_story)

# ----------------------------------------------------------------------------
# 2. GIF animato
# ----------------------------------------------------------------------------
LARGHEZZA = 760
FPS_OUT = 24
imgs = []
for f in files:
    im = Image.open(f).convert("RGB")
    h = int(round(im.height * LARGHEZZA / im.width))
    imgs.append(im.resize((LARGHEZZA, h), Image.LANCZOS))

p_gif = os.path.join(OUT, "animazione_sistema_combinato.gif")
imgs[0].save(p_gif, save_all=True, append_images=imgs[1:],
             duration=int(round(1000.0 / FPS_OUT)), loop=0, optimize=True)
print("GIF    :", p_gif, "%.2f MB" % (os.path.getsize(p_gif) / 1e6))
