"""Figura 7 - Studio qualitativo dei flussi fotonici e vettore di spinta.

Modello geometrico illustrativo coerente con il paper "Modellazione geometrica 3D
del sistema Spin / Contro-spin per il Motore Fotonico v4": non e' una validazione
fisica, ma la visualizzazione dei flussi e del vettore di spinta risultante.

Output: Figura_7_studio_flussi_spinta.png nella cartella del progetto.
"""

import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = r"C:\Users\acese\Desktop\motore fotonico del vuoto\versione5"

rng = np.random.default_rng(7)
t = np.linspace(0.0, 1.0, 220)

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection="3d")

u = np.linspace(0, 2 * np.pi, 160)

# --- Sistema Spin: 3 anelli intersecati, raggio 1.0 (rame) ---
spin_rot = [(0.0, 0.0, 0.0), (np.pi / 2, 0.0, 0.0), (0.0, np.pi / 2, 0.0)]
for rx, ry, rz in spin_rot:
    Rx = np.array([[1, 0, 0], [0, np.cos(rx), -np.sin(rx)], [0, np.sin(rx), np.cos(rx)]])
    Ry = np.array([[np.cos(ry), 0, np.sin(ry)], [0, 1, 0], [-np.sin(ry), 0, np.cos(ry)]])
    Rz = np.array([[np.cos(rz), -np.sin(rz), 0], [np.sin(rz), np.cos(rz), 0], [0, 0, 1]])
    pts = np.vstack((np.cos(u), np.sin(u), np.zeros_like(u)))
    pts = (Rz @ Ry @ Rx) @ pts
    ax.plot(pts[0], pts[1], pts[2], color="#e67e22", lw=2.4, alpha=0.95)

# --- Sistema Contro-spin: 3 anelli intersecati, raggio 1.25 (giallo) ---
for rx, ry, rz in spin_rot:
    Rx = np.array([[1, 0, 0], [0, np.cos(rx), -np.sin(rx)], [0, np.sin(rx), np.cos(rx)]])
    Ry = np.array([[np.cos(ry), 0, np.sin(ry)], [0, 1, 0], [-np.sin(ry), 0, np.cos(ry)]])
    Rz = np.array([[np.cos(rz), -np.sin(rz), 0], [np.sin(rz), np.cos(rz), 0], [0, 0, 1]])
    pts = np.vstack((1.25 * np.cos(u), 1.25 * np.sin(u), np.zeros_like(u)))
    pts = (Rz @ Ry @ Rx) @ pts
    ax.plot(pts[0], pts[1], pts[2], color="#f1c40f", lw=1.6, alpha=0.65)

# --- Nucleo centrale ---
ax.scatter([0], [0], [0], color="#ffffff", s=110, edgecolors="black", zorder=5)

# --- Flussi fotonici: spirali che entrano dal basso e vengono canalizzate in +Z ---
for k in range(18):
    ang0 = 2 * np.pi * k / 18
    r0 = 1.9 + 0.25 * rng.random()
    z0 = -2.0
    swirl = 2.2 * t
    rad = r0 * (1.0 - 0.78 * t) + 0.15
    x = rad * np.cos(ang0 + swirl)
    y = rad * np.sin(ang0 + swirl)
    z = z0 + (t ** 1.6) * 4.6
    ax.plot(x, y, z, color="#00d4ff", lw=0.9, alpha=0.55)

# --- Vettore di spinta risultante (+Z) ---
ax.quiver(0, 0, 0.15, 0, 0, 2.3, color="#2ecc71", lw=4.0,
          arrow_length_ratio=0.16)

ax.text(0.05, 0.05, 2.55, "F$_z$ (spinta risultante)", color="#2ecc71", fontsize=10)
ax.text(1.45, -0.35, 0.35, "Contro-spin", color="#b7950b", fontsize=10)
ax.text(-1.75, 0.55, -0.35, "Spin 720", color="#a04000", fontsize=10)

ax.set_title(
    "Motore Fotonico v4 - Flussi fotonici e vettore di spinta (modello geometrico illustrativo)",
    fontsize=11,
)
ax.set_xlabel("X (u.a.)")
ax.set_ylabel("Y (u.a.)")
ax.set_zlabel("Z - asse di propulsione (u.a.)")
ax.set_xlim(-2.2, 2.2)
ax.set_ylim(-2.2, 2.2)
ax.set_zlim(-2.2, 3.0)
ax.view_init(elev=24, azim=-58)
ax.grid(alpha=0.25)

plt.tight_layout()
path = os.path.join(OUT, "Figura_7_studio_flussi_spinta.png")
plt.savefig(path, dpi=220)
plt.close(fig)

# --- Stima d'ordine di grandezza, dichiarata come tale ---
P = 1.0e-6
c = 299792458.0
eta_spin, R_cont, alpha = 0.70, 0.85, np.radians(15.0)
F = (P / c) * (eta_spin + R_cont * np.cos(alpha))
print("Figura 7:", path)
print(f"Stima illustrativa F_z = {F * 1e15:.3f} fN per 1 uW (non e' una validazione fisica)")
