# Motore Fotonico del Vuoto — V7.2

Appendice di revisione della V7.1. La V7.1 resta intatta in
`../versione 7.1/`; questa cartella contiene i tre passi richiesti, le
correzioni e i dati che li sostengono.

## I tre passi

1. **1064 nm nella griglia.** La griglia passa da 532/800 a
   532/800/1064 nm. La profondità balistica passa da 22.0 a 28.4 µm. Il buco
   termico a 2.34 µm (800 nm) si sposta a 3.11 µm e scende da 27.2 a 15.7 K,
   sotto il limite di 20 K. Il bordo a basso R resta termico: a
   R = 0.3 µm `P_target` cresce con λ (644 → 836 µW), quindi 1064 nm non
   sposta la soglia di potenza.
2. **Costo energetico del contrappeso.** I cresce di 6.31, β di 3.88.
   A pari Pe_r, `P_target` non dipende da β: la potenza assorbita è la
   stessa. Il ciclo però dura di più. Con il latch a metà ciclo:
   `T_720 = T_360_no + T_360_cw` → 2.44× un 720° senza contrappeso e 4.88×
   un 360° senza. Il 3.88 della V7.1 è il limite conservativo con
   contrappeso sempre innestato (7.76× un 360° senza). A pari angolo il
   costo è +144 % (metà ciclo) o +288 % (limite conservativo).
3. **Natura fisica del latch.** Il latch non è un cavo: è un vincolo
   reversibile. Braccio = linker di lunghezza 2R; contrappeso = bead di
   volume V/2; latch = spot ottico secondario, linker DNA, bead magnetica o
   foto-switch. La scena Blender dichiara che i cavi sono un diagramma del
   collegamento topologico.

## Correzioni alla V7.1

- La modulazione di ω è **74.2 %** (`1 - β_rot/β_tot`), non 84.2 %.
  L'84.2 % è `1 - I_rot/I_tot`: è la frazione di inerzia, conta nel
  transitorio (`τ_spin ≈ 7 ns`), non nella velocità di regime.
- Il subarmonico a `f_rot/2` è nello **stato A/B** del latch, non
  nell'orientazione. L'orientazione è 360°-periodica; solo un osservabile
  dispari sotto 360° ha periodo 720°.
- I picchi di ΔT vs R sono zeri di `<1-cos δ>` (forma di ritardo), non
  risonanze Mie del modello. Con dispersione di taglia del 15 %, il picco a
  2.34 µm scende da 27.2 a 12.5 K; con 20 % a 7.3 K. Il buco termico è un
  artefatto di monodispersità.
- La soglia di 20 K è 57 °C, coagulazione. La soglia prudente è 5 K
  (42 °C). Il punto sicuro più veloce è R ≈ 0.5–0.6 µm,
  T_720 ≈ 7–12 s, E_720 ≈ 1.4–1.8 mJ.

## File principali

| file | cosa fa |
|---|---|
| `codice/v72_contrappeso.py` | griglia, costo energetico, rumore, viscosità, latch |
| `codice/v72_figure.py` | figure `figure/fig_V72_*.png` |
| `codice/v72_rumore.py` | Langevin rotazionale e figura del rumore |
| `codice/v72_tabelle.py` | tabelle LaTeX da `dati/V72_sintesi.json` |
| `codice/v72_blender.py` | scena e render `blender/v72_*.png` |
| `documento/Preprint_V72_appendice.tex` | appendice compilata in PDF |
| `dati/V72_griglia.json` | griglia completa (1920 celle) |
| `dati/V72_sintesi.json` | sintesi, costo, finestre, rumore, latch |
| `dati/V72_rumore.json` | simulazione di Langevin |

## Come si esegue

```powershell
cd codice
python v72_contrappeso.py
python v72_figure.py
python v72_rumore.py
python v72_tabelle.py

# Blender (Windows)
& 'C:\Program Files\Blender Foundation\Blender 5.2\blender.exe' -b -P v72_blender.py

# Preprint
cd ..\documento
pdflatex Preprint_V72_appendice.tex
pdflatex Preprint_V72_appendice.tex
```

## Cosa non dimostra

La V7.2 non dimostra che il latch esista. Calcola quanto costerebbe se
esistesse e come si falsificherebbe. I limiti principali sono:

- sangue non newtoniano: 4 mPa s è il valore ad alto shear; a basso shear
  può essere 10 mPa s, e T_720/E_720 raddoppiano;
- R = 0.3 µm è più piccolo di un globulo rosso: il continuo di Stokes è
  un'approssimazione;
- la termica è una sorgente lineare in mezzo infinito; la perfusione
  (Pennes) e la parete del capillare non sono modellate;
- il canale a ritardo usa `δ = 2π Δn (2R)/λ`, non un calcolo Mie completo;
- il tempo di risposta del latch non è misurato.
