# Motore Fotonico del Vuoto — V7.3

Data: 10 novembre 2026

La V7.3 non aggiunge un modello: chiude i conti della V7.2, dichiara che cosa è
calcolato e a quali ipotesi, e indica quale misura potrebbe falsificarlo. Non
introduce nuove lunghezze d'onda e non introduce varianti del modello.

La conclusione della revisione è una sola, e non è quella che la V7.2 lasciava
intendere:

> **Non esiste un punto di funzionamento contemporaneamente dentro il dominio
> del continuo e lontano dallo snervamento.** Il punto submicrometrico è fragile,
> il latch non è realizzato e il flusso fisiologico è chiuso dal cancello.

Questa cartella contiene tutto ciò che serve per rifare i numeri, le tabelle, le
figure e i render. Contiene anche il preprint, in italiano e in inglese:
`paper/it/Preprint_V73_it.tex` è la versione di riferimento,
`paper/en/Preprint_V73_en.tex` è il mirror. Il testo italiano è quello da
citare; l'inglese non è un paper diverso.

---

## 1. Che cosa c'è qui

| Cartella | Contenuto |
|----------|-----------|
| `codice/` | Gli script di calcolo, le tabelle, le figure, i render e la suite di test |
| `dati/` | I JSON prodotti dagli script; nessun numero è scritto a mano |
| `figure/` | Anteprime PNG delle figure, in italiano e in inglese |
| `blender/` | Il file `.blend` e i cinque render della V7.3, in italiano e in inglese |



I numeri citati nel seguito vengono tutti dai JSON in `dati/`. Se cambia il
modello, cambia il JSON e cambiano tabelle, macro e figure.

---

## 2. Riproduzione

### 2.1 Requisiti

- Python 3.13 (verificato con 3.13.13)
- NumPy 2.2.6
- Matplotlib 3.10.1
- Blender 5.2 LTS, solo per i render
- pdfTeX / TeX Live 2026, solo per compilare le tabelle in un documento

### 2.2 Pipeline completa

L'ordine conta: il cancello legge i risultati del motore, gli ambienti leggono il
cancello, le tabelle leggono tutti i JSON.

```powershell
cd "versione 7.3"

python codice\v73_contrappeso.py    # motore: griglia, ciclo 720, costo, viscosità
python codice\v73_taglio.py         # cancello sul flusso: quattro uscite
python codice\v73_ambienti.py       # tabella ambienti A–D
python codice\v73_termica.py        # termica a mezzo uniforme, sink/mirror
python codice\v73_rumore.py         # rumore esatto (inverse Gaussian), latch
python codice\v73_tabelle.py        # tabelle LaTeX e macro
python codice\v73_figure.py         # figure italiane
$env:V72_FIG_LANG="en"; python codice\v73_figure.py; Remove-Item Env:V72_FIG_LANG
python codice\v73_test.py           # 222 controlli di regressione
```

I render Blender si rifanno a parte. Il default è l'italiano; la versione
inglese si genera con `V73_LINGUA=en` e scrive file con suffisso `_en`:

```powershell
& "C:\Program Files\Blender Foundation\Blender 5.2\blender.exe" `
    -b -P codice\v73_blender.py
$env:V73_LINGUA="en"
& "C:\Program Files\Blender Foundation\Blender 5.2\blender.exe" `
    -b -P codice\v73_blender.py
Remove-Item Env:V73_LINGUA
```

`v73_termica.py` e `v73_rumore.py` scrivono da soli le proprie figure in
`paper/it/figures` e `paper/en/figures`; `v73_figure.py` scrive le nove figure
del modello.

### 2.3 Che cosa produce ogni script

| Script | Output |
|--------|--------|
| `v73_contrappeso.py` | `V73_griglia.json` (540 celle), `V73_sintesi.json` |
| `v73_taglio.py` | `V73_taglio.json` (quattro uscite, due rami, correzione idrodinamica) |
| `v73_ambienti.py` | `V73_ambienti.json`, `tab_v73_ambienti.tex` IT/EN |
| `v73_termica.py` | `V73_termica.json`, tabella e figura termica IT/EN |
| `v73_rumore.py` | `V73_rumore.json`, tabella, figura e specifica del latch |
| `v73_tabelle.py` | 14 tabelle IT/EN e `paper/numeri_v73.tex` |
| `v73_figure.py` | 9 figure PDF/PNG IT/EN |
| `v73_blender.py` | 5 PNG IT, 5 PNG EN con `V73_LINGUA=en`, e `motore_fotonico_v73.blend` |
| `v73_test.py` | Nessun file: verifica che tutto quanto sopra sia coerente |

`v7_mappa.py` e `v7_blender.py` sono helper importati dagli script V7.3; non
vanno eseguiti da soli.

---

## 3. Il cancello sul flusso

Il cancello non dipende dalle costanti ottiche ancora da calibrare
(`lambda`, `Delta_n`, `Q`, `eta_int`): dipende da `tau_ott`, `beta_cw` e `mu`.
Per questo la sua conclusione è robusta.

Quattro uscite, tutte in `s^-1`, ma solo due sono soglie: la coppia e la forza.
Le altre due sono una verifica di coerenza (Jeffery) e una barriera in kT (fuga).

| Uscita | Formula | Valore |
|--------|---------|--------|
| 1. Coppia, corpo vincolato | `(tau_ott - M_y) / beta_cw` | **2.700** |
| 1b. Coppia, sfera libera | `2 tau_ott / beta_no` | 20.96 (solo convenzione) |
| 2. Forza di trappola | `F_trap / (6 pi mu R_eff z_c)`, `Q = 0.1–0.3` | 5.35 / 10.69 / 16.04 |
| 3. Jeffery | `gamma` a `T720` | 5.20 (solo coerenza) |
| 4. Fuga | `U0 = F w0 / 4`, `tau_esc = tau_D exp(U0/kT)` | 8.6 / 17.2 / 25.8 kT |

La soglia è il minimo delle due uscite che sono soglie, coppia e forza:

```
gamma_crit = min(2.700, 5.35, 10.69, 16.04) = 2.700 s^-1
```

Il minimo è la scelta conservativa. L'intervallo 2.7–5.4 `s^-1` è una
**coincidenza** (fattore 1.98), non una finestra: Jeffery a 5.20 cade dentro
quell'intervallo, ma non definisce la soglia. La V7.3 ritratta la "finestra
2.5–8.7 `s^-1`" della V7.2.

### 3.1 Sangue e Bingham

Nel sangue lo snervamento entra come **scostamento additivo** della coppia:

```
tau_app = tau_ott + M_y,     M_y = 4 pi R^3 sigma_y
```

La convenzione `4 pi R^3` è conservativa: la soluzione esatta della sfera di
Bingham dà `(8/3) pi R^3 sigma_y`, un fattore 1.5 più bassa. La soglia nel
sangue scende a:

| `mu` (mPa s) | `sigma_y` (mPa) | `gamma_crit` (s^-1) |
|--------------|-----------------|---------------------|
| 7 | 5 | 0.409 |
| 8 | 10 | 0.278 |
| 9 | 20 | 0.104 |

Il valore ufficiale per l'ambiente B è `mu_B = 8.0 mPa s`, con banda 7–9. Il
valore 4.0 mPa s non è l'ambiente B: è il riferimento nominale del continuo usato
solo nella tabella a due zone del CFL. A `Pe_r` fissato la coppia è
`tau_ott = Pe_r kT`, indipendente da `mu`: la viscosità cambia durata ed energia
ottica, non la coppia istantanea.

Il rapporto 3.881 tra i due stati (vincolato e libero) resta invariato con lo
scostamento costante.

### 3.2 Perché il flusso è chiuso

Lo shear fisiologico è 99–741 volte sopra la soglia:

| Distretto | `gamma` (s^-1) | Rapporto |
|-----------|----------------|----------|
| Venula 30 µm | 267 | ×98.9 |
| Capillare 6 µm | 667 | ×247 |
| Arteriola | 2000 | ×741 |

Il verdetto è nel JSON: *il flusso fisiologico è 99–741 volte sopra la soglia; il
dispositivo vive in stasi ex vivo o in stasi in vivo.*

### 3.3 Correzione idrodinamica

Il rapporto 3.881 è una geometria ideale. La correzione a due sfere vale
`(r_cw/d)^2 = 0.1575` e porta il rapporto a:

| Correzione | Rapporto | Variazione | Limite CV appaiato |
|------------|----------|------------|--------------------|
| Nessuna (conservativa) | 3.881 | — | 17.83% |
| Termine traslazionale della bead | 4.256 | +9.7% | 19.10% |
| Termine completo della bead | 4.335 | +11.7% | 19.36% |
| Limite superiore (tutto il drag cw) | 4.492 | +15.7% | 19.86% |

La correzione non è validata; nel calcolo del rumore si usa il valore
conservativo 3.881.

---

## 4. Gli ambienti

Il CFL è un **limite del dominio, non un ambiente**: esiste solo con flusso, e il
flusso è chiuso dal cancello. "Plasma + stasi" si realizza solo in vitro.

| ID | Ambiente | `mu` (mPa s) | Stato |
|----|----------|--------------|-------|
| A | Plasma/tampone in vitro, stasi | 1.3–1.5 | Ammesso |
| B | Sangue intero in vivo, stasi | 8.0 (7–9) | Ammesso con riserva |
| C | Sangue intero in vivo, flusso | 4–10 | **Escluso** |
| D | Microfluidica ex vivo, shear controllato | 1.3–10 | Ammesso solo sotto soglia |

Il ramo A ha `T720 = 2.93 s` ed `E720 = 1.94 mJ`. Il ramo B, con lo snervamento,
sale a `T720 = 15.8–20.3 s` ed `E720 = 12.4–23.3 mJ`.

---

## 5. La termica

Modello a mezzo uniforme, senza Pennes. La giustificazione è dichiarata: il tempo
di perfusione è circa 60 s contro `T720 = 2.93 s` (peso 4.8%), e la lunghezza di
diffusione termica a 3 s è 640 µm: il calore esce dal vaso.

Quattro contributi:

- `dT_part`: assorbimento della sfera (Mie, con tetto geometrico)
- `dT_pre`: assorbimento del fascio dall'interfaccia alla sfera
- `dT_post`: assorbimento oltre il fuoco, fino a `z_R = pi w0^2 / lambda`,
  con prefattore `2 pi` da sorgente lineare
- `dT_diff`: luce diffusa nei pressi della sfera,
  `mu_a P_sca / (4 pi kappa)`, con `P_sca` limitata da `eta_int P_inc`

La parete ha un segno esplicito:

| Parete | Significato | Formula |
|--------|-------------|---------|
| `sink` | Isoterma (vetro, endotelio) | `dT_0 - dT_pre - img` |
| `0` | Nessuna correzione | `dT_0` |
| `mirror` | Adiabatica | `dT_0 + dT_pre + img` |

con `img = dT_part R / (2 z)`. Il `sink` è il limite inferiore, il `mirror` il
limite superiore: la forcella è dichiarata, non scelta a posteriori.

**Il criterio di 5 K è un criterio termico scelto, non una prova di sicurezza
biologica.** Le tabelle dicono "sotto 5 K", non "sicuro".

Caso sangue, `z = 2 µm`:

| `lambda` (nm) | `R` (µm) | `sink` (K) | nessuna correzione (K) | `mirror` (K) |
|---------------|----------|------------|------------------------|--------------|
| 532 | 0.5 | 5.72 | 10.87 | 16.02 |
| 532 | 0.3 | 22.25 | 45.06 | 67.88 |
| 800 | 0.5 | 0.35 | 0.51 | 0.67 |
| 800 | 0.3 | 1.13 | 1.78 | 2.42 |
| 1064 | 0.5 | 0.28 | 0.39 | 0.50 |
| 1064 | 0.3 | 0.89 | 1.29 | 1.70 |

A 532 nm il sangue fallisce **anche nel pozzo al valore centrale** `Delta n =
0.08`: non c'è una scelta di parete che lo salvi. Per `Delta n = 0.09–0.10` il
pozzo passa (4.61–3.82 K) ma lo specchio resta sopra 10 K: la conclusione cambia
con la parete, e un caso che cambia con la parete non è robusto. A 800 e 1064 nm
passa **anche nello specchio**: la conclusione non dipende dal segno della parete.

Il plasma non usa più l'ottica del sangue: `mu_a = 0.044 / 2.3 / 15 m^-1` a
532/800/1064 nm. Nel sangue `mu_a(532) = 2.35e4 m^-1`, non 7500.

### 5.1 La calibrazione di `Delta n_eff`

`Delta n = 0.08` non è una scelta ereditaria: è una grandezza misurabile. Parkin
e collaboratori danno `Delta n = 0.06`; Modin, Ben Zion e Chaikin danno
`0.075 ± 0.015`, con letteratura 0.06–0.09. Il valore 0.10 è il monocristallo,
non la sfera efficace; 0.172 è la calcite. Il valore 0.08 è dentro `1 sigma`
della misura di Modin e il 33% sopra quella di Parkin.

L'effetto sui numeri (a `R = 0.5 µm`, `z = 2 µm`):

| `Delta n` | `P_inc` 800 (mW) | `Delta T_0` 800 (K) | `P_inc` 532 (mW) | `Delta T_0` 532 (K) | sink 532 (K) | mirror 532 (K) |
|-----------|------------------|---------------------|------------------|---------------------|--------------|----------------|
| 0.06 | 1.161 | 0.89 | 0.888 | 18.66 | 9.82 | 27.51 |
| 0.075 | 0.751 | 0.58 | 0.583 | 12.25 | 6.44 | 18.05 |
| 0.08 | 0.663 | 0.51 | 0.517 | 10.87 | 5.72 | 16.02 |
| 0.09 | 0.529 | 0.41 | 0.417 | 8.77 | 4.61 | 12.93 |
| 0.10 | 0.433 | 0.33 | 0.346 | 7.27 | 3.82 | 10.72 |

La conclusione non cambia: 532 nm non è utilizzabile, 800 e 1064 nm sì. Cambia la
lettura del margine: il bordo conservativo non è 0.10, è 0.06; e a 532 nm il
fallimento nel pozzo vale per `Delta n <= 0.08`, non per tutta la banda.

---

## 6. Rumore e latch

Il rumore non è gaussiano. La durata di un giro è la **prima passa esatta** di un
processo diffusivo, cioè una inverse Gaussian:

```
T ~ IG(m, lambda),    CV^2 = 1/(pi Pe_r) per 360 gradi
```

La somma di `N` giri resta inverse Gaussian, con shape `N lambda`. Da qui una
conseguenza che la V7.2 sbagliava: **l'errore non migliora con `N`**, perché la
shape e la soglia scalano insieme.

Il test obbligatorio è **appaiato**: la somma dei 3 mezzi-giri veloci contro la
somma dei 3 mezzi-giri lenti. Il test sulla durata totale del ciclo non è
utilizzabile.

| Grandezza | Valore |
|-----------|--------|
| Errore appaiato a CV ideale 17.84% | `1.01e-4` |
| Errore appaiato geometrico | `1.07e-4` |
| Errore a durata totale | `1.14e-2` |
| Limite CV appaiato per `10^-4` | **17.83%** |
| Limite CV a durata totale | 11.58% |
| Errore appaiato a CV 30.2% | `2.07e-2` |

Il punto ideale è **già fuori budget**: CV 17.84% contro un limite di 17.83%,
cioè `1.01e-4` contro `1e-4`. Non è "senza margine": è dalla parte sbagliata del
limite già nello scenario ideale. Il target del 30% non è supportato da nessun
calcolo esatto.

La specifica del latch:

| Voce | Valore | Stato |
|------|--------|-------|
| Rapporto | 3.881 (fino a 4.492 corretto) | Geometria ideale verificata; correzione non validata |
| CV | 17.8% ideale | Saturato; ogni rumore aggiuntivo lo porta sopra |
| Tempo di commutazione | `< 0.23 s` richiesto, `5.4e-8 s` il rilassamento viscoso | Attuazione non modellata |
| Energia | 10–60 kT dentro 62.8 kT per giro | Non dimostrata |
| Fatica | `10^6` cicli (~34 giorni a `T720`) | Non misurata |

**Il latch non è dimostrato e non ha margine.** Il JSON lo dice esplicitamente.

---

## 7. I sette tagli della V7.3

1. **Soglia, non finestra.** `gamma_crit = 2.700 s^-1`; l'intervallo 2.7–5.4 è
   una coincidenza, non una finestra.
2. **Jeffery retrocesso.** 5.20 `s^-1` è una verifica di coerenza del tumbling,
   non la soglia; sta sopra il minimo e non lo definisce.
3. **`F_trap` nella banda `Q` dichiarata.** 0.1–0.3, non un parametro scelto per
   far passare il caso.
4. **CFL e stasi mutuamente esclusivi.** Il CFL esiste solo con flusso; "plasma +
   stasi" solo in vitro.
5. **Etichette termiche corrette.** 5 K è un criterio scelto, non una prova di
   sicurezza; le tabelle dicono "sotto 5 K".
6. **`mu_a` del plasma corretto.** 0.044/2.3/15 `m^-1`, separato dall'ottica del
   sangue.
7. **Nessun punto fisso per `mu(gamma_dot)` nel plasma.** Il plasma è newtoniano:
   non c'è un punto fisso da risolvere, c'è una banda dichiarata (1.3–1.5 mPa s)
   e la soglia scala linearmente con `mu`.

---

## 8. Limiti dichiarati

- **Nessun punto di funzionamento validato.** La V7.3 non dice "funziona"; dice
  che cosa sarebbe necessario per funzionare e perché non è stato dimostrato.
- **Il cancello è conservativo.** Usa `4 pi R^3 sigma_y`; la soluzione esatta
  della sfera di Bingham darebbe un fattore 1.5 in meno.
- **La correzione idrodinamica non è validata.** Il rapporto sta tra 3.881 e
  4.492; il calcolo del rumore usa il valore conservativo.
- **La termica è a mezzo uniforme.** Niente Pennes, niente convezione, niente
  fotoforesi. La parete è un caso dichiarato, non una misura.
- **Il rumore è la prima passa esatta di un giro.** Non include deriva del
  laser, attrito statico, né urti con cellule.
- **Il latch non è realizzato.** L'attuazione non è modellata; la fatica non è
  misurata.
- **L'aggregazione nel sangue intero non è modellata.** Rouleaux ed ematocrito
  pieno entrano solo come `mu` e `sigma_y` efficaci.
- **Il piano dell'orbita della bead non è misurato.** La correzione a due sfere
  assume la geometria dichiarata.
- **`Q` non è misurato.** Il cancello non ne dipende; l'uscita 2 sì.
- **Il 5 K non è una soglia biologica.** È un criterio termico scelto.

---

## 9. Ritrattazioni

Le ritrattazioni sono conservate nei JSON e riprodotte dagli script, non
cancellate:

- La "finestra 2.5–8.7 `s^-1`" non esiste: è una soglia a 2.700 `s^-1`.
- I valori 1.74 / 0.64 / 0.90–0.57 pN e i fattori ×19 / ×194 sono ritrattati.
- La "finestra CFL verde, non IR" è ritrattata: il guadagno è 5.3e5 / 190 / 13,
  ma il CFL esiste solo con flusso e non entra negli ambienti A, B, D.
- L'asserzione che la V7.2 avesse un errore in `mu(gamma_dot)` nel plasma è
  ritrattata: il plasma è newtoniano.
- Il crossover a `R = 0.422 µm` è un artefatto di Rayleigh, ritrattato.
- Il limite CV del 30.2% e l'errore `8.8e-7` su tre cicli sono artefatti
  gaussiani, ritrattati: il limite esatto è 17.83% e l'errore non dipende da `N`.
- L'eq. (7) della V7.2 dà `lambda`, non `lambda^2`.

---

## 10. Domande aperte

1. **Misura di `Q`.** È l'unico input del cancello che non è stato misurato, e
   governa l'uscita 2.
2. **Piano dell'orbita della bead.** La correzione idrodinamica dipende dalla
   geometria; serve una misura, non un'assunzione.
3. **Cinetica di aggregazione nel sangue intero.** `mu` e `sigma_y` efficaci
   sono una parametrizzazione, non una reologia misurata.
4. **Attuazione del latch.** Il requisito è `< 0.23 s`; il rilassamento viscoso è
   `5.4e-8 s`. Manca il meccanismo, non il tempo.
5. **Fatica.** `10^6` cicli non sono stati misurati.
6. **Validazione termica.** La forcella sink/mirror è un calcolo, non un
   termometro.

---

## 11. Test

```powershell
python codice\v73_test.py
```

La suite esegue **222 controlli**: identità del cancello, ambienti, forcella
termica (incluso il termine immagine), limiti esatti della inverse Gaussian,
scaling del modello, conteggi dei JSON, presenza e contenuto di tabelle, macro,
figure e render. Non è una validazione fisica: è una regressione. Se un numero
cambia senza che cambi il modello, il test lo segnala.

---

## 12. Il preprint

Il preprint è incluso, in due lingue: https://doi.org/10.5281/zenodo.24137503 


