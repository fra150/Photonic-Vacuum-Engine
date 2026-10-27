# Nota tecnica — Formalizzazione matematica dei rapporti tra Spin e Contro-spin

**Motore Fotonico v4** — sviluppo futuro n. 4 del paper
*"Modellazione geometrica 3D del sistema Spin / Contro-spin per il Motore Fotonico v4"*
Dr. Francesco Bulla — Catania, 25 ottobre 2026

> **Perimetro.** Questa nota formalizza le **relazioni geometriche e cinematiche** già
> contenute nel modello 3D (Figure 1–8) e ne verifica la coerenza dimensionale. Non è una
> validazione fisica del dispositivo: i valori numerici derivano da un modello a corpi
> rigidi con densità uniforme e non includono dissipazione, attrito viscoso, effetti
> plasmonici o accoppiamento con il vuoto. Ogni numero va letto come **stima di ordine di
> grandezza coerente con la geometria**, non come previsione.
>
> **File collegati** (stessa cartella)
> - `Figura_11_formalizzazione.png` — i quattro pannelli (a)–(d) discussi qui sotto
> - `parametri_formalizzazione.json` — tutti i valori usati e ottenuti
> - `formalizzazione_spin_controspin.py` — codice che genera figura e parametri
> - `Figura_9_interazione_flussi.png`, `Figura_10_scale_dimensionali.png` — sviluppi 2 e 3

---

## 1. Notazione

| Simbolo | Significato | Unità |
|---|---|---|
| $1\ \mathrm{u.a.}$ | unità di scena Blender = 100 nm | $10^{-7}$ m |
| $I_1$ | momento d'inerzia del sistema **Spin** (3 orbite ellittiche) | kg·m² |
| $I_2$ | momento d'inerzia del sistema **Contro-spin** (disco + 2 lenti) | kg·m² |
| $\omega_1,\ \omega_2$ | velocità angolari dei due sistemi | rad/s |
| $L_1,\ L_2,\ L_{tot}$ | momenti angolari e loro somma | kg·m²/s |
| $\tau_{net}$ | coppia netta sul sistema accoppiato | N·m |
| $P$ | potenza ottica in ingresso | W |
| $\eta_s,\ \eta_c$ | efficienze di accoppiamento Spin e Contro-spin | – |
| $\alpha$ | angolo di disallineamento tra i due assi di spinta | rad |
| $F_z$ | spinta assiale risultante | N |

**Convenzione di segno.** $\omega_1 > 0$ (Spin antiorario visto da $+z$), $\omega_2 < 0$
(Contro-spin orario). Tutti i rapporti che seguono sono **adimensionali** e quindi
invarianti rispetto al fattore di scala $1\ \mathrm{u.a.} = 100$ nm: cambiare la scala
non cambia la fisica del modello, cambia solo i valori assoluti di $I$, $L$, $E$.

---

## 2. Lo stato di spin è uno spinore: copertura doppia di SU(2)

Il sistema Spin compie una rotazione di $720°$ per ciclo (due giri completi), non $360°$.
La ragione formale è che lo stato non è un vettore ma uno **spinore**, cioè un elemento
della rappresentazione fondamentale di SU(2), che è una copertura doppia di SO(3):

$$
U(\theta) = \exp\!\left(-i\,\frac{\theta}{2}\,\sigma_z\right)
= \cos\frac{\theta}{2}\,I - i\,\sin\frac{\theta}{2}\,\sigma_z
$$

da cui, per $\sigma_z$ arbitrario,

$$
U(2\pi) = -I \qquad\qquad U(4\pi) = +I .
$$

**Pannello (a) di Figura 11.** Re $\psi = \cos(\theta/2)$ e Im $\psi = \sin(\theta/2)$:
a $\theta = 360°$ lo spinore cambia segno ($\psi \to -\psi$), a $\theta = 720°$ torna
esattamente allo stato iniziale. Questo è il motivo per cui la timeline è costruita su
**240 frame = 2 giri** e non su 120: un solo giro lascerebbe il sistema in uno stato
fisicamente distinto (segno opposto), e il ciclo non sarebbe chiuso.

Conseguenza operativa: la periodicità del sistema combinato è $\theta_1 \in [0, 4\pi]$,
cioè il **doppio** della periodicità geometrica visibile.

---

## 3. Momenti di inerzia ricavati dalla geometria

I momenti di inerzia non sono parametri liberi: si calcolano dai fili della scena
trattati come tubi di rame ($\rho_{Cu} = 8960$ kg/m³) con le sezioni dichiarate negli
script di modellazione.

$$
V = L_{filo}\,\pi\,r_{sezione}^2 \qquad
m = \rho_{Cu}\,V \qquad
I = m\,\langle x^2 + y^2\rangle
$$

dove $\langle x^2+y^2\rangle$ è la media quadratica delle distanze dall'asse, misurata
**direttamente sui vertici delle curve** della scena Blender (non stimata a occhio).

| Grandezza | Spin | Contro-spin |
|---|---|---|
| Lunghezza totale dei fili | $43.240$ u.a. | $42.986$ u.a. |
| Volume di rame | $0.1391$ u.a.³ | $0.5009$ u.a.³ |
| Massa | $1.246\times10^{-18}$ kg | $4.488\times10^{-18}$ kg |
| $\langle x^2+y^2\rangle$ | $2.4807$ u.a.² | $5.7567$ u.a.² |
| **Momento d'inerzia $I$** | $\mathbf{3.092\times10^{-32}}$ kg·m² | $\mathbf{2.583\times10^{-31}}$ kg·m² |

Il Contro-spin ha momento d'inerzia **8.35 volte** maggiore: massa 3.6× e raggio
quadratico medio 2.3×. È questa asimmetria, e non una scelta arbitraria, a fissare il
rapporto di rotazione.

---

## 4. Bilancio del momento angolare → legge di contro-rotazione

Se il sistema combinato è **isolato** rispetto alle coppie esterne attorno all'asse $z$,
il momento angolare totale si conserva:

$$
L_{tot} = I_1\omega_1 + I_2\omega_2 = 0
\quad\Longrightarrow\quad
\boxed{\ \frac{\omega_2}{\omega_1} = -\frac{I_1}{I_2}\ }
$$

Il segno meno è la **legge di contro-rotazione**: i due sistemi non possono girare nello
stesso verso se il momento angolare totale deve restare nullo.

Valore del modello:

$$
\frac{I_1}{I_2} = \frac{3.092\times10^{-32}}{2.583\times10^{-31}} = 0.11968
\quad\Longrightarrow\quad
\omega_2 = -\,0.11968\;\omega_1
$$

ovvero: **il Contro-spin compie 1 giro ogni 8.36 giri dello Spin**.

**Pannello (b) di Figura 11.** Con $\omega_1 = 2\pi$ rad/s (1 Hz):

| Termine | Valore |
|---|---|
| $L_1 = I_1\omega_1$ | $+1.943\times10^{-31}$ kg·m²/s |
| $L_2 = -I_2\omega_2$ con $\omega_2/\omega_1 = -I_1/I_2$ | $-1.943\times10^{-31}$ kg·m²/s |
| $L_{tot}$ | $0$ (esatto, per costruzione) |
| $L_{tot}$ residuo se si imponesse $\omega_2/\omega_1 = -0.85$ (valore usato nell'animazione) | $2.91\times10^{-32}$ kg·m²/s |

Il residuo della curva grigia vale il **15% di $L_1$**: è la misura di quanto
l'animazione (dove il rapporto è fissato a $0.85$ per ragioni di leggibilità visiva) si
discosti dal vincolo dinamico. Il rapporto $0.11968$ è quello fisicamente coerente; il
rapporto $0.85$ serve solo a rendere percepibile il moto relativo a occhio.

---

## 5. Condizione di coppia nulla

Derivando il bilancio rispetto al tempo si ottiene la condizione di **coppia netta
nulla**, cioè di funzionamento senza accumulo di momento angolare reattivo:

$$
\tau_{net} = I_1\alpha_1 + I_2\alpha_2 = 0
\qquad\text{con}\qquad
\alpha_i = \dot\omega_i .
$$

Poiché $\omega_2/\omega_1$ è costante, $\alpha_2/\alpha_1 = -I_1/I_2$ e la relazione è
soddisfatta **a ogni istante**, anche in transitorio (accelerazione e decelerazione
simmetriche). Ne segue la cancellazione della componente trasversale della spinta: i due
getti (Figura 9) hanno quantità di moto trasversale uguale e opposta e si annullano nel
punto di stagnazione, e resta solo la componente assiale.

**Limite dichiarato.** La coppia del campo ottico non è nulla:
$\tau_{ph} = P/\omega$, che per $P = 1$ µW a 1 Hz vale $1.59\times10^{-7}$ N·m — circa
24 ordini di grandezza sopra $\tau_{net}=0$. Il vincolo di contro-rotazione va quindi
letto come **condizione di progetto (phase-locking imposto dal campo)**, non come
conseguenza spontanea della dinamica libera: un sistema lasciato libero verrebbe
accelerato dal campo, non si auto-bilancierebbe. Questa è la discontinuità principale tra
modello e dispositivo.

---

## 6. Spinta assiale e coppia fotonica

La spinta risultante è la somma vettoriale dei due contributi, proiettata sull'asse di
propulsione:

$$
\boxed{\ F_z = \frac{P}{c}\bigl(\eta_s + \eta_c\cos\alpha\bigr)\ }
$$

con $c = 299\,792\,458$ m/s, $\eta_s = 0.70$, $\eta_c = 0.85$.

| $\alpha$ | $\eta_s + \eta_c\cos\alpha$ | $F_z$ per 1 µW |
|---|---|---|
| $0°$ | $1.550$ | $5.170$ fN |
| $15°$ (nominale) | $1.521$ | $\mathbf{5.074}$ fN |
| $30°$ | $1.436$ | $4.790$ fN |

**Pannello (c) di Figura 11**: la dipendenza da $P$ è rigorosamente lineare, mentre la
dipendenza da $\alpha$ è un coseno: il disallineamento costa solo il $6\%$ della spinta
fino a $30°$, il che rende il progetto tollerante agli errori di allineamento.

La coppia fotonica associata è

$$
\tau_{ph} = \frac{P}{\omega}
$$

e cresce al diminuire della frequenza di rotazione: a parità di potenza conviene lavorare
a bassa frequenza per massimizzare la coppia, ma a frequenza più alta si riduce il tempo
di ciclo. Da qui la tabella dei regimi di `Figura_10_scale_dimensionali.png`.

> Confronto puramente dimensionale, **non** una rivendicazione di spinta netta nel vuoto:
> $5.07$ fN equivalgono a circa $90$ volte il peso dell'assieme a $1$ g. Il confronto è
> riportato solo per mostrare che, a queste scale, le forze radiative non sono
> trascurabili rispetto alla gravità.

---

## 7. Struttura quasi-periodica del moto combinato

Le due rotazioni sono angoli su un toro $\mathbb{T}^2 = S^1_{\theta_1}\times S^1_{\theta_2}$,
con

$$
\theta_2(t) = -\frac{I_1}{I_2}\,\theta_1(t) + \theta_2(0).
$$

Il moto si chiude (cioè torna esattamente alla configurazione iniziale) **solo se**
$I_1/I_2$ è razionale. Con $I_1/I_2 = 0.1196797$ lo sviluppo in frazione continua dà

$$
0.1196797 \approx [0;\,8,\,2,\,1,\,4,\,3] = \frac{45}{376} = 0.1196809
\qquad (\text{errore relativo} \approx 1\times10^{-5})
$$

cioè il ciclo completo richiede **376 giri dello Spin e 45 giri del Contro-spin**
($376 \times 360° = 135\,360°$ di $\theta_1$). A 1 Hz di rotazione sono $\sim 6$ min 16 s;
a $0.2$ Hz (velocità dell'animazione) $\sim 31$ min.

Poiché il valore reale non è esattamente razionale, il moto è **quasi-periodico**: le
traiettorie sul toro sono dense e il sistema non ritorna mai esattamente allo stato
iniziale, ma ne passa arbitrariamente vicino. È la firma matematica della copertura
doppia SU(2) combinata con un rapporto inerziale non razionale — **pannello (d) di
Figura 11**.

Energia rotazionale del sistema combinato, in forma chiusa:

$$
E_{rot} = \tfrac12 I_1\omega_1^2 + \tfrac12 I_2\omega_2^2
= \tfrac12 I_1\omega_1^2\left(1 + \frac{I_1}{I_2}\right)
$$

A $\omega_1 = 2\pi$ rad/s: $E_{rot} = 6.83\times10^{-31}$ J, di cui il $10.7\%$ nel
Contro-spin. Il fattore $(1 + I_1/I_2)$ mostra che l'asimmetria inerziale **aggiunge**
energia al sistema: non è un vincolo gratuito.

---

## 8. Riepilogo delle relazioni formalizzate

| # | Relazione | Tipo | Verifica |
|---|---|---|---|
| 1 | $U(2\pi) = -I$, $U(4\pi) = +I$ | cinematica (SU(2)) | esatta per definizione |
| 2 | $I = m\langle x^2+y^2\rangle$, $m = \rho L \pi r^2$ | geometria | $I_1 = 3.09\times10^{-32}$, $I_2 = 2.58\times10^{-31}$ kg·m² |
| 3 | $\omega_2/\omega_1 = -I_1/I_2 = -0.1197$ | dinamica (bilancio $L$) | $L_{tot} = 0$ esatto |
| 4 | $\tau_{net} = I_1\alpha_1 + I_2\alpha_2 = 0$ | dinamica | valida anche in transitorio |
| 5 | $F_z = (P/c)(\eta_s + \eta_c\cos\alpha) = 5.074$ fN/µW | accoppiamento radiativo | lineare in $P$, $\cos$ in $\alpha$ |
| 6 | $\tau_{ph} = P/\omega$ | coppia fotonica | $1.59\times10^{-7}$ N·m a 1 µW / 1 Hz |
| 7 | $\theta_2 = -(I_1/I_2)\theta_1$ su $\mathbb{T}^2$ | topologia del moto | chiusura a $45/376$, altrimenti quasi-periodico |
| 8 | $E_{rot} = \tfrac12 I_1\omega_1^2(1 + I_1/I_2)$ | energia | $6.83\times10^{-31}$ J a 1 Hz |

---

## 9. Limiti di questa formalizzazione

1. **Corpi rigidi.** I fili sono trattati come rigidi e a densità uniforme: non si tiene
   conto della nanostruttura né della distribuzione reale della massa.
2. **Vincolo di fase imposto.** Il rapporto $\omega_2/\omega_1 = -I_1/I_2$ è una
   *condizione di progetto*, non un equilibrio dinamico spontaneo (§5).
3. **Nessuna dissipazione.** Non ci sono attrito viscoso, smorzamento strutturale o
   perdite per irraggiamento; i regimi di Reynolds di `Figura_10` mostrano però che a
   $1$ kHz–$1$ MHz il regime non è più Stokes e la dissipazione non è trascurabile.
4. **Accoppiamento col vuoto non modellato.** $F_z$ è un termine di bilancio di quantità
   di moto del campo, non una forza di propulsione nel vuoto dimostrata.
5. **Scala.** Tutto è ancorato a $1\ \mathrm{u.a.} = 100$ nm tramite il vincolo della
   nanoparticella di rame ($r = 30$ nm $\leftrightarrow 0.30$ u.a.): il modello non
   predice la scala, la assume.

## 10. Sviluppi coerenti con questa nota

- Sostituire il rapporto imposto $0.85$ con $0.1197$ in un secondo render di
  animazione, per mostrare il comportamento dinamico coerente (richiede tempi diversi per
  i due assi e quindi un nuovo rig).
- Introdurre un termine dissipativo $\tau_{diss} = -\beta\omega$ e studiare la
  condizione di regime stazionario $\tau_{ph} = \tau_{diss}$, che darebbe una frequenza
  di rotazione *prevista* invece che imposta.
- Verificare se $I_1/I_2$ può essere reso razionale per progetto (per esempio
  aggiungendo una massa di bilanciamento al Contro-spin) così da chiudere esattamente il
  ciclo a 376 giri.
