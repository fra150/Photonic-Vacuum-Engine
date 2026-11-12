# Protocollo di misura V7.3
Sezione 13 · 12 novembre 2026 · una pagina
Ogni numero viene dai JSON in `dati/`: nessun valore è scritto a mano.
Riferimento: R = 0,5 µm, λ = 800 nm, z = 2,0 µm, plasma 1,3 mPa·s, n = 1,33, T = 310 K.
Obiettivo: misurare latch, termica e cancello; non dimostrare il motore.

---

## Setup 1 — Latch: rapporto e rumore appaiato
**Misura.** 3 mezzi-giri veloci contro 3 mezzi-giri lenti sullo stesso bead (r = 0,3969 µm, linker 1,00 µm), alla stessa potenza. Un mezzo-giro è metà del ciclo a 720°, cioè un giro di 360°. Il mezzo-giro veloce (contrappeso disinnestato) vale 0,5996 s, il lento (contrappeso innestato) 2,327 s; la somma è 2,926 s, il rapporto è 3,881, con la correzione idrodinamica non validata fino a 4,492.
**Atteso.** Rapporto 3,881 per la geometria ideale; 4,256 con il termine traslazionale della bead, 4,335 con il termine completo, 4,492 come limite superiore.
**Errore.** CV = 17,84%; σR/R = CV·√(2/3) = 14,57%; ε per giro = 1,07×10⁻⁴ (soglia geometrica) e 1,01×10⁻⁴ (soglia ottimizzata); ε sulle somme di 3 = 1,95×10⁻¹¹ (geometrica).
**Falsificazione.** Rapporto fuori da 3,881–4,492 di 2σR (sotto 2: nessuna asimmetria CW/No); un CV misurato sopra 17,76% rende irraggiungibile il budget 1×10⁻⁴ per giro, e sopra 30,77% lo rende irraggiungibile anche sulle somme.
**Nota.** Il criterio per giro non migliora con N; le somme migliorano con la radice di N. A CV = 30,2% il per giro vale 2,21×10⁻² (fuori budget) e le somme 7,38×10⁻⁵ (dentro). Il test sulla durata totale darebbe 1,14×10⁻² e non si usa.

## Setup 2 — Termica: forcella pozzo-specchio
**Misura.** ΔT allo stato stazionario in sangue, R = 0,5 µm, z = 2,0 µm, a 800 nm e 1064 nm; potenza incidente 0,6628 mW e 0,8507 mW.
**Atteso.** A 800 nm: pozzo 0,3471 K, nessuna correzione 0,5086 K, specchio 0,6701 K. A 1064 nm: 0,2827 K, 0,3906 K, 0,4986 K.
**Errore.** Forcella di modello ±0,161 K a 800 nm e ±0,108 K a 1064 nm; risoluzione dichiarata 0,05 K; il criterio è 5,0 K.
**Falsificazione.** ΔT ≥ 5,0 K a una delle due lunghezze d'onda, oppure ΔT sopra lo specchio di 2σ. Il criterio 5,0 K è scelto, non una soglia biologica.
**Margine.** Lo specchio resta sotto 5,0 K di un fattore 7,46 a 800 nm e 10,03 a 1064 nm.

## Setup 3 — Cancello: soglia di shear e fuga
**Misura.** Rampa di shear in plasma (μ = 1,3 mPa·s–1,5 mPa·s) e in sangue (μ = 8,0 mPa·s, σy = 10,0 mPa); fuga a Q = 0,1.
**Atteso.** In plasma γ⋆ = 2,700 s⁻¹, con 2,340 s⁻¹ a μ = 1,5 mPa·s; in sangue 0,2777 s⁻¹ con M_y = 4πR³σy, e banda 0,1037 s⁻¹–0,4094 s⁻¹. A Q = 0,1: U₀ = 8,59 kT, τesc = 200,7 s.
**Errore.** In plasma ±0,180 s⁻¹ (6,67%), interamente dalla banda di μ; in sangue la banda 0,1037 s⁻¹–0,4094 s⁻¹ copre μ = 7,0 mPa·s–9,0 mPa·s e σy = 5,0 mPa–20,0 mPa. La convenzione di fuga è τD = ξ/krad = 0,0374 s, non w₀²/(4D) = 0,321 s.
**Falsificazione.** Rotazione persistente a 267,0 s⁻¹ (venula, 98,88× la soglia) o γ⋆ fuori banda di 2σ; nessuna fuga in 2007 s a Q = 0,1 falsifica la barriera.
**Nota.** Il prefattore Bingham è 4πR³σy, conservativo; la soluzione esatta della sfera darebbe un fattore 1,5 più basso. Il sangue intero resta una reologia efficace, non una curva misurata.

## Criteri comuni
- **Soglia prima dei dati.** La soglia geometrica si fissa prima di guardare i tempi; l'ottimizzata si riporta solo come controllo.
- **Appaiamento.** Ogni ciclo veloce è appaiato al suo lento; il test sulla durata totale non si usa.
- **Budget.** Il budget è 1×10⁻⁴ per confronto; il CV massimo è 17,76% per giro e 30,77% sulle somme.
- **Referto.** Si riportano insieme valore, errore e falsificazione; un risultato fuori banda si dichiara, non si scarta.

## Che cosa questo protocollo non prova
- **Non prova la sicurezza biologica.** 5,0 K è un criterio termico scelto, non una soglia di danno.
- **Non prova il latch.** Misura rapporto e rumore; non misura l'attuazione, l'energia o la fatica.
- **Non prova il funzionamento in flusso.** Misura la soglia in stasi; il sangue intero resta una reologia efficace, non una curva misurata.
