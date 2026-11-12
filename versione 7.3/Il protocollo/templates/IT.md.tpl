# Protocollo di misura ${VERSIONE}
Sezione ${SEZIONE} · ${DATA} · una pagina
Ogni numero viene dai JSON in `dati/`: nessun valore è scritto a mano.
Riferimento: R = ${R_REF}, λ = ${LAMBDA_800}, z = ${Z_REF}, plasma ${MU_PLASMA}, n = ${N_MED}, T = ${T0}.
Obiettivo: misurare latch, termica e cancello; non dimostrare il motore.

---

## Setup ${N_SETUP_1} — Latch: rapporto e rumore appaiato
**Misura.** ${N_CICLI} mezzi-giri veloci contro ${N_CICLI} mezzi-giri lenti sullo stesso bead (r = ${R_BEAD}, linker ${D_LINKER}), alla stessa potenza. Un mezzo-giro è metà del ciclo a ${GRADI_CICLO}, cioè un giro di ${GRADI_MEZZO}. Il mezzo-giro veloce (contrappeso disinnestato) vale ${T_MEZZO_VELOCE}, il lento (contrappeso innestato) ${T_MEZZO_LENTO}; la somma è ${T720_S}, il rapporto è ${RATIO_LATCH}, con la correzione idrodinamica non validata fino a ${RATIO_UP}.
**Atteso.** Rapporto ${RATIO_LATCH} per la geometria ideale; ${RATIO_BEAD} con il termine traslazionale della bead, ${RATIO_BEAD_TOT} con il termine completo, ${RATIO_UP} come limite superiore.
**Errore.** CV = ${CV_IDEALE}; σR/R = CV·√(${DUE}/${N_CICLI}) = ${SIGMA_RATIO_REL}; ε per giro = ${EPS_1_GEO} (soglia geometrica) e ${EPS_1_OPT} (soglia ottimizzata); ε sulle somme di ${N_CICLI} = ${EPS_3_GEO} (geometrica).
**Falsificazione.** Rapporto fuori da ${RATIO_LATCH}–${RATIO_UP} di ${SIGMA_2}σR (sotto ${RATIO_FALSIF_LOW}: nessuna asimmetria CW/No); un CV misurato sopra ${CV_MAX_1_GEO} rende irraggiungibile il budget ${BUDGET} per giro, e sopra ${CV_MAX_3_GEO} lo rende irraggiungibile anche sulle somme.
**Nota.** Il criterio per giro non migliora con N; le somme migliorano con la radice di N. A CV = ${CV_TARGET} il per giro vale ${EPS_AT_302_N1} (fuori budget) e le somme ${EPS_AT_302_N3} (dentro). Il test sulla durata totale darebbe ${EPS_TOTALE} e non si usa.

## Setup ${N_SETUP_2} — Termica: forcella pozzo-specchio
**Misura.** ΔT allo stato stazionario in sangue, R = ${R_REF}, z = ${Z_REF}, a ${LAMBDA_800} e ${LAMBDA_1064}; potenza incidente ${P_INC_800} e ${P_INC_1064}.
**Atteso.** A ${LAMBDA_800}: pozzo ${DT_SINK_800}, nessuna correzione ${DT_0_800}, specchio ${DT_MIRROR_800}. A ${LAMBDA_1064}: ${DT_SINK_1064}, ${DT_0_1064}, ${DT_MIRROR_1064}.
**Errore.** Forcella di modello ±${FORCELLA_800} a ${LAMBDA_800} e ±${FORCELLA_1064} a ${LAMBDA_1064}; risoluzione dichiarata ${RISOLUZIONE_TERMICA}; il criterio è ${T_CRITERIO}.
**Falsificazione.** ΔT ≥ ${T_CRITERIO} a una delle due lunghezze d'onda, oppure ΔT sopra lo specchio di ${SIGMA_2}σ. Il criterio ${T_CRITERIO} è scelto, non una soglia biologica.
**Margine.** Lo specchio resta sotto ${T_CRITERIO} di un fattore ${MARGINE_800} a ${LAMBDA_800} e ${MARGINE_1064} a ${LAMBDA_1064}.

## Setup ${N_SETUP_3} — Cancello: soglia di shear e fuga
**Misura.** Rampa di shear in plasma (μ = ${MU_PLASMA}–${MU_PLASMA_ALT}) e in sangue (μ = ${MU_BLOOD}, σy = ${SIGMA_Y_MID}); fuga a Q = ${Q_01}.
**Atteso.** In plasma γ⋆ = ${GAMMA_LIM}, con ${GAMMA_LIM_15} a μ = ${MU_PLASMA_ALT}; in sangue ${GAMMA_BLOOD_MID} con M_y = ${PREF_BINGHAM}πR³σy, e banda ${GAMMA_BLOOD_LOW}–${GAMMA_BLOOD_HIGH}. A Q = ${Q_01}: U₀ = ${U0_Q01}, τesc = ${TAU_ESC_Q01}.
**Errore.** In plasma ±${GAMMA_BAND_HALF} (${GAMMA_BAND_REL}), interamente dalla banda di μ; in sangue la banda ${GAMMA_BLOOD_LOW}–${GAMMA_BLOOD_HIGH} copre μ = ${MU_BLOOD_LOW}–${MU_BLOOD_HIGH} e σy = ${SIGMA_Y_LOW}–${SIGMA_Y_HIGH}. La convenzione di fuga è τD = ξ/krad = ${TAU_D_Q01}, non w₀²/(${QUATTRO}D) = ${W0_SQ_4D}.
**Falsificazione.** Rotazione persistente a ${GAMMA_VENULA} (venula, ${RAPPORTO_VENULA}× la soglia) o γ⋆ fuori banda di ${SIGMA_2}σ; nessuna fuga in ${FINESTRA_FUGA} a Q = ${Q_01} falsifica la barriera.
**Nota.** Il prefattore Bingham è ${PREF_BINGHAM}πR³σy, conservativo; la soluzione esatta della sfera darebbe un fattore ${BINGHAM_RATIO} più basso. Il sangue intero resta una reologia efficace, non una curva misurata.

## Criteri comuni
- **Soglia prima dei dati.** La soglia geometrica si fissa prima di guardare i tempi; l'ottimizzata si riporta solo come controllo.
- **Appaiamento.** Ogni ciclo veloce è appaiato al suo lento; il test sulla durata totale non si usa.
- **Budget.** Il budget è ${BUDGET} per confronto; il CV massimo è ${CV_MAX_1_GEO} per giro e ${CV_MAX_3_GEO} sulle somme.
- **Referto.** Si riportano insieme valore, errore e falsificazione; un risultato fuori banda si dichiara, non si scarta.

## Che cosa questo protocollo non prova
- **Non prova la sicurezza biologica.** ${T_CRITERIO} è un criterio termico scelto, non una soglia di danno.
- **Non prova il latch.** Misura rapporto e rumore; non misura l'attuazione, l'energia o la fatica.
- **Non prova il funzionamento in flusso.** Misura la soglia in stasi; il sangue intero resta una reologia efficace, non una curva misurata.
