# Il protocollo di misura V7.3

Questa cartella contiene la scheda da banco del protocollo V7.3: una pagina,
in italiano e in inglese, con tre setup di misura. Non aggiunge fisica e non
modifica il preprint; mette in fila le misure che cambierebbero lo stato della
V7.3 e dichiara, per ognuna, valore atteso, errore atteso e criterio di
falsificazione.

## File

| File | Contenuto |
|------|-----------|
| `IT.md`, `EN.md` | Protocollo renderizzato in Markdown |
| `IT.tex`, `EN.tex` | Sorgente LaTeX del PDF |
| `protocollo_IT.pdf`, `protocollo_EN.pdf` | PDF di una pagina |
| `templates/IT.md.tpl`, `templates/EN.md.tpl` | Template con i segnaposto `${MACRO}` |

Il codice che produce tutto questo sta in `../codice/`:

| Script | Ruolo |
|--------|-------|
| `v73_chiavi.py` | Inventario delle chiavi numeriche e registro delle macro |
| `v73_protocollo.py` | Payload del protocollo e specifica del latch |
| `v73_render_protocollo.py` | Markdown, LaTeX, PDF e controllo di una pagina |

I dati stanno in `../dati/`:

| File | Contenuto |
|------|-----------|
| `V73_chiavi.json` | Inventario delle foglie numeriche e registro delle macro |
| `V73_ig_verifica.json` | Verifica preliminare della inverse Gaussian e Monte Carlo |
| `V73_protocollo.json` | Macro usate dai template e sezioni strutturate |
| `V73_latch_spec.json` | CV massimo per i quattro rapporti, per giro e sulle somme |

## I tre setup

1. **Latch.** Rapporto tra mezzo-giro lento e mezzo-giro veloce, con test
   appaiato 3+3. Il valore atteso è 3.881 (fino a 4.492 con la correzione
   idrodinamica non validata). L'errore per giro è 1.07e-4; quello sulle somme
   di 3 è 1.95e-11. La falsificazione è un rapporto sotto 2 o fuori banda, o un
   CV misurato sopra il limite del budget.
2. **Termica.** ΔT in sangue a 800 e 1064 nm. La forcella pozzo-specchio è
   0.347–0.670 K a 800 nm e 0.283–0.499 K a 1064 nm. Il criterio è 5 K, scelto
   e non biologico. La falsificazione è ΔT ≥ 5 K.
3. **Cancello.** Soglia di shear in plasma (2.700 s⁻¹, banda 2.340–2.700) e in
   sangue (0.104–0.409 s⁻¹ con Bingham). Comprende la fuga a Q = 0.1:
   U₀ = 8.59 kT, τesc = 200.7 s, con la convenzione τD = ξ/krad e non
   w₀²/(4D). La falsificazione è rotazione persistente a 267 s⁻¹ o nessuna
   fuga in 10τesc (2007 s) a Q = 0.1.

## Il ponte 17.76 / 17.83

Il limite 17.76% si riferisce alla soglia geometrica; il 17.83% del preprint è
quello della soglia ottimizzata, qui riportata solo come controllo. La
geometrica è più severa e si fissa prima dei dati: per questo il protocollo usa
17.76% nella falsificazione, mentre 17.83% serve solo al confronto con il
preprint. La distinzione è sviluppata nell'addendum al preprint:
`../paper/it/Addendum_V73_it.pdf` (IT) e `../paper/en/Addendum_V73_en.pdf` (EN).

## Regola dei numeri

Ogni numero stampato nel protocollo viene da una macro in `../dati/V73_*.json`.
I template non contengono cifre: contengono solo prosa e segnaposto
`${MACRO}`. Una macro mancante solleva `KeyError` con il nome; non esiste un
valore di default. Il test di regressione verifica che i template IT/EN abbiano
la stessa sequenza di macro, lo stesso numero di righe e nessuna cifra fuori
dai segnaposto.

## Rigenerare

Dalla cartella `versione 7.3`:

```powershell
python codice\v73_protocollo.py
python codice\v73_render_protocollo.py
```

`v73_protocollo.py` richiama `v73_chiavi.py`, quindi l'inventario e la verifica
IG si aggiornano prima del payload. Il renderer compila con `pdflatex` e
verifica con `pypdf` che ogni PDF abbia esattamente una pagina.
