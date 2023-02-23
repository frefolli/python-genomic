# Guida all'Uso

## Esecuzione

A seconda di come hai installato il software:

  - [esegui il modulo python](glossary.md#running-python-module)
  - [esegui il modulo pip](glossary.md#running-pip-package)

Questo progetto ha usato come esempio due file, che non sono provvisti nel repository. Se vuoi usarli devi scaricarli:

| file | descrizione |
| :---: | :---------- |
| sample.bam    | Allineamenti del cromosoma X della Drosophila Melanogaster |
| BDGP6.X.fasta | DNA del cromosoma X della Drosophila Melanogaster |

Di default il software usera' tutti i threads disponibili (= ```os.cpu_count()```).
Se vuoi limitare questo numero usa `-j <job-number>` opzione da linea di comando.