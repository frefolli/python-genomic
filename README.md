# bioinf-progetto

## Istruzioni

### Scaricare il sorgente

`git clone git@github.com:frefolli/bioinf-progetto.git` o tramite l'opzione di github per esportare lo zip.

### Installare i requisiti

E' richiesto `python >= 3.10`, perche' vengono usate alcune feature delle ultime versioni.

Per evitare conflitti con l'installazione globale di pip e' preferibile usare un virtual environment.
In questo caso per crearlo e' sufficiente usare `python3 -m venv .env` e quindi attivare l'environment con `source .env/bin/activate`.
Python chiedera' di aggiornare pip, in questo caso occorre usare `python3 -m pip install --upgrade pip`.

Quindi bisogna installare i pacchetti richiesti da python. La lista e' contenuta in `requirements.txt`. Si possono installare con il comando `pip3 install -r requirements.txt`.

### Eseguire il progetto

Si dispongano nella cartella di lavoro i seguenti file:

| file | descrizione |
| :---: | :---------- |
| sample.bam    | Allineamenti del cromosoma X della Drosophila Melanogaster |
| BDGP6.X.fasta | DNA del cromosoma X della Drosophila Melanogaster |

Di default il software usera' tutti i threads disponibili (= ```os.cpu_count()```).
Se si vuole limitare questo numero, si usi `-j <job-number>` opzione da linea di comando.

Per otterere informazioni circa le opzioni da linea di comando digitare `-h`.

Una volta eseguito il software con `python3 -m lib`, esso produrra' un file in formato CSV con le informazioni richieste dal task progetto. Produrra' anche un file intermedio sempre in formato CSV con informazioni accessorie prodotte durante l'analisi dei file.
