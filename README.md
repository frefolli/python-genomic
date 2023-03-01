# Progetto di Bioinformatica

## Task

Prendere in input un file in formato BAM che contiene gli allineamenti rispetto a un reference e il file in formato FASTA che contiene il reference stesso. Si richiede di produrre in output i reads spliced, cioè i reads che inducono un introne sul reference. In output, per ogni read spliced (supponendo che ogni read allineato supporti un solo introne) devono essere prodotti:

1. il prefisso del read che si allinea al suffisso dell'esone al 5'
2. il suffisso del read che si allinea al prefisso dell'esone al 3'
3. le prime 20 basi dell'introne
4. le ultime 20 basi dell'introne. Si devono specificare in output l'identificatore del read e se l'introne supportato è canonico (cioé inizia con GT e finisce con AG)

## Sviluppo

Il progetto e' stato sviluppato da

- Refolli Francesco 865955
- Terzi Telemaco 865981

Abbiamo sviluppato una libreria (che abbiamo chiamato `genomic`) in grado di estrarre da un file BAM le reads che contengono un introne e quindi estrarre i dati su quell'introne e le relative porzioni di esoni che lo circondano nel read.

## Istruzioni

### Scaricare il sorgente

`git clone git@github.com:frefolli/bioinf-progetto.git` o tramite l'opzione di github per esportare lo zip.

### Installare i requisiti

E' richiesto `python >= 3.10`, perche' vengono usate alcune feature delle ultime versioni.

Per evitare conflitti con l'installazione globale di pip e' preferibile usare un virtual environment.
In questo caso per crearlo e' sufficiente usare `python3 -m venv .env` e quindi attivare l'environment con `source .env/bin/activate`.
Python chiedera' di aggiornare pip, in questo caso occorre usare `python3 -m pip install --upgrade pip`.

Quindi bisogna installare i pacchetti richiesti da python. La lista e' contenuta in `requirements.txt`.
Si possono installare con il comando `pip3 install -r requirements.txt`.

### Eseguire il progetto

Si dispongano nella cartella di lavoro i seguenti file:

| file | descrizione |
| :---: | :---------- |
| sample.bam    | Allineamenti del cromosoma X della Drosophila Melanogaster |
| BDGP6.X.fasta | DNA del cromosoma X della Drosophila Melanogaster |

Una volta eseguito il software con `python3 -m lib`, esso produrra' un file `output.csv` in formato CSV con le informazioni richieste dal task progetto. Produrra' anche un file intermedio `reads.csv` sempre in formato CSV con informazioni accessorie prodotte durante l'analisi dei file.

Di default il software usera' tutti i threads disponibili (= ```os.cpu_count()```).
Se si vuole limitare questo numero, si usi `-j <job-number>` come opzione da linea di comando.

Per otterere informazioni circa le opzioni da linea di comando digitare `-h`.
