# Guida allo Sviluppo

## Creazione dell'Ambiente di Sviluppo

Dopo aver [clonato il repository](glossary.md#cloning-repository),
non resta che [creare una nuova branch](glossary.md#creating-new-branch) e
quindi [installare i requisiti del modulo](glossary.md#installing-module-requirements)

## Struttura del Progetto

| file | descrizione |
| :---: | :---------- |
| mkdocs.yml    | file di configurazione per mkdocs |
| docs/         | file di documentazione |
| lib           | codice sorgente del modulo |
| tests         | vari test sul modulo |
| samples       | vari file semplici bam/fasta per testing e esempi |
| sample.bam    | Allineamenti del cromosoma X della Drosophila Melanogaster |
| BDGP6.X.fasta | DNA del cromosoma X della Drosophila Melanogaster |

### Contenuto of lib/

| file | descrizione |
| :---: | :---------- |
| lib/__init__.py | indice di modulo |
| lib/__main__.py | punto d'entrata del modulo |
| lib/aligned_segment.py | classe Aligned Segment |
| lib/alignment_worker.py | classe Alignment Worker |
| lib/bam_file.py | astrazione BAM file |
| lib/cigar_operation.py | enumerazione Cigar Operation |
| lib/cli.py | astrazione Command Line Interface |
| lib/fasta_file.py | astrazione FASTA file |
| lib/csv_file.py | astrazione CSV file |
| lib/intron.py | classe Intron |
| lib/reference.py | classe Reference |

### Contenuto of tests/

| file | descrizione |
| :---: | :---------- |
| tests/__init__.py | indice di modulo |
| tests/alignment_test.py | test su AlignedSegment e AlignmentWorker |
| tests/bam_file_test.py | test su BamFile |
| tests/cli_test.py | test su CLI |
| tests/fasta_file_test.py | test su FastaFile |
| tests/csv_file_test.py | test su CsvFile |
| tests/intron_test.py | test su Intron |

## Esecuzione dei Test

Dopo aver [installato i requisiti per i test](glossary.md#installing-test-requirements),
continua [eseguendo i test](glossary.md#running-tests)

## Esecuzione dei Lint

Dopo aver [installato i requisiti per i test](glossary.md#installing-test-requirements),
il prossimo passo e' [eseguire i lint](glossary.md#running-lints)

## Esecuzione di SonarQube

Segui le istruzioni di [eseguire sonarqube](glossary.md#running-sonarqube)