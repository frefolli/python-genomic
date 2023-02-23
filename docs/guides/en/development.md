# Development Guide

## Project layout

| file | description |
| :---: | :---------- |
| mkdocs.yml    | The mkdocs configuration file. |
| docs/         | mkdocs documentation files |
| lib/          | source code of library (and in fact of whole module) |
| tests/        | tests for lib |
| samples/      | sample bam/fasta files for both examples and testing |
| sample.bam    | Alignments of Drosophila Melanogaster X chromosome |
| BDGP6.X.fasta | DNA of Drosophila Melanogaster X chromosome |

### Content of lib/

| file | description |
| :---: | :---------- |
| lib/__init__.py | top module |
| lib/__main__.py | entry point |
| lib/aligned_segment.py | Aligned Segment datatype |
| lib/alignment_worker.py | Alignment Worker |
| lib/bam_file.py | BAM file abstraction |
| lib/cigar_operation.py | Cigar Operation enum |
| lib/cli.py | Command Line Interface abstraction |
| lib/fasta_file.py | FASTA file abstraction |
| lib/csv_file.py | CSV file abstraction |
| lib/intron.py | Intron datatype |
| lib/reference.py | Reference datatype |

### Content of tests/

| file | description |
| :---: | :---------- |
| tests/__init__.py | top module |
| tests/alignment_test.py | tests on AlignedSegment and AlignmentWorker |
| tests/bam_file_test.py | tests on BamFile |
| tests/cli_test.py | tests on CLI |
| tests/fasta_file_test.py | tests on FastaFile |
| tests/csv_file_test.py | tests on CsvFile |
| tests/intron_test.py | tests on Intron |

## Running Tests

After [installing test requirements](glossary.md#installing-test-requirements),
continue by [running tests](glossary.md#running-tests)

## Running Lints

After [installing test requirements](glossary.md#installing-test-requirements),
next step is [running lints](glossary.md#running-lints)

## Running SonarQube

Follow [running-sonarqube](glossary.md#running-sonarqube)