# Usage Guide

## Run

Depending on how you installed this module, either:

  - [run python module](running-python-module)
  - [run pip package](running-pip-package)

This project used as example two files, which aren't provided in repository.
You'll have to download them by yourself if you'll want to use them:

| file | description |
| :---: | :---------- |
| sample.bam    | Alignments of Drosophila Melanogaster X chromosome |
| BDGP6.X.fasta | DNA of Drosophila Melanogaster X chromosome |

By default this software will use all threads available (= ```os.cpu_count()```).
If you you wish to limit that number use `-j <job-number>` as command line option.