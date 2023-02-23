# Glossary

## Instructions

### Cloning Repository

| method | command |
| ------ | ------- |
| https | `git clone https://github.com/frefolli/bioinf-progetto.git` |
| ssh | `git clone git@github.com:frefolli/bioinf-progetto.git` |
| Github CLI | `gh repo clone frefolli/bioinf-progetto` |

### Downloading Release Source Code

After [finding release](#finding-release), simply download `Source Code` with your preferred file format: either `.tar.gz` or `.zip`.

### Downloading Release PyPi Package 

After [finding release](#finding-release), simply download `genomic-*.tar.gz`.

### Finding Release

Goto [project github page](https://github.com/frefolli/bioinf-progetto) with your browser, you will see the tags section:

![Github Page Tags](../../images/github_page_tags.png)

Then click on it and open `latest-master` tag, it will have a section called "Assets".

![Github Page Tags](../../images/github_page_latest_master_tag.png)

### Installing Module Requirements

Use `pip install -r requirements.txt`.

### Installing Build Requirements

Use `pip install build`.

### Installing Documentation Requirements

Use `pip install mkdocs "mkdocstrings[python]" mkdocs-material`.

### Installing Test Requirements

Use `pip install flake8 pylint pytest coverage`.

### Building PyPi Package

Use `python -m build` or `./actions.sh build`

### Installing PyPi Package

Use `pip install genomic-*.tar.gz`.
It will install module dependency for you.

### Running Python Module

Use `python -m lib -h`, you will see help screen:

![Github Page Tags](../../images/running_python_module_help.png)

### Running PyPi Package

Use `genomic -h`, you will see help screen:

![Github Page Tags](../../images/running_pip_package_help.png)

### Running Tests

In order to run only tests, use `pytest`.

If you wish to see coverage, `./actions.sh CI` gets the job done.

### Running Lints

With flake8:

  - `flake8 ./lib/*.py ./tests/*.py --count --select=E9,F63,F7,F82 --show-source --statistics`
  - `flake8 ./lib/*.py ./tests/*.py --count --exit-zero --max-complexity=10 --max-line-length=127 --statistics`

With pylint:

  - `pylint ./lib/*.py ./tests/*.py --exit-zero`

### Running Sonarqube

Use `./actions.sh CI` then `./actions.sh sonarqube`.

### Creating New Branch

Your new branch will be called `dev-{username}-{feature}`.
Where username is your username or optionally an acronym of it,
and feature is a short description of your contribution.

First ensure you are on master.

typing `git branch` should return something like:

```git
  ...
* master
  ...
```

If it isn't the case:

  - if `git status` says that there are changes use `git stash`
  - finally use `git checkout master`

Now pull from master to ensure you are forking from the latest version of it: 

use ```
  git pull origin master
```.

Then use `git checkout -b dev-{username}-{feature}` to create the new branch.