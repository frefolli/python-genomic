#!/bin/bash
ENV=.env

function install_python_env() {
    echo -e "installing python enviroment"
    python3 -m venv $ENV
}

function install_python_requirements() {
    echo -e "asserting python dependencies"
    $ENV/bin/python3 -m pip install --upgrade pip
    $ENV/bin/pip3 install -r requirements.txt
}

function verify_python_env() {
    if [ ! -f $ENV/bin/python3 ]; then
        install_python_env
    fi
    install_python_requirements
}

function run_main() {
    $ENV/bin/python3 -m lib
}

function run_sonarqube() {
    sonar-scanner
}

function run_tests() {
    $ENV/bin/coverage run -m pytest
    $ENV/bin/coverage xml
}

function run_lints() {
    $ENV/bin/pylint ./lib/*.py ./tests/*.py --exit-zero
    $ENV/bin/flake8 ./lib/*.py ./tests/*.py --count --select=E9,F63,F7,F82 --show-source --statistics
    $ENV/bin/flake8 ./lib/*.py ./tests/*.py --count --exit-zero --max-complexity=10 --max-line-length=127 --statistics
}

case "$1" in
    'run')
        run_main
        ;;
    'verify')
        verify_python_env
        ;;
    'CI')
        run_tests
        run_sonarqube
        ;;
    'sonarqube')
        run_sonarqube
        ;;
    'test')
        run_tests
        ;;
    'lint')
        run_lints
        ;;
    *)
        echo "Usage: $0 { run | verify | test | sonarqube | CI | lint }"
        exit 1
        ;;
esac

exit 0
