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
    $ENV/bin/python3 -m main
}

function run_sonarqube() {
    sonar-scanner
}

function run_tests() {
    $ENV/bin/pytest --cov=. --junit-xml=coverage.xml tests/
}

case "$1" in
    'run')
        run_main
        ;;
    'verify')
        verify_python_env
        ;;
    'sample')
        run_tests
        ;;
    'sonarqube')
        run_sonarqube
        ;;
    'test')
        run_tests
        ;;
    *)
        echo "Usage: $0 { run | verify | test | sonarqube | sample }"
        exit 1
        ;;
esac

exit 0
