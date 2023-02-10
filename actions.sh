#!/bin/bash
ENV=env

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

function main() {
    verify_python_env
    $ENV/bin/python3 -m main
}

main