#!/bin/bash

# Get directory of script (rather than directory the script is run from)
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export zzt="${DIR}"

if [ -d "$zzt"/recipe/vpython ]; then
    echo "Activating python virtual environment"
    source "$zzt"/recipe/vpython/bin/activate
fi

export PYTHONPATH="${PYTHONPATH}":"${zzt}"

echo "ZZTools setup complete"

if [ -d "$zzt"/RooUnfold-1.1.1 ]; then
    export LD_LIBRARY_PATH="$LD_LIBRARY_PATH":"$zzt"/RooUnfold-1.1.1
fi
