#!/bin/bash

# Get directory of script (rather than directory the script is run from)
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export zzt="${DIR}"

if [ ! -d "$zzt"/recipe/cmssw ]; then
    pushd "$zzt"/recipe
    scram pro -n cmssw CMSSW CMSSW_8_0_26
    popd
fi
pushd "$zzt"/recipe/cmssw/src
cmsenv
popd

if [ -d "$zzt"/recipe/vpython ]; then
    echo "Activating python virtual environment"
    source "$zzt"/recipe/vpython/bin/activate
fi

export PYTHONPATH="${zzt}":"${PYTHONPATH}"

echo "ZZTools setup complete"

if [ -d "$zzt"/RooUnfold-1.1.1 ]; then
    export LD_LIBRARY_PATH="$LD_LIBRARY_PATH":"$zzt"/RooUnfold-1.1.1
fi
if [ -d "$zzt"/RooUnfold ]; then
    export LD_LIBRARY_PATH="$LD_LIBRARY_PATH":"$zzt"/RooUnfold
fi

