#!/bin/bash
####
####    Set up a virtual environment with Rootpy, IPython, and Jupyter.
####    Installs things the first time it runs, just activates the virtualenv
####    thereafter.
####
####    Rootpy: Pythonic wrappers for PyROOT (using v0.8.3 because 0.9.0 is broken somehow)
####    IPython: Nicer interactive Python
####    Jupyter: IPython notebooks
####    BRIL Work Suite: tools for luminosity (etc.)
####
####    Author: Nate Woods, U. Wisconsin
####


if [ "$zzt" == "" ]; then
    echo "Please set up this package with setup.sh before running this script."
    exit 2
fi

export vpython="$zzt"/recipe/vpython

if [ ! -d "$vpython" ]; then
    python "$zzt"/recipe/virtualenv/virtualenv.py "$vpython"

    source "$vpython"/bin/activate

    pip install -U pip
    pip install -U rootpy==0.8.3
    pip install -U ipython
    pip install -U jupyter
    if [ "$CMSSW_BASE" != "" ]; then
        pip install brilws
    fi
else
    source "$vpython"/bin/activate
fi

export PYTHONPATH="$vpython"/lib/python2.7/site-packages/:"$PYTHONPATH"

