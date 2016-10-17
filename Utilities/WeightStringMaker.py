'''

WeightStringMaker.py

Some utilities to make strings to do event weighting
from histograms, JSONs, or C++/ROOT macros.

Author: Nate Woods, U. Wisconsin

'''

import logging
from rootpy import log as rlog; rlog = rlog["/WeightStringMaker"]
# don't show most silly ROOT messages
logging.basicConfig(level=logging.WARNING)
rlog["/rootpy.compiled"].setLevel(rlog.WARNING)

import rootpy.compiled as _comp
from rootpy.ROOT import gROOT as _gROOT
try:
    gROOT = _gROOT._gROOT
except AttributeError:
    gROOT = _gROOT

class _WeightStringSingleton(type):
    '''
    Have to do some kind of singleton thing to avoid duplicate functions.
    Will need to change if ever used multithreaded.
    '''
    _instances = {}
    def __call__(cls, fName, *args, **kwargs):
        try:
            return cls._instances[fName]
        except KeyError:
            cls._instances[fName] = super(_WeightStringSingleton, cls).__call__(fName, *args, **kwargs)
            return cls._instances[fName]

class WeightStringMaker(object):
    __metaclass__ = _WeightStringSingleton

    _counter = 0
    _hists = []
    _functions = []

    def __init__(self, fName='weightFun'):
        self.fName = fName

        self.codeBase = '''
            #include "TH1.h"
            #include "TROOT.h"
            #include <iostream>

            double {}{{0}}({{1}})
            {{{{
              TH1* h = (TH1*)gROOT->Get("{{2}}");
              if(h == 0)
                {{{{
                  std::cout << "Can't find {{2}}!" << std::endl;
                  return 0.;
                }}}}
              return h->GetBinContent(h->FindBin({{3}}));
            }}}}'''.format(self.fName)

    def makeWeightStringFromHist(self, h, *variables, **kwargs):
        '''
        Return a string that weights an event by the value of histogram h in
        the bin that would be filled by variables.
        '''
        # make a copy so we can change directory, save it in global scope
        hCopy = h.clone()
        hCopy.SetDirectory(gROOT)
        self._hists.append(hCopy)

        iName = "{0}{1}".format(self.fName, self._counter)

        _comp.register_code(
            self.codeBase.format(self._counter,
                                 ', '.join('double x%d'%i for i in range(len(variables))),
                                 hCopy.GetName(),
                                 ', '.join("x%d"%i for i in range(len(variables)))),
            [iName,])
        out = '{0}({1})'.format(iName, ', '.join(variables))

        # Save function. Also forces system to compile code now
        self._functions.append(getattr(_comp, iName))

        self._counter += 1

        return out

    def getWeightFunction(self, i):
        '''
        Return a function object for the ith weight function made.
        '''
        return self._functions[i]
