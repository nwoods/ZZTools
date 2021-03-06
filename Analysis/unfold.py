
import logging
from rootpy import log as rlog; rlog = rlog["/unfold"]
# don't show most silly ROOT messages
logging.basicConfig(level=logging.WARNING)
rlog["/ROOT.TUnixSystem.SetDisplay"].setLevel(rlog.ERROR)

from rootpy import asrootpy
from rootpy.io import root_open
from rootpy.plotting import Canvas, Legend, Hist, Hist2D, HistStack, Graph
from rootpy.plotting.utils import draw
from rootpy.ROOT import cout, TDecompSVD, TBox, TLatex
from rootpy.ROOT import RooUnfoldBayes as RooUnfoldIter # it's frequentist!

from SampleTools import MCSample, DataSample, SampleGroup, SampleStack
from PlotTools import PlotStyle as _Style
from PlotTools import makeLegend, addPadsBelow, makeRatio, fixRatioAxes, makeErrorBand
from Utilities import WeightStringMaker, identityFunction, Z_MASS, \
    deltaRFunction, deltaRString, deltaPhiFunction, deltaPhiString
from Analysis.setupStandardSamples import *
from Analysis.unfoldingHelpers import getResponse, getResponsePDFErrors, \
    getResponseScaleErrors, getResponseAlphaSErrors
from Analysis.weightHelpers import puWeight, baseMCWeight
from Metadata.metadata import sampleInfo

from os import environ as _env
from os import makedirs as _mkdir
from os.path import join as _join
from os.path import isdir as _isdir
from os.path import exists as _exists
from math import sqrt


_channels = ['eeee','eemm', 'mmmm']

def _makeReturner(var, doAbs=False):
    '''
    Make a little function to return var (or its absolute value) from an
    ntuple row
    '''
    newVar = str(var) # insurance against stupid scoping issues
    if doAbs:
        return lambda row: abs(getattr(row, newVar))
    else:
        return lambda row: getattr(row, newVar)

def _makeVectorReturner(var, n, doAbs):
    if doAbs:
        return lambda row: abs(getattr(row, var).at(n-1))

    return lambda row: getattr(row, var).at(n-1)

def _makeComparator(var, compareTo):
    '''
    Make a little function to check whether var is greater than compareTo
    '''
    # insurance against stupid scoping issues
    newVar = str(var)
    newComp = type(compareTo)(compareTo)

    return lambda row: getattr(row, newVar) > newComp

_wrongZRejectionTemp = ('(({0}1Charge=={0}3Charge||abs(_getInvariantMass({0}1Pt,{0}1Eta,{0}1Phi,{0}3Pt,{0}3Eta,{0}3Phi)-{1}) > abs({0}1_{0}2_Mass-{1})) && '
                        ' ({0}1Charge=={0}4Charge||abs(_getInvariantMass({0}1Pt,{0}1Eta,{0}1Phi,{0}4Pt,{0}4Eta,{0}4Phi)-{1}) > abs({0}1_{0}2_Mass-{1})) && '
                        ' ({0}2Charge=={0}3Charge||abs(_getInvariantMass({0}2Pt,{0}2Eta,{0}2Phi,{0}3Pt,{0}3Eta,{0}3Phi)-{1}) > abs({0}1_{0}2_Mass-{1})) && '
                        ' ({0}2Charge=={0}4Charge||abs(_getInvariantMass({0}2Pt,{0}2Eta,{0}2Phi,{0}4Pt,{0}4Eta,{0}4Phi)-{1}) > abs({0}1_{0}2_Mass-{1})))')
_wrongZRejectionStr = {
    'eeee' : _wrongZRejectionTemp.format('e', Z_MASS),
    'mmmm' : _wrongZRejectionTemp.format('m', Z_MASS),
    'eemm' : '1'
    }


# set up variables, selection, binnings etc.
# (jet-related variables and selections done later)
_variables = {
    'pt' : {c:'Pt' for c in _channels},
    'mass' : {c:'Mass' for c in _channels},
    'z1Mass' : {'eeee':'e1_e2_Mass', 'mmmm':'m1_m2_Mass',
                'eemm':['e1_e2_Mass','m1_m2_Mass']},
    'z2Mass' : {'eeee':'e3_e4_Mass', 'mmmm':'m3_m4_Mass',
                'eemm':['m1_m2_Mass','e1_e2_Mass']},
    'z1Pt' : {'eeee':'e1_e2_Pt', 'mmmm':'m1_m2_Pt',
              'eemm':['e1_e2_Pt','m1_m2_Pt']},
    'z2Pt' : {'eeee':'e3_e4_Pt', 'mmmm':'m3_m4_Pt',
              'eemm':['m1_m2_Pt','e1_e2_Pt']},
    'zHigherPt' : {
        'eeee' : ['e1_e2_Pt', 'e3_e4_Pt'],
        'mmmm' : ['m1_m2_Pt', 'm3_m4_Pt'],
        'eemm' : ['e1_e2_Pt','m1_m2_Pt']
        },
    'zLowerPt' : {
        'eeee' : ['e1_e2_Pt', 'e3_e4_Pt'],
        'mmmm' : ['m1_m2_Pt', 'm3_m4_Pt'],
        'eemm' : ['e1_e2_Pt','m1_m2_Pt']
        },
    'zPt' : {
        'eeee' : ['e1_e2_Pt', 'e3_e4_Pt'],
        'mmmm' : ['m1_m2_Pt', 'm3_m4_Pt'],
        'eemm' : ['e1_e2_Pt','m1_m2_Pt']
        },
    'deltaPhiZZ' : {
        'eeee' : 'abs({}(e1_e2_Phi, e3_e4_Phi))'.format(deltaPhiString()),
        'eemm' : 'abs({}(e1_e2_Phi, m1_m2_Phi))'.format(deltaPhiString()),
        'mmmm' : 'abs({}(m1_m2_Phi, m3_m4_Phi))'.format(deltaPhiString()),
        },
    'deltaRZZ' : {
        'eeee' : '{}(e1_e2_Eta, e1_e2_Phi, e3_e4_Eta, e3_e4_Phi)'.format(deltaRString()),
        'eemm' : '{}(e1_e2_Eta, e1_e2_Phi, m1_m2_Eta, m1_m2_Phi)'.format(deltaRString()),
        'mmmm' : '{}(m1_m2_Eta, m1_m2_Phi, m3_m4_Eta, m3_m4_Phi)'.format(deltaRString()),
        },
    'lPt' : {
        'eeee' : ['e1Pt', 'e2Pt', 'e3Pt', 'e4Pt'],
        'eemm' : ['e1Pt', 'e2Pt', 'm1Pt', 'm2Pt'],
        'mmmm' : ['m1Pt', 'm2Pt', 'm3Pt', 'm4Pt'],
        },
    'l1Pt' : {
        'eeee' : 'max(e1Pt, e3Pt)',
        'eemm' : 'max(e1Pt, m1Pt)',
        'mmmm' : 'max(m1Pt, m3Pt)',
        },
    }
_fDPhi = deltaPhiFunction()
_fDR = deltaRFunction()
_varFunctions = {
    'pt' : {c:lambda row: row.Pt for c in _channels},
    'mass' : {c:lambda row: row.Mass for c in _channels},
    'z1Mass' : {
        'eeee' : lambda row: row.e1_e2_Mass,
        'mmmm' : lambda row: row.m1_m2_Mass,
        'eemm' : lambda row: row.e1_e2_Mass if abs(row.e1_e2_Mass - Z_MASS) < abs(row.m1_m2_Mass - Z_MASS) else row.m1_m2_Mass,
        },
    'z2Mass' : {
        'eeee' : lambda row: row.e3_e4_Mass,
        'mmmm' : lambda row: row.m3_m4_Mass,
        'eemm' : lambda row: row.e1_e2_Mass if abs(row.e1_e2_Mass - Z_MASS) > abs(row.m1_m2_Mass - Z_MASS) else row.m1_m2_Mass,
        },
    'z1Pt' : {
        'eeee' : lambda row: row.e1_e2_Pt,
        'mmmm' : lambda row: row.m1_m2_Pt,
        'eemm' : lambda row: row.e1_e2_Pt if abs(row.e1_e2_Mass - Z_MASS) < abs(row.m1_m2_Mass - Z_MASS) else row.m1_m2_Pt,
        },
    'z2Pt' : {
        'eeee' : lambda row: row.e3_e4_Pt,
        'mmmm' : lambda row: row.m3_m4_Pt,
        'eemm' : lambda row: row.e1_e2_Pt if abs(row.e1_e2_Mass - Z_MASS) > abs(row.m1_m2_Mass - Z_MASS) else row.m1_m2_Pt,
        },
    'zPt' : {
        'eeee' : [lambda row: row.e1_e2_Pt, lambda row: row.e3_e4_Pt],
        'eemm' : [lambda row: row.e1_e2_Pt, lambda row: row.m1_m2_Pt],
        'mmmm' : [lambda row: row.m1_m2_Pt, lambda row: row.m3_m4_Pt],
        },
    'zHigherPt' : {
        'eeee' : lambda row: max(row.e1_e2_Pt, row.e3_e4_Pt),
        'eemm' : lambda row: max(row.e1_e2_Pt, row.m1_m2_Pt),
        'mmmm' : lambda row: max(row.m1_m2_Pt, row.m3_m4_Pt),
        },
    'zLowerPt' : {
        'eeee' : lambda row: min(row.e1_e2_Pt, row.e3_e4_Pt),
        'eemm' : lambda row: min(row.e1_e2_Pt, row.m1_m2_Pt),
        'mmmm' : lambda row: min(row.m1_m2_Pt, row.m3_m4_Pt),
        },
    'deltaPhiZZ' : {
        'eeee' : lambda row: abs(_fDPhi(row.e1_e2_Phi, row.e3_e4_Phi)),
        'eemm' : lambda row: abs(_fDPhi(row.e1_e2_Phi, row.m1_m2_Phi)),
        'mmmm' : lambda row: abs(_fDPhi(row.m1_m2_Phi, row.m3_m4_Phi)),
        },
    'deltaRZZ' : {
        'eeee' : lambda row: _fDR(row.e1_e2_Eta, row.e1_e2_Phi, row.e3_e4_Eta, row.e3_e4_Phi),
        'eemm' : lambda row: _fDR(row.e1_e2_Eta, row.e1_e2_Phi, row.m1_m2_Eta, row.m1_m2_Phi),
        'mmmm' : lambda row: _fDR(row.m1_m2_Eta, row.m1_m2_Phi, row.m3_m4_Eta, row.m3_m4_Phi),
        },
    'lPt' : {
        'eeee' : [lambda row: row.e1Pt, lambda row: row.e2Pt, lambda row: row.e3Pt, lambda row: row.e4Pt],
        'eemm' : [lambda row: row.e1Pt, lambda row: row.e2Pt, lambda row: row.m1Pt, lambda row: row.m2Pt],
        'mmmm' : [lambda row: row.m1Pt, lambda row: row.m2Pt, lambda row: row.m3Pt, lambda row: row.m4Pt],
        },
    'l1Pt' : {
        'eeee' : lambda row: max(row.e1Pt, row.e3Pt),
        'eemm' : lambda row: max(row.e1Pt, row.m1Pt),
        'mmmm' : lambda row: max(row.m1Pt, row.m3Pt),
        },
    }

_binning = {
    'pt' : [25.*i for i in range(4)] + [100., 150., 200., 300.],
    'nJets' : [6,-0.5,5.5],
    'mass' : [100.] + [200.+50.*i for i in range(5)] + [500.,600.,800.],
    'jet1Pt' : [0., 50., 100., 200., 300., 500.],
    'jet1Eta' : [0., 1.5, 3., 4.7],
    'jet2Pt' : [30., 100., 200., 500.],
    'jet2Eta' : [0., 1.5, 3., 4.7],
    'mjj' : [0., 100., 300., 800.],
    'deltaEtajj' : [6, 0.,6.],
    'z1Mass' : [60., 80., 84., 86.] + [87.+i for i in range(10)] + [98., 102., 120.], #[12, 60., 120.],
    'z2Mass' : [60., 75., 83.] + [84.+i for i in range(14)] + [105., 120.],#[12, 60., 120.],
    'z1Pt' : [i * 25. for i in range(7)] + [200., 300.],
    'z2Pt' : [i * 25. for i in range(7)] + [200., 300.],
    'zPt' : [i * 25. for i in range(7)] + [200., 300.],
    'zHigherPt' : [i * 25. for i in range(7)] + [200., 300.],
    'zLowerPt' : [i * 25. for i in range(7)] + [200., 300.],
    'deltaPhiZZ' : [0., 1.5] + [2.+.25*i for i in range(6)],
    'deltaRZZ' : [6, 0., 6.],
    'lPt' : [15, 0., 150.],
    'l1Pt' : [15, 0., 150.],
    }

_units = {
    'pt' : 'GeV',
    'nJets' : '',
    'mass' : 'GeV',
    'jet1Pt' : 'GeV',
    'jet1Eta' : '',
    'jet2Pt' : 'GeV',
    'jet2Eta' : '',
    'mjj' : 'GeV',
    'deltaEtajj' : '',
    'z1Mass' : 'GeV',
    'z2Mass' : 'GeV',
    'zPt' : 'GeV',
    'z1Pt' : 'GeV',
    'z2Pt' : 'GeV',
    'zHigherPt' : 'GeV',
    'zLowerPt' : 'GeV',
    'deltaPhiZZ' : '',
    'deltaRZZ' : '',
    'lPt' : 'GeV',
    'l1Pt' : 'GeV',
    }

_prettyVars = {
    'pt' : 'p_T^{4\\ell}',
    'nJets' : '# jets',
    'mass' : 'm_{4\\ell}',
    'jet1Pt' : 'p_T^\\text{j1}',
    'jet1Eta' : '\\eta_\\text{j1}',
    'jet2Pt' : 'p_T^\\text{j2}',
    'jet2Eta' : '\\eta_\\text{j2}',
    'mjj' : 'm_\\text{jj}',
    'deltaEtajj' : '|\\Delta \\eta_{\\text{jj}}|',
    'z1Mass' : 'm_{\\text{Z}_{1}}',
    'z2Mass' : 'm_{\\text{Z}_{2}}',
    'z1Pt' : 'p_T^{\\text{Z}_{1}}',
    'z2Pt' : 'p_T^{\\text{Z}_{2}}',
    'zPt' : 'p_T^{\\text{Z}}',
    'zHigherPt' : 'p_T^{\\text{Z}_{\\text{lead}}}',
    'zLowerPt' : 'p_T^{\\text{Z}_{\\text{sublead}}}',
    'deltaPhiZZ' : '\\Delta \\phi_{\\text{Z}_1,\\text{Z}_2}',
    'deltaRZZ' : '\\Delta \\text{R}_{\\text{Z}_1,\\text{Z}_2}',
    'lPt' : 'p_{T}^{\\ell}',
    'l1Pt' : 'p_{T}^{\\ell_1}',
    }

_xTitle = {}
_yTitle = {}
_yTitleTemp = '\\frac{{1}}{{\\sigma_{{\\text{{fid}}}}}} \\frac{{d\\sigma_{{\\text{{fid}}}}}}{{d{xvar}}} {units}'
for var, prettyVar in _prettyVars.iteritems():
    xt = prettyVar
    if _units[var]:
        xt += ' \\, \\text{{[{}]}}'.format(_units[var])
        yt = _yTitleTemp.format(xvar=prettyVar,
                                units='\\, \\left[ \\frac{{1}}{{\\text{{{unit}}}}} \\right]'.format(unit=_units[var]))
    else:
        yt = _yTitleTemp.format(xvar=prettyVar, units='')

    _xTitle[var] = xt
    _yTitle[var] = yt

_selections = {
    'pt' : {c:'' for c in _channels},
    'mass' : {c:'' for c in _channels},
    'z1Mass' : {'eeee':'','mmmm':'',
                'eemm':['abs(e1_e2_Mass - {0}) < abs(m1_m2_Mass - {0})'.format(Z_MASS),
                        'abs(e1_e2_Mass - {0}) > abs(m1_m2_Mass - {0})'.format(Z_MASS)],
                },
    'z2Mass' : {'eeee':'','mmmm':'',
                'eemm':['abs(e1_e2_Mass - {0}) < abs(m1_m2_Mass - {0})'.format(Z_MASS),
                        'abs(e1_e2_Mass - {0}) > abs(m1_m2_Mass - {0})'.format(Z_MASS)],
                },
    'z1Pt' : {'eeee':'','mmmm':'',
              'eemm':['abs(e1_e2_Mass - {0}) < abs(m1_m2_Mass - {0})'.format(Z_MASS),
                      'abs(e1_e2_Mass - {0}) > abs(m1_m2_Mass - {0})'.format(Z_MASS)],
              },
    'z2Pt' : {'eeee':'','mmmm':'',
              'eemm':['abs(e1_e2_Mass - {0}) < abs(m1_m2_Mass - {0})'.format(Z_MASS),
                      'abs(e1_e2_Mass - {0}) > abs(m1_m2_Mass - {0})'.format(Z_MASS)],
              },
    'zPt' : {c:'' for c in _channels},
    'zHigherPt' : {
        'eeee' : ['e1_e2_Pt > e3_e4_Pt', 'e1_e2_Pt < e3_e4_Pt'],
        'mmmm' : ['m1_m2_Pt > m3_m4_Pt', 'm1_m2_Pt < m3_m4_Pt'],
        'eemm' : ['e1_e2_Pt > m1_m2_Pt', 'e1_e2_Pt < m1_m2_Pt'],
        },
    'zLowerPt' : {
        'eeee' : ['e1_e2_Pt < e3_e4_Pt', 'e1_e2_Pt > e3_e4_Pt'],
        'mmmm' : ['m1_m2_Pt < m3_m4_Pt', 'm1_m2_Pt > m3_m4_Pt'],
        'eemm' : ['e1_e2_Pt < m1_m2_Pt', 'e1_e2_Pt > m1_m2_Pt'],
        },
    'deltaPhiZZ' : {c:'' for c in _channels},
    'deltaRZZ' : {c:'' for c in _channels},
    'lPt' : {c:'' for c in _channels},
    'l1Pt' : {c:'' for c in _channels},
    }
_selFunctions = {
    'pt' : identityFunction,
    'mass' : identityFunction,
    'z1Mass' : identityFunction, # Z choice done in var function
    'z2Mass' : identityFunction,
    'z1Pt' : identityFunction, # Z choice done in var function
    'z2Pt' : identityFunction,
    'zPt' : identityFunction,
    'zHigherPt' : identityFunction, # Z choice done in var function
    'zLowerPt' : identityFunction,
    'deltaPhiZZ' : identityFunction,
    'deltaRZZ' : identityFunction,
    'lPt' : identityFunction,
    'l1Pt' : identityFunction,
    }

# do jet variables separately because we have to deal with systematics
for sys in ['', '_jerUp', '_jerDown', '_jesUp','_jesDown']:
    for varName in ['nJets', 'mjj', 'deltaEtajj']:
        doAbs = 'eta' in varName.lower()

        varName += sys
        var = varName
        if doAbs:
            var = 'abs({})'.format(var)

        _variables[varName] = {c:var for c in _channels}
        _varFunctions[varName] = {c:_makeReturner(varName, doAbs) for c in _channels}

        if 'jj' in varName.lower():
            _selections[varName] = {c:'nJets{} > 1'.format(sys) for c in _channels}
            _selFunctions[varName] = _makeComparator('nJets'+sys, 1)
        else:
            _selections[varName] = {c:'' for c in _channels}
            _selFunctions[varName] = identityFunction

    for baseVar in ['jetPt'+sys, 'jetEta'+sys]:
        for j in [1,2]:
            doAbs = 'eta' in baseVar.lower()

            var = baseVar+'[{}]'.format(j-1)
            if doAbs:
                var = 'abs({})'.format(var)
            varName = baseVar.replace('jet','jet{}'.format(j))

            _variables[varName] = {c:var for c in _channels}
            _varFunctions[varName] = {c:_makeVectorReturner(baseVar, j, doAbs) for c in _channels}
            _selections[varName] = {c:'nJets{} >= {}'.format(sys,j) for c in _channels}
            _selFunctions[varName] = _makeComparator('nJets'+sys, j-1)

# some variables need to be blinded
_blind = {
    'mass' : 500.,
    'z1Pt' : 200.,
    'z2Pt' : 200.,
    'zPt' : 200.,
    'zHigherPt' : 200.,
    'zLowerPt' : 200.,
    }

# list of variables not counting systematic shifts
_varList = [v for v in _variables if 'Up' not in v and 'Down' not in v]

# Sometimes need to more or resize legend
_legDefaults = {
    'textsize' : .023,
    'leftmargin' : 0.3,
    }
_legParams = {v:_legDefaults.copy() for v in _varList}
_legParams['z1Mass'] = {
    'textsize' : .015,
    'leftmargin' : .03,
    'rightmargin' : .5,
    'entryheight' : .023,
    'entrysep' : .007,
    }
_legParams['z2Mass'] = _legParams['z1Mass'].copy()
_legParams['deltaRZZ'] = _legParams['z1Mass'].copy()
_legParams['deltaPhiZZ']['leftmargin'] = 0.05
_legParams['deltaPhiZZ']['rightmargin'] = 0.32
_legParams['deltaEtajj'] = _legParams['z1Mass'].copy()
_legParams['deltaEtajj']['leftmargin'] = .5
_legParams['deltaEtajj']['rightmargin'] = .03
_legParams['deltaEtajj']['topmargin'] = .05
_legParams['lPt']['topmargin'] = 0.05
_legParams['l1Pt']['topmargin'] = 0.05

_systSaveFileTemplate = _join(_env['zzt'], 'Analysis', 'savedResults',
                              'unfoldSave_{}Iter.root')

def _normalizeBins(h):
    binUnit = 1 # min(h.GetBinWidth(b) for b in range(1,len(h)+1))
    for ib in xrange(1,len(h)+1):
        w = h.GetBinWidth(ib)
        h.SetBinContent(ib, h.GetBinContent(ib) * binUnit / w)
        h.SetBinError(ib, h.GetBinError(ib) * binUnit / w)
        if h.GetBinError(ib) > h.GetBinContent(ib):
            h.SetBinError(ib, h.GetBinContent(ib))
    h.sumw2()


def _compatibleHistOrNone(fileOrDir, name, desiredBinning):
    try:
        h = getattr(fileOrDir, name)
        h.check_compatibility(desiredBinning, True)
        return h
    except:
        # lots of ways to raise errors, but all of them should return None
        return None


def main(inData, inMC, plotDir, fakeRateFile, puWeightFile, lumi, nIter, redo,
         systSaveFile, amcatnlo=False, *varNames, **kwargs):

    style = _Style()

    systSaveFile.cd()

    channels = _channels

    puWeightStr, puWt = puWeight(puWeightFile, '')
    puWeightStrUp, puWtUp = puWeight(puWeightFile, 'up')
    puWeightStrDn, puWtDn = puWeight(puWeightFile, 'dn')
    fPUWt = {'':puWt,'up':puWtUp,'dn':puWtDn}

    true = genZZSamples('zz', inMC, 'smp', lumi, amcatnlo=amcatnlo)
    reco = zzStackSignalOnly('zz', inMC, 'smp', puWeightFile,
                             lumi, amcatnlo=amcatnlo, asGroup=True)
    bkg = standardZZBkg('zz', inData, inMC, 'smp', puWeightFile,
                        fakeRateFile, lumi)
    bkgSyst = {
        'eup' : standardZZBkg('zz', inData, inMC, 'smp', puWeightFile,
                              fakeRateFile, lumi, eFakeRateSyst='up'),
        'edn' : standardZZBkg('zz', inData, inMC, 'smp', puWeightFile,
                              fakeRateFile, lumi, eFakeRateSyst='dn'),
        'mup' : standardZZBkg('zz', inData, inMC, 'smp', puWeightFile,
                              fakeRateFile, lumi, mFakeRateSyst='up'),
        'mdn' : standardZZBkg('zz', inData, inMC, 'smp', puWeightFile,
                              fakeRateFile, lumi, mFakeRateSyst='dn'),
        }

    data = standardZZData('zz', inData, 'smp')

    altReco = zzStackSignalOnly('zz', inMC, 'smp', puWeightFile, lumi,
                                amcatnlo=(not amcatnlo), asGroup=True)
    altTrue = genZZSamples('zz', inMC, 'smp', lumi,
                           amcatnlo=(not amcatnlo))

    if amcatnlo:
        signalName = 'MG5_aMC@NLO'
        signalNameAlt = 'POWHEG'
    else:
        signalName = 'POWHEG'
        signalNameAlt = 'MG5_aMC@NLO'
    signalName += '+MCFM'
    signalNameAlt += '+MCFM'

    recoSyst = {}
    for syst in ['eScaleUp', 'eScaleDn', 'eRhoResUp',
                 'eRhoResDn', 'ePhiResUp']:
        recoSyst[syst] = zzStackSignalOnly('eeee,eemm', inMC.replace('mc_','mc_{}_'.format(syst)),
                                           'smp', puWeightFile, lumi, amcatnlo=amcatnlo,
                                           asGroup=True)
    for syst in ['mClosureUp','mClosureDn']:
        recoSyst[syst] = zzStackSignalOnly('eemm,mmmm', inMC.replace('mc_','mc_{}_'.format(syst)),
                                           'smp', puWeightFile, lumi, amcatnlo=amcatnlo,
                                           asGroup=True)

    for varName in varNames:
        if not hasattr(systSaveFile, varName):
            d = systSaveFile.mkdir(varName)
            d.write()
        varDir = getattr(systSaveFile, varName)

        binning = _binning[varName]

        hErr = {}
        hUnfoldedByChan = {}
        #hUnfoldedAltByChan = {}
        for chan in channels:
            if not hasattr(varDir, chan):
                d = varDir.mkdir(chan)
                d.write()
            systSaveDir = getattr(varDir, chan)

            var = _variables[varName][chan]
            fVar = _varFunctions[varName][chan]
            sel = _selections[varName][chan]
            fSel = _selFunctions[varName]

            # regular weight, no systematics. Apply just in case.
            nominalWeight = baseMCWeight(chan, puWeightFile)
            reco[chan].applyWeight(nominalWeight, True)
            altReco[chan].applyWeight(nominalWeight, True)
            for s in recoSyst.values():
                try:
                    s[chan].applyWeight(nominalWeight, True)
                except KeyError:
                    pass

            response = getResponse(chan, true[chan], reco[chan], bkg[chan], var,
                                   fVar, binning, puWt, selectionStr=sel,
                                   selectionFunction=fSel)

            svd = TDecompSVD(response.Mresponse())
            sig = svd.GetSig()
            try:
                condition = sig.Max() / max(0., sig.Min())
            except ZeroDivisionError:
                condition = float('inf')
            print ''
            print '{}: {}'.format(chan, varName)
            print 'condition: {}'.format(condition)
            print ''

            hBkg = bkg[chan].makeHist(var, sel, binning, perUnitWidth=False)
            hData = data[chan].makeHist(var, sel, binning, perUnitWidth=False)
            hData -= hBkg

            unfolder = RooUnfoldIter(response, hData, nIter)

            #print chan
            #unfolders[chan].PrintTable(cout)

            hUnfolded = asrootpy(unfolder.Hreco())
            hUnfoldedByChan[chan] = hUnfolded.clone()

            hUnfolded /= hUnfolded.Integral(0,hUnfolded.GetNbinsX()+1)

            hFrame = hUnfolded.empty_clone()

            ### Unfold amc@NLO "data" as a sanity check
            hAlt = altReco[chan].makeHist(var, sel, binning,
                                          perUnitWidth=False)
            #unfolderAlt = RooUnfoldIter(response, hAlt, nIter)
            #hUnfoldedAlt = asrootpy(unfolderAlt.Hreco())
            #hUnfoldedAltByChan[chan] = hUnfoldedAlt.clone()
            #
            #hUnfoldedAlt /= hUnfoldedAlt.Integral(0,hUnfoldedAlt.GetNbinsX()+1)

            #### represent systematic errors as histograms where the bin content
            #### is the systematic error from that source
            hErr[chan] = {'up':{},'dn':{}}


            # PU reweight uncertainty
            for sys in ['up','dn']:
                hUnf = _compatibleHistOrNone(systSaveDir, 'pu_'+sys, hFrame)
                if hUnf is None or redo:
                    wtStr = baseMCWeight(chan, puWeightFile, puSyst=sys)
                    reco[chan].applyWeight(wtStr, True)

                    res = getResponse(chan, true[chan], reco[chan], bkg[chan],
                                      var, fVar, binning, fPUWt[sys],
                                      selectionStr=sel,
                                      selectionFunction=fSel)
                    unf = RooUnfoldIter(res, hData, nIter)

                    systSaveDir.cd()
                    hUnf = asrootpy(unf.Hreco())
                    hUnf /= hUnf.Integral(0,hUnf.GetNbinsX()+1)
                    hUnf.SetName('pu_'+sys)
                    hUnf.write()

                    reco[chan].applyWeight(nominalWeight, True)

                hErr[chan][sys]['pu'] = hUnf - hUnfolded
                hErr[chan][sys]['pu'].title = "PU"
                hErr[chan][sys]['pu'].fillstyle = 'solid'
                hErr[chan][sys]['pu'].drawstyle = "hist"
                hErr[chan][sys]['pu'].color = "green"
                hErr[chan][sys]['pu'].legendstyle = 'F'

            # lepton efficiency uncertainty
            for lep in ['e','m']:
                if lep not in chan:
                    continue

                for sys in ['up','dn']:
                    hUnf = _compatibleHistOrNone(systSaveDir, lep+'Eff_'+sys, hFrame)
                    if hUnf is None or redo:
                        wtArg = {lep+'Syst':sys}
                        wtStr = baseMCWeight(chan, puWeightFile, **wtArg)
                        reco[chan].applyWeight(wtStr, True)

                        res = getResponse(chan, true[chan], reco[chan], bkg[chan], var,
                                          fVar, binning, fPUWt[''], sys,
                                          selectionStr=sel,
                                          selectionFunction=fSel)
                        unf = RooUnfoldIter(res, hData, nIter)

                        systSaveDir.cd()
                        hUnf = asrootpy(unf.Hreco())
                        hUnf /= hUnf.Integral(0,hUnf.GetNbinsX()+1)
                        hUnf.SetName(lep+'Eff_'+sys)
                        hUnf.write()

                        reco[chan].applyWeight(nominalWeight, True)

                    hErr[chan][sys][lep+'Eff'] = hUnf - hUnfolded
                    hErr[chan][sys][lep+'Eff'].title = "{}on eff.".format('Electr' if lep == 'e' else 'Mu')
                    hErr[chan][sys][lep+'Eff'].fillstyle = 'solid'
                    hErr[chan][sys][lep+'Eff'].drawstyle = "hist"
                    hErr[chan][sys][lep+'Eff'].color = ("blue" if lep == 'e' else '#002db3')
                    hErr[chan][sys][lep+'Eff'].legendstyle = 'F'

            # unfolding uncertainty by checking difference with alternate generator
            hUnf = _compatibleHistOrNone(systSaveDir, 'generator', hFrame)
            if hUnf is None or redo:
                res = getResponse(chan, altTrue[chan], altReco[chan], bkg[chan], var,
                                  fVar, binning, fPUWt[''],
                                  selectionStr=sel,
                                  selectionFunction=fSel)
                unf = RooUnfoldIter(res, hData, nIter)

                systSaveDir.cd()
                hUnf = asrootpy(unf.Hreco())
                hUnf /= hUnf.Integral(0,hUnf.GetNbinsX()+1)
                hUnf.SetName('generator')
                hUnf.write()

            hErr[chan]['up']['generator'] = hUnf - hUnfolded
            hErr[chan]['dn']['generator'] = hErr[chan]['up']['generator'].clone()
            # make symmetrical
            for bUp, bDn in zip(hErr[chan]['up']['generator'],hErr[chan]['dn']['generator']):
                bUp.value = abs(bUp.value)
                bDn.value = -1 * abs(bDn.value)

            for sys in ['up','dn']:
                hErr[chan][sys]['generator'].title = "Generator choice"
                hErr[chan][sys]['generator'].fillstyle = 'solid'
                hErr[chan][sys]['generator'].drawstyle = "hist"
                hErr[chan][sys]['generator'].color = "magenta"
                hErr[chan][sys]['generator'].legendstyle = 'F'


            # luminosity uncertainty
            lumiUncBase = .062
            for sys in ['up','dn']:
                lumiUnc = lumiUncBase
                if sys == 'dn':
                    lumiUnc *= -1.

                hUnf = _compatibleHistOrNone(systSaveDir, 'lumi_'+sys, hFrame)
                if hUnf is None or redo:
                    reco[chan].applyWeight(str(1.+lumiUnc))
                    true[chan].applyWeight(str(1.+lumiUnc), True)
                    res = getResponse(chan, true[chan], reco[chan], bkg[chan], var,
                                      fVar, binning, fPUWt[''],
                                      selectionStr=sel,
                                      selectionFunction=fSel)
                    unf = RooUnfoldIter(res, hData, nIter)

                    systSaveDir.cd()
                    hUnf = asrootpy(unf.Hreco())
                    hUnf /= hUnf.Integral(0,hUnf.GetNbinsX()+1)
                    hUnf.SetName('lumi_'+sys)
                    hUnf.write()

                    reco[chan].applyWeight(nominalWeight, True)
                    true[chan].applyWeight('', True)

                hErr[chan][sys]['lumi'] = hUnf - hUnfolded
                hErr[chan][sys]['lumi'].title = "Luminosity"
                hErr[chan][sys]['lumi'].fillstyle = 'solid'
                hErr[chan][sys]['lumi'].drawstyle = "hist"
                hErr[chan][sys]['lumi'].color = "orange"
                hErr[chan][sys]['lumi'].legendstyle = 'F'


            # Fake rate uncertainty
            for sys in ['up','dn']:
                for lep in ['e','m']:
                    if lep not in chan:
                        continue

                    hUnf = _compatibleHistOrNone(systSaveDir, lep+'FR_'+sys, hFrame)
                    if hUnf is None or redo:
                        res = getResponse(chan, true[chan], reco[chan],
                                          bkgSyst[lep+sys][chan], var,
                                          fVar, binning, fPUWt[''],
                                          selectionStr=sel,
                                          selectionFunction=fSel)
                        unf = RooUnfoldIter(res, hData, nIter)

                        systSaveDir.cd()
                        hUnf = asrootpy(unf.Hreco())
                        hUnf /= hUnf.Integral(0,hUnf.GetNbinsX()+1)
                        hUnf.SetName(lep+'FR_'+sys)
                        hUnf.write()

                    hErr[chan][sys][lep+'FR'] = hUnf - hUnfolded
                    if lep == 'e':
                        hErr[chan][sys][lep+'FR'].title = "Electron Fake Rate"
                        hErr[chan][sys][lep+'FR'].color = "#00cc99"
                    elif lep == 'm':
                        hErr[chan][sys][lep+'FR'].title = "Muon Fake Rate"
                        hErr[chan][sys][lep+'FR'].color = "#00ff00"
                    hErr[chan][sys][lep+'FR'].fillstyle = 'solid'
                    hErr[chan][sys][lep+'FR'].drawstyle = "hist"
                    hErr[chan][sys][lep+'FR'].legendstyle = 'F'


            # Jet energy scale and resolution uncertainties (jet variables only)
            if 'jet' in varName.lower() or 'jj' in varName.lower():
                for shift in ['up','dn']:
                    sysStr = 'Up' if shift == 'up' else 'Down'

                    for sys in ['jer', 'jes']:
                        shiftedName = varName + '_' + sys + sysStr
                        varShifted = _variables[shiftedName][chan]
                        fVarShifted = _varFunctions[shiftedName][chan]
                        selShifted = _selections[shiftedName][chan]
                        fSelShifted = _selFunctions[shiftedName]

                        hUnf = _compatibleHistOrNone(systSaveDir, sys+'_'+shift, hFrame)
                        if hUnf is None or redo:
                            res = getResponse(chan, true[chan], reco[chan], bkg[chan],
                                              varShifted, fVarShifted, binning, fPUWt[''],
                                              altVar=var,
                                              selectionStr=selShifted,
                                              selectionFunction=fSelShifted,
                                              selectionStrAlt=sel,
                                              varFunctionAlt=fVar,
                                              selectionFunctionAlt=fSel)
                            unf = RooUnfoldIter(res, hData, nIter)

                            systSaveDir.cd()
                            hUnf = asrootpy(unf.Hreco())
                            hUnf /= hUnf.Integral(0,hUnf.GetNbinsX()+1)
                            hUnf.SetName(sys+'_'+shift)
                            hUnf.write()

                        hErr[chan][shift][sys] = hUnf - hUnfolded
                        hErr[chan][shift][sys].fillstyle = 'solid'
                        hErr[chan][shift][sys].drawstyle = "hist"
                        hErr[chan][shift][sys].legendstyle = 'F'


                    hErr[chan][shift]['jer'].title = "Jet energy resolution"
                    hErr[chan][shift]['jer'].color = "cyan"
                    hErr[chan][shift]['jes'].title = "Jet energy scale"
                    hErr[chan][shift]['jes'].color = "darkblue"


            if 'e' in chan:
                for sys in ['up', 'dn']:
                    # EES
                    hUnf = _compatibleHistOrNone(systSaveDir, 'ees_'+sys, hFrame)
                    if hUnf is None or redo:
                        res = getResponse(chan, true[chan],
                                          recoSyst['eScale'+sys[0].upper()+sys[1:]][chan],
                                          bkg[chan], var, fVar, binning,
                                          fPUWt[''],
                                          selectionStr=sel,
                                          selectionFunction=fSel)
                        unf = RooUnfoldIter(res, hData, nIter)

                        systSaveDir.cd()
                        hUnf = asrootpy(unf.Hreco())
                        hUnf /= hUnf.Integral(0,hUnf.GetNbinsX()+1)
                        hUnf.SetName('ees_'+sys)
                        hUnf.write()

                    hErr[chan][sys]['ees'] = hUnf - hUnfolded
                    hErr[chan][sys]['ees'].title = "Electron energy scale"
                    hErr[chan][sys]['ees'].fillstyle = 'solid'
                    hErr[chan][sys]['ees'].drawstyle = "hist"
                    hErr[chan][sys]['ees'].color = "purple"
                    hErr[chan][sys]['ees'].legendstyle = 'F'

                    # EER (rho)
                    hUnf = _compatibleHistOrNone(systSaveDir, 'eerRho_'+sys, hFrame)
                    if hUnf is None or redo:
                        res = getResponse(chan, true[chan],
                                          recoSyst['eRhoRes'+sys[0].upper()+sys[1:]][chan],
                                          bkg[chan], var, fVar, binning,
                                          fPUWt[''],
                                          selectionStr=sel,
                                          selectionFunction=fSel)
                        unf = RooUnfoldIter(res, hData, nIter)

                        systSaveDir.cd()
                        hUnf = asrootpy(unf.Hreco())
                        hUnf /= hUnf.Integral(0,hUnf.GetNbinsX()+1)
                        hUnf.SetName('eerRho_'+sys)
                        hUnf.write()

                    hErr[chan][sys]['eerRho'] = hUnf - hUnfolded
                    hErr[chan][sys]['eerRho'].title = "Electron energy res. (rho)"
                    hErr[chan][sys]['eerRho'].fillstyle = 'solid'
                    hErr[chan][sys]['eerRho'].drawstyle = "hist"
                    hErr[chan][sys]['eerRho'].color = "lavender"
                    hErr[chan][sys]['eerRho'].legendstyle = 'F'

                # EER (phi)
                hUnf = _compatibleHistOrNone(systSaveDir, 'eerPhi_up', hFrame)
                if hUnf is None or redo:
                    res = getResponse(chan, true[chan],
                                      recoSyst['ePhiResUp'][chan],
                                      bkg[chan], var, fVar, binning,
                                      fPUWt[''],
                                      selectionStr=sel,
                                      selectionFunction=fSel)
                    unf = RooUnfoldIter(res, hData, nIter)

                    systSaveDir.cd()
                    hUnf = asrootpy(unf.Hreco())
                    hUnf /= hUnf.Integral(0,hUnf.GetNbinsX()+1)
                    hUnf.SetName('eerPhi_up')
                    hUnf.write()

                hErr[chan]['up']['eerPhi'] = hUnf - hUnfolded
                hErr[chan]['up']['eerPhi'].title = "Electron energy res. (phi)"
                hErr[chan]['up']['eerPhi'].fillstyle = 'solid'
                hErr[chan]['up']['eerPhi'].drawstyle = "hist"
                hErr[chan]['up']['eerPhi'].color = "violet"
                hErr[chan]['up']['eerPhi'].legendstyle = 'F'

                hErr[chan]['dn']['eerPhi'] = hErr[chan]['up']['eerPhi'].clone()
                hErr[chan]['dn']['eerPhi'].title = "Electron energy res. (phi)"
                hErr[chan]['dn']['eerPhi'].fillstyle = 'solid'
                hErr[chan]['dn']['eerPhi'].drawstyle = "hist"
                hErr[chan]['dn']['eerPhi'].color = "violet"
                hErr[chan]['dn']['eerPhi'].legendstyle = 'F'


            if 'm' in chan:
                for sys in ['up', 'dn']:
                    # MES/MER (closure)
                    hUnf = _compatibleHistOrNone(systSaveDir, 'mClosure_'+sys, hFrame)
                    if hUnf is None or redo:
                        res = getResponse(chan, true[chan],
                                          recoSyst['mClosure'+sys[0].upper()+sys[1:]][chan],
                                          bkg[chan], var, fVar, binning,
                                          fPUWt[''],
                                          selectionStr=sel,
                                          selectionFunction=fSel)
                        unf = RooUnfoldIter(res, hData, nIter)

                        systSaveDir.cd()
                        hUnf = asrootpy(unf.Hreco())
                        hUnf /= hUnf.Integral(0,hUnf.GetNbinsX()+1)
                        hUnf.SetName('mClosure_'+sys)
                        hUnf.write()

                    hErr[chan][sys]['mClosure'] = hUnf - hUnfolded
                    hErr[chan][sys]['mClosure'].title = "Muon calibration"
                    hErr[chan][sys]['mClosure'].fillstyle = 'solid'
                    hErr[chan][sys]['mClosure'].drawstyle = "hist"
                    hErr[chan][sys]['mClosure'].color = "#c61aff"
                    hErr[chan][sys]['mClosure'].legendstyle = 'F'


            # PDF uncertainties
            hUnfUp = _compatibleHistOrNone(systSaveDir, 'pdf_up', hFrame)
            hUnfDn = _compatibleHistOrNone(systSaveDir, 'pdf_dn', hFrame)
            if hUnfUp is None or hUnfDn is None or redo:
                resDn, resUp = getResponsePDFErrors(chan, true[chan],
                                                    reco[chan], bkg[chan], var,
                                                    fVar, binning, fPUWt[''],
                                                    selectionStr = sel,
                                                    selectionFunction=fSel)
                unfUp = RooUnfoldIter(resUp, hData, nIter)
                unfDn = RooUnfoldIter(resDn, hData, nIter)

                systSaveDir.cd()
                hUnfUp = asrootpy(unfUp.Hreco())
                hUnfUp /= hUnfUp.Integral(0,hUnfUp.GetNbinsX()+1)
                hUnfUp.SetName('pdf_up')
                hUnfUp.write()
                hUnfDn = asrootpy(unfDn.Hreco())
                hUnfDn /= hUnfDn.Integral(0,hUnfDn.GetNbinsX()+1)
                hUnfDn.SetName('pdf_dn')
                hUnfDn.write()

            hErr[chan]['up']['pdf'] = hUnfUp - hUnfolded
            hErr[chan]['up']['pdf'].title = "PDF"
            hErr[chan]['up']['pdf'].fillstyle = 'solid'
            hErr[chan]['up']['pdf'].drawstyle = "hist"
            hErr[chan]['up']['pdf'].color = "#80aaff"
            hErr[chan]['up']['pdf'].legendstyle = 'F'
            hErr[chan]['dn']['pdf'] = hUnfDn - hUnfolded
            hErr[chan]['dn']['pdf'].title = "PDF"
            hErr[chan]['dn']['pdf'].fillstyle = 'solid'
            hErr[chan]['dn']['pdf'].drawstyle = "hist"
            hErr[chan]['dn']['pdf'].color = "#80aaff"
            hErr[chan]['dn']['pdf'].legendstyle = 'F'


            # alpha_s uncertainty
            hUnf1 = _compatibleHistOrNone(systSaveDir, 'alphas_1', hFrame)
            hUnf2 = _compatibleHistOrNone(systSaveDir, 'alphas_2', hFrame)
            if hUnf1 is None or hUnf2 is None or redo:
                res1, res2 = getResponseAlphaSErrors(chan, true[chan],
                                                     reco[chan], bkg[chan], var,
                                                     fVar, binning, fPUWt[''],
                                                     selectionStr = sel,
                                                     selectionFunction=fSel)
                unf1 = RooUnfoldIter(res1, hData, nIter)
                unf2 = RooUnfoldIter(res2, hData, nIter)

                systSaveDir.cd()
                hUnf1 = asrootpy(unf1.Hreco())
                hUnf1 /= hUnf1.Integral(0,hUnf1.GetNbinsX()+1)
                hUnf1.SetName('alphas_1')
                hUnf1.write()
                hUnf2 = asrootpy(unf2.Hreco())
                hUnf2 /= hUnf2.Integral(0,hUnf2.GetNbinsX()+1)
                hUnf2.SetName('alphas_2')
                hUnf2.write()

            # PDF4LHC recommendation
            hAlphaSErr = (hUnf2 - hUnf1) / 2.
            # Not sure which is up and which is down, so just do the
            # absolute value up and down.
            # Factor of 1.5 comes from slide 14 of
            # https://indico.cern.ch/event/459797/contributions/1961581/attachments/1181555/1800214/mcaod-Feb15-2016.pdf
            hErr[chan]['up']['alphas'] = hAlphaSErr.empty_clone()
            hErr[chan]['dn']['alphas'] = hAlphaSErr.empty_clone()
            for b, bUp, bDn in zip(hAlphaSErr, hErr[chan]['up']['alphas'], hErr[chan]['dn']['alphas']):
                absNum = 1.5 * abs(b.value)
                bUp.value = absNum
                bDn.value = -1 * absNum

            for sys in ['up','dn']:
                hErr[chan][sys]['alphas'].title = "#alpha_{s}"
                hErr[chan][sys]['alphas'].fillstyle = 'solid'
                hErr[chan][sys]['alphas'].drawstyle = "hist"
                hErr[chan][sys]['alphas'].color = "#4e72ba"
                hErr[chan][sys]['alphas'].legendstyle = 'F'


            # QCD scale uncertainties
            hUnfUp = _compatibleHistOrNone(systSaveDir, 'scale_up', hFrame)
            hUnfDn = _compatibleHistOrNone(systSaveDir, 'scale_dn', hFrame)
            if hUnfUp is None or hUnfDn is None or redo:
                res = getResponseScaleErrors(chan, true[chan], reco[chan],
                                             bkg[chan], var, fVar, binning,
                                             fPUWt[''], selectionStr = sel,
                                             selectionFunction=fSel)

                unf = [RooUnfoldIter(r, hData, nIter) for r in res]

                hUnf = [asrootpy(u.Hreco()) for u in unf]

                hUnfUp = hUnfolded.empty_clone()
                hUnfDn = hUnfolded.empty_clone()
                for bUp, bDn, allUnfoldedBins in zip(hUnfUp, hUnfDn, zip(*hUnf)):
                    bUp.value = max(b.value for b in allUnfoldedBins)
                    bDn.value = min(b.value for b in allUnfoldedBins)

                systSaveDir.cd()
                hUnfUp.SetName('scale_up')
                hUnfUp /= hUnfUp.Integral(0,hUnfUp.GetNbinsX()+1)
                hUnfUp.write()
                hUnfDn.SetName('scale_dn')
                hUnfDn /= hUnfDn.Integral(0,hUnfDn.GetNbinsX()+1)
                hUnfDn.write()

            hErr[chan]['up']['scale'] = hUnfUp - hUnfolded
            hErr[chan]['up']['scale'].title = "QCD Scale"
            hErr[chan]['up']['scale'].fillstyle = 'solid'
            hErr[chan]['up']['scale'].drawstyle = "hist"
            hErr[chan]['up']['scale'].color = "#800000"
            hErr[chan]['up']['scale'].legendstyle = 'F'
            hErr[chan]['dn']['scale'] = hUnfDn - hUnfolded
            hErr[chan]['dn']['scale'].title = "QCD Scale"
            hErr[chan]['dn']['scale'].fillstyle = 'solid'
            hErr[chan]['dn']['scale'].drawstyle = "hist"
            hErr[chan]['dn']['scale'].color = "#800000"
            hErr[chan]['dn']['scale'].legendstyle = 'F'


            ### Since we have no LHE info for MCFM samples, we just vary the
            ### normalization uncertainty up and down by the cross section's
            ### PDF and scale uncertainties
            # gg LO uncertainty: +18%/-15% (MCFM, AN-2016-029)
            mcfmUnc = {'up':.18,'dn':-.15}
            for sys, shift in mcfmUnc.iteritems():
                hUnf = _compatibleHistOrNone(systSaveDir, 'mcfmxsec_'+sys, hFrame)
                if hUnf is None or redo:
                    reco[chan].applyWeight(nominalWeight, True)
                    for n,s in reco[chan].itersamples():
                        if 'GluGluZZ' in n:
                            s.applyWeight(str(1.+shift))
                    for n,s in true[chan].itersamples():
                        if 'GluGluZZ' in n:
                            s.applyWeight(str(1.+shift))
                    res = getResponse(chan, true[chan], reco[chan], bkg[chan], var,
                                      fVar, binning, fPUWt[''],
                                      selectionStr=sel,
                                      selectionFunction=fSel)
                    unf = RooUnfoldIter(res, hData, nIter)

                    systSaveDir.cd()
                    hUnf = asrootpy(unf.Hreco())
                    hUnf /= hUnf.Integral(0,hUnf.GetNbinsX()+1)
                    hUnf.SetName('mcfmxsec_'+sys)
                    hUnf.write()

                    reco[chan].applyWeight(nominalWeight, True)
                    true[chan].applyWeight('',True)

                hErr[chan][sys]['mcfmxsec'] = hUnf - hUnfolded
                hErr[chan][sys]['mcfmxsec'].title = "MCFM PDF/scale"
                hErr[chan][sys]['mcfmxsec'].fillstyle = 'solid'
                hErr[chan][sys]['mcfmxsec'].drawstyle = "hist"
                hErr[chan][sys]['mcfmxsec'].color = "red"
                hErr[chan][sys]['mcfmxsec'].legendstyle = 'F'


            ### Get the total uncertainties
            # this method of assigning up vs down should be revisited
            hUncUp = hUnfolded.empty_clone()
            hUncDn = hUnfolded.empty_clone()
            uncTypes = hErr[chan]['up'].keys()
            for bUncUp, bUncDn, allUncUp, allUncDn in zip(hUncUp, hUncDn,
                                                          zip(*hErr[chan]['up'].values()),
                                                          zip(*hErr[chan]['dn'].values())):
                for b1,b2 in zip(allUncUp,allUncDn):
                    if b1.value > 0. and b2.value < 0.:
                        bUncUp.value += b1.value**2
                        bUncDn.value += b2.value**2
                    elif b2.value > 0. and b1.value < 0.:
                        bUncUp.value += b2.value**2
                        bUncDn.value += b1.value**2
                    else:
                        bUncUp.value += max(b1.value,b2.value,0.)**2
                        bUncDn.value += min(b1.value,b2.value,0.)**2

                bUncUp.value = sqrt(bUncUp.value)
                bUncDn.value = sqrt(bUncDn.value)


            # Make all error histograms positive
            # and make uncertainties fractional (as a percentage)
            for sys in hErr[chan].values():
                for h in sys.values():
                    h /= hUnfolded
                    h *= 100.
                    for b in h:
                        b.value = abs(b.value)

            # Make plots of uncertainties (added linearly)
            cErrUp = Canvas(1000,1000)
            errStackUp = HistStack(hErr[chan]['up'].values(), drawstyle = 'histnoclear')
            draw(errStackUp, cErrUp, xtitle=_xTitle[varName], ytitle="+Error (%)",yerror_in_padding=False)
            leg = makeLegend(cErrUp, *hErr[chan]['up'].values(), leftmargin=0.25,
                             entryheight=.02, entrysep=.007, textsize=.022,
                             rightmargin=.25)
            leg.Draw('same')
            style.setCMSStyle(cErrUp, '', dataType='Preliminary', intLumi=lumi)
            cErrUp.Print(_join(plotDir, 'errUp_{}_{}.png'.format(varName, chan)))
            cErrUp.Print(_join(plotDir, 'errUp_{}_{}.C'.format(varName, chan)))

            cErrDn = Canvas(1000,1000)
            errStackDn = HistStack(hErr[chan]['dn'].values(), drawstyle = 'histnoclear')
            draw(errStackDn, cErrDn, xtitle=_xTitle[varName], ytitle="-Error (%)",yerror_in_padding=False)
            leg = makeLegend(cErrDn, *hErr[chan]['dn'].values(), leftmargin=0.25,
                             entryheight=.02, entrysep=.007, textsize=.022,
                             rightmargin=.25)
            leg.Draw('same')
            style.setCMSStyle(cErrDn, '', dataType='Preliminary', intLumi=lumi)
            cErrDn.Print(_join(plotDir, 'errDown_{}_{}.png'.format(varName, chan)))
            cErrDn.Print(_join(plotDir, 'errDown_{}_{}.C'.format(varName, chan)))

            # Errors should no longer be fractional or a percentage
            for sys in hErr[chan].values():
                for h in sys.values():
                    h *= hUnfolded
                    h /= 100.

            ### plot
            hUnfolded.color = 'black'
            hUnfolded.drawstyle = 'PE'
            hUnfolded.legendstyle = 'LPE'
            hUnfolded.title = 'Data + stat. unc.'
            _normalizeBins(hUnfolded)

            #hUnfoldedAlt.color = 'magenta'
            #hUnfoldedAlt.drawstyle = 'hist'
            #hUnfoldedAlt.fillstyle = 'hollow'
            #hUnfoldedAlt.legendstyle = 'L'
            #hUnfoldedAlt.title = 'Unfolded {}'.format(signalNameAlt)
            #_normalizeBins(hUnfoldedAlt)

            if not sel:
                selGen = _wrongZRejectionStr[chan]
            elif isinstance(sel, str):
                selGen = sel + ' && ' + _wrongZRejectionStr[chan]
            else:
                # better be an iterable of strings
                selGen = [(s + ' && ' if s else '') + _wrongZRejectionStr[chan] for s in sel]

            hTrue = true[chan].makeHist(var, selGen, binning, perUnitWidth=False)
            hTrue.fillcolor = '#99ccff'
            hTrue.linecolor = '#000099'
            hTrue.drawstyle = 'hist'
            hTrue.fillstyle = 'solid'
            hTrue.legendstyle = 'F'
            hTrue.title = '{}'.format(signalName)
            hTrue /= hTrue.Integral(0,hTrue.GetNbinsX()+1)
            _normalizeBins(hTrue)

            hTrueAlt = altTrue[chan].makeHist(var, selGen, binning,
                                              perUnitWidth=False)
            hTrueAlt.color = 'red'
            hTrueAlt.drawstyle = 'hist'
            hTrueAlt.fillstyle = 'hollow'
            hTrueAlt.legendstyle = 'L'
            hTrueAlt.SetLineWidth(hTrueAlt.GetLineWidth()*2)
            hTrueAlt.title = '{}'.format(signalNameAlt)
            hTrueAlt /= hTrueAlt.Integral(0,hTrueAlt.GetNbinsX()+1)
            _normalizeBins(hTrueAlt)

            _normalizeBins(hUncUp)
            _normalizeBins(hUncDn)

            if varName in _blind:
                for b, bUp, bDn in zip(hUnfolded, hUncUp, hUncDn):
                    if hUnfolded.xaxis.GetBinLowEdge(b.idx) >= _blind[varName]:
                        b.value = 0
                        b.error = 0
                        bUp.value = 0
                        bUp.error = 0
                        bDn.value = 0
                        bDn.error = 0

            errorBand = makeErrorBand(hUnfolded, hUncUp, hUncDn)

            cUnf = Canvas(1000,1200)
            mainPad, ratioPad1, ratioPad2 = addPadsBelow(cUnf, 0.15, 0.15, bottomMargin=0.35)

            mainPad.cd()
            (xaxis, yaxis), (xmin,xmax,ymin,ymax) = draw([hTrue, hTrueAlt,
                                                          #hUnfoldedAlt,
                                                          hUnfolded,
                                                          errorBand], mainPad,
                                                         xtitle=_xTitle[varName],
                                                         ytitle=_yTitle[varName])#,
                                                         #yerror_in_padding=False)
            yaxis.SetTitleSize(0.75*yaxis.GetTitleSize())
            yaxis.SetTitleOffset(1.25*yaxis.GetTitleOffset())
            yaxis.SetLabelSize(0.82*yaxis.GetLabelSize())
            xaxis.SetLabelSize(0.82*xaxis.GetLabelSize())

            leg = makeLegend(cUnf, hTrueAlt, #hUnfoldedAlt,
                             hTrue, errorBand,
                             hUnfolded, **_legParams[varName])

            if varName in _blind and _blind[varName] < xmax:
                box = TBox(max(xmin,_blind[varName]), ymin, xmax, ymax)
                box.SetFillColor(1)
                box.SetFillStyle(3002)
                box.Draw("same")
                leg.SetFillStyle(1001)

            leg.Draw("same")

            latex = TLatex()
            latex.SetNDC()
            latex.SetTextSize(.13)
            latex.SetTextFont(62)
            latex.SetTextAlign(11)

            ratioPad1.cd()
            ratio1, unity1 = makeRatio(hUnfolded, hTrue)
            ratioError1 = makeErrorBand(hUnfolded/hTrue, hUncUp/hTrue,
                                        hUncDn/hTrue)
            (ratio1X, ratio1Y), ratio1Limits = draw([ratio1,ratioError1],
                                                    ratioPad1,
                                                    ytitle='Data / MC',
                                                    xlimits=(xmin,xmax),
                                                    ylimits=(0.50001,1.9999),
                                                    ydivisions=5)
            unity1.Draw("same")
            latex.DrawLatex(0.15, 0.8, signalName)


            ratioPad2.cd()
            ratio2, unity2 = makeRatio(hUnfolded, hTrueAlt)
            ratioError2 = makeErrorBand(hUnfolded/hTrueAlt, hUncUp/hTrueAlt,
                                        hUncDn/hTrueAlt)
            (ratio2X, ratio2Y), ratio2Limits = draw([ratio2,ratioError2],
                                                    ratioPad2,
                                                    ytitle='Data / MC',
                                                    xlimits=(xmin,xmax),
                                                    ylimits=(0.5,1.9999),
                                                    ydivisions=5)
            unity2.Draw("same")
            latex.SetTextSize(latex.GetTextSize() * ratioPad1.height / ratioPad2.height)
            latex.DrawLatex(0.15, 1.-.2*ratioPad1.height/ratioPad2.height,
                            signalNameAlt)

            cUnf.cd()
            ratioPad2.Draw()
            ratioPad1.Draw()
            mainPad.Draw()

            fixRatioAxes(xaxis,yaxis,ratio1X,ratio1Y, mainPad.height, ratioPad1.height)
            fixRatioAxes(ratio1X,ratio1Y,ratio2X,ratio2Y, ratioPad1.height, ratioPad2.height)

            style.setCMSStyle(cUnf, '', dataType='Preliminary', intLumi=lumi)
            cUnf.Print(_join(plotDir, "unfold_{}_{}.png".format(varName, chan)))
            cUnf.Print(_join(plotDir, "unfold_{}_{}.C".format(varName, chan)))

            cRes = Canvas(1000,1000)
            hRes = asrootpy(response.Hresponse())
            hRes.drawstyle = 'colztext'
            hRes.xaxis.title = '\\text{Reco} '+_xTitle[varName]
            hRes.yaxis.title = '\\text{True} '+_xTitle[varName]
            hRes.draw()
            style.setCMSStyle(cRes, '', dataType='Preliminary Simulation', intLumi=lumi)
            cRes.Print(_join(plotDir, "response_{}_{}.png".format(varName, chan)))
            cRes.Print(_join(plotDir, "response_{}_{}.C".format(varName, chan)))

            cCov = Canvas(1000,1000)
            covariance = unfolder.Ereco(2) # TMatrixD
            covariance.Draw("colztext")
            style.setCMSStyle(cCov, '', dataType='Preliminary', intLumi=lumi)
            cCov.Print(_join(plotDir, "covariance_{}_{}.png".format(varName, chan)))
            cCov.Print(_join(plotDir, "covariance_{}_{}.C".format(varName, chan)))

        hTot = sum(hUnfoldedByChan.values())
        hTot.color = 'black'
        hTot.drawstyle = 'PE'
        hTot.legendstyle = 'LPE'
        hTot.title = 'Data + stat. unc.'
        # total normalization
        hTotNoNorm = hTot.clone()
        hTot /= hTot.Integral(0,hTot.GetNbinsX()+1)
        _normalizeBins(hTot)

        #hTotAlt = sum(hUnfoldedAltByChan.values())
        #hTotAlt.color = 'magenta'
        #hTotAlt.drawstyle = 'hist'
        #hTotAlt.fillstyle = 'hollow'
        #hTotAlt.legendstyle = 'L'
        #hTotAlt.title = 'Unfolded {}'.format(signalNameAlt)
        #hTotAlt /= hTotAlt.Integral(0,hTotAlt.GetNbinsX()+1)
        #_normalizeBins(hTotAlt)

        selGen = {}
        for chan in channels:
            if not _selections[varName][chan]:
                selGen[chan] = _wrongZRejectionStr[chan]
            elif isinstance(_selections[varName][chan], str):
                selGen[chan] = _selections[varName][chan] + ' && ' + _wrongZRejectionStr[chan]
            else:
                # better be an iterable of strings
                selGen[chan] = [(s + ' && ' if s else '') + _wrongZRejectionStr[chan] for s in _selections[varName][chan]]

        hTrue = true.makeHist(_variables[varName], selGen,
                              binning, perUnitWidth=False)
        hTrue.fillcolor = '#99ccff'
        hTrue.linecolor = '#000099'
        hTrue.drawstyle = 'hist'
        hTrue.fillstyle = 'solid'
        hTrue.legendstyle = 'F'
        hTrue.title = '{}'.format(signalName)
        hTrue /= hTrue.Integral(0,hTrue.GetNbinsX()+1)
        _normalizeBins(hTrue)

        hTrueAlt = altTrue.makeHist(_variables[varName], selGen,
                                    binning, perUnitWidth=False)
        hTrueAlt.color = 'r'
        hTrueAlt.drawstyle = 'hist'
        hTrueAlt.fillstyle = 'hollow'
        hTrueAlt.legendstyle = 'L'
        hTrueAlt.SetLineWidth(hTrueAlt.GetLineWidth()*2)
        hTrueAlt.title = '{}'.format(signalNameAlt)
        hTrueAlt /= hTrueAlt.Integral(0,hTrueAlt.GetNbinsX()+1)
        _normalizeBins(hTrueAlt)

        hUncTot = {}
        uncList = []
        for chan in channels:
            for sys in ['up','dn']:
                uncList += hErr[chan][sys].keys()
        uncList = set(uncList)
        for sys in ['up','dn']:
            hUncTot[sys] = {}
            for unc in uncList:
                hUncTot[sys][unc] = hTot.empty_clone()
                for chan in _channels:
                    try:
                        hThis = hErr[chan][sys][unc].clone()
                    except KeyError:
                        continue

                    # weight by size of channel, re-normalize to total
                    hThis *= hUnfoldedByChan[chan] / hTotNoNorm
                    # for bErr, bChan, bTot in zip(hThis, hUnfoldedByChan[chan], hTotNoNorm):
                    #     try:
                    #         bErr.value *= bChan.value / bTot.value
                    #     except ZeroDivisionError:
                    #         pass
                    hUncTot[sys][unc] += hThis


        hUncUp = hTot.clone()
        for bTot, allbs in zip(hUncUp, zip(*hUncTot['up'].values())):
            bTot.value = sqrt(sum(b.value**2 for b in allbs))
        hUncDn = hTot.clone()
        for bTot, allbs in zip(hUncDn, zip(*hUncTot['dn'].values())):
            bTot.value = sqrt(sum(b.value**2 for b in allbs))

        _normalizeBins(hUncUp)
        _normalizeBins(hUncDn)

        if varName in _blind:
            for b, bUp, bDn in zip(hTot, hUncUp, hUncDn):
                if hUnfolded.xaxis.GetBinLowEdge(b.idx) >= _blind[varName]:
                    b.value = 0
                    b.error = 0
                    bUp.value = 0
                    bUp.error = 0
                    bDn.value = 0
                    bDn.error = 0

        errorBand = makeErrorBand(hTot, hUncUp, hUncDn)

        cUnf = Canvas(1000,1200)
        mainPad, ratioPad1, ratioPad2 = addPadsBelow(cUnf, 0.15, 0.15, bottomMargin=.35)

        mainPad.cd()

        (xaxis, yaxis), (xmin,xmax,ymin,ymax) = draw([hTrue, hTrueAlt,
                                                      #hTotAlt,
                                                      hTot, errorBand],
                                                     mainPad, xtitle=_xTitle[varName],
                                                     ytitle=_yTitle[varName])#,
                                                     #yerror_in_padding=False)
        yaxis.SetTitleSize(0.75*yaxis.GetTitleSize())
        yaxis.SetTitleOffset(1.25*yaxis.GetTitleOffset())
        yaxis.SetLabelSize(0.82*yaxis.GetLabelSize())
        xaxis.SetLabelSize(0.82*xaxis.GetLabelSize())

        leg = makeLegend(mainPad, hTrueAlt, #hTotAlt,
                         hTrue, errorBand,
                         hTot, **_legParams[varName])

        if varName in _blind and _blind[varName] < xmax:
            box = TBox(max(xmin,_blind[varName]), ymin, xmax, ymax)
            box.SetFillColor(1)
            box.SetFillStyle(3002)
            box.Draw("same")
            leg.SetFillStyle(1001)

        leg.Draw("same")

        latex = TLatex()
        latex.SetNDC()
        latex.SetTextSize(.13)
        latex.SetTextFont(62)
        latex.SetTextAlign(11)

        ratioPad1.cd()
        ratio1, unity1 = makeRatio(hTot, hTrue)
        ratioError1 = makeErrorBand(hTot/hTrue, hUncUp/hTrue, hUncDn/hTrue)
        (ratio1X, ratio1Y), ratio1Limits = draw([ratio1,ratioError1], ratioPad1, ytitle='Data / MC',
                                                xlimits=(xmin,xmax),
                                                ylimits=(0.5,1.99999), ydivisions=5)
        unity1.Draw("same")
        latex.DrawLatex(0.15, 0.8, signalName)

        ratioPad2.cd()
        ratio2, unity2 = makeRatio(hTot, hTrueAlt)
        ratioError2 = makeErrorBand(hTot/hTrueAlt, hUncUp/hTrueAlt, hUncDn/hTrueAlt)
        (ratio2X, ratio2Y), ratio2Limits = draw([ratio2,ratioError2], ratioPad2, ytitle='Data / MC',
                                                xlimits=(xmin,xmax),
                                                ylimits=(0.5,1.99999), ydivisions=5)
        unity2.Draw("same")
        latex.SetTextSize(latex.GetTextSize() * ratioPad1.height / ratioPad2.height)
        latex.DrawLatex(0.15, 1.-.2*ratioPad1.height/ratioPad2.height, signalNameAlt)

        cUnf.cd()
        ratioPad2.Draw()
        ratioPad1.Draw()
        mainPad.Draw()

        fixRatioAxes(xaxis,yaxis,ratio1X,ratio1Y, mainPad.height, ratioPad1.height)
        fixRatioAxes(ratio1X,ratio1Y,ratio2X,ratio2Y, ratioPad1.height, ratioPad2.height)

        style.setCMSStyle(cUnf, '', dataType='Preliminary', intLumi=lumi)
        cUnf.Print(_join(plotDir, "unfold_{}.png".format(varName)))
        cUnf.Print(_join(plotDir, "unfold_{}.C".format(varName)))


if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser(description="Do full unfolding with all systematics")
    parser.add_argument('--dataDir', type=str, nargs='?',
                        default='uwvvNtuples_data_12oct2016',
                        help='Directory where data ntuples live')
    parser.add_argument('--mcDir', type=str, nargs='?',
                        default='uwvvNtuples_mc_12oct2016',
                        help='Directory where MC ntuples live')
    parser.add_argument('--plotDir', type=str, nargs='?',
                        default='/afs/cern.ch/user/n/nawoods/www/UWVVPlots/unfold',
                        help='Directory to put plots in')
    parser.add_argument('--fakeRateFile', type=str, nargs='?',
                        default='fakeRate_08sep2016',
                        help=('Name of fake rate file (assumed to be in usual '
                              'data directory unless full path is specified)'))
    parser.add_argument('--puWeightFile', type=str, nargs='?',
                        default='puWeight_69200_08sep2016.root',
                        help=('Name of pileup weight file (assumed to be in usual '
                              'data directory unless full path is specified)'))
    parser.add_argument('--lumi', type=float, nargs='?', default=15937.,
                        help='Integrated luminosity of sample (in pb^-1)')
    parser.add_argument('--nIter', type=int, nargs='?', default=8,
                        help='Number of iterations for D\'Agostini method')
    parser.add_argument('--redo', action='store_true',
                        help='Redo everything, ignoring stored results')
    parser.add_argument('--amcatnlo', action='store_true',
                        help='Use MadGraph5_aMC@NLO as the primary MC and '
                        'Powheg as the cross-check, instead of the other way '
                        'around.')
    parser.add_argument('--variables', type=str, nargs='*',
                        default=_varList,
                        help=('Names of variables to use. If not specified, '
                              'all are used ({})').format(', '.join(_varList)))

    args=parser.parse_args()

    systSaveFileName = _systSaveFileTemplate.format(args.nIter)
    if args.amcatnlo:
        systSaveFileName = systSaveFileName.replace('.root','_amcatnlo.root')

    if not _exists(args.plotDir):
        _mkdir(args.plotDir)
    elif not _isdir(args.plotDir):
        raise IOError("There is already some non-directory object called {}.".format(args.plotDir))

    with root_open(systSaveFileName, 'UPDATE') as systSaveFile:
        main(args.dataDir, args.mcDir, args.plotDir, args.fakeRateFile,
             args.puWeightFile, args.lumi, args.nIter, args.redo,
             systSaveFile, args.amcatnlo, *args.variables)
