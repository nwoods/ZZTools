import logging
from rootpy import log as rlog; rlog = rlog["/unfoldFast"]
# don't show most silly ROOT messages
logging.basicConfig(level=logging.WARNING)
rlog["/ROOT.TUnixSystem.SetDisplay"].setLevel(rlog.ERROR)
rlog["/ROOT.TROOT.Append"].setLevel(rlog.ERROR)

from rootpy.ROOT import gROOT, gStyle, TAttFill
gROOT.SetBatch(True)

from rootpy import asrootpy
from rootpy.io import root_open
from rootpy.plotting import Canvas, Legend, Hist, HistStack, Graph
from rootpy.plotting.utils import draw
from rootpy.ROOT import cout, TDecompSVD, TBox, TLatex, TLegend, TAttFill
from rootpy.ROOT import RooUnfoldResponse as Response
from rootpy.ROOT import RooUnfoldBayes as RooUnfoldIter # it's frequentist!
from rootpy.context import preserve_current_directory
import rootpy.compiled as _rootComp
from rootpy.stl import vector as _Vec
_VFloat = _Vec('float')

from PlotTools import PlotStyle as _Style, pdfViaTex as _pdfViaTex
from PlotTools import makeLegend, addPadsBelow, makeRatio, fixRatioAxes, makeErrorBand
from Utilities import WeightStringMaker, Z_MASS, deltaRString, deltaPhiString, zeroNegativeBins, combineWeights
from Analysis.setupStandardSamples import *
# from Analysis.unfoldingHelpers import getResponse, getResponsePDFErrors, \
#     getResponseScaleErrors, getResponseAlphaSErrors
from Analysis.weightHelpers import puWeight, baseMCWeight
from Metadata.metadata import sampleInfo

from os import environ as _env
from os import makedirs as _mkdir
from os import system as _bash
from os.path import join as _join
from os.path import isdir as _isdir
from os.path import exists as _exists
from math import sqrt

# need to load RooUnfold libraries for cling
try:
    _zztBaseDir = _env['zzt']
except KeyError:
    rlog.error("Can't find ZZTools base directory. Is your area set up properly?")
    raise

_style = _Style()
gStyle.SetLineScalePS(1.8)

_channels = ['eeee','eemm', 'mmmm']

# set up variables, selection, binnings etc.
# (jet-related variables and selections done later)
_variables = {
    'pt' : {c:'Pt' for c in _channels},
    'mass' : {c:'Mass' for c in _channels},
    'massFull' : {c:'Mass' for c in _channels},
    'eta' : {c:'abs(Eta)' for c in _channels},
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

_blind = {}

_binning = {
    'pt' : [25.*i for i in range(4)] + [100., 150., 200., 300.],
    'nJets' : [5,-0.5,4.5],
    'mass' : [100.] + [200.+50.*i for i in range(5)] + [500.,600.,800.],
    'massFull' : [80.,100.,120.,130.,150.,180.,200.,240.,300.,400.,1000],
    'eta' : [6,0.,6.],
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
    'l1Pt' : [0.,15.,30.,40.,50.]+[60.+15.*i for i in range(9)]+[195.,225.],#[14,0.,210.],#[15, 0., 150.],
    }

_units = {
    'pt' : 'GeV',
    'nJets' : '',
    'mass' : 'GeV',
    'massFull' : 'GeV',
    'eta' : '',
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
    'pt' : 'p_T^{\\text{ZZ}}',
    'nJets' : '# jets',
    'mass' : 'm_{\\text{ZZ}}',
    'massFull' : 'm_{4\\ell}',
    'eta' : '\\eta_{\\text{ZZ}}',
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
_yTitleNoNorm = {}
_yTitleTemp = '{prefix} \\frac{{d\\sigma_{{\\text{{fid}}}}}}{{d{xvar}}} {units}'
for var, prettyVar in _prettyVars.iteritems():
    xt = prettyVar
    if _units[var]:
        xt += ' \\, \\left(\\text{{{}}}\\right)'.format(_units[var])
        yt = _yTitleTemp.format(xvar=prettyVar,
                                prefix='\\frac{1}{\\sigma_{\\text{fid}}}',
                                units='\\, \\left( \\frac{{1}}{{\\text{{{unit}}}}} \\right)'.format(unit=_units[var]))
        ytnn = _yTitleTemp.format(xvar=prettyVar, prefix='',
                                  units='\\, \\left( \\frac{{\\text{{fb}}}}{{\\text{{{unit}}}}} \\right)'.format(unit=_units[var]))
    else:
        yt = _yTitleTemp.format(prefix='\\frac{1}{\\sigma_{\\text{fid}}}',
                                xvar=prettyVar, units='')
        ytnn = _yTitleTemp.format(prefix='', xvar=prettyVar, units='\\left( \\text{fb} \\right)')

    _xTitle[var] = xt
    _yTitle[var] = yt
    _yTitleNoNorm[var] = ytnn

_selections = {
    'pt' : {c:'' for c in _channels},
    'mass' : {c:'' for c in _channels},
    'massFull' : {c:'' for c in _channels},
    'eta' : {c:'' for c in _channels},
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


# do jet variables separately because we have to deal with systematics
for sys in ['', '_jerUp', '_jerDown', '_jesUp','_jesDown']:
    for varName in ['nJets', 'mjj', 'deltaEtajj']:
        doAbs = 'eta' in varName.lower()

        varName += sys
        var = varName
        if doAbs:
            var = 'abs({})'.format(var)

        _variables[varName] = {c:var for c in _channels}

        if 'jj' in varName.lower():
            _selections[varName] = {c:'nJets{} > 1'.format(sys) for c in _channels}
        else:
            _selections[varName] = {c:'' for c in _channels}

    for baseVar in ['jetPt'+sys, 'jetEta'+sys]:
        for j in [1,2]:
            doAbs = 'eta' in baseVar.lower()

            var = baseVar+'[{}]'.format(j-1)
            if doAbs:
                var = 'abs({})'.format(var)
            varName = baseVar.replace('jet','jet{}'.format(j))

            _variables[varName] = {c:var for c in _channels}
            _selections[varName] = {c:'nJets{} >= {}'.format(sys,j) for c in _channels}

_trueSelections = {
    v : {
        'eeee' : 'e1_e2_Mass > 60. && e3_e4_Mass > 60.',
        'eemm' : 'e1_e2_Mass > 60. && m1_m2_Mass > 60.',
        'mmmm' : 'm1_m2_Mass > 60. && m3_m4_Mass > 60.',
        } for v in _selections
    }
_trueSelections['massFull'] = {c:'' for c in _channels}

# Names of compiled C++ classes to make response matrices fast
# (this is extremely slow in Python because it requires a combination of
# information from multiple trees, which can't be done with TTree::Draw())
_responseClassNames = {
    'mass' : {c:'FloatBranchResponseMatrixMaker' for c in _channels},
    'massFull' : {c:'FullSpectrumFloatResponseMatrixMaker' for c in _channels},
    'pt' : {c:'FloatBranchResponseMatrixMaker' for c in _channels},
    'eta' : {c:'AbsFloatBranchResponseMatrixMaker' for c in _channels},
    'nJets' : {c:'JetUIntBranchResponseMatrixMaker' for c in _channels},
    'mjj' : {c:'DijetBranchResponseMatrixMaker' for c in _channels},
    'deltaEtajj' : {c:'AbsDijetBranchResponseMatrixMaker' for c in _channels},
    'z1Mass' : {'eeee':'FloatBranchResponseMatrixMaker',
                'mmmm':'FloatBranchResponseMatrixMaker',
                'eemm':'Z1ByMassResponseMatrixMaker',},
    'z2Mass' : {'eeee':'FloatBranchResponseMatrixMaker',
                'mmmm':'FloatBranchResponseMatrixMaker',
                'eemm':'Z2ByMassResponseMatrixMaker',},
    'z1Pt' : {'eeee':'FloatBranchResponseMatrixMaker',
              'mmmm':'FloatBranchResponseMatrixMaker',
              'eemm':'Z1ByMassResponseMatrixMaker',},
    'z2Pt' : {'eeee':'FloatBranchResponseMatrixMaker',
              'mmmm':'FloatBranchResponseMatrixMaker',
              'eemm':'Z2ByMassResponseMatrixMaker',},
    'zHigherPt' : {c:'Z1ByPtResponseMatrixMaker' for c in _channels},
    'zLowerPt' : {c:'Z2ByPtResponseMatrixMaker' for c in _channels},
    'deltaPhiZZ' : {c:'ZZAbsDeltaPhiResponseMatrixMaker' for c in _channels},
    'deltaRZZ' : {c:'ZZDeltaRResponseMatrixMaker' for c in _channels},
    'lPt' : {c:'AllLeptonBranchResponseMatrixMaker' for c in _channels},
    'l1Pt' : {c:'LeptonMaxBranchResponseMatrixMaker' for c in _channels},
    'zPt' : {c:'BothZsBranchResponseMatrixMaker' for c in _channels},
    'jet1Pt' : {c:'FirstJetFloatResponseMatrixMaker' for c in _channels},
    'jet2Pt' : {c:'SecondJetFloatResponseMatrixMaker' for c in _channels},
    'jet1Eta' : {c:'FirstJetAbsFloatResponseMatrixMaker' for c in _channels},
    'jet2Eta' : {c:'SecondJetAbsFloatResponseMatrixMaker' for c in _channels},
    }

# Variable names usable by response maker classes
_varNamesForResponseMaker = {
    'mass' : {c:'Mass' for c in _channels},
    'massFull' : {c:'Mass' for c in _channels},
    'pt' : {c:'Pt' for c in _channels},
    'eta' : {c:'Eta' for c in _channels},
    'nJets' : {c:'nJets' for c in _channels},
    'mjj' : {c:'mjj' for c in _channels},
    'deltaEtajj' : {c:'deltaEtajj' for c in _channels},
    'z1Mass' : {'eeee':'e1_e2_Mass','mmmm':'m1_m2_Mass','eemm':'Mass'}, # 4e/4mu just use 1 variable because that's easy
    'z2Mass' : {'eeee':'e3_e4_Mass','mmmm':'m3_m4_Mass','eemm':'Mass'}, # for 2e2mu, the response maker class will figure it out
    'z1Pt' : {'eeee':'e1_e2_Pt','mmmm':'m1_m2_Pt','eemm':'Pt'}, # 4e/4mu just use 1 variable because that's easy
    'z2Pt' : {'eeee':'e3_e4_Pt','mmmm':'m3_m4_Pt','eemm':'Pt'}, # for 2e2mu, the response maker class will figure it out
    'zPt' : {c:'Pt' for c in _channels},
    'zHigherPt' : {c:'Pt' for c in _channels},
    'zLowerPt' : {c:'Pt' for c in _channels},
    'deltaPhiZZ' : {c:'deltaPhiZZ' for c in _channels},
    'deltaRZZ' : {c:'deltaRZZ' for c in _channels},
    'lPt' : {c:'Pt' for c in _channels},
    'l1Pt' : {c:'Pt' for c in _channels},
    'zPt' : {c:'Pt' for c in _channels},
    'jet1Pt' : {c:'jetPt' for c in _channels},
    'jet2Pt' : {c:'jetPt' for c in _channels},
    'jet1Eta' : {c:'jetEta' for c in _channels},
    'jet2Eta' : {c:'jetEta' for c in _channels},
}


# list of variables not counting systematic shifts
_varList = [v for v in _variables if 'Up' not in v and 'Down' not in v]
_varListNoFull = _varList[:]
_varListNoFull.remove("massFull")

# Sometimes need to more or resize legend
_legDefaults = {
    'textsize' : .027,#2,
    'leftmargin' : 0.35,
    'entryheight' : 0.037,
    'rightmargin' : 0.03,
    }
_legParams = {v:_legDefaults.copy() for v in _varList}
_legParams['z1Mass'] = {
    'textsize' : .02,
    'leftmargin' : .03,
    'rightmargin' : .46,
    'entryheight' : .034,#23
    'entrysep' : .007,
    }
_legParams['z2Mass'] = _legParams['z1Mass'].copy()
_legParams['deltaRZZ'] = _legParams['z1Mass'].copy()
#_legParams['deltaRZZ']['topmargin'] = 0.7
#_legParams['deltaRZZ']['leftmargin'] = 0.3
#_legParams['deltaRZZ']['rightmargin'] = 0.23
#_legParams['deltaRZZ']['solid'] = True
_legParams['deltaPhiZZ']['leftmargin'] = 0.05
_legParams['deltaPhiZZ']['rightmargin'] = 0.32
_legParams['deltaEtajj'] = _legParams['z1Mass'].copy()
_legParams['deltaEtajj']['leftmargin'] = .5
_legParams['deltaEtajj']['rightmargin'] = .03
_legParams['deltaEtajj']['topmargin'] = .05
_legParams['eta'] = _legParams['deltaEtajj'].copy()
_legParams['jet1Eta']['topmargin'] = 0.058
_legParams['jet2Eta']['topmargin'] = 0.058
_legParams['nJets']['topmargin'] = 0.058
_legParams['massFull']['leftmargin'] = 0.25

_legParamsLogy = {v:p.copy() for v,p in _legParams.iteritems()}
_legParamsLogy['deltaRZZ']['textsize'] = .025
_legParamsLogy['deltaRZZ']['entryheight'] = .03
_legParamsLogy['deltaRZZ']['topmargin'] = .7
_legParamsLogy['deltaRZZ']['leftmargin'] = .18
_legParamsLogy['deltaRZZ']['rightmargin'] = .22
_legParamsLogy['l1Pt']['topmargin'] = 0.05
_legParamsLogy['l1Pt']['leftmargin'] = 0.425
_legParamsLogy['l1Pt']['rightmargin'] = 0.01
_legParamsLogy['l1Pt']['textsize'] = 0.025
_legParamsLogy['mass']['topmargin'] = 0.05
_legParamsLogy['mass']['leftmargin'] = 0.39
_legParamsLogy['mass']['rightmargin'] = 0.025
_legParamsLogy['mass']['textsize'] = 0.025
_legParamsLogy['lPt']['topmargin'] = 0.05
_legParamsLogy['zHigherPt']['topmargin'] = 0.045
_legParamsLogy['massFull']['topmargin'] = 0.035


_uncertaintyTitles = {
    'pu' : 'PU',
    'eEff' : 'Electron eff.',
    'mEff' : 'Muon eff.',
    'generator' : 'Generator choice',
    'lumi' : 'Luminosity',
    'eFR' : 'Electron fake rate',
    'mFR' : 'Muon fake rate',
    'jer' : 'Jet energy res.',
    'jes' : 'Jet energy scale',
    'eScale' : 'Electron energy scale',
    'ePhiRes' : 'Electron energy res. (phi)',
    'eRhoRes' : 'Electron energy res. (rho)',
    'mClosure' : 'Muon calibration',
    'pdf' : 'PDF',
    'alphaS' : '#alpha_{s}',
    'scale' : 'QCD scale',
    'mcfmxsec' : 'MCFM PDF/scale',
    }

_uncertaintyColors = {
    'pu' : 'green',
    'eEff' : 'blue',
    'mEff' : '#002db3',
    'generator' : 'magenta',
    'lumi' : 'orange',
    'eFR' : '#00cc99',
    'mFR' : '#00ff00',
    'jer' : 'cyan',
    'jes' : 'darkblue',
    'eScale' : 'purple',
    'ePhiRes' : 'violet',
    'eRhoRes' : 'lavender',
    'mClosure' : '#c61aff',
    'pdf' : '#80aaff',
    'alphaS' : '#4e72ba',
    'scale' : '#800000',
    'mcfmxsec' : 'red',
    }

# info for adding the MATRIX NNLO curves (has nothign to do with reponse matrix)
_matrixNames = {
    'deltaRZZ' : 'dR.Z_0.25__NNLO_QCD',
    'deltaPhiZZ' : 'dphi.Z_0.25__NNLO_QCD',
    'mass' : 'm.ZZ_5.0__NNLO_QCD',
    'nJets' : 'n_jets__NNLO_QCD',
    'pt' : 'pT.ZZ_2.5__NNLO_QCD',
    'zPt' : 'pT.Z_5.0__NNLO_QCD',
    'zHigherPt' : 'pT.Z_5.0__NNLO_QCD',
    'l1Pt' : 'pT.lep_5.0__NNLO_QCD',
    }
_matrixXSecs = {
    '' : 17.5413,
    'dn' : 17.1878,
    'up' : 17.9555,
    }
_matrixPath = '/data/nawoods/ZZMatrixDistributions'#'/afs/cern.ch/user/k/kelong/www/ZZMatrixDistributions'

_cacheFileTemplate = _join(_env['zzt'], 'Analysis', 'savedResults',
                           'unfoldCache_{}Iter.root')

def _normalizeBins(h):
    binUnit = 1 # min(h.GetBinWidth(b) for b in range(1,len(h)+1))
    for ib in xrange(1,len(h)+1):
        w = h.GetBinWidth(ib)
        h.SetBinContent(ib, h.GetBinContent(ib) * binUnit / w)
        h.SetBinError(ib, h.GetBinError(ib) * binUnit / w)
        if h.GetBinError(ib) > h.GetBinContent(ib):
            h.SetBinError(ib, h.GetBinContent(ib))
    h.sumw2()

def _unnormalizeBins(h):
    binUnit = 1 # min(h.GetBinWidth(b) for b in range(1,len(h)+1))
    for ib in xrange(1,len(h)+1):
        w = h.GetBinWidth(ib)
        h.SetBinContent(ib, h.GetBinContent(ib) * w / binUnit)
        h.SetBinError(ib, h.GetBinError(ib) * w / binUnit)
        if h.GetBinError(ib) > h.GetBinContent(ib):
            h.SetBinError(ib, h.GetBinContent(ib))
    h.sumw2()

_printNext = False
_printCounter = 0
def _getUnfolded(hSig, hBkg, hTrue, hResponse, hData, nIter,
                 withRespAndCov=False, printIt=False):
    global _printNext
    global _printCounter

    response = Response(hSig, hTrue.clone(), hResponse.clone())

    try:
        svd = TDecompSVD(response.Mresponse())
        sig = svd.GetSig()
        try:
            condition = sig.Max() / max(0., sig.Min())
        except ZeroDivisionError:
            condition = float('inf')
            raise

        print ''
        print 'condition: {}'.format(condition)
        print ''

    except:
        print "It broke! Printing debug info"
        print "Sig: {}, bkg: {}, true: {}, response: {}".format(hSig.Integral(), hBkg.Integral(), hTrue.Integral(), hResponse.Integral())
        print "Generating debug plots sig{0}.png, bkg{0}.png, true{0}.png, data{0}.png, resp{0}.png, resp{0}.root".format(_printCounter)
        c = Canvas(1000,1000)
        hSig.draw()
        _style.setCMSStyle(c, '', dataType='Debug', intLumi=35860.)
        c.Print("sig{}.png".format(_printCounter))
        hBkg.draw()
        _style.setCMSStyle(c, '', dataType='Debug', intLumi=35860.)
        c.Print("bkg{}.png".format(_printCounter))
        hTrue.draw()
        _style.setCMSStyle(c, '', dataType='Debug', intLumi=35860.)
        c.Print("true{}.png".format(_printCounter))
        hData.draw()
        _style.setCMSStyle(c, '', dataType='Debug', intLumi=35860.)
        c.Print("data{}.png".format(_printCounter))
        hResponse.drawstyle = 'colz'
        hResponse.draw()
        _style.setCMSStyle(c, '', dataType='Debug', intLumi=35860.)
        c.Print("resp{}.png".format(_printCounter))
        c.Print("resp{}.root".format(_printCounter))
        _printCounter += 1


    hDataMinusBkg = hData - hBkg
    zeroNegativeBins(hDataMinusBkg)

    unf = RooUnfoldIter(response, hDataMinusBkg, nIter)

    if _printNext or printIt:
        _printNext = False

        print "Generating debug plots sig{0}.png, bkg{0}.png, true{0}.png, data{0}.png, resp{0}.png, resp{0}.root".format(_printCounter)
        c = Canvas(1000,1000)
        hSig.draw()
        _style.setCMSStyle(c, '', dataType='Debug', intLumi=35860.)
        c.Print("sig{}.png".format(_printCounter))
        hBkg.draw()
        _style.setCMSStyle(c, '', dataType='Debug', intLumi=35860.)
        c.Print("bkg{}.png".format(_printCounter))
        hTrue.draw()
        _style.setCMSStyle(c, '', dataType='Debug', intLumi=35860.)
        c.Print("true{}.png".format(_printCounter))
        hData.draw()
        _style.setCMSStyle(c, '', dataType='Debug', intLumi=35860.)
        c.Print("data{}.png".format(_printCounter))
        hResponse.drawstyle = 'colz'
        hResponse.draw()
        _style.setCMSStyle(c, '', dataType='Debug', intLumi=35860.)
        c.Print("resp{}.png".format(_printCounter))
        c.Print("resp{}.root".format(_printCounter))
        _printCounter += 1

    hOut = unf.Hreco()
    if not hOut:
        print hOut
        raise ValueError("The unfolded histogram got screwed up somehow!")

    if withRespAndCov:
        return asrootpy(hOut), unf.Ereco(2).Clone(), asrootpy(response.Hresponse()).clone()

    return asrootpy(hOut)


def _generateAnalysisInputs(puWeightFile, looseSIP=False, noSIP=False,
                            sfRemake=False):
    sfFiles = {}
    sipForBkg = 4.
    if noSIP:
        sfFiles['eSelSFFile'] = 'eleSelectionSF_HZZ_NWRemake_NoSIP'
        sfFiles['eSelSFFileGap'] = 'eleSelectionSFGap_HZZ_NWRemake_NoSIP'
        sfFiles['eRecoSFFile'] = 'eleRecoSF_HZZ_Moriond17'
        sfFiles['mSFFile'] = 'muSelectionAndRecoSF_HZZ_Moriond17_NoSIP'
        sipForBkg = -1
        if looseSIP:
            raise ValueError("You can use scale factors for loose SIP cut or "
                             "no SIP cut, but not both.")
    elif looseSIP:
        sfFiles['eSelSFFile'] = 'eleSelectionSF_HZZ_NWRemake_LooseSIP'
        sfFiles['eSelSFFileGap'] = 'eleSelectionSFGap_HZZ_NWRemake_LooseSIP'
        sfFiles['mSFFile'] = 'muSelectionAndRecoSF_HZZ_Moriond17_LooseSIP'
        sfFiles['eRecoSFFile'] = 'eleRecoSF_HZZ_Moriond17'
        sipForBkg = 10.
    elif sfRemake:
        sfFiles['eSelSFFile'] = 'eleSelectionSF_HZZ_NWRemake'
        sfFiles['eSelSFFileGap'] = 'eleSelectionSFGap_HZZ_NWRemake'
        sfFiles['eRecoSFFile'] = 'eleRecoSF_HZZ_Moriond17'
        sfFiles['mSFFile'] = 'muSelectionAndRecoSF_HZZ_Moriond17'

    puWeightFileFull = _join(_env['zzt'],'data','pileup',puWeightFile+'.root')
    with preserve_current_directory():
        with root_open(puWeightFileFull) as fPU:
            hPUWt = {
                '' : asrootpy(fPU.puScaleFactor),
                'up' : asrootpy(fPU.puScaleFactor_up),
                'dn' : asrootpy(fPU.puScaleFactor_down),
                }
            for hpu in hPUWt.values():
                hpu.SetDirectory(0)

    hSF = {}
    if sfFiles:
        with preserve_current_directory():
            with root_open(_join(_env['zzt'],'data','leptonScaleFactors',
                                      sfFiles['eSelSFFile']+'.root')) as fEleSel:
                hSF['eSel'] = asrootpy(fEleSel.EGamma_SF2D).clone()
                hSF['eSel'].SetDirectory(0)
            with root_open(_join(_env['zzt'],'data','leptonScaleFactors',
                                 sfFiles['eSelSFFileGap']+'.root')) as fEleSelGap:
                hSF['eSelGap'] = asrootpy(fEleSelGap.EGamma_SF2D).clone()
                hSF['eSelGap'].SetDirectory(0)
            with root_open(_join(_env['zzt'],'data','leptonScaleFactors',
                                 sfFiles['eRecoSFFile']+'.root')) as fEleReco:
                hSF['eReco'] = asrootpy(fEleReco.EGamma_SF2D).clone()
                hSF['eReco'].SetDirectory(0)
            with root_open(_join(_env['zzt'],'data','leptonScaleFactors',
                                 sfFiles['mSFFile']+'.root')) as fMuSF:
                hSF['m'] = asrootpy(fMuSF.FINAL).clone()
                hSF['m'].SetDirectory(0)
                hSF['mErr'] = asrootpy(fMuSF.ERROR).clone()
                hSF['mErr'].SetDirectory(0)

    return sfFiles, hPUWt, hSF, sipForBkg


def _generateSamples(inData, inMC, ana, fakeRateFile, puWeightFile, lumi,
                     amcatnlo=False, sipForBkg=4., sfFiles={}):
    puWeightStr, puWt = puWeight(puWeightFile, '')
    puWeightStrUp, puWtUp = puWeight(puWeightFile, 'up')
    puWeightStrDn, puWtDn = puWeight(puWeightFile, 'dn')

    allSamples = {}

    allSamples['true'] = genZZSamples('zz', inMC, ana, lumi, amcatnlo=amcatnlo,
                                      higgs=(ana=='full'))
    allSamples['reco'] = zzStackSignalOnly('zz', inMC, ana, puWeightFile,
                                           lumi, amcatnlo=amcatnlo,
                                           asGroup=True, higgs=(ana=='full'),
                                           **sfFiles)

    allSamples['bkgMC'] = zzIrreducibleBkg('zz', inMC, ana, puWeightFile, lumi,
                                           **sfFiles)
    allSamples['bkg'] = standardZZBkg('zz', inData, inMC, ana, puWeightFile,
                                      fakeRateFile, lumi,
                                      sipCut=sipForBkg)
    allSamples['bkgSyst'] = {
        'eup' : standardZZBkg('zz', inData, inMC, ana, puWeightFile,
                              fakeRateFile, lumi, eFakeRateSyst='up',
                              sipCut=sipForBkg),
        'edn' : standardZZBkg('zz', inData, inMC, ana, puWeightFile,
                              fakeRateFile, lumi, eFakeRateSyst='dn',
                              sipCut=sipForBkg),
        'mup' : standardZZBkg('zz', inData, inMC, ana, puWeightFile,
                              fakeRateFile, lumi, mFakeRateSyst='up',
                              sipCut=sipForBkg),
        'mdn' : standardZZBkg('zz', inData, inMC, ana, puWeightFile,
                              fakeRateFile, lumi, mFakeRateSyst='dn',
                              sipCut=sipForBkg),
        }

    allSamples['data'] = standardZZData('zz', inData, ana)

    allSamples['altReco'] = zzStackSignalOnly('zz', inMC, ana, puWeightFile,
                                              lumi, amcatnlo=(not amcatnlo),
                                              asGroup=True,
                                              higgs=(ana=='full'), **sfFiles)

    allSamples['altTrue'] = genZZSamples('zz', inMC, ana, lumi,
                                         amcatnlo=(not amcatnlo),
                                         higgs=(ana=='full'))
    allSamples['recoSyst'] = {}
    allSamples['bkgMCSyst'] = {}
    sigFileNamesSyst = {}
    for syst in ['eScaleUp', 'eScaleDn', 'eRhoResUp',
                 'eRhoResDn', 'ePhiResUp']:
        allSamples['recoSyst'][syst] = zzStackSignalOnly('eeee,eemm',
                                                         inMC.replace('mc_','mc_{}_'.format(syst)),
                                                         ana, puWeightFile,
                                                         lumi,
                                                         amcatnlo=amcatnlo,
                                                         asGroup=True,
                                                         higgs=(ana=='full'),
                                                         **sfFiles)
        allSamples['bkgMCSyst'][syst] = zzIrreducibleBkg('eeee,eemm', inMC.replace('mc_','mc_{}_'.format(syst)),
                                                         ana, puWeightFile,
                                                         lumi, **sfFiles)

    for syst in ['mClosureUp','mClosureDn']:
        allSamples['recoSyst'][syst] = zzStackSignalOnly('eemm,mmmm',
                                                         inMC.replace('mc_','mc_{}_'.format(syst)),
                                                         ana, puWeightFile,
                                                         lumi,
                                                         amcatnlo=amcatnlo,
                                                         asGroup=True,
                                                         higgs=(ana=='full'),
                                                         **sfFiles)
        sigFileNamesSyst[syst] = {s.name : [f for f in s.getFileNames()]
                                  for s in allSamples['recoSyst'][syst].values()[0].getBaseSamples()}
        allSamples['bkgMCSyst'][syst] = zzIrreducibleBkg('eemm,mmmm',
                                                         inMC.replace('mc_','mc_{}_'.format(syst)),
                                                         ana, puWeightFile,
                                                         lumi, **sfFiles)

    return allSamples


def _generateResponseClass(varName, channel, samples, hPUWt, hSF={}):
    className = _responseClassNames[varName][channel]
    if hSF:
        className = 'SFHist'+className

    if not hasattr(_rootComp, className):
        # compile the code and register all the classes
        classesNeeded = set([cn for cv in _responseClassNames.values() for cn in cv.values()])
        classesNeeded |= set(['SFHist'+cn for cn in classesNeeded])

        _rootComp.register_file(_join(_zztBaseDir, 'Utilities',
                                      'ResponseMatrixMaker.cxx'),
                                list(classesNeeded))

    C = getattr(_rootComp, className)

    sigFileNames = {s.name : [f for f in s.getFileNames()]
                    for s in samples['reco'].values()[0].getBaseSamples()}
    sigConstWeights = {s.name : s.xsec * s.intLumi * float(s.kFactor) / s.sumW
                       for s in samples['reco'].values()[0].getBaseSamples()}
    altSigFileNames = {s.name : [f for f in s.getFileNames()]
                       for s in samples['altReco'].values()[0].getBaseSamples()}
    altSigConstWeights = {s.name : s.xsec * s.intLumi * float(s.kFactor) / s.sumW
                          for s in samples['altReco'].values()[0].getBaseSamples()}
    sigFileNamesSyst = {syst : {s.name : [f for f in s.getFileNames()]
                                for s in samples['recoSyst'][syst].values()[0].getBaseSamples()}
                        for syst in samples['recoSyst']}

    binning = _binning[varName]
    vBinning = _VFloat()
    if len(binning) == 3:
        binningTemp = [binning[1] + i * (binning[2] - binning[1])/float(binning[0]) for i in xrange(binning[0]+1)]
        for b in binningTemp:
            vBinning.push_back(b)
    else:
        for b in binning:
            vBinning.push_back(b)

    responseMakers = {}
    for sample, fNameList in sigFileNames.iteritems():
        resp = C(channel, _varNamesForResponseMaker[varName][channel], vBinning)

        for fName in fNameList:
            resp.registerFile(fName)
        for syst in sigFileNamesSyst:
            for fName in sigFileNamesSyst[syst][sample]:
                resp.registerFile(fName, syst)

        resp.registerPUWeights(hPUWt[''])
        resp.registerPUWeights(hPUWt['up'], 'up')
        resp.registerPUWeights(hPUWt['dn'], 'dn')
        resp.setConstantScale(sigConstWeights[sample])
        if hSF:
            resp.registerElectronSelectionSFHist(hSF['eSel'])
            resp.registerElectronSelectionGapSFHist(hSF['eSelGap'])
            resp.registerElectronRecoSFHist(hSF['eReco'])
            resp.registerMuonSFHist(hSF['m'])
            resp.registerMuonSFErrorHist(hSF['mErr'])

        responseMakers[sample] = resp

    altResponseMakers = {}
    for sample, fNameList in altSigFileNames.iteritems():
        if sample in responseMakers:
            continue
        resp = C(channel, _varNamesForResponseMaker[varName][channel], vBinning)

        for fName in fNameList:
            resp.registerFile(fName)
        resp.registerPUWeights(hPUWt[''])
        resp.setConstantScale(altSigConstWeights[sample])
        resp.setSkipSystematics()
        if hSF:
            resp.registerElectronSelectionSFHist(hSF['eSel'])
            resp.registerElectronSelectionGapSFHist(hSF['eSelGap'])
            resp.registerElectronRecoSFHist(hSF['eReco'])
            resp.registerMuonSFHist(hSF['m'])
            resp.registerMuonSFErrorHist(hSF['mErr'])

        altResponseMakers[sample] = resp

    return responseMakers, altResponseMakers


def _unfold(varName, chan, samples, puWeightFile, sfFiles,
            responseMakers, altResponseMakers, nIter, plotDir=''):
    # outputs
    hUnfolded = {}
    hTrue = {}
    hTrueAlt = {}

    var = _variables[varName][chan]
    sel = _selections[varName][chan]
    if isinstance(sel,str):
        selTrue = combineWeights(sel, _trueSelections[varName][chan], selections=True)
    else:
        selTrue = [combineWeights(s, _trueSelections[varName][chan], selections=True) for s in sel]

    binning = _binning[varName]
    hData = samples['data'][chan].makeHist(var, sel, binning, perUnitWidth=False)

    # regular weight, no systematics. Apply just in case.
    nominalWeight = baseMCWeight(chan, puWeightFile,
                                 **sfFiles)
    samples['reco'][chan].applyWeight(nominalWeight, True)
    samples['altReco'][chan].applyWeight(nominalWeight, True)
    samples['bkgMC'][chan].applyWeight(nominalWeight, True)
    for s in samples['recoSyst'].values():
        try:
            s[chan].applyWeight(nominalWeight, True)
        except KeyError:
            pass

    hTrue[''] = samples['true'][chan].makeHist(var, selTrue, binning,
                                               perUnitWidth=False)
    hSigNominal = samples['reco'][chan].makeHist(var, sel, binning, perUnitWidth=False)
    hBkgMCNominal = samples['bkgMC'][chan].makeHist(var, sel, binning, perUnitWidth=False)
    hBkgNominal = samples['bkg'][chan].makeHist(var, sel, binning, perUnitWidth=False,
                                                postprocess=True)
    hResponseNominal = {s:asrootpy(resp()) for s,resp in responseMakers.iteritems()}
    hResponseNominalTotal = sum(resp for resp in hResponseNominal.values())


    hUnfolded[''], hCov, hResp = _getUnfolded(hSigNominal,
                                              hBkgMCNominal+hBkgNominal,
                                              hTrue[''],
                                              hResponseNominalTotal,
                                              hData, nIter, True)

    # plot covariance and response
    if plotDir:
        cRes = Canvas(1000,1000)
        if varName == 'massFull':
            cRes.SetLogx()
            cRes.SetLogy()
        hResp.drawstyle = 'colztext'
        hResp.xaxis.title = '\\text{Reco} '+_xTitle[varName]
        hResp.yaxis.title = '\\text{True} '+_xTitle[varName]
        hResp.draw()
        _style.setCMSStyle(cRes, '', dataType='Internal', intLumi=35860.)
        cRes.Print(_join(plotDir, 'pngs', "response_{}_{}.png".format(varName, chan)))
        cRes.Print(_join(plotDir, 'Cs', "response_{}_{}.C".format(varName, chan)))

        cCov = Canvas(1000,1000)
        if varName == 'massFull':
            cCov.SetLogx()
            cCov.SetLogy()
        hCov.Draw("colztext")
        _style.setCMSStyle(cCov, '', dataType='Internal', intLumi=35860.)
        cCov.Print(_join(plotDir, 'pngs', "covariance_{}_{}.png".format(varName, chan)))
        cCov.Print(_join(plotDir, 'Cs', "covariance_{}_{}.C".format(varName, chan)))


    # PU reweight uncertainty
    for sys in ['up','dn']:
        wtStr = baseMCWeight(chan, puWeightFile, puSyst=sys,
                             **sfFiles)
        samples['reco'][chan].applyWeight(wtStr, True)
        samples['bkgMC'][chan].applyWeight(wtStr, True)

        hSig = samples['reco'][chan].makeHist(var, sel, binning, perUnitWidth=False)
        hBkgMC = samples['bkgMC'][chan].makeHist(var, sel, binning, perUnitWidth=False)

        hResponse = sum(asrootpy(resp('pu_'+sys)) for resp in responseMakers.values())

        hUnfolded['pu_'+sys] = _getUnfolded(hSig, hBkgMC+hBkgNominal,
                                            hTrue[''], hResponse,
                                            hData, nIter)

        samples['reco'][chan].applyWeight(nominalWeight, True)
        samples['bkgMC'][chan].applyWeight(nominalWeight, True)


    # lepton efficiency uncertainty
    for lep in set(chan):
        for sys in ['up','dn']:
            wtArg = {lep+'Syst':sys}
            wtArg.update(sfFiles)
            wtStr = baseMCWeight(chan, puWeightFile, **wtArg)
            samples['reco'][chan].applyWeight(wtStr, True)
            samples['bkgMC'][chan].applyWeight(wtStr, True)

            hSig = samples['reco'][chan].makeHist(var, sel, binning, perUnitWidth=False)
            hBkgMC = samples['bkgMC'][chan].makeHist(var, sel, binning, perUnitWidth=False)

            hResponse = sum(asrootpy(resp(lep+'Eff_'+sys)) for resp in responseMakers.values())

            hUnfolded[lep+'Eff_'+sys] = _getUnfolded(hSig,
                                                     hBkgMC+hBkgNominal,
                                                     hTrue[''],
                                                     hResponse,
                                                     hData, nIter)

            samples['reco'][chan].applyWeight(nominalWeight, True)
            samples['bkgMC'][chan].applyWeight(nominalWeight, True)

    # alternate generator
    hSig = samples['altReco'][chan].makeHist(var, sel, binning,
                                             perUnitWidth=False)
    hTrueAlt[''] = samples['altTrue'][chan].makeHist(var, selTrue, binning,
                                                     perUnitWidth=False)
    hResponses = []
    altSigFileNames = {s.name : [f for f in s.getFileNames()]
                       for s in samples['altReco'].values()[0].getBaseSamples()}
    for s in altSigFileNames.keys():
        try:
            hResponses.append(asrootpy(altResponseMakers[s]()))
        except KeyError:
            hResponses.append(hResponseNominal[s])
    hResponse = sum(h for h in hResponses)

    hUnfolded['generator'] = _getUnfolded(hSig,
                                          hBkgMCNominal+hBkgNominal,
                                          hTrueAlt[''],
                                          hResponse,
                                          hData, nIter)

    # luminosity
    lumiUnc = 0.025
    lumiScale = {'up':1.+lumiUnc,'dn':1.-lumiUnc}
    for sys, scale in lumiScale.iteritems():
        hSig = hSigNominal * scale
        hBkgMC = hBkgMCNominal * scale
        hTrueLumiShift = hTrue[''] * scale
        hResponse = hResponseNominalTotal * scale

        hUnfolded['lumi_'+sys] = _getUnfolded(hSig,
                                              hBkgMC+hBkgNominal,
                                              hTrueLumiShift,
                                              hResponse,
                                              hData, nIter)

    # lepton fake rate uncertainty
    for lep in set(chan):
        for sys in ['up','dn']:
            hBkg = samples['bkgSyst'][lep+sys][chan].makeHist(var, sel, binning,
                                                              perUnitWidth=False,
                                                              postprocess=True)

            hUnfolded[lep+'FR_'+sys] = _getUnfolded(hSigNominal,
                                                    hBkgMCNominal+hBkg,
                                                    hTrue[''],
                                                    hResponseNominalTotal,
                                                    hData, nIter)

    # jet stuff
    if 'jet' in varName.lower() or 'jj' in varName.lower():
        for shift in ['up','dn']:
            sysStr = 'Up' if shift == 'up' else 'Down'

            for sys in ['jer','jes']:
                shiftedVarName = varName + '_' + sys + sysStr
                varShifted = _variables[shiftedVarName][chan]
                selShifted = _selections[shiftedVarName][chan]

                hSig = samples['reco'][chan].makeHist(varShifted, selShifted,
                                                      binning,
                                                      perUnitWidth=False)
                hBkgMC = samples['bkgMC'][chan].makeHist(varShifted, selShifted,
                                                         binning,
                                                         perUnitWidth=False)

                hResponse = sum(asrootpy(resp(sys+'_'+shift)) for resp in responseMakers.values())

                hUnfolded[sys+'_'+shift] = _getUnfolded(hSig,
                                                        hBkgMC+hBkgNominal,
                                                        hTrue[''],
                                                        hResponse,
                                                        hData, nIter)

    # lepton momentum uncertainties
    if 'e' in chan:
        for sys in ['eScale', 'eRhoRes', 'ePhiRes']:
            for shift in ['up','dn']:
                if sys == 'ePhiRes' and shift == 'dn':
                    continue
                sysStr = 'Up' if shift == 'up' else 'Dn'

                hSig = samples['recoSyst'][sys+sysStr][chan].makeHist(var, sel,
                                                                      binning,
                                                                      perUnitWidth=False)
                hBkgMC = samples['bkgMCSyst'][sys+sysStr][chan].makeHist(var, sel,
                                                                         binning,
                                                                         perUnitWidth=False)
                hResponse = sum(asrootpy(resp(sys+'_'+shift)) for resp in responseMakers.values())

                storeAs = sys+'_'+shift
                if sys == 'ePhiRes':
                    storeAs = sys
                hUnfolded[storeAs] = _getUnfolded(hSig,
                                                  hBkgMC+hBkgNominal,
                                                  hTrue[''],
                                                  hResponse,
                                                  hData, nIter)
    if 'm' in chan:
        sys = 'mClosure'
        for shift in ['up','dn']:
            sysStr = 'Up' if shift == 'up' else 'Dn'

            hSig = samples['recoSyst'][sys+sysStr][chan].makeHist(var, sel,
                                                                  binning,
                                                                  perUnitWidth=False)
            hBkgMC = samples['bkgMCSyst'][sys+sysStr][chan].makeHist(var, sel,
                                                                     binning,
                                                                     perUnitWidth=False)
            hResponse = sum(asrootpy(resp(sys+'_'+shift)) for resp in responseMakers.values())

            hUnfolded[sys+'_'+shift] = _getUnfolded(hSig,
                                                    hBkgMC+hBkgNominal,
                                                    hTrue[''],
                                                    hResponse,
                                                    hData, nIter)

    # PDF uncertainties
    hSigVariations = []
    for s in samples['reco'][chan].values():
        if 'GluGluZZ' not in s.name and 'phantom' not in s.name:
            hSigVariations.append(s.makeHist2(var, 'Iteration$', sel, binning,
                                              [100,0.,100.], 'pdfWeights/pdfWeights[0]', False))
    hResponseVariations = []
    for s, resp in responseMakers.iteritems():
        if "GluGluZZ" not in s and 'phantom' not in s:
            hResponseVariations.append(asrootpy(resp.getPDFResponses()))

    # for each var bin in each sample, get the RMS across all the variations
    allSigRMSes = [[Graph(h.ProjectionY('slice{}'.format(i), i+1,i+1)).GetRMS(2) for i in xrange(h.GetNbinsX())] for h in hSigVariations]
    allResponseRMSes = [[[Graph(h.ProjectionZ('slice_{}_{}'.format(x,y), x+1, x+1, y+1, y+1)).GetRMS(2)
                          for y in xrange(h.GetNbinsY())]
                         for x in xrange(h.GetNbinsX())]
                        for h in hResponseVariations]
    # for each var bin, add variations for all samples
    sigBinRMSes = [sum(rmses) for rmses in zip(*allSigRMSes)]
    responseBinRMSes = [[sum(rmses) for rmses in zip(*colForAllSamples)] for colForAllSamples in zip(*allResponseRMSes)]

    hSigUp = hSigNominal.clone()
    hSigDn = hSigNominal.clone()
    hResponseUp = hResponseNominalTotal.clone()
    hResponseDn = hResponseNominalTotal.clone()

    # apply variations
    for i in xrange(hSigUp.GetNbinsX()):
        hSigUp[i+1].value += sigBinRMSes[i]
        hSigDn[i+1].value = max(0.,hSigDn[i+1].value - sigBinRMSes[i])
    for x in xrange(hResponseUp.GetNbinsX()):
        for y in xrange(hResponseUp.GetNbinsY()):
            hResponseUp[x+1,y+1].value += responseBinRMSes[x][y]
            hResponseDn[x+1,y+1].value = max(0., hResponseDn[x+1,y+1].value - responseBinRMSes[x][y])

    hTrue['pdf_up'] = hTrue[''].clone()
    hTrue['pdf_dn'] = hTrue[''].clone()
    hTrueVariations = []
    for s in samples['true'][chan].values():
        if 'GluGluZZ' not in s.name and 'phantom' not in s.name:
            hTrueVariations.append(s.makeHist2(var, 'Iteration$', selTrue, binning,
                                               [100,0.,100.], 'pdfWeights/pdfWeights[0]', False))
    allTrueRMSes = [[Graph(h.ProjectionY('slice{}'.format(i), i+1,i+1)).GetRMS(2) for i in xrange(h.GetNbinsX())] for h in hTrueVariations]
    binTrueRMSes = [sum(rmses) for rmses in zip(*allTrueRMSes)]

    #hTruePDFErr[chan] = hTrue.empty_clone() # save true variation for later
    for i in xrange(hTrue['pdf_up'].GetNbinsX()):
        hTrue['pdf_up'][i+1].value += binTrueRMSes[i]
        hTrue['pdf_dn'][i+1].value = max(0.,hTrue['pdf_dn'][i+1].value - binTrueRMSes[i])
        #hTruePDFErr[chan][i+1].value = binTrueRMSes[i]

    # get the other sample's uncertainty too as long as we're at it
    #hTruePDFErrAlt[chan] = hTrueAlt.empty_clone()
    hTrueAlt['pdf_up'] = hTrueAlt[''].clone()
    hTrueAlt['pdf_dn'] = hTrueAlt[''].clone()
    hTrueVariationsAlt = []
    for s in samples['altTrue'][chan].values():
        if 'GluGluZZ' not in s.name and 'phantom' not in s.name:
            hTrueVariationsAlt.append(s.makeHist2(var, 'Iteration$', selTrue, binning,
                                                  [100,0.,100.], 'pdfWeights/pdfWeights[0]', False))
    allTrueRMSesAlt = [[Graph(h.ProjectionY('slice{}'.format(i), i+1,i+1)).GetRMS(2) for i in xrange(h.GetNbinsX())] for h in hTrueVariationsAlt]
    binTrueRMSesAlt = [sum(rmses) for rmses in zip(*allTrueRMSesAlt)]
    for i in xrange(hTrueAlt['pdf_up'].GetNbinsX()):
        hTrueAlt['pdf_up'][i+1].value += binTrueRMSesAlt[i]
        hTrueAlt['pdf_dn'][i+1].value = max(0.,hTrueAlt['pdf_dn'][i+1].value - binTrueRMSesAlt[i])
        #hTruePDFErrAlt[chan][i+1].value = binTrueRMSesAlt[i]


    hUnfolded['pdf_up'] = _getUnfolded(hSigUp,
                                       hBkgMCNominal+hBkgNominal,
                                       hTrue['pdf_up'], hResponseUp,
                                       hData, nIter)
    hUnfolded['pdf_dn'] = _getUnfolded(hSigDn,
                                       hBkgMCNominal+hBkgNominal,
                                       hTrue['pdf_dn'], hResponseDn,
                                       hData, nIter)

    # QCD scale uncertainties
    variationIndices = [1,2,3,4,6,8]
    hSigs = [samples['reco'][chan].makeHist(var, sel, binning,
                                            {
                'ZZTo4L':'scaleWeights[{}]/scaleWeights[0]'.format(i),
                'ZZTo4L-amcatnlo':'scaleWeights[{}]/scaleWeights[0]'.format(i),
                'ZZJJTo4L_EWK':'scaleWeights[{}]/scaleWeights[0]'.format(i),
                },
                                            perUnitWidth=False)
             for i in variationIndices]

    hTrues = [samples['true'][chan].makeHist(var, selTrue, binning,
                                             {
                'ZZTo4L':'scaleWeights[{}]/scaleWeights[0]'.format(i),
                'ZZTo4L-amcatnlo':'scaleWeights[{}]/scaleWeights[0]'.format(i),
                'ZZJJTo4L_EWK':'scaleWeights[{}]/scaleWeights[0]'.format(i),
                },
                                             perUnitWidth=False)
              for i in variationIndices]
    nominalArea = hTrue[''].Integral(0,hTrue[''].GetNbinsX()+1)
    for h in hTrues:
        h *= nominalArea / h.Integral(0,h.GetNbinsX()+1)

    # save true-level uncertainty for later
    hTrue['scale_up'] = hTrue[''].empty_clone()
    hTrue['scale_dn'] = hTrue[''].empty_clone()
    for bUp, bDn, variations in zip(hTrue['scale_up'],
                                    hTrue['scale_dn'],
                                    zip(*hTrues)):
        bUp.value = max(b.value for b in variations)
        bDn.value = min(b.value for b in variations)

    # get the true-level uncertainty too while we're at it
    hTruesAlt = [samples['altTrue'][chan].makeHist(var, selTrue, binning,
                                        {
                'ZZTo4L':'scaleWeights[{}]/scaleWeights[0]'.format(i),
                'ZZTo4L-amcatnlo':'scaleWeights[{}]/scaleWeights[0]'.format(i),
                'ZZJJTo4L_EWK':'scaleWeights[{}]/scaleWeights[0]'.format(i),
                },
                                        perUnitWidth=False)
                 for i in variationIndices]
    nominalAreaAlt = hTrueAlt[''].Integral(0,hTrueAlt[''].GetNbinsX()+1)
    for h in hTruesAlt:
        h *= nominalAreaAlt / h.Integral(0,h.GetNbinsX()+1)

    hTrueAlt['scale_up'] = hTrueAlt[''].empty_clone()
    hTrueAlt['scale_dn'] = hTrueAlt[''].empty_clone()
    for bUp, bDn, variations in zip(hTrueAlt['scale_up'],
                                    hTrueAlt['scale_dn'],
                                    zip(*hTruesAlt)):
        bUp.value = max(b.value for b in variations)
        bDn.value = min(b.value for b in variations)

    hResponseVariations = [hResponseNominalTotal.empty_clone() for v in variationIndices]
    for s, resp in responseMakers.iteritems():
        vResponses = resp.getScaleResponses()
        if vResponses.size() == len(hResponseVariations):
            for iResp in xrange(vResponses.size()):
                hResponseVariations[iResp] += asrootpy(vResponses.at(iResp))
        else:
            for hrv in hResponseVariations:
                hrv += hResponseNominal[s]

    hUnfoldedVariations = []
    for hSig, hTr, hResponse in zip(hSigs, hTrues, hResponseVariations):
        hUnfoldedVariations.append(_getUnfolded(hSig,
                                                hBkgMCNominal+hBkgNominal,
                                                hTr, hResponse,
                                                hData, nIter))

    hUnfoldedUp = hUnfoldedVariations[0].empty_clone()
    hUnfoldedDn = hUnfoldedVariations[0].empty_clone()
    for bUp, bDn, bVars in zip(hUnfoldedUp, hUnfoldedDn, zip(*hUnfoldedVariations)):
        bUp.value = max(b.value for b in bVars)
        bDn.value = min(b.value for b in bVars)

    hUnfolded['scale_up'] = hUnfoldedUp
    hUnfolded['scale_dn'] = hUnfoldedDn

    # alpha_s uncertainties
    alphaSIndices = [100,101]
    hSigs = [samples['reco'][chan].makeHist(var, sel, binning,
                         {
                'ZZTo4L':'pdfWeights[{}]/pdfWeights[0]'.format(i),
                'ZZTo4L-amcatnlo':'pdfWeights[{}]/pdfWeights[0]'.format(i),
                'ZZJJTo4L_EWK':'pdfWeights[{}]/pdfWeights[0]'.format(i),
                },
                         perUnitWidth=False)
            for i in alphaSIndices]
    hTrues = [samples['true'][chan].makeHist(var, selTrue, binning,
                          {
                'ZZTo4L':'pdfWeights[{}]/pdfWeights[0]'.format(i),
                'ZZTo4L-amcatnlo':'pdfWeights[{}]/pdfWeights[0]'.format(i),
                'ZZJJTo4L_EWK':'pdfWeights[{}]/pdfWeights[0]'.format(i),
                },
                          perUnitWidth=False)
             for i in alphaSIndices]

    hTrue['alphaS_up'] = hTrues[0]
    hTrue['alphaS_dn'] = hTrues[1]

    hTruesAlt = [samples['altTrue'][chan].makeHist(var, selTrue, binning,
                                                   {
                'ZZTo4L':'pdfWeights[{}]/pdfWeights[0]'.format(i),
                'ZZTo4L-amcatnlo':'pdfWeights[{}]/pdfWeights[0]'.format(i),
                'ZZJJTo4L_EWK':'pdfWeights[{}]/pdfWeights[0]'.format(i),
                },
                                                   perUnitWidth=False)
                 for i in alphaSIndices]

    hTrueAlt['alphaS_up'] = hTruesAlt[0]
    hTrueAlt['alphaS_dn'] = hTruesAlt[1]

    hResponses = [hResponseNominalTotal.empty_clone(),
                  hResponseNominalTotal.empty_clone()]
    for s, resp in responseMakers.iteritems():
        if resp.hasSystematic('alphaS_up'):
            hResponses[0] += asrootpy(resp('alphaS_up'))
            hResponses[1] += asrootpy(resp('alphaS_dn'))
        else:
            hResponses[0] += hResponseNominal[s]
            hResponses[1] += hResponseNominal[s]

    hUnfUp = _getUnfolded(hSigs[0], hBkgNominal+hBkgMCNominal,
                          hTrues[0], hResponses[0],
                          hData, nIter)
    hUnfDn = _getUnfolded(hSigs[1], hBkgNominal+hBkgMCNominal,
                          hTrues[1], hResponses[1],
                          hData, nIter)

    unc = hUnfUp - hUnfDn
    unc /= 2.
    for b in unc:
        b.value = abs(b.value)

    hUnfolded['alphaS_up'] = hUnfolded[''] + unc
    hUnfolded['alphaS_dn'] = hUnfolded[''] - unc

    # since MCFM samples don't have LHE information, we just vary by
    # the cross section uncertainties
    mcfmUnc = {'up':.18,'dn':-.15}
    for sys, shift in mcfmUnc.iteritems():
        hSig = samples['reco'][chan].makeHist(var, sel, binning,
                                              {'GluGluZZ':str(1.+shift)},
                                              perUnitWidth=False)
        hTrue['mcfmxsec_'+sys] = samples['true'][chan].makeHist(var, selTrue, binning,
                                                                {'GluGluZZ':str(1.+shift)},
                                                                perUnitWidth=False)
        hTrueAlt['mcfmxsec_'+sys] = samples['altTrue'][chan].makeHist(var, selTrue, binning,
                                                                      {'GluGluZZ':str(1.+shift)},
                                                                      perUnitWidth=False)
        hResponse = hResponseNominalTotal.empty_clone()
        for s, h in hResponseNominal.iteritems():
            if 'GluGluZZ' in s:
                hResponse += h * (1.+shift)
            else:
                hResponse += h

        hUnfolded['mcfmxsec_'+sys] = _getUnfolded(hSig,
                                                  hBkgNominal+hBkgMCNominal,
                                                  hTrue['mcfmxsec_'+sys],
                                                  hResponse,
                                                  hData, nIter)

    # make everything local (we'll cache copies)
    for h in hUnfolded.values()+hTrue.values()+hTrueAlt.values():
        h.SetDirectory(0)

    return hUnfolded, hTrue, hTrueAlt


def _sumUncertainties(errDict):
    hUncUp = errDict['up'].values()[0].empty_clone()
    hUncDn = errDict['dn'].values()[0].empty_clone()
    sysList = errDict['up'].keys()
    for bUncUp, bUncDn, allUncUp, allUncDn in zip(hUncUp, hUncDn,
                                                  zip(*[errDict['up'][sys] for sys in sysList]),
                                                  zip(*[errDict['dn'][sys] for sys in sysList])):
        for b1,b2 in zip(allUncUp,allUncDn):
            bUncUp.value += max(b1.value,b2.value)**2
            bUncDn.value += min(b1.value,b2.value)**2

        bUncUp.value = sqrt(bUncUp.value)
        bUncDn.value = sqrt(bUncDn.value)

    return hUncUp, hUncDn


def _combineChannelUncertainties(*errDicts):
    hUncTot = {}
    uncList = []
    for errDict in errDicts:
        for sys in ['up','dn']:
            uncList += errDict[sys].keys()
    uncList = set(uncList)

    for sys in ['up','dn']:
        hUncTot[sys] = {}
        for unc in uncList:
            hUncTot[sys][unc] = errDicts[0][sys].values()[0].empty_clone()
            for errDict in errDicts:
                try:
                    hUncTot[sys][unc] += errDict[sys][unc]
                except KeyError:
                    continue

    return hUncTot


def _generateUncertainties(hDict, norm, **plotArgs):
    plot = bool(plotArgs) # only plot if told to
    if plot:
        lumi = plotArgs.get('lumi', 35860.)
        varName = plotArgs['varName']
        ana = plotArgs['ana']
        plotDir = plotArgs['plotDir']
        chan = plotArgs['chan']

    nominalArea = hDict[''].Integral(0,hDict[''].GetNbinsX()+1)
    hErr = {'up':{},'dn':{}}
    for sys, h in hDict.iteritems():
        if not sys:
            continue

        he = h.clone()

        if norm:
            he *= nominalArea / he.Integral(0,he.GetNbinsX()+1)
        # we already shifted the response matrix for lumi, but not the
        # final normalization
        elif sys == 'lumi_up':
            he /= 1.025
        elif sys == 'lumi_dn':
            he /= 0.975

        he -= hDict['']
        sysName = sys.replace('_up','').replace('_dn','')
        if plot:
            he.title = _uncertaintyTitles[sysName]
            he.color = _uncertaintyColors[sysName]
            he.fillstyle = 'solid'
            he.drawstyle = 'hist'
            he.legendstyle = 'F'

        if '_up' in sys:
            hErr['up'][sysName] = he
        elif '_dn' in sys:
            hErr['dn'][sysName] = he
        else:
            hErr['up'][sysName] = he
            he2 = he.clone()
            if plot:
                he2.title = _uncertaintyTitles[sysName]
                he2.color = _uncertaintyColors[sysName]
                he2.fillstyle = 'solid'
                he2.drawstyle = 'hist'
                he2.legendstyle = 'F'
            hErr['dn'][sysName] = he2


    # Make all error histograms positive
    # and make uncertainties fractional (as a percentage)
    if plot:
        for sys in hErr.values():
            for h in sys.values():
                h /= hDict['']
                h *= 100.
                for b in h:
                    b.value = abs(b.value)

        statErr = hDict[''].empty_clone()
        for bUnf, bStat in zip(hDict[''], statErr):
            bStat.value = bUnf.error
        statErr /= hDict['']
        statErr *= 100.
        statErr.color = 'lightgrey'
        statErr.fillstyle = 'solid'
        statErr.legendstyle = 'F'
        statErr.drawstyle = 'hist'
        statErr.title = 'Stat/unfolding'
        errListUp = [statErr]+list(hErr['up'].values())
        errListDn = [statErr]+list(hErr['up'].values())

        # quadrature sum of errors to put on top
        totErrUp = statErr.empty_clone()
        totErrDn = statErr.empty_clone()
        for bTot, bs in zip(totErrUp, zip(*errListUp)):
            bTot.value = sqrt(sum(b.value**2 for b in bs))
        for bTot, bs in zip(totErrDn, zip(*errListDn)):
            bTot.value = sqrt(sum(b.value**2 for b in bs))
        totErrUp.title = 'Total (quadrature sum)'
        totErrUp.color = 'black'
        totErrUp.fillstyle = 'hollow'
        totErrUp.drawstyle = 'hist'
        totErrUp.legendstyle = 'L'
        totErrUp.SetLineWidth(3*totErrUp.GetLineWidth())
        totErrDn.title = 'Total (sum of squares)'
        totErrDn.color = 'black'
        totErrDn.fillstyle = 'hollow'
        totErrDn.drawstyle = 'hist'
        totErrDn.legendstyle = 'L'
        totErrDn.SetLineWidth(3*totErrDn.GetLineWidth())

        # Make plots of uncertainties (added linearly)
        cErrUp = Canvas(1000,1000)
        drawOpts = {
            'xtitle' : _xTitle[varName],
            'ytitle' : "+Error (%)",
            'yerror_in_padding' : False,
            }
        if ana == 'full':
            drawOpts['logx'] = True
        if varName == 'deltaRZZ':
            if chan == 'eeee':
                drawOpts['ylimits'] = (0.,140.)
            elif chan == 'mmmm':
                drawOpts['ylimits'] = (0.,75.)
            else:
                drawOpts['ylimits'] = (0.,58.)
        errStackUp = HistStack(errListUp, drawstyle = 'histnoclear')
        draw([errStackUp,totErrUp], cErrUp, **drawOpts)
        leg = makeLegend(cErrUp, *(errListUp+[totErrUp]), leftmargin=0.25,
                         entryheight=.02, entrysep=.007, textsize=.022,
                         rightmargin=.25)
        leg.Draw('same')
        _style.setCMSStyle(cErrUp, '', dataType='Internal', intLumi=lumi)
        cErrUp.Print(_join(plotDir, 'pngs', 'errUp_{}_{}.png'.format(varName, chan)))
        cErrUp.Print(_join(plotDir, 'Cs', 'errUp_{}_{}.C'.format(varName, chan)))

        cErrDn = Canvas(1000,1000)
        drawOpts = {
            'xtitle' : _xTitle[varName],
            'ytitle' : "-Error (%)",
            'yerror_in_padding' : False,
            }
        if ana == 'full':
            drawOpts['logx'] = True
        if varName == 'deltaRZZ' and chan == 'eeee':
            drawOpts['ylimits'] = (0.,140.)
        errStackDn = HistStack(errListDn, drawstyle = 'histnoclear')
        draw([errStackDn, totErrDn], cErrDn, **drawOpts)
        leg = makeLegend(cErrDn, *(errListDn+[totErrDn]), leftmargin=0.25,
                         entryheight=.02, entrysep=.007, textsize=.022,
                         rightmargin=.25)
        leg.Draw('same')
        _style.setCMSStyle(cErrDn, '', dataType='Internal', intLumi=lumi)
        cErrDn.Print(_join(plotDir, 'pngs', 'errDown_{}_{}.png'.format(varName, chan)))
        cErrDn.Print(_join(plotDir, 'Cs', 'errDown_{}_{}.C'.format(varName, chan)))

        # Errors should no longer be fractional or a percentage
        for sys in hErr.values():
            for h in sys.values():
                h *= hDict['']
                h /= 100.

    return hErr


def _generatePlots(hUnfolded, hUncUp, hUncDn,
                   hTrue, hTrueUncUp, hTrueUncDn,
                   hTrueAlt, hTrueUncUpAlt, hTrueUncDnAlt,
                   varName, chan, ana, plotDir, plotType, lumi, norm, logy,
                   amcatnlo):
    # for normalization if needed
    nominalArea = hUnfolded.Integral(0,hUnfolded.GetNbinsX()+1)
    # Make uncertainties out of the unfolded histos
    ### plot
    hUnf = hUnfolded.clone()

    lumifb = lumi / 1000.

    if norm:
        hUnf /= hUnf.Integral(0,hUnf.GetNbinsX()+1)
    else:
        hUnf /= lumifb

    if logy:
        legParams = _legParamsLogy
    else:
        legParams = _legParams

    hUnf.color = 'black'
    hUnf.drawstyle = 'PE1'
    hUnf.legendstyle = 'LPE1'
    hUnf.title = '\\textbf{Data + stat.\ unc.}'
    if not norm:
        print "Inclusive {} fiducial cross section = {} fb".format(chan, hUnf.Integral(0,hUnf.GetNbinsX()+1))
    _normalizeBins(hUnf)

    signalName = 'POWHEG+MCFM+Pythia8'
    signalNameAlt = 'MG5\_aMC@NLO+MCFM'
    if ana == 'full':
        signalNameAlt += '+POWHEG'
    signalNameAlt += '+Pythia8'
    if amcatnlo:
        signalName, signalNameAlt = signalNameAlt, signalName

    signalName = r'\textbf{'+signalName+'}'
    signalNameAlt = r'\textbf{'+signalNameAlt+'}'

    hTrueDraw = hTrue.clone()
    hTrueDraw.fillcolor = '#99ccff'
    hTrueDraw.linecolor = '#000099'
    hTrueDraw.drawstyle = 'hist'
    TAttFill.SetFillStyle(hTrueDraw, 3004)
    hTrueDraw.SetLineWidth(2*hTrueDraw.GetLineWidth())
    #hTrue.fillstyle = '/'#'solid'
    hTrueDraw.title = '{}'.format(signalName)

    if norm:
        trueInt = hTrueDraw.Integral(0,hTrueDraw.GetNbinsX()+1)
        hTrueDraw /= trueInt
        hTrueUncUp /= trueInt # (trueInt + hTrueUncUp.Integral(0,hTrueUncUp.GetNbinsX()+1))
        hTrueUncDn /= trueInt # (trueInt - hTrueUncDn.Integral(0,hTrueUncDn.GetNbinsX()+1))
    else:
        hTrueDraw /= lumifb
        hTrueUncUp /= lumifb
        hTrueUncDn /= lumifb

    _normalizeBins(hTrueDraw)
    _normalizeBins(hTrueUncUp)
    _normalizeBins(hTrueUncDn)

    # true uncertainty band
    errorBandTrue = makeErrorBand(hTrueDraw, hTrueUncUp, hTrueUncDn)
    errorBandTrue.fillstyle = 'solid'
    errorBandTrue.SetFillColorAlpha(600, 0.5)

    hTrueLeg = Hist(1, 0, 1)
    hTrueLeg.fillcolor = hTrueDraw.fillcolor
    hTrueLeg.linecolor = '#7f7fff' # color transparent error band looks like on white background, found with GIMP eyedropper
    TAttFill.SetFillStyle(hTrueLeg, 3004)
    hTrueLeg.SetLineWidth(5*hTrueDraw.GetLineWidth())
    hTrueLeg.title = hTrueDraw.title
    hTrueLeg.legendstyle = 'FL'

    # use later to put fill and error bars on same legend entry
    trueLegOverlay = Hist(1, 0, 1, title='', legendstyle='L',
                          linecolor=hTrueDraw.linecolor)
    trueLegOverlay.SetLineWidth(hTrueDraw.GetLineWidth())

    toPlot = [hTrueDraw, errorBandTrue]
    forLegend = [hTrueLeg]#Leg]


    hTrueDrawAlt = hTrueAlt.clone()
    hTrueDrawAlt.color = 'red'
    hTrueDrawAlt.drawstyle = 'hist'
    hTrueDrawAlt.fillstyle = 'hollow'
    hTrueDrawAlt.linestyle = 'verylongdash'
    hTrueDrawAlt.SetLineWidth(hTrueDrawAlt.GetLineWidth()*2)
    hTrueDrawAlt.title = '{}'.format(signalNameAlt)

    if norm:
        trueIntAlt = hTrueDrawAlt.Integral(0,hTrueDrawAlt.GetNbinsX()+1)
        hTrueDrawAlt /= trueIntAlt
        hTrueUncUpAlt /= trueIntAlt # (trueIntAlt + hTrueUncUpAlt.Integral(0,hTrueUncUpAlt.GetNbinsX()+1))
        hTrueUncDnAlt /= trueIntAlt # (trueIntAlt - hTrueUncDnAlt.Integral(0,hTrueUncDnAlt.GetNbinsX()+1))
    else:
        hTrueDrawAlt /= lumifb
        hTrueUncUpAlt /= lumifb
        hTrueUncDnAlt /= lumifb
    _normalizeBins(hTrueDrawAlt)
    _normalizeBins(hTrueUncUpAlt)
    _normalizeBins(hTrueUncDnAlt)

    errorBandTrueAlt = makeErrorBand(hTrueDrawAlt, hTrueUncUpAlt, hTrueUncDnAlt)
    errorBandTrueAlt.fillstyle = 'solid'
    errorBandTrueAlt.SetFillColorAlpha(628, 0.5)

    hTrueLegAlt = Hist(1, 0, 1)
    hTrueLegAlt.linecolor = '#ff9898'
    hTrueLegAlt.SetLineWidth(5*hTrueDrawAlt.GetLineWidth())
    hTrueLegAlt.title = hTrueDrawAlt.title
    hTrueLegAlt.legendstyle = 'L'

    trueLegOverlayAlt = Hist(1, 0, 1, linecolor=hTrueDrawAlt.linecolor,
                             legendstyle='L', title='',
                             linestyle=hTrueDrawAlt.linestyle)
    trueLegOverlayAlt.SetLineWidth(hTrueDrawAlt.GetLineWidth())

    toPlot += [errorBandTrueAlt, hTrueDrawAlt]
    forLegend.append(hTrueLegAlt)

    if varName in _matrixNames:
        with root_open(_join(_matrixPath,_matrixNames[varName]+'.root')) as fMat:
            cMat = asrootpy(fMat.canvas)
            matDist = asrootpy(cMat.FindObject(_matrixNames[varName])).clone()
            # scaleDown is for the UPPER error, scaleUp is for the LOWER
            matDistUp = asrootpy(cMat.FindObject(_matrixNames[varName]+'__scaleDown')).clone()
            matDistDn = asrootpy(cMat.FindObject(_matrixNames[varName]+'__scaleUp')).clone()
            matDist.SetDirectory(0)
            matDistUp.SetDirectory(0)
            matDistDn.SetDirectory(0)

        # un-normalize the bins, rebin, renormalize
        _unnormalizeBins(matDist)
        matDist = matDist.rebinned([e for e in hUnf._edges(0)])
        _unnormalizeBins(matDistUp)
        matDistUp = matDistUp.rebinned([e for e in hUnf._edges(0)])
        _unnormalizeBins(matDistDn)
        matDistDn = matDistDn.rebinned([e for e in hUnf._edges(0)])
        if norm:
            matDist /= _matrixXSecs['']
            matDistUp /= _matrixXSecs['up']
            matDistDn /= _matrixXSecs['dn']
        elif chan.lower() in ['eeee','mmmm']:
            matDist /= 2
            matDistUp /= 2
            matDistDn /= 2
        elif chan.lower() == 'zz':
            matDist *= 2
            matDistUp *= 2
            matDistDn *= 2

        matDist.title = r'\textbf{MATRIX}'
        matDist.color = 'forestgreen'
        matDist.drawstyle = 'hist'
        matDist.fillstyle = 'hollow'
        matDist.linestyle = 'longdash'
        matDist.SetLineWidth(matDist.GetLineWidth()*2)

        matDistUp -= matDist
        matDistDn -= matDist

        _normalizeBins(matDist)
        _normalizeBins(matDistUp)
        _normalizeBins(matDistDn)

        errorBandMat = makeErrorBand(matDist, matDistUp, matDistDn)
        errorBandMat.fillstyle = 'solid'
        errorBandMat.SetFillColorAlpha(819, 0.5)

        matDistLeg = Hist(1,0,1)
        matDistLeg.linecolor = '#98e57f'
        matDistLeg.SetLineWidth(5*matDist.GetLineWidth())
        matDistLeg.title = matDist.title
        matDistLeg.legendstyle = 'L'

        matDistLegOverlay = Hist(1, 0, 1, linecolor=matDist.linecolor,
                                 legendstyle='L', title='',
                                 linestyle=matDist.linestyle,)
        matDistLegOverlay.SetLineWidth(matDist.GetLineWidth())

        toPlot += [errorBandMat,matDist]
        forLegend.append(matDistLeg)


    if norm:
        hUncUp /= hUnfolded.Integral(0,hUnfolded.GetNbinsX()+1)
        hUncDn /= hUnfolded.Integral(0,hUnfolded.GetNbinsX()+1)
    else:
        hUncUp /= lumifb
        hUncDn /= lumifb

    _normalizeBins(hUncUp)
    _normalizeBins(hUncDn)

    if varName in _blind:
        for b, bUp, bDn in zip(hUnf, hUncUp, hUncDn):
            if hUnf.xaxis.GetBinLowEdge(b.idx) >= _blind[varName]:
                b.value = 0
                b.error = 0
                bUp.value = 0
                bUp.error = 0
                bDn.value = 0
                bDn.error = 0

    errorBand = makeErrorBand(hUnf, hUncUp, hUncDn)

    toPlot += [errorBand,hUnf]
    forLegend += [errorBand,hUnf]

    drawOpts = {
        'xtitle' : _xTitle[varName],
        }
    if norm:
        drawOpts['ytitle'] = _yTitle[varName]
    else:
        drawOpts['ytitle'] = _yTitleNoNorm[varName]
    if ana == 'full':
        drawOpts['logx'] = True
    if logy:
        drawOpts['logy'] = True
        drawOpts['yerror_in_padding'] = True
        drawOpts['ypadding'] = 0.04
    if varName == 'l1Pt':
        drawOpts['logy_crop_value'] = 4e-5

    if varName in _matrixNames:
        cUnf = Canvas(1000,1400)
        mainPad, ratioPadMat, ratioPadMain, ratioPadAlt = addPadsBelow(cUnf, 0.12, 0.12, 0.12, bottomMargin=0.35)
    else:
        cUnf = Canvas(1000,1200)
        mainPad, ratioPadMain, ratioPadAlt = addPadsBelow(cUnf, 0.15, 0.15, bottomMargin=0.35)

    mainPad.cd()
    (xaxis, yaxis), (xmin,xmax,ymin,ymax) = draw(toPlot, mainPad,
                                                 **drawOpts)
    yaxis.SetTitleSize(0.75*yaxis.GetTitleSize())
    yaxis.SetTitleOffset(1.25*yaxis.GetTitleOffset())
    #yaxis.SetLabelSize(0.82*yaxis.GetLabelSize()*1.2)
    #yaxis.SetMoreLogLabels()
    #xaxis.SetLabelSize(1.2*xaxis.GetLabelSize())#(0.82*xaxis.GetLabelSize()*1.2)
    xaxis.SetTitleSize(1.2*xaxis.GetTitleSize())
    xaxis.SetTitleOffset(1.5)#*xaxis.GetTitleOffset())#(0.95*xaxis.GetTitleOffset())

    leg = makeLegend(cUnf, *forLegend, **legParams[varName])
    leg.SetFillStyle(1001)

    if varName in _blind and _blind[varName] < xmax:
        box = TBox(max(xmin,_blind[varName]), ymin, xmax, ymax)
        box.SetFillColor(1)
        box.SetFillStyle(3002)
        box.Draw("same")
        leg.SetFillStyle(1001)

    leg.Draw("same")

    # want shaded error band and hatched fill for true legend entry,
    # so cheat and overlay a second legend
    legOverlay = asrootpy(TLegend(leg))
    legOverlay.SetFillStyle(0)
    for entry in legOverlay.primitives:
        if entry.GetLabel() == hTrueDraw.title:
            overlay = trueLegOverlay
        elif entry.GetLabel() == hTrueDrawAlt.title:
            overlay = trueLegOverlayAlt
        elif varName in _matrixNames and entry.GetLabel() == matDist.title:
            overlay = matDistLegOverlay
        else:
            entry.SetObject(None)
            entry.SetOption('')
            entry.SetLabel('')
            continue
        entry.SetObject(overlay)
        entry.SetOption(overlay.legendstyle)
        entry.SetLabel('')
    legOverlay.Draw("same")

    latex = TLatex()
    latex.SetNDC()
    latex.SetTextAlign(11)
    latex.SetTextSize(.16)#3)
    latex.SetTextFont(62)
    latexXMargin = 0.15
    if varName in ['deltaRZZ','massFull']:
        latex.SetTextAlign(31)
        latexXMargin = 1.-latexXMargin
    elif varName == 'l1Pt':
        latexXMargin = 0.35

    drawOpts = {
        'ytitle' : 'Data / Theo.',
        'xlimits' : (xmin,xmax),
        'ylimits' : (0.50001,1.9999),
        'ydivisions' : 505,
        }
    if ana == 'full':
        drawOpts['logx'] = True
    if varName in ['pt','deltaRZZ','mass']:
        drawOpts['ylimits'] = (0.2500001, 1.9999)

    if varName in _matrixNames:
        ratioPadMat.cd()

        matNoErrs = matDist.clone() # need central value only to keep ratio uncertainties consistent
        for b in matNoErrs: b.error = 0
        ratioMat, unityMat = makeRatio(hUnf, matNoErrs)

        matUnity = matDist.clone() # need errors only as baseline for ratio theory uncertainties
        for b in matUnity: b.value = 1.
        ratioTheoryErrorMat = makeErrorBand(matUnity,
                                            matDistUp/matNoErrs,
                                            matDistDn/matNoErrs)
        ratioTheoryErrorMat.fillstyle = 'solid'
        ratioTheoryErrorMat.SetFillColorAlpha(819, 0.5)

        ratioErrorMat = makeErrorBand(hUnf/matNoErrs, hUncUp/matNoErrs,
                                      hUncDn/matNoErrs)
        (ratioMatX, ratioMatY), ratioMatLimits = draw([ratioTheoryErrorMat,
                                                       ratioErrorMat,
                                                       ratioMat],
                                                      ratioPadMat,
                                                      **drawOpts)
        ratioMatY.CenterTitle()
        unityMat.Draw("same")
        latex.DrawLatex(latexXMargin, 0.8, r"\textbf{MATRIX}")

    ratioPadMain.cd()

    hTrueNoErrs = hTrueDraw.clone() # need central value only to keep ratio uncertainties consistent
    for b in hTrueNoErrs: b.error = 0

    ratioMain, unityMain = makeRatio(hUnf, hTrueNoErrs)
    ratioErrorMain = makeErrorBand(hUnf/hTrueNoErrs, hUncUp/hTrueNoErrs,
                                   hUncDn/hTrueNoErrs)

    hTrueUnity = hTrueDraw.clone() # need errors only as baseline for ratio theory uncertainties
    for b in hTrueUnity: b.value = 1.
    ratioTheoryError = makeErrorBand(hTrueUnity,
                                     hTrueUncUp/hTrueNoErrs,
                                     hTrueUncDn/hTrueNoErrs)
    ratioTheoryError.fillstyle = 'solid'
    ratioTheoryError.SetFillColorAlpha(600, 0.5)

    (ratioMainX, ratioMainY), ratioMainLimits = draw([ratioTheoryError,
                                                      ratioErrorMain,
                                                      ratioMain],
                                                     ratioPadMain,
                                                     **drawOpts)
    ratioMainY.CenterTitle()
    unityMain.Draw("same")
    latex.DrawLatex(latexXMargin, 0.8, signalName)

    ratioPadAlt.cd()

    hTrueNoErrsAlt = hTrueDrawAlt.clone() # need central value only to keep ratio uncertainties consistent
    for b in hTrueNoErrsAlt: b.error = 0

    ratioAlt, unityAlt = makeRatio(hUnf, hTrueNoErrsAlt)
    ratioErrorAlt = makeErrorBand(hUnf/hTrueNoErrsAlt, hUncUp/hTrueNoErrsAlt,
                                  hUncDn/hTrueNoErrsAlt)

    hTrueUnityAlt = hTrueDrawAlt.clone() # need errors only as baseline for ratio theory uncertainties
    for b in hTrueUnityAlt: b.value = 1.
    ratioTheoryErrorAlt = makeErrorBand(hTrueUnityAlt,
                                        hTrueUncUpAlt/hTrueNoErrsAlt,
                                        hTrueUncDnAlt/hTrueNoErrsAlt)
    ratioTheoryErrorAlt.fillstyle = 'solid'
    ratioTheoryErrorAlt.SetFillColorAlpha(628, 0.5)

    (ratioAltX, ratioAltY), ratioAltLimits = draw([ratioTheoryErrorAlt,
                                                   ratioErrorAlt,
                                                   ratioAlt],
                                                  ratioPadAlt,
                                                  **drawOpts)
    ratioAltY.CenterTitle()
    unityAlt.Draw("same")
    latex.SetTextSize(latex.GetTextSize() * ratioPadMain.height / ratioPadAlt.height)
    latex.DrawLatex(latexXMargin, 1.-.2*ratioPadMain.height/ratioPadAlt.height,
                    signalNameAlt)

    cUnf.cd()
    ratioPadAlt.Draw()
    ratioPadMain.Draw()
    if varName in _matrixNames:
        ratioPadMat.Draw()
    mainPad.Draw()

    if varName in _matrixNames:
        fixRatioAxes(xaxis,yaxis,ratioMatX,ratioMatY, mainPad.height, ratioPadMat.height)
        fixRatioAxes(ratioMatX,ratioMatY,ratioMainX,ratioMainY, ratioPadMat.height, ratioPadMain.height)
    else:
        fixRatioAxes(xaxis,yaxis,ratioMainX,ratioMainY, mainPad.height, ratioPadMain.height)
    fixRatioAxes(ratioMainX,ratioMainY,ratioAltX,ratioAltY, ratioPadMain.height, ratioPadAlt.height)

    yaxis.SetTitleSize(0.05)#0.042)
    yaxis.SetTitleOffset(1.05)

    # raster formats apparently need different fill styles?
    errorBand.fillstyle = 'x'
    errorBand.SetFillColorAlpha(1,0.6)
    ratioErrorMain.fillstyle = 'x'
    ratioErrorMain.SetFillColorAlpha(1,0.6)
    ratioErrorAlt.fillstyle = 'x'
    ratioErrorAlt.SetFillColorAlpha(1,0.6)
    if varName in _matrixNames:
        ratioErrorMat.fillstyle = 'x'
    errorBand.drawstyle = '2'
    cUnf.Update()

    suffix = ''
    if len(chan) > 2:
        suffix += '_'+chan


    _style.setCMSStyle(cUnf, '', dataType=plotType, intLumi=lumi, forLatex=True)
    cUnf.Print(_join(plotDir, 'pngs', "unfold_{}{}.png".format(varName, suffix)))
    cUnf.Print(_join(plotDir, 'Cs', "unfold_{}{}.C".format(varName, suffix)))
    _pdfViaTex(cUnf, 'unfold_{}{}'.format(varName, suffix),
               _join(plotDir, 'texs'), _join(plotDir, 'pdfs'))

    # change formatting for the postscript formats, then change it back
    #hatchWidth = gStyle.GetHatchesLineWidth()
    #hatchSpace = gStyle.GetHatchesSpacing()
    #gStyle.SetHatchesLineWidth(1)
    #gStyle.SetHatchesSpacing(1.01)
    #TAttFill.SetFillStyle(errorBand, 3244)
    #TAttFill.SetFillStyle(ratioErrorMain, 3244)
    #TAttFill.SetFillStyle(ratioErrorAlt, 3244)
    #if varName in _matrixNames:
    #    TAttFill.SetFillStyle(ratioErrorMat, 3244)
    #cUnf.Update()
    #cUnf.Print(_join(plotDir, 'epses', "unfold_{}_{}.eps".format(varName, chan)))
    #_bash('epstopdf {basedir}/epses/{filename}.eps --outfile={basedir}/pdfs/{filename}.pdf'.format(basedir=plotDir,
    #                                                                                               filename='_'.join(['unfold',varName,chan])))
    #gStyle.SetHatchesLineWidth(hatchWidth)
    #gStyle.SetHatchesSpacing(hatchSpace)


def main(inData, inMC, plotDir, fakeRateFile, puWeightFile, lumi, nIter,
         amcatnlo=False, norm=True, logy=False, looseSIP=False, noSIP=False,
         sfRemake=False, forceRedo=False, *varNames, **kwargs):
    channels = _channels[:]

    plotType = '' #'Preliminary'

    for fileType in 'Cs', 'pngs', 'texs', 'pdfs',:
        subDir = _join(plotDir, fileType)
        if not _exists(subDir):
            _mkdir(subDir)
        elif not _isdir(subDir):
            raise IOError("There is already some non-directory object called {}.".format(subDir))

    ana = 'smp'
    if 'massFull' in varNames:
        # massFull is totally different from everything else. Instead of
        # trying to make them play nicely, just do massFull totally separately
        if len(varNames) != 1:
            raise ValueError("massFull must be run separately from other variables")
        else:
            ana = 'full'

    # variables we will need if we can't just use cached histos
    samples = None
    sfFiles = None
    hPUWt = None
    hSF = None

    for varName in varNames:

        # save unfolded distributions by channel, then systematic
        hUnfolded = {}
        hTrue = {}
        hTrueAlt = {}
        hErr = {}
        hErrTrue = {}
        hErrTrueAlt = {}

        for chan in channels:
            print ""
            print "**************************************************"
            print "**** " + varName
            print "**** " + chan
            print "**************************************************"
            print ""

            # if the histograms are cached, get them (unless we don't want to)
            needToCreate = forceRedo

            if not needToCreate:
                try:
                    with root_open(_cacheFileTemplate.format(nIter), 'UPDATE') as cacheFile:
                        chanDir = getattr(cacheFile, chan)
                        varDir = getattr(chanDir, varName)

                        hUnfolded[chan] = {h.name:h.clone() for h in varDir.unfolded.objects()}
                        hUnfolded[chan][''] = varDir.hUnfolded.clone()
                        for h in hUnfolded[chan].values():
                            h.SetDirectory(0)

                        hTrue[chan] = {h.name:h.clone() for h in varDir.true.objects()}
                        hTrue[chan][''] = varDir.hTrue.clone()
                        for h in hTrue[chan].values():
                            h.SetDirectory(0)

                        hTrueAlt[chan] = {h.name:h.clone() for h in varDir.trueAlt.objects()}
                        hTrueAlt[chan][''] = varDir.hTrueAlt.clone()
                        for h in hTrueAlt[chan].values():
                            h.SetDirectory(0)
                except Exception as e:
                    rlog.warning("Rebuilding everything. Reason:")
                    rlog.warning('    '+str(e))
                    needToCreate = True

            # otherwise, we have to make everything
            # this can't be an else with the previous if
            if needToCreate:
                if samples is None:
                    sfFiles, hPUWt, hSF, sipForBkg = _generateAnalysisInputs(puWeightFile,
                                                                             looseSIP,
                                                                             noSIP,
                                                                             sfRemake)
                    samples = _generateSamples(inData, inMC, ana, fakeRateFile,
                                               puWeightFile, lumi, amcatnlo,
                                               sipForBkg, sfFiles)

                responseMakers, altResponseMakers = _generateResponseClass(
                    varName, chan, samples, hPUWt, hSF
                    )

                hUnfolded[chan], hTrue[chan], hTrueAlt[chan] = _unfold(
                    varName, chan, samples, puWeightFile, sfFiles,
                    responseMakers, altResponseMakers, nIter, plotDir
                    )

                with root_open(_cacheFileTemplate.format(nIter), 'UPDATE') as cacheFile:
                    if not hasattr(cacheFile, chan):
                        cacheFile.mkdir(chan)
                    chanDir = getattr(cacheFile, chan)
                    if hasattr(chanDir, varName):
                        chanDir.rm(varName)
                    varDir = chanDir.mkdir(varName)

                    varDir['hTrue'] = hTrue[chan][''].clone(name='hTrue')
                    if hasattr(varDir, 'true'):
                        varDir.rm('true')
                    trueDir = varDir.mkdir('true')
                    for syst, hist in hTrue[chan].iteritems():
                        if syst:
                            hist.name = syst
                            trueDir[syst] = hist.clone(name=syst)

                    varDir['hTrueAlt'] = hTrueAlt[chan][''].clone('hTrueAlt')
                    if hasattr(varDir, 'trueAlt'):
                        varDir.rm('trueAlt')
                    trueDirAlt = varDir.mkdir('trueAlt')
                    for syst, hist in hTrueAlt[chan].iteritems():
                        if syst:
                            hist.name = syst
                            trueDirAlt[syst] = hist.clone(name=syst)

                    varDir['hUnfolded'] = hUnfolded[chan][''].clone(name='hUnfolded')
                    if hasattr(varDir, 'unfolded'):
                        varDir.rm('unfolded')
                    unfDir = varDir.mkdir('unfolded')
                    for syst, hist in hUnfolded[chan].iteritems():
                        if syst:
                            hist.name = syst
                            unfDir[syst] = hist.clone(name=syst)

            hErr[chan]= _generateUncertainties(hUnfolded[chan], norm,
                                               lumi=lumi, varName=varName,
                                               ana=ana, plotDir=plotDir,
                                               chan=chan)
            (hUncUp, hUncDn) = _sumUncertainties(hErr[chan])
            hErrTrue[chan] = _generateUncertainties(hTrue[chan], norm)
            (hTrueUncUp, hTrueUncDn) = _sumUncertainties(hErrTrue[chan])
            hErrTrueAlt[chan] = _generateUncertainties(hTrueAlt[chan], norm)
            (hTrueUncUpAlt, hTrueUncDnAlt) = _sumUncertainties(hErrTrueAlt[chan])

            _generatePlots(hUnfolded[chan][''], hUncUp, hUncDn,
                           hTrue[chan][''], hTrueUncUp, hTrueUncDn,
                           hTrueAlt[chan][''], hTrueUncUpAlt, hTrueUncDnAlt,
                           varName, chan, ana, plotDir, plotType, lumi, norm, logy,
                           amcatnlo)


        hTot = sum(hUnfolded[c][''] for c in channels)
        hTrueTot = sum(hTrue[c][''] for c in channels)
        hTrueTotAlt = sum(hTrueAlt[c][''] for c in channels)

        hErrTot = _combineChannelUncertainties(*hErr.values())
        hUncUp, hUncDn = _sumUncertainties(hErrTot)
        hErrTotTrue = _combineChannelUncertainties(*hErrTrue.values())
        hTrueUncUp, hTrueUncDn = _sumUncertainties(hErrTotTrue)
        hErrTotTrueAlt = _combineChannelUncertainties(*hErrTrueAlt.values())
        hTrueUncUpAlt, hTrueUncDnAlt = _sumUncertainties(hErrTotTrueAlt)

        _generatePlots(hTot, hUncUp, hUncDn,
                       hTrueTot, hTrueUncUp, hTrueUncDn,
                       hTrueTotAlt, hTrueUncUpAlt, hTrueUncDnAlt,
                       varName, 'zz', ana, plotDir, plotType, lumi, norm, logy,
                       amcatnlo)



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
    parser.add_argument('--nIter', type=int, nargs='?', default=4,
                        help='Number of iterations for D\'Agostini method')
    parser.add_argument('--amcatnlo', action='store_true',
                        help='Use MadGraph5_aMC@NLO as the primary MC and '
                        'Powheg as the cross-check, instead of the other way '
                        'around.')
    parser.add_argument('--analysis', '--ana', type=str, nargs='?',
                        default='smp',
                        help='Which set of cuts to use (full, smp, etc.).')
    parser.add_argument('--variables', type=str, nargs='*',
                        default=_varListNoFull,
                        help=('Names of variables to use. Options are: {}. '
                              'If not specified, all are used except massFull'
                              ).format(', '.join(_varList)))
    parser.add_argument('--noNorm', action='store_true',
                        help='Leave differential cross sections in abolute normalization rather than normalizing to unity area.')
    parser.add_argument('--logy', '--logY', '--log', action='store_true',
                        help='Put vertical axis on a log scale.')
    parser.add_argument('--looseSIP', action='store_true',
                        help='Use scale factors for SIP<10 with no extra IP cuts.')
    parser.add_argument('--noSIP', action='store_true',
                        help='Use scale factors for no SIP cut and no extra IP cuts.')
    parser.add_argument('--sfRemake', action='store_true',
                        help='Use homebrewed scale factors for electrons.')
    parser.add_argument('--redo', action='store_true',
                        help='Make new histograms even if some are cached.')

    args=parser.parse_args()

    if not _exists(args.plotDir):
        _mkdir(args.plotDir)
    elif not _isdir(args.plotDir):
        raise IOError("There is already some non-directory object called {}.".format(args.plotDir))

    main(args.dataDir, args.mcDir, args.plotDir, args.fakeRateFile,
         args.puWeightFile, args.lumi, args.nIter, args.amcatnlo,
         not args.noNorm, args.logy, args.looseSIP, args.noSIP, args.sfRemake,
         args.redo, *args.variables)

