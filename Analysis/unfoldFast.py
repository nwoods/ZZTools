import logging
from rootpy import log as rlog; rlog = rlog["/unfoldFast"]
# don't show most silly ROOT messages
logging.basicConfig(level=logging.WARNING)
rlog["/ROOT.TUnixSystem.SetDisplay"].setLevel(rlog.ERROR)
rlog["/ROOT.TROOT.Append"].setLevel(rlog.ERROR)

from rootpy import asrootpy
from rootpy.io import root_open
from rootpy.plotting import Canvas, Legend, Hist, HistStack, Graph
from rootpy.plotting.utils import draw
from rootpy.ROOT import cout, TDecompSVD, TBox, TLatex
from rootpy.ROOT import gSystem, TString, gCling
from rootpy.ROOT import RooUnfoldResponse as Response
from rootpy.ROOT import RooUnfoldBayes as RooUnfoldIter # it's frequentist!
from rootpy.context import preserve_current_directory
import rootpy.compiled as _rootComp
from rootpy.stl import vector as _Vec
_VFloat = _Vec('float')

from PlotTools import PlotStyle as _Style
from PlotTools import makeLegend, addPadsBelow, makeRatio, fixRatioAxes, makeErrorBand
from Utilities import WeightStringMaker, Z_MASS, deltaRString, deltaPhiString, zeroNegativeBins
from Analysis.setupStandardSamples import *
# from Analysis.unfoldingHelpers import getResponse, getResponsePDFErrors, \
#     getResponseScaleErrors, getResponseAlphaSErrors
from Analysis.weightHelpers import puWeight, baseMCWeight
from Metadata.metadata import sampleInfo

from os import environ as _env
from os import makedirs as _mkdir
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
#gCling.AddIncludePath(_join(_zztBaseDir, 'RooUnfold-1.1.1'))
#gSystem.Load(_join(_zztBaseDir, 'RooUnfold-1.1.1', 'libRooUnfold'))


_channels = ['eeee','eemm', 'mmmm']

# set up variables, selection, binnings etc.
# (jet-related variables and selections done later)
_variables = {
    'pt' : {c:'Pt' for c in _channels},
    'mass' : {c:'Mass' for c in _channels},
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
    'nJets' : [6,-0.5,5.5],
    'mass' : [100.] + [200.+50.*i for i in range(5)] + [500.,600.,800.],
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
    'l1Pt' : [15, 0., 150.],
    }

_binningFull = _binning.copy()
_binningFull['mass'] = [80.,100.,120.,130.,150.,180.,200.,240.,300.,400.,1000]
_binningFull['z1Mass'] = [40,40.,120.]
_binningFull['z1Mass'] = [60,0.,120.]

_units = {
    'pt' : 'GeV',
    'nJets' : '',
    'mass' : 'GeV',
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
    'pt' : 'p_T^{4\\ell}',
    'nJets' : '# jets',
    'mass' : 'm_{4\\ell}',
    'eta' : '\\eta_{4\\ell}',
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

# Names of compiled C++ classes to make response matrices fast
# (this is extremely slow in Python because it requires a combination of
# information from multiple trees, which can't be done with TTree::Draw())
_responseClassNames = {
    'mass' : {c:'FloatBranchResponseMatrixMaker' for c in _channels},
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
    'deltaPhiZZ' : {c:'ZZDeltaPhiResponseMatrixMaker' for c in _channels},
    'deltaRZZ' : {c:'ZZDeltaRResponseMatrixMaker' for c in _channels},
    'lPt' : {c:'AllLeptonBranchResponseMatrixMaker' for c in _channels},
    'l1Pt' : {c:'LeptonMaxBranchResponseMatrixMaker' for c in _channels},
    'zPt' : {c:'BothZsBranchResponseMatrixMaker' for c in _channels},
    }

# Variable names usable by response maker classes
_varNamesForResponseMaker = {
    'mass' : {c:'Mass' for c in _channels},
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
_legParams['eta'] = _legParams['deltaEtajj'].copy()
_legParams['lPt']['topmargin'] = 0.05
_legParams['l1Pt']['topmargin'] = 0.05


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


def _normalizeBins(h):
    binUnit = 1 # min(h.GetBinWidth(b) for b in range(1,len(h)+1))
    for ib in xrange(1,len(h)+1):
        w = h.GetBinWidth(ib)
        h.SetBinContent(ib, h.GetBinContent(ib) * binUnit / w)
        h.SetBinError(ib, h.GetBinError(ib) * binUnit / w)
        if h.GetBinError(ib) > h.GetBinContent(ib):
            h.SetBinError(ib, h.GetBinContent(ib))
    h.sumw2()

_printNext = False

def _getUnfolded(hSig, hBkg, hTrue, hResponse, hData, nIter,
                 withRespAndCov=False, printIt=False):
    global _printNext

    hSigPlusBkg = hSig+hBkg
    zeroNegativeBins(hSigPlusBkg)

    response = Response(hSigPlusBkg, hTrue.clone(), hResponse.clone())

    try:
        svd = TDecompSVD(response.Mresponse())
        sig = svd.GetSig()
        try:
            condition = sig.Max() / max(0., sig.Min())
        except ZeroDivisionError:
            condition = float('inf')
            print hSig.Integral(), hBkg.Integral(), hTrue.Integral(), hResponse.Integral()

        print ''
        print 'condition: {}'.format(condition)
        print ''

    except:
        print "It broke!"
        print hSig.Integral(), hBkg.Integral(), hTrue.Integral(), hResponse.Integral()
        c = Canvas(1000,1000)
        hSig.draw()
        c.Print("sig.png")
        hBkg.draw()
        c.Print("bkg.png")
        hTrue.draw()
        c.Print("true.png")
        hData.draw()
        c.Print("data.png")
        hResponse.drawstyle = 'colz'
        hResponse.draw()
        c.Print("resp.png")


    hDataMinusBkg = hData - hBkg
    zeroNegativeBins(hDataMinusBkg)

    unf = RooUnfoldIter(response, hDataMinusBkg, nIter)

    if _printNext or printIt:
        _printNext = False

        c = Canvas(1000,1000)
        hSig.draw()
        c.Print("sig.png")
        hBkg.draw()
        c.Print("bkg.png")
        hTrue.draw()
        c.Print("true.png")
        hData.draw()
        c.Print("data.png")
        hResponse.drawstyle = 'colz'
        hResponse.draw()
        c.Print("resp.png")

    hOut = unf.Hreco()
    if not hOut:
        print hOut
        raise ValueError("The unfolded histogram got screwed up somehow!")

    if withRespAndCov:
        return asrootpy(hOut), unf.Ereco(2).Clone(), asrootpy(response.Hresponse()).clone()

    return asrootpy(hOut)



def main(inData, inMC, plotDir, fakeRateFile, puWeightFile, lumi, nIter,
         amcatnlo=False, ana='smp', *varNames, **kwargs):

    classesNeeded = list(set(cls  for v in varNames for cls in _responseClassNames[v].values()))
    _rootComp.register_file(_join(_zztBaseDir, 'Utilities',
                                  'ResponseMatrixMaker.cxx'),
                            classesNeeded)

    # force compilation
    _C = getattr(_rootComp, classesNeeded[0])

    style = _Style()

    channels = _channels

    puWeightStr, puWt = puWeight(puWeightFile, '')
    puWeightStrUp, puWtUp = puWeight(puWeightFile, 'up')
    puWeightStrDn, puWtDn = puWeight(puWeightFile, 'dn')

    puWeightFileFull = _join(_env['zzt'],'data','pileup',puWeightFile+'.root')
    with preserve_current_directory():
        with root_open(puWeightFileFull) as fPU:
            hPUWt = asrootpy(fPU.puScaleFactor)
            hPUWtUp = asrootpy(fPU.puScaleFactor_up)
            hPUWtDn = asrootpy(fPU.puScaleFactor_down)
            hPUWt.SetDirectory(0)
            hPUWtUp.SetDirectory(0)
            hPUWtDn.SetDirectory(0)

    true = genZZSamples('zz', inMC, ana, lumi, amcatnlo=amcatnlo,
                        higgs=(ana=='full'))
    reco = zzStackSignalOnly('zz', inMC, ana, puWeightFile,
                             lumi, amcatnlo=amcatnlo, asGroup=True,
                             higgs=(ana=='full'))
    sigFileNames = {s.name : [f for f in s.getFileNames()]
                    for s in reco.values()[0].getBaseSamples()}
    sigConstWeights = {s.name : s.xsec * s.intLumi * float(s.kFactor) / s.sumW
                       for s in reco.values()[0].getBaseSamples()}

    bkgMC = zzIrreducibleBkg('zz', inMC, ana, puWeightFile, lumi)
    bkg = standardZZBkg('zz', inData, inMC, ana, puWeightFile,
                        fakeRateFile, lumi)
    bkgSyst = {
        'eup' : standardZZBkg('zz', inData, inMC, ana, puWeightFile,
                              fakeRateFile, lumi, eFakeRateSyst='up'),
        'edn' : standardZZBkg('zz', inData, inMC, ana, puWeightFile,
                              fakeRateFile, lumi, eFakeRateSyst='dn'),
        'mup' : standardZZBkg('zz', inData, inMC, ana, puWeightFile,
                              fakeRateFile, lumi, mFakeRateSyst='up'),
        'mdn' : standardZZBkg('zz', inData, inMC, ana, puWeightFile,
                              fakeRateFile, lumi, mFakeRateSyst='dn'),
        }

    data = standardZZData('zz', inData, ana)

    altReco = zzStackSignalOnly('zz', inMC, ana, puWeightFile, lumi,
                                amcatnlo=(not amcatnlo), asGroup=True,
                                higgs=(ana=='full'))
    altSigFileNames = {s.name : [f for f in s.getFileNames()]
                       for s in altReco.values()[0].getBaseSamples()}
    altSigConstWeights = {s.name : s.xsec * s.intLumi * float(s.kFactor) / s.sumW
                          for s in altReco.values()[0].getBaseSamples()}

    altTrue = genZZSamples('zz', inMC, ana, lumi,
                           amcatnlo=(not amcatnlo), higgs=(ana=='full'))

    if amcatnlo:
        signalName = 'MG5_aMC@NLO'
        signalNameAlt = 'POWHEG'
        if ana=='full':
            signalName += '+POWHEG'
    else:
        signalName = 'POWHEG'
        signalNameAlt = 'MG5_aMC@NLO'
        if ana=='full':
            signalNameAlt += '+POWHEG'
    signalName += '+MCFM'
    signalNameAlt += '+MCFM'

    recoSyst = {}
    bkgMCSyst = {}
    sigFileNamesSyst = {}
    for syst in ['eScaleUp', 'eScaleDn', 'eRhoResUp',
                 'eRhoResDn', 'ePhiResUp']:
        recoSyst[syst] = zzStackSignalOnly('eeee,eemm', inMC.replace('mc_','mc_{}_'.format(syst)),
                                           ana, puWeightFile, lumi, amcatnlo=amcatnlo,
                                           asGroup=True, higgs=(ana=='full'))
        sigFileNamesSyst[syst] = {s.name : [f for f in s.getFileNames()]
                        for s in recoSyst[syst].values()[0].getBaseSamples()}
        bkgMCSyst[syst] = zzIrreducibleBkg('eeee,eemm', inMC.replace('mc_','mc_{}_'.format(syst)),
                                           ana, puWeightFile, lumi)

    for syst in ['mClosureUp','mClosureDn']:
        recoSyst[syst] = zzStackSignalOnly('eemm,mmmm', inMC.replace('mc_','mc_{}_'.format(syst)),
                                           ana, puWeightFile, lumi, amcatnlo=amcatnlo,
                                           asGroup=True, higgs=(ana=='full'))
        sigFileNamesSyst[syst] = {s.name : [f for f in s.getFileNames()]
                                  for s in recoSyst[syst].values()[0].getBaseSamples()}
        bkgMCSyst[syst] = zzIrreducibleBkg('eemm,mmmm', inMC.replace('mc_','mc_{}_'.format(syst)),
                                           ana, puWeightFile, lumi)


    for varName in varNames:

        if ana == 'full':
            binning = _binningFull[varName]
        else:
            binning = _binning[varName]
        vBinning = _VFloat()
        if len(binning) == 3:
            binningTemp = [binning[1] + i * (binning[2] - binning[1])/float(binning[0]) for i in xrange(binning[0]+1)]
            for b in binningTemp:
                vBinning.push_back(b)
        else:
            for b in binning:
                vBinning.push_back(b)

        # save unfolded distributions by channel, then systematic
        hErr = {}
        hUnfolded = {c:{} for c in channels}

        for chan in channels:
            print ""
            print "**************************************************"
            print "**** " + varName
            print "**** " + chan
            print "**************************************************"
            print ""

            var = _variables[varName][chan]
            sel = _selections[varName][chan]
            selTrue = sel #combineWeights(sel, _wrongZRejectionStr[chan], selections=True)

            ResponseMakerClass = getattr(_rootComp, _responseClassNames[varName][chan])

            responseMakers = {}
            for sample, fNameList in sigFileNames.iteritems():
                resp = ResponseMakerClass(chan,
                                          _varNamesForResponseMaker[varName][chan],
                                          vBinning)

                for fName in fNameList:
                    resp.registerFile(fName)
                for syst in sigFileNamesSyst:
                    for fName in sigFileNamesSyst[syst][sample]:
                        resp.registerFile(fName, syst)

                resp.registerPUWeights(hPUWt)
                resp.registerPUWeights(hPUWtUp, 'up')
                resp.registerPUWeights(hPUWtDn, 'dn')
                resp.setConstantScale(sigConstWeights[sample])
                responseMakers[sample] = resp

            altResponseMakers = {}
            for sample, fNameList in altSigFileNames.iteritems():
                if sample in responseMakers:
                    continue
                resp = ResponseMakerClass(chan,
                                          _varNamesForResponseMaker[varName][chan],
                                          vBinning)
                for fName in fNameList:
                    resp.registerFile(fName)
                resp.registerPUWeights(hPUWt)
                resp.setConstantScale(altSigConstWeights[sample])
                resp.setSkipSystematics()
                altResponseMakers[sample] = resp

            ### Do all the unfolding

            hData = data[chan].makeHist(var, sel, binning, perUnitWidth=False)

            # regular weight, no systematics. Apply just in case.
            nominalWeight = baseMCWeight(chan, puWeightFile)
            reco[chan].applyWeight(nominalWeight, True)
            altReco[chan].applyWeight(nominalWeight, True)
            bkgMC[chan].applyWeight(nominalWeight, True)
            for s in recoSyst.values():
                try:
                    s[chan].applyWeight(nominalWeight, True)
                except KeyError:
                    pass

            hTrueNominal = true[chan].makeHist(var, selTrue, binning,
                                               perUnitWidth=False)
            hSigNominal = reco[chan].makeHist(var, sel, binning, perUnitWidth=False)
            hBkgMCNominal = bkgMC[chan].makeHist(var, sel, binning, perUnitWidth=False)
            hBkgNominal = bkg[chan].makeHist(var, sel, binning, perUnitWidth=False)
            hResponseNominal = {s:asrootpy(resp()) for s,resp in responseMakers.iteritems()}
            hResponseNominalTotal = sum(resp for resp in hResponseNominal.values())
            hDataNoBkgNominal = hData - hBkgMCNominal - hBkgNominal

            # also get covariance and response matrices for later plotting
            hUnfolded[chan][''], hCov, hResp = _getUnfolded(hSigNominal,
                                                            hBkgMCNominal+hBkgNominal,
                                                            hTrueNominal,
                                                            hResponseNominalTotal,
                                                            hData, nIter, True)


            # PU reweight uncertainty
            for sys in ['up','dn']:
                wtStr = baseMCWeight(chan, puWeightFile, puSyst=sys)
                reco[chan].applyWeight(wtStr, True)
                bkgMC[chan].applyWeight(wtStr, True)

                hSig = reco[chan].makeHist(var, sel, binning, perUnitWidth=False)
                hBkgMC = bkgMC[chan].makeHist(var, sel, binning, perUnitWidth=False)

                hResponse = sum(asrootpy(resp('pu_'+sys)) for resp in responseMakers.values())

                hUnfolded[chan]['pu_'+sys] = _getUnfolded(hSig,
                                                          hBkgMC+hBkgNominal,
                                                          hTrueNominal,
                                                          hResponse,
                                                          hData, nIter)

                reco[chan].applyWeight(nominalWeight, True)
                bkgMC[chan].applyWeight(nominalWeight, True)


            # lepton efficiency uncertainty
            for lep in set(chan):
                for sys in ['up','dn']:
                    wtArg = {lep+'Syst':sys}
                    wtStr = baseMCWeight(chan, puWeightFile, **wtArg)
                    reco[chan].applyWeight(wtStr, True)
                    bkgMC[chan].applyWeight(wtStr, True)

                    hSig = reco[chan].makeHist(var, sel, binning, perUnitWidth=False)
                    hBkgMC = bkgMC[chan].makeHist(var, sel, binning, perUnitWidth=False)

                    hResponse = sum(asrootpy(resp(lep+'Eff_'+sys)) for resp in responseMakers.values())

                    hUnfolded[chan][lep+'Eff_'+sys] = _getUnfolded(hSig,
                                                                   hBkgMC+hBkgNominal,
                                                                   hTrueNominal,
                                                                   hResponse,
                                                                   hData, nIter)

                    reco[chan].applyWeight(nominalWeight, True)
                    bkgMC[chan].applyWeight(nominalWeight, True)

            # alternate generator
            hSig = altReco[chan].makeHist(var, sel, binning,
                                          perUnitWidth=False)
            hTrue = altTrue[chan].makeHist(var, selTrue, binning,
                                           perUnitWidth=False)
            hResponses = []
            for s in altSigFileNames.keys():
                try:
                    hResponses.append(asrootpy(altResponseMakers[s]()))
                except KeyError:
                    hResponses.append(hResponseNominal[s])
            hResponse = sum(h for h in hResponses)

            hUnfolded[chan]['generator'] = _getUnfolded(hSig,
                                                        hBkgMCNominal+hBkgNominal,
                                                        hTrue,
                                                        hResponse,
                                                        hData, nIter)

            # luminosity
            lumiUnc = 0.062
            lumiScale = {'up':1.+lumiUnc,'dn':1.-lumiUnc}
            for sys, scale in lumiScale.iteritems():
                hSig = hSigNominal * scale
                hBkgMC = hBkgMCNominal * scale
                hTrue = hTrueNominal * scale
                hResponse = hResponseNominalTotal * scale

                hUnfolded[chan]['lumi_'+sys] = _getUnfolded(hSig,
                                                            hBkgMC+hBkgNominal,
                                                            hTrue,
                                                            hResponse,
                                                            hData, nIter)

            # lepton fake rate uncertainty
            for lep in set(chan):
                for sys in ['up','dn']:
                    hBkg = bkgSyst[lep+sys][chan].makeHist(var, sel, binning, perUnitWidth=False)

                    hUnfolded[chan][lep+'FR_'+sys] = _getUnfolded(hSigNominal,
                                                                  hBkgMCNominal+hBkg,
                                                                  hTrueNominal,
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

                        hSig = reco[chan].makeHist(varShifted, selShifted,
                                                   binning, perUnitWidth=False)
                        hBkgMC = bkgMC[chan].makeHist(varShifted, selShifted,
                                                     binning, perUnitWidth=False)

                        hResponse = sum(asrootpy(resp(sys+'_'+shift)) for resp in responseMakers.values())

                        hUnfolded[chan][sys+'_'+shift] = _getUnfolded(hSig,
                                                                      hBkgMC+hBkgNominal,
                                                                      hTrueNominal,
                                                                      hResponse,
                                                                      hData, nIter)

            # lepton momentum uncertainties
            if 'e' in chan:
                for sys in ['eScale', 'eRhoRes', 'ePhiRes']:
                    for shift in ['up','dn']:
                        if sys == 'ePhiRes' and shift == 'dn':
                            continue
                        sysStr = 'Up' if shift == 'up' else 'Dn'

                        hSig = recoSyst[sys+sysStr][chan].makeHist(var, sel,
                                                                   binning,
                                                                   perUnitWidth=False)
                        hBkgMC = bkgMCSyst[sys+sysStr][chan].makeHist(var, sel,
                                                                      binning,
                                                                      perUnitWidth=False)
                        hResponse = sum(asrootpy(resp(sys+'_'+shift)) for resp in responseMakers.values())

                        storeAs = sys+'_'+shift
                        if sys == 'ePhiRes':
                            storeAs = sys
                        hUnfolded[chan][storeAs] = _getUnfolded(hSig,
                                                                      hBkgMC+hBkgNominal,
                                                                      hTrueNominal,
                                                                      hResponse,
                                                                      hData, nIter)
            if 'm' in chan:
                sys = 'mClosure'
                for shift in ['up','dn']:
                    sysStr = 'Up' if shift == 'up' else 'Dn'

                    hSig = recoSyst[sys+sysStr][chan].makeHist(var, sel,
                                                               binning,
                                                               perUnitWidth=False)
                    hBkgMC = bkgMCSyst[sys+sysStr][chan].makeHist(var, sel,
                                                                  binning,
                                                                  perUnitWidth=False)
                    hResponse = sum(asrootpy(resp(sys+'_'+shift)) for resp in responseMakers.values())

                    hUnfolded[chan][sys+'_'+shift] = _getUnfolded(hSig,
                                                                  hBkgMC+hBkgNominal,
                                                                  hTrueNominal,
                                                                  hResponse,
                                                                  hData, nIter)

            # PDF uncertainties
            hSigVariations = []
            for s in reco[chan].values():
                if 'GluGluZZ' not in s.name:
                    hSigVariations.append(s.makeHist2(var, 'Iteration$', sel, binning,
                                                      [100,0.,100.], 'pdfWeights', False))
            hResponseVariations = []
            for s, resp in responseMakers.iteritems():
                if "GluGluZZ" not in s:
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
                hSigDn[i+1].value = max(0.,hSigDn[i+1].value + sigBinRMSes[i])
            for x in xrange(hResponseUp.GetNbinsX()):
                for y in xrange(hResponseUp.GetNbinsY()):
                    hResponseUp[x+1,y+1].value += responseBinRMSes[x][y]
                    hResponseDn[x+1,y+1].value = max(0., hResponseDn[x+1,y+1].value - responseBinRMSes[x][y])

            hTrueUp = hTrueNominal.clone()
            hTrueDn = hTrueNominal.clone()
            hTrueVariations = []
            for s in true[chan].values():
                if 'GluGluZZ' not in s.name:
                    hTrueVariations.append(s.makeHist2(var, 'Iteration$', sel, binning,
                                                       [100,0.,100.], 'pdfWeights', False))
            allTrueRMSes = [[Graph(h.ProjectionY('slice{}'.format(i), i+1,i+1)).GetRMS(2) for i in xrange(h.GetNbinsX())] for h in hTrueVariations]
            binTrueRMSes = [sum(rmses) for rmses in zip(*allTrueRMSes)]

            for i in xrange(hTrueUp.GetNbinsX()):
                hTrueUp[i+1].value += binTrueRMSes[i]
                hTrueDn[i+1].value = max(0.,hTrueDn[i+1].value + binTrueRMSes[i])

            hUnfolded[chan]['pdf_up'] = _getUnfolded(hSigUp,
                                                     hBkgMCNominal+hBkgNominal,
                                                     hTrueUp, hResponseUp,
                                                     hData, nIter)
            hUnfolded[chan]['pdf_dn'] = _getUnfolded(hSigDn,
                                                     hBkgMCNominal+hBkgNominal,
                                                     hTrueDn, hResponseDn,
                                                     hData, nIter)

            # QCD scale uncertainties
            variationIndices = [1,2,3,4,6,8]
            hSigs = [reco[chan].makeHist(var, sel, binning,
                                 {
                        'ZZTo4L':'scaleWeights[{}]'.format(i),
                        'ZZTo4L-amcatnlo':'scaleWeights[{}]'.format(i),
                        'ZZJJTo4L_EWK':'scaleWeights[{}]'.format(i),
                        },
                                 perUnitWidth=False)
                     for i in variationIndices]

            hTrues = [true[chan].makeHist(var, selTrue, binning,
                                          {
                        'ZZTo4L':'scaleWeights[{}]'.format(i),
                        'ZZTo4L-amcatnlo':'scaleWeights[{}]'.format(i),
                        'ZZJJTo4L_EWK':'scaleWeights[{}]'.format(i),
                        },
                                          perUnitWidth=False)
                      for i in variationIndices]

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
            for hSig, hTrue, hResponse in zip(hSigs, hTrues, hResponseVariations):
                hUnfoldedVariations.append(_getUnfolded(hSig,
                                                        hBkgMCNominal+hBkgNominal,
                                                        hTrue, hResponse,
                                                        hData, nIter))

            hUnfoldedUp = hUnfoldedVariations[0].empty_clone()
            hUnfoldedDn = hUnfoldedVariations[0].empty_clone()
            for bUp, bDn, bVars in zip(hUnfoldedUp, hUnfoldedDn, zip(*hUnfoldedVariations)):
                bUp.value = max(b.value for b in bVars)
                bDn.value = min(b.value for b in bVars)

            hUnfolded[chan]['scale_up'] = hUnfoldedUp
            hUnfolded[chan]['scale_dn'] = hUnfoldedDn

            # alpha_s uncertainties
            alphaSIndices = [100,101]
            hSigs = [reco[chan].makeHist(var, sel, binning,
                                 {
                        'ZZTo4L':'pdfWeights[{}]'.format(i),
                        'ZZTo4L-amcatnlo':'pdfWeights[{}]'.format(i),
                        'ZZJJTo4L_EWK':'pdfWeights[{}]'.format(i),
                        },
                                 perUnitWidth=False)
                    for i in alphaSIndices]
            hTrues = [reco[chan].makeHist(var, sel, binning,
                                  {
                        'ZZTo4L':'pdfWeights[{}]'.format(i),
                        'ZZTo4L-amcatnlo':'pdfWeights[{}]'.format(i),
                        'ZZJJTo4L_EWK':'pdfWeights[{}]'.format(i),
                        },
                                  perUnitWidth=False)
                     for i in alphaSIndices]

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

            hUnfolded[chan]['alphaS_up'] = hUnfolded[chan][''] + unc
            hUnfolded[chan]['alphaS_dn'] = hUnfolded[chan][''] - unc

            # since MCFM samples don't have LHE information, we just vary by
            # the cross section uncertainties
            mcfmUnc = {'up':.18,'dn':-.15}
            for sys, shift in mcfmUnc.iteritems():
                hSig = reco[chan].makeHist(var, sel, binning,
                                           {'GluGluZZ':str(1.+shift)},
                                           perUnitWidth=False)
                hTrue = true[chan].makeHist(var, selTrue, binning,
                                           {'GluGluZZ':str(1.+shift)},
                                           perUnitWidth=False)
                hResponse = hResponseNominalTotal.empty_clone()
                for s, h in hResponseNominal.iteritems():
                    if 'GluGluZZ' in s:
                        hResponse += h * (1.+shift)
                    else:
                        hResponse += h

                hUnfolded[chan]['mcfmxsec_'+sys] = _getUnfolded(hSig,
                                                                hBkgNominal+hBkgMCNominal,
                                                                hTrue,
                                                                hResponse,
                                                                hData, nIter)

            # Normalize all the things
            for sys, hUnf in hUnfolded[chan].iteritems():
                hUnf /= hUnf.Integral(0,hUnf.GetNbinsX()+1)

            # Make uncertainties out of the unfolded histos
            hErr[chan] = {'up':{},'dn':{}}
            for sys, hUnf in hUnfolded[chan].iteritems():
                if not sys:
                    continue

                he = hUnf - hUnfolded[chan]['']
                sysName = sys.replace('_up','').replace('_dn','')
                he.title = _uncertaintyTitles[sysName]
                he.color = _uncertaintyColors[sysName]
                he.fillstyle = 'solid'
                he.drawstyle = 'hist'
                he.legendstyle = 'F'

                if '_up' in sys:
                    hErr[chan]['up'][sysName] = he
                elif '_dn' in sys:
                    hErr[chan]['dn'][sysName] = he
                else:
                    hErr[chan]['up'][sysName] = he
                    he2 = he.clone()
                    he2.title = _uncertaintyTitles[sysName]
                    he2.color = _uncertaintyColors[sysName]
                    he2.fillstyle = 'solid'
                    he2.drawstyle = 'hist'
                    he2.legendstyle = 'F'
                    hErr[chan]['dn'][sysName] = he2



            ### Get the total uncertainties
            # this method of assigning up vs down should be revisited
            hUncUp = hUnfolded[chan][''].empty_clone()
            hUncDn = hUnfolded[chan][''].empty_clone()
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
                    h /= hUnfolded[chan]['']
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
                    h *= hUnfolded[chan]['']
                    h /= 100.

            ### plot
            hUnf = hUnfolded[chan][''].clone()
            hUnf.color = 'black'
            hUnf.drawstyle = 'PE'
            hUnf.legendstyle = 'LPE'
            hUnf.title = 'Data + stat. unc.'
            _normalizeBins(hUnf)

            hTrue = hTrueNominal.clone()
            hTrue.fillcolor = '#99ccff'
            hTrue.linecolor = '#000099'
            hTrue.drawstyle = 'hist'
            hTrue.fillstyle = 'solid'
            hTrue.legendstyle = 'F'
            hTrue.title = '{}'.format(signalName)
            hTrue /= hTrue.Integral(0,hTrue.GetNbinsX()+1)
            _normalizeBins(hTrue)

            hTrueAlt = altTrue[chan].makeHist(var, selTrue, binning,
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
                for b, bUp, bDn in zip(hUnf, hUncUp, hUncDn):
                    if hUnf.xaxis.GetBinLowEdge(b.idx) >= _blind[varName]:
                        b.value = 0
                        b.error = 0
                        bUp.value = 0
                        bUp.error = 0
                        bDn.value = 0
                        bDn.error = 0

            errorBand = makeErrorBand(hUnf, hUncUp, hUncDn)

            cUnf = Canvas(1000,1200)
            mainPad, ratioPad1, ratioPad2 = addPadsBelow(cUnf, 0.15, 0.15, bottomMargin=0.35)

            mainPad.cd()
            (xaxis, yaxis), (xmin,xmax,ymin,ymax) = draw([hTrue, hTrueAlt,
                                                          #hUnfoldedAlt,
                                                          hUnf,
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
                             hUnf, **_legParams[varName])

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
            ratio1, unity1 = makeRatio(hUnf, hTrue)
            ratioError1 = makeErrorBand(hUnf/hTrue, hUncUp/hTrue,
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
            ratio2, unity2 = makeRatio(hUnf, hTrueAlt)
            ratioError2 = makeErrorBand(hUnf/hTrueAlt, hUncUp/hTrueAlt,
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
            hResp.drawstyle = 'colztext'
            hResp.xaxis.title = '\\text{Reco} '+_xTitle[varName]
            hResp.yaxis.title = '\\text{True} '+_xTitle[varName]
            hResp.draw()
            style.setCMSStyle(cRes, '', dataType='Preliminary Simulation', intLumi=lumi)
            cRes.Print(_join(plotDir, "response_{}_{}.png".format(varName, chan)))
            cRes.Print(_join(plotDir, "response_{}_{}.C".format(varName, chan)))

            cCov = Canvas(1000,1000)
            hCov.Draw("colztext")
            style.setCMSStyle(cCov, '', dataType='Preliminary', intLumi=lumi)
            cCov.Print(_join(plotDir, "covariance_{}_{}.png".format(varName, chan)))
            cCov.Print(_join(plotDir, "covariance_{}_{}.C".format(varName, chan)))

        hTot = sum(hUnfolded[c][''] for c in channels)
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

        hTrue = true.makeHist(_variables[varName], _selections[varName],
                              binning, perUnitWidth=False)
        hTrue.fillcolor = '#99ccff'
        hTrue.linecolor = '#000099'
        hTrue.drawstyle = 'hist'
        hTrue.fillstyle = 'solid'
        hTrue.legendstyle = 'F'
        hTrue.title = '{}'.format(signalName)
        hTrue /= hTrue.Integral(0,hTrue.GetNbinsX()+1)
        _normalizeBins(hTrue)

        hTrueAlt = altTrue.makeHist(_variables[varName], _selections[varName],
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
                    hThis *= hUnfolded[chan][''] / hTotNoNorm
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
                if hTot.xaxis.GetBinLowEdge(b.idx) >= _blind[varName]:
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
    parser.add_argument('--amcatnlo', action='store_true',
                        help='Use MadGraph5_aMC@NLO as the primary MC and '
                        'Powheg as the cross-check, instead of the other way '
                        'around.')
    parser.add_argument('--analysis', '--ana', type=str, nargs='?',
                        default='smp',
                        help='Which set of cuts to use (full, smp, etc.).')
    parser.add_argument('--variables', type=str, nargs='*',
                        default=_varList,
                        help=('Names of variables to use. If not specified, '
                              'all are used ({})').format(', '.join(_varList)))

    args=parser.parse_args()

    if not _exists(args.plotDir):
        _mkdir(args.plotDir)
    elif not _isdir(args.plotDir):
        raise IOError("There is already some non-directory object called {}.".format(args.plotDir))

    main(args.dataDir, args.mcDir, args.plotDir, args.fakeRateFile,
         args.puWeightFile, args.lumi, args.nIter, args.amcatnlo,
         args.analysis, *args.variables)

