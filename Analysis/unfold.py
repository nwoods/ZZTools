
import logging
from rootpy import log as rlog; rlog = rlog["/unfold"]
# don't show most silly ROOT messages
logging.basicConfig(level=logging.WARNING)
rlog["/ROOT.TUnixSystem.SetDisplay"].setLevel(rlog.ERROR)

from rootpy import asrootpy
from rootpy.io import root_open
from rootpy.plotting import Canvas, Legend, Hist, Hist2D, HistStack
from rootpy.plotting.utils import draw
from rootpy.ROOT import cout, TDecompSVD
from rootpy.ROOT import RooUnfoldBayes as RooUnfoldIter # it's frequentist!

from SampleTools import MCSample, DataSample, SampleGroup, SampleStack
from PlotTools import PlotStyle as _Style
from PlotTools import makeLegend, addPadBelow, makeRatio, fixRatioAxes, makeErrorBand
from Utilities import WeightStringMaker, identityFunction, Z_MASS, \
    deltaRFunction, deltaRString, deltaPhiFunction, deltaPhiString
from Analysis.setupStandardSamples import *
from Analysis.unfoldingHelpers import getResponse
from Analysis.weightHelpers import puWeight, baseMCWeight
from Metadata.metadata import sampleInfo

from os import environ as _env
from os.path import join as _join
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

def _makeComparator(var, compareTo):
    '''
    Make a little function to check whether var is greater than compareTo
    '''
    # insurance against stupid scoping issues
    newVar = str(var)
    newComp = type(compareTo)(compareTo)

    return lambda row: getattr(row, newVar) > newComp

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
    'deltaPhiZZ' : {
        'eeee' : '{}(e1_e2_Phi, e3_e4_Phi)'.format(deltaPhiString()),
        'eemm' : '{}(e1_e2_Phi, m1_m2_Phi)'.format(deltaPhiString()),
        'mmmm' : '{}(m1_m2_Phi, m3_m4_Phi)'.format(deltaPhiString()),
        },
    'deltaRZZ' : {
        'eeee' : '{}(e1_e2_Eta, e1_e2_Phi, e3_e4_Eta, e3_e4_Phi)'.format(deltaRString()),
        'eemm' : '{}(e1_e2_Eta, e1_e2_Phi, m1_m2_Eta, m1_m2_Phi)'.format(deltaRString()),
        'mmmm' : '{}(m1_m2_Eta, m1_m2_Phi, m3_m4_Eta, m3_m4_Phi)'.format(deltaRString()),
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
    'z1Mass' : {'eeee':lambda row: row.e1_e2_Mass,
                'mmmm':lambda row: row.m1_m2_Mass,
                'eemm':lambda row: row.e1_e2_Mass if abs(row.e1_e2_Mass - Z_MASS) < abs(row.m1_m2_Mass - Z_MASS) else row.m1_m2_Mass,
                },
    'z2Mass' : {'eeee':lambda row: row.e3_e4_Mass,
                'mmmm':lambda row: row.m3_m4_Mass,
                'eemm':lambda row: row.e1_e2_Mass if abs(row.e1_e2_Mass - Z_MASS) > abs(row.m1_m2_Mass - Z_MASS) else row.m1_m2_Mass,
                },
    'z1Pt' : {'eeee':lambda row: row.e1_e2_Pt,
              'mmmm':lambda row: row.m1_m2_Pt,
              'eemm':lambda row: row.e1_e2_Pt if abs(row.e1_e2_Mass - Z_MASS) < abs(row.m1_m2_Mass - Z_MASS) else row.m1_m2_Pt,
              },
    'z2Pt' : {'eeee':lambda row: row.e3_e4_Pt,
              'mmmm':lambda row: row.m3_m4_Pt,
              'eemm':lambda row: row.e1_e2_Pt if abs(row.e1_e2_Mass - Z_MASS) > abs(row.m1_m2_Mass - Z_MASS) else row.m1_m2_Pt,
              },
    'deltaPhiZZ' : {
        'eeee' : lambda row: _fDPhi(row.e1_e2_Phi, row.e3_e4_Phi),
        'eemm' : lambda row: _fDPhi(row.e1_e2_Phi, row.m1_m2_Phi),
        'mmmm' : lambda row: _fDPhi(row.m1_m2_Phi, row.m3_m4_Phi),
        },
    'deltaRZZ' : {
        'eeee' : lambda row: _fDR(row.e1_e2_Eta, row.e1_e2_Phi, row.e3_e4_Eta, row.e3_e4_Phi),
        'eemm' : lambda row: _fDR(row.e1_e2_Eta, row.e1_e2_Phi, row.m1_m2_Eta, row.m1_m2_Phi),
        'mmmm' : lambda row: _fDR(row.m1_m2_Eta, row.m1_m2_Phi, row.m3_m4_Eta, row.m3_m4_Phi),
        },
    'l1Pt' : {
        'eeee' : lambda row: max(row.e1Pt, row.e3Pt),
        'eemm' : lambda row: max(row.e1Pt, row.m1Pt),
        'mmmm' : lambda row: max(row.m1Pt, row.m3Pt),
        },
    }

_binning = {
    'pt' : [20.*i for i in range(4)] + [100., 140., 200., 300.],
    'nJets' : [6,-0.5,5.5],
    'mass' : [100.] + [200.+50.*i for i in range(5)] + [500.,600.,800.],
    'jet1Pt' : [0., 50., 100., 200., 300., 500.],
    'jet1Eta' : [0., 1.5, 3., 4.7],
    'jet2Pt' : [30., 100., 200., 500.],
    'jet2Eta' : [0., 1.5, 3., 4.7],
    'mjj' : [0., 100., 300., 800.],
    'deltaEtajj' : [6, 0.,6.],
    'z1Mass' : [12, 60., 120.],
    'z2Mass' : [12, 60., 120.],
    'z1Pt' : [i * 25. for i in range(7)] + [200., 300.],
    'z2Pt' : [i * 25. for i in range(7)] + [200., 300.],
    'deltaPhiZZ' : [0., 1.5] + [2.+.25*i for i in range(6)],
    'deltaRZZ' : [6, 0., 6.],
    'l1Pt' : [15, 0., 150.],
    }

_xTitle = {
    'pt' : '4\\ell p_T \\, \\text{(GeV)}',
    'nJets' : '# jets',
    'mass' : 'm_{4\\ell} \\, \\text{(GeV)}',
    'jet1Pt' : 'p_T^\\text{j1} \\, \\text{(GeV)}',
    'jet1Eta' : '\\eta_\\text{j1}',
    'jet2Pt' : 'p_T^\\text{j2} \\, \\text{(GeV)}',
    'jet2Eta' : '\\eta_\\text{j2}',
    'mjj' : 'm_\\text{jj} \\, \\text{(GeV)}',
    'deltaEtajj' : '|\\Delta \\eta_{\\text{jj}}|',
    'z1Mass' : 'm_{\\text{Z}_{1}} \\, \\text{(GeV)}',
    'z2Mass' : 'm_{\\text{Z}_{2}} \\, \\text{(GeV)}',
    'z1Pt' : '\\text{Z}_{1} \\, p_T \\, \\text{(GeV)}',
    'z2Pt' : '\\text{Z}_{2} \\, p_T \\, \\text{(GeV)}',
    'deltaPhiZZ' : '\\Delta \\phi (\\text{Z}_1, \\text{Z}_2)',
    'deltaRZZ' : '\\Delta \\text{R} (\\text{Z}_1, \\text{Z}_2)',
    'l1Pt' : '\\text{Lead. lep. } \\, p_{T} \\, \\text{(GeV)}',
    }

_selections = {
    'pt' : {c:'' for c in _channels},
    'mass' : {c:'' for c in _channels},
    'z1Mass' : {'eeee':'','mmmm':'',
                'eemm':['abs(e1_e2_Mass - {0}) < abs(m1_m2_Mass - {0})'.format(Z_MASS),
                        'abs(e1_e2_Mass - {0}) > abs(m1_m2_Mass - {0})'.format(Z_MASS)],
                },
    'z2Mass' : {'eeee':'','mmmm':'',
                'eemm':['abs(e1_e2_Mass - {0}) > abs(m1_m2_Mass - {0})'.format(Z_MASS),
                        'abs(e1_e2_Mass - {0}) < abs(m1_m2_Mass - {0})'.format(Z_MASS)],
                },
    'z1Pt' : {'eeee':'','mmmm':'',
              'eemm':['abs(e1_e2_Mass - {0}) < abs(m1_m2_Mass - {0})'.format(Z_MASS),
                      'abs(e1_e2_Mass - {0}) > abs(m1_m2_Mass - {0})'.format(Z_MASS)],
              },
    'z2Pt' : {'eeee':'','mmmm':'',
              'eemm':['abs(e1_e2_Mass - {0}) > abs(m1_m2_Mass - {0})'.format(Z_MASS),
                      'abs(e1_e2_Mass - {0}) < abs(m1_m2_Mass - {0})'.format(Z_MASS)],
              },
    'deltaPhiZZ' : {c:'' for c in _channels},
    'deltaRZZ' : {c:'' for c in _channels},
    'l1Pt' : {c:'' for c in _channels},
    }
_selFunctions = {
    'pt' : identityFunction,
    'mass' : identityFunction,
    'z1Mass' : identityFunction, # Z choice done in var function
    'z2Mass' : identityFunction,
    'z1Pt' : identityFunction,
    'z2Pt' : identityFunction,
    'deltaPhiZZ' : identityFunction,
    'deltaRZZ' : identityFunction,
    'l1Pt' : identityFunction,
    }

# do jet variables separately because we have to deal with systematics
for sys in ['', '_jerUp', '_jerDown', '_jesUp','_jesDown']:
    for varName in ['nJets', 'jet1Pt', 'jet1Eta', 'jet2Pt', 'jet2Eta', 'mjj', 'deltaEtajj']:
        doAbs = 'eta' in varName.lower()

        varName += sys
        var = varName
        if doAbs:
            var = 'abs({})'.format(var)

        _variables[varName] = {c:var for c in _channels}
        _varFunctions[varName] = {c:_makeReturner(varName, doAbs) for c in _channels}

        if 'jet1' in varName.lower():
            _selections[varName] = {c:'nJets{} > 0'.format(sys) for c in _channels}
            _selFunctions[varName] = _makeComparator('nJets'+sys, 0)
        elif 'jj' in varName.lower() or 'jet2' in varName.lower():
            _selections[varName] = {c:'nJets{} > 1'.format(sys) for c in _channels}
            _selFunctions[varName] = _makeComparator('nJets'+sys, 1)
        else:
            _selections[varName] = {c:'' for c in _channels}
            _selFunctions[varName] = identityFunction


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
_legParams['l1Pt']['topmargin'] = 0.05

_systSaveFileTemplate = _join(_env['zzt'], 'Analysis', 'savedResults',
                              'unfoldSave_{}Iter.root')

def _normalizeBins(h):
    binUnit = min(h.GetBinWidth(b) for b in range(1,len(h)+1))
    for ib in xrange(1,len(h)+1):
        w = h.GetBinWidth(ib)
        h.SetBinContent(ib, h.GetBinContent(ib) * binUnit / w)
        h.SetBinError(ib, h.GetBinError(ib) * sqrt(binUnit / w))
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
         *varNames, **kwargs):

    style = _Style()

    systSaveFileName = _systSaveFileTemplate.format(nIter)

    systSaveFile = root_open(systSaveFileName, 'UPDATE')

    channels = _channels

    puWeightStr, puWt = puWeight(puWeightFile, '')
    puWeightStrUp, puWtUp = puWeight(puWeightFile, 'up')
    puWeightStrDn, puWtDn = puWeight(puWeightFile, 'dn')
    fPUWt = {'':puWt,'up':puWtUp,'dn':puWtDn}


    true = genZZSamples('zz', inMC, 'smp', lumi)
    reco = {c:zzStackMCOnly(c, inMC, 'smp', puWeightFile, lumi) for c in channels}
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

    altReco = {c:zzStackMCOnly(c, inMC, 'smp', puWeightFile, lumi, amcatnlo=True) for c in channels}
    altTrue = genZZSamples('zz', inMC, 'smp', lumi, amcatnlo=True)

    recoSyst = {}
    for syst in ['eScaleUp', 'eScaleDn', 'eRhoResUp',
                 'eRhoResDn', 'ePhiResUp']:
        recoSyst[syst] = {
            c:zzStackMCOnly(c, inMC.replace('mc_','mc_{}_'.format(syst)),
                            'smp', puWeightFile, lumi)
            for c in ['eeee','eemm']
            }
    for syst in ['mClosureUp','mClosureDn']:
        recoSyst[syst] = {
            c:zzStackMCOnly(c, inMC.replace('mc_','mc_{}_'.format(syst)),
                            'smp', puWeightFile, lumi)
            for c in ['eemm','mmmm']
            }

    for varName in varNames:
        if not hasattr(systSaveFile, varName):
            d = systSaveFile.mkdir(varName)
            d.write()
        varDir = getattr(systSaveFile, varName)

        binning = _binning[varName]

        hTot = {}
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

            hFrame = hUnfolded.empty_clone()

            ### Unfold amc@NLO "data" as a sanity check
            hAlt = sum(altReco[chan].makeHist(var, sel, binning,
                                              perUnitWidth=False).hists)
            unfolderAlt = RooUnfoldIter(response, hAlt, nIter)
            hUnfoldedAlt = asrootpy(unfolderAlt.Hreco())

            #### represent systematic errors as histograms where the bin content
            #### is the systematic error from that source
            hErr = {'up':{},'dn':{}}


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
                    hUnf.SetName('pu_'+sys)
                    hUnf.write()

                    reco[chan].applyWeight(nominalWeight, True)

                hErr[sys]['pu'] = hUnf - hUnfolded
                hErr[sys]['pu'].title = "PU"
                hErr[sys]['pu'].fillstyle = 'solid'
                hErr[sys]['pu'].drawstyle = "hist"
                hErr[sys]['pu'].color = "green"
                hErr[sys]['pu'].legendstyle = 'F'

            # lepton efficiency uncertainty
            for sys in ['up','dn']:
                hUnf = _compatibleHistOrNone(systSaveDir, 'lep_'+sys, hFrame)
                if hUnf is None or redo:
                    wtStr = baseMCWeight(chan, puWeightFile, lepSyst=sys)
                    reco[chan].applyWeight(wtStr, True)

                    res = getResponse(chan, true[chan], reco[chan], bkg[chan], var,
                                      fVar, binning, fPUWt[''], sys,
                                      selectionStr=sel,
                                      selectionFunction=fSel)
                    unf = RooUnfoldIter(res, hData, nIter)

                    systSaveDir.cd()
                    hUnf = asrootpy(unf.Hreco())
                    hUnf.SetName('lep_'+sys)
                    hUnf.write()

                    reco[chan].applyWeight(nominalWeight, True)

                hErr[sys]['lep'] = hUnf - hUnfolded
                hErr[sys]['lep'].title = "Lepton eff."
                hErr[sys]['lep'].fillstyle = 'solid'
                hErr[sys]['lep'].drawstyle = "hist"
                hErr[sys]['lep'].color = "blue"
                hErr[sys]['lep'].legendstyle = 'F'

            # unfolding uncertainty by checking difference with alternate generator
            for sys in ['up','dn']:
                hUnf = _compatibleHistOrNone(systSaveDir, 'generator_'+sys, hFrame)
                if hUnf is None or redo:
                    res = getResponse(chan, altTrue[chan], altReco[chan], bkg[chan], var,
                                      fVar, binning, fPUWt[''],
                                      selectionStr=sel,
                                      selectionFunction=fSel)
                    unf = RooUnfoldIter(res, hData, nIter)

                    systSaveDir.cd()
                    hUnf = asrootpy(unf.Hreco())
                    hUnf.SetName('generator_'+sys)
                    hUnf.write()

                hErr[sys]['generator'] = hUnf - hUnfolded
                hErr[sys]['generator'].title = "Generator choice"
                hErr[sys]['generator'].fillstyle = 'solid'
                hErr[sys]['generator'].drawstyle = "hist"
                hErr[sys]['generator'].color = "magenta"
                hErr[sys]['generator'].legendstyle = 'F'


            # qq NLO uncertainty: +/-1.4% (POWHEG, AN-2016-029)
            # gg LO uncertainty: +18%/-15% (MCFM, AN-2016-029)
            hUnf = _compatibleHistOrNone(systSaveDir, 'xsec_dn', hFrame)
            if hUnf is None or redo:
                reco[chan].applyWeight(nominalWeight, True)
                for s in reco[chan]:
                    if 'GluGlu' in s.name:
                        s.applyWeight('.85')
                    elif 'ZZTo4L' in s.name:
                        s.applyWeight('1.014')
                for n,s in true[chan].itersamples():
                    if 'GluGlu' in n:
                        s.applyWeight('.85')
                    elif 'ZZTo4L' in n:
                        s.applyWeight('1.014')
                res = getResponse(chan, true[chan], reco[chan], bkg[chan], var,
                                  fVar, binning, fPUWt[''],
                                  selectionStr=sel,
                                  selectionFunction=fSel)
                unf = RooUnfoldIter(res, hData, nIter)

                systSaveDir.cd()
                hUnf = asrootpy(unf.Hreco())
                hUnf.SetName('xsec_dn')
                hUnf.write()

                reco[chan].applyWeight(nominalWeight, True)
                true[chan].applyWeight('',True)

            hErr['dn']['xsec'] = hUnf - hUnfolded
            hErr['dn']['xsec'].title = "Cross section"
            hErr['dn']['xsec'].fillstyle = 'solid'
            hErr['dn']['xsec'].drawstyle = "hist"
            hErr['dn']['xsec'].color = "red"
            hErr['dn']['xsec'].legendstyle = 'F'

            hUnf = _compatibleHistOrNone(systSaveDir, 'xsec_up', hFrame)
            if hUnf is None or redo:
                for s in reco[chan]:
                    if 'GluGlu' in s.name:
                        s.applyWeight('1.18')
                    elif 'ZZTo4L' in s.name:
                        s.applyWeight('0.986')
                for n,s in true[chan].itersamples():
                    if 'GluGlu' in n:
                        s.applyWeight('1.18', True)
                    elif 'ZZTo4L' in n:
                        s.applyWeight('0.986', True)
                res = getResponse(chan, true[chan], reco[chan], bkg[chan], var,
                                  fVar, binning, fPUWt[''],
                                  selectionStr=sel,
                                  selectionFunction=fSel)
                unf = RooUnfoldIter(res, hData, nIter)
                systSaveDir.cd()
                hUnf = asrootpy(unf.Hreco())
                hUnf.SetName('xsec_up')
                hUnf.write()

                reco[chan].applyWeight(nominalWeight, True)
                true[chan].applyWeight('',True)

            hErr['up']['xsec'] = hUnf - hUnfolded
            hErr['up']['xsec'].title = "Cross section"
            hErr['up']['xsec'].fillstyle = 'solid'
            hErr['up']['xsec'].drawstyle = "hist"
            hErr['up']['xsec'].color = "red"
            hErr['up']['xsec'].legendstyle = 'F'


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
                    hUnf.SetName('lumi_'+sys)
                    hUnf.write()

                    reco[chan].applyWeight(nominalWeight, True)

                hErr[sys]['lumi'] = hUnf - hUnfolded
                hErr[sys]['lumi'].title = "Luminosity"
                hErr[sys]['lumi'].fillstyle = 'solid'
                hErr[sys]['lumi'].drawstyle = "hist"
                hErr[sys]['lumi'].color = "orange"
                hErr[sys]['lumi'].legendstyle = 'F'


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
                        hUnf.SetName(lep+'FR_'+sys)
                        hUnf.write()

                    hErr[sys][lep+'FR'] = hUnf - hUnfolded
                    if lep == 'e':
                        hErr[sys][lep+'FR'].title = "Electron Fake Rate"
                        hErr[sys][lep+'FR'].color = "#00cc99"
                    elif lep == 'm':
                        hErr[sys][lep+'FR'].title = "Muon Fake Rate"
                        hErr[sys][lep+'FR'].color = "#00ff00"
                    hErr[sys][lep+'FR'].fillstyle = 'solid'
                    hErr[sys][lep+'FR'].drawstyle = "hist"
                    hErr[sys][lep+'FR'].legendstyle = 'F'


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
                            hUnf.SetName(sys+'_'+shift)
                            hUnf.write()

                        hErr[shift][sys] = hUnf - hUnfolded
                        hErr[shift][sys].fillstyle = 'solid'
                        hErr[shift][sys].drawstyle = "hist"
                        hErr[shift][sys].legendstyle = 'F'


                    hErr[shift]['jer'].title = "Jet energy resolution"
                    hErr[shift]['jer'].color = "cyan"
                    hErr[shift]['jes'].title = "Jet energy scale"
                    hErr[shift]['jes'].color = "darkblue"


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
                        hUnf.SetName('ees_'+sys)
                        hUnf.write()

                    hErr[sys]['ees'] = hUnf - hUnfolded
                    hErr[sys]['ees'].title = "Electron energy scale"
                    hErr[sys]['ees'].fillstyle = 'solid'
                    hErr[sys]['ees'].drawstyle = "hist"
                    hErr[sys]['ees'].color = "purple"
                    hErr[sys]['ees'].legendstyle = 'F'

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
                        hUnf.SetName('eerRho_'+sys)
                        hUnf.write()

                    hErr[sys]['eerRho'] = hUnf - hUnfolded
                    hErr[sys]['eerRho'].title = "Electron energy res. (rho)"
                    hErr[sys]['eerRho'].fillstyle = 'solid'
                    hErr[sys]['eerRho'].drawstyle = "hist"
                    hErr[sys]['eerRho'].color = "lavender"
                    hErr[sys]['eerRho'].legendstyle = 'F'

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
                    hUnf.SetName('eerPhi_up')
                    hUnf.write()

                hErr['up']['eerPhi'] = hUnf - hUnfolded
                hErr['up']['eerPhi'].title = "Electron energy res. (phi)"
                hErr['up']['eerPhi'].fillstyle = 'solid'
                hErr['up']['eerPhi'].drawstyle = "hist"
                hErr['up']['eerPhi'].color = "violet"
                hErr['up']['eerPhi'].legendstyle = 'F'

                hErr['dn']['eerPhi'] = hErr['up']['eerPhi'].clone()
                hErr['dn']['eerPhi'].title = "Electron energy res. (phi)"
                hErr['dn']['eerPhi'].fillstyle = 'solid'
                hErr['dn']['eerPhi'].drawstyle = "hist"
                hErr['dn']['eerPhi'].color = "violet"
                hErr['dn']['eerPhi'].legendstyle = 'F'


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
                        hUnf.SetName('mClosure_'+sys)
                        hUnf.write()

                    hErr[sys]['mClosure'] = hUnf - hUnfolded
                    hErr[sys]['mClosure'].title = "Muon calibration"
                    hErr[sys]['mClosure'].fillstyle = 'solid'
                    hErr[sys]['mClosure'].drawstyle = "hist"
                    hErr[sys]['mClosure'].color = "#c61aff"
                    hErr[sys]['mClosure'].legendstyle = 'F'


            ### Get the total uncertainties
            # this method of assigning up vs down should be revisited
            hUncUp = hUnfolded.empty_clone()
            hUncDn = hUnfolded.empty_clone()
            uncTypes = hErr['up'].keys()
            for bUncUp, bUncDn, allUncUp, allUncDn in zip(hUncUp, hUncDn,
                                                          zip(*(hErr['up'][t] for t in uncTypes)),
                                                          zip(*(hErr['dn'][t] for t in uncTypes))):
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

            # Note: hUncUp and hUncDn are both squared at this point

            if 'uncUpSqr' not in hTot:
                hTot['uncUpSqr'] = hUncUp.empty_clone()
            hTot['uncUpSqr'] += hUncUp
            if 'uncDnSqr' not in hTot:
                hTot['uncDnSqr'] = hUncDn.empty_clone()
            hTot['uncDnSqr'] += hUncDn

            for bUp, bDn in zip(hUncUp, hUncDn):
                bUp.value = sqrt(bUp.value)
                bDn.value = sqrt(bDn.value)


            # Make all error histograms positive (only matters for error plot)
            # and make uncertainties fractional (as a percent)
            for sys in hErr.values():
                for h in sys.values():
                    h /= hUnfolded
                    h *= 100.
                    for b in h:
                        b.value = abs(b.value)

            # Make plots of uncertainties (added linearly)
            cErrUp = Canvas(1000,1000)
            errStackUp = HistStack(hErr['up'].values(), drawstyle = 'histnoclear')
            draw(errStackUp, cErrUp, xtitle=_xTitle[varName], ytitle="+Error (%)",yerror_in_padding=False)
            leg = makeLegend(cErrUp, *hErr['up'].values(), leftmargin=0.25,
                             entryheight=.02, entrysep=.007, textsize=.022,
                             rightmargin=.25)
            leg.Draw('same')
            style.setCMSStyle(cErrUp, '', dataType='Preliminary', intLumi=lumi)
            cErrUp.Print(_join(plotDir, 'errUp_{}_{}.png'.format(varName, chan)))

            cErrDn = Canvas(1000,1000)
            errStackDn = HistStack(hErr['dn'].values(), drawstyle = 'histnoclear')
            draw(errStackDn, cErrDn, xtitle=_xTitle[varName], ytitle="-Error (%)",yerror_in_padding=False)
            leg = makeLegend(cErrDn, *hErr['dn'].values(), leftmargin=0.25,
                             entryheight=.02, entrysep=.007, textsize=.022,
                             rightmargin=.25)
            leg.Draw('same')
            style.setCMSStyle(cErrDn, '', dataType='Preliminary', intLumi=lumi)
            cErrDn.Print(_join(plotDir, 'errDown_{}_{}.png'.format(varName, chan)))


            ### plot
            hUnfolded.color = 'black'
            hUnfolded.drawstyle = 'PE'
            hUnfolded.legendstyle = 'LPE'
            hUnfolded.title = 'Unfolded data + stat. unc.'
            if 'unfolded' not in hTot:
                hTot['unfolded'] = hUnfolded.empty_clone()
            hTot['unfolded'] += hUnfolded
            _normalizeBins(hUnfolded)

            hUnfoldedAlt.color = 'r'
            hUnfoldedAlt.drawstyle = 'hist'
            hUnfoldedAlt.fillstyle = 'hollow'
            hUnfoldedAlt.legendstyle = 'L'
            hUnfoldedAlt.title = 'Unfolded MG5_aMC@NLO+MCFM'
            if 'unfoldedAlt' not in hTot:
                hTot['unfoldedAlt'] = hUnfoldedAlt.empty_clone()
            hTot['unfoldedAlt'] += hUnfoldedAlt
            _normalizeBins(hUnfoldedAlt)

            hTrue = true[chan].makeHist(var, sel, binning)
            hTrue.color = 'blue'
            hTrue.drawstyle = 'hist'
            hTrue.fillstyle = 'hollow'
            hTrue.legendstyle = 'L'
            hTrue.title = 'POWHEG (true) [training sample]'

            hTrueAlt = altTrue[chan].makeHist(var, sel, binning)
            hTrueAlt.color = 'magenta'
            hTrueAlt.drawstyle = 'hist'
            hTrueAlt.fillstyle = 'hollow'
            hTrueAlt.legendstyle = 'L'
            hTrueAlt.title = 'MG5_aMC@NLO (true)'

            _normalizeBins(hUncUp)
            _normalizeBins(hUncDn)
            errorBand = makeErrorBand(hUnfolded, hUncUp, hUncDn)

            cUnf = Canvas(1000,1000)
            draw([hTrue, hTrueAlt, hUnfoldedAlt, hUnfolded, errorBand], cUnf,
                 xtitle=_xTitle[varName], ytitle='Events',yerror_in_padding=False)
            leg = makeLegend(cUnf, hTrueAlt, hUnfoldedAlt, hTrue, errorBand,
                             hUnfolded, **_legParams[varName])
            leg.Draw("same")

            style.setCMSStyle(cUnf, '', dataType='Preliminary', intLumi=lumi)
            cUnf.Print(_join(plotDir, "unfold_{}_{}.png".format(varName, chan)))

            cRes = Canvas(1000,1000)
            hRes = asrootpy(response.Hresponse())
            hRes.drawstyle = 'colztext'
            hRes.draw()
            style.setCMSStyle(cRes, '', dataType='Preliminary Simulation', intLumi=lumi)
            cRes.Print(_join(plotDir, "response_{}_{}.png".format(varName, chan)))

            cCov = Canvas(1000,1000)
            covariance = unfolder.Ereco(2) # TMatrixD
            covariance.Draw("colztext")
            style.setCMSStyle(cCov, '', dataType='Preliminary', intLumi=lumi)
            cCov.Print(_join(plotDir, "covariance_{}_{}.png".format(varName, chan)))

        hTot['unfolded'].color = 'black'
        hTot['unfolded'].drawstyle = 'PE'
        hTot['unfolded'].legendstyle = 'LPE'
        hTot['unfolded'].title = 'Unfolded data + stat. unc.'
        _normalizeBins(hTot['unfolded'])

        hTot['unfoldedAlt'].color = 'r'
        hTot['unfoldedAlt'].drawstyle = 'hist'
        hTot['unfoldedAlt'].fillstyle = 'hollow'
        hTot['unfoldedAlt'].legendstyle = 'L'
        hTot['unfoldedAlt'].title = 'Unfolded MG5_aMC@NLO+MCFM'
        _normalizeBins(hTot['unfoldedAlt'])

        hTrue = true.makeHist(_variables[varName], _selections[varName], binning)
        hTrue.color = 'blue'
        hTrue.drawstyle = 'hist'
        hTrue.fillstyle = 'hollow'
        hTrue.legendstyle = 'L'
        hTrue.title = 'POWHEG (true) [training sample]'

        hTrueAlt = altTrue.makeHist(_variables[varName], _selections[varName], binning)
        hTrueAlt.color = 'magenta'
        hTrueAlt.drawstyle = 'hist'
        hTrueAlt.fillstyle = 'hollow'
        hTrueAlt.legendstyle = 'L'
        hTrueAlt.title = 'MG5_aMC@NLO (true)'

        hUncUp = hTot['uncUpSqr'].clone()
        for b in hUncUp:
            b.value = sqrt(b.value)
        hUncDn = hTot['uncDnSqr'].clone()
        for b in hUncDn:
            b.value = sqrt(b.value)

        _normalizeBins(hUncUp)
        _normalizeBins(hUncDn)
        errorBand = makeErrorBand(hTot['unfolded'], hUncUp, hUncDn)

        cUnf = Canvas(1000,1000)
        draw([hTrue, hTrueAlt, hTot['unfoldedAlt'], errorBand, hTot['unfolded']],
             cUnf, xtitle=_xTitle[varName], ytitle='Events',yerror_in_padding=False)
        leg = makeLegend(cUnf, hTrueAlt, hTot['unfoldedAlt'], hTrue, errorBand,
                         hTot['unfolded'], **_legParams[varName])
        leg.Draw("same")

        style.setCMSStyle(cUnf, '', dataType='Preliminary', intLumi=lumi)
        cUnf.Print(_join(plotDir, "unfold_{}.png".format(varName)))

    systSaveFile.close()

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
    parser.add_argument('--variables', type=str, nargs='*',
                        default=_varList,
                        help=('Names of variables to use. If not specified, '
                              'all are used ({})').format(', '.join(_varList)))

    args=parser.parse_args()
    main(args.dataDir, args.mcDir, args.plotDir, args.fakeRateFile,
         args.puWeightFile, args.lumi, args.nIter, args.redo, *args.variables)
