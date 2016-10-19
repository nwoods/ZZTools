
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
from Utilities import WeightStringMaker
from Analysis.setupStandardSamples import *
from Analysis.unfoldingHelpers import getResponse
from Analysis.weightHelpers import puWeight, baseMCWeight
from Metadata.metadata import sampleInfo

from os import environ as _env
from os.path import join as _join
from math import sqrt


_variables = {
    'pt' : 'Pt',
    'nJets' : 'nJets',
    'mass' : 'Mass',
    'jet1Pt' : 'jet1Pt',
    'jet1Eta' : 'abs(jet1Eta)',
    'mjj' : 'mjj',
    }

_binning = {
    'pt' : [20.*i for i in range(4)] + [100., 140., 200., 300.],
    'nJets' : [6,-0.5,5.5],
    'mass' : [100.] + [200.+50.*i for i in range(5)] + [500.,600.,800.],
    'jet1Pt' : [0., 50., 100., 200., 300., 500.],
    'jet1Eta' : [0., 1.5, 3., 4.7],
    'mjj' : [0., 100., 300., 800.],
    }

_xTitle = {
    'pt' : '4\\ell p_T',
    'nJets' : '# jets',
    'mass' : 'm_{4\\ell}',
    'jet1Pt' : 'p_T^\\text{j1}',
    'jet1Pt' : '|\\eta^\\text{j1}|',
    'jet1Pt' : 'm_\\text{jj}',
    }

_selections = {
    'pt' : ('', None),
    'nJets' : ('', None),
    'mass' : ('', None),
    'jet1Pt' : ('nJets > 0', lambda row: row.nJets > 0),
    'mjj' : ('nJets >= 2', lambda row: row.nJets >= 2),
    }

_systSaveFileName = _join(_env['zzt'], 'Analysis', 'savedResults',
                          'unfoldSave.root')

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

    systSaveFile = root_open(_systSaveFileName, 'UPDATE')

    channels = ['eeee','eemm', 'mmmm']

    puWeightStr, puWt = puWeight(puWeightFile, '')
    puWeightStrUp, puWtUp = puWeight(puWeightFile, 'up')
    puWeightStrDn, puWtDn = puWeight(puWeightFile, 'dn')
    fPUWt = {'':puWt,'up':puWtUp,'dn':puWtDn}


    true = genZZSamples('zz', inMC, 'smp', lumi)
    reco = zzStackMCOnly('zz', inMC, 'smp', puWeightFile, lumi)
    recoByChan = {c:SampleStack('stack', c, reco.getSamplesForChannel(c)) for c in channels}
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

    altReco = zzStackMCOnly('zz', inMC, 'smp', puWeightFile, lumi, amcatnlo=True)
    altRecoByChan = {c:SampleStack('stack', c, altReco.getSamplesForChannel(c)) for c in channels}
    altTrue = genZZSamples('zz', inMC, 'smp', lumi, amcatnlo=True)

    recoSyst = {}
    recoSystByChan = {}
    for syst in ['eScaleUp', 'eScaleDn', 'eRhoResUp',
                 'eRhoResDn', 'ePhiResUp']:
        recoSyst[syst] = zzStackMCOnly('eeee,eemm',
                                       inMC.replace('mc_','mc_{}_'.format(syst)),
                                       'smp', puWeightFile, lumi)
        recoSystByChan[syst] = {
            c:SampleStack('stack', c,
                          recoSyst[syst].getSamplesForChannel(c)
                          ) for c in ['eeee','eemm']
            }
    for syst in ['mClosureUp','mClosureDn']:
        recoSyst[syst] = zzStackMCOnly('eemm,mmmm',
                                       inMC.replace('mc_','mc_{}_'.format(syst)),
                                       'smp', puWeightFile, lumi)
        recoSystByChan[syst] = {
            c:SampleStack('stack', c,
                          recoSyst[syst].getSamplesForChannel(c)
                          ) for c in ['eemm','mmmm']
            }

    for varName in varNames:
        if not hasattr(systSaveFile, varName):
            d = systSaveFile.mkdir(varName)
            d.write()
        varDir = getattr(systSaveFile, varName)

        var = _variables[varName]
        hTot = {}
        for chan in channels:
            if not hasattr(varDir, chan):
                d = varDir.mkdir(chan)
                d.write()
            systSaveDir = getattr(varDir, chan)

            # regular weight, no systematics. Apply just in case.
            nominalWeight = baseMCWeight(chan, puWeightFile)
            recoByChan[chan].applyWeight(nominalWeight, True)

            response = getResponse(chan, true[chan], recoByChan[chan], bkg[chan], var,
                                   _binning[varName], puWt, selectionStr=_selections[varName][0],
                                   selectionFunction=_selections[varName][1])

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

            hBkg = bkg[chan].makeHist(var, _selections[varName][0], _binning[varName],
                                      perUnitWidth=False)
            hData = data[chan].makeHist(var, _selections[varName][0], _binning[varName],
                                        perUnitWidth=False)
            hData -= hBkg

            unfolder = RooUnfoldIter(response, hData, nIter)

            #print chan
            #unfolders[chan].PrintTable(cout)

            hUnfolded = asrootpy(unfolder.Hreco())

            hFrame = hUnfolded.empty_clone()

            ### Unfold amc@NLO "data" as a sanity check
            hAlt = sum(altRecoByChan[chan].makeHist(var, _selections[varName][0],
                                                    _binning[varName],
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
                    recoByChan[chan].applyWeight(wtStr, True)

                    res = getResponse(chan, true[chan], recoByChan[chan], bkg[chan], var,
                                      _binning[varName], fPUWt[sys],
                                      selectionStr=_selections[varName][0],
                                      selectionFunction=_selections[varName][1])
                    unf = RooUnfoldIter(res, hData, nIter)

                    systSaveDir.cd()
                    hUnf = asrootpy(unf.Hreco())
                    hUnf.SetName('pu_'+sys)
                    hUnf.write()

                    recoByChan[chan].applyWeight(nominalWeight, True)

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
                    recoByChan[chan].applyWeight(wtStr, True)

                    res = getResponse(chan, true[chan], recoByChan[chan], bkg[chan], var,
                                      _binning[varName], fPUWt[''], sys,
                                      selectionStr=_selections[varName][0],
                                      selectionFunction=_selections[varName][1])
                    unf = RooUnfoldIter(res, hData, nIter)

                    systSaveDir.cd()
                    hUnf = asrootpy(unf.Hreco())
                    hUnf.SetName('lep_'+sys)
                    hUnf.write()

                    recoByChan[chan].applyWeight(nominalWeight, True)

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
                    res = getResponse(chan, altTrue[chan], altRecoByChan[chan], bkg[chan], var,
                                      _binning[varName], fPUWt[''],
                                      selectionStr=_selections[varName][0],
                                      selectionFunction=_selections[varName][1])
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
                recoByChan[chan].applyWeight(nominalWeight, True)
                for s in recoByChan[chan]:
                    if 'GluGlu' in s.name:
                        s.applyWeight('.85')
                    elif 'ZZTo4L' in s.name:
                        s.applyWeight('1.014')
                for n,s in true[chan].itersamples():
                    if 'GluGlu' in n:
                        s.applyWeight('.85')
                    elif 'ZZTo4L' in n:
                        s.applyWeight('1.014')
                res = getResponse(chan, true[chan], recoByChan[chan], bkg[chan], var,
                                  _binning[varName], fPUWt[''],
                                  selectionStr=_selections[varName][0],
                                  selectionFunction=_selections[varName][1])
                unf = RooUnfoldIter(res, hData, nIter)

                systSaveDir.cd()
                hUnf = asrootpy(unf.Hreco())
                hUnf.SetName('xsec_dn')
                hUnf.write()

                recoByChan[chan].applyWeight(nominalWeight, True)
                true[chan].applyWeight('',True)

            hErr['dn']['xsec'] = hUnf - hUnfolded
            hErr['dn']['xsec'].title = "Cross section"
            hErr['dn']['xsec'].fillstyle = 'solid'
            hErr['dn']['xsec'].drawstyle = "hist"
            hErr['dn']['xsec'].color = "red"
            hErr['dn']['xsec'].legendstyle = 'F'

            hUnf = _compatibleHistOrNone(systSaveDir, 'xsec_up', hFrame)
            if hUnf is None or redo:
                for s in recoByChan[chan]:
                    if 'GluGlu' in s.name:
                        s.applyWeight('1.18')
                    elif 'ZZTo4L' in s.name:
                        s.applyWeight('0.986')
                for n,s in true[chan].itersamples():
                    if 'GluGlu' in n:
                        s.applyWeight('1.18', True)
                    elif 'ZZTo4L' in n:
                        s.applyWeight('0.986', True)
                res = getResponse(chan, true[chan], recoByChan[chan], bkg[chan], var,
                                  _binning[varName], fPUWt[''],
                                  selectionStr=_selections[varName][0],
                                  selectionFunction=_selections[varName][1])
                unf = RooUnfoldIter(res, hData, nIter)
                systSaveDir.cd()
                hUnf = asrootpy(unf.Hreco())
                hUnf.SetName('xsec_up')
                hUnf.write()

                recoByChan[chan].applyWeight(nominalWeight, True)
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
                    recoByChan[chan].applyWeight(str(1.+lumiUnc))
                    true[chan].applyWeight(str(1.+lumiUnc), True)
                    res = getResponse(chan, true[chan], recoByChan[chan], bkg[chan], var,
                                      _binning[varName], fPUWt[''],
                                      selectionStr=_selections[varName][0],
                                      selectionFunction=_selections[varName][1])
                    unf = RooUnfoldIter(res, hData, nIter)

                    systSaveDir.cd()
                    hUnf = asrootpy(unf.Hreco())
                    hUnf.SetName('lumi_'+sys)
                    hUnf.write()

                    recoByChan[chan].applyWeight(nominalWeight, True)

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
                        res = getResponse(chan, true[chan], recoByChan[chan],
                                          bkgSyst[lep+sys][chan], var,
                                          _binning[varName], fPUWt[''],
                                          selectionStr=_selections[varName][0],
                                          selectionFunction=_selections[varName][1])
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


            # Jet energy scale and resolution uncertainties (nJets only)
            if varName == 'nJets':
                for sys in ['up','dn']:
                    sysStr = 'Up' if sys == 'up' else 'Down'

                    # JER
                    hUnf = _compatibleHistOrNone(systSaveDir, 'jer_'+sys, hFrame)
                    if hUnf is None or redo:
                        res = getResponse(chan, true[chan], recoByChan[chan], bkg[chan],
                                          'nJets_jer'+sysStr, _binning[varName], fPUWt[''],
                                          altVar='nJets',
                                          selectionStr=_selections[varName][0],
                                          selectionFunction=_selections[varName][1])
                        unf = RooUnfoldIter(res, hData, nIter)

                        systSaveDir.cd()
                        hUnf = asrootpy(unf.Hreco())
                        hUnf.SetName('jer_'+sys)
                        hUnf.write()

                    hErr[sys]['jer'] = hUnf - hUnfolded
                    hErr[sys]['jer'].title = "Jet energy resolution"
                    hErr[sys]['jer'].fillstyle = 'solid'
                    hErr[sys]['jer'].drawstyle = "hist"
                    hErr[sys]['jer'].color = "cyan"
                    hErr[sys]['jer'].legendstyle = 'F'

                    # JES
                    hUnf = _compatibleHistOrNone(systSaveDir, 'jes_'+sys, hFrame)
                    if hUnf is None or redo:
                        res = getResponse(chan, true[chan], recoByChan[chan], bkg[chan],
                                          'nJets_jes'+sysStr, _binning[varName], fPUWt[''],
                                          altVar='nJets',
                                          selectionStr=_selections[varName][0],
                                          selectionFunction=_selections[varName][1])
                        unf = RooUnfoldIter(res, hData, nIter)

                        systSaveDir.cd()
                        hUnf = asrootpy(unf.Hreco())
                        hUnf.SetName('jes_'+sys)
                        hUnf.write()

                    hErr[sys]['jes'] = hUnf - hUnfolded
                    hErr[sys]['jes'].title = "Jet energy scale"
                    hErr[sys]['jes'].fillstyle = 'solid'
                    hErr[sys]['jes'].drawstyle = "hist"
                    hErr[sys]['jes'].color = "darkblue"
                    hErr[sys]['jes'].legendstyle = 'F'

            if 'e' in chan:
                for sys in ['up', 'dn']:
                    # EES
                    hUnf = _compatibleHistOrNone(systSaveDir, 'ees_'+sys, hFrame)
                    if hUnf is None or redo:
                        res = getResponse(chan, true[chan],
                                          recoSystByChan['eScale'+sys[0].upper()+sys[1:]][chan],
                                          bkg[chan], var, _binning[varName],
                                          fPUWt[''],
                                          selectionStr=_selections[varName][0],
                                          selectionFunction=_selections[varName][1])
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
                                          recoSystByChan['eRhoRes'+sys[0].upper()+sys[1:]][chan],
                                          bkg[chan], var, _binning[varName],
                                          fPUWt[''],
                                          selectionStr=_selections[varName][0],
                                          selectionFunction=_selections[varName][1])
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
                                      recoSystByChan['ePhiResUp'][chan],
                                      bkg[chan], var, _binning[varName],
                                      fPUWt[''],
                                      selectionStr=_selections[varName][0],
                                      selectionFunction=_selections[varName][1])
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
                                          recoSystByChan['mClosure'+sys[0].upper()+sys[1:]][chan],
                                          bkg[chan], var, _binning[varName],
                                          fPUWt[''],
                                          selectionStr=_selections[varName][0],
                                          selectionFunction=_selections[varName][1])
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

            hTrue = true[chan].makeHist(var, _selections[varName][0], _binning[varName])
            hTrue.color = 'blue'
            hTrue.drawstyle = 'hist'
            hTrue.fillstyle = 'hollow'
            hTrue.legendstyle = 'L'
            hTrue.title = 'POWHEG (true) [training sample]'

            hTrueAlt = altTrue[chan].makeHist(var, _selections[varName][0], _binning[varName])
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
            leg = makeLegend(cUnf, hUnfolded, errorBand, hUnfoldedAlt, hTrue,
                             hTrueAlt, textsize=.023, leftmargin=0.3)
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

        hTrue = true.makeHist(var, _selections[varName][0], _binning[varName])
        hTrue.color = 'blue'
        hTrue.drawstyle = 'hist'
        hTrue.fillstyle = 'hollow'
        hTrue.legendstyle = 'L'
        hTrue.title = 'POWHEG (true) [training sample]'

        hTrueAlt = altTrue.makeHist(var, _selections[varName][0], _binning[varName])
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
        draw([hTrue, hTrueAlt, hTot['unfoldedAlt'], hTot['unfolded'], errorBand],
             cUnf, xtitle=_xTitle[varName], ytitle='Events',yerror_in_padding=False)
        leg = makeLegend(cUnf, hTot['unfolded'], errorBand, hTot['unfoldedAlt'], hTrue,
                         hTrueAlt, textsize=.023, leftmargin=0.3)
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
                        default=list(_variables.keys()),
                        help=('Names of variables to use. If not specified, '
                              'all are used ({})').format(', '.join(_variables.keys())))

    args=parser.parse_args()
    main(args.dataDir, args.mcDir, args.plotDir, args.fakeRateFile,
         args.puWeightFile, args.lumi, args.nIter, args.redo, *args.variables)
