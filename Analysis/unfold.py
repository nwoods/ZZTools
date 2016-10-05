
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

from os.path import join as _join
from math import sqrt

inData = 'uwvvNtuples_data_08sep2016'
inMC = 'uwvvNtuples_mc_08sep2016'

plotDir = '/afs/cern.ch/user/n/nawoods/www/UWVVPlots/unfold'

fakeRateFile = 'fakeRate_08sep2016'

#fakeDataFile = '/afs/cern.ch/user/n/nawoods/ZZTools/ZZTo4LTest_amcatnlo.root'

puWeightFile = 'puWeight_69200_08sep2016.root'

style = _Style()

lumi = 15937.

nIter = 8

channels = ['eeee','eemm', 'mmmm']

variables = {
    'pt' : 'Pt',
    #'nJets' : 'nJets',
    'mass' : 'Mass',
    }

binning = {
    'pt' : [20.*i for i in range(4)] + [100., 140., 200., 300.],
    'nJets' : [6,-0.5,5.5],
    'mass' : [100.] + [200.+50.*i for i in range(5)] + [500.,600.,800.],
    }

xTitle = {
    'pt' : '4\\ell p_T',
    'nJets' : '# jets',
    'mass' : 'm_{4\\ell}',
    }

# Making PU weight strings makes compiled functions; get them
puWeightStr = puWeight(puWeightFile, '')
from rootpy.ROOT import puWeight0 as puWt
puWeightStrUp = puWeight(puWeightFile, 'up')
from rootpy.ROOT import puWeight1 as puWtUp
puWeightStrDn = puWeight(puWeightFile, 'dn')
from rootpy.ROOT import puWeight2 as puWtDn
fPUWt = {'':puWt,'up':puWtUp,'dn':puWtDn}


def normalizeBins(h):
    binUnit = min(h.GetBinWidth(b) for b in range(1,len(h)+1))
    for ib in xrange(1,len(h)+1):
        w = h.GetBinWidth(ib)
        h.SetBinContent(ib, h.GetBinContent(ib) * binUnit / w)
        h.SetBinError(ib, h.GetBinError(ib) * sqrt(binUnit / w))
        if h.GetBinError(ib) > h.GetBinContent(ib):
            h.SetBinError(ib, h.GetBinContent(ib))
    h.sumw2()


true = genZZSamples('zz', inMC, 'smp', lumi)
reco = zzStackMCOnly('zz', inMC, 'smp', puWeightFile, lumi)
recoByChan = {c:SampleStack('stack', c, reco.getSamplesForChannel(c)) for c in channels}
bkg = standardZZBkg('zz', inData, inMC, 'smp', puWeightFile,
                    fakeRateFile, lumi)
data = standardZZData('zz', inData, 'smp')

altReco = zzStackMCOnly('zz', inMC, 'smp', puWeightFile, lumi, amcatnlo=True)
altRecoByChan = {c:SampleStack('stack', c, altReco.getSamplesForChannel(c)) for c in channels}
altTrue = genZZSamples('zz', inMC, 'smp', lumi, amcatnlo=True)


for varName, var in variables.iteritems():
    hTot = {}
    for chan in channels:

        # regular weight, no systematics. Apply just in case.
        nominalWeight = baseMCWeight(chan, puWeightFile)
        recoByChan[chan].applyWeight(nominalWeight, True)

        response = getResponse(chan, true[chan], recoByChan[chan], bkg[chan], var, 
                               binning[varName], puWt)

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

        hBkg = bkg[chan].makeHist(var, '', binning[varName], 
                                  perUnitWidth=False)
        hData = data[chan].makeHist(var, '', binning[varName], 
                                    perUnitWidth=False)
        hData -= hBkg

        unfolder = RooUnfoldIter(response, hData, nIter)
        
        #print chan
        #unfolders[chan].PrintTable(cout)

        hUnfolded = asrootpy(unfolder.Hreco())

        
        ### Unfold amc@NLO "data" as a sanity check
        hAlt = sum(altRecoByChan[chan].makeHist(var, '', 
                                                binning[varName], 
                                                perUnitWidth=False).hists)
        unfolderAlt = RooUnfoldIter(response, hAlt, nIter)
        hUnfoldedAlt = asrootpy(unfolderAlt.Hreco())

        #### represent systematic errors as histograms where the bin content
        #### is the systematic error from that source
        hErr = {'up':{},'dn':{}}


        # PU reweight uncertainty
        for sys in ['up','dn']:
            wtStr = baseMCWeight(chan, puWeightFile, puSyst=sys)
            recoByChan[chan].applyWeight(wtStr, True)

            res = getResponse(chan, true[chan], recoByChan[chan], bkg[chan], var,
                              binning[varName], fPUWt[sys])
            unf = RooUnfoldIter(res, hData, nIter)
            hUnf = asrootpy(unf.Hreco())
            hErr[sys]['pu'] = hUnf - hUnfolded
            hErr[sys]['pu'].title = "PU"
            hErr[sys]['pu'].fillstyle = 'solid'
            hErr[sys]['pu'].drawstyle = "hist"
            hErr[sys]['pu'].color = "green"
            hErr[sys]['pu'].legendstyle = 'F'

        # lepton efficiency uncertainty
        # for sys in ['up','dn']:
        #     wtStr = baseMCWeight(chan, puWeightFile, lepSyst=sys)
        #     recoByChan[chan].applyWeight(wtStr, True)
        # 
        #     res = getResponse(chan, true[chan], recoByChan[chan], bkg[chan], var,
        #                       binning[varName], fPUWt[''], sys)
        #     unf = RooUnfoldIter(res, hData, nIter)
        #     hUnf = asrootpy(unf.Hreco())
        #     hErr[sys]['lep'] = hUnf - hUnfolded
        #     hErr[sys]['lep'].title = "Lepton eff."
        #     hErr[sys]['lep'].fillstyle = 'solid'
        #     hErr[sys]['lep'].drawstyle = "hist"
        #     hErr[sys]['lep'].color = "blue"
        #     hErr[sys]['lep'].legendstyle = 'F'

        # unfolding uncertainty by checking difference with alternate generator
        res = getResponse(chan, altTrue[chan], altRecoByChan[chan], bkg[chan], var,
                          binning[varName], fPUWt[''])
        unf = RooUnfoldIter(res, hData, nIter)
        hUnf = asrootpy(unf.Hreco())
        hErr['up']['generator'] = hUnf - hUnfolded
        hErr['up']['generator'].title = "Generator choice"
        hErr['up']['generator'].fillstyle = 'solid'
        hErr['up']['generator'].drawstyle = "hist"
        hErr['up']['generator'].color = "magenta"
        hErr['up']['generator'].legendstyle = 'F'

        hErr['dn']['generator'] = hErr['up']['generator'].clone()
        hErr['dn']['generator'].title = "Generator choice"
        hErr['dn']['generator'].fillstyle = 'solid'
        hErr['dn']['generator'].drawstyle = "hist"
        hErr['dn']['generator'].color = "magenta"
        hErr['dn']['generator'].legendstyle = 'F'

        # qq NLO uncertainty: +/-1.4% (POWHEG, AN-2016-029)
        # gg LO uncertainty: +18%/-15% (MCFM, AN-2016-029)
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
                          binning[varName], fPUWt[''])
        unf = RooUnfoldIter(res, hData, nIter)
        hUnf = asrootpy(unf.Hreco())
        hErr['dn']['xsec'] = hUnf - hUnfolded
        hErr['dn']['xsec'].title = "Cross section"
        hErr['dn']['xsec'].fillstyle = 'solid'
        hErr['dn']['xsec'].drawstyle = "hist"
        hErr['dn']['xsec'].color = "red"
        hErr['dn']['xsec'].legendstyle = 'F'

        recoByChan[chan].applyWeight(nominalWeight, True)
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
                          binning[varName], fPUWt[''])
        unf = RooUnfoldIter(res, hData, nIter)
        hUnf = asrootpy(unf.Hreco())
        hErr['up']['xsec'] = hUnf - hUnfolded
        hErr['up']['xsec'].title = "Cross section"
        hErr['up']['xsec'].fillstyle = 'solid'
        hErr['up']['xsec'].drawstyle = "hist"
        hErr['up']['xsec'].color = "red"
        hErr['up']['xsec'].legendstyle = 'F'

        recoByChan[chan].applyWeight(nominalWeight, True)
        true[chan].applyWeight('',True)

        # luminosity uncertainty
        lumiUnc = .062
        recoByChan[chan].applyWeight(str(1.+lumiUnc))
        true[chan].applyWeight(str(1.+lumiUnc), True)
        res = getResponse(chan, true[chan], recoByChan[chan], bkg[chan], var,
                          binning[varName], fPUWt[''])
        unf = RooUnfoldIter(res, hData, nIter)
        hUnf = asrootpy(unf.Hreco())
        hErr['up']['lumi'] = hUnf - hUnfolded
        hErr['up']['lumi'].title = "Luminosity"
        hErr['up']['lumi'].fillstyle = 'solid'
        hErr['up']['lumi'].drawstyle = "hist"
        hErr['up']['lumi'].color = "orange"
        hErr['up']['lumi'].legendstyle = 'F'

        recoByChan[chan].applyWeight(nominalWeight, True)
        recoByChan[chan].applyWeight(str(1.-lumiUnc))
        true[chan].applyWeight(str(1.-lumiUnc), True)
        res = getResponse(chan, true[chan], recoByChan[chan], bkg[chan], var,
                          binning[varName], fPUWt[''])
        unf = RooUnfoldIter(res, hData, nIter)
        hUnf = asrootpy(unf.Hreco())
        hErr['dn']['lumi'] = hUnf - hUnfolded
        hErr['dn']['lumi'].title = "Luminosity"
        hErr['dn']['lumi'].fillstyle = 'solid'
        hErr['dn']['lumi'].drawstyle = "hist"
        hErr['dn']['lumi'].color = "orange"
        hErr['dn']['lumi'].legendstyle = 'F'

        recoByChan[chan].applyWeight(nominalWeight, True)
        true[chan].applyWeight('',True)

        # Make all error histograms positive (only matters for error plot)
        for sys in hErr.values():
            for h in sys.values():
                for b in h:
                    b.value = abs(b.value)

        # Make plots of uncertainties (added linearly)
        cErrUp = Canvas(1000,1000)
        errStackUp = HistStack(hErr['up'].values(), drawstyle = 'histnoclear')
        draw(errStackUp, cErrUp, xtitle=xTitle[varName], ytitle="Error (+)",yerror_in_padding=False)
        leg = makeLegend(cErrUp, *hErr['up'].values())
        leg.Draw('same')
        style.setCMSStyle(cErrUp, '', dataType='Preliminary Simulation', intLumi=lumi)
        cErrUp.Print(_join(plotDir, 'errUp_{}_{}.png'.format(varName, chan)))

        cErrDn = Canvas(1000,1000)
        errStackDn = HistStack(hErr['dn'].values(), drawstyle = 'histnoclear')
        draw(errStackDn, cErrDn, xtitle=xTitle[varName], ytitle="Error (-)",yerror_in_padding=False)
        leg = makeLegend(cErrDn, *hErr['dn'].values())
        leg.Draw('same')
        style.setCMSStyle(cErrDn, '', dataType='Preliminary Simulation', intLumi=lumi)
        cErrDn.Print(_join(plotDir, 'errDown_{}_{}.png'.format(varName, chan)))


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

            bUncUp.value = sqrt(bUncUp.value)
            bUncDn.value = sqrt(bUncDn.value)
                                                      

        errorBand = makeErrorBand(hUnfolded, hUncUp, hUncDn)
        if 'uncUpSqr' not in hTot:
            hTot['uncUpSqr'] = hUncUp.empty_clone()
        for b, bTot in zip(hUncUp, hTot['uncUpSqr']):
            bTot.value += b.value**2
        if 'uncDnSqr' not in hTot:
            hTot['uncDnSqr'] = hUncDn.empty_clone()
        for b, bTot in zip(hUncDn, hTot['uncDnSqr']):
            bTot.value += b.value**2


        ### plot
        hUnfolded.color = 'black'
        hUnfolded.drawstyle = 'PE'
        hUnfolded.legendstyle = 'LPE'
        hUnfolded.title = 'Unfolded data + stat. unc.'
        if 'unfolded' not in hTot:
            hTot['unfolded'] = hUnfolded.empty_clone()
        hTot['unfolded'] += hUnfolded
        normalizeBins(hUnfolded)

        hUnfoldedAlt.color = 'r'
        hUnfoldedAlt.drawstyle = 'hist'
        hUnfoldedAlt.fillstyle = 'hollow'
        hUnfoldedAlt.legendstyle = 'L'
        hUnfoldedAlt.title = 'Unfolded fullsim aMC@NLO+MCFM'
        if 'unfoldedAlt' not in hTot:
            hTot['unfoldedAlt'] = hUnfoldedAlt.empty_clone()
        hTot['unfoldedAlt'] += hUnfoldedAlt
        normalizeBins(hUnfoldedAlt)

        hTrue = true[chan].makeHist(var, '', binning[varName])
        hTrue.color = 'blue'
        hTrue.drawstyle = 'hist'
        hTrue.fillstyle = 'hollow'
        hTrue.legendstyle = 'L'
        hTrue.title = 'POWHEG (training sample)'

        hTrueAlt = altTrue[chan].makeHist(var, '', binning[varName])
        hTrueAlt.color = 'magenta'
        hTrueAlt.drawstyle = 'hist'
        hTrueAlt.fillstyle = 'hollow'
        hTrueAlt.legendstyle = 'L'
        hTrueAlt.title = 'aMC@NLO'

        cUnf = Canvas(1000,1000)
        draw([hTrue, hTrueAlt, hUnfoldedAlt, hUnfolded, errorBand], cUnf, 
             xtitle=xTitle[varName], ytitle='Events',yerror_in_padding=False)
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
    normalizeBins(hTot['unfolded'])

    hTot['unfoldedAlt'].color = 'r'
    hTot['unfoldedAlt'].drawstyle = 'hist'
    hTot['unfoldedAlt'].fillstyle = 'hollow'
    hTot['unfoldedAlt'].legendstyle = 'L'
    hTot['unfoldedAlt'].title = 'Unfolded fullsim aMC@NLO+MCFM'
    normalizeBins(hTot['unfoldedAlt'])

    hTrue = true.makeHist(var, '', binning[varName])
    hTrue.color = 'blue'
    hTrue.drawstyle = 'hist'
    hTrue.fillstyle = 'hollow'
    hTrue.legendstyle = 'L'
    hTrue.title = 'POWHEG (training sample)'

    hTrueAlt = altTrue.makeHist(var, '', binning[varName])
    hTrueAlt.color = 'magenta'
    hTrueAlt.drawstyle = 'hist'
    hTrueAlt.fillstyle = 'hollow'
    hTrueAlt.legendstyle = 'L'
    hTrueAlt.title = 'aMC@NLO'

    hUncUp = hTot['uncUpSqr'].clone()
    for b in hUncUp:
        b.value = sqrt(b.value)
    hUncDn = hTot['uncDnSqr'].clone()
    for b in hUncDn:
        b.value = sqrt(b.value)

    errorBand = makeErrorBand(hTot['unfolded'], hUncUp, hUncDn)

    cUnf = Canvas(1000,1000)
    draw([hTrue, hTrueAlt, hTot['unfoldedAlt'], hTot['unfolded'], errorBand], 
         cUnf, xtitle=xTitle[varName], ytitle='Events',yerror_in_padding=False)
    leg = makeLegend(cUnf, hTot['unfolded'], errorBand, hTot['unfoldedAlt'], hTrue, 
                     hTrueAlt, textsize=.023, leftmargin=0.3)
    leg.Draw("same")

    style.setCMSStyle(cUnf, '', dataType='Preliminary', intLumi=lumi)
    cUnf.Print(_join(plotDir, "unfold_{}.png".format(varName)))


