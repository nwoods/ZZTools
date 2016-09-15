
import logging
from rootpy import log as rlog; rlog = rlog["/unfold"]
# don't show most silly ROOT messages
logging.basicConfig(level=logging.WARNING)
rlog["/ROOT.TUnixSystem.SetDisplay"].setLevel(rlog.ERROR)

from rootpy import asrootpy
from rootpy.io import root_open
from rootpy.plotting import Canvas, Legend, Hist, Hist2D
from rootpy.plotting.utils import draw
from rootpy.ROOT import RooUnfoldResponse, cout, TDecompSVD
from rootpy.ROOT import RooUnfoldBayes as RooUnfoldIter # it's frequentist!

from SampleTools import MCSample, DataSample, SampleGroup, SampleStack
from PlotTools import PlotStyle as _Style
from PlotTools import makeLegend, addPadBelow, makeRatio, fixRatioAxes
from Utilities import WeightStringMaker
from Analysis.setupStandardSamples import standardZZSamples, genZZSamples
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

channels = ['eeee','eemm', 'mmmm']

data, stack = standardZZSamples(inData, inMC, 'smp', puWeightFile, 
                                fakeRateFile, lumi)
true = genZZSamples(inMC, 'smp', lumi)

fakeDataFile = _join('/data/nawoods/ntuples', inMC, 'results_smp', 
                     'ZZTo4L-amcatnlo*.root')
fakeData = SampleGroup('ZZTo4L-amcatnlo', 'zz', {
        c : MCSample('ZZTo4L-amcatnlo', c, 
                     fakeDataFile, True, lumi) for c in channels
        }, True)
for c in channels:
    fakeData[c].applyWeight(stack[2][c]._weight)

fakeDataTruth = SampleGroup('ZZTo4L-amcatnlo', 'zzGen', {
        c : MCSample('ZZTo4L-amcatnlo', c+'Gen', 
                     fakeDataFile, True, lumi) for c in channels
        }, True)

# the PU weight function was created during standardZZSamples(); get it.
from rootpy.ROOT import puWeight0 as puWt

binning = [20.*i for i in range(4)] + [100., 140., 200., 300.] #[15, 0., 300.]
response = {}
unfolders = {}
hUnfolded = {}

def normalizeBins(h):
    binUnit = min(h.GetBinWidth(b) for b in range(1,len(h)+1))
    for ib in xrange(1,len(h)+1):
        w = h.GetBinWidth(ib)
        h.SetBinContent(ib, h.GetBinContent(ib) * binUnit / w)
        h.SetBinError(ib, h.GetBinError(ib) * sqrt(binUnit / w))
        if h.GetBinError(ib) > h.GetBinContent(ib):
            h.SetBinError(ib, h.GetBinContent(ib))
    h.sumw2()


hTot = {}
genPt = {}
for c in channels:
    chanGenPt = {}

    for name, sample in true[c].itersamples():
        sampleGenPt = {}

        for row in sample:
            sampleGenPt[(row.run,row.lumi,row.evt)] = row.Pt

        chanGenPt[name] = sampleGenPt

    genPt[c] = chanGenPt

    hReco = sum(stack.makeHist({c:'Pt'}, '', binning, perUnitWidth=False).hists)

    if len(binning) == 3:
        hResponse = Hist2D(*(binning+binning))
    else:
        hResponse = Hist2D(binning, binning)

    stackMC = {}
    for sample in stack:
        if sample.name in ['Z+X', 'bkg']:
            hBkg = sample[c].makeHist('Pt', '', binning, perUnitWidth=False)
            continue
        s = sample[c]
        if isinstance(s, SampleGroup):
            for name, ss in s.itersamples():
                stackMC[ss.name] = ss
        else:
            stackMC[s.name] = s

    for name, sample in stackMC.iteritems():
        xsec = sample.xsec
        sumW = sample.sumW
        for row in sample:
            evtID = (row.run, row.lumi, row.evt)
            weight = puWt(row.nTruePU) * row.genWeight * xsec * lumi / sumW
            try:
                hResponse.Fill(row.Pt, chanGenPt[name][evtID], weight)
            except KeyError:
                pass

    hTrue = true[c].makeHist('Pt', '', binning, perUnitWidth=False)

    response[c] = RooUnfoldResponse(hReco, hTrue, hResponse)

    svd = TDecompSVD(response[c].Mresponse())
    sig = svd.GetSig()

    hData = data[c].makeHist('Pt', '', binning, perUnitWidth=False)
    hData -= hBkg
    unfolders[c] = RooUnfoldIter(response[c], hData, 8)

    #print c
    #unfolders[c].PrintTable(cout)

    hUnfolded[c] = asrootpy(unfolders[c].Hreco())
    hUnfolded[c].color = 'black'
    hUnfolded[c].drawstyle = 'PE'
    hUnfolded[c].legendstyle = 'LPE'
    hUnfolded[c].title = 'Unfolded data'
    if 'unfolded' not in hTot:
        hTot['unfolded'] = hUnfolded[c].empty_clone()
    hTot['unfolded'] += hUnfolded[c]
    normalizeBins(hUnfolded[c])

    hFake = fakeData[c].makeHist('Pt', '', binning, perUnitWidth=False)
    fakeUnfolder = RooUnfoldIter(response[c], hFake, 8)
    hUnfoldedFake = asrootpy(fakeUnfolder.Hreco())
    hUnfoldedFake.color = 'red'
    hUnfoldedFake.title = 'Unfolded aMC@NLO'
    hUnfoldedFake.drawstyle = 'PE'
    hUnfoldedFake.legendstyle = 'LPE'
    if 'unfoldedFake' not in hTot:
        hTot['unfoldedFake'] = hUnfoldedFake.empty_clone()
    hTot['unfoldedFake'] += hUnfoldedFake
    normalizeBins(hUnfoldedFake)

    hTrue.title = 'POWHEG + MCFM'
    hTrue.drawstyle = 'hist'
    hTrue.color = 'blue'
    hTrue.legendstyle = 'L'
    hTrue.fillstyle = 'hollow'
    if 'true' not in hTot:
        hTot['true'] = hTrue.empty_clone()
    hTot['true'] += hTrue
    normalizeBins(hTrue)

    hFakeTruth = fakeDataTruth[c].makeHist('Pt', '', binning, 
                                           perUnitWidth=False)
    hFakeTruth.title = 'aMC@NLO'
    hFakeTruth.drawstyle = 'hist'
    hFakeTruth.color = 'orange'
    hFakeTruth.legendstyle = 'L'
    hFakeTruth.fillstyle = 'hollow'
    if 'fakeTruth' not in hTot:
        hTot['fakeTruth'] = hFakeTruth.empty_clone()
    hTot['fakeTruth'] += hFakeTruth
    normalizeBins(hFakeTruth)

    cUnf = Canvas(1000,1200)

    draw([hTrue, hFakeTruth, hUnfoldedFake, hUnfolded[c]], cUnf, xtitle='4\\ell p_T', ytitle='Events / 20 GeV')
    leg = makeLegend(cUnf, hUnfoldedFake, hFakeTruth, hUnfolded[c], hTrue, 
                     textsize=.023)
    leg.Draw("same")

    style.setCMSStyle(cUnf, '', dataType='Preliminary', intLumi=lumi)
    cUnf.Print(_join(plotDir, "unfold_Pt_{}.png".format(c)))

    cRes = Canvas(1000,1000)
    hRes = asrootpy(response[c].Hresponse())
    hRes.drawstyle = 'colztext'
    hRes.draw()
    style.setCMSStyle(cRes, '', dataType='Preliminary Simulation', intLumi=lumi)
    cRes.Print(_join(plotDir, "response_Pt_{}.png".format(c)))


cUnf = Canvas(1000,1200)

hTot['unfolded'].color = 'black'
hTot['unfolded'].drawstyle = 'PE'
hTot['unfolded'].legendstyle = 'LPE'
hTot['unfolded'].title = 'Unfolded data'

hTot['unfoldedFake'].color = 'red'
hTot['unfoldedFake'].title = 'Unfolded aMC@NLO'
hTot['unfoldedFake'].drawstyle = 'PE'
hTot['unfoldedFake'].legendstyle = 'LPE'
 
hTot['true'].title = 'POWHEG + MCFM'
hTot['true'].drawstyle = 'hist'
hTot['true'].color = 'blue'
hTot['true'].legendstyle = 'L'
hTot['true'].fillstyle = 'hollow'

hTot['fakeTruth'].title = 'aMC@NLO'
hTot['fakeTruth'].drawstyle = 'hist'
hTot['fakeTruth'].color = 'orange'
hTot['fakeTruth'].legendstyle = 'L'
hTot['fakeTruth'].fillstyle = 'hollow'

for n,h in hTot.iteritems():
    normalizeBins(h)

draw([hTot['true'], hTot['fakeTruth'], hTot['unfoldedFake'], 
      hTot['unfolded']], cUnf, xtitle='4\\ell p_T', ytitle='Events / 20 GeV')
leg = makeLegend(cUnf, hTot['unfoldedFake'], hTot['fakeTruth'], 
                 hTot['unfolded'], hTot['true'], 
                 textsize=.023)
leg.Draw("same")
    
style.setCMSStyle(cUnf, '', dataType='Preliminary', intLumi=lumi)
cUnf.Print(_join(plotDir, "unfold_Pt.png"))
