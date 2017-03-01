
import logging
from rootpy import log as rlog; rlog = rlog["/FakeRateCalculator"]
# don't show most silly ROOT messages
logging.basicConfig(level=logging.WARNING)
rlog["/ROOT.TUnixSystem.SetDisplay"].setLevel(rlog.ERROR)

from SampleTools import MCSample, DataSample, SampleGroup, SampleStack
from PlotTools import PlotStyle as _Style
from PlotTools import makeLegend
from Utilities import WeightStringMaker

from rootpy import asrootpy
from rootpy.io import root_open
from rootpy.plotting import Canvas
from rootpy.plotting.utils import draw

from os import path as _path
from os import makedirs as _mkdir
from os import environ



def calculateFakeRate(sampleID, outFile, puFile, lumi, plot=True,
                      plotDir='/afs/cern.ch/user/n/nawoods/www/UWVVPlots/fakeRate'):
    '''
    Calculate fake rates/factors and put them in a root file, stored in
    histograms binned in (abs(eta),pt).
    Note: fake factor = 1/(1-fake rate)

    sampleID (str): e.g. '08sep2016'
    outFile (str): Location for output fake rate/fake factor histograms,
        relative to ZZTools/data/fakeRate
    puFile (str): Location of PU weights, relative to ZZTools/data/pileup
    lumi (float): integrated luminosity of data sample
    plot (bool): if True, draw plots
    plotDir (str): absolute path of directory to put plots in, if applicable
    '''
    if plot:
        style = _Style()
        if not _path.exists(plotDir):
            _mkdir(plotDir)
        elif not _path.isdir(plotDir):
            raise IOError("There is already some non-directory object called {}.".format(plotDir))


    fStr = '/data/nawoods/ntuples/uwvvZPlusl_{{}}_{}/results{{}}/{{}}*.root'.format(sampleID)
    outFile = _path.join(_path.join(environ['zzt'], 'data', 'fakeRate'), outFile)

    channels = ['eee','eem','emm','mmm']

    puWeight = WeightStringMaker('puWeight')
    fPU = root_open(_path.join(environ['zzt'], 'data', 'pileup',
                               puFile))
    hPU = fPU.puScaleFactor
    strPU = puWeight.makeWeightStringFromHist(hPU, 'nTruePU')


    mcWeight = {
        'eee' : 'e1EffScaleFactor * e2EffScaleFactor * e3EffScaleFactor',
        'eem' : 'e1EffScaleFactor * e2EffScaleFactor * mEffScaleFactor',
        'emm' : 'eEffScaleFactor * m1EffScaleFactor * m2EffScaleFactor',
        'mmm' : 'm1EffScaleFactor * m2EffScaleFactor * m3EffScaleFactor',
    }

    mcSamples = ['DYToLL-0J', 'DYToLL-1J', 'DYToLL-2J', 'TTTo2L2Nu']
    signalSamples = ['WZTo3LNu', 'ZZTo4L',
                     'GluGluZZTo4e','GluGluZZTo4mu','GluGluZZTo2e2mu']

    mcLooseByChan = {}
    mcTightByChan = {}
    for s in mcSamples:
        mcLooseByChan[s] = {}
        mcTightByChan[s] = {}
        for c in channels:
            mcLooseByChan[s][c] = MCSample(s, c, fStr.format('mc','Loose',s),
                                           True, lumi)
            mcTightByChan[s][c] = MCSample(s, c, fStr.format('mc','Tight',s),
                                           True, lumi)

    mcLoose = {s:SampleGroup(s,'zl',mcLooseByChan[s],True) for s in mcSamples}
    mcStackLoose = SampleStack('mcLoose', 'zl', list(mcLoose.values()))
    mcTight = {s:SampleGroup(s,'zl',mcTightByChan[s],True) for s in mcSamples}
    mcStackTight = SampleStack('mcTight', 'zl', list(mcTight.values()))

    signalLooseByChan = {}
    signalTightByChan = {}
    for s in signalSamples:
        signalLooseByChan[s] = {}
        signalTightByChan[s] = {}
        for c in channels:
            signalLooseByChan[s][c] = MCSample(s, c, fStr.format('mc','Loose',s),
                                               True, lumi)
            signalTightByChan[s][c] = MCSample(s, c, fStr.format('mc','Tight',s),
                                               True, lumi)

    signalLoose = {s:SampleGroup(s,'zl',
                                 signalLooseByChan[s],True)
                   for s in signalSamples}
    signalStackLoose = SampleStack('signalLoose', 'zl', list(signalLoose.values()))
    signalTight = {s:SampleGroup(s,'zl',
                                 signalTightByChan[s],True)
                   for s in signalSamples}
    signalStackTight = SampleStack('signalTight', 'zl', list(signalTight.values()))

    dataLooseByChan = {}
    dataTightByChan = {}
    for c in channels:
        samplesByEraLoose = {}
        samplesByEraTight = {}
        for era in ['2016'+let for let in 'BCDEFGH']:
            samplesByEraLoose[era] = DataSample('data{}'.format(era), c,
                                                fStr.format('data', 'Loose',
                                                            'Run{}'.format(era)))
            samplesByEraTight[era] = DataSample('data{}'.format(era), c,
                                                fStr.format('data', 'Tight',
                                                            'Run{}'.format(era)))
        dataLooseByChan[c] = SampleGroup('dataLoose', c, samplesByEraLoose)
        dataTightByChan[c] = SampleGroup('dataTight', c, samplesByEraTight)

    dataLoose = SampleGroup('data', 'zl', dataLooseByChan)
    dataTight = SampleGroup('data', 'zl', dataTightByChan)

    dataLoose.format(color='black',drawstyle='PE',legendstyle='LPE')
    dataTight.format(color='black',drawstyle='PE',legendstyle='LPE')

    ptBinning = {
        'e' : [7.,]+[8.+2.*i for i in range(6)]+[20.+5.*i for i in range(6)]+[50.+10.*i for i in range(4)]+[200.,],
        'm' : [5.+1.*i for i in range(5)]+[10.,20.,30.,70.,200.],#[10.+5.*i for i in range(4)]+[30.,50.,60.,80.,200.],
        }
    etaBinning = {
        'e' : [0.,0.8,1.47,2.,2.5],
        'm' : [0.,1.2,2.4,2.5],
        }

    ptVars = {
        'e' : {
            'eee' : 'e3Pt',
            'emm' : 'ePt',
            },
        'm' : {
            'eem' : 'mPt',
            'mmm' : 'm3Pt',
            },
        }
    etaVars = {
        'e' : {
            'eee' : 'abs(e3Eta)',
            'emm' : 'abs(eEta)',
            },
        'm' : {
            'eem' : 'abs(mEta)',
            'mmm' : 'abs(m3Eta)',
            },
        }

    writeOut = []
    for lep in ptVars:
        # MC
        mcNumStack = mcStackTight.makeHist2(etaVars[lep], ptVars[lep], '',
                                            etaBinning[lep], ptBinning[lep],
                                            mcWeight)
        mcNum = asrootpy(mcNumStack.GetStack().Last())

        mcDenomStack = mcStackLoose.makeHist2(etaVars[lep], ptVars[lep], '',
                                              etaBinning[lep], ptBinning[lep],
                                              mcWeight)
        mcDenom = asrootpy(mcDenomStack.GetStack().Last())

        fMC = mcNum.clone(name='fakeRateMC_{}'.format(lep))
        fMC.Divide(mcDenom)

        writeOut.append(fMC)

        mcFakeFactor = mcNum.empty_clone(name="fakeFactorMC_{}".format(lep))
        for fb,b in zip(fMC, mcFakeFactor):
            if b.overflow:
                continue
            b.value = fb.value / (1. - fb.value)

        writeOut.append(mcFakeFactor)

        # ZZ/WZ MC to subtract from data
        signalNumStack = signalStackTight.makeHist2(etaVars[lep], ptVars[lep],
                                                    '', etaBinning[lep],
                                                    ptBinning[lep],
                                                    mcWeight)
        signalNum = asrootpy(signalNumStack.GetStack().Last())

        signalDenomStack = signalStackLoose.makeHist2(etaVars[lep],
                                                      ptVars[lep], '',
                                                      etaBinning[lep],
                                                      ptBinning[lep],
                                                      mcWeight)
        signalDenom = asrootpy(signalDenomStack.GetStack().Last())

        # data
        num = dataTight.makeHist2(etaVars[lep], ptVars[lep], '',
                                  etaBinning[lep], ptBinning[lep])
        denom = dataLoose.makeHist2(etaVars[lep], ptVars[lep], '',
                                    etaBinning[lep], ptBinning[lep])

        print "{} Loose: {}".format(lep, denom.GetEntries())
        print "{} Tight: {}".format(lep, num.GetEntries())

        num -= signalNum
        denom -= signalDenom

        for bNum,bDenom in zip(num, denom):
            if bNum.value < 0. or bDenom.value <= 0.:
                bNum.value = 0.
                bNum.error = 0.
                bDenom.value = 0.0000001
                bDenom.error = 0.

        f = num.clone(name='fakeRate_{}'.format(lep))
        f.Divide(denom)

        writeOut.append(f)

        fakeFactor = num.empty_clone(name="fakeFactor_{}".format(lep))
        for fb,b in zip(f, fakeFactor):
            if b.overflow:
                continue
            b.value = fb.value / (1. - fb.value)

        writeOut.append(fakeFactor)


        if plot:
            # MC fake rate 2D
            c2DMC = Canvas(1200,1000)
            c2DMC.SetLogy()
            fMC.xaxis.title = '|#eta|'
            fMC.yaxis.title = 'p_{T}'
            fMC.drawstyle = 'colztexte'
            fMC.draw()
            style.setCMSStyle(c2DMC, '', dataType='Preliminary', intLumi=lumi)
            c2DMC.Print('{}/{}fakeRateMC.png'.format(plotDir, lep))

            # data fake rate 2D
            c2D = Canvas(1200,1000)
            c2D.SetLogy()
            f.xaxis.title = '|#eta|'
            f.yaxis.title = 'p_{T}'
            f.drawstyle = 'colztexte'
            f.draw()
            style.setCMSStyle(c2D, '', dataType='Preliminary', intLumi=lumi)
            c2D.Print('{}/{}fakeRate.png'.format(plotDir, lep))

            # denominator vs pt
            cPtDenom = Canvas(1000,1000)
            ptMC = mcStackLoose.makeHist(ptVars[lep], '', ptBinning[lep],
                                         mcWeight)
            ptMCTot = asrootpy(ptMC.GetStack().Last()).clone() # for ratio
            ptSig = signalStackLoose.makeHist(ptVars[lep], '', ptBinning[lep],
                                              mcWeight)
            ptSigTot = asrootpy(ptSig.GetStack().Last()).clone() # for ratio
            ptMC += ptSig

            ptData = dataLoose.makeHist(ptVars[lep], '', ptBinning[lep],
                                        poissonErrors=True)

            ptDataTot = dataLoose.makeHist(ptVars[lep], '', ptBinning[lep]) # for ratio

            legPtDenom = makeLegend(cPtDenom, ptMC, ptData)

            draw([ptMC, ptData], cPtDenom, xtitle='p_{T} (GeV)', ytitle='Leptons')
            legPtDenom.Draw("same")
            style.setCMSStyle(cPtDenom, '', dataType='Preliminary', intLumi=lumi)
            cPtDenom.Print('{}/{}PtLoose.png'.format(plotDir, lep))

            # numerator vs pt
            cPtNum = Canvas(1000,1000)
            ptMCTight = mcStackTight.makeHist(ptVars[lep], '', ptBinning[lep],
                                              mcWeight)
            ptMCTotTight = asrootpy(ptMCTight.GetStack().Last()).clone() # for ratio
            ptSigTight = signalStackTight.makeHist(ptVars[lep], '',
                                                   ptBinning[lep], mcWeight)
            ptSigTotTight = asrootpy(ptSigTight.GetStack().Last()).clone() # for ratio
            ptMCTight += ptSigTight

            ptDataTight = dataTight.makeHist(ptVars[lep], '', ptBinning[lep],
                                             poissonErrors=True)
            ptDataTotTight = dataTight.makeHist(ptVars[lep], '', # for ratio
                                                ptBinning[lep])

            legPtNum = makeLegend(cPtNum, ptMCTight, ptDataTight)

            draw([ptMCTight, ptDataTight], cPtNum, xtitle='p_{T} (GeV)',
                 ytitle='Tight+Iso Leptons')
            legPtNum.Draw("same")
            style.setCMSStyle(cPtNum, '', dataType='Preliminary', intLumi=lumi)
            cPtNum.Print('{}/{}PtTight.png'.format(plotDir, lep))

            # fake rate vs pt
            cPt = Canvas(1000,1000)
            fPtMC = ptMCTotTight.clone(name='MC', title='MC')
            fPtMC.Divide(ptMCTot)
            fPtMC.color = 'r'
            fPtMC.drawstyle = 'LE'

            ptDataTotTight -= ptSigTotTight
            ptDataTot -= ptSigTot

            for bNum,bDenom in zip(ptDataTotTight, ptDataTot):
                if bNum.value < 0. or bDenom.value <= 0.:
                    bNum.value = 0.
                    bNum.error = 0.
                    bDenom.value = 0.0000001
                    bDenom.error = 0.

            fPt = ptDataTotTight.clone(name='Data')
            fPt.Divide(ptDataTot)
            fPt.drawstyle = 'PE'

            fPtMC.legendstyle = 'L'
            fPt.legendstyle = 'LPE'
            legPt = makeLegend(cPt, fPtMC, fPt)

            draw([fPtMC, fPt], cPt, xtitle='p_{T} (GeV)', ytitle='Fake Rate',
                 ylimits=(0.,.25))
            legPt.Draw("same")
            style.setCMSStyle(cPt, '', dataType='Preliminary', intLumi=lumi)
            cPt.Print('{}/{}FakeRate_pt.png'.format(plotDir, lep))

            # denominator vs eta
            cEtaDenom = Canvas(1000,1000)
            etaMC = mcStackLoose.makeHist(etaVars[lep], '', etaBinning[lep],
                                          mcWeight)
            etaMCTot = asrootpy(etaMC.GetStack().Last()).clone() # for ratio
            etaSig = signalStackLoose.makeHist(etaVars[lep], '',
                                               etaBinning[lep], mcWeight)
            etaSigTot = asrootpy(etaSig.GetStack().Last()).clone() # for ratio
            etaMC += etaSig

            etaData = dataLoose.makeHist(etaVars[lep], '', etaBinning[lep],
                                        poissonErrors=True)
            etaDataTot = dataLoose.makeHist(etaVars[lep], '', etaBinning[lep]) # for ratio

            legEtaDenom = makeLegend(cEtaDenom, etaMC, etaData)

            draw([etaMC, etaData], cEtaDenom, xtitle='p_{T} (GeV)',
                 ytitle='Leptons')
            legEtaDenom.Draw("same")
            style.setCMSStyle(cEtaDenom, '', dataType='Preliminary',
                              intLumi=lumi)
            cEtaDenom.Print('{}/{}EtaLoose.png'.format(plotDir, lep))

            # numerator vs eta
            cEtaNum = Canvas(1000,1000)
            etaMCTight = mcStackTight.makeHist(etaVars[lep], '',
                                               etaBinning[lep], mcWeight)
            etaMCTotTight = asrootpy(etaMCTight.GetStack().Last()).clone() # for ratio
            etaSigTight = signalStackTight.makeHist(etaVars[lep], '',
                                                    etaBinning[lep], mcWeight)
            etaSigTotTight = asrootpy(etaSigTight.GetStack().Last()).clone() # for ratio
            etaMCTight += etaSigTight

            etaDataTight = dataTight.makeHist(etaVars[lep], '',
                                              etaBinning[lep],
                                              poissonErrors=True)
            etaDataTotTight = dataTight.makeHist(etaVars[lep], '', # for ratio
                                                 etaBinning[lep])

            legEtaNum = makeLegend(cEtaNum, etaMCTight, etaDataTight)

            draw([etaMCTight, etaDataTight], cEtaNum, xtitle='p_{T} (GeV)',
                 ytitle='Tight+Iso Leptons')
            legEtaNum.Draw("same")
            style.setCMSStyle(cEtaNum, '', dataType='Preliminary', intLumi=lumi)
            cEtaNum.Print('{}/{}EtaTight.png'.format(plotDir, lep))

            # fake rate vs eta
            cEta = Canvas(1000,1000)
            fEtaMC = etaMCTotTight.clone(name='MC', title='MC')
            fEtaMC.Divide(etaMCTot)
            fEtaMC.color = 'r'
            fEtaMC.drawstyle = 'LE'

            etaDataTotTight -= etaSigTotTight
            etaDataTot -= etaSigTot

            for bNum,bDenom in zip(etaDataTotTight, etaDataTot):
                if bNum.value < 0. or bDenom.value <= 0.:
                    bNum.value = 0.
                    bNum.error = 0.
                    bDenom.value = 0.0000001
                    bDenom.error = 0.

            fEta = etaDataTotTight.clone(name='Data', title='Data')
            fEta.Divide(etaDataTot)
            fEta.drawstyle = 'PE'

            fEtaMC.legendstyle = 'L'
            fEta.legendstyle = 'LPE'
            legEta = makeLegend(cEta, fEtaMC, fEta)

            draw([fEtaMC, fEta], cEta, xtitle='p_{T} (GeV)', ytitle='Fake Rate',
                 ylimits=(0.,.25))
            legEta.Draw("same")
            style.setCMSStyle(cEta, '', dataType='Preliminary', intLumi=lumi)
            cEta.Print('{}/{}FakeRate_eta.png'.format(plotDir, lep))


    with root_open(outFile, 'recreate') as fOut:
        for h in writeOut:
            h.write()


if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Script for calculating fake rates')
    parser.add_argument('sampleID', nargs=1, type=str,
                        help='ID string for samples, e.g. "08sep2016".')
    parser.add_argument('outFile', nargs='?', type=str, default='fakeRate.root',
                        help=('File name to store output histograms in, '
                              'relative to ZZTools/data/fakeRate.'))
    parser.add_argument('puFile', nargs='?', type=str,
                        default='puWeight_69200_08sep2016.root',
                        help=('File name for pileup weights, '
                              'relative to ZZTools/data/pileup.'))
    parser.add_argument('lumi', nargs='?', type=float, default=15937.,
                        help='Integrated luminosity of sample, in pb^-1.')
    parser.add_argument('--plot', action='store_true',
                        help='Make plots of fake rates and inputs.')
    parser.add_argument('--plotDir', nargs='?', type=str,
                        default='/afs/cern.ch/user/n/nawoods/www/UWVVPlots/fakeRate',
                        help='Location for plots, if applicable.')

    args = parser.parse_args()

    calculateFakeRate(args.sampleID[0], args.outFile, args.puFile, args.lumi,
                      args.plot, args.plotDir)

