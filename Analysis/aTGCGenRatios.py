import logging
from rootpy import log as rlog; rlog = rlog["/aTGCGenRatios"]
# don't show most silly ROOT messages
logging.basicConfig(level=logging.WARNING)
rlog["/ROOT.TUnixSystem.SetDisplay"].setLevel(rlog.ERROR)

from rootpy import asrootpy
from rootpy.io import root_open
from rootpy.plotting import Graph, Graph2D, Canvas
from rootpy.ROOT import TF2

from Utilities import WeightStringMaker, identityFunction
from Analysis.setupStandardSamples import *
from SampleTools import SampleGroup
from Analysis.weightHelpers import puWeight, baseMCWeight

from os.path import join, exists
from os import makedirs


var = 'Mass'
binning = [100.+100.*i for i in range(6)] + [800.,1000.,1200.]
aTGCParams = ['f4','f5']

channels = ['eeee','eemm','mmmm']

fileTemp = 'histo_{fg}I{fz}_{param}_file.root'

fgs = [-.0038, 0., .0038]
fzs = [-0.003, 0., 0.003]

paramStrs = {
    -.0038 : 'm0p0038',
    -0.003 : 'm0p003',
    0.     : '0',
    0.003  : '0p003',
    .0038  : '0p0038',
    }

fileNames = {
    p:{(fg,fz):fileTemp.format(fg=paramStrs[fg],fz=paramStrs[fz],param=p)
       for fg in fgs
       for fz in fzs}
    for p in aTGCParams
    }

systDirNames = {
    'ees_up' : 'eScaleUp',
    'ees_dn' : 'eScaleDn',
    'eerPhi_up' : 'ePhiResUp',
    'eerRho_up' : 'eRhoResUp',
    'eerRho_dn' : 'eRhoResDn',
    'mClosure_up' : 'mClosureUp',
    'mClosure_dn' : 'mClosureDn',
    }

def main(inData, inMC, inATGC, outDir, fakeRateFile, puWeightFile, lumi,
         lumiATGC=40000):

    true = standardZZGen('zz', inMC, 'ZZTo4L', 'smp', lumiATGC)

    hTrue = true.makeHist(var, '', binning, perUnitWidth=False)
    cTrue = Canvas(1000,1000)
    hTrue.draw('hist')
    cTrue.Print('~/www/UWVVPlots/trueYield.png')
    # merge overflow
    hTrue = hTrue.merge_bins([(-2,-1)])

    # get sherpa/powheg ratio to use for reweighting
    qqWeights = {}
    for param in aTGCParams:
        qqWeights[param] = {}
        hSM = hTrue.empty_clone()
        with root_open(join(inATGC, fileTemp.format(fg='0',fz='0',param=param))) as f:
            h = asrootpy(f.h_ratio_MZZ_wt)
            hSM += h
        wtMaker = WeightStringMaker(param)
        for iWP, ((fg,fz), fileName) in enumerate(fileNames[param].iteritems()):
            ratio = hTrue.empty_clone()
            with root_open(join(inATGC,fileName)) as f:
                hATGC = asrootpy(f.h_ratio_MZZ_wt)
                ratio += hATGC
            ratio /= hSM
            ratio[-1].value = ratio[-2].value
            cRatio = Canvas(1000,1000)
            ratio.Draw('hist')
            cRatio.Print("~/www/UWVVPlots/ratio_{}_{}_{}.png".format(param, fg, fz))
            qqWeights[param][(fg,fz)] = wtMaker.makeWeightStringFromHist(ratio, 'Mass')

    nominalWeight = baseMCWeight('zz', puWeightFile)

    # prepare samples
    qqZZ = standardZZMC('zz', inMC, 'ZZTo4L', 'smp', puWeightFile, lumi)

    # gg->ZZ
    ggZZByChan = {}
    for c in channels:
        ggZZByFS = {
            fs : standardZZMC(c, inMC, 'GluGluZZTo{}'.format(fs),
                              'smp', puWeightFile, lumi)
            for fs in ['4e', '4mu', '2e2mu']
            }
        ggZZByChan[c] = SampleGroup('GluGluZZ', c, ggZZByFS, True)
    ggZZ = SampleGroup('GluGluZZ', 'zz', ggZZByChan, True)

    bkgMC = zzIrreducibleBkg('zz', inMC, 'smp', puWeightFile, lumi)

    qqZZSyst = {}
    ggZZSyst = {}
    bkgMCSyst = {}
    for s in ['ees_up', 'ees_dn', 'eerRho_up',
              'eerRho_dn', 'eerPhi_up']:
        syst = systDirNames[s]
        qqZZSyst[s] = standardZZMC('eeee,eemm',
                                   inMC.replace('mc_','mc_{}_'.format(syst)),
                                   'ZZTo4L', 'smp', puWeightFile, lumi)
        ggZZSystByChan = {}
        for c in ['eeee','eemm']:
            ggZZSystByFS = {
                fs : standardZZMC(c, inMC.replace('mc_','mc_{}_'.format(syst)),
                                  'GluGluZZTo{}'.format(fs),
                                  'smp', puWeightFile, lumi)
                for fs in ['4e', '2e2mu']
                }
            ggZZSystByChan[c] = SampleGroup('GluGluZZ', c, ggZZSystByFS, True)
        ggZZSyst[s] = SampleGroup('GluGluZZ', 'eeee,eemm', ggZZSystByChan, True)
        bkgMCSyst[s] = zzIrreducibleBkg('eeee,eemm', inMC.replace('mc_','mc_{}_'.format(syst)),
                                        'smp', puWeightFile, lumi)
    for s in ['mClosure_up','mClosure_dn']:
        syst = systDirNames[s]
        qqZZSyst[s] = standardZZMC('eemm,mmmm',
                                   inMC.replace('mc_','mc_{}_'.format(syst)),
                                   'ZZTo4L', 'smp', puWeightFile, lumi)
        ggZZSystByChan = {}
        for c in ['eemm','mmmm']:
            ggZZSystByFS = {
                fs : standardZZMC(c, inMC.replace('mc_','mc_{}_'.format(syst)),
                                  'GluGluZZTo{}'.format(fs),
                                  'smp', puWeightFile, lumi)
                for fs in ['4mu', '2e2mu']
                }
            ggZZSystByChan[c] = SampleGroup('GluGluZZ', c, ggZZSystByFS, True)
        ggZZSyst[s] = SampleGroup('GluGluZZ', 'eemm,mmmm', ggZZSystByChan, True)
        bkgMCSyst[s] = zzIrreducibleBkg('eemm,mmmm', inMC.replace('mc_','mc_{}_'.format(syst)),
                                        'smp', puWeightFile, lumi)

    bkgData = standardZZBkg('zz', inData, inMC, 'smp', puWeightFile,
                            fakeRateFile, lumi)
    bkgDataSyst = {
        'eup' : standardZZBkg('eeee,eemm', inData, inMC, 'smp', puWeightFile,
                              fakeRateFile, lumi, eFakeRateSyst='up'),
        'edn' : standardZZBkg('eeee,eemm', inData, inMC, 'smp', puWeightFile,
                              fakeRateFile, lumi, eFakeRateSyst='dn'),
        'mup' : standardZZBkg('eemm,mmmm', inData, inMC, 'smp', puWeightFile,
                              fakeRateFile, lumi, mFakeRateSyst='up'),
        'mdn' : standardZZBkg('eemm,mmmm', inData, inMC, 'smp', puWeightFile,
                              fakeRateFile, lumi, mFakeRateSyst='dn'),
        }

    # weigh MC
    qqZZ.applyWeight(nominalWeight, True)
    for s in qqZZSyst.values():
        s.applyWeight(nominalWeight, True)
    ggZZ.applyWeight(nominalWeight, True)
    for s in ggZZSyst.values():
        s.applyWeight(nominalWeight, True)
    bkgMC.applyWeight(nominalWeight, True)
    for s in bkgMCSyst.values():
        s.applyWeight(nominalWeight, True)

    # base histograms for reuse
    qqZZNominalSM = qqZZ.makeHist(var, '', binning, perUnitWidth=False,
                                  mergeOverflow=True)

    qqZZHists = {
        '' : {
            p:{fgfz:qqZZ.makeHist(var, '', binning, w, perUnitWidth=False,
                                  mergeOverflow=True)
               for fgfz, w in wts.iteritems()}
            for p, wts in qqWeights.iteritems()
            }
        }

    ggZZHists = {
        '' : ggZZ.makeHist(var, '', binning, perUnitWidth=False,
                           mergeOverflow=True)
        }
    bkgDataHists = {
        '' : bkgData.makeHist(var, '', binning, perUnitWidth=False,
                              mergeOverflow=True)
        }
    bkgMCHists = {
        '' : bkgMC.makeHist(var, '', binning, perUnitWidth=False,
                            mergeOverflow=True)
        }

    sigSM = qqZZNominalSM + ggZZHists['']


    ### Systematics

    # PU weight systematic
    for sys in ['up','dn']:
        puReweight = baseMCWeight('zz', puWeightFile, puSyst=sys)
        qqZZ.applyWeight(puReweight, True)
        ggZZ.applyWeight(puReweight, True)
        bkgMC.applyWeight(puReweight, True)

        qqZZHists['pu_'+sys] = {
            p:{fgfz:qqZZ.makeHist(var, '', binning, w, perUnitWidth=False,
                                  mergeOverflow=True)
               for fgfz, w in wts.iteritems()}
            for p, wts in qqWeights.iteritems()
            }

        ggZZHists['pu_'+sys] = ggZZ.makeHist(var, '', binning, perUnitWidth=False,
                                             mergeOverflow=True)

        bkgMCHists['pu_'+sys] = bkgMC.makeHist(var, '', binning, perUnitWidth=False,
                                               mergeOverflow=True)

        qqZZ.applyWeight(nominalWeight, True)
        ggZZ.applyWeight(nominalWeight, True)
        bkgMC.applyWeight(nominalWeight, True)


    # Lepton efficiency systematics
    for lep in ['e','m']:
        for sys in ['up','dn']:
            wtArg = {lep+'Syst':sys}
            lepEffReweight = baseMCWeight('zz', puWeightFile, **wtArg)
            qqZZ.applyWeight(lepEffReweight, True)
            ggZZ.applyWeight(lepEffReweight, True)
            bkgMC.applyWeight(lepEffReweight, True)

            qqZZHists[lep+'Eff_'+sys] = {
                p:{fgfz:qqZZ.makeHist(var, '', binning, w, perUnitWidth=False,
                                      mergeOverflow=True)
                   for fgfz, w in wts.iteritems()}
                for p, wts in qqWeights.iteritems()
                }

            ggZZHists[lep+'Eff_'+sys] = ggZZ.makeHist(var, '', binning, perUnitWidth=False,
                                                      mergeOverflow=True)

            bkgMCHists[lep+'Eff_'+sys] = bkgMC.makeHist(var, '', binning, perUnitWidth=False,
                                                        mergeOverflow=True)

            qqZZ.applyWeight(nominalWeight, True)
            ggZZ.applyWeight(nominalWeight, True)
            bkgMC.applyWeight(nominalWeight, True)


    # Lepton fake rate systematics
    for sys in ['up','dn']:
        for lep in ['e','m']:
            bkgDataHists[lep+'FakeRate_'+sys] = bkgDataSyst[lep+sys].makeHist(var, '', binning, perUnitWidth=False,
                                                                              mergeOverflow=True)
            if lep == 'e':
                otherChannel = 'mmmm'
            else:
                otherChannel = 'eeee'
            bkgDataHists[lep+'FakeRate_'+sys] += bkgData[otherChannel].makeHist(var, '', binning, perUnitWidth=False,
                                                                                mergeOverflow=True)


    # electron energy scale
    for sys in ['up','dn']:
        qqZZHists['eEnergyScale_'+sys] = {
            p:{fgfz:qqZZSyst['ees_'+sys].makeHist(var, '', binning, w, perUnitWidth=False,
                                                  mergeOverflow=True) + \
                   qqZZ['mmmm'].makeHist(var, '', binning, w, perUnitWidth=False,
                                         mergeOverflow=True)
               for fgfz, w in wts.iteritems()}
            for p, wts in qqWeights.iteritems()
            }

        ggZZHists['eEnergyScale_'+sys] = ggZZSyst['ees_'+sys].makeHist(var, '', binning, perUnitWidth=False,
                                                                       mergeOverflow=True)
        ggZZHists['eEnergyScale_'+sys] += ggZZ['mmmm'].makeHist(var, '', binning, perUnitWidth=False,
                                                                mergeOverflow=True)

        bkgMCHists['eEnergyScale_'+sys] = bkgMCSyst['ees_'+sys].makeHist(var, '', binning, perUnitWidth=False,
                                                                          mergeOverflow=True)
        bkgMCHists['eEnergyScale_'+sys] += bkgMC['mmmm'].makeHist(var, '', binning, perUnitWidth=False,
                                                                  mergeOverflow=True)


    # electron energy resolution
    for sys in ['up','dn']:
        qqZZHists['eEnergyResolutionRho_'+sys] = {
            p:{fgfz:qqZZSyst['eerRho_'+sys].makeHist(var, '', binning, w, perUnitWidth=False,
                                                     mergeOverflow=True) + \
                   qqZZ['mmmm'].makeHist(var, '', binning, w, perUnitWidth=False,
                                         mergeOverflow=True)
               for fgfz, w in wts.iteritems()}
            for p, wts in qqWeights.iteritems()
            }

        ggZZHists['eEnergyResolutionRho_'+sys] = ggZZSyst['eerRho_'+sys].makeHist(var, '', binning, perUnitWidth=False,
                                                                                  mergeOverflow=True)
        ggZZHists['eEnergyResolutionRho_'+sys] += ggZZ['mmmm'].makeHist(var, '', binning, perUnitWidth=False,
                                                                        mergeOverflow=True)

        bkgMCHists['eEnergyResolutionRho_'+sys] = bkgMCSyst['eerRho_'+sys].makeHist(var, '', binning, perUnitWidth=False,
                                                                                    mergeOverflow=True)
        bkgMCHists['eEnergyResolutionRho_'+sys] += bkgMC['mmmm'].makeHist(var, '', binning, perUnitWidth=False,
                                                                          mergeOverflow=True)


    qqZZHists['eEnergyResolutionPhi'] = {
        p:{fgfz:qqZZSyst['eerPhi_up'].makeHist(var, '', binning, w, perUnitWidth=False,
                                               mergeOverflow=True) + \
               qqZZ['mmmm'].makeHist(var, '', binning, w, perUnitWidth=False,
                                     mergeOverflow=True)
           for fgfz, w in wts.iteritems()}
        for p, wts in qqWeights.iteritems()
        }

    ggZZHists['eEnergyResolutionPhi'] = ggZZSyst['eerPhi_up'].makeHist(var, '', binning, perUnitWidth=False,
                                                                       mergeOverflow=True)
    ggZZHists['eEnergyResolutionPhi'] += ggZZ['mmmm'].makeHist(var, '', binning, perUnitWidth=False,
                                                               mergeOverflow=True)

    bkgMCHists['eEnergyResolutionPhi'] = bkgMCSyst['eerPhi_up'].makeHist(var, '', binning, perUnitWidth=False,
                                                                         mergeOverflow=True)
    bkgMCHists['eEnergyResolutionPhi'] += bkgMC['mmmm'].makeHist(var, '', binning, perUnitWidth=False,
                                                                 mergeOverflow=True)


    # muon energy scale/resolution
    for sys in ['up','dn']:
        qqZZHists['mEnergy_'+sys] = {
            p:{fgfz:qqZZSyst['mClosure_'+sys].makeHist(var, '', binning, w, perUnitWidth=False,
                                                       mergeOverflow=True) + \
                   qqZZ['eeee'].makeHist(var, '', binning, w, perUnitWidth=False,
                                         mergeOverflow=True)
               for fgfz, w in wts.iteritems()}
            for p, wts in qqWeights.iteritems()
            }

        ggZZHists['mEnergy_'+sys] = ggZZSyst['mClosure_'+sys].makeHist(var, '', binning, perUnitWidth=False,
                                                                       mergeOverflow=True)
        ggZZHists['mEnergy_'+sys] += ggZZ['eeee'].makeHist(var, '', binning, perUnitWidth=False,
                                                           mergeOverflow=True)

        bkgMCHists['mEnergy_'+sys] = bkgMCSyst['mClosure_'+sys].makeHist(var, '', binning, perUnitWidth=False,
                                                                         mergeOverflow=True)
        bkgMCHists['mEnergy_'+sys] += bkgMC['eeee'].makeHist(var, '', binning, perUnitWidth=False,
                                                             mergeOverflow=True)


    # PDF
    qqZZHists['pdf'] = {}
    for param, wts in qqWeights.iteritems():
        qqZZHists['pdf'][param] = {}
        for fgfz, wt in wts.iteritems():
            variations = qqZZ.makeHist2(var, 'Iteration$', '', binning,
                                        [100, 0, 100], 'pdfWeights * '+wt, False,
                                        mergeOverflowX=True)
            binRMSes = [Graph(variations.ProjectionY('slice{}'.format(i),
                                                     i+1,i+1)).GetRMS(2)
                        for i in xrange(variations.GetNbinsX())]

            qqZZHists['pdf'][param][fgfz] = qqZZHists[''][param][fgfz].clone()
            for b, rms in zip(qqZZHists['pdf'][param][fgfz].bins(), binRMSes):
                b.value += rms


    # alpha_s
    qqZZHists['alphas'] = {}
    for param, wts in qqWeights.iteritems():
        qqZZHists['alphas'][param] = {}
        for fgfz, wt in wts.iteritems():
            qqZZAlphaS1 = qqZZ.makeHist(var, '', binning,
                                        'scaleWeights[100] * '+wt,
                                        perUnitWidth=False,
                                        mergeOverflow=True)
            qqZZAlphaS2 = qqZZ.makeHist(var, '', binning,
                                        'scaleWeights[101] * '+wt,
                                        perUnitWidth=False,
                                        mergeOverflow=True)

            # PDF4LHC recommendation
            qqZZAlphaSErr = (qqZZAlphaS1 - qqZZAlphaS2) / 2.
            # Not sure which is up and which is down, so just do the
            # absolute value up and down.
            # Factor of 1.5 comes from slide 14 of
            # https://indico.cern.ch/event/459797/contributions/1961581/attachments/1181555/1800214/mcaod-Feb15-2016.pdf
            qqZZHists['alphas'][param][fgfz] = qqZZHists[''][param][fgfz].clone()
            for b, bShift in zip(qqZZHists['alphas'][param][fgfz], qqZZAlphaSErr):
                b.value += 1.5 * abs(bShift.value)


    # QCD scale
    scaleVarIndices = [1,2,3,4,6,8]

    qqZZHists['scale_up'] = {}
    qqZZHists['scale_dn'] = {}
    for param, wts in qqWeights.iteritems():
        qqZZHists['scale_up'][param] = {}
        qqZZHists['scale_dn'][param] = {}
        for fgfz, wt in wts.iteritems():
            qqZZVariations = [qqZZ.makeHist(var, '', binning,
                                            ('scaleWeights[{}] * '+wt).format(i),
                                            perUnitWidth=False,
                                            mergeOverflow=True)
                      for i in scaleVarIndices]

            qqZZHists['scale_up'][param][fgfz] = qqZZHists[''][param][fgfz].empty_clone()
            qqZZHists['scale_dn'][param][fgfz] = qqZZHists[''][param][fgfz].empty_clone()

            for bUp, bDn, thisBinAllHists in zip(qqZZHists['scale_up'][param][fgfz],
                                                 qqZZHists['scale_dn'][param][fgfz],
                                                 zip(*qqZZVariations)):
                bUp.value = max(b.value for b in thisBinAllHists)
                bDn.value = min(b.value for b in thisBinAllHists)

    bkgMCVariations = [bkgMC.makeHist(var, '', binning,
                                      'scaleWeights[{}]'.format(i),
                                      perUnitWidth=False,
                                      mergeOverflow=True)
                       for i in scaleVarIndices]

    bkgMCHists['scale_up'] = bkgMCHists[''].empty_clone()
    bkgMCHists['scale_dn'] = bkgMCHists[''].empty_clone()

    for bUp, bDn, thisBinAllHists in zip(bkgMCHists['scale_up'],
                                         bkgMCHists['scale_dn'],
                                         zip(*bkgMCVariations)):
        bUp.value = max(b.value for b in thisBinAllHists)
        bDn.value = min(b.value for b in thisBinAllHists)


    ### Since we have no LHE info for MCFM samples, we just vary the
    ### normalization uncertainty up and down by the cross section's
    ### PDF and scale uncertainties
    # gg LO uncertainty: +18%/-15% (MCFM, AN-2016-029)
    mcfmUnc = {'up':.18,'dn':-.15}
    for sys, shift in mcfmUnc.iteritems():
        ggZZHists['mcfmxsec_'+sys] = ggZZHists[''].clone() * (1.+shift)


    # save this all in a file
    for param in qqWeights:
        for (fg, fz) in qqWeights[param]:
            with root_open(join(outDir,
                                'mZZ_signal_aTGC-{}I{}_{}.root'.format(paramStrs[fg],
                                                                       paramStrs[fz],
                                                                       param)),
                           'recreate') as f:
                for syst in set(qqZZHists.keys()+ggZZHists.keys()):
                    h = hTrue.empty_clone(name='mZZ_reco_signal{}'.format('_'+syst if syst else ''))
                    try:
                        h += qqZZHists[syst][param][(fg,fz)]
                    except KeyError:
                        if syst in qqZZHists:
                            raise
                        h += qqZZHists[''][param][(fg,fz)]
                    try:
                        h += ggZZHists[syst]
                    except KeyError:
                        h += ggZZHists['']

                    h.write()


        # Use the points we got to parameterize N(aTGC)/N(SM) in each bin
        # as a function of fg and fz
        ratioGraphs = [Graph2D(len(qqZZHists[''][param])-2, # don't include overflow
                               type='default') for b in hTrue.bins()]
        for iWP, ((fg,fz), h) in enumerate(qqZZHists[''][param].iteritems()):
            yieldRatio = (h+ggZZHists['']) / sigSM
            for g, b in zip(ratioGraphs, yieldRatio.bins()):
                g.SetPoint(iWP, fg, fz, b.value)

        # fit ratios to a function of fg and fz and save the fits
        with root_open(join(outDir, 'mZZ_aTGCFits_{}.root'.format(param)),
                       'recreate') as f:
            ratioFits = []
            for i, bin in enumerate(hTrue.bins()):
                binfo = bin.x
                f = TF2('aTGCToSMRatioFit_M{}to{}'.format(binfo.low, binfo.high),
                        "[0]+[1]*x+[2]*y+[3]*x*x+[4]*y*y+[5]*x*y",
                        fgs[0], fgs[-1], fzs[0], fzs[-1])
                ratioGraphs[i].Fit(f, "WRN0EX0")
                ratioFits.append(f)
                f.write()


    # Data and background are totally independent of aTGC stuff, so they can
    # go in a single file
    data = standardZZData('zz', inData, 'smp')

    with root_open(join(outDir, 'mZZ_data-bkg.root'), 'recreate') as f:
        hData = data.makeHist(var, '', binning, perUnitWidth=False,
                              mergeOverflow=True)
        hData.name = 'mZZ_data_15p9ifb'
        hData.write()

        for syst in set(bkgMCHists.keys()+bkgDataHists.keys()):
            h = bkgMCHists[''].empty_clone(name='mZZ_bkg{}'.format('_'+syst if syst else ''))
            try:
                h += bkgMCHists[syst]
            except KeyError:
                h += bkgMCHists['']
            try:
                h += bkgDataHists[syst]
            except KeyError:
                h += bkgDataHists['']

            h.write()



if __name__ == '__main__':
    inData = 'uwvvNtuples_data_25nov2016'
    inMC = 'uwvvNtuples_mc_25nov2016'
    inATGC = '/data/nawoods/aTGCSherpaHistos'
    outDir = '/afs/cern.ch/user/n/nawoods/public/aTGC'
    if not exists(outDir):
        makedirs(outDir)
    fakeRateFile = 'fakeRate_08sep2016'
    puWeightFile = 'puWeight_69200_08sep2016.root'
    lumi = 15937.

    main(inData, inMC, inATGC, outDir, fakeRateFile, puWeightFile, lumi)
