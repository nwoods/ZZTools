import logging
from rootpy import log as rlog; rlog = rlog["/aTGCFullShape"]
# don't show most silly ROOT messages
logging.basicConfig(level=logging.WARNING)
rlog["/ROOT.TUnixSystem.SetDisplay"].setLevel(rlog.ERROR)

from rootpy.io import root_open
from rootpy.plotting import Graph, Graph2D, Canvas
from rootpy.ROOT import TF2

from Analysis.setupStandardSamples import standardZZMC, zzIrreducibleBkg, \
    standardZZBkg, standardZZData
from SampleTools import SampleGroup
from Analysis.weightHelpers import baseMCWeight

from os.path import join, exists
from os import makedirs


var = 'Mass'
binning = [100.*(i+1) for i in range(13)]+[2000.,2600]# range([100.+100.*i for i in range(6)] + [800.,1000.,1100.,1200.,1300.]
aTGCParams = ['f4','f5']

channels = ['eeee','eemm','mmmm']

sampleTemp = 'ZZTo4L-aTGC-{param}-fg{fg}-fz{fz}'

fgs = [-.0038, 0.0019, 0., 0.0019, 0.0038]
fzs = [-0.003, -0.0015, 0., 0.0015, 0.003]

paramStrs = {
    -0.0038 : 'm0p0038',
    -0.003  : 'm0p003',
    -0.0019 : 'm0p0019',
    -0.0015 : 'm0p0015',
    0.      : '0',
    0.0015  : '0p0015',
    0.0019  : '0p0019',
    0.003   : '0p003',
    0.0038  : '0p0038',
    }

fileNames = {
    p:{(fg,fz):sampleTemp.format(param=p,fg=paramStrs[fg],fz=paramStrs[fz])
       for fg in fgs
       for fz in fzs}
    for p in aTGCParams
    }

for p in aTGCParams:
    fileNames[p][(0.,0.)] = 'ZZTo4L-sherpa'

systDirNames = {
    'ees_up' : 'eScaleUp',
    'ees_dn' : 'eScaleDn',
    'eerPhi_up' : 'ePhiResUp',
    'eerRho_up' : 'eRhoResUp',
    'eerRho_dn' : 'eRhoResDn',
    'mClosure_up' : 'mClosureUp',
    'mClosure_dn' : 'mClosureDn',
    }

def main(inData, inMC, inATGC, outDir, fakeRateFile, puWeightFile, lumi):

    nominalWeight = baseMCWeight('zz', puWeightFile)

    # prepare non-aTGC samples

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

    ggZZSyst = {}
    bkgMCSyst = {}
    for s in ['ees_up', 'ees_dn', 'eerRho_up',
              'eerRho_dn', 'eerPhi_up']:
        syst = systDirNames[s]
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
                            fakeRateFile, lumi, sipCut=10.)
    bkgDataSyst = {
        'eup' : standardZZBkg('eeee,eemm', inData, inMC, 'smp', puWeightFile,
                              fakeRateFile, lumi, sipCut=10.,
                              eFakeRateSyst='up'),
        'edn' : standardZZBkg('eeee,eemm', inData, inMC, 'smp', puWeightFile,
                              fakeRateFile, lumi, sipCut=10.,
                              eFakeRateSyst='dn'),
        'mup' : standardZZBkg('eemm,mmmm', inData, inMC, 'smp', puWeightFile,
                              fakeRateFile, lumi, sipCut=10.,
                              mFakeRateSyst='up'),
        'mdn' : standardZZBkg('eemm,mmmm', inData, inMC, 'smp', puWeightFile,
                              fakeRateFile, lumi, sipCut=10.,
                              mFakeRateSyst='dn'),
        }

    # weight MC
    ggZZ.applyWeight(nominalWeight, True)
    for s in ggZZSyst.values():
        s.applyWeight(nominalWeight, True)
    bkgMC.applyWeight(nominalWeight, True)
    for s in bkgMCSyst.values():
        s.applyWeight(nominalWeight, True)

    # and make other weights to use later
    lepEffReweight = {}
    for lep in ['e','m']:
        lepEffReweight[lep] = {}
        for sys in ['up','dn']:
            wtArg = {lep+'Syst':sys}
            lepEffReweight[lep][sys] = baseMCWeight('zz', puWeightFile, **wtArg)
    puReweight = {}
    for sys in ['up','dn']:
        puReweight[sys] = baseMCWeight('zz', puWeightFile, puSyst=sys)


    ### Get all non-aTGC histograms now

    print "Making non-aTGC histograms"

    # Nominal
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

    # Systematics

    # PU weight systematic
    for sys in ['up','dn']:
        ggZZ.applyWeight(puReweight[sys], True)
        bkgMC.applyWeight(puReweight[sys], True)

        ggZZHists['pu_'+sys] = ggZZ.makeHist(var, '', binning, perUnitWidth=False,
                                             mergeOverflow=True)

        bkgMCHists['pu_'+sys] = bkgMC.makeHist(var, '', binning, perUnitWidth=False,
                                               mergeOverflow=True)

        ggZZ.applyWeight(nominalWeight, True)
        bkgMC.applyWeight(nominalWeight, True)


    # Lepton efficiency systematics
    for lep in ['e','m']:
        for sys in ['up','dn']:
            wtArg = {lep+'Syst':sys}
            ggZZ.applyWeight(lepEffReweight[lep][sys], True)
            bkgMC.applyWeight(lepEffReweight[lep][sys], True)

            ggZZHists[lep+'Eff_'+sys] = ggZZ.makeHist(var, '', binning, perUnitWidth=False,
                                                      mergeOverflow=True)

            bkgMCHists[lep+'Eff_'+sys] = bkgMC.makeHist(var, '', binning, perUnitWidth=False,
                                                        mergeOverflow=True)

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
        ggZZHists['eEnergyResolutionRho_'+sys] = ggZZSyst['eerRho_'+sys].makeHist(var, '', binning, perUnitWidth=False,
                                                                                  mergeOverflow=True)
        ggZZHists['eEnergyResolutionRho_'+sys] += ggZZ['mmmm'].makeHist(var, '', binning, perUnitWidth=False,
                                                                        mergeOverflow=True)

        bkgMCHists['eEnergyResolutionRho_'+sys] = bkgMCSyst['eerRho_'+sys].makeHist(var, '', binning, perUnitWidth=False,
                                                                                    mergeOverflow=True)
        bkgMCHists['eEnergyResolutionRho_'+sys] += bkgMC['mmmm'].makeHist(var, '', binning, perUnitWidth=False,
                                                                          mergeOverflow=True)


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
        ggZZHists['mEnergy_'+sys] = ggZZSyst['mClosure_'+sys].makeHist(var, '', binning, perUnitWidth=False,
                                                                       mergeOverflow=True)
        ggZZHists['mEnergy_'+sys] += ggZZ['eeee'].makeHist(var, '', binning, perUnitWidth=False,
                                                           mergeOverflow=True)

        bkgMCHists['mEnergy_'+sys] = bkgMCSyst['mClosure_'+sys].makeHist(var, '', binning, perUnitWidth=False,
                                                                         mergeOverflow=True)
        bkgMCHists['mEnergy_'+sys] += bkgMC['eeee'].makeHist(var, '', binning, perUnitWidth=False,
                                                             mergeOverflow=True)



    ### Since we have no LHE info for MCFM samples, we just vary the
    ### normalization uncertainty up and down by the cross section's
    ### PDF and scale uncertainties
    # gg LO uncertainty: +18%/-15% (MCFM, AN-2016-029)
    mcfmUnc = {'up':.18,'dn':-.15}
    for sys, shift in mcfmUnc.iteritems():
        ggZZHists['mcfmxsec_'+sys] = ggZZHists[''] * (1.+shift)


    ### Make the same histograms for each aTGC sample
    aTGCHists = {}
    for param in aTGCParams:
        aTGCHists[param] = {}
        for (fg, fz), fName in fileNames[param].iteritems():

            print "making histograms for {}, fg={}, fz={}".format(param,fg,fz)

            try:
                aTGCByChan = {}
                for c in channels:
                    aTGCByChan[c] = standardZZMC(c, inATGC, fName, 'smp',
                                                 puWeightFile, lumi)
            except ValueError:
                rlog.warning("Can't find files for {}: fg={}, fz={} -- skipping.".format(param,fg,fz))
                continue

            aTGC = SampleGroup(fName, 'zz', aTGCByChan, True)

            aTGCSyst = {}
            for s in ['ees_up', 'ees_dn', 'eerRho_up',
                      'eerRho_dn', 'eerPhi_up']:
                syst = systDirNames[s]

                aTGCSystByChan = {}
                for c in ['eeee','eemm']:
                    aTGCSystByChan[c] = standardZZMC(c, inATGC.replace('mc_','mc_{}_'.format(syst)),
                                                     fName, 'smp',
                                                     puWeightFile, lumi)

                aTGCSyst[s] = SampleGroup(fName, 'eeee,eemm', aTGCSystByChan, True)

            for s in ['mClosure_up','mClosure_dn']:
                syst = systDirNames[s]
                aTGCSystByChan = {}
                for c in ['eemm','mmmm']:
                    aTGCSystByChan[c] = standardZZMC(c, inATGC.replace('mc_','mc_{}_'.format(syst)),
                                                     fName, 'smp',
                                                     puWeightFile, lumi)

                aTGCSyst[s] = SampleGroup('GluGluZZ', 'eemm,mmmm', aTGCSystByChan, True)


            # Nominal
            aTGCHists[param][(fg,fz)] = {
                '' : aTGC.makeHist(var, '', binning, perUnitWidth=False,
                                   mergeOverflow=True)
                }

            # Systematics

            # PU weight systematic
            for sys in ['up','dn']:
                aTGC.applyWeight(puReweight[sys], True)

                aTGCHists[param][(fg,fz)]['pu_'+sys] = aTGC.makeHist(var, '',
                                                                     binning,
                                                                     perUnitWidth=False,
                                                                     mergeOverflow=True)

                aTGC.applyWeight(nominalWeight, True)


            # Lepton efficiency systematics
            for lep in ['e','m']:
                for sys in ['up','dn']:
                    wtArg = {lep+'Syst':sys}
                    aTGC.applyWeight(lepEffReweight[lep][sys], True)

                    aTGCHists[param][(fg,fz)][lep+'Eff_'+sys] = aTGC.makeHist(var, '', binning, perUnitWidth=False,
                                                                              mergeOverflow=True)

                    aTGC.applyWeight(nominalWeight, True)


            # electron energy scale
            for sys in ['up','dn']:
                aTGCHists[param][(fg,fz)]['eEnergyScale_'+sys] = aTGCSyst['ees_'+sys].makeHist(var, '', binning, perUnitWidth=False,
                                                                                               mergeOverflow=True)
                aTGCHists[param][(fg,fz)]['eEnergyScale_'+sys] += aTGC['mmmm'].makeHist(var, '', binning, perUnitWidth=False,
                                                                                        mergeOverflow=True)


            # electron energy resolution
            for sys in ['up','dn']:
                aTGCHists[param][(fg,fz)]['eEnergyResolutionRho_'+sys] = aTGCSyst['eerRho_'+sys].makeHist(var, '', binning, perUnitWidth=False,
                                                                                                          mergeOverflow=True)
                aTGCHists[param][(fg,fz)]['eEnergyResolutionRho_'+sys] += aTGC['mmmm'].makeHist(var, '', binning, perUnitWidth=False,
                                                                                                mergeOverflow=True)

            aTGCHists[param][(fg,fz)]['eEnergyResolutionPhi'] = aTGCSyst['eerPhi_up'].makeHist(var, '', binning, perUnitWidth=False,
                                                                                               mergeOverflow=True)
            aTGCHists[param][(fg,fz)]['eEnergyResolutionPhi'] += aTGC['mmmm'].makeHist(var, '', binning, perUnitWidth=False,
                                                                                       mergeOverflow=True)


            # muon energy scale/resolution
            for sys in ['up','dn']:
                aTGCHists[param][(fg,fz)]['mEnergy_'+sys] = aTGCSyst['mClosure_'+sys].makeHist(var, '', binning, perUnitWidth=False,
                                                                                               mergeOverflow=True)
                aTGCHists[param][(fg,fz)]['mEnergy_'+sys] += aTGC['eeee'].makeHist(var, '', binning, perUnitWidth=False,
                                                                                   mergeOverflow=True)


    # save this all in a file
    hSigSM = aTGCHists[aTGCParams[0]][(0.,0.)][''] + ggZZHists['']

    for param in aTGCHists:
        for (fg, fz) in aTGCHists[param]:
            with root_open(join(outDir,
                                'mZZ_signal_aTGC-{}I{}_{}.root'.format(paramStrs[fg],
                                                                       paramStrs[fz],
                                                                       param)),
                           'recreate') as f:
                for syst in set(aTGCHists[param][(fg,fz)].keys()+ggZZHists.keys()):
                    h = hSigSM.empty_clone(name='mZZ_reco_signal{}'.format('_'+syst if syst else ''))
                    try:
                        h += aTGCHists[param][(fg,fz)][syst]
                    except KeyError:
                        if syst in aTGCHists[param][(fg,fz)]:
                            raise
                        h += aTGCHists[param][(fg,fz)]['']
                    try:
                        h += ggZZHists[syst]
                    except KeyError:
                        h += ggZZHists['']

                    h.write()

        # Use the points we got to parameterize N(aTGC)/N(SM) in each bin
        # as a function of fg and fz
        ratioGraphs = [Graph2D(len(aTGCHists[param]),
                               type='default') for b in hSigSM.bins()]
        for iWP, ((fg,fz), hs) in enumerate(aTGCHists[param].iteritems()):
            h = hs['']
            yieldRatio = (h+ggZZHists['']) / hSigSM
            for g, b in zip(ratioGraphs, yieldRatio.bins()):
                g.SetPoint(iWP, fg, fz, b.value)



        # fit ratios to a function of fg and fz and save the fits
        with root_open(join(outDir, 'mZZ_aTGCFits_{}.root'.format(param)),
                       'recreate') as f:
            ratioFits = []
            for i, bin in enumerate(hSigSM.bins()):
                binfo = bin.x
                f = TF2('aTGCToSMRatioFit_M{}to{}'.format(int(binfo.low), int(binfo.high)),
                        "[0]+[1]*x+[2]*y+[3]*x*x+[4]*y*y+[5]*x*y",
                        min(fgs), max(fgs), min(fzs), max(fzs))
                ratioGraphs[i].Fit(f, "WRN0EX0")
                ratioFits.append(f)
                if binfo.low < 1100:
                    print f.Eval(0.,0.)
                f.write()
                if binfo.low < 1100:
                    print f.Eval(0.,0.)
                ratioGraphs[i].name = 'ratioGraph_M{}to{}'.format(int(binfo.low),int(binfo.high))
                ratioGraphs[i].write()


    # Data and background are totally independent of aTGC stuff, so they can
    # go in a single file
    data = standardZZData('zz', inData, 'smp')

    with root_open(join(outDir, 'mZZ_data-bkg.root'), 'recreate') as f:
        hData = data.makeHist(var, '', binning, perUnitWidth=False,
                              mergeOverflow=True)
        hData.name = 'mZZ_data_35p9ifb'
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

        print "{:.2f} + {:.2f} = {:.2f}".format(bkgDataHists[''].Integral(),bkgMCHists[''].Integral(),bkgDataHists[''].Integral()+bkgMCHists[''].Integral())

if __name__ == '__main__':
    inData = 'uwvvNtuples_data_10mar2017_LooseSIPLooseVtx'
    inMC = 'uwvvNtuples_mc_23mar2017_LooseSIPLooseVtx'
    inATGC = 'uwvvNtuples_mc_21feb2017_aTGC_LooseSIP'
    outDir = '/afs/cern.ch/user/n/nawoods/public/aTGC_reminiAOD_LooseSIP'
    if not exists(outDir):
        makedirs(outDir)
    fakeRateFile = 'fakeRate_20feb2017_LooseSIP'
    puWeightFile = 'puWeight_69200_24jan2017.root'
    lumi = 35860.

    main(inData, inMC, inATGC, outDir, fakeRateFile, puWeightFile, lumi)
