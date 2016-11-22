import logging
from rootpy import log as rlog; rlog = rlog["/makeCards"]
# don't show most silly ROOT messages
logging.basicConfig(level=logging.WARNING)
rlog["/ROOT.TUnixSystem.SetDisplay"].setLevel(rlog.ERROR)

from rootpy.ROOT import Double

from SampleTools import MCSample, DataSample, SampleGroup, SampleStack
from Utilities import WeightStringMaker, Z_MASS
from Analysis.setupStandardSamples import *
from Analysis.weightHelpers import puWeight, baseMCWeight
from Metadata.metadata import sampleInfo

from os import environ as _env
from os.path import join as _join
from math import sqrt


_channels = ['eeee','eemm', 'mmmm']


def _getYield(sample, errorToo=False):
    hYield = sample.makeHist('1', '', [1,0.,2.], perUnitWidth=False)

    if errorToo:
        e = Double(0)
        y = hYield.IntegralAndError(0,hYield.GetNbinsX(),e)
        return y, e

    return hYield.Integral()


def main(inData, inMC, ana, cardName, fakeRateFile, puWeightFile, lumi):

    reco = {c:zzStackMCOnly(c, inMC, ana, puWeightFile, lumi) for c in _channels}
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

    recoSyst = {}
    for syst in ['eScaleUp', 'eScaleDn', 'eRhoResUp',
                 'eRhoResDn', 'ePhiResUp']:
        recoSyst[syst.lower()] = {
            c:zzStackMCOnly(c, inMC.replace('mc_','mc_{}_'.format(syst)),
                            ana, puWeightFile, lumi)
            for c in ['eeee','eemm']
            }
    for syst in ['mClosureUp','mClosureDn']:
        recoSyst[syst.lower()] = {
            c:zzStackMCOnly(c, inMC.replace('mc_','mc_{}_'.format(syst)),
                            ana, puWeightFile, lumi)
            for c in ['eemm','mmmm']
            }

    sigYield = {c:0. for c in _channels}
    bkgYield = {c:0. for c in _channels}
    sigSystSqr = {c:0. for c in _channels}
    bkgSystSqr = {c:0. for c in _channels}
    sigStatSqr = {c:0. for c in _channels}
    bkgStatSqr = {c:0. for c in _channels}
    nObs = {}

    for chan in _channels:
        cardNameChan = cardName + '_' + chan + '.txt'
        nominalWeight = baseMCWeight(chan, puWeightFile)

        with open(cardNameChan, 'w') as f:
            nSigSamples = len(reco[chan])
            colWidths = [max(len(s.name)+3, 8) for s in reco[chan]]

            lines = []
            lines.append('imax 1  number of channels')
            lines.append('jmax {}  number of backgrounds plus signals minus 1'.format(nSigSamples))
            lines.append('kmax *  number of nuisance parameters')
            lines.append('------------')

            lines.append('bin            1')

            hObs = data[chan].makeHist('1', '', [1,0.,2.], perUnitWidth=False)
            nObs[chan] = int(hObs.GetEntries())
            lines.append('observation    {}'.format(nObs[chan]))
            lines.append('------------')
            lines.append('bin                     '+''.join(['1'+' '*(wid-1) for i, wid in zip(range(nSigSamples+1),colWidths+[1])]))
            lines.append('process                 '+''.join([s.name+' '*(wid-len(s.name)) for s, wid in zip(reco[chan], colWidths)] + ['fake']))
            lines.append('process                 '+''.join([str(i+1)+' '*(wid-len(str(i+1))) for i,wid in zip(range(-1*nSigSamples, 1), colWidths+[1])]))

            yields, mcStatErrs = zip(*(_getYield(s,True) for s in reco[chan]))

            sigYield[chan] = sum(y for y in yields)

            bkgYield[chan], bkgStatErr = _getYield(bkg[chan], True)
            bkgStatSqr[chan] = bkgStatErr ** 2

            yields = [y for y in yields] + [bkgYield[chan]]

            lines.append('rate                    '+''.join([('{:.3f}'.format(y)+' '*(wid-len('{:.3f}'.format(y)))) if y else '1.e-8'+' '*(wid-5) for y,wid in zip(yields, colWidths+[7])]))
            lines.append('------------')

            errs = {}

            errs['mcStat_'+chan] = [e/y if y else 0 for e,y in zip(mcStatErrs,yields)] + [0]
            sigStatSqr[chan] = sum(e**2 for e in mcStatErrs)

            errs['bkgStat'] = [0.]*len(yields)
            errs['bkgStat'][-1] = bkgStatErr / bkgYield[chan]

            # PU uncertainty
            yield_puShift = {}
            for sys in ['up','dn']:
                wtStr = baseMCWeight(chan, puWeightFile, puSyst=sys)
                reco[chan].applyWeight(wtStr, True)
                yield_puShift[sys] = [_getYield(s) for s in reco[chan]]
                reco[chan].applyWeight(nominalWeight, True)

            errs['pu'] = [(abs(y-yUp)+abs(y-yDn))/(2.*y) for y, yUp, yDn in zip(yields[:-1],
                                                                                yield_puShift['up'],
                                                                                yield_puShift['dn'])] + [0]

            sigSystSqr[chan] += sum(((abs(y-yUp)+abs(y-yDn))/2.) ** 2 for y, yUp, yDn in zip(yields[:-1],
                                                                                         yield_puShift['up'],
                                                                                         yield_puShift['dn']))

            # Lepton efficiency uncertainty
            for lep in ['e','m']:
                if lep not in chan:
                    continue
                yield_lepEff = {}
                for sys in ['up','dn']:
                    wtOpts = {lep+'Syst':sys}
                    wtStr = baseMCWeight(chan, puWeightFile, **wtOpts)
                    reco[chan].applyWeight(wtStr, True)
                    yield_lepEff[sys] = [_getYield(s) for s in reco[chan]]
                    reco[chan].applyWeight(nominalWeight, True)

                errs[lep+'Eff'] = [(abs(y-yUp)+abs(y-yDn))/(2.*y) for y, yUp, yDn in zip(yields[:-1],
                                                                                         yield_lepEff['up'],
                                                                                         yield_lepEff['dn'])] + [0]

                sigSystSqr[chan] += sum(((abs(y-yUp)+abs(y-yDn))/2.) ** 2 for y, yUp, yDn in zip(yields[:-1],
                                                                                                 yield_lepEff['up'],
                                                                                                 yield_lepEff['dn']))

            # Lepton fake rate uncertainty
            for lep in ['e','m']:
                if lep not in chan:
                    continue

                yield_fr = {}
                for sys in ['up','dn']:
                    yield_fr[sys] = _getYield(bkgSyst[lep+sys][chan])

                errs[lep+'FR'] = [0]*len(yields)
                errs[lep+'FR'][-1] = (abs(yields[-1] - yield_fr['up']) + abs(yields[-1] - yield_fr['dn'])) / (2. * yields[-1])

                bkgSystSqr[chan] += (abs(yields[-1] - yield_fr['up']) + abs(yields[-1] - yield_fr['dn']) / (2.)) ** 2


            # Lepton scale/resolution
            if 'e' in chan:
                yield_ees = {}
                yield_eerRho = {}
                for sys in ['up','dn']:
                    yield_ees[sys] = [_getYield(s) for s in recoSyst['escale'+sys][chan]]
                    yield_eerRho[sys] = [_getYield(s) for s in recoSyst['erhores'+sys][chan]]

                errs['eScale'] = [(abs(y-yUp)+abs(y-yDn))/(2.*y) for y, yUp, yDn in zip(yields[:-1],
                                                                                        yield_ees['up'],
                                                                                        yield_ees['dn'])] + [0]
                errs['eRhoRes'] = [(abs(y-yUp)+abs(y-yDn))/(2.*y) for y, yUp, yDn in zip(yields[:-1],
                                                                                         yield_eerRho['up'],
                                                                                         yield_eerRho['dn'])] + [0]

                sigSystSqr[chan] += sum(((abs(y-yUp)+abs(y-yDn))/2.) ** 2 for y, yUp, yDn in zip(yields[:-1],
                                                                                                 yield_ees['up'],
                                                                                                 yield_ees['dn']))
                sigSystSqr[chan] += sum(((abs(y-yUp)+abs(y-yDn))/2.) ** 2 for y, yUp, yDn in zip(yields[:-1],
                                                                                                 yield_eerRho['up'],
                                                                                                 yield_eerRho['dn']))

                yield_eerPhi = [_getYield(s) for s in recoSyst['ephiresup'][chan]]
                errs['ePhiRes'] = [abs(y-yUp)/y for y, yUp in zip(yields[:-1],yield_eerPhi)] + [0]

                sigSystSqr[chan] += sum((abs(y-yUp)/y) ** 2. for y, yUp in zip(yields[:-1],yield_eerPhi))

            if 'm' in chan:
                yield_mEnergy = {}
                for sys in ['up','dn']:
                    yield_mEnergy[sys] = [_getYield(s) for s in recoSyst['mclosure'+sys][chan]]

                errs['mEnergy'] = [(abs(y-yUp)+abs(y-yDn))/(2.*y) for y, yUp, yDn in zip(yields[:-1],
                                                                                         yield_mEnergy['up'],
                                                                                         yield_mEnergy['dn'])] + [0]

                sigSystSqr[chan] += sum(((abs(y-yUp)+abs(y-yDn))/2.) ** 2 for y, yUp, yDn in zip(yields[:-1],
                                                                                                 yield_mEnergy['up'],
                                                                                                 yield_mEnergy['dn']))

            # few more by hand
            errs['pdf'] = [.01]*(len(yields)-1) + [0]
            errs['scale'] = [.01]*(len(yields)-1) + [0]
            errs['met'] = [0]*(len(yields)-1) + [0.01]
            errs['trigger'] = [.02]*(len(yields)-1) + [0]

            sigSystSqr[chan] += sum((e*y) ** 2 for e,y in zip(errs['pdf'],yields))
            sigSystSqr[chan] += sum((e*y) ** 2 for e,y in zip(errs['scale'],yields))
            sigSystSqr[chan] += sum((e*y) ** 2 for e,y in zip(errs['met'],yields))
            sigSystSqr[chan] += sum((e*y) ** 2 for e,y in zip(errs['trigger'],yields))


            # Make the rest of the card
            col1Width = 18
            lnN = 'lnN   '

            for errName, errors in errs.iteritems():
                colStrs = ['{:.3f}'.format(1.+e) if e else '-' for e in errors]
                otherCols = [s+' '*(wid-len(s)) for s, wid in zip(colStrs, colWidths+[8])]
                lines.append(''.join([errName,' '*(col1Width-len(errName)),lnN]+otherCols))


            f.write('\n'.join(lines))


    # print these numbers as a LaTeX table
    table = '''
\\begin{{table}}[htbp]
  \\begin{{center}}
    \\begin{{small}}
      \\begin{{tabular}}{{|c|c|c|c|c|}}
        \\hline
         & & & & \\\\ [-1.6ex]
        Decay & $N_{{\\mathrm{{ZZ}}}}^{{\\mathrm{{exp}}}}$ & Background & Total & Observed \\\\
        channel & & & expected & \\\\
        \\hline
        \\hline
        $4\Pgm$       & $ {mmmmSig:.2f} \pm {mmmmSigStat:.2f} \pm {mmmmSigSyst:.2f} $ & $ {mmmmBkg:.2f} \pm {mmmmBkgStat:.2f} \pm {mmmmBkgSyst:.2f} $ & $ {mmmmTot:.2f} \pm {mmmmTotStat:.2f} \pm {mmmmTotSyst:.2f} $ & $ {mmmmObs} $ \\\\
        $2\Pe 2\Pgm$  & $ {eemmSig:.2f} \pm {eemmSigStat:.2f} \pm {eemmSigSyst:.2f} $ & $ {eemmBkg:.2f} \pm {eemmBkgStat:.2f} \pm {eemmBkgSyst:.2f} $ & $ {eemmTot:.2f} \pm {eemmTotStat:.2f} \pm {eemmTotSyst:.2f} $ & $ {eemmObs} $ \\\\
        $4\Pe$        & $ {eeeeSig:.2f} \pm {eeeeSigStat:.2f} \pm {eeeeSigSyst:.2f} $ & $ {eeeeBkg:.2f} \pm {eeeeBkgStat:.2f} \pm {eeeeBkgSyst:.2f} $ & $ {eeeeTot:.2f} \pm {eeeeTotStat:.2f} \pm {eeeeTotSyst:.2f} $ & $ {eeeeObs} $ \\\\
        \\hline
        Total         & $ {totSig:.2f}  \pm {totSigStat:.2f}  \pm {totSigSyst:.2f} $  & $ {totBkg:.2f}  \pm {totBkgStat:.2f}  \pm {totBkgSyst:.2f} $  & $ {totTot:.2f}  \pm {totTotStat:.2f}  \pm {totTotSyst:.2f} $  & $ {totObs} $  \\\\
        \\hline
      \\end{{tabular}}
    \\end{{small}}
  \\vspace{{0.5cm}}
  \\caption{{ The observed and expected yield of events and estimated background
             yield evaluated from data, shown for each decay channel and summed
             to the total expected yield (``Total expected''). The first
             uncertainty is statistical, the second is systematic.
            }}
  \\end{{center}}
  \\label{{table:results}}
\\end{{table}}
'''

    toFormat = {
        'totSig' : 0.,
        'totSigStat' : 0.,
        'totSigSyst' : 0.,
        'totBkg' : 0.,
        'totBkgStat' : 0.,
        'totBkgSyst' : 0.,
        'totTot' : 0.,
        'totTotStat' : 0.,
        'totTotSyst' : 0.,
        'totObs' : 0,
        }

    for chan in _channels:
        tot = 0.
        totStatSqr = 0.
        totSystSqr = 0.
        toFormat[chan+'Sig'] = sigYield[chan]
        tot += sigYield[chan]
        toFormat[chan+'SigStat'] = sqrt(sigStatSqr[chan])
        totStatSqr += sigStatSqr[chan]
        toFormat[chan+'SigSyst'] = sqrt(sigSystSqr[chan])
        totSystSqr += sigStatSqr[chan]

        toFormat[chan+'Bkg'] = bkgYield[chan]
        tot += bkgStatSqr[chan]
        toFormat[chan+'BkgStat'] = sqrt(bkgStatSqr[chan])
        totStatSqr += bkgStatSqr[chan]
        toFormat[chan+'BkgSyst'] = sqrt(bkgSystSqr[chan])
        totSystSqr += bkgSystSqr[chan]

        toFormat[chan+'Tot'] = tot
        toFormat[chan+'TotStat'] = sqrt(totStatSqr)
        toFormat[chan+'TotSyst'] = sqrt(totSystSqr)

        toFormat[chan+'Obs'] = nObs[chan]

        toFormat['totSig'] += sigYield[chan]
        toFormat['totSigStat'] += sigStatSqr[chan]
        toFormat['totSigSyst'] += sigSystSqr[chan]
        toFormat['totBkg'] += bkgYield[chan]
        toFormat['totBkgStat'] += bkgStatSqr[chan]
        toFormat['totBkgSyst'] += bkgSystSqr[chan]
        toFormat['totTot'] += tot
        toFormat['totTotStat'] += totStatSqr
        toFormat['totTotSyst'] += totSystSqr

        toFormat['totObs'] += nObs[chan]

    for sOrB in ['Sig', 'Bkg', 'Tot']:
        for unc in ['Stat', 'Syst']:
            toFormat['tot'+sOrB+unc] = sqrt(toFormat['tot'+sOrB+unc])

    print table.format(**toFormat)


if __name__ == "__main__":

    from argparse import ArgumentParser

    defaultPath = _join(_env['zzt'], 'Analysis', 'savedResults')

    parser = ArgumentParser(description="Make combine cards to find an inclusive cross section")
    parser.add_argument('--dataDir', type=str, nargs='?',
                        default='uwvvNtuples_data_12oct2016',
                        help='Directory where data ntuples live')
    parser.add_argument('--mcDir', type=str, nargs='?',
                        default='uwvvNtuples_mc_12oct2016',
                        help='Directory where MC ntuples live')
    parser.add_argument('--analysis', type=str, nargs='?', default='smp',
                        help=('Which set of cuts to use. Options are "smp" '
                              '(default) and "z4l"'))
    parser.add_argument('--cardName', type=str, nargs='?',
                        default='combineCard',
                        help=('Start (i.e. everything but channel and ".txt") '
                              'of names for output files, including path if '
                              'different from $zzt/Analysis/savedResults,'))
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

    args=parser.parse_args()

    cardName = _join(defaultPath, args.cardName)

    main(args.dataDir, args.mcDir, args.analysis, cardName, args.fakeRateFile,
         args.puWeightFile, args.lumi)