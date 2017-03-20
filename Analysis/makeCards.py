import logging
from rootpy import log as rlog; rlog = rlog["/makeCards"]
# don't show most silly ROOT messages
logging.basicConfig(level=logging.WARNING)
rlog["/ROOT.TUnixSystem.SetDisplay"].setLevel(rlog.ERROR)

from rootpy.ROOT import Double
from rootpy.plotting import Graph

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


def main(inData, inMC, ana, cardName, fakeRateFile, puWeightFile, lumi,
         sfFromHists=False):

    reco = {c:zzStackSignalOnly(c, inMC, ana, puWeightFile,
                                lumi, scaleFactorsFromHists=sfFromHists)
            for c in _channels}
    bkg = standardZZBkg('zz', inData, inMC, ana, puWeightFile,
                        fakeRateFile, lumi)
    irreducible = zzIrreducibleBkg('zz', inMC, ana, puWeightFile, lumi,
                                   scaleFactorsFromHists=sfFromHists)

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
    irreducibleSyst = {}
    for syst in ['eScaleUp', 'eScaleDn', 'eRhoResUp',
                 'eRhoResDn', 'ePhiResUp']:
        recoSyst[syst.lower()] = {
            c:zzStackSignalOnly(c, inMC.replace('mc_','mc_{}_'.format(syst)),
                                ana, puWeightFile, lumi,
                                scaleFactorsFromHists=sfFromHists)
            for c in ['eeee','eemm']
            }
        irreducibleSyst[syst.lower()] = zzIrreducibleBkg('eeee,eemm',
                                                         inMC.replace('mc_','mc_{}_'.format(syst)),
                                                         ana, puWeightFile,
                                                         lumi,
                                                         scaleFactorsFromHists=sfFromHists)
    for syst in ['mClosureUp','mClosureDn']:
        recoSyst[syst.lower()] = {
            c:zzStackSignalOnly(c, inMC.replace('mc_','mc_{}_'.format(syst)),
                                ana, puWeightFile, lumi,
                                scaleFactorsFromHists=sfFromHists)
            for c in ['eemm','mmmm']
            }
        irreducibleSyst[syst.lower()] = zzIrreducibleBkg('eemm,mmmm',
                                                         inMC.replace('mc_','mc_{}_'.format(syst)),
                                                         ana, puWeightFile,
                                                         lumi,
                                                         scaleFactorsFromHists=sfFromHists)

    sigYield = {c:0. for c in _channels}
    bkgYield = {c:0. for c in _channels}
    irrYield = {c:0. for c in _channels}
    sigSystByChan = {}
    bkgSystByChan = {}
    irrSystByChan = {}
    sigStatSqr = {c:0. for c in _channels}
    bkgStatSqr = {c:0. for c in _channels}
    irrStatSqr = {c:0. for c in _channels}
    nObs = {}

    for chan in _channels:
        cardNameChan = cardName + '_' + chan + '.txt'
        nominalWeight = baseMCWeight(chan, puWeightFile,
                                     scaleFactorsFromHists=sfFromHists)

        with open(cardNameChan, 'w') as f:
            nSigSamples = len(reco[chan])
            nBkgSamples = 2
            nMCSamples = nSigSamples + nBkgSamples - 1
            names = [s.name for s in reco[chan]] + ['irreducible','fake']
            colWidths = [max(len(n)+3, 8) for n in names]

            lines = []
            lines.append('imax 1  number of channels')
            lines.append('jmax {}  number of backgrounds plus signals minus 1'.format(len(names)-1))
            lines.append('kmax *  number of nuisance parameters')
            lines.append('------------')

            lines.append('bin            1')

            hObs = data[chan].makeHist('1', '', [1,0.,2.], perUnitWidth=False)
            nObs[chan] = int(hObs.GetEntries())
            lines.append('observation    {}'.format(nObs[chan]))
            lines.append('------------')
            lines.append('bin                     '+''.join(['1'+' '*(wid-1) for wid in colWidths]))
            lines.append('process                 '+''.join([n+' '*(wid-len(n)) for n, wid in zip(names, colWidths)]))
            lines.append('process                 '+''.join([str(i+1)+' '*(wid-len(str(i+1))) for i,wid in zip(range(-1*nSigSamples, nBkgSamples),
                                                                                                               colWidths)]))

            yields, sigStatErrs = zip(*(_getYield(s,True) for s in reco[chan]))

            yields = [y for y in yields]
            sigStatErrs = [e for e in sigStatErrs]

            sigYield[chan] = sum(y for y in yields)

            bkgYield[chan], bkgStatErr = _getYield(bkg[chan], True)
            bkgStatSqr[chan] = bkgStatErr ** 2

            irrYield[chan], irrStatErr = _getYield(irreducible[chan], True)
            irrStatSqr[chan] = irrStatErr ** 2

            yields += [irrYield[chan],bkgYield[chan]]

            lines.append('rate                    '+''.join([('{:.3f}'.format(y)+' '*(wid-len('{:.3f}'.format(y)))) if y else '1.e-8'+' '*(wid-5) for y,wid in zip(yields, colWidths)]))
            lines.append('------------')

            errs = {}

            for i, (e, y) in enumerate(zip(sigStatErrs+[irrStatErr],yields[:-1])):
                mcStatErrName = 'mcStat_{}_{}'.format(chan,i)
                errs[mcStatErrName] = [0.]*len(yields)
                errs[mcStatErrName][i] = e/y

            sigStatSqr[chan] = sum(e**2 for e in sigStatErrs)

            errs['bkgStat'] = [0.]*len(yields)
            errs['bkgStat'][-1] = bkgStatErr / bkgYield[chan]

            # PU uncertainty
            yield_puShift = {}
            for sys in ['up','dn']:
                wtStr = baseMCWeight(chan, puWeightFile, puSyst=sys,
                                     scaleFactorsFromHists=sfFromHists)
                reco[chan].applyWeight(wtStr, True)
                irreducible[chan].applyWeight(wtStr, True)
                yield_puShift[sys] = [_getYield(s) for s in reco[chan]]+[_getYield(irreducible[chan])]
                reco[chan].applyWeight(nominalWeight, True)
                irreducible[chan].applyWeight(nominalWeight, True)

            errs['pu'] = [(abs(y-yUp)+abs(y-yDn))/(2.*y) for y, yUp, yDn in zip(yields[:-1],
                                                                                yield_puShift['up'],
                                                                                yield_puShift['dn'])] + [0]

            if 'pu' not in sigSystByChan:
                sigSystByChan['pu'] = {}
            sigSystByChan['pu'][chan] = sum((abs(y-yUp)+abs(y-yDn))/2. for y, yUp, yDn in zip(yields[:nSigSamples],
                                                                                              yield_puShift['up'][:nSigSamples],
                                                                                              yield_puShift['dn'][:nSigSamples])
                                            )

            if 'pu' not in irrSystByChan:
                irrSystByChan['pu'] = {}
            irrSystByChan['pu'][chan] = ((abs(yields[-2]-yield_puShift['up'][-1])+abs(yields[-2]-yield_puShift['dn'][-1]))/2)


            # Lepton efficiency uncertainty
            for lep in ['e','m']:
                if lep not in chan:
                    continue
                yield_lepEff = {}
                for sys in ['up','dn']:
                    wtOpts = {lep+'Syst':sys}
                    wtStr = baseMCWeight(chan, puWeightFile,
                                         scaleFactorsFromHists=sfFromHists,
                                         **wtOpts)
                    reco[chan].applyWeight(wtStr, True)
                    irreducible[chan].applyWeight(wtStr, True)
                    yield_lepEff[sys] = [_getYield(s) for s in reco[chan]]+[_getYield(irreducible[chan])]
                    reco[chan].applyWeight(nominalWeight, True)
                    irreducible[chan].applyWeight(nominalWeight, True)

                errs[lep+'Eff'] = [(abs(y-yUp)+abs(y-yDn))/(2.*y) for y, yUp, yDn in zip(yields[:-1],
                                                                                         yield_lepEff['up'],
                                                                                         yield_lepEff['dn'])] + [0]

                if lep+'Eff' not in sigSystByChan:
                    sigSystByChan[lep+'Eff'] = {}
                sigSystByChan[lep+'Eff'][chan] = sum((abs(y-yUp)+abs(y-yDn))/2. for y, yUp, yDn in zip(yields[:nSigSamples],
                                                                                                       yield_lepEff['up'][:nSigSamples],
                                                                                                       yield_lepEff['dn'][:nSigSamples])
                                                     )

                if lep+'Eff' not in irrSystByChan:
                    irrSystByChan[lep+'Eff'] = {}
                irrSystByChan[lep+'Eff'][chan] = ((abs(yields[-2]-yield_lepEff['up'][-1])+abs(yields[-2]-yield_lepEff['dn'][-1]))/2)


            # Lepton fake rate uncertainty
            for lep in ['e','m']:
                if lep not in chan:
                    continue

                yield_fr = {}
                for sys in ['up','dn']:
                    yield_fr[sys] = _getYield(bkgSyst[lep+sys][chan])

                errs[lep+'FR'] = [0]*len(yields)
                errs[lep+'FR'][-1] = (abs(yields[-1] - yield_fr['up']) + abs(yields[-1] - yield_fr['dn'])) / (2. * yields[-1])

                if lep+'FR' not in bkgSystByChan:
                    bkgSystByChan[lep+'FR'] = {}
                bkgSystByChan[lep+'FR'][chan] = (abs(yields[-1] - yield_fr['up']) + abs(yields[-1] - yield_fr['dn']) / (2.))


            # Lepton scale/resolution
            if 'e' in chan:
                yield_ees = {}
                yield_eerRho = {}
                for sys in ['up','dn']:
                    yield_ees[sys] = [_getYield(s) for s in recoSyst['escale'+sys][chan]]
                    yield_ees[sys] += [_getYield(irreducibleSyst['escale'+sys][chan])]
                    yield_eerRho[sys] = [_getYield(s) for s in recoSyst['erhores'+sys][chan]]
                    yield_eerRho[sys] += [_getYield(irreducibleSyst['erhores'+sys][chan])]

                errs['eScale'] = [(abs(y-yUp)+abs(y-yDn))/(2.*y) for y, yUp, yDn in zip(yields[:-1],
                                                                                        yield_ees['up'],
                                                                                        yield_ees['dn'])] + [0]
                errs['eRhoRes'] = [(abs(y-yUp)+abs(y-yDn))/(2.*y) for y, yUp, yDn in zip(yields[:-1],
                                                                                         yield_eerRho['up'],
                                                                                         yield_eerRho['dn'])] + [0]

                if 'eScale' not in sigSystByChan:
                    sigSystByChan['eScale'] = {}
                sigSystByChan['eScale'][chan] = sum((abs(y-yUp)+abs(y-yDn))/2. for y, yUp, yDn in zip(yields[:nSigSamples],
                                                                                                      yield_ees['up'][:nSigSamples],
                                                                                                      yield_ees['dn'][:nSigSamples])
                                                    )
                if 'eRhoRes' not in sigSystByChan:
                    sigSystByChan['eRhoRes'] = {}
                sigSystByChan['eRhoRes'][chan] = sum((abs(y-yUp)+abs(y-yDn))/2. for y, yUp, yDn in zip(yields[:-1][:nSigSamples],
                                                                                                       yield_eerRho['up'][:nSigSamples],
                                                                                                       yield_eerRho['dn'][:nSigSamples])
                                                     )

                if 'eScale' not in irrSystByChan:
                    irrSystByChan['eScale'] = {}
                irrSystByChan['eScale'][chan] = ((abs(yields[-2]-yield_ees['up'][-1])+abs(yields[-2]-yield_ees['dn'][-1]))/2)
                if 'eRhoRes' not in irrSystByChan:
                    irrSystByChan['eRhoRes'] = {}
                irrSystByChan['eRhoRes'][chan] = ((abs(yields[-2]-yield_eerRho['up'][-1])+abs(yields[-2]-yield_eerRho['dn'][-1]))/2)

                yield_eerPhi = [_getYield(s) for s in recoSyst['ephiresup'][chan]]
                yield_eerPhi += [_getYield(irreducibleSyst['ephiresup'][chan])]
                errs['ePhiRes'] = [abs(y-yUp)/y for y, yUp in zip(yields[:-1],yield_eerPhi)] + [0]

                if 'ePhiRes' not in sigSystByChan:
                    sigSystByChan['ePhiRes'] = {}
                sigSystByChan['ePhiRes'][chan] = sum(abs(y-yUp)/y for y, yUp in zip(yields[:nSigSamples],yield_eerPhi[:nSigSamples]))

                if 'ePhiRes' not in irrSystByChan:
                    irrSystByChan['ePhiRes'] = {}
                irrSystByChan['ePhiRes'][chan] = ((abs(yields[-2]-yield_eerPhi[-1])+abs(yields[-2]-yield_eerPhi[-1]))/2)


            if 'm' in chan:
                yield_mEnergy = {}
                for sys in ['up','dn']:
                    yield_mEnergy[sys] = [_getYield(s) for s in recoSyst['mclosure'+sys][chan]]
                    yield_mEnergy[sys] += [_getYield(irreducibleSyst['mclosure'+sys][chan])]

                errs['mEnergy'] = [(abs(y-yUp)+abs(y-yDn))/(2.*y) for y, yUp, yDn in zip(yields[:-1],
                                                                                         yield_mEnergy['up'],
                                                                                         yield_mEnergy['dn'])] + [0]

                if 'mEnergy' not in sigSystByChan:
                    sigSystByChan['mEnergy'] = {}
                sigSystByChan['mEnergy'][chan] = sum((abs(y-yUp)+abs(y-yDn))/2. for y, yUp, yDn in zip(yields[:nSigSamples],
                                                                                                       yield_mEnergy['up'][:nSigSamples],
                                                                                                       yield_mEnergy['dn'][:nSigSamples])
                                                     )

                if 'mEnergy' not in irrSystByChan:
                    irrSystByChan['mEnergy'] = {}
                irrSystByChan['mEnergy'][chan] = ((abs(yields[-2]-yield_mEnergy['up'][-1])+abs(yields[-2]-yield_mEnergy['dn'][-1]))/2)


            ###############################
            # Scale errors matter and are calculated fully for the
            # irreducible backgrounds, because here the resulting normalization
            # difference matters. For signal, they will cancel out when the
            # signal strength is multiplied by the cross section, so we just
            # use a flat 1% uncertainty on the acceptance.
            ###############################

            irrScaleVariations = []
            variationIndices = [1,2,3,4,6,8]
            for varInd in variationIndices:
                irreducible[chan].applyWeight('scaleWeights[{}]'.format(varInd))
                irrScaleVariations.append(_getYield(irreducible[chan]))
                irreducible[chan].applyWeight(nominalWeight, True)

            irrScaleErrUp = max(irrScaleVariations) - yields[-2]
            irrScaleErrDn = yields[-2] - min(irrScaleVariations)

            if 'scale_norm' not in irrSystByChan:
                irrSystByChan['scale_norm'] = {}
            irrSystByChan['scale_norm'][chan] = max(abs(irrScaleErrUp),abs(irrScaleErrDn))

            errs['scale_norm'] = [0]*len(yields)
            errs['scale_norm'][-2] = irrSystByChan['scale_norm'][chan] / yields[-2]


            # few more by hand
            errs['pdf_acc'] = [0.01]*nSigSamples + [0,0]
            errs['scale_acc'] = [0.01]*nSigSamples + [0,0]
            errs['met'] = [0]*(nMCSamples) + [0.01]
            errs['trigger'] = [.04 if ana == 'z4l' and chan == 'eeee' else .02]*(nMCSamples) + [0]

            if 'pdf_acc' not in sigSystByChan:
                sigSystByChan['pdf_acc'] = {}
            sigSystByChan['pdf_acc'][chan] = sum(0.01*y for y in yields[:nSigSamples])

            if 'scale_acc' not in sigSystByChan:
                sigSystByChan['scale_acc'] = {}
            sigSystByChan['scale_acc'][chan] = sum(0.01*y for y in yields[:nSigSamples])

            if 'trigger' not in sigSystByChan:
                sigSystByChan['trigger'] = {}
            sigSystByChan['trigger'][chan] = sum((e*y) for e,y in zip(errs['trigger'],yields[:nSigSamples]))

            if 'met' not in bkgSystByChan:
                bkgSystByChan['met'] = {}
            bkgSystByChan['met'][chan] = (errs['met'][-1]*yields[-1])

            if 'trigger' not in irrSystByChan:
                irrSystByChan['trigger'] = {}
            irrSystByChan['trigger'][chan] = (errs['trigger'][-2]*yields[-2])


            # Make the rest of the card
            col1Width = 18
            lnN = 'lnN   '

            for errName, errors in errs.iteritems():
                colStrs = ['{:.3f}'.format(1.+e) if e else '-' for e in errors]
                otherCols = [s+' '*(wid-len(s)) for s, wid in zip(colStrs, colWidths)]
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

    toFormat = {}

    allSystNames = set(sigSystByChan.keys() + bkgSystByChan.keys() + irrSystByChan.keys())

    for chan in _channels:

        toFormat[chan+'Sig'] = sigYield[chan]
        toFormat[chan+'SigStat'] = sqrt(sigStatSqr[chan])
        toFormat[chan+'SigSyst'] = sqrt(sum(sysVals.get(chan, 0.) ** 2 for sysVals in sigSystByChan.values()))

        toFormat[chan+'Bkg'] = bkgYield[chan] + irrYield[chan]
        toFormat[chan+'BkgStat'] = sqrt(bkgStatSqr[chan] + irrStatSqr[chan])
        toFormat[chan+'BkgSyst'] = sqrt(sum(sysVals.get(chan, 0.) ** 2 for sysVals in bkgSystByChan.values()) + \
                                        sum(sysVals.get(chan, 0.) ** 2 for sysVals in irrSystByChan.values()))

        toFormat[chan+'Tot'] = toFormat[chan+'Sig'] + toFormat[chan+'Bkg']
        toFormat[chan+'TotStat'] = sqrt(sigStatSqr[chan] + bkgStatSqr[chan] + irrStatSqr[chan])
        chanSystTot = {sys : sigSystByChan.get(sys, {}).get(chan, 0.) + bkgSystByChan.get(sys, {}).get(chan, 0.) + irrSystByChan.get(sys, {}).get(chan, 0.) for sys in allSystNames}
        toFormat[chan+'TotSyst'] = sqrt(sum(v ** 2 for v in chanSystTot.values()))

        toFormat[chan+'Obs'] = nObs[chan]

    toFormat['totSig'] = sum(sigYield.values())
    toFormat['totSigStat'] = sqrt(sum(sigStatSqr.values()))
    sigSystTot = {sysName : sum(sys.get(c, 0.) for c in _channels) for sysName, sys in sigSystByChan.iteritems()}
    toFormat['totSigSyst'] = sqrt(sum(v ** 2 for v in sigSystTot.values()))

    toFormat['totBkg'] = sum(bkgYield.values()) + sum(irrYield.values())
    toFormat['totBkgStat'] = sqrt(sum(bkgStatSqr.values()) + sum(irrStatSqr.values()))
    bkgSystTot = {sysName : sum(sys.get(c, 0.) for c in _channels) for sysName, sys in bkgSystByChan.iteritems()}
    irrSystTot = {sysName : sum(sys.get(c, 0.) for c in _channels) for sysName, sys in irrSystByChan.iteritems()}
    toFormat['totBkgSyst'] = sqrt(sum(v ** 2 for v in bkgSystTot.values()) + sum(v ** 2 for v in irrSystTot.values()))

    toFormat['totTot'] = toFormat['totSig'] + toFormat['totBkg']
    toFormat['totTotStat'] = sqrt(sum(sigStatSqr.values()) + sum(bkgStatSqr.values()) + sum(irrStatSqr.values()))
    allSystNames = set(sigSystTot.keys() + bkgSystTot.keys() + irrSystTot.keys())
    totSystTot = {sys : sigSystTot.get(sys, 0.) + bkgSystTot.get(sys, 0.) + irrSystTot.get(sys, 0.) for sys in allSystNames}
    toFormat['totTotSyst'] = sqrt(sum(v ** 2 for v in totSystTot.values()))

    toFormat['totObs'] = sum(nObs.values())


    print table.format(**toFormat) # for AN
    print ''
    print table.replace(':.2f',':.1f').format(**toFormat) # for PAS/paper


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
    parser.add_argument('--sfHists', action='store_true',
                        help='Get lepton scale factors from files instead of directly from the ntuples.')

    args=parser.parse_args()

    cardName = _join(defaultPath, args.cardName)

    main(args.dataDir, args.mcDir, args.analysis, cardName, args.fakeRateFile,
         args.puWeightFile, args.lumi, args.sfHists)
