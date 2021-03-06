
import logging as _logging
from rootpy import log as _log; _log = _log["/setupStandardSamples"]
# don't show most silly ROOT messages
_logging.basicConfig(level=_logging.WARNING)
_log["/ROOT.TUnixSystem.SetDisplay"].setLevel(_log.ERROR)

from SampleTools import MCSample as _MC
from SampleTools import DataSample as _Data
from SampleTools import SampleGroup as _Group
from SampleTools import SampleStack as _Stack
from Utilities import WeightStringMaker as _Weight
from Analysis.weightHelpers import leptonEfficiencyWeights as _leptonEfficiencyWeights
from Analysis.weightHelpers import puWeight as _puWeight
from Analysis.weightHelpers import baseMCWeight as _baseMCWeight
from Utilities.helpers import parseChannels as _parseChannels
from Utilities.helpers import mapObjects as _mapObjects

from rootpy.io import root_open as _open

from os import environ as _env
from os import path as _path
from collections import OrderedDict as _ODict
from re import compile as _reComp


def _ensureNonneg(h):
    '''
    Remove any bin of h with value < 0. Any nonnegative bin whose error bar
    would go below 0 has the error set to the bin value.
    '''
    for b in h:
        if b.value < 0.:
            b.value = 0.
            b.error = 0.
        elif b.error > b.value:
            b.error = b.value

    return h



def standardZZMC(channel, inDir, sampleName, resultType, puWeightFile, lumi,
                 eEfficiencySyst='', mEfficiencySyst='', puSyst='',
                 **kwargs):
    channels = _parseChannels(channel)

    mcWeight = _baseMCWeight(channel, puWeightFile, eEfficiencySyst,
                             mEfficiencySyst, puSyst,
                             **kwargs)

    mcFiles = _path.join('/data/nawoods/ntuples', inDir, # if inDir is absolute, first argument is ignored
                         'results_{}'.format(resultType),
                         '{}_*.root'.format(sampleName))

    byChan = {
        c : _MC(sampleName, c,
                mcFiles, True, lumi) for c in channels
        }

    if len(channels) == 1:
        mc = _MC(sampleName, channel, mcFiles, True, lumi)
        mc.applyWeight(mcWeight)
    else:
        mc = _Group(sampleName, channel, byChan, True)
        mc.applyWeight(mcWeight)

    return mc


def standardZZData(channel, inDir, resultType, eras='BCDEFGH'):
    channels = _parseChannels(channel)
    dataFileTemp = _path.join('/data/nawoods/ntuples', inDir, # if mcDir is absolute, first argument is ignored
                              'results_{}'.format(resultType), 'Run2016{}_*.root')
    byChan = {}
    for c in channels:
        samplesByEra = {}
        for era in eras:
            samplesByEra['2016{}'.format(era)] = _Data('2016{}_{}'.format(era, c),
                                                       c, dataFileTemp.format(era)
                                                       )
        byChan[c] = _Group('Data', c, samplesByEra)
        byChan[c].format(color='black',drawstyle='PE',legendstyle='LPE')

    if len(channels) == 1:
        return byChan[channels[0]]

    data = _Group('Data', channel, byChan)
    data.format(color='black',drawstyle='PE',legendstyle='LPE')

    return data


def zzStackSignalOnly(channel, inDir, resultType, puWeightFile, lumi,
                      eEfficiencySyst='', mEfficiencySyst='', puSyst='',
                      amcatnlo=False, higgs=False, asGroup=False,
                      skipEWK=False, madgraphEWK=False, jetBinnedZZ=False,
                      *extraSamples, **kwargs):
    channels = _parseChannels(channel)

    samplesByChan = _ODict()

    if not jetBinnedZZ:
        qqZZSampleName = 'ZZTo4L'
        if amcatnlo:
            qqZZSampleName += '-amcatnlo'

        # qq->ZZ
        samplesByChan[qqZZSampleName] = {
            c : standardZZMC(c, inDir, qqZZSampleName, resultType, puWeightFile,
                             lumi, eEfficiencySyst, mEfficiencySyst,
                             puSyst,
                             **kwargs) for c in channels
            }
    else:
        qqZZByChan = {}
        for c in channels:
            qqZZByJetBin = {}
            for nJ in 0, 1, 2:
                sName = 'ZZTo4L-{}J'.format(nJ)
                qqZZByJetBin[sName] = standardZZMC(c, inDir, sName, resultType,
                                                   puWeightFile, lumi,
                                                   eEfficiencySyst,
                                                   mEfficiencySyst,
                                                   puSyst,
                                                   **kwargs)
            qqZZByChan[c] = _Group('ZZTo4L-nJ', c, qqZZByJetBin, True)
        samplesByChan['ZZTo4L-nJ'] = qqZZByChan


    # gg->ZZ
    ggZZByChan = {}
    for c in channels:
        ggZZByFS = {
            fs : standardZZMC(c, inDir, 'GluGluZZTo{}'.format(fs),
                              resultType, puWeightFile, lumi,
                              eEfficiencySyst, mEfficiencySyst,
                              puSyst,
                              **kwargs) for fs in ['4e', '4mu', '2e2mu']
            }
        ggZZByChan[c] = _Group('GluGluZZ', c, ggZZByFS, True)
    samplesByChan['GluGluZZ'] = ggZZByChan

    if not skipEWK:
        # EWK ZZ+2jets (if it's not there, just skip it)
        try:
            if madgraphEWK:
                samplesByChan['ZZJJTo4L-EWK'] = {
                    c : standardZZMC(c, inDir, 'ZZJJTo4L-EWK', resultType,
                                     puWeightFile, lumi, eEfficiencySyst,
                                     mEfficiencySyst, puSyst,
                                     **kwargs) for c in channels
                }
            else:
                ewkByChan = {}
                for c in channels:
                    ewkByFS = {
                        fs : standardZZMC(c, inDir, 'ZZJJTo{}-EWK-phantom'.format(fs),
                                          resultType, puWeightFile, lumi,
                                          eEfficiencySyst, mEfficiencySyst,
                                          puSyst,
                                          **kwargs) for fs in ['4e', '4mu', '2e2mu']
                        }
                    ewkByChan[c] = _Group('ZZJJTo4L-EWK-phantom', c, ewkByFS, True)
                samplesByChan['ZZJJTo4L-EWK-phantom'] = ewkByChan
        except IOError:
            pass

    # Higgs if desired
    if higgs:
        samplesByChan['ggHZZ'] = {
            c : standardZZMC(c, inDir, 'ggHZZ', resultType, puWeightFile,
                             lumi, eEfficiencySyst, mEfficiencySyst,
                             puSyst,
                             **kwargs) for c in channels
            }

    # anything else
    for name in extraSamples:
        samplesByChan[name] = {
            c : standardZZMC(c, inDir, name, resultType, puWeightFile, lumi,
                             eEfficiencySyst, mEfficiencySyst,
                             puSyst,
                             **kwargs) for c in channels
            }

    if asGroup:
        if len(channels) == 1:
            samples = {n:s[channels[0]]
                       for n,s in samplesByChan.iteritems()}
        else:
            samples = {c:_Group('signal',c,{n:s[c] for n,s in samplesByChan.iteritems()})
                       for c in channels}
        return _Group('signal', channel, samples)

    if len(channels) == 1:
        samples = [s[channels[0]] for s in samplesByChan.values()]
    else:
        samples = [_Group(n, channel, s, True) for n,s in samplesByChan.iteritems()]

    return _Stack('stack', channel, samples)


def zzIrreducibleBkg(channel, inDir, resultType, puWeightFile, lumi,
                     eEfficiencySyst='', mEfficiencySyst='', puSyst='',
                     **kwargs):
    channels = _parseChannels(channel)

    samplesByChan = _ODict()

    samplesByChan['TTZ'] = {
        c : standardZZMC(c, inDir, 'TTZ', resultType, puWeightFile,
                         lumi, eEfficiencySyst, mEfficiencySyst,
                         puSyst,
                         **kwargs) for c in channels
        }

    samplesByChan['WWZ'] = {
        c : standardZZMC(c, inDir, 'WWZ', resultType, puWeightFile,
                         lumi, eEfficiencySyst, mEfficiencySyst,
                         puSyst,
                         **kwargs) for c in channels
        }

    groupByChan = {
        c : _Group('irreducible', c,
                   {n:samplesByChan[n][c] for n in samplesByChan},
                   True) for c in channels
        }

    if len(channels) == 1:
        return groupByChan[channels[0]]
    else:
        return _Group('irreducible', channel, groupByChan, True)


def standardZZBkg(channel, dataDir, mcDir, resultType, puWeightFile,
                  fakeRateFile, lumi, eEfficiencySyst='', mEfficiencySyst='',
                  puSyst='', eFakeRateSyst='', mFakeRateSyst='',
                  eras='BCDEFGH', sipCut=4., **kwargs):
    channels = _parseChannels(channel)

    data2P2F = standardZZData(channel, dataDir, resultType+'_2P2F',
                              eras=eras)
    data3P1F = standardZZData(channel, dataDir, resultType+'_3P1F',
                              eras=eras)
    mc2P2F = zzStackSignalOnly(channel, mcDir, resultType+'_2P2F', puWeightFile,
                               lumi, eEfficiencySyst, mEfficiencySyst, puSyst,
                               skipEWK=True, **kwargs)
    mc3P1F = zzStackSignalOnly(channel, mcDir, resultType+'_3P1F', puWeightFile,
                               lumi, eEfficiencySyst, mEfficiencySyst, puSyst,
                               skipEWK=True, **kwargs)

    if  '.root' not in fakeRateFile:
        fakeRateFile = fakeRateFile + '.root'
    wCR = _Weight('fakeFactor')
    with _open(_path.join(_env['zzt'], 'data', 'fakeRate',
                          fakeRateFile)) as f:
        fakeFactorStrE = wCR.makeWeightStringFromHist(f.fakeFactor_e, 'abs({lep}Eta)', '{lep}Pt')
        fakeFactorStrM = wCR.makeWeightStringFromHist(f.fakeFactor_m, 'abs({lep}Eta)', '{lep}Pt')
        if eFakeRateSyst.lower() == 'up':
            fakeFactorStrE = '1.4 * {}'.format(fakeFactorStrE)
        elif eFakeRateSyst.lower() in ['dn', 'down']:
            fakeFactorStrE = '0.6 * {}'.format(fakeFactorStrE)
        if mFakeRateSyst.lower() == 'up':
            fakeFactorStrM = '1.4 * {}'.format(fakeFactorStrM)
        elif mFakeRateSyst.lower() in ['dn', 'down']:
            fakeFactorStrM = '0.6 * {}'.format(fakeFactorStrM)

    # CR samples weighted by fake factor
    # 2P2F is subtracted from 3P1F so weight it by -1
    # ... but MC CRs are subtracted from data CRs, so give them all opposite sign
    data2P2F.applyWeight('-1.')
    mc3P1F.applyWeight('-1.')

    zCRWeightTemp = ('({{lep1}}ZZTightID && {{lep1}}ZZIsoPass ? 1. : {fr1}) * '
                     '({{lep2}}ZZTightID && {{lep2}}ZZIsoPass ? 1. : {fr2})')

    if abs(sipCut-4.) > 0.001:
        idPattern = _reComp(r'(?<=(?P<lep>\{\{lep\d\}\})ZZTightID)')
        idRepl = r'NoVtx && abs(\g<lep>PVDXY) < 0.5 && abs(\g<lep>PVDZ) < 1.'
        if sipCut > 0:
            idRepl += ' && \g<lep>SIP3D < {}'.format(sipCut)
        zCRWeightTemp = idPattern.sub(idRepl, zCRWeightTemp,2)

    zeCRWeight = zCRWeightTemp.format(fr1=fakeFactorStrE.format(lep='{lep1}'),
                                      fr2=fakeFactorStrE.format(lep='{lep2}'))
    zmCRWeight = zCRWeightTemp.format(fr1=fakeFactorStrM.format(lep='{lep1}'),
                                      fr2=fakeFactorStrM.format(lep='{lep2}'))

    crWeight = {
        'eeee' : zeCRWeight.format(lep1='e1',lep2='e2') + ' * ' + zeCRWeight.format(lep1='e3',lep2='e4'),
        'eemm' : zeCRWeight.format(lep1='e1',lep2='e2') + ' * ' + zmCRWeight.format(lep1='m1',lep2='m2'),
        'mmmm' : zmCRWeight.format(lep1='m1',lep2='m2') + ' * ' + zmCRWeight.format(lep1='m3',lep2='m4'),
        }

    if len(channels) == 1:
        crWt = crWeight[channels[0]]
    else:
        crWt = {c:crWeight[c] for c in channels}

    data2P2F.applyWeight(crWt)
    data3P1F.applyWeight(crWt)
    mc2P2F.applyWeight(crWt)
    mc3P1F.applyWeight(crWt)

    qqZZ2P2F = mc2P2F[0]
    qqZZ3P1F = mc3P1F[0]
    ggZZ2P2F = mc2P2F[1]
    ggZZ3P1F = mc3P1F[1]

    if len(channels) == 1:
        data2P2F = {channels[0] : data2P2F}
        data3P1F = {channels[0] : data3P1F}
        qqZZ2P2F = {channels[0] : qqZZ2P2F}
        qqZZ3P1F = {channels[0] : qqZZ3P1F}
        ggZZ2P2F = {channels[0] : ggZZ2P2F}
        ggZZ3P1F = {channels[0] : ggZZ3P1F}

    bkgByChan = {
        c : _Group('Z+X', c, {'2P2F':data2P2F[c],
                              '3P1F':data3P1F[c],
                              'qqZZ2P2F':qqZZ2P2F[c],
                              'qqZZ3P1F':qqZZ3P1F[c],
                              'ggZZ2P2F':ggZZ2P2F[c],
                              'ggZZ3P1F':ggZZ3P1F[c]
                              }
                        ) for c in channels
        }

    if len(channels) == 1:
        bkg = bkgByChan[channels[0]]
    else:
        bkg = _Group('Z+X', channel, bkgByChan, True)

    bkg.setPostprocessor(_ensureNonneg, True)

    return bkg


def zzStackBkgOnly(channel, dataDir, mcDir, resultType, puWeightFile,
                   fakeRateFile, lumi, eEfficiencySyst='', mEfficiencySyst='',
                   puSyst='', eFakeRateSyst='', mFakeRateSyst='',
                   amcatnlo=False, higgs=False, eras='BCDEFGH',
                   sipCut=4., *extraSamples, **kwargs):
    irreducible = zzIrreducibleBkg(channel, mcDir, resultType, puWeightFile,
                                   lumi, eEfficiencySyst, mEfficiencySyst,
                                   puSyst, **kwargs)
    reducible = standardZZBkg(channel, dataDir, mcDir, resultType, puWeightFile,
                              fakeRateFile, lumi, eEfficiencySyst,
                              mEfficiencySyst, puSyst, eFakeRateSyst,
                              mFakeRateSyst, eras=eras, sipCut=sipCut,
                              **kwargs)

    return _Stack('bkg', channel, [irreducible, reducible])


def standardZZStack(channel, dataDir, mcDir, resultType, puWeightFile,
                    fakeRateFile, lumi, eEfficiencySyst='', mEfficiencySyst='',
                    puSyst='', eFakeRateSyst='', mFakeRateSyst='',
                    amcatnlo=False, higgs=False, skipEWK=False,
                    madgraphEWK=False, eras='BCDEFGH', jetBinnedZZ=False,
                    sipCut=4., *extraSamples, **kwargs):
    stack = zzStackSignalOnly(channel, mcDir, resultType, puWeightFile, lumi,
                              eEfficiencySyst, mEfficiencySyst, puSyst, amcatnlo,
                              higgs=higgs, skipEWK=skipEWK,
                              madgraphEWK=madgraphEWK, jetBinnedZZ=jetBinnedZZ,
                              *extraSamples, **kwargs)
    irreducible = zzIrreducibleBkg(channel, mcDir, resultType, puWeightFile,
                                   lumi, eEfficiencySyst, mEfficiencySyst,
                                   puSyst, **kwargs)
    reducible = standardZZBkg(channel, dataDir, mcDir, resultType, puWeightFile,
                              fakeRateFile, lumi, eEfficiencySyst,
                              mEfficiencySyst, puSyst, eFakeRateSyst,
                              mFakeRateSyst, eras=eras, sipCut=sipCut,
                              **kwargs)

    stack.addSample(irreducible)
    stack.addSample(reducible)
    return stack


def standardZZSamples(channel, dataDir, mcDir, resultType, puWeightFile,
                      fakeRateFile, lumi, eEfficiencySyst='',
                      mEfficiencySyst='', puSyst='', eFakeRateSyst='',
                      mFakeRateSyst='', amcatnlo=False, higgs=False,
                      eras='BCDEFGH', skipEWK=False, madgraphEWK=False,
                      jetBinnedZZ=False, sipCut=4., *extraSamples, **kwargs):
    '''
    Return dataSampleGroup, bkgAndMCStack for data files in
    [dataDir]/results_[resultType] and MC in [mcDir]/results[resultType],
    using pileup weights from data/pileup/[puWeightFile].
    If the data and MC directories start with /data/nawoods/ntuples, this can
    be omitted.
    Background is constructed from samples in
    [dataDir]/results_[resultType]_2P2F and ..._3P1F, using fake factor from
    data/fakeRate/[puWeightFile].
    MC scaled to integrated luminosity lumi.

    efficiencySyst, puSyst, and fakeRateSyst give systematic shifts for those
    things with 'up' or 'down'.

    If amcatnlo is True, the aMC@NLO qq->ZZ sample will be used instead of the
    standard POWHEG sample. Ignored if jetBinnedZZ is True.

    if jetBinnedZZ is True, the jet-binned aMC@NLO qq->ZZ samples will be used
    instead of the standard POWHEG sample. Note that this sample contains
    only events with two on-shell Zs.

    To put extra histograms in the stack, add extra keyword arguments of the
    for sampleName='fileGlob*.root'. The files are assumed to be in the MC
    directory.

    kwargs are eventually passed to Analysis.weightHelpers.baseMCWeight,
    see it for options
    '''
    data = standardZZData(channel, dataDir, resultType, eras=eras)

    stack = standardZZStack(channel, dataDir, mcDir, resultType,
                            puWeightFile, fakeRateFile,
                            lumi, eEfficiencySyst, mEfficiencySyst, puSyst,
                            eFakeRateSyst, mFakeRateSyst, amcatnlo,
                            higgs=higgs, eras=eras,
                            skipEWK=skipEWK, madgraphEWK=madgraphEWK,
                            jetBinnedZZ=jetBinnedZZ, sipCut=sipCut,
                            *extraSamples, **kwargs)

    return data, stack


def standardZZGen(channel, inDir, sampleName, resultType, lumi):
    channels = _parseChannels(channel)

    files = _path.join('/data/nawoods/ntuples', inDir, # if inDir is absolute, first argument is ignored
                       'results_{}'.format(resultType),
                       '{}_*.root'.format(sampleName))

    byChan = {
        c : _MC(sampleName, c+'Gen',
                files, True, lumi) for c in channels
        }

    if len(channels) == 1:
        return _MC(sampleName, channel+'Gen', files, True, lumi)

    return _Group(sampleName, channel+'Gen', byChan, True)


def genZZSamples(channel, fileDir, resultType, lumi, amcatnlo=False,
                 higgs=False, madgraphEWK=False, jetBinnedZZ=False):
    qqZZSampleName = 'ZZTo4L'
    if amcatnlo:
        qqZZSampleName += '-amcatnlo'
    elif jetBinnedZZ:
        qqZZSampleName += '-{}J'

    channels = _parseChannels(channel)

    samplesByChan = {c:{} for c in channels}

    for c in channels:
        theseSamples = {}
        if jetBinnedZZ:
            for nj in 0, 1, 2:
                sample = qqZZSampleName.format(nj)
                theseSamples[sample] = standardZZGen(c, fileDir, sample,
                                                     resultType, lumi)
        else:
            theseSamples[qqZZSampleName] = standardZZGen(c, fileDir,
                                                         qqZZSampleName,
                                                         resultType,
                                                         lumi)

        for fs in ['4e', '4mu', '2e2mu']:
            sample = 'GluGluZZTo{}'.format(fs)
            theseSamples[sample] = standardZZGen(c, fileDir, sample,
                                                 resultType,
                                                 lumi)

        if madgraphEWK:
            theseSamples['ZZJJTo4L-EWK'] = standardZZGen(c, fileDir,
                                                         'ZZJJTo4L-EWK',
                                                         resultType,
                                                         lumi)
        else:
            for fs in ['4e', '4mu', '2e2mu']:
                sample = 'ZZJJTo{}-EWK-phantom'.format(fs)
                theseSamples[sample] = standardZZGen(c, fileDir, sample,
                                                     resultType,
                                                     lumi)

        if higgs:
            theseSamples['ggHZZ'] = standardZZGen(c, fileDir, 'ggHZZ',
                                                  resultType, lumi)

        samplesByChan[c] = _Group("GenZZ", c+'Gen', theseSamples, False)

    if len(channels) == 1:
        return samplesByChan[channels[0]]

    return _Group('GenZZ', channel+'Gen', samplesByChan, False)

