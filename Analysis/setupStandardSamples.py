
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
                 efficiencySyst='', puSyst=''):
    channels = _parseChannels(channel)

    mcWeight = _baseMCWeight(channel, puWeightFile, efficiencySyst, puSyst)

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


def standardZZData(channel, inDir, resultType):
    channels = _parseChannels(channel)
    dataFileTemp = _path.join('/data/nawoods/ntuples', inDir, # if mcDir is absolute, first argument is ignored
                              'results_{}'.format(resultType), 'Run2016{}_*.root')

    byChan = {}
    for c in channels:
        samplesByEra = {}
        for era in ['B','C','D','E']:
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


def zzStackMCOnly(channel, inDir, resultType, puWeightFile, lumi, 
                  efficiencySyst='', puSyst='', amcatnlo=False, 
                  **extraSamples):
    qqZZSampleName = 'ZZTo4L'
    if amcatnlo:
        qqZZSampleName += '-amcatnlo'

    channels = _parseChannels(channel)

    qqZZByChan = {
        c : standardZZMC(c, inDir, qqZZSampleName, resultType, puWeightFile,
                         lumi, efficiencySyst, puSyst) for c in channels
        }

    ggZZByChan = {}
    for c in channels:
        ggZZByFS = {
            fs : standardZZMC(c, inDir, 'GluGluZZTo{}'.format(fs), 
                              resultType, puWeightFile, lumi, 
                              efficiencySyst, 
                              puSyst) for fs in ['4e', '4mu', '2e2mu']
            }
        ggZZByChan[c] = _Group('GluGluZZ', c, ggZZByFS, True)

    otherSamplesByChan = {
        name : {
            c : standardZZMC(c, inDir, name, resultType, puWeightFile, lumi,
                             efficiencySyst, puSyst) for c in channels
            } for name, files in extraSamples.iteritems()
        }

    if len(channels) == 1:
        qqZZ = qqZZByChan[channels[0]]
        ggZZ = ggZZByChan[channels[0]]
        otherMC = [s[channels[0]] for s in otherSamplesByChan.values()]
    else:
        qqZZ = _Group(qqZZSampleName, channel, qqZZByChan, True)
        ggZZ = _Group('GluGluZZ', channel, ggZZByChan, True)
        otherMC = [_Group(n, channel, s, True) for n,s in otherSamplesByChan.iteritems()]

    return _Stack('stack', channel, [qqZZ, ggZZ]+otherMC)


def standardZZBkg(channel, dataDir, mcDir, resultType, puWeightFile,
                  fakeRateFile, lumi, efficiencySyst='', puSyst='', 
                  fakeRateSyst='', **extraSamples):
    channels = _parseChannels(channel)

    data2P2F = standardZZData(channel, dataDir, resultType+'_2P2F')
    data3P1F = standardZZData(channel, dataDir, resultType+'_3P1F')
    mc2P2F = zzStackMCOnly(channel, mcDir, resultType+'_2P2F', puWeightFile, 
                           lumi, efficiencySyst, puSyst)
    mc3P1F = zzStackMCOnly(channel, mcDir, resultType+'_3P1F', puWeightFile, 
                           lumi, efficiencySyst, puSyst)

    if  '.root' not in fakeRateFile:
        fakeRateFile = fakeRateFile + '.root'
    wCR = _Weight('fakeFactor')
    with _open(_path.join(_env['zzt'], 'data', 'fakeRate',
                          fakeRateFile)) as f:
        fakeFactorStrE = wCR.makeWeightStringFromHist(f.fakeFactor_e, 'abs({lep}Eta)', '{lep}Pt')
        fakeFactorStrM = wCR.makeWeightStringFromHist(f.fakeFactor_m, 'abs({lep}Eta)', '{lep}Pt')
        if fakeRateSyst.lower() == 'up':
            fakeFactorStrE = '1.4 * {}'.format(fakeFactorStrE)
            fakeFactorStrM = '1.4 * {}'.format(fakeFactorStrM)
        elif fakeRateSyst.lower() in ['dn', 'down']:
            fakeFactorStrE = '0.6 * {}'.format(fakeFactorStrE)
            fakeFactorStrM = '0.6 * {}'.format(fakeFactorStrM)

    # CR samples weighted by fake factor
    # 2P2F is subtracted from 3P1F so weight it by -1
    # ... but MC CRs are subtracted from data CRs, so give them all opposite sign
    data2P2F.applyWeight('-1.')
    mc3P1F.applyWeight('-1.')

    zCRWeightTemp = ('({{lep1}}ZZTightID && {{lep1}}ZZIsoPass ? 1. : {fr1}) * '
                     '({{lep2}}ZZTightID && {{lep2}}ZZIsoPass ? 1. : {fr2})')
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
    
    # don't allow negative backgrounds
    bkg.setPostprocessor(_ensureNonneg)

    return bkg


def standardZZStack(channel, dataDir, mcDir, resultType, puWeightFile, 
                    fakeRateFile, lumi, efficiencySyst='', puSyst='', 
                    fakeRateSyst='', amcatnlo=False, **extraSamples):
    stack = zzStackMCOnly(channel, mcDir, resultType, puWeightFile, lumi, 
                          efficiencySyst, puSyst, amcatnlo)
    bkg = standardZZBkg(channel, dataDir, mcDir, resultType, puWeightFile, 
                        fakeRateFile, lumi, efficiencySyst, puSyst, fakeRateSyst)

    stack.addSample(bkg)
    return stack


def standardZZSamples(channel, dataDir, mcDir, resultType, puWeightFile, 
                      fakeRateFile, lumi, efficiencySyst='', puSyst='', 
                      fakeRateSyst='', amcatnlo=False, **extraSamples):
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

    If applyWeights is False, PU weights and data/MC scale factors are not used

    efficiencySyst, puSyst, and fakeRateSyst give systematic shifts for those
    things with 'up' or 'down'.

    If amcatnlo is True, the aMC@NLO qq->ZZ sample will be used instead of the
    standard POWHEG sample.

    To put extra histograms in the stack, add extra keyword arguments of the 
    for sampleName='fileGlob*.root'. The files are assumed to be in the MC 
    directory.
    '''
    data = standardZZData(channel, dataDir, resultType)

    stack = standardZZStack(channel, dataDir, mcDir, resultType, 
                            puWeightFile, fakeRateFile,
                            lumi, efficiencySyst, puSyst, fakeRateSyst, 
                            amcatnlo, **extraSamples)

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


def genZZSamples(channel, fileDir, resultType, lumi, amcatnlo=False):
    qqZZSampleName = 'ZZTo4L'
    if amcatnlo:
        qqZZSampleName += '-amcatnlo'

    channels = _parseChannels(channel)

    samplesByChan = {c:{} for c in channels}

    for c in channels:
        theseSamples = {}
        theseSamples[qqZZSampleName] = standardZZGen(c, fileDir,  
                                                     qqZZSampleName, 'smp', 
                                                     lumi)

        for fs in ['4e', '4mu', '2e2mu']:
            sample = 'GluGluZZTo{}'.format(fs)
            theseSamples[sample] = standardZZGen(c, fileDir, sample, 'smp',
                                                 lumi)

        samplesByChan[c] = _Group("GenZZ", c+'Gen', theseSamples, False)
    
    if len(channels) == 1:
        return samplesByChan[channels[0]]

    return _Group('GenZZ', channel+'Gen', samplesByChan, False)

