
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



def standardZZSamples(dataDir, mcDir, resultType, puWeightFile, fakeRateFile,
                      lumi, **extraSamples):
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

    To put extra histograms in the stack, add extra keyword arguments of the 
    for sampleName='fileGlob*.root'. The files are assumed to be in the MC 
    directory.
    '''
    channels = ['eeee','eemm','mmmm']

    mcFileTemp = _path.join('/data/nawoods/ntuples', mcDir, # if mcDir is absolute, first argument is ignored
                            'results_{}'.format(resultType), '{}.root')
    mc2P2FTemp = _path.join('/data/nawoods/ntuples', mcDir, # if mcDir is absolute, first argument is ignored
                            'results_{}_2P2F'.format(resultType), '{}.root')
    mc3P1FTemp = _path.join('/data/nawoods/ntuples', mcDir, # if mcDir is absolute, first argument is ignored
                            'results_{}_3P1F'.format(resultType), '{}.root')

    qqZZByChan = {
        c : _MC('ZZTo4L', c, 
                     mcFileTemp.format('ZZTo4L_*'), 
                     True, lumi) for c in channels
        }
    qqZZByChan2P2F = {
        c : _MC('ZZTo4L', c, 
                     mc2P2FTemp.format('ZZTo4L_*'), 
                     True, lumi) for c in channels
        }
    qqZZByChan3P1F = {
        c : _MC('ZZTo4L', c, 
                     mc3P1FTemp.format('ZZTo4L_*'), 
                     True, lumi) for c in channels
        }

    ggZZByChan = {}
    ggZZByChan2P2F = {}
    ggZZByChan3P1F = {}
    for c in channels:
        ggZZByFS = {
            fs : _MC('GluGluZZTo{}'.format(fs), c, 
                          mcFileTemp.format('GluGluZZTo{}*'.format(fs)),
                          True, lumi) for fs in ['4e', '4mu', '2e2mu']
            }
        ggZZByChan[c] = _Group('GluGluZZ', c, ggZZByFS, True)
        ggZZByFS2P2F = {
            fs : _MC('GluGluZZTo{}'.format(fs), c, 
                          mc2P2FTemp.format('GluGluZZTo{}*'.format(fs)),
                          True, lumi) for fs in ['4e', '4mu', '2e2mu']
            }
        ggZZByChan2P2F[c] = _Group('GluGluZZ', c, ggZZByFS2P2F, True)
        ggZZByFS3P1F = {
            fs : _MC('GluGluZZTo{}'.format(fs), c, 
                          mc3P1FTemp.format('GluGluZZTo{}*'.format(fs)),
                          True, lumi) for fs in ['4e', '4mu', '2e2mu']
            }
        ggZZByChan3P1F[c] = _Group('GluGluZZ', c, ggZZByFS3P1F, True)

    dataFileTemp = _path.join('/data/nawoods/ntuples', dataDir, # if mcDir is absolute, first argument is ignored
                              'results_{}'.format(resultType), 'Run2016{}_*.root')
    data2P2FTemp = _path.join('/data/nawoods/ntuples', dataDir, # if mcDir is absolute, first argument is ignored
                              'results_{}_2P2F'.format(resultType), 'Run2016{}_*.root')
    data3P1FTemp = _path.join('/data/nawoods/ntuples', dataDir, # if mcDir is absolute, first argument is ignored
                              'results_{}_3P1F'.format(resultType), 'Run2016{}_*.root')

    dataByChan = {}
    dataByChan2P2F = {}
    dataByChan3P1F = {}
    for c in channels:
        samplesByEra = {}
        samplesByEra2P2F = {}
        samplesByEra3P1F = {}
        for era in ['B','C','D']:
            samplesByEra['2016{}'.format(era)] = _Data('data2016{}_{}'.format(era, c), c, 
                                                            dataFileTemp.format(era)
                                                            )
            samplesByEra2P2F['2016{}'.format(era)] = _Data('data2016{}_{}'.format(era, c), c, 
                                                                data2P2FTemp.format(era)
                                                                )
            samplesByEra3P1F['2016{}'.format(era)] = _Data('data2016{}_{}'.format(era, c), c, 
                                                                data3P1FTemp.format(era)
                                                                )

        dataByChan[c] = _Group('data_{}'.format(c), c, samplesByEra)
        dataByChan2P2F[c] = _Group('data_{}'.format(c), c, samplesByEra2P2F)
        dataByChan3P1F[c] = _Group('data_{}'.format(c), c, samplesByEra3P1F)

    if '.root' not in puWeightFile:
        puWeightFile = puWeightFile + '.root'
    puWeight = _Weight('puWeight')
    with _open(_path.join(_env['zzt'], 'data', 'pileup', 
                          puWeightFile)) as f:
        strPU = puWeight.makeWeightStringFromHist(f.puScaleFactor, 'nTruePU')
    
    mcWeight = {
        'eeee' : 'e1EffScaleFactor * e2EffScaleFactor * e3EffScaleFactor * e4EffScaleFactor * {}'.format(strPU),
        'eemm' : 'e1EffScaleFactor * e2EffScaleFactor * m1EffScaleFactor * m2EffScaleFactor * {}'.format(strPU),
        'mmmm' : 'm1EffScaleFactor * m2EffScaleFactor * m3EffScaleFactor * m4EffScaleFactor * {}'.format(strPU),
        }

    otherSamplesByChan = {
        name : {
            c : _MC('ZZTo4L', c, 
                    mcFileTemp.format(files), 
                    True, lumi) for c in channels
            } for name, files in extraSamples.iteritems()
        }
    
    if  '.root' not in fakeRateFile:
        fakeRateFile = fakeRateFile + '.root'
    wCR = _Weight('fakeFactor')
    with _open(_path.join(_env['zzt'], 'data', 'fakeRate',
                          fakeRateFile)) as f:
        fakeFactorStrE = wCR.makeWeightStringFromHist(f.fakeFactor_e, 'abs({lep}Eta)', '{lep}Pt')
        fakeFactorStrM = wCR.makeWeightStringFromHist(f.fakeFactor_m, 'abs({lep}Eta)', '{lep}Pt')
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

    for c in channels:
        # all MC samples need standard MC weight
        qqZZByChan[c].applyWeight(mcWeight[c])
        qqZZByChan2P2F[c].applyWeight(mcWeight[c])
        qqZZByChan3P1F[c].applyWeight(mcWeight[c])
        ggZZByChan[c].applyWeight(mcWeight[c])
        ggZZByChan2P2F[c].applyWeight(mcWeight[c])
        ggZZByChan3P1F[c].applyWeight(mcWeight[c])
    
        for s in otherSamplesByChan:
            otherSamplesByChan[s][c].applyWeight(mcWeight[c])

        # CR samples weighted by fake factor
        # 2P2F is subtracted from 3P1F so weight it by -1
        # ... but MC CRs are subtracted from data CRs, so give them all opposite sign
        dataByChan2P2F[c].applyWeight(crWeight[c])
        dataByChan2P2F[c].applyWeight('-1.')
        dataByChan3P1F[c].applyWeight(crWeight[c])
    
        qqZZByChan2P2F[c].applyWeight(crWeight[c])
        qqZZByChan3P1F[c].applyWeight(crWeight[c])
        qqZZByChan3P1F[c].applyWeight('-1.')
        ggZZByChan2P2F[c].applyWeight(crWeight[c])
        ggZZByChan3P1F[c].applyWeight(crWeight[c])
        ggZZByChan3P1F[c].applyWeight('-1.')

        print c
        print qqZZByChan[c].weight
        for n,s in ggZZByChan[c].itersamples():
            print n, s.weight
        print ''

    ### Make the stack and data points we'll actually use
    
    data = _Group('Data', 'zz', dataByChan)
    data.format(color='black',drawstyle='PE',legendstyle='LPE')
    
    qqZZ = _Group('ZZTo4L', 'zz', qqZZByChan, True)
    ggZZ = _Group('GluGluZZ', 'zz', ggZZByChan, True)
    otherMC = [_Group(n, 'zz', s, True) for n,s in otherSamplesByChan.iteritems()]

    bkgByChan = {
        c : _Group('Z+X', c, {'2P2F':dataByChan2P2F[c], 
                              '3P1F':dataByChan3P1F[c], 
                              'qqZZ2P2F':qqZZByChan2P2F[c], 
                              'qqZZ3P1F':qqZZByChan3P1F[c], 
                              'ggZZ2P2F':ggZZByChan2P2F[c], 
                              'ggZZ3P1F':ggZZByChan3P1F[c]
                              }
                        ) for c in channels
        }
    
    bkg = _Group('Z+X', 'zz', bkgByChan, True)
    
    # don't allow negative backgrounds
    bkg.setPostprocessor(_ensureNonneg)
    
    stack = _Stack('stack', 'zz', [bkg, ggZZ, qqZZ]+otherMC)

    return data, stack



def genZZSamples(fileDir, resultType, lumi):
    channels = ['eeee','eemm','mmmm']

    fileTemp = _path.join('/data/nawoods/ntuples', fileDir, # if fileDir is absolute, first argument is ignored
                          'results_{}'.format(resultType), '{}.root')

    samplesByChan = {c:{} for c in channels}

    for c in channels:
        theseSamples = {}
        theseSamples['ZZTo4L'] = _MC('ZZTo4L', c+'Gen', 
                                     fileTemp.format('ZZTo4L_*'), 
                                     True, lumi)

        for fs in ['4e', '4mu', '2e2mu']:
            sample = 'GluGluZZTo{}'.format(fs)
            theseSamples[sample] = _MC(sample, c+'Gen', 
                                       fileTemp.format('{}*'.format(sample)),
                                       True, lumi)

        samplesByChan[c] = _Group("GenZZ", c+'Gen', theseSamples, False)
    
    out = _Group('GenZZ', 'zzGen', samplesByChan, False)

    return out
