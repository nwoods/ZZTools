'''

Helpers to get weight strings, possibly with systematic shifts.

Nate Woods, U. Wisconsin

'''

from Utilities.helpers import parseChannels as _parseChannels
from Utilities.helpers import mapObjects as _mapObjects
from Utilities import WeightStringMaker as _Weight

from rootpy import asrootpy as _asRP
from rootpy.io import root_open as _open

from os import environ as _env
from os import path as _path


def leptonEfficiencyWeights(channel, eSystematic='', mSystematic=''):
    channels = _parseChannels(channel)

    sfTemp = {'e':'{lep}EffScaleFactor','m':'{lep}EffScaleFactor'}
    if eSystematic.lower() == 'up':
        sfTemp['e'] = '({0} + {0}Error)'.format(sfTemp['e'])
    elif eSystematic.lower() in ['dn', 'down']:
        sfTemp['e'] = '({0} - {0}Error)'.format(sfTemp['e'])
    elif eSystematic:
        raise ValueError("Unknown electron efficiency systematic {}".format(eSystematic))
    if mSystematic.lower() == 'up':
        sfTemp['m'] = '({0} + {0}Error)'.format(sfTemp['m'])
    elif mSystematic.lower() in ['dn', 'down']:
        sfTemp['m'] = '({0} - {0}Error)'.format(sfTemp['m'])
    elif mSystematic:
        raise ValueError("Unknown muon efficiency systematic {}".format(mSystematic))

    out = {
        c : ' * '.join([sfTemp[obj[0]].format(lep=obj) for obj in _mapObjects(c)]) for c in channels
        }

    if len(channels) == 1:
        return out[channels[0]]

    return out


_sfStrings = {'e':{},'m':{}}

def leptonEfficiencyWeightsFromHists(channel, eSystematic='', mSystematic='',
                                     eSelSFFile='eleSelectionSF_HZZ_Moriond17',
                                     eSelSFFileGap='',
                                     eRecoSFFile='eleRecoSF_HZZ_Moriond17',
                                     mSFFile='muSelectionAndRecoSF_HZZ_Moriond17'):
    global _sfStrings

    if not eSelSFFileGap:
        eSelSFFileGap = eSelSFFile.replace('SF', 'SFGap')

    if not eSelSFFile.endswith('.root'):
        eSelSFFile += '.root'
    if not eSelSFFileGap.endswith('.root'):
        eSelSFFileGap += '.root'
    if not eRecoSFFile.endswith('.root'):
        eRecoSFFile += '.root'
    if not mSFFile.endswith('.root'):
        mSFFile += '.root'

    channels = _parseChannels(channel)

    sfTemp = {}

    if any('e' in chan for chan in channels):
        if eSystematic not in _sfStrings['e']:
            with _open(_path.join(_env['zzt'],'data','leptonScaleFactors',
                                  eSelSFFile)) as fEleSel:
                hEleSel = _asRP(fEleSel.EGamma_SF2D).clone()
                hEleSel.SetDirectory(0)
            with _open(_path.join(_env['zzt'],'data','leptonScaleFactors',
                                  eSelSFFileGap)) as fEleSelGap:
                hEleSelGap = _asRP(fEleSelGap.EGamma_SF2D).clone()
                hEleSelGap.SetDirectory(0)
            with _open(_path.join(_env['zzt'],'data','leptonScaleFactors',
                                  eRecoSFFile)) as fEleReco:
                hEleReco = _asRP(fEleReco.EGamma_SF2D).clone()
                hEleReco.SetDirectory(0)

            wtMaker = _Weight('leptonSFs')

            eSelTemp = wtMaker.makeWeightStringFromHist(hEleSel, '{obj}SCEta', '{obj}Pt')
            eSelGapTemp = wtMaker.makeWeightStringFromHist(hEleSelGap, '{obj}SCEta', '{obj}Pt')
            eRecoTemp = wtMaker.makeWeightStringFromHist(hEleReco, '{obj}SCEta', '100.')
            eSFTemp = '({} * (({{obj}}IsGap)*{} + (!{{obj}}IsGap)*{}))'.format(eRecoTemp,
                                                                               eSelGapTemp,
                                                                               eSelTemp)
            if eSystematic:
                eSelErrTemp = wtMaker.makeWeightStringFromHist(hEleSel, '{obj}SCEta', '{obj}Pt', getError=True)
                eSelErrGapTemp = wtMaker.makeWeightStringFromHist(hEleSelGap, '{obj}SCEta', '{obj}Pt', getError=True)
                eRecoErrTemp = wtMaker.makeWeightStringFromHist(hEleReco, '{obj}SCEta', '100', getError=True)
                eSFErrTemp = ('sqrt(({} + ({{obj}}Pt < 20. || {{obj}}Pt > 75.)*0.01)^2 + '
                              '(({{obj}}IsGap)*{} + (!{{obj}}IsGap)*{})^2)').format(eRecoErrTemp,
                                                                                    eSelErrGapTemp,
                                                                                    eSelErrTemp)
                if eSystematic.lower() == 'up':
                    sign = '+'
                elif eSystematic.lower() in ['dn','down']:
                    sign = '-'
                else:
                    raise ValueError("Unknown electron efficiency systematic {}".format(eSystematic))

                eSFTemp = '({} {} {})'.format(eSFTemp, sign, eSFErrTemp)

            _sfStrings['e'][eSystematic] = eSFTemp

        sfTemp['e'] = _sfStrings['e'][eSystematic]

    if any('m' in chan for chan in channels):
        if mSystematic not in _sfStrings['m']:
            with _open(_path.join(_env['zzt'],'data','leptonScaleFactors',
                                  mSFFile)) as fMuSF:
                hMuSF = _asRP(fMuSF.FINAL).clone()
                hMuSF.SetDirectory(0)
                if mSystematic:
                    hMuSFErr = _asRP(fMuSF.ERROR).clone()
                    hMuSFErr.SetDirectory(0)

            wtMaker = _Weight('leptonSFs')

            muSFTemp = wtMaker.makeWeightStringFromHist(hMuSF, '{obj}Eta', '{obj}Pt')

            if mSystematic:
                muSFErrTemp = wtMaker.makeWeightStringFromHist(hMuSFErr,
                                                               '{obj}Eta',
                                                               '{obj}Pt')

                if mSystematic.lower() == 'up':
                    sign = '+'
                elif mSystematic.lower() in ['dn','down']:
                    sign = '-'
                else:
                    raise ValueError("Unknown muon efficiency systematic {}".format(mSystematic))

                muSFTemp = '({} {} {})'.format(muSFTemp, sign, muSFErrTemp)

            _sfStrings['m'][mSystematic] = muSFTemp

        sfTemp['m'] = _sfStrings['m'][mSystematic]

    out = {
        c : ' * '.join([sfTemp[obj[0]].format(obj=obj) for obj in _mapObjects(c)]) for c in channels
        }

    if len(channels) == 1:
        return out[channels[0]]

    return out


_puFuns = {}
def puWeight(weightFile, systematic=''):
    global _puFuns

    try:
        return _puFuns[weightFile + '_' + systematic]
    except KeyError:
        pass

    weightFileName = weightFile
    if '.root' not in weightFile:
        weightFileName = weightFile + '.root'
    weightFileName = _path.join(_env['zzt'], 'data', 'pileup', weightFileName)

    sfName = 'puScaleFactor'
    if systematic.lower() == 'up':
        sfName += '_up'
    elif systematic.lower() in ['dn', 'down']:
        sfName += '_down'
    elif systematic:
        raise ValueError("Unknown pileup systematic {}".format(systematic))

    puWeight = _Weight('puWeight')
    with _open(weightFileName) as f:
        wtStr = puWeight.makeWeightStringFromHist(f.Get(sfName), 'nTruePU')
    wtFunc = puWeight.getWeightFunction(-1)
    _puFuns[weightFile + '_' + systematic] = wtStr, wtFunc
    return wtStr, wtFunc


def baseMCWeight(channel, puWeightFile, eSyst='', mSyst='', puSyst='',
                 **kwargs):
    '''
    kwargs can be used to specify files from which to retrieve scale factor
    histograms (instead of using the ones in the ntuples). Options are:
        eSelSFFile : non-gap electron selections
        eSelSFFileGap : gap electron selections
        eRecoSFFile : electron GSF track reconstruction
        mSFFile : muon overall efficiency
        scaleFactorsFromHists : if False, SF hists are not used, if True
            SF hists are definitely used, possibly default ones (for
            backwards compatibility)
    If the full path is not specified, they're assumed to be in
    ZZTools/data/leptonScaleFactors. The '.root' suffix is optional.
    If some but not all are specified, the others use defaults except
    eSelSFFileGap which is guessed from eSelSFFile. If none are specified,
    everything is taken from the ntuple rows as usual.
    '''
    useSFHists = kwargs.pop('scaleFactorsFromHists', bool(kwargs))
    if useSFHists:
        lepWeights = leptonEfficiencyWeightsFromHists(channel, eSyst, mSyst,
                                                      **kwargs)
    else:
        lepWeights = leptonEfficiencyWeights(channel, eSyst, mSyst)
    puWtStr, puWtFun = puWeight(puWeightFile, puSyst)

    if isinstance(lepWeights, str):
        return '{} * {}'.format(lepWeights, puWtStr)

    return {
        c : '{} * {}'.format(w,puWtStr) for c,w in lepWeights.iteritems()
        }
