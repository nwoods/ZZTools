'''

Helpers to get weight strings, possibly with systematic shifts.

Nate Woods, U. Wisconsin

'''

from Utilities.helpers import parseChannels as _parseChannels
from Utilities.helpers import mapObjects as _mapObjects
from Utilities import WeightStringMaker as _Weight

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


def baseMCWeight(channel, puWeightFile, eSyst='', mSyst='', puSyst=''):
    lepWeights = leptonEfficiencyWeights(channel, eSyst, mSyst)
    puWtStr, puWtFun = puWeight(puWeightFile, puSyst)

    if isinstance(lepWeights, str):
        return '{} * {}'.format(lepWeights, puWtStr)

    return {
        c : '{} * {}'.format(w,puWtStr) for c,w in lepWeights.iteritems()
        }
