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


def leptonEfficiencyWeights(channel, systematic=''):
    channels = _parseChannels(channel)

    sfTemp = '{lep}EffScaleFactor'
    if systematic.lower() == 'up':
        sfTemp = '({0} + {0}Error)'.format(sfTemp)
    elif systematic.lower() in ['dn', 'down']:
        sfTemp = '({0} - {0}Error)'.format(sfTemp)
    elif systematic:
        raise ValueError("Unknown lepton efficiency systematic {}".format(systematic))

    out = {
        c : ' * '.join([sfTemp.format(lep=obj) for obj in _mapObjects(c)]) for c in channels
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


def baseMCWeight(channel, puWeightFile, lepSyst='', puSyst=''):
    lepWeights = leptonEfficiencyWeights(channel, lepSyst)
    puWtStr, puWtFun = puWeight(puWeightFile, puSyst)

    if isinstance(lepWeights, str):
        return '{} * {}'.format(lepWeights, puWtStr)

    return {
        c : '{} * {}'.format(w,puWtStr) for c,w in lepWeights.iteritems()
        }
