'''

Some useful functions

Nate Woods, U. Wisconsin

'''


import logging
from rootpy import log as rlog; rlog = rlog["/helpers"]
# don't show most silly ROOT messages
logging.basicConfig(level=logging.WARNING)
rlog["/rootpy.compiled"].setLevel(rlog.WARNING)

import rootpy.compiled as _rootComp
from rootpy.ROOT import TGraphAsymmErrors as _TGAE

from numbers import Number as _NumType

def makeNumberPretty(n, maxDigits=10):
    '''
    Take a number, return a string of it with the right number of digits
    and whatnot.
    Cuts off at maxDigits places after the decimal point.
    No trailing zeros, or decimal points if there are no digits after it.
    '''
    if maxDigits < 0:
        raise ValueError("Invalid number of digits: {}".format(maxDigits))

    digitsBeforeDecimal = len('{:.0f}'.format(n))
    digits = digitsBeforeDecimal + maxDigits
    return '{{:.{}g}}'.format(digits).format(n)


def parseChannels(channels):
    '''
    Take a string or list of strings and return a list of channels.
    '''
    if type(channels) == str:
        if channels in ['4l', 'zz', 'ZZ']:
            return ['eeee', 'eemm', 'mmmm']
        elif channels in ['3l', 'zl', 'Zl', 'z+l', 'Z+l']:
            return ['eee', 'eem', 'emm', 'mmm']
        elif channels in ['z', '2l', 'll', 'Z']:
            return ['ee', 'mm']
        else:
            chanList = channels.split(',')
            assert all(all(letter in ['e','m','t','g','j'] for letter in ch) and len(ch) <= 4 for ch in chanList),\
                'Invalid channel ' + channels
            return chanList
    else:
        assert isinstance(channels, list), 'Channels must be a list or a string'
        out = []
        for ch in channels:
            out += parseChannels(ch)
        return out


_zzhelpers_object_maps_ = {}
def mapObjects(channel):
    '''
    Return a list of objects of the form ['e1','e2','m1','m2'] or ['e1','e2','m']
    Objects are in alphabetical/numerical order order
    '''
    global _zzhelpers_object_maps_

    try:
        return _zzhelpers_object_maps_[channel]
    except KeyError:
        nObjects = {}
        objects = []

        for obj in channel:
            if obj not in nObjects:
                nObjects[obj] = 1
            else:
                nObjects[obj] += 1

        for obj, num in nObjects.iteritems():
            if num == 1:
                objects.append(obj)
            else:
                for i in range(num):
                    objects.append(obj+str(i+1))

        objects.sort()

        _zzhelpers_object_maps_[channel] = objects
        return objects


def dictFromJSONFile(fName):
    '''
    Take a JSON file, return a dict with the information.
    '''
    with open(fName) as f:
        return loadJSON(f)


def combineWeights(*wts, **kwargs):
    '''
    Combine all non-null items in wts into a string that multiplies them all
    together (' * ' between items), unless certain keyword arguments
    are present:
    kwargs['selections'] evaluates to True -> join with ' && ' instead of ' * '
    kwargs['joinwith'] can be used to replace ' * ' with an arbitrary string
    '''
    if not all(isinstance(w,str) or isinstance(w,_NumType) for w in wts):
        raise TypeError("You can only combine weights made of numbers or strings")
    goodWeights = [str(w) for w in wts if w]
    if not goodWeights:
        return '1'
    joiner = ' * '
    if kwargs.get('selections', False):
        joiner = ' && '
    joiner = '){}('.format(kwargs.get('joinwith', joiner))
    return '('+joiner.join(goodWeights)+')'


identityFunction = lambda *args, **kwargs: True


Z_MASS = 91.1876


_strDeltaPhi = ''
_fDeltaPhi = None
_strDeltaR = ''
_fDeltaR = None

def _setupDeltaR():
    global _strDeltaPhi
    global _fDeltaPhi
    global _strDeltaR
    global _fDeltaR

    if _strDeltaR and _fDeltaR is not None:
        return

    _rootComp.register_code(
        '''
        #include<cmath> // std::sqrt

        float _f_Delta_Phi(float phi1, float phi2)
        {
          const float pi = 3.14159265;
          float out = phi1 - phi2;
          while(out > pi)
            out -= 2. * pi;
          while(out < -1. * pi)
            out += 2. * pi;

          return out;
        }

        float _f_Delta_R(float eta1, float phi1, float eta2, float phi2)
        {
          float dPhi = _f_Delta_Phi(phi1, phi2);
          float dEta = eta1 - eta2;
          return std::sqrt(dPhi * dPhi + dEta * dEta);
        }
        ''', ['_f_Delta_Phi', '_f_Delta_R'])

    _strDeltaPhi = '_f_Delta_Phi'
    _fDeltaPhi = _rootComp._f_Delta_Phi
    _strDeltaR = '_f_Delta_R'
    _fDeltaR = _rootComp._f_Delta_R


def deltaPhiFunction():
    if _fDeltaPhi is None:
        _setupDeltaR()

    return _fDeltaPhi


def deltaPhiString():
    if not _strDeltaPhi:
        _setupDeltaR()

    return _strDeltaPhi


def deltaRFunction():
    if _fDeltaR is None:
        _setupDeltaR()

    return _fDeltaR


def deltaRString():
    if not _strDeltaR:
        _setupDeltaR()

    return _strDeltaR


def zMassDist(m):
    return abs(m - Z_MASS)


def zeroNegativeBins(h):
    for b in h:
        b.value = max(0., b.value)


def removeXErrors(p):
    '''
    Set all x-axis errors on a TGraphAsymmErrors to 0, or set a TH1 to not
    draw its x errors.
    '''
    try:
        for i in xrange(len(p)):
            p.SetPointEXhigh(i,0)
            p.SetPointEXlow(i,0)
    except: # histogram
        p.drawstyle = p.drawstyle + 'X0'

    if 'l' in p.legendstyle.lower():
        p.legendstyle = p.legendstyle.replace('L','').replace('l','')

