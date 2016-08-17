'''

Some useful functions

Nate Woods, U. Wisconsin

'''



def makeNumberPretty(n, maxDigits=10):
    '''
    Take a number, return a string of it with the right number of digits
    and whatnot.
    Cuts off at maxDigits places after the decimal point.
    Assumes you want all digits before the decimal point no matter what.
    '''
    if int(n) == n: # integer
        return "%d"%n

    nDecimals = 0
    m = n
    while nDecimals < maxDigits:
        nDecimals += 1
        m *= 10
        if int(m) == m:
            break
    
    preFormat = "%%.%df"%nDecimals
    return preFormat%n


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
    goodWeights = [str(w) for w in wts if w]
    if not goodWeights:
        return '1'
    joiner = ' * '
    if kwargs.get('selections', False):
        joiner = ' && '
    joiner = '){}('.format(kwargs.get('joinwith', joiner))
    return '('+joiner.join(goodWeights)+')'

