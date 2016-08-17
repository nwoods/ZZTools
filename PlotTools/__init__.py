from PlotStyle import PlotStyle

from rootpy.plotting import Legend as _Legend
from rootpy.plotting import HistStack as _HistStack


_defaultLegParams = {
    'entryheight' : 0.03,
    'entrysep' : 0.01,
    'leftmargin' : 0.5,
    'topmargin' : 0.05,
    'rightmargin' : 0.05,
    'textsize' : 0.033,
    }
def makeLegend(pad, *objects, **params):
    '''
    Make a legend initialized with parameters params, containing objects, 
    on pad.
    '''
    legParams = _defaultLegParams.copy()
    legParams.update(params)

    obs = []
    for ob in objects:
        if isinstance(ob, _HistStack):
            for h in ob:
                if h.Integral() > 0.:
                    obs.append(h)
        else:
            obs.append(ob)

    return _Legend(obs[::-1], pad, **legParams)

