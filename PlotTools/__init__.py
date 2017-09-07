from PlotStyle import PlotStyle
from pdfViaTex import pdfViaTex

from rootpy.plotting import Legend as _Legend
from rootpy.plotting import HistStack as _HistStack
from rootpy.plotting import Pad as _Pad
from rootpy.plotting import Graph as _Graph
from rootpy.plotting.utils import get_band as _band
from rootpy.plotting.base import Color as _Color
from rootpy.ROOT import TLine as _Line, TAttFill as _Fill

from numbers import Number
from math import sqrt

_defaultLegParams = {
    'entryheight' : 0.03,
    'entrysep' : 0.01,
    'leftmargin' : 0.5,
    'topmargin' : 0.08,
    'rightmargin' : 0.05,
    'textsize' : 0.033,
    }
def makeLegend(pad, *objects, **params):
    '''
    Make a legend initialized with parameters params, containing objects,
    on pad.
    If 'solid' is in params and evaluates True, the legend will be opaque.
    '''
    legParams = _defaultLegParams.copy()
    solid = params.pop('solid', False)
    legParams.update(params)

    obs = []
    for ob in objects:
        if isinstance(ob, _HistStack):
            for h in ob:
                if h.Integral() > 0.:
                    obs.append(h)
        else:
            obs.append(ob)

    out = _Legend(obs[::-1], pad, **legParams)
    if solid:
        out.SetFillStyle(1001)

    return out

def addPadsBelow(p, *heights, **kwargs):
    '''
    Split pad p into multiple pads, and return a tuple containing all pads
    (ordered from top of the canvas to bottom).

    p (Pad): Pad or Canvas to split.
    heights (float(s)): heights of the lower pads, as a fraction of the canvas
        height, ordered top to bottom. Should not include top pad.

    Recognized keyword arguments:
    bottomMargin (float, default=0.3): Height of the margin at the bottom of
        the pad (where the x-axis title goes), as a fraction of the height of
        the last pad.
    topMargin (float, default=0.085): Height of the margin above the top pad,
        as a fraction of the height of the top pad.
    '''
    bottomMargin = kwargs.get('bottomMargin', 0.3)
    topMargin = kwargs.get('topMargin', 0.085)

    p.cd()

    heights = list(heights)

    # make room for bottom margin
    heights[-1] *= 1.+bottomMargin

    top = 1.
    bottom = sum(heights)
    pads = [_Pad(0.,bottom, 1., top)]
    pads[0].SetTopMargin(topMargin)
    pads[0].SetBottomMargin(0.005)
    for h in heights:
        top = bottom
        bottom = max(bottom - h, 0) # max to avoid "-0" errors
        pads.append(_Pad(0., bottom, 1., top))
        pads[-1].SetTopMargin(0.)
        pads[-1].SetBottomMargin(0)

    pads[-1].SetBottomMargin(bottomMargin)

    return tuple(pads)

def addPadBelow(p, height, bottomMargin=0.3, topMargin=0.085):
    '''
    Split pad p into two pads, and return (upper, lower). Just calls
    addPadsBelow(); mostly here for backwards compatibility.
    '''
    return addPadsBelow(p, height, bottomMargin=bottomMargin,
                        topMargin=topMargin)

def makeRatio(numerator, denominator):
    '''
    Return the graph of numerator/denominator and a line at y=1.
    '''
    if isinstance(numerator, _HistStack):
        for h in numerator.hists:
            h.sumw2()
        hNum = sum(numerator.hists)
    else:
        hNum = numerator.clone()
    try:
        hNum.sumw2()
    except AttributeError:
        pass
    num = _Graph(hNum)

    if isinstance(denominator, _HistStack):
        for h in denominator.hists:
            h.sumw2()
        hDenom = sum(denominator.hists)
    else:
        hDenom = denominator.clone()
    try:
        hDenom.sumw2()
    except AttributeError:
        pass
    denom = _Graph(hDenom)

    nRemoved = 0
    for i in range(num.GetN()):
        if hDenom[i+1].value <= 0. or hNum[i+1].value <= 0:
            num.RemovePoint(i - nRemoved)
            denom.RemovePoint(i - nRemoved)
            nRemoved += 1

    ratio = num / denom

    ratio.drawstyle = 'PE'
    ratio.color = 'black'

    unity = _Line(hNum.lowerbound(), 1, hNum.upperbound(), 1)
    unity.SetLineStyle(7)
    unity.SetLineWidth(2*unity.GetLineWidth())

    return ratio, unity

def fixRatioAxes(mainXAxis, mainYAxis, ratioXAxis, ratioYAxis,
                 mainPadHeight, ratioPadHeight):
    '''
    Remove the x axis title and labels from the main pad and recreate them by
    modifying the ratio plot axes.
    Resizes the y axis title and labels so they're the same size
    on both (the size they are on the main pad).
    '''
    ratioXAxis.title = mainXAxis.title
    ratioXAxis.SetTitleSize(mainXAxis.GetTitleSize() * mainPadHeight / ratioPadHeight)
    ratioXAxis.SetLabelSize(mainXAxis.GetLabelSize() * mainPadHeight / ratioPadHeight)

    ratioYAxis.SetTitleSize(mainYAxis.GetTitleSize() * mainPadHeight / ratioPadHeight)
    ratioYAxis.SetLabelSize(mainYAxis.GetLabelSize() * mainPadHeight / ratioPadHeight)
    ratioYAxis.SetTitleOffset(mainYAxis.GetTitleOffset() * ratioPadHeight / mainPadHeight)

    # mainXAxis.SetLabelOffset(20)
    # mainXAxis.SetTitleOffset(20)
    mainXAxis.SetTitle("")
    # don't print axis labels for the top plot
    for ib in range(mainXAxis.GetNbins()):
        mainXAxis.SetLabelSize(0)#Color(0,0.)

def makeErrorBand(hMean, errUp, errDn=None):
    '''
    Make a hatched black band to represent error bars around hMean.

    errUp and errDn are taken to be systematic errors, to be added in
    quadrature to the statistical error bars already on hMean. If they are
    floats, they are taken to be the fractional error on
    all bins. If they are histograms, the value of their bins is taken to be
    the absolute error on the corresponding bin of hMean.

    If errDn is not specified, errors are taken to be symmetric.
    '''
    if isinstance(hMean, _HistStack):
        hMean = sum(h for h in hMean)

    hUp = hMean.clone()
    hDn = hMean.clone()
    if isinstance(errUp, Number):
        if errDn is None:
            errDn = errUp

        for bMean, bUp, bDn in zip(hMean, hUp, hDn):
            bUp.value = bMean.value + sqrt(bMean.error**2 + (errUp * bMean.value)**2)
            bDn.value = max(bMean.value - sqrt(bMean.error**2 + (errDn * bMean.value)**2), 0.)
    else:
        if errDn is None:
            errDn = errUp.clone()

        for bMean, bUp, bDn, bErrUp, bErrDn in zip(hMean, hUp,
                                                   hDn, errUp, errDn):
            bUp.value = bMean.value + sqrt(bMean.error**2 + bErrUp.value**2)
            bDn.value = max(bMean.value - sqrt(bMean.error**2 + bErrDn.value**2), 0.)

    for bMean, bUp, bDn in zip(hMean, hUp, hDn):
        if (bUp.value < bMean.value or bDn.value > bMean.value) and not bMean.overflow:
            print "problem in bin {} ({:.2f} +{:.2f}/-{:.2f})".format(bMean.idx, bMean.value, bUp.value, bDn.value)

    err = _band(hDn, hUp, hMean)
    err.SetLineWidth(0)
    err.SetFillColorAlpha(1,0.3)
    err.fillstyle = 'solid'
    err.drawstyle = '2'
    err.title = r'\textbf{Stat.} \ \! \boldsymbol{\oplus} \ \! \textbf{syst.\ unc.}'
    err.legendstyle = 'F'
    try:
        err.SetLineColorAlpha(_Color(err.GetLineColor())('root'),0.)
    except TypeError:
        err.SetLineColorAlpha(err.GetLineColor(),0.)

    return err
