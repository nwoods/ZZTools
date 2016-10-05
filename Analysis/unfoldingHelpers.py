from SampleTools import SampleGroup as _Group
from SampleTools import SampleStack as _Stack
from Utilities import mapObjects as _mapObjects

from rootpy.ROOT import RooUnfoldResponse as _Response
from rootpy.plotting import Hist2D as _Hist2D

from operator import mul as _times

def _normalizeBins(h):
    binUnit = min(h.GetBinWidth(b) for b in range(1,len(h)+1))
    for ib in xrange(1,len(h)+1):
        w = h.GetBinWidth(ib)
        h.SetBinContent(ib, h.GetBinContent(ib) * binUnit / w)
        h.SetBinError(ib, h.GetBinError(ib) * sqrt(binUnit / w))
        if h.GetBinError(ib) > h.GetBinContent(ib):
            h.SetBinError(ib, h.GetBinContent(ib))
    h.sumw2()


def _baseSamples(sample):
    '''
    From a sample group or stack, gets a list of the the lowest-level samples
    in the composite.
    '''
    out = []
    if isinstance(sample, _Group):
        for s in sample.values():
            out += _baseSamples(s)
    elif isinstance(sample, _Stack):
        for s in sample:
            out += _baseSamples(s)
    else:
        out.append(sample)

    return out

def _makeScaleFactorFunction(channel, syst=''):
    objects = _mapObjects(channel)
    
    if syst.lower() == 'up':
        return lambda row: reduce(_times,(getattr(row, ob+'EffScaleFactor')+getattr(row, ob+'EffScaleFactorError') for ob in objects))
    if syst.lower() in ['dn','down']:
        return lambda row: reduce(_times,(getattr(row, ob+'EffScaleFactor')-getattr(row, ob+'EffScaleFactorError') for ob in objects))
    return lambda row: reduce(_times,(getattr(row, ob+'EffScaleFactor') for ob in objects))

_genVars = {}

def getResponse(channel, truth, mc, bkg, var, binning, fPUWeight, lepSyst=''):
    '''
    Get the unfolding response matrix as a RooUnfoldResponse object.
    channel (str): single channel to use (matters for lepton scale factors)
    truth (SampleGroup): gen-level information
    stack (SampleStack): MC+background stack
    var (str): distribution to unfold
    binning (list): binning for var
    fPUWeight (callable): function for calculating pileup weight from nTruePU
    lepSyst (str): if 'up' or 'dn'/'down', will shift the lepton efficiency
        scale factors up/down by 1 sigma
    '''
    # cache gen vars because they're slow to collect and don't change much
    global _genVars

    try:
        varDict = _genVars[var]
        print "Already have something for variable {}".format(var)
    except KeyError:
        varDict = {}
        _genVars[var] = varDict
    try:
        channelDict = varDict[channel]
        print "Already have something there for channel {}".format(channel)
    except KeyError:
        channelDict = {}
        varDict[channel] = channelDict

    genVar = {}
    for name, sample in truth.itersamples():
        try:
            genVar[name] = channelDict[name]
            print 'Reusing gen events for gen {}, in {} sample {}!'.format(var,channel,name)
        except KeyError:
            sampleGenVar = {}

            for row in sample:
                sampleGenVar[(row.run,row.lumi,row.evt)] = getattr(row, var)

            genVar[name] = sampleGenVar
            channelDict[name] = genVar[name]

    hReco = sum(mc.makeHist(var, '', binning, perUnitWidth=False).hists)
    hReco += bkg.makeHist(var, '', binning, perUnitWidth=False)

    if len(binning) == 3:
        hResponse = _Hist2D(*(binning+binning))
    else:
        hResponse = _Hist2D(binning, binning)

    mcList = _baseSamples(mc)

    fLepSF = _makeScaleFactorFunction(channel, lepSyst)

    for sample in mcList:
        name = sample.name
        wConst = sample.xsec*sample.intLumi/sample.sumW
        for row in sample:
            evtID = (row.run, row.lumi, row.evt)
            weight = fPUWeight(row.nTruePU) * row.genWeight * fLepSF(row) * wConst
            try:
                hResponse.Fill(getattr(row, var), genVar[name][evtID], weight)
            except KeyError:
                pass

    hTrue = truth.makeHist(var, '', binning, perUnitWidth=False)

    return _Response(hReco, hTrue, hResponse)

