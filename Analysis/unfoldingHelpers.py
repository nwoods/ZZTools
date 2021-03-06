import logging
from rootpy import log as rlog; rlog = rlog["/unfoldingHelpers"]
# don't show most silly ROOT messages
logging.basicConfig(level=logging.WARNING)
rlog["/ROOT.TUnixSystem.SetDisplay"].setLevel(rlog.ERROR)

from SampleTools import SampleGroup as _Group
from SampleTools import SampleStack as _Stack
from Utilities import mapObjects as _mapObjects, identityFunction as _identityFunction, \
    zMassDist as _zMassDist, Z_MASS as _MZ

from rootpy import ROOTError as _RootError
from rootpy.ROOT import RooUnfoldResponse as _Response
from rootpy.plotting import Hist2D as _Hist2D, Graph as _Graph
import rootpy.compiled as _RootComp

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


def _objVar(row, var, obj):
    return getattr(row, obj+var)
def _nObjVar(row, var, *objects):
    toGet = '_'.join(list(objects)+[var])
    return getattr(row, toGet)

_RootComp.register_code('''
#include "TLorentzVector.h"

float _getInvariantMass(float pt1, float eta1, float phi1,
                        float pt2, float eta2, float phi2)
{
  TLorentzVector p1, p2;
  p1.SetPtEtaPhiM(pt1, eta1, phi1, 0.);
  p2.SetPtEtaPhiM(pt2, eta2, phi2, 0.);

  return (p1+p2).M();
}
''', ['_getInvariantMass'])
_asdf = _RootComp._getInvariantMass(0,0,0,0,0,0) # force to compile
_fInvMass = _RootComp._getInvariantMass
# note no leading '_'
fInvMassString = '_getInvariantMass({0}Pt, {0}Eta, {0}Phi, {1}Pt, {1}Eta, {1}Phi)'

def _getInvMass(row, lep1, lep2):
    return _fInvMass(_objVar(row, 'Pt', lep1), _objVar(row, 'Eta', lep1),
                     _objVar(row, 'Phi', lep1), _objVar(row, 'Pt', lep2),
                     _objVar(row, 'Eta', lep2), _objVar(row, 'Phi', lep2))

def _makeWrongZRejecter(channel):
    # only needed for 4e and 4mu
    if not all(lep == channel[0] for lep in channel):
        return _identityFunction

    obj = _mapObjects(channel)
    def f(row):
        dZ1 = _zMassDist(_nObjVar(row, 'Mass', obj[0], obj[1]))
        return ((_objVar(row, 'Charge', obj[0]) == _objVar(row, 'Charge', obj[2]) or _zMassDist(_getInvMass(row, obj[0], obj[2])) > dZ1) and
                (_objVar(row, 'Charge', obj[0]) == _objVar(row, 'Charge', obj[3]) or _zMassDist(_getInvMass(row, obj[0], obj[3])) > dZ1) and
                (_objVar(row, 'Charge', obj[1]) == _objVar(row, 'Charge', obj[2]) or _zMassDist(_getInvMass(row, obj[1], obj[2])) > dZ1) and
                (_objVar(row, 'Charge', obj[1]) == _objVar(row, 'Charge', obj[3]) or _zMassDist(_getInvMass(row, obj[1], obj[3])) > dZ1))
    return f

_wrongZRejectionTemp = ('(({0}1Charge=={0}3Charge||abs(_getInvariantMass({0}1Pt,{0}1Eta,{0}1Phi,{0}3Pt,{0}3Eta,{0}3Phi)-{1}) > abs({0}1_{0}2_Mass-{1})) && '
                        ' ({0}1Charge=={0}4Charge||abs(_getInvariantMass({0}1Pt,{0}1Eta,{0}1Phi,{0}4Pt,{0}4Eta,{0}4Phi)-{1}) > abs({0}1_{0}2_Mass-{1})) && '
                        ' ({0}2Charge=={0}3Charge||abs(_getInvariantMass({0}2Pt,{0}2Eta,{0}2Phi,{0}3Pt,{0}3Eta,{0}3Phi)-{1}) > abs({0}1_{0}2_Mass-{1})) && '
                        ' ({0}2Charge=={0}4Charge||abs(_getInvariantMass({0}2Pt,{0}2Eta,{0}2Phi,{0}4Pt,{0}4Eta,{0}4Phi)-{1}) > abs({0}1_{0}2_Mass-{1})))')
_wrongZRejectionStr = {
    'eeee' : _wrongZRejectionTemp.format('e', _MZ),
    'mmmm' : _wrongZRejectionTemp.format('m', _MZ),
    'eemm' : '1'
    }


_genCache = {}
_cacheVar = ''
_cacheChannel = ''
def _getGenVarDict(channel, truth, var, varFunction, selectionFunction):
    '''
    Get a dictionary of gen info for this channel, variable, and selection,
    with caching for performance.

    channel (str): single channel to use
    truth (SampleGroup): gen-level sample
    var (str or iterable of str): variable in question (used only for caching)
    varFunction (callable or iterable of callable): function to extract the
        variable from an ntuple row
    selectionFunction (callable): function that takes an ntuple row and returns
        a boolean to select events we care about

    Returns a variable- and channel-specific dict whose keys are sample names.
        Each corresponding value is a dict containing the variable values keyed
        to the event ID in the form of a tuple containing (run, lumi, evt)
    '''
    # cache gen vars because they're slow to collect and don't change much
    global _genCache
    global _cacheVar
    global _cacheChannel

    # variables given as lists need a name for the cache
    if not isinstance(var, str) and hasattr(var, '__iter__'):
        varName = ''.join(var)
    else:
        varName = var

    if _cacheVar != varName or _cacheChannel != channel:
        _genCache.clear()
        _cacheVar = varName
        _cacheChannel = channel

    for name, sample in truth.itersamples():
        if name in _genCache:
            print "Reusing {} {}".format(channel, varName)
            continue
        fBestZ = _makeWrongZRejecter(channel)
        if hasattr(varFunction, '__iter__'):
            _genCache[name] = {
                (row.run, row.lumi, row.evt) : [v(row) for v in varFunction]
                for row in sample if selectionFunction(row) and fBestZ(row)
                }
        else:
            _genCache[name] = {
                (row.run, row.lumi, row.evt) : varFunction(row)
                for row in sample if selectionFunction(row) and fBestZ(row)
                }

    return _genCache


def getResponse(channel, truth, mc, bkg, var, varFunction, binning, fPUWeight,
                lepSyst='', altVar='', selectionStr='',
                selectionFunction=_identityFunction,
                selectionStrAlt='', varFunctionAlt=None,
                selectionFunctionAlt=None):
    '''
    Get the unfolding response matrix as a RooUnfoldResponse object.
    channel (str): single channel to use
    truth (SampleGroup): gen-level information
    stack (SampleStack): MC+background stack
    var (str or iterable of str): distribution to unfold
    varFunction (callable or iterable of callable): function(s) that take(s)
        an ntuple row and return(s) the same quantity(ies) as var. If an
        iterable, all are used with the same selection; if only one of several
        values should be used (as with the variable strings), the selection
        should be done inside the function.
    binning (list): binning for var
    fPUWeight (callable): function for calculating pileup weight from nTruePU
    lepSyst (str): if 'up' or 'dn'/'down', will shift the lepton efficiency
        scale factors up/down by 1 sigma
    altVar (str or iterable of str): If non-empty, background and truth
        ntuples use this instead of var (for systematic shifts and other
        variables that only make sense for MC reco)
    selectionStr (str or iterable of str): selection to apply when using draw
        strings
    selectionFunction (callable): function that takes an ntuple row and
        returns a boolean to select events we care about
    selectionStrAlt (str or iterable or str): Selection string for background
        and truth
    varFunctionAlt (callable or None): function equivalent to selectionFunction
        but for true-level ntuples only
    selectionFunctionAlt (callable, list of callable, or None): function(s)
        corresponding to selectionStrAlt
    '''
    if not altVar:
        altVar = var
    if not selectionStrAlt:
        selectionStrAlt = selectionStr

    if not selectionStrAlt:
        selectionStrGen = _wrongZRejectionStr[channel]
    elif isinstance(selectionStrAlt, str):
        selectionStrGen = selectionStrAlt + ' && ' + _wrongZRejectionStr[channel]
    else:
        # better be an iterable of strings
        selectionStrGen = [(s + ' && ' if s else '') + _wrongZRejectionStr[channel] for s in selectionStrAlt]

    if varFunctionAlt is None:
        varFunctionAlt = varFunction
    if selectionFunctionAlt is None:
        selectionFunctionAlt = selectionFunction

    hReco = mc.makeHist(var, selectionStr, binning, perUnitWidth=False)
    hReco += bkg.makeHist(altVar, selectionStrAlt, binning, perUnitWidth=False)

    if len(binning) == 3:
        hResponse = _Hist2D(*(binning+binning))
    else:
        hResponse = _Hist2D(binning, binning)

    mcList = _baseSamples(mc)

    fLepSF = _makeScaleFactorFunction(channel, lepSyst)

    genVar = _getGenVarDict(channel, truth, altVar, varFunctionAlt,
                            selectionFunctionAlt)

    if hasattr(varFunction, '__iter__'):
        def _fill(h, rw, gen, wt):
            for f,g in zip(varFunction, gen):
                h.Fill(f(rw), g, wt)
    else:
        def _fill(h, rw, gen, wt):
            h.Fill(varFunction(rw), gen, wt)

    for sample in mcList:
        name = sample.name
        wConst = sample.xsec*sample.intLumi/sample.sumW
        for row in sample:
            if selectionFunction(row):
                evtID = (row.run, row.lumi, row.evt)
                try:
                    genInfo = genVar[name][evtID]
                except KeyError:
                    # it's ok to miss an event, not a whole sample
                    if name not in genVar:
                        raise
                    continue
                weight = fPUWeight(row.nTruePU) * row.genWeight * fLepSF(row) * wConst
                _fill(hResponse, row, genInfo, weight)

    hTrue = truth.makeHist(altVar,
                           selectionStrGen,
                           binning, perUnitWidth=False)

    return _Response(hReco, hTrue, hResponse)


class _NoPDFWeightLogging(logging.Filter):
    '''
    Suppress ROOT messages about 'pdfWeights' not appearing in trees, since
    some samples don't have them anyway.
    It's ok to remove these messages for other samples because they will throw
    an exception anyway.
    '''
    def filter(self, record):
        return 'Bad numerical expression : "pdfWeights"' not in record.msg

rlog["/ROOT.TTreeFormula.Compile"].addFilter(_NoPDFWeightLogging())

def getResponsePDFErrors(channel, truth, mc, bkg, var, varFunction,
                         binning, fPUWeight, lepSyst='', altVar='',
                         selectionStr='', selectionFunction=_identityFunction,
                         selectionStrAlt='', varFunctionAlt=None,
                         selectionFunctionAlt=None):
    '''
    As getResponse() above, but with a systematic shift corresponding to a
    1sigma shift up and down in the PDF.

    Inputs and outputs as in getResponse() except that the return value is
    a 2-tuple (responseDown, responseUp)
    '''
    if not altVar:
        altVar = var
    if not selectionStrAlt:
        selectionStrAlt = selectionStr

    if not selectionStrAlt:
        selectionStrGen = _wrongZRejectionStr[channel]
    elif isinstance(selectionStrAlt, str):
        selectionStrGen = selectionStrAlt + ' && ' + _wrongZRejectionStr[channel]
    else:
        # better be an iterable of strings
        selectionStrGen = [(s + ' && ' if s else '') + _wrongZRejectionStr[channel] for s in selectionStrAlt]

    if varFunctionAlt is None:
        varFunctionAlt = varFunction
    if selectionFunctionAlt is None:
        selectionFunctionAlt = selectionFunction

    hRecoUp = mc.makeHist(var, selectionStr, binning, perUnitWidth=False)
    hRecoUp += bkg.makeHist(altVar, selectionStrAlt, binning, perUnitWidth=False)
    hRecoDn = hRecoUp.clone()

    # get 2D hist of var vs PDF variation index for each sample
    hRecoVariations = []
    for s in mc.values():
        try:
            hRecoVariations.append(s.makeHist2(var, 'Iteration$', selectionStr, binning,
                                   [100, 0, 100], 'pdfWeights', False))
        except _RootError:
            # MCFM samples don't have LHE info
            if 'GluGluZZ' not in s.name:
                raise

    # for each var bin in each sample, get the RMS across all the variations
    allRMSes = [[_Graph(h.ProjectionY('slice{}'.format(i), i+1,i+1)).GetRMS(2) for i in xrange(h.GetNbinsX())] for h in hRecoVariations]
    # for each var bin, add variations for all samples
    binRMSes = [sum(rmses) for rmses in zip(*allRMSes)]

    # apply variations
    for i in xrange(hRecoUp.GetNbinsX()):
        hRecoUp[i+1].value += binRMSes[i]
        hRecoDn[i+1].value = max(0.,hRecoDn[i+1].value + binRMSes[i])

    # same for gen info
    hTrueUp = truth.makeHist(altVar, selectionStrGen, binning,
                             perUnitWidth=False)
    hTrueDn = hTrueUp.clone()

    hTrueVariations = []
    for s in truth.values():
        try:
            hTrueVariations.append(s.makeHist2(var, 'Iteration$',
                                               selectionStrGen, binning,
                                               [100, 0, 100], 'pdfWeights',
                                               False))
        except _RootError:
            # MCFM samples don't have LHE info
            if 'GluGluZZ' not in s.name:
                raise

    allTrueRMSes = [[_Graph(h.ProjectionY('slice{}'.format(i), i+1,i+1)).GetRMS(2) for i in xrange(h.GetNbinsX())] for h in hTrueVariations]
    binTrueRMSes = [sum(rmses) for rmses in zip(*allTrueRMSes)]
    # get RMS as fraction of bin size for later use
    trueRMSFrac = [rms / hTrueUp[i+1].value if hTrueUp[i+1].value else 0. for i,rms in enumerate(binTrueRMSes)]

    for i in xrange(hTrueUp.GetNbinsX()):
        hTrueUp[i+1].value += binTrueRMSes[i]
        hTrueDn[i+1].value = max(0.,hTrueDn[i+1].value + binTrueRMSes[i])


    if len(binning) == 3:
        hResponseUp = _Hist2D(*(binning+binning))
    else:
        hResponseUp = _Hist2D(binning, binning)

    mcList = _baseSamples(mc)

    fLepSF = _makeScaleFactorFunction(channel, lepSyst)

    genVar = _getGenVarDict(channel, truth, altVar, varFunctionAlt,
                            selectionFunctionAlt)

    if hasattr(varFunction, '__iter__'):
        def _fill(h, rw, gen, wt):
            for f,g in zip(varFunction, gen):
                h.Fill(f(rw), g, wt)
    else:
        def _fill(h, rw, gen, wt):
            h.Fill(varFunction(rw), gen, wt)

    for sample in mcList:
        name = sample.name
        wConst = sample.xsec*sample.intLumi/sample.sumW
        for row in sample:
            if selectionFunction(row):
                evtID = (row.run, row.lumi, row.evt)
                try:
                    genInfo = genVar[name][evtID]
                except KeyError:
                    if name not in genVar:
                        raise
                    continue
                weight = fPUWeight(row.nTruePU) * row.genWeight * fLepSF(row) * wConst
                _fill(hResponseUp, row, genInfo, weight)

    hResponseDn = hResponseUp.clone()
    # assume matrix is diagonal enough that we can just use the variance of the
    # gen bin
    for bUp,bDn in zip(hResponseUp, hResponseDn):
        if bUp.overflow:
            continue
        bUp.value *= (1.+trueRMSFrac[bUp.xyz[1]-1])
        bDn.value *= (1.-trueRMSFrac[bDn.xyz[1]-1])

    return _Response(hRecoDn, hTrueDn, hResponseDn), _Response(hRecoUp, hTrueUp, hResponseUp)


_variationIndices = [1,2,3,4,6,8] # indices of the scale variation sets we care about
def getResponseScaleErrors(channel, truth, mc, bkg, var, varFunction, binning, fPUWeight,
                           lepSyst='', altVar='', selectionStr='',
                           selectionFunction=_identityFunction,
                           selectionStrAlt='', varFunctionAlt=None,
                           selectionFunctionAlt=None):
    '''
    Get the unfolding response matrices for all interesting scale variations.
    channel (str): single channel to use
    truth (SampleGroup): gen-level information
    stack (SampleStack): MC+background stack
    var (str or iterable of str): distribution to unfold
    varFunction (callable): function that takes an ntuple row and returns the
        same quantity as var
    binning (list): binning for var
    fPUWeight (callable): function for calculating pileup weight from nTruePU
    lepSyst (str): if 'up' or 'dn'/'down', will shift the lepton efficiency
        scale factors up/down by 1 sigma
    altVar (str or iterable of str): If non-empty, background and truth
        ntuples use this instead of var (for systematic shifts and other
        variables that only make sense for MC reco)
    selectionStr (str or iterable of str): selection to apply when using draw
        strings
    selectionFunction (callable): function that takes an ntuple row
        and returns a boolean to select the same events as selectionStr
    selectionStrAlt (str or iterable or str): Selection string for background
        and truth
    varFunctionAlt (callable or None): function corresponding to altVar
    selectionFunctionAlt (callable or None): function corresponding to
        selectionStrAlt
    '''
    if not altVar:
        altVar = var
    if not selectionStrAlt:
        selectionStrAlt = selectionStr

    if not selectionStrAlt:
        selectionStrGen = _wrongZRejectionStr[channel]
    elif isinstance(selectionStrAlt, str):
        selectionStrGen = selectionStrAlt + ' && ' + _wrongZRejectionStr[channel]
    else:
        # better be an iterable of strings
        selectionStrGen = [(s + ' && ' if s else '') + _wrongZRejectionStr[channel] for s in selectionStrAlt]

    if varFunctionAlt is None:
        varFunctionAlt = varFunction
    if selectionFunctionAlt is None:
        selectionFunctionAlt = selectionFunction

    hReco = [mc.makeHist(var, selectionStr, binning,
                         {
                'ZZTo4L':'scaleWeights[{}]'.format(i),
                'ZZTo4L-amcatnlo':'scaleWeights[{}]'.format(i),
                'ZZJJTo4L_EWK':'scaleWeights[{}]'.format(i),
                },
                         perUnitWidth=False)
             for i in _variationIndices]
    for h in hReco:
        h += bkg.makeHist(altVar, selectionStrAlt, binning, perUnitWidth=False)

    mcList = _baseSamples(mc)

    fLepSF = _makeScaleFactorFunction(channel, lepSyst)

    genVar = _getGenVarDict(channel, truth, altVar, varFunctionAlt,
                            selectionFunctionAlt)

    if hasattr(varFunction, '__iter__'):
        def _fill(h, rw, gen, wt):
            for f,g in zip(varFunction, gen):
                h.Fill(f(rw), g, wt)
    else:
        def _fill(h, rw, gen, wt):
            h.Fill(varFunction(rw), gen, wt)

    hResponse = []
    for iVar in _variationIndices:
        if len(binning) == 3:
            hResponse.append(_Hist2D(*(binning+binning)))
        else:
            hResponse.append(_Hist2D(binning, binning))

        for sample in mcList:
            name = sample.name
            wConst = sample.xsec*sample.intLumi/sample.sumW
            if 'GluGluZZ' in name:
                fScaleWt = lambda *args: 1.
            else:
                fScaleWt = lambda row: row.scaleWeights.at(iVar)

            for row in sample:
                if selectionFunction(row):
                    evtID = (row.run, row.lumi, row.evt)
                    try:
                        genInfo = genVar[name][evtID]
                    except KeyError:
                        if name not in genVar:
                            raise
                        continue
                    weight = fPUWeight(row.nTruePU) * row.genWeight * fLepSF(row) * wConst * fScaleWt(row)
                    _fill(hResponse[-1], row, genInfo, weight)

    hTrue = [truth.makeHist(altVar, selectionStrGen, binning,
                            {
                'ZZTo4L':'scaleWeights[{}]'.format(i),
                'ZZTo4L-amcatnlo':'scaleWeights[{}]'.format(i),
                'ZZJJTo4L_EWK':'scaleWeights[{}]'.format(i),
                },
                            perUnitWidth=False)
             for i in _variationIndices]

    return [_Response(hR, hT, hResp) for hR, hT, hResp in zip(hReco, hTrue, hResponse)]


_alphaSIndices = [100,101]
def getResponseAlphaSErrors(channel, truth, mc, bkg, var, varFunction, binning, fPUWeight,
                            lepSyst='', altVar='', selectionStr='',
                            selectionFunction=_identityFunction,
                            selectionStrAlt='', varFunctionAlt=None,
                            selectionFunctionAlt=None):
    '''
    Get the unfolding response matrices for alpha_s varied up and down.
    truth (SampleGroup): gen-level information
    stack (SampleStack): MC+background stack
    var (str or iterable of str): distribution to unfold
    varFunction (callable): function that takes an ntuple row and returns the
        same quantity as var
    binning (list): binning for var
    fPUWeight (callable): function for calculating pileup weight from nTruePU
    lepSyst (str): if 'up' or 'dn'/'down', will shift the lepton efficiency
        scale factors up/down by 1 sigma
    altVar (str or iterable of str): If non-empty, background and truth
        ntuples use this instead of var (for systematic shifts and other
        variables that only make sense for MC reco)
    selectionStr (str or iterable of str): selection to apply when using draw
        strings
    selectionFunction (callable): function that takes an ntuple row
        and returns a boolean to select the same events as selectionStr
    selectionStrAlt (str or iterable or str): Selection string for background
        and truth
    varFunctionAlt (callable or None): function corresponding to altVar
    selectionFunctionAlt (callable or None): function corresponding to
        selectionStrAlt
    '''
    if not altVar:
        altVar = var
    if not selectionStrAlt:
        selectionStrAlt = selectionStr

    if not selectionStrAlt:
        selectionStrGen = _wrongZRejectionStr[channel]
    elif isinstance(selectionStrAlt, str):
        selectionStrGen = selectionStrAlt + ' && ' + _wrongZRejectionStr[channel]
    else:
        # better be an iterable of strings
        selectionStrGen = [(s + ' && ' if s else '') + _wrongZRejectionStr[channel] for s in selectionStrAlt]

    if varFunctionAlt is None:
        varFunctionAlt = varFunction
    if selectionFunctionAlt is None:
        selectionFunctionAlt = selectionFunction

    hReco = [mc.makeHist(var, selectionStr, binning,
                         {
                'ZZTo4L':'pdfWeights[{}]'.format(i),
                'ZZTo4L-amcatnlo':'pdfWeights[{}]'.format(i),
                'ZZJJTo4L_EWK':'pdfWeights[{}]'.format(i),
                },
                         perUnitWidth=False)
             for i in _alphaSIndices]
    for h in hReco:
        h += bkg.makeHist(altVar, selectionStrAlt, binning, perUnitWidth=False)

    mcList = _baseSamples(mc)

    fLepSF = _makeScaleFactorFunction(channel, lepSyst)

    genVar = _getGenVarDict(channel, truth, altVar, varFunctionAlt,
                            selectionFunctionAlt)

    if hasattr(varFunction, '__iter__'):
        def _fill(h, rw, gen, wt):
            for f,g in zip(varFunction, gen):
                h.Fill(f(rw), g, wt)
    else:
        def _fill(h, rw, gen, wt):
            h.Fill(varFunction(rw), gen, wt)

    hResponse = []
    for iVar in _alphaSIndices:
        if len(binning) == 3:
            hResponse.append(_Hist2D(*(binning+binning)))
        else:
            hResponse.append(_Hist2D(binning, binning))

        for sample in mcList:
            name = sample.name
            wConst = sample.xsec*sample.intLumi/sample.sumW
            if 'GluGluZZ' in name:
                fScaleWt = lambda *args: 1.
            else:
                fScaleWt = lambda row: row.pdfWeights.at(iVar)

            for row in sample:
                if selectionFunction(row):
                    evtID = (row.run, row.lumi, row.evt)
                    try:
                        genInfo = genVar[name][evtID]
                    except KeyError:
                        if name not in genVar:
                            raise
                        continue
                    weight = fPUWeight(row.nTruePU) * row.genWeight * fLepSF(row) * wConst * fScaleWt(row)
                    _fill(hResponse[-1], row, genInfo, weight)

    hTrue = [truth.makeHist(altVar, selectionStrGen, binning,
                            {
                'ZZTo4L':'pdfWeights[{}]'.format(i),
                'ZZTo4L-amcatnlo':'pdfWeights[{}]'.format(i),
                'ZZJJTo4L_EWK':'pdfWeights[{}]'.format(i),
                },
                            perUnitWidth=False)
             for i in _alphaSIndices]

    return tuple(_Response(hR, hT, hResp) for hR, hT, hResp in zip(hReco, hTrue, hResponse))


