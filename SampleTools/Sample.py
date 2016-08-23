#############################################################################
#                                                                           #
#    Classes for objects that contain, manipulate, and use physics          #
#    samples.                                                               #
#                                                                           #
#    Author: Nate Woods, U. Wisconsin                                       #
#                                                                           #
#############################################################################


import logging
from rootpy import log as rlog; rlog = rlog["/Sample"]
# don't show most silly ROOT messages
logging.basicConfig(level=logging.WARNING)
rlog["/ROOT.TCanvas.Print"].setLevel(rlog.WARNING)
rlog["/ROOT.TUnixSystem.SetDisplay"].setLevel(rlog.ERROR)
rlog["/rootpy.tree.chain"].setLevel(rlog.WARNING)

from Metadata.metadata import sampleInfo as _samples
from Utilities import combineWeights

from rootpy.io import root_open, TemporaryFile
from rootpy.io import DoesNotExist as _RootpyDNE
from rootpy.tree import Tree, TreeChain
from rootpy.plotting import Hist, Canvas

from glob import glob
from math import sqrt
from os import remove as _rm
from os import close as _close

# Workaround for weird ROOT bug
_dummy = Hist(1,0,1)
_cDummy = Canvas(10,10)
_dummy.xaxis.SetTitle('\\mu\\mu')
_dummy.draw()


class _SampleBase(object):
    def __init__(self, name, channel, dataIn, initFromMetadata=False, *args, **kwargs):
        self.name = name
        self.channel = channel
        self.shortName = name
        self.prettyName = name

        self.formatDefault()

        if initFromMetadata:
            self.initFromMetadata()

        self.storeInputs(dataIn)


    def initFromMetadata(self):
        if self.name not in _samples:
            raise ValueError("Sample {} is not in the metadata".format(self.name))

        info = _samples[self.name]

        self.shortName = info.get('shortName', self.name)
        self.prettyName = info.get('prettyName', self.shortName)

        self.format(False, **info.get('format', {}))


    def storeInputs(self, inputs):
        pass


    def format(self, reset=False, **newFormat):
        if reset:
            self._format = newFormat.copy()
        else:
            self._format.update(newFormat)


    def formatDefault(self):
        self._format = {}


    @property
    def channel(self):
        return self._channel
    @channel.setter
    def channel(self, ch):
        self._channel = ch

    @property
    def name(self):
        return self._name
    @name.setter
    def name(self, name):
        self._name = name

    @property
    def shortName(self):
        return self._shortName
    @shortName.setter
    def shortName(self, name):
        self._shortName = name

    @property
    def prettyName(self):
        return self._prettyName
    @prettyName.setter
    def prettyName(self, name):
        self._prettyName = name



class NtupleSample(_SampleBase):
    def __init__(self, name, channel, dataIn, initFromMetadata=False, *args, **kwargs):
        self.weight = ''

        super(NtupleSample, self).__init__(name, channel, dataIn, initFromMetadata, *args, **kwargs)

            
    def storeInputs(self, inputs):
        if isinstance(inputs, str):
            inputs = [inputs]

        self.files = []
        for inp in inputs:
            self.files += glob(inp)

        if not self.files:
            raise IOError("No files found for {}".format(self.name))

        self.ntuple = self.combineNtuples(self.files, self.channel)

    def combineNtuples(self, files, chan):
        '''
        Makes a TreeChain of the files (or a bare TTree if there's just one 
        file). No redundant row culling or anything.
        '''
        if len(files) == 1:
            self.ntupleFile = root_open(files[0])
            return self.ntupleFile.Get('{}/ntuple'.format(chan))

        return TreeChain('{}/ntuple'.format(chan), files)


    def addToHist(self, hist, var, selection):
        self.ntuple.Draw(var, selection, 'goff', hist)
        hist.sumw2()

    def makeHist(self, var, selection, binning, weight='', perUnitWidth=True, **kwargs):
        '''
        If var, selection, and/or weight are iterables (which must be of the 
        same length), the hists from the resulting var/selection/weight sets
        will be added together. If one or two are strings, they are reused the
        appropriate number of times.
        '''
        if len(binning) != 3:
            binning = [binning]
        
        # use TH1D instead of TH1F because some datasets are now big enough for
        # floating point stuff to matter (!!)
        h = Hist(*binning, type='D', title=self.prettyName, **self._format)

        nToAdd = max(1 if isinstance(var,str) else len(var),
                     1 if isinstance(selection,str) else len(selection),
                     1 if isinstance(weight,str) else len(weight))
        if isinstance(var, str):
            var = nToAdd * [var]
        if isinstance(selection, str):
            selection = nToAdd * [selection]
        if isinstance(weight, str):
            weight = nToAdd * [weight]

        assert len(var) == len(selection) and len(selection) == len(weight), \
            "Invalid plotting parameters! Variable: {}, selection: {}, weight: {}.".format(var,selection,weight)

        for v, s, w in zip(var, selection, weight):
            w = combineWeights(w, self.weight)

            s = combineWeights(w, s)

            self.addToHist(h, v, s)

        h.sumw2()

        if perUnitWidth and not h.uniform():
            binUnit = min(h.GetBinWidth(b) for b in range(1,len(h)+1))
            for ib in xrange(1,len(h)+1):
                w = h.GetBinWidth(ib)
                h.SetBinContent(ib, h.GetBinContent(ib) * binUnit / w)
                h.SetBinError(ib, h.GetBinContent(ib) * sqrt(binUnit / w))
            h.sumw2()

        return h


    def implicitWeight(self):
        return ''
        

    def applyWeight(self, w, reset=False):
        '''
        Factor by which to always weight the sample. If reset evaluates to
        True, previous weights are removed (does not effect default weights
        like MC cross sections and so on if these are enabled).
        '''
        if reset:
            self.weight = str(w)
        else:
            self.weight = combineWeights(self.weight, w)

    @property
    def weight(self):
        return combineWeights(self._weight, self.implicitWeight())
    @weight.setter
    def weight(self, w):
        if not isinstance(w, str):
            w = str(w)
        self._weight = w



class MCSample(NtupleSample):
    def __init__(self, name, channel, dataIn, initFromMetadata=False, intLumi=1000, *args, **kwargs):
        self.isSignal = False
        self.sumW = -1
        self.xsec = -1
        self.kFactor = '1'
        self.intLumi = intLumi
        super(MCSample, self).__init__(name, channel, dataIn, initFromMetadata, *args, **kwargs)


    def initFromMetadata(self):
        super(MCSample, self).initFromMetadata()

        info = _samples[self.name]

        self.isSignal = info.get('isSignal', False)
        self.sumW = info.get('sumW', 0)
        self.xsec = info.get('xsec', -1)
        self.kFactor = info.get('kFactor', '1.')


    def storeInputs(self, inputs):
        super(MCSample, self).storeInputs(inputs)

        # get the sum of the weights if applicable
        if len(self.files) == 1:
            try:
                with root_open(self.files[0]) as f:
                    metaTree = f.Get('metaInfo/metaInfo')
                    self.sumW = metaTree.Draw('1', 'summedWeights').Integral()
            except _RootpyDNE:
                pass
        else:
            try:
                metaChain = TreeChain('metaInfo/metaInfo', self.files)
                self.sumW = metaChain.Draw('1', 'summedWeights').Integral()
            except _RootpyDNE:
                pass


    def implicitWeight(self):
        if self.xsec > 0:
            w = 'genWeight * {} * {} * {} / {}'.format(self.intLumi, self.xsec, 
                                                       self.kFactor, self.sumW)
        else:
            w = ''

        return combineWeights(super(MCSample, self).implicitWeight(), w)
            

    def formatDefault(self):
        self.format(True, drawstyle='hist', fillstyle='solid', legendstyle='F')


    @property
    def isSignal(self):
        return self._signal
    @isSignal.setter
    def isSignal(self, sig):
        self._signal = sig

    @property
    def sumW(self):
        return self._sumW
    @sumW.setter
    def sumW(self, val):
        self._sumW = val

    @property
    def xsec(self):
        return self._xsec
    @xsec.setter
    def xsec(self, val):
        self._xsec = val

    @property
    def kFactor(self):
        return self._kFactor
    @kFactor.setter
    def kFactor(self, k):
        self._kFactor = k

    @property
    def intLumi(self):
        return self._intLumi
    @intLumi.setter
    def intLumi(self, il):
        self._intLumi = il



class _FileWrapper(object):
    '''
    Hold a temporary file and delete it when this object is deleted.
    *** FIXME *** __del__ is fragile
    '''
    def __init__(self, f):
        self.f = f
        
        # *** FIXME *** there must be a better way to do this
        self.descriptor = getattr(self.f, '_{}__fd'.format(self.f.__class__.__name__))
        self.name = getattr(self.f, '_{}__tmp_path'.format(self.f.__class__.__name__))
    def __del__(self):
        if self.f.IsOpen():
            self.f.Close()
        try:
            _close(self.descriptor)
        except OSError:
            pass
        try:
            _rm(self.name)
        except OSError:
            rlog.warning("Failed to delete temporary file {} -- you may need to delete it manually")

class DataSample(NtupleSample):
    def __init__(self, name, channel, dataIn, *args, **kwargs):
        super(DataSample, self).__init__(name, channel, dataIn, *args, **kwargs)


    def combineNtuples(self, files, chan):
        '''
        Combines a number of ntuples into a new ntuple in a temporary file, 
        ensuring that each event only appears once. For now, the version
        of the event chosen is just the first one seen.
        '''
        chain = TreeChain('{}/ntuple'.format(chan), files)
        if not hasattr(self, 'tempFile'):
            self.tempFile = _FileWrapper(TemporaryFile())
        self.tempFile.f.cd()
        
        out = Tree('{}_{}_ntuple'.format(self.name, chan))
        out.set_buffer(chain._buffer, create_branches=True)

        found = set()
        
        for i,ev in enumerate(chain):
            evID = (ev.run, ev.lumi, ev.evt)
            if evID in found:
                continue

            out.fill()
            found.add(evID)
    
        return out


    def formatDefault(self):
        self.format(True, drawstyle='PE', legendstyle='LPE')


    def makeHist(self, var, selection, binning, weight='', perUnitWidth=True, 
                 poissonErrors=False, **kwargs):
        h = super(DataSample, self).makeHist(var, selection, binning, weight,
                                             perUnitWidth, **kwargs)
        
        if poissonErrors:
            pois = h.poisson_errors()
            pois.SetTitle(self.prettyName)
            for a,b in self._format.iteritems():
                setattr(pois, a, b)
            return pois

        return h
