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
from rootpy.plotting import Hist, Hist2D, Canvas
from rootpy.context import preserve_current_directory

from glob import glob
from math import sqrt
from os import remove as _rm
from os import close as _close
from os.path import isfile as _isfile
from numbers import Number as _Number

# Workaround for weird ROOT bug
_dummy = Hist(1,0,1)
_cDummy = Canvas(10,10)
_dummy.xaxis.SetTitle('\\mu\\mu')
_dummy.draw()


class _SampleBase(object):
    def __init__(self, name, channel, dataIn,
                 initFromMetadata=False, *args, **kwargs):
        self.name = name
        self.channel = channel
        self.shortName = name
        self.prettyName = name

        self.formatDefault()

        if initFromMetadata:
            self.initFromMetadata()

        self.storeInputs(dataIn)

        self._postprocessor = lambda *args, **kwargs: None


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

    def setPostprocessor(self, f, *args):
        '''
        Set a function that is always run on a histogram made by this object.
        '''
        if not hasattr(f, '__call__'):
            raise ValueError('Failed to set postprocessor: '
                             '{} is not callable'.format(f))

        self._postprocessor = f


class NtupleSample(_SampleBase):
    def __init__(self, name, channel, dataIn, initFromMetadata=False,
                 *args, **kwargs):
        self.weight = ''

        super(NtupleSample, self).__init__(name, channel, dataIn, initFromMetadata, *args, **kwargs)

        self.oldNtuples = []

    def storeInputs(self, inputs):
        if isinstance(inputs, Tree):
            self.ntuple = inputs
            self.files = [self.ntuple.GetCurrentFile().GetName()]
            return self.ntuple

        if isinstance(inputs, str):
            inputs = [inputs]

        self.files = []
        for inp in inputs:
            self.files += glob(inp)

        if not self.files:
            raise IOError(("No files found for {}.\n"
                           "Paths checked:\n{}").format(self.name,
                                                        '\n'.join(inputs))
                          )

        self.ntuple = self.combineNtuples(self.files, self.channel)
        return self.ntuple

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

    def makeHist(self, var, selection, binning, weight='', perUnitWidth=True,
                 postprocess=False, mergeOverflow=False, **kwargs):
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
            w = combineWeights(w, self.fullWeight())

            s = combineWeights(w, s)

            self.addToHist(h, v, s)

        h.sumw2()

        if mergeOverflow:
            h = h.merge_bins([(-2,-1)])
            h.sumw2()
            # have to reformat
            h.SetTitle(self.prettyName)
            for a,b in self._format.iteritems():
                setattr(h, a, b)

        if perUnitWidth:
            if not isinstance(perUnitWidth, _Number):
                perUnitWidth = 1.
            for ib in xrange(1,len(h)+1):
                w = h.GetBinWidth(ib) / perUnitWidth
                h.SetBinContent(ib, h.GetBinContent(ib) / w)
                h.SetBinError(ib, h.GetBinError(ib) / w)
            h.sumw2()

        if postprocess:
            self._postprocessor(h)

        return h


    def makeHist2(self, varX, varY, selection, binningX, binningY, weight='',
                  postprocess=False, mergeOverflowX=False,
                  mergeOverflowY=False, **kwargs):
        '''
        If var[XY], selection, and/or weight are iterables (which must be of the
        same length), the hists from the resulting var/selection/weight sets
        will be added together. If some are strings, they are reused the
        appropriate number of times.
        '''
        if len(binningX) != 3:
            binningX = [binningX]
        if len(binningY) != 3:
            binningY = [binningY]
        binning = binningX + binningY

        # use TH2D instead of TH2F because some datasets are now big enough for
        # floating point stuff to matter (!!)
        h = Hist2D(*binning, type='D', title=self.prettyName, drawstyle='colz')

        nToAdd = max(1 if isinstance(varX,str) else len(varX),
                     1 if isinstance(varY,str) else len(varY),
                     1 if isinstance(selection,str) else len(selection),
                     1 if isinstance(weight,str) else len(weight))
        if isinstance(varX, str):
            varX = nToAdd * [varX]
        if isinstance(varY, str):
            varY = nToAdd * [varY]
        if isinstance(selection, str):
            selection = nToAdd * [selection]
        if isinstance(weight, str):
            weight = nToAdd * [weight]

        assert len(varX) == len(varY) and len(varY) == len(selection) and len(selection) == len(weight), \
            "Invalid plotting parameters! Variables: {} and {}, selection: {}, weight: {}.".format(varX,varY,selection,weight)

        for vx, vy, s, w in zip(varX, varY, selection, weight):
            v = '{}:{}'.format(vx,vy)

            w = combineWeights(w, self.fullWeight())

            s = combineWeights(w, s)

            self.addToHist(h, v, s)

        h.sumw2()

        if mergeOverflowX:
            h = h.merge_bins([(-2,-1)])
            h.sumw2()
            if not mergeOverflowY: # this will happen there anyway
                h.SetTitle(self.prettyName)
                for a,b in self._format.iteritems():
                    setattr(h, a, b)
        if mergeOverflowY:
            h = h.merge_bins([(-2,-1)], 1)
            h.sumw2()
            # have to reformat
            h.SetTitle(self.prettyName)
            for a,b in self._format.iteritems():
                setattr(h, a, b)


        if postprocess:
            self._postprocessor(h)

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
        '''
        Get the weight not including the automatically set weights (cross
        section etc.). I.e., only the weights set by the user.
        '''
        return self._weight
    @weight.setter
    def weight(self, w):
        if not isinstance(w, str):
            w = str(w)
        self._weight = w
    def fullWeight(self):
        '''
        Get the weight including the automatically set weights (cross
        section etc.).
        '''
        return combineWeights(self.weight, self.implicitWeight())

    def applyCut(self, cut='', name=''):
        '''
        cut (str, Cut, function, or other callable): selection for rows to copy.
            If a string or Cut, the Tree will be copied with CopyTree using
            that cut. If a function or other callable, each row will be copied
            if cut(row) evaluates to True.
        name (str): name for new tree (defaults to name of old tree. Only used
            if cut is callable).
        '''
        self.oldNtuples.append(self.ntuple)
        if hasattr(cut, '__call__'):
            if not name:
                name = self.oldNtuples[-1].GetName()
            self.ntuple = Tree(name)
            self.ntuple.set_buffer(self.oldNtuples[-1]._buffer,
                                   create_branches=True)
            for row in self.oldNtuples[-1]:
                if cut(row):
                    self.ntuple.fill()
        else:
            self.ntuple = asrootpy(self.oldNtuples[-1].CopyTree(cut))

        return self.ntuple


    def __iter__(self):
        for row in self.ntuple:
            yield row


    def rows(self):
        for row in self:
            yield row


    def __len__(self):
        return self.ntuple.GetEntries()


    def getFileNames(self):
        for fn in self.files:
            yield fn



class MCSample(NtupleSample):
    def __init__(self, name, channel, dataIn, initFromMetadata=False, intLumi=1000, *args, **kwargs):
        self.isSignal = 0
        self.sumW = -1
        self.xsec = -1
        self.kFactor = '1'
        self.intLumi = intLumi
        super(MCSample, self).__init__(name, channel, dataIn, initFromMetadata, *args, **kwargs)


    def initFromMetadata(self):
        super(MCSample, self).initFromMetadata()

        info = _samples[self.name]

        self.isSignal = info.get('isSignal', 0)
        self.sumW = info.get('sumW', 0)
        self.xsec = info.get('xsec', -1)
        self.kFactor = info.get('kFactor', '1.')


    def storeInputs(self, inputs):
        stored = super(MCSample, self).storeInputs(inputs)

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

        return stored


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
            if _isfile(self.name):
                _rm(self.name)
        except OSError:
            rlog.warning("Failed to delete temporary file {} -- you may "
                         "need to delete it manually".format(self.name))

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
                 poissonErrors=False, postprocess=False, mergeOverflow=False,
                 **kwargs):
        h = super(DataSample, self).makeHist(var, selection, binning, weight,
                                             perUnitWidth=(perUnitWidth if not poissonErrors else False),
                                             **kwargs)

        if mergeOverflow:
            h = h.merge_bins([(-2,-1)])
            h.sumw2()
            if not poissonErrors: # this will happen anyway
                h.SetTitle(self.prettyName)
                for a,b in self._format.iteritems():
                    setattr(h, a, b)

        if poissonErrors:
            pois = h.poisson_errors()
            # have to reformat
            pois.SetTitle(self.prettyName)
            for a,b in self._format.iteritems():
                setattr(pois, a, b)

            if perUnitWidth:
                if not isinstance(perUnitWidth, _Number):
                    perUnitWidth = 1.
                x = pois.GetX()
                y = pois.GetY()
                for i in xrange(pois.GetN()):
                    width = (pois.GetErrorXlow(i) + pois.GetErrorXhigh(i)) / perUnitWidth
                    pois.SetPoint(i, x[i], y[i] / width)
                    pois.SetPointEYlow(i, pois.GetErrorYlow(i) / width)
                    pois.SetPointEYhigh(i, pois.GetErrorYhigh(i) / width)

            if postprocess:
                self._postprocessor(pois)

            return pois

        if postprocess:
            self._postprocessor(h)

        return h

