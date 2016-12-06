#############################################################################
#                                                                           #
#    Classes for objects that contain, manipulate, and use physics          #
#    samples in groups or other combinations.                               #
#                                                                           #
#    Author: Nate Woods, U. Wisconsin                                       #
#                                                                           #
#############################################################################


from Metadata.metadata import groupInfo as _groups
from Metadata.metadata import sampleInfo as _samples

from Sample import _SampleBase
from . import DataSample as _DataSample

from rootpy.plotting import Hist, Hist2D, HistStack

from numbers import Number as _Number

class SampleGroup(_SampleBase):
    '''
    A sample group where the samples are simply added together.
    '''
    def __init__(self, name, channel, samplesIn, initFromMetadata=False,
                 *args, **kwargs):
        '''
        samplesIn should be a dict where each sample has an identifying key.
        '''
        super(SampleGroup, self).__init__(name, channel, samplesIn,
                                          initFromMetadata, *args, **kwargs)

        self._recursePostprocessor = False


    def initFromMetadata(self):
        if self.name in _groups:
            info = _groups[self.name]
        elif self.name in _samples:
            info = _samples[self.name]
        else:
            raise ValueError("Sample {} is not in the metadata".format(self.name))

        self.shortName = info.get('shortName', self.name)
        self.prettyName = info.get('prettyName', self.shortName)

        self._format.update(info.get('format', {}))

        self.isSignal = info.get('isSignal',False)


    def storeInputs(self, inputs):
        self._samples = {k:s for k,s in inputs.iteritems()}


    def addSample(self, sampleName, sample):
        self._samples[sampleName] = sample


    def makeHist(self, var, selection, binning, weight='', perUnitWidth=True,
                 poissonErrors=False, postprocess=False,
                 mergeOverflow=False, **kwargs):
        '''
        If var is a dictionary, only samples with those keys will be used. If
            it is a string or appropriate iterable, it will be applied to all
            samples.
        If selection and/or weight are dictionaries, they will be applied only
            to samples of their keys. Any sample not in the dict will be given
            a null selection/weight (i.e., an empty string). If they are
            strings or appropriate iterables, they will be applied to all
            samples.
        If poissonErrors evaluates to True, all samples are required to be
            data samples, and a TGraphAsymmErrors is returned instead of a Hist
        If perUnitWidth is a positive number, bins are normalized to that
            width. If it is a non-number that evaluates to True, bins are
            normalized to their width in the units of the x-axis.
        Note: if the postprocessor was added recursively, it is run only here,
            not on the histograms made by the child samples
        '''
        if not len(self):
            raise KeyError(("Group {} can't be drawn because it contains no "
                            "samples.").format(self.name))

        if isinstance(var, dict):
            samplesToUse = [k for k in var.keys() if k in self._samples]
            if not len(samplesToUse):
                errMsg = ('No samples to draw!\n'
                          'You tried to draw ({}),\n'
                          'but group {} contains only ({}).')
                errMsg = errMsg.format(', '.join(var.keys()),
                                       self.name,
                                       ', '.join(self.keys()))
                raise KeyError(errMsg)
        else:
            samplesToUse = self.keys()
            var = {k:var for k in samplesToUse}

        # Do it this way to allow extra samples to be in the dict without
        # without modifying dictionaries passed in from elsewhere
        if isinstance(selection, dict):
            selections = {k:'' for k in samplesToUse} # note the "s"!
            selections.update(selection)
        else:
            selections = {k:selection for k in samplesToUse}
        if isinstance(weight, dict):
            weights = {k:'' for k in samplesToUse} # note the "s"!
            weights.update(weight)
        else:
            weights = {k:weight for k in samplesToUse}

        if poissonErrors:
            assert all(isinstance(self._samples[s], _DataSample) or isinstance(self._samples[s], SampleGroup) for s in samplesToUse), \
                "Poisson errors only make sense with data."
            h = sum(
                self._samples[s].makeHist(var[s], selections[s], binning,
                                          weights[s],
                                          perUnitWidth=False,
                                          poissonErrors=False,
                                          postprocess=(postprocess and not self._recursePostprocessor),
                                          mergeOverflow=mergeOverflow,
                                          **kwargs) for s in samplesToUse)
            out = h.poisson_errors()
            out.title = self.prettyName
            for a,b in self._format.iteritems():
                setattr(out,a,b)

            if perUnitWidth:
                if not isinstance(perUnitWidth, _Number):
                    perUnitWidth = 1.
                x = out.GetX()
                y = out.GetY()
                for i in xrange(out.GetN()):
                    width = (out.GetErrorXlow(i) + out.GetErrorXhigh(i)) / perUnitWidth
                    out.SetPoint(i, x[i], y[i] / width)
                    out.SetPointEYlow(i, out.GetErrorYlow(i) / width)
                    out.SetPointEYhigh(i, out.GetErrorYhigh(i) / width)

            if postprocess:
                self._postprocessor(out)

            return out

        bins = binning[:]
        if len(bins) != 3:
            bins = [bins]
        # use TH1D instead of TH1F because some datasets are now big enough for
        # floating point stuff to matter (!!)
        h = Hist(*bins, type='D', title=self.prettyName, **self._format)

        for s in samplesToUse:
            h += self._samples[s].makeHist(var[s], selections[s], binning,
                                           weights[s], perUnitWidth,
                                           mergeOverflow=mergeOverflow,
                                           **kwargs)

        if postprocess:
            self._postprocessor(h)

        return h


    def makeHist2(self, varX, varY, selection, binningX, binningY,
                  weight='', postprocess=False, mergeOverflowX=False,
                  mergeOverflowY=False, **kwargs):
        '''
        var[XY], selection and weight work the same as SampleGroup.makeHist()
        '''
        if not len(self):
            raise KeyError(("Group {} can't be drawn because it contains no "
                            "samples.").format(self.name))

        if isinstance(varX, dict):
            samplesToUse = [k for k in varX.keys() if k in self._samples]
            if not len(samplesToUse):
                errMsg = ('No samples to draw!\n'
                          'You tried to draw ({}),\n'
                          'but group {} contains only ({}).')
                errMsg = errMsg.format(', '.join(varX.keys()),
                                       self.name,
                                       ', '.join(self.keys()))
                raise KeyError(errMsg)
            if isinstance(varY, str):
                varY = {s:varY for s in samplesToUse}
            elif not all(k in varY for k in samplesToUse):
                raise KeyError('X and Y variables must be defined for the same'
                               ' samples')
        elif isinstance(varY, dict):
            samplesToUse = [k for k in varY.keys() if k in self._samples]
            if not len(samplesToUse):
                errMsg = ('No samples to draw!\n'
                          'You tried to draw ({}),\n'
                          'but group {} contains only ({}).')
                errMsg = errMsg.format(', '.join(varY.keys()),
                                       self.name,
                                       ', '.join(self.keys()))
                raise KeyError(errMsg)
            varX = {s:varX for s in samplesToUse} # better be a str!
        else:
            samplesToUse = self.keys()
            varX = {k:varX for k in samplesToUse}
            varY = {k:varY for k in samplesToUse}

        # Do it this way to allow extra samples to be in the dict without
        # without modifying dictionaries passed in from elsewhere
        if isinstance(selection, dict):
            selections = {k:'' for k in samplesToUse} # note the "s"!
            selections.update(selection)
        else:
            selections = {k:selection for k in samplesToUse}
        if isinstance(weight, dict):
            weights = {k:'' for k in samplesToUse} # note the "s"!
            weights.update(weight)
        else:
            weights = {k:weight for k in samplesToUse}

        binsX = binningX[:]
        if len(binsX) != 3:
            binsX = [binsX]
        binsY = binningY[:]
        if len(binsY) != 3:
            binsY = [binsY]
        binning = binsX + binsY

        # use TH2D instead of TH2F because some datasets are now big enough for
        # floating point stuff to matter (!!)
        h = Hist2D(*binning, type='D', title=self.prettyName, **self._format)

        for s in samplesToUse:
            h += self._samples[s].makeHist2(varX[s], varY[s], selections[s],
                                            binningX, binningY, weights[s],
                                            postprocess=postprocess,
                                            mergeOverflowX=mergeOverflowX,
                                            mergeOverflowY=mergeOverflowY,
                                            **kwargs)

        if postprocess:
            self._postprocessor(h)

        return h


    def formatDefault(self):
        self.format(True, drawstyle='hist', fillstyle='solid', legendstyle='F')


    def applyWeight(self, w, reset=False):
        '''
        Apply weight w to samples in the group.
        If w is a dict, each value is taken to be the weight for the sample in
        the group with the same key. Otherwise, w is applied to all samples.
        '''
        if isinstance(w, dict):
            for name, wt in w.iteritems():
                if name in self._samples:
                    self[name].applyWeight(wt, reset)
        else:
            for s in self.values():
                s.applyWeight(w, reset)


    def __iter__(self):
        for s in self._samples:
            yield s


    def keys(self):
        return self._samples.keys()


    def values(self):
        return self._samples.values()


    def itersamples(self):
        for x in self._samples.iteritems():
            yield x


    def __getitem__(self, x):
        return self._samples[x]


    def __len__(self):
        return len(self._samples)


    def setPostprocessor(self, f, recursive=False):
        '''
        If recursive is True, this is also passed down to child samples,
        though this.makeHist won't use the postprocessor on the children
        to avoid running it twice.
        '''
        super(SampleGroup, self).setPostprocessor(f)

        self._recursePostprocessor = recursive

        if recursive:
            for s in self._samples.values():
                s.setPostprocessor(f, True)



class SampleStack(_SampleBase):
    '''
    A sample group where the samples are plotted in a stack
    '''
    def __init__(self, name, channel, samplesIn, *args, **kwargs):
        '''
        samplesIn should be an iterable of samples
        '''
        super(SampleStack, self).__init__(name, channel, samplesIn,
                                          initFromMetadata=False,
                                          *args, **kwargs)
        self._format['drawstyle'] = 'histnoclear'

        self._recursePostprocessor = False


    def storeInputs(self, inputs):
        self._samples = inputs[:]


    def addSample(self, sample):
        self._samples.append(sample)


    def makeHist(self, var, selection, binning, weight='', perUnitWidth=True,
                 postprocess=False, mergeOverflow=False,
                 *extraHists, **kwargs):
        '''
        Note: if the postprocessor was added recursively, it is run only here,
        not on the histograms made by the child samples.
        '''
        sortByMax = kwargs.pop('sortByMax', True)

        sig = []
        bkg = []
        for s in self._samples:
            h = s.makeHist(var, selection, binning, weight,
                           perUnitWidth,
                           postprocess=(postprocess and not self._recursePostprocessor),
                           mergeOverflow=mergeOverflow,
                           **kwargs)
            try:
                isSignal = s.isSignal
            except AttributeError:
                isSignal = False
            if isSignal:
                sig.append(h)
            else:
                bkg.append(h)

        hists = SampleStack.orderForStack(list(extraHists)+bkg)
        hists += SampleStack.orderForStack(sig)

        stack = HistStack(hists, drawstyle='histnoclear')

        if postprocess:
            self._postprocessor(h)

        return stack


    def makeHist2(self, varX, varY, selection, binningX, binningY,
                  weight='', postprocess=False, mergeOverflowX=False,
                  mergeOverflowY=False, *extraHists, **kwargs):
        sortByMax = kwargs.pop('sortByMax', True)

        sig = []
        bkg = []
        for s in self._samples:
            h = s.makeHist2(varX, varY, selection, binningX, binningY, weight,
                            postprocess=(postprocess and not self._recursePostprocessor),
                            mergeOverflowX=mergeOverflowX,
                            mergeOverflowY=mergeOverflowY,
                            **kwargs)
            try:
                isSignal = s.isSignal
            except AttributeError:
                isSignal = False
            if isSignal:
                sig.append(h)
            else:
                bkg.append(h)

        hists = SampleStack.orderForStack(list(extraHists)+bkg)
        hists += SampleStack.orderForStack(sig)

        stack = HistStack(hists, drawstyle='histnoclear')

        if postprocess:
            self._postprocessor(h)

        return stack


    def applyWeight(self, w, reset=False):
        '''
        Apply weight w to all samples in the stack.
        '''
        for s in self._samples:
            s.applyWeight(w, reset)


    @staticmethod
    def orderForStack(hists, sortByMax=True):
        '''
        Sort a list of histograms so that the largest is listed last.
        If sortByMax determines whether the sort variable is the histogram's
        largest or smallest bin.
        '''
        if sortByMax:
            binGetter = Hist.GetMaximumBin
        else:
            binGetter = Hist.GetMinimumBin

        key = lambda h: (h.GetBinContent(binGetter(h)))

        hists.sort(key=key)

        return hists


    def __iter__(self):
        for s in self._samples:
            yield s


    def __getitem__(self, n):
        return self._samples[n]


    def getSamplesForChannel(self, channel):
        '''
        Get a list of all samples in the stack for one channel.
        All samples in the stack must be groups which include a sample for
        this channel.
        '''
        return [s[channel] for s in self]


    def __len__(self):
        return len(self._samples)


    def setPostprocessor(self, f, recursive=False):
        '''
        If recursive is True, this is also passed down to child samples,
        though this.makeHist won't use the postprocessor on the children
        to avoid running it twice.
        '''
        super(SampleStack, self).setPostprocessor(f)

        self._recursePostprocessor = recursive

        if recursive:
            for s in self:
                s.setPostprocessor(f, True)
