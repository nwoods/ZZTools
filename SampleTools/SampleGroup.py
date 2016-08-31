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


    def storeInputs(self, inputs):
        self._samples = {k:s for k,s in inputs.iteritems()}


    def addSample(self, sampleName, sample):
        self._samples[sampleName] = sample


    def makeHist(self, var, selection, binning, weight='', perUnitWidth=True, 
                 poissonErrors=False, postprocess=False, **kwargs):
        '''
        If var, selection, and/or weight are dictionaries, they will be used to
        add histograms from the sample of the same key. If they are strings or 
        appropriate iterables, they will be used to add histograms from all 
        samples.
        If poissonErrors evaluates to True, all samples are required to be 
        data samples, and a TGraphAsymmErrors is returned instead of a Hist
        '''
        if isinstance(var, dict):
            samplesToUse = var.keys()
        elif isinstance(selection, dict):
            samplesToUse = selection.keys()
        elif isinstance(weight, dict):
            samplesToUse = weight.keys()
        else:
            samplesToUse = self._samples.keys()

        if not isinstance(var, dict):
            var = {k:var for k in samplesToUse}
        if not isinstance(selection, dict):
            selection = {k:selection for k in samplesToUse}
        if not isinstance(weight, dict):
            weight = {k:weight for k in samplesToUse}

        if poissonErrors:
            assert all(isinstance(self._samples[s], _DataSample) or isinstance(self._samples[s], SampleGroup) for s in samplesToUse), \
                "Poisson errors only make sense with data."
            h = sum(self._samples[s].makeHist(var[s], selection[s], binning,
                                              weight[s], perUnitWidth, 
                                              poissonErrors=False,
                                              postprocess=postprocess,
                                              **kwargs) for s in samplesToUse)
            out = h.poisson_errors()
            out.title = self.prettyName
            for a,b in self._format.iteritems():
                setattr(out,a,b)

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
            h += self._samples[s].makeHist(var[s], selection[s], binning,
                                           weight[s], perUnitWidth, **kwargs)

        if postprocess:
            self._postprocessor(h)

        return h


    def makeHist2(self, varX, varY, selection, binningX, binningY, 
                  weight='', postprocess=False, **kwargs):
        '''
        If var[XY], selection, and/or weight are dictionaries, they will be 
        used to add histograms from the sample of the same key. If they are 
        strings or appropriate iterables, they will be used to add histograms 
        from all samples.
        '''
        if isinstance(varX, dict):
            samplesToUse = varX.keys()
            if not isinstance(varY, str):
                assert all(k in varY for k in varX), "X and Y variables must match"
        elif isinstance(varY, dict):
            samplesToUse = varY.keys()
        elif isinstance(selection, dict):
            samplesToUse = selection.keys()
        elif isinstance(weight, dict):
            samplesToUse = weight.keys()
        else:
            samplesToUse = self._samples.keys()

        if not isinstance(varX, dict):
            varX = {k:varX for k in samplesToUse}
        if not isinstance(varY, dict):
            varY = {k:varY for k in samplesToUse}
        if not isinstance(selection, dict):
            selection = {k:selection for k in samplesToUse}
        if not isinstance(weight, dict):
            weight = {k:weight for k in samplesToUse}

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
            h += self._samples[s].makeHist2(varX[s], varY[s], selection[s], 
                                            binningX, binningY, weight[s], 
                                            postprocess=postprocess, 
                                            **kwargs)

        if postprocess:
            self._postprocessor(h)

        return h


    def formatDefault(self):
        self.format(True, drawstyle='hist', fillstyle='solid', legendstyle='F')


    def applyWeight(self, w, reset=False):
        '''
        Apply weight w to all samples in the group.
        '''
        for s in self._samples.values():
            s.applyWeight(w, reset)


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


    def storeInputs(self, inputs):
        self._samples = [s for s in inputs]


    def addSample(self, sample):
        self._samples.append(sample)


    def makeHist(self, var, selection, binning, weight='', perUnitWidth=True, 
                 postprocess=False, *extraHists, **kwargs):
        sortByMax = kwargs.pop('sortByMax', True)

        sig = []
        bkg = []
        for s in self._samples:
            h = s.makeHist(var, selection, binning, weight, 
                           perUnitWidth, postprocess=postprocess, **kwargs)
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
                  weight='', postprocess=False, *extraHists, **kwargs):
        sortByMax = kwargs.pop('sortByMax', True)

        sig = []
        bkg = []
        for s in self._samples:
            h = s.makeHist2(varX, varY, selection, binningX, binningY, weight, 
                            postprocess=postprocess, **kwargs)
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

        # lazy way to prevent data hist from reaching sampleInfo query
        key = lambda h: (h.GetBinContent(binGetter(h)))

        hists.sort(key=key)

        return hists
