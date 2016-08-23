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

from rootpy.plotting import Hist, HistStack


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
                 poissonErrors=False, **kwargs):
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
            out = sum(self._samples[s].makeHist(var[s], selection[s], binning,
                                                weight[s], perUnitWidth, 
                                                poissonErrors=True,
                                                **kwargs) for s in samplesToUse)
            out.title = self.prettyName
            return out

        # use TH1D instead of TH1F because some datasets are now big enough for
        # floating point stuff to matter (!!)
        h = Hist(*binning, type='D', title=self.prettyName, **self._format)

        for s in samplesToUse:
            h += self._samples[s].makeHist(var[s], selection[s], binning,
                                           weight[s], perUnitWidth, **kwargs)

        return h


    def formatDefault(self):
        self.format(True, drawstyle='hist', fillstyle='solid', legendstyle='F')



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
                 *extraHists, **kwargs):
        sortByMax = kwargs.pop('sortByMax', True)

        sig = []
        bkg = []
        for s in self._samples:
            h = s.makeHist(var, selection, binning, weight, 
                           perUnitWidth, **kwargs)
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

        # # workaround for weird root bug
        # emptyHist = hists[0].empty_clone()
        # hists.append(emptyHist)

        stack = HistStack(hists, drawstyle='histnoclear')

        return stack


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
