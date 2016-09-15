
import logging
from rootpy import log as rlog; rlog = rlog["/zzPlots"]
# don't show most silly ROOT messages
logging.basicConfig(level=logging.WARNING)
rlog["/ROOT.TUnixSystem.SetDisplay"].setLevel(rlog.ERROR)

from rootpy.io import root_open
from rootpy.plotting import Canvas, Legend
from rootpy.plotting.utils import draw
from rootpy.ROOT import TBox

from SampleTools import MCSample, DataSample, SampleGroup, SampleStack
from PlotTools import PlotStyle as _Style
from PlotTools import makeLegend, addPadBelow, makeRatio, fixRatioAxes
from Utilities import WeightStringMaker
from Analysis import standardZZSamples

from os import environ
from os import path as _path
from os import makedirs as _mkdir



inData = 'uwvvNtuples_data_08sep2016'
inMC = 'uwvvNtuples_mc_08sep2016'

puWeightFile = 'puWeight_69200_08sep2016'
fakeRateFile = 'fakeRate_08sep2016'

ana = 'z4l' #'smp'

outdir = '/afs/cern.ch/user/n/nawoods/www/UWVVPlots/zz_{}'.format(ana)
try:
    _mkdir(outdir)
except OSError: # already exists
    pass

style = _Style()

lumi = 15937.

data, stack = standardZZSamples(inData, inMC, ana, puWeightFile, 
                                fakeRateFile, lumi)

### Set up variable specific info

units = {
    'Mass' : 'GeV',
    'Eta' : '',
    'Phi' : '',
    'Pt' : 'GeV',
    'nJets' : '',
    'Iso' : '',
    'PVDXY' : 'cm',
    'PVDZ' : 'cm',
    'nvtx' : '',
    'SIP3D' : '',
    }

binning4l = {
    'Mass'  : [26, 150., 700.],#[35, 250., 2000.],
    'Pt'    : [20.*i for i in range(4)] + [100., 140., 200., 300.], #[40, 0., 200.],
    'Eta'   : [16, -5., 5.],
    'Phi'   : [12, -3.15, 3.15],
    'nvtx'  : [40, 0., 40.],
    'nJets' : [6, -0.5, 5.5],
    }

if ana == 'z4l':
    binning4l['Mass'] = [20, 80., 100.]

for chan in ['zz', 'eeee', 'eemm', 'mmmm']:
    for varName, binning in binning4l.iteritems():
        print "Plotting {} {}".format(chan, varName)

        if chan == 'zz':
            var = varName
        else:
            var = {chan:varName}

        # blinding
        dataSelection = ''
        if varName == 'Mass':
            dataSelection = 'Mass < 700.'

        hStack = stack.makeHist(var, '', binning, postprocess=True)
        dataPts = data.makeHist(var, dataSelection, binning, poissonErrors=True)
        
        c = Canvas(1000,1000)

        leg = makeLegend(c, hStack, dataPts)

        (xaxis, yaxis), (xmin,xmax,ymin,ymax) = draw([hStack, dataPts], c, 
                                                     xtitle='{}'.format(varName)+(' ({})'.format(units[varName]) if units[varName] else ''), 
                                                     ytitle='Events')
        # blinding box
        if varName == 'Mass' and binning4l['Mass'][-1] > 700.:
            box = TBox(max(xmin,700.), ymin, min(binning4l['Mass'][-1], xmax), ymax)
            box.SetFillColor(1)
            box.SetFillStyle(3002)
            box.Draw("same")
            leg.SetFillStyle(1001)

        leg.Draw("same")

        style.setCMSStyle(c, '', dataType='Preliminary', intLumi=lumi)
        c.Print('{}/{}{}.png'.format(outdir, chan, varName))



binning2l = {
    'Mass' : [60, 60., 120.],
    'Pt' : [60, 0., 300.],
    'Eta' : [48,-6.,6.],
    'Phi' : [24, -3.15,3.15],
    }

if ana != 'smp':
    binning2l['Mass'] = [60, 0., 120.]

ze1VarTemp = 'e1_e2_{var}'
ze2VarTemp = 'e3_e4_{var}'
zm1VarTemp = 'm1_m2_{var}'
zm2VarTemp = 'm3_m4_{var}'
varTemplates2l = {
    'z' : {
        'eeee' : [ze1VarTemp, ze2VarTemp],
        'eemm' : [ze1VarTemp, zm1VarTemp],
        'mmmm' : [zm1VarTemp, zm2VarTemp],
        },
    'z1' : {
        'eeee' : [ze1VarTemp],
        'eemm' : [ze1VarTemp, zm1VarTemp],
        'mmmm' : [zm1VarTemp],
        },
    'z2' : {
        'eeee' : [ze2VarTemp],
        'eemm' : [ze1VarTemp, zm1VarTemp],
        'mmmm' : [zm2VarTemp],
        },
    'ze' : {
        'eeee' : [ze1VarTemp, ze2VarTemp],
        'eemm' : [ze1VarTemp],
        },
    'zm' : {
        'mmmm' : [zm1VarTemp, zm2VarTemp],
        'eemm' : [zm1VarTemp],
        },
    }

selections2l = {z:'' for z in varTemplates2l}
selections2l['z1'] = {
    'eeee' : '',
    'mmmm' : '',
    'eemm' : ['abs(e1_e2_Mass - 91.1876) < abs(m1_m2_Mass - 91.1876)',
              'abs(e1_e2_Mass - 91.1876) > abs(m1_m2_Mass - 91.1876)']
    }
selections2l['z2'] = {
    'eeee' : '',
    'mmmm' : '',
    'eemm' : ['abs(e1_e2_Mass - 91.1876) > abs(m1_m2_Mass - 91.1876)',
              'abs(e1_e2_Mass - 91.1876) < abs(m1_m2_Mass - 91.1876)']
    }


for z in ['z', 'ze', 'zm', 'z1', 'z2']:
    for varName, binning in binning2l.iteritems():
        print "Plotting {} {}".format(z, varName)

        var = {c:[vt.format(var=varName) for vt in varTemplates2l[z][c]] for c in varTemplates2l[z]}

        hStack = stack.makeHist(var, selections2l[z], binning2l[varName], postprocess=True)
        dataPts = data.makeHist(var, selections2l[z], binning2l[varName], poissonErrors=True)
        
        # for ratio
        dataHist = data.makeHist(var, '', binning2l[varName])

        c = Canvas(1000,1000)

        leg = makeLegend(c, hStack, dataPts)

        (xaxis, yaxis), (xmin,xmax,ymin,ymax) = draw([hStack, dataPts], c, 
                                                     xtitle='{}'.format(varName)+(' ({})'.format(units[varName]) if units[varName] else ''), 
                                                     ytitle='Events')
        leg.Draw("same")

        style.setCMSStyle(c, '', dataType='Preliminary', intLumi=lumi)
        c.Print('{}/{}{}.png'.format(outdir, z, varName))

                       
