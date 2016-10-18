
import logging
from rootpy import log as rlog; rlog = rlog["/zzPlots"]
# don't show most silly ROOT messages
logging.basicConfig(level=logging.WARNING)
rlog["/ROOT.TUnixSystem.SetDisplay"].setLevel(rlog.ERROR)

from rootpy.io import root_open
from rootpy.plotting import Canvas, Legend
from rootpy.plotting.utils import draw
from rootpy.ROOT import TBox, Double

from SampleTools import MCSample, DataSample, SampleGroup, SampleStack
from PlotTools import PlotStyle as _Style
from PlotTools import makeLegend, addPadBelow, makeRatio, fixRatioAxes
from Utilities import WeightStringMaker
from Analysis import standardZZSamples

from os import environ
from os import path as _path
from os import makedirs as _mkdir
from collections import OrderedDict
from math import sqrt


inData = 'uwvvNtuples_data_12oct2016'
inMC = 'uwvvNtuples_mc_12oct2016'

puWeightFile = 'puWeight_69200_08sep2016'
fakeRateFile = 'fakeRate_08sep2016'

ana = 'smp'

outdir = '/afs/cern.ch/user/n/nawoods/www/UWVVPlots/zz_{}'.format(ana)
try:
    _mkdir(outdir)
except OSError: # already exists
    pass

style = _Style()

lumi = 15937.

amcatnlo=False
if amcatnlo:
    outdir += '_amcatnlo'

data, stack = standardZZSamples('zz', inData, inMC, ana, puWeightFile,
                                fakeRateFile, lumi, amcatnlo=amcatnlo)

# count events
tot = OrderedDict()
totErrSqr = {}
for c in ['eeee','eemm','mmmm']:
    print c + ':'

    expected = 0.
    expErrSqr = 0.
    for s in stack:
        h = s[c].makeHist('1.', '', [1,0.,2.])
        yErr = Double(0)
        y = h.IntegralAndError(0, h.GetNbinsX(), yErr)
        if s.name == 'Z+X':
            yErr = sqrt(yErr ** 2 + (.4 * y) ** 2)
        print '    {}: {} +/- {}'.format(s.name, y, yErr)
        if s.name not in tot:
            tot[s.name] = 0.
            totErrSqr[s.name] = 0.
        tot[s.name] += y
        totErrSqr[s.name] += yErr ** 2

        expected += y
        expErrSqr += yErr ** 2
    print '    Total expected: {} +/- {}'.format(expected, sqrt(expErrSqr))
    if 'expected' not in tot:
        tot['expected'] = 0.
        totErrSqr['expected'] = 0.
    tot['expected'] += expected
    totErrSqr['expected'] += expErrSqr

    hData = data[c].makeHist('1.', '', [1,0.,2.])
    dataErr = Double(0)
    yieldData = hData.IntegralAndError(0, hData.GetNbinsX(), dataErr)
    print '    Data: {} +/- {}'.format(yieldData, dataErr)
    if 'data' not in tot:
        tot['data'] = 0.
        totErrSqr['data'] = 0.
    tot['data'] += yieldData
    totErrSqr['data'] += dataErr ** 2

    print ''

print 'Total:'
for n,t in tot.iteritems():
    print '    {}: {} +/- {}'.format(n,t,sqrt(totErrSqr[n]))


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
    'Mass'  : [100.] + [200.+50.*i for i in range(5)] + [500.,600.,800.], #[26, 150., 700.],#[35, 250., 2000.],
    'Pt'    : [20.*i for i in range(4)] + [100., 140., 200., 300.],#[20.*i for i in range(4)] + [100., 140., 200., 300.], #[40, 0., 200.],
    'Eta'   : [16, -5., 5.],
    'Phi'   : [12, -3.15, 3.15],
    'nvtx'  : [40, 0., 40.],
    'nJets' : [6, -0.5, 5.5],
    }

if amcatnlo:
    binning4l = {'Mass' : binning4l['Mass'], 'Pt' : binning4l['Pt']}

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
            dataSelection = 'Mass < 800.'

        hStack = stack.makeHist(var, '', binning, postprocess=True)
        dataPts = data.makeHist(var, dataSelection, binning, poissonErrors=True)

        c = Canvas(1000,1000)

        leg = makeLegend(c, hStack, dataPts)

        (xaxis, yaxis), (xmin,xmax,ymin,ymax) = draw([hStack, dataPts], c,
                                                     xtitle='{}'.format(varName)+(' ({})'.format(units[varName]) if units[varName] else ''),
                                                     ytitle='Events')
        # blinding box
        if varName == 'Mass' and binning4l['Mass'][-1] > 800.:
            box = TBox(max(xmin,800.), ymin, min(binning4l['Mass'][-1], xmax), ymax)
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

if amcatnlo:
    binning2l = {}

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


binning1l = {
    'Pt' : [20, 0., 200.],
    'Eta' : [20, -2.5, 2.5],
    'Phi' : [24, -3.15, 3.15],
    'Iso' : [8, 0., .4],
    'PVDXY' : [20, -.1, .1],
    'PVDZ' : [20, -.2, .2],
    'SIP3D' : [20, 0., 5.],
    }

if amcatnlo:
    binning1l = {}

ze1LepVarTemp = ['e1{var}','e2{var}']
ze2LepVarTemp = ['e3{var}','e4{var}']
zm1LepVarTemp = ['m1{var}','m2{var}']
zm2LepVarTemp = ['m3{var}','m4{var}']
varTemplates1l = {
    'l' : {
        'eeee' : ze1LepVarTemp + ze2LepVarTemp,
        'eemm' : ze1LepVarTemp + zm1LepVarTemp,
        'mmmm' : zm1LepVarTemp + zm2LepVarTemp,
        },
    'l1' : {
        'eeee' : [ze1LepVarTemp[0], ze2LepVarTemp[0]],
        'eemm' : [ze1LepVarTemp[0], zm1LepVarTemp[0]],
        'mmmm' : [zm1LepVarTemp[0], zm2LepVarTemp[0]],
        },
    'e' : {
        'eeee' : ze1LepVarTemp + ze2LepVarTemp,
        'eemm' : ze1LepVarTemp,
        },
    'm' : {
        'mmmm' : zm1LepVarTemp + zm2LepVarTemp,
        'eemm' : zm1LepVarTemp,
        },
    }

#vars1l = {l:{v:{c:[vt.format(v) for vt in varTemplates1l[l][c]] for c in varTemplates1l[l]} for v in binning1l} for l in varTemplates1l}

selections1l = {l:'' for l in varTemplates1l}
selections1l['l1'] = {
    'eeee' : ['e1Pt > e3Pt', 'e3Pt > e1Pt'],
    'eemm' : ['e1Pt > m1Pt', 'm1Pt > e1Pt'],
    'mmmm' : ['m1Pt > m3Pt', 'm3Pt > m1Pt'],
    }

for lep in varTemplates1l:
    for varName, binning in binning1l.iteritems():
        print "Plotting {} {}".format(lep, varName)

        if varName == 'Iso':
            varStr = 'ZZIso'
        else:
            varStr = varName

        var = {c:[vt.format(var=varStr) for vt in varTemplates1l[lep][c]] for c in varTemplates1l[lep]}

        hStack = stack.makeHist(var, selections1l[lep], binning1l[varName], postprocess=True)
        dataPts = data.makeHist(var, selections1l[lep], binning1l[varName], poissonErrors=True)

        # for ratio
        dataHist = data.makeHist(var, selections1l[lep], binning1l[varName])

        c = Canvas(1000,1000)

        leg = makeLegend(c, hStack, dataPts)

        (xaxis, yaxis), (xmin,xmax,ymin,ymax) = draw([hStack, dataPts], c,
                                                     xtitle='{}'.format(varName)+(' ({})'.format(units[varName]) if units[varName] else ''),
                                                     ytitle='Leptons')
        leg.Draw("same")

        style.setCMSStyle(c, '', dataType='Preliminary', intLumi=lumi)
        c.Print('{}/{}{}.png'.format(outdir, lep, varName))
