
import logging
from rootpy import log as rlog; rlog = rlog["/singleZPlots"]
# don't show most silly ROOT messages
logging.basicConfig(level=logging.WARNING)
rlog["/ROOT.TUnixSystem.SetDisplay"].setLevel(rlog.ERROR)

from rootpy.io import root_open
from rootpy.plotting import Canvas, Legend
from rootpy.plotting.utils import draw

from SampleTools import MCSample, DataSample, SampleGroup, SampleStack
from PlotTools import PlotStyle as _Style
from PlotTools import makeLegend, addPadBelow, makeRatio, fixRatioAxes
from Utilities import WeightStringMaker
from Analysis.weightHelpers import baseMCWeight

from os import environ
from os import path as _path



test = False

sfFromHists = True

outdir = '/afs/cern.ch/user/n/nawoods/www/UWVVPlots/singleZ_24jan2017_sfFromHists'
if test:
    outdir = '/afs/cern.ch/user/n/nawoods/www/UWVVPlots/test'

style = _Style()

lumi = 35860.

channels = ['ee','mm']

puWeightFile = 'puWeight_69200_24jan2017.root'

mcWeight = baseMCWeight('ee,mm', puWeightFile,
                        scaleFactorsFromHists=sfFromHists)

if test:
    dyByChan = {
        c : MCSample('DYJets', c,
                     '/data/nawoods/ntuples/uwvvSingleZ_mc_test/DYJets*.root',
                     True, lumi) for c in channels
        }
else:
    dyByChan = {}
    for c in channels:
        byNJets = {
            '{}J'.format(nJ) : MCSample('DYToLL-{}J'.format(nJ), c,
                                        '/data/nawoods/ntuples/uwvvSingleZ_mc_24jan2017/results/DYToLL-{}J_*.root'.format(nJ),
                                        True, lumi) for nJ in range(3)}
        dyByChan[c] = SampleGroup('DYJets', c, byNJets, True)

    # dyByChan = {
    #     c : MCSample('DYJets', c,
    #                  '/data/nawoods/ntuples/uwvvSingleZ_mc_24jan2017/results/DYJets*.root',
    #                  True, lumi) for c in channels
    #     }

dy = SampleGroup('DYJets', 'z', dyByChan, True)

if test:
    ttByChan = {
        c : MCSample('TTJets', c,
                     '/data/nawoods/ntuples/uwvvSingleZ_mc_test/TTJets*.root',
                     True, lumi) for c in channels
        }
else:
    ttByChan = {
        c : MCSample('TTTo2L2Nu', c,
                     '/data/nawoods/ntuples/uwvvSingleZ_mc_24jan2017/results/TTTo2L2Nu*.root',
                     True, lumi) for c in channels
        }
    # ttByChan = {
    #     c : MCSample('TTJets', c,
    #                  '/data/nawoods/ntuples/uwvvSingleZ_mc_24jan2017/results/TTJets*.root',
    #                  True, lumi) for c in channels
    #     }

tt = SampleGroup('TTJets', 'z', ttByChan, True)

stack = SampleStack('stack', 'z', [dy,tt])


dataByChan = {}
for c in channels:
    samplesByEra = {}
    if test:
        samplesByEra['test'] = DataSample('test', c,
                                          '/data/nawoods/ntuples/uwvvSingleZ_data_test/*.root',
                                          )
    else:
        for era in 'BCDEFGH':
            samplesByEra['2016{}'.format(era)] = DataSample('2016{}'.format(era), c,
                                                            '/data/nawoods/ntuples/uwvvSingleZ_data_24jan2017/results/Run2016{}*.root'.format(era)
                                                            )


    dataByChan[c] = SampleGroup('data_{}'.format(c), c, samplesByEra)
    dataByChan[c].format(drawstyle='PE', legendstyle='LPE', color='black')


data = SampleGroup('Data', 'z', dataByChan)
data.format(drawstyle='PE', legendstyle='LPE', color='black')


binning2l = {
    'Mass' : [60, 60., 120.],
    'Pt' : [60, 0., 300.],
    'Eta' : [48,-6.,6.],
    'Phi' : [24, -3.15,3.15],
    'nJets' : [6,-0.5,5.5],
    }

binning1l = {
    'Pt' : [50, 0., 150.],
    'Eta' : [20, -2.5, 2.5],
    'Phi' : [24, -3.15, 3.15],
    'Iso' : [8, 0., .4],
    'PVDXY' : [20, -.1, .1],
    'PVDZ' : [20, -.2, .2],
    'SIP3D' : [20, 0., 5.],
    }

xTitles = {
    'Mass' : 'm_{{{obj}}} \\, (\\text{{GeV}})',
    'Eta' : '\\eta_{{{obj}}}',
    'Phi' : '\\phi_{{{obj}}}',
    'Pt' : '{obj} \\, p_{{T}} \\, (\\text{{GeV}})',
    'nJets' : 'N_{{\\text{{jets}}}} \\, (\\text{{single-}}{obj} \\, \\text{{events}})',
    'Iso' : 'R_{{\\text{{Iso}}}} \\, ({obj})',
    'PVDXY' : '#Delta_{{xy}} \\, ({obj}) \\, (\\text{{cm}})',
    'PVDZ' : '#Delta_{{z}} ({obj}) \\, (\\text{{cm}})',
    'nvtx' : 'N_{{vtx}} \\, (\\text{{single-}}{obj} \, \\text{{events}})',
    'SIP3D' : 'SIP_{{3D}} \\, ({obj})',
    }

objectNames = {
    'z' : 'Z',
    'ze' : 'Z \\rightarrow e^{+}e^{-}',
    'zm' : 'Z \\rightarrow \\mu^{+}\\mu^{-}',
    'l' : '\\ell',
    'e' : 'e',
    'm' : '\\mu',
    }

legParams = {v:{} for v in xTitles}
legParams['Eta'] = {'topmargin':.25,'leftmargin':.3,'rightmargin':.3}
legParams['Phi'] = {'topmargin':.25,'leftmargin':.3,'rightmargin':.3}



vars2l = {}
vars2l['ze'] = {v : {'ee':v} for v in binning2l}
vars2l['zm'] = {v : {'mm':v} for v in binning2l}

vars2l['z'] = {v:{} for v in binning2l}
for v in binning2l:
    vars2l['z'][v].update(vars2l['ze'][v].copy())
    vars2l['z'][v].update(vars2l['zm'][v].copy())

vars1l = {}
vars1l['e'] = {v : {'ee':['e1'+v, 'e2'+v]} for v in binning1l}
vars1l['e']['Iso'] = {'ee':['e1ZZIso','e2ZZIso']}
vars1l['m'] = {v : {'mm':['m1'+v, 'm2'+v]} for v in binning1l}
vars1l['m']['Iso'] = {'mm':['m1ZZIso','m2ZZIso']}

vars1l['l'] = {v:{} for v in binning1l}
for v in binning1l:
    vars1l['l'][v].update(vars1l['e'][v].copy())
    vars1l['l'][v].update(vars1l['m'][v].copy())


for chan in ['z', 'ze', 'zm']:
    for varName, var in vars2l[chan].iteritems():
        print "Plotting {} {}".format(chan, varName)

        hStack = stack.makeHist(var, '', binning2l[varName], mcWeight)
        dataPts = data.makeHist(var, '', binning2l[varName], poissonErrors=True)

        # for ratio
        dataHist = data.makeHist(var, '', binning2l[varName])

        print "    data: {} MC: {:.2f}".format(dataHist.Integral(), hStack.Integral())
        print ''

        c = Canvas(1000,1000)

        leg = makeLegend(c, hStack, dataPts, solid=True, **legParams[varName])

        pad1, pad2 = addPadBelow(c, .23)

        pad1.cd()
        (xaxis, yaxis), (xmin,xmax,ymin,ymax) = draw([hStack, dataPts], pad1,
                                                     xtitle=xTitles[varName].format(obj=objectNames[chan]),
                                                     ytitle='Events')
        leg.Draw("same")

        pad2.cd()
        ratio, unity = makeRatio(dataHist, hStack)
        (ratioX, ratioY), ratioLimits = draw(ratio, pad2, ytitle='Data / MC',
                                             xlimits=(xmin,xmax),
                                             ylimits=(0.7,1.3), ydivisions=5)
        unity.Draw("same")

        c.cd()
        pad1.Draw()
        pad2.Draw()

        fixRatioAxes(xaxis,yaxis,ratioX,ratioY, pad1.height, pad2.height)

        style.setCMSStyle(c, '', dataType='   Work in Progress', intLumi=lumi)
        c.Print('{}/{}{}.png'.format(outdir, chan, varName))

for chan in ['l', 'e', 'm']:
    for varName, var in vars1l[chan].iteritems():
        print "Plotting {} {}".format(chan, varName)

        hStack = stack.makeHist(var, '', binning1l[varName], mcWeight)
        dataPts = data.makeHist(var, '', binning1l[varName], poissonErrors=True)

        # for ratio
        dataHist = data.makeHist(var, '', binning1l[varName])

        print "    data: {} MC: {:.2f}".format(dataHist.Integral(), hStack.Integral())
        print ''

        c = Canvas(1000,1200)

        leg = makeLegend(c, hStack, dataPts, solid=True, **legParams[varName])

        pad1, pad2 = addPadBelow(c, .23)

        pad1.cd()
        (xaxis, yaxis), (xmin,xmax,ymin,ymax) = draw([hStack, dataPts], pad1,
                                                     xtitle=xTitles[varName].format(obj=objectNames[chan]),
                                                     ytitle='Leptons')
        leg.Draw("same")

        pad2.cd()
        ratio, unity = makeRatio(dataHist, hStack)
        (ratioX, ratioY), ratioLimits = draw(ratio, pad2, ytitle='Data / MC',
                                             xlimits=(xmin,xmax),
                                             ylimits=(0.7,1.3), ydivisions=5)
        unity.Draw("same")

        c.cd()
        pad1.Draw()
        pad2.Draw()

        fixRatioAxes(xaxis,yaxis,ratioX,ratioY, pad1.height, pad2.height)

        style.setCMSStyle(c, '', dataType='   Work in Progress', intLumi=lumi)
        c.Print('{}/{}{}.png'.format(outdir, chan, varName))

