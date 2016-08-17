

from rootpy.plotting import Canvas, Legend
from rootpy.plotting.utils import draw

from SampleTools import MCSample, DataSample, SampleGroup, SampleStack
from PlotTools import PlotStyle as _Style
from PlotTools import makeLegend


outdir = '/afs/cern.ch/user/n/nawoods/www/UWVVPlots/singleZ'

style = _Style()

lumi = 12900.

channels = ['ee','mm']

mcWeight = {
    'ee' : 'e1EffScaleFactor * e2EffScaleFactor',
    'mm' : 'm1EffScaleFactor * m2EffScaleFactor',
    }
dyByChan = {
    c : MCSample('DYJets', c, 
                 '/data/nawoods/ntuples/singleZ_mc_11aug2016/results/DYJets*.root', 
                 True, lumi) for c in channels
    }
dy = SampleGroup('DYJets', 'z', dyByChan, True)

ttByChan = {
    c : MCSample('TTJets', c, 
                 '/data/nawoods/ntuples/singleZ_mc_11aug2016/results/TTJets*.root', 
                 True, lumi) for c in channels
    }
tt = SampleGroup('TTJets', 'z', ttByChan, True)

stack = SampleStack('stack', 'z', [dy,tt])


dataByChan = {}
for c in channels:
    samplesByEra = {}
    for era in ['b','c','d']:
        samplesByEra['2016{}'.format(era)] = DataSample('data2016{}_{}'.format(era, c), c, 
                                                        '/data/nawoods/ntuples/singleZ_data2016{}_11aug2016/results/*.root'.format(era)
                                                        )

    dataByChan[c] = SampleGroup('data_{}'.format(c), c, samplesByEra)

data = SampleGroup('Data', 'z', dataByChan)


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

binning2l = {
    'Mass' : [60, 60., 120.],
    'Pt' : [60, 0., 300.],
    'Eta' : [40,-5.,5.],
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


for chan in ['z', 'ze','zm']:
    for varName, var in vars2l[chan].iteritems():
        print "Plotting {} {}".format(chan, varName)

        hStack = stack.makeHist(var, '', binning2l[varName], mcWeight)
        dataPts = data.makeHist(var, '', binning2l[varName], poissonErrors=True)

        c = Canvas(1000,1000)

        leg = makeLegend(c, hStack, dataPts)

        draw([hStack, dataPts], c, 
             xtitle='{}'.format(varName)+(' ({})'.format(units[varName]) if units[varName] else ''), 
             ytitle='Events')
        leg.Draw("same")

        style.setCMSStyle(c, '', dataType='Preliminary', intLumi=lumi)
        c.Print('{}/{}{}.png'.format(outdir, chan, varName))

for chan in ['l', 'e', 'm']:
    for varName, var in vars1l[chan].iteritems():
        print "Plotting {} {}".format(chan, varName)
        hStack = stack.makeHist(var, '', binning1l[varName], mcWeight)
        dataPts = data.makeHist(var, '', binning1l[varName], poissonErrors=True)

        c = Canvas(1000,1000)

        leg = makeLegend(c, hStack, dataPts)

        draw([hStack, dataPts], c, 
             xtitle='{}'.format(varName)+(' ({})'.format(units[varName]) if units[varName] else ''), 
             ytitle='Leptons')
        leg.Draw("same")

        style.setCMSStyle(c, '', dataType='Preliminary', intLumi=lumi)
        c.Print('{}/{}{}.png'.format(outdir, chan, varName))
                       
