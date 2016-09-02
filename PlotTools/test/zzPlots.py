
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

from os import environ
from os import path as _path



test = False
inData = '/data/nawoods/ntuples/uwvvNtuples_data_23aug2016'
sampleIDMC = '/data/nawoods/ntuples/uwvvNtuples_data_23aug2016'

puWeightFile = 'puWeight_69200_23aug2016.root'

outdir = '/afs/cern.ch/user/n/nawoods/www/UWVVPlots/zz'
if test:
    outdir = '/afs/cern.ch/user/n/nawoods/www/UWVVPlots/test'

style = _Style()

lumi = 12900.

channels = ['eeee','eemm', 'mmmm']

### Set up base samples

if test:
    qqZZByChan = {
        c : MCSample('ZZTo4L', c, 
                     '/data/nawoods/ntuples/uwvvNtuples_mc_test/ZZTo4L*.root', 
                     True, lumi) for c in channels
        }
    qqZZByChan2P2F = {
        c : MCSample('ZZTo4L', c, 
                     '/data/nawoods/ntuples/uwvvNtuples_mc_test/2P2F/ZZTo4L*.root', 
                     True, lumi) for c in channels
        }
    qqZZByChan3P1F = {
        c : MCSample('ZZTo4L', c, 
                     '/data/nawoods/ntuples/uwvvNtuples_mc_test/3P1F/ZZTo4L*.root', 
                     True, lumi) for c in channels
        }
else:
    qqZZByChan = {
        c : MCSample('ZZTo4L', c, 
                     '/data/nawoods/ntuples/uwvvNtuples_mc_23aug2016/results_smp/ZZTo4L*.root', 
                     True, lumi) for c in channels
        }
    qqZZByChan2P2F = {
        c : MCSample('ZZTo4L', c, 
                     '/data/nawoods/ntuples/uwvvNtuples_mc_23aug2016/results_smp_2P2F/ZZTo4L*.root', 
                     True, lumi) for c in channels
        }
    qqZZByChan3P1F = {
        c : MCSample('ZZTo4L', c, 
                     '/data/nawoods/ntuples/uwvvNtuples_mc_23aug2016/results_smp_3P1F/ZZTo4L*.root', 
                     True, lumi) for c in channels
        }

ggZZByChan = {}
ggZZByChan2P2F = {}
ggZZByChan3P1F = {}
if test:
    for c in channels:
        ggZZByFS = {
            fs : MCSample('GluGluZZTo{}'.format(fs), c, 
                          '/data/nawoods/ntuples/uwvvNtuples_mc_test/GluGluZZTo{}*.root'.format(fs),
                          True, lumi) for fs in ['4e', '4mu', '2e2mu']
            }
        ggZZByChan[c] = SampleGroup('GluGluZZ', c, ggZZByFS, True)
        ggZZByFS2P2F = {
            fs : MCSample('GluGluZZTo{}'.format(fs), c, 
                          '/data/nawoods/ntuples/uwvvNtuples_mc_test/2P2F/GluGluZZTo{}*.root'.format(fs),
                          True, lumi) for fs in ['4e', '4mu', '2e2mu']
            }
        ggZZByChan3P1F[c] = SampleGroup('GluGluZZ', c, ggZZByFS3P1F, True)
        ggZZByFS3P1F = {
            fs : MCSample('GluGluZZTo{}'.format(fs), c, 
                          '/data/nawoods/ntuples/uwvvNtuples_mc_test/3P1F/GluGluZZTo{}*.root'.format(fs),
                          True, lumi) for fs in ['4e', '4mu', '2e2mu']
            }
        ggZZByChan3P1F[c] = SampleGroup('GluGluZZ', c, ggZZByFS3P1F, True)
else:
    for c in channels:
        ggZZByFS = {
            fs : MCSample('GluGluZZTo{}'.format(fs), c, 
                          '/data/nawoods/ntuples/uwvvNtuples_mc_23aug2016/results_smp/GluGluZZTo{}*.root'.format(fs),
                          True, lumi) for fs in ['4e', '4mu', '2e2mu']
            }
        ggZZByChan[c] = SampleGroup('GluGluZZ', c, ggZZByFS, True)
        ggZZByFS2P2F = {
            fs : MCSample('GluGluZZTo{}'.format(fs), c, 
                          '/data/nawoods/ntuples/uwvvNtuples_mc_23aug2016/results_smp_2P2F/GluGluZZTo{}*.root'.format(fs),
                          True, lumi) for fs in ['4e', '4mu', '2e2mu']
            }
        ggZZByChan2P2F[c] = SampleGroup('GluGluZZ', c, ggZZByFS2P2F, True)
        ggZZByFS3P1F = {
            fs : MCSample('GluGluZZTo{}'.format(fs), c, 
                          '/data/nawoods/ntuples/uwvvNtuples_mc_23aug2016/results_smp_3P1F/GluGluZZTo{}*.root'.format(fs),
                          True, lumi) for fs in ['4e', '4mu', '2e2mu']
            }
        ggZZByChan3P1F[c] = SampleGroup('GluGluZZ', c, ggZZByFS3P1F, True)


dataByChan = {}
dataByChan2P2F = {}
dataByChan3P1F = {}
for c in channels:
    samplesByEra = {}
    samplesByEra2P2F = {}
    samplesByEra3P1F = {}
    if test:
        samplesByEra['test'] = DataSample('test', c, 
                                          '/data/nawoods/ntuples/uwvvNtuples_data_test/*.root',
                                          )
        samplesByEra2P2F['test'] = DataSample('test', c, 
                                              '/data/nawoods/ntuples/uwvvNtuples_data_test/2P2F/*.root',
                                          )
        samplesByEra3P1F['test'] = DataSample('test', c, 
                                          '/data/nawoods/ntuples/uwvvNtuples_data_test/3P1F/*.root',
                                          )
    else:
        for era in ['B','C','D']:
            samplesByEra['2016{}'.format(era)] = DataSample('data2016{}_{}'.format(era, c), c, 
                                                            '/data/nawoods/ntuples/uwvvNtuples_data_23aug2016/results_smp/Run2016{}*.root'.format(era)
                                                            )
            samplesByEra2P2F['2016{}'.format(era)] = DataSample('data2016{}_{}'.format(era, c), c, 
                                                                '/data/nawoods/ntuples/uwvvNtuples_data_23aug2016/'
                                                                'results_smp_2P2F/Run2016{}*.root'.format(era)
                                                                )
            samplesByEra3P1F['2016{}'.format(era)] = DataSample('data2016{}_{}'.format(era, c), c, 
                                                                '/data/nawoods/ntuples/uwvvNtuples_data_23aug2016/'
                                                                'results_smp_3P1F/Run2016{}*.root'.format(era)
                                                                )

    dataByChan[c] = SampleGroup('data_{}'.format(c), c, samplesByEra)
    dataByChan2P2F[c] = SampleGroup('data_{}'.format(c), c, samplesByEra2P2F)
    dataByChan3P1F[c] = SampleGroup('data_{}'.format(c), c, samplesByEra3P1F)


### Apply weights to samples

puWeight = WeightStringMaker('puWeight')
with root_open(_path.join(environ['zzt'], 'data', 'pileup', 
                          puWeightFile)) as f:
    strPU = puWeight.makeWeightStringFromHist(f.puScaleFactor, 'nTruePU')

mcWeight = {
    'eeee' : 'e1EffScaleFactor * e2EffScaleFactor * e3EffScaleFactor * e4EffScaleFactor * {}'.format(strPU),
    'eemm' : 'e1EffScaleFactor * e2EffScaleFactor * m1EffScaleFactor * m2EffScaleFactor * {}'.format(strPU),
    'mmmm' : 'm1EffScaleFactor * m2EffScaleFactor * m3EffScaleFactor * m4EffScaleFactor * {}'.format(strPU),
    }


wCR = WeightStringMaker('fakeFactor')
with root_open(_path.join(environ['zzt'], 'data', 'fakeRate',
                          'fakeRate_24aug2016.root')) as f:
    fakeFactorStrE = wCR.makeWeightStringFromHist(f.fakeFactor_e, 'abs({lep}Eta)', '{lep}Pt')
    fakeFactorStrM = wCR.makeWeightStringFromHist(f.fakeFactor_m, 'abs({lep}Eta)', '{lep}Pt')
zCRWeightTemp = ('({{lep1}}ZZTightID && {{lep1}}ZZIsoPass ? 1. : {fr1}) * '
                 '({{lep2}}ZZTightID && {{lep2}}ZZIsoPass ? 1. : {fr2})')
zeCRWeight = zCRWeightTemp.format(fr1=fakeFactorStrE.format(lep='{lep1}'),
                                  fr2=fakeFactorStrE.format(lep='{lep2}'))
zmCRWeight = zCRWeightTemp.format(fr1=fakeFactorStrM.format(lep='{lep1}'),
                                  fr2=fakeFactorStrM.format(lep='{lep2}'))
crWeight = {
    'eeee' : zeCRWeight.format(lep1='e1',lep2='e2') + ' * ' + zeCRWeight.format(lep1='e3',lep2='e4'), 
    'eemm' : zeCRWeight.format(lep1='e1',lep2='e2') + ' * ' + zmCRWeight.format(lep1='m1',lep2='m2'), 
    'mmmm' : zmCRWeight.format(lep1='m1',lep2='m2') + ' * ' + zmCRWeight.format(lep1='m3',lep2='m4'), 
    }

def ensureNonneg(h):
    '''
    Remove any bin of h with value < 0. Any nonnegative bin whose error bar
    would go below 0 has the error set to the bin value.
    '''
    for b in h:
        if b.value < 0.:
            b.value = 0.
            b.error = 0.
        elif b.error > b.value:
            b.error = b.value

    return h

for c in channels:
    # all MC samples need standard MC weight
    qqZZByChan[c].applyWeight(mcWeight[c])
    qqZZByChan2P2F[c].applyWeight(mcWeight[c])
    qqZZByChan3P1F[c].applyWeight(mcWeight[c])
    ggZZByChan[c].applyWeight(mcWeight[c])
    ggZZByChan2P2F[c].applyWeight(mcWeight[c])
    ggZZByChan3P1F[c].applyWeight(mcWeight[c])

    # CR samples weighted by fake factor
    # 2P2F is subtracted from 3P1F so weight it by -1
    # ... but MC CRs are subtracted from data CRs, so give them all opposite sign
    dataByChan2P2F[c].applyWeight(crWeight[c])
    dataByChan2P2F[c].applyWeight('-1.')
    dataByChan3P1F[c].applyWeight(crWeight[c])

    qqZZByChan2P2F[c].applyWeight(crWeight[c])
    qqZZByChan3P1F[c].applyWeight(crWeight[c])
    qqZZByChan3P1F[c].applyWeight('-1.')
    ggZZByChan2P2F[c].applyWeight(crWeight[c])
    ggZZByChan3P1F[c].applyWeight(crWeight[c])
    ggZZByChan3P1F[c].applyWeight('-1.')

### Make the stack and data points we'll actually use

data = SampleGroup('Data', 'zz', dataByChan)
data.format(color='black',drawstyle='PE',legendstyle='LPE')

qqZZ = SampleGroup('ZZTo4L', 'zz', qqZZByChan, True)
ggZZ = SampleGroup('GluGluZZ', 'zz', ggZZByChan, True)
bkgByChan = {
    c : SampleGroup('Z+X', c, {'2P2F':dataByChan2P2F[c], 
                               '3P1F':dataByChan3P1F[c], 
                               'qqZZ2P2F':qqZZByChan2P2F[c], 
                               'qqZZ3P1F':qqZZByChan3P1F[c], 
                               'ggZZ2P2F':ggZZByChan2P2F[c], 
                               'ggZZ3P1F':ggZZByChan3P1F[c]
                               }
                    ) for c in channels
    }

bkg = SampleGroup('Z+X', 'zz', bkgByChan, True)

# don't allow negative backgrounds
bkg.setPostprocessor(ensureNonneg)

stack = SampleStack('stack', 'zz', [bkg, ggZZ, qqZZ])


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
    'Mass'  : [35, 250., 2000.],
    'Pt'    : [40, 0., 200.],
    'Eta'   : [16, -5., 5.],
    'Phi'   : [12, -3.15, 3.15],
    'nvtx'  : [40, 0., 40.],
    'nJets' : [6, -0.5, 5.5],
    }

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
        if varName == 'Mass':
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

                       
