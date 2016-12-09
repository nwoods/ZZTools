
import logging
from rootpy import log as rlog; rlog = rlog["/zzPlots"]
# don't show most silly ROOT messages
logging.basicConfig(level=logging.WARNING)
rlog["/ROOT.TUnixSystem.SetDisplay"].setLevel(rlog.ERROR)

from rootpy import asrootpy
from rootpy.io import root_open
from rootpy.plotting import Canvas, Legend
from rootpy.plotting.utils import draw
from rootpy.ROOT import TBox, Double

from SampleTools import MCSample, DataSample, SampleGroup, SampleStack
from PlotTools import PlotStyle as _Style
from PlotTools import makeLegend, addPadBelow, makeRatio, fixRatioAxes
from Utilities import WeightStringMaker, deltaRString, deltaPhiString, \
    makeNumberPretty
from Analysis import standardZZSamples

from os import environ
from os import path as _path
from os import makedirs as _mkdir
from collections import OrderedDict
from math import sqrt


inData = 'uwvvNtuples_data_25nov2016'
inMC = 'uwvvNtuples_mc_25nov2016'

puWeightFile = 'puWeight_69200_08sep2016'
fakeRateFile = 'fakeRate_08sep2016'

ana = 'full'

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
                                fakeRateFile, lumi, amcatnlo=amcatnlo,
                                higgs=(ana=='full'))

# # count events
# tot = OrderedDict()
# totErrSqr = {}
# for c in ['eeee','eemm','mmmm']:
#     print c + ':'
#
#     expected = 0.
#     expErrSqr = 0.
#     for s in stack:
#         h = s[c].makeHist('1.', '', [1,0.,2.])
#         yErr = Double(0)
#         y = h.IntegralAndError(0, h.GetNbinsX(), yErr)
#         if s.name == 'Z+X':
#             yErr = sqrt(yErr ** 2 + (.4 * y) ** 2)
#         print '    {}: {} +/- {}'.format(s.name, y, yErr)
#         if s.name not in tot:
#             tot[s.name] = 0.
#             totErrSqr[s.name] = 0.
#         tot[s.name] += y
#         totErrSqr[s.name] += yErr ** 2
#
#         expected += y
#         expErrSqr += yErr ** 2
#     print '    Total expected: {} +/- {}'.format(expected, sqrt(expErrSqr))
#     if 'expected' not in tot:
#         tot['expected'] = 0.
#         totErrSqr['expected'] = 0.
#     tot['expected'] += expected
#     totErrSqr['expected'] += expErrSqr
#
#     hData = data[c].makeHist('1.', '', [1,0.,2.])
#     dataErr = Double(0)
#     yieldData = hData.IntegralAndError(0, hData.GetNbinsX(), dataErr)
#     print '    Data: {} +/- {}'.format(yieldData, dataErr)
#     if 'data' not in tot:
#         tot['data'] = 0.
#         totErrSqr['data'] = 0.
#     tot['data'] += yieldData
#     totErrSqr['data'] += dataErr ** 2
#
#     print ''
#
# print 'Total:'
# for n,t in tot.iteritems():
#     print '    {}: {} +/- {}'.format(n,t,sqrt(totErrSqr[n]))


objNames = {
    'zz' : 'ZZ',
    'eeee' : '4e',
    'eemm' : '2e2\\mu',
    'mmmm' : '4\\mu',
    'z' : 'Z',
    'z1' : 'Z_{1}',
    'z2' : 'Z_{2}',
    'ze' : 'Z \\rightarrow e^{+}e^{-}',
    'zm' : 'Z \\rightarrow \\mu^{+}\\mu^{-}',
    'l1' : '\\ell_{1}',
    'l' : '\\ell',
    'e' : 'e',
    'm' : '\\mu',
    }

### Set up variable specific info

units = {
    'Pt' : 'GeV',
    'Eta' : '',
    'Phi' : '',
    'nJets' : '',
    'Mass' : 'GeV',
    'jet1Pt' : 'GeV',
    'jet1Eta' : '',
    'jet2Pt' : 'GeV',
    'jet2Eta' : '',
    'mjj' : 'GeV',
    'deltaEtajj' : '',
    'deltaPhiZZ' : '',
    'deltaRZZ' : '',
    'Iso' : '',
    'PVDXY' : 'cm',
    'PVDZ' : 'cm',
    'nvtx' : '',
    'SIP3D' : '',
    }

xTitles = {
    'Mass' : 'm_{{{obj}}} \\, (\\text{{GeV}})',
    'Eta' : '\\eta_{{{obj}}}',
    'Phi' : '\\phi_{{{obj}}}',
    'Pt' : '{obj} \\, p_{{T}} \\, (\\text{{GeV}})',
    'nJets' : 'N_{\\text{jets}}',
    'Iso' : 'R_{{Iso}} \\, ({obj})',
    'PVDXY' : '\\Delta_{{xy}} \\, ({obj}) \\, (\\text{{cm}})',
    'PVDZ' : '\\Delta_{{z}} \\, ({obj}) \\, (\\text{{cm}})',
    'nvtx' : 'N_{\\text{vtx}}',
    'SIP3D' : 'SIP_{{3D}} \\, ({obj})',
    'jet1Pt' : 'p_T^\\text{j1} \\, \\text{(GeV)}',
    'jet1Eta' : '\\eta_\\text{j1}',
    'jet2Pt' : 'p_T^\\text{j2} \\, \\text{(GeV)}',
    'jet2Eta' : '\\eta_\\text{j2}',
    'mjj' : 'm_\\text{jj} \\, \\text{(GeV)}',
    'deltaEtajj' : '|\\Delta \\eta_{\\text{jj}}}|',
    'deltaPhiZZ' : '\\Delta \\phi (\\text{Z}_1, \\text{Z}_2)',
    'deltaRZZ' : '\\Delta \\text{R} (\\text{Z}_1, \\text{Z}_2)',
    }

for v,t in xTitles.iteritems():
    if units[v]:
        t += ' \\, [\\text{{{{{}}}}}]'.format(units[v])

# some distributions need the legend moved
legParamsLeft = {
    'leftmargin' : 0.05,
    'rightmargin' : 0.5,
    }


binning4l = {
    'Mass'  : [100.] + [200.+50.*i for i in range(5)] + [500.,600.,800., 1000., 1200.],
    'Pt'    : [25.*i for i in range(4)] + [100., 150., 200., 300.],
    'Eta'   : [16, -5., 5.],
    'Phi'   : [12, -3.15, 3.15],
    'nvtx'  : [40, 0., 40.],
    'nJets' : [6, -0.5, 5.5],
    'jet1Pt' : [0., 50., 100., 200., 300., 500.],
    'jet1Eta' : [0., 1.5, 3., 4.7],
    'jet2Pt' : [30., 100., 200., 500.],
    'jet2Eta' : [0., 1.5, 3., 4.7],
    'mjj' : [0., 100., 300., 800.],
    'deltaEtajj' : [6, 0.,6.],
    'deltaPhiZZ' : [0., 1.5] + [2.+.25*i for i in range(6)],
    'deltaRZZ' : [6, 0., 6.],
    }

binNormWidth4l = {
    'Mass' : 50.,
    'Pt' : 25.,
    'Eta' : 1.,
    'Phi' : 1.,
    'nvtx' : 1.,
    'nJets' : 1.,
    'jet1Pt' : 50.,
    'jet2Pt' : 50.,
    'jet1Eta' : 1.,
    'jet2Eta' : 1.,
    'mjj' : 100.,
    'deltaEtajj' : 1.,
    'deltaPhiZZ' : 1.,
    'deltaRZZ' : 1.,
    }

vars4l = {v:v for v in binning4l}
vars4l['jet1Pt'] = 'jetPt[0]'
vars4l['jet2Pt'] = 'jetPt[1]'
vars4l['jet1Eta'] = 'abs(jetEta[0])'
vars4l['jet2Eta'] = 'abs(jetEta[1])'
vars4l = {v:{c:vars4l[v] for c in ['eeee','eemm','mmmm']} for v in vars4l}
vars4l['deltaPhiZZ'] = {
    'eeee' : '{}(e1_e2_Phi, e3_e4_Phi)'.format(deltaPhiString()),
    'eemm' : '{}(e1_e2_Phi, m1_m2_Phi)'.format(deltaPhiString()),
    'mmmm' : '{}(m1_m2_Phi, m3_m4_Phi)'.format(deltaPhiString()),
    }
vars4l['deltaRZZ'] = {
    'eeee' : '{}(e1_e2_Eta, e1_e2_Phi, e3_e4_Eta, e3_e4_Phi)'.format(deltaRString()),
    'eemm' : '{}(e1_e2_Eta, e1_e2_Phi, m1_m2_Eta, m1_m2_Phi)'.format(deltaRString()),
    'mmmm' : '{}(m1_m2_Eta, m1_m2_Phi, m3_m4_Eta, m3_m4_Phi)'.format(deltaRString()),
    }

selections4l = {v:'' for v in vars4l}
selections4l['jet1Pt'] = 'nJets >= 1'
selections4l['jet2Pt'] = 'nJets >= 2'
selections4l['jet1Eta'] = 'nJets >= 1'
selections4l['jet2Eta'] = 'nJets >= 2'

if amcatnlo:
    binning4l = {'Mass' : binning4l['Mass'], 'Pt' : binning4l['Pt']}

if ana == 'z4l':
    binning4l['Mass'] = [20, 80., 100.]
    binNormWidth4l['Mass'] = 1.
elif ana == 'full':
    binning4l['Mass'] = [25.*i for i in range(17)] + [500.,600.,800.]
    binNormWidth4l['Mass'] = 25.

for chan in ['zz', 'eeee', 'eemm', 'mmmm']:
    for varName, binning in binning4l.iteritems():
        print "Plotting {} {}".format(chan, varName)

        var = vars4l[varName]
        if chan != 'zz':
            var = {chan:var[chan]}

        # blinding
        dataSelection = ''
        if varName == 'Mass':
            dataSelection = 'Mass < 500.'
            if selections4l[varName]:
                dataSelection += ' && ' + selections4l[varName]

        hStack = stack.makeHist(var, selections4l[varName], binning,
                                postprocess=True,
                                perUnitWidth=binNormWidth4l[varName])
        dataPts = data.makeHist(var, dataSelection, binning,
                                poissonErrors=True,
                                perUnitWidth=binNormWidth4l[varName])
        toPlot = [hStack, dataPts]

        c = Canvas(1000,1000)

        legParams = {}
        if ana == 'z4l' and varName == 'Mass' or ana == 'smp' and varName == 'deltaRZZ':
            legParams = legParamsLeft.copy()
        leg = makeLegend(c, *toPlot, **legParams)

        xTitle = xTitles[varName]
        if 'obj' in xTitle:
            xTitle = xTitle.format(obj=objNames[chan])

        yTitle = 'Events / {} {}'.format(makeNumberPretty(binNormWidth4l[varName], 2),
                                         units[varName])

        (xaxis, yaxis), (xmin,xmax,ymin,ymax) = draw(toPlot, c,
                                                     xtitle=xTitle,
                                                     ytitle=yTitle)
        # blinding box
        if varName == 'Mass' and binning4l['Mass'][-1] > 500.:
            box = TBox(max(xmin,500.), ymin, min(binning4l['Mass'][-1], xmax), ymax)
            box.SetFillColor(1)
            box.SetFillStyle(3002)
            box.Draw("same")
            leg.SetFillStyle(1001)

        leg.Draw("same")

        style.setCMSStyle(c, '', dataType='Preliminary', intLumi=lumi)
        c.Print('{}/{}{}.png'.format(outdir, chan, varName))


        if varName == 'Mass' and ana == 'smp' and chan == 'zz':
            aTGCFileTemp = '/data/nawoods/aTGCSherpaHistos/histo_{fg}I{fz}_{param}_file.root'

            try:
                with root_open(aTGCFileTemp.format(fg='0',fz='0',param='f4')) as f:
                    hSherpaSM = asrootpy(f.h_ratio_MZZ_wt)
                    hSherpaSM.SetDirectory(0)
                with root_open(aTGCFileTemp.format(fg='0p0038',fz='0p003',param='f5')) as f:
                    hSherpa_f5_maxTGC = asrootpy(f.h_ratio_MZZ_wt)
                    hSherpa_f5_maxTGC.SetDirectory(0)
                with root_open(aTGCFileTemp.format(fg='0p0038',fz='0p003',param='f4')) as f:
                    hSherpa_f4_maxTGC = asrootpy(f.h_ratio_MZZ_wt)
                    hSherpa_f4_maxTGC.SetDirectory(0)

                ratio_f5_maxTGC = hSherpa_f5_maxTGC / hSherpaSM
                ratio_f5_maxTGC[-1].value = ratio_f5_maxTGC[-2].value
                ratio_f4_maxTGC = hSherpa_f4_maxTGC / hSherpaSM
                ratio_f4_maxTGC[-1].value = ratio_f4_maxTGC[-2].value
                wtMaker = WeightStringMaker('aTGC')
                wt_f5_maxTGC = wtMaker.makeWeightStringFromHist(ratio_f5_maxTGC, 'Mass')
                wt_f4_maxTGC = wtMaker.makeWeightStringFromHist(ratio_f4_maxTGC, 'Mass')

                h_f5_maxTGC = hStack[0].empty_clone()
                h_f4_maxTGC = hStack[0].empty_clone()
                for s in stack:
                    if 'ZZTo4L' in s.name:
                        h_f4_maxTGC += s.makeHist(var, selections4l[varName],
                                                  binning,
                                                  wt_f4_maxTGC,
                                                  postprocess=True,
                                                  perUnitWidth=binNormWidth4l[varName],
                                                  mergeOverflow=True)
                        h_f5_maxTGC += s.makeHist(var, selections4l[varName],
                                                  binning,
                                                  wt_f5_maxTGC,
                                                  postprocess=True,
                                                  perUnitWidth=binNormWidth4l[varName],
                                                  mergeOverflow=True)
                    else:
                        hTemp = s.makeHist(var, selections4l[varName],
                                           binning,
                                           postprocess=True,
                                           perUnitWidth=binNormWidth4l[varName],
                                           mergeOverflow=True)
                        h_f4_maxTGC += hTemp
                        h_f5_maxTGC += hTemp

                h_f5_maxTGC.color = 'r'
                h_f5_maxTGC.drawstyle = 'hist'
                h_f5_maxTGC.fillstyle = 'hollow'
                h_f5_maxTGC.linestyle = 'dashed'
                h_f5_maxTGC.SetLineWidth(h_f5_maxTGC.GetLineWidth() * 2)
                h_f5_maxTGC.legendstyle = 'L'
                h_f5_maxTGC.title = 'f_{5}^{#gamma} = 0.0038, f_{5}^{Z} = 0.003'

                h_f4_maxTGC.color = 'magenta'
                h_f4_maxTGC.drawstyle = 'hist'
                h_f4_maxTGC.fillstyle = 'hollow'
                h_f4_maxTGC.linestyle = 'dashed'
                h_f4_maxTGC.SetLineWidth(h_f4_maxTGC.GetLineWidth() * 2)
                h_f4_maxTGC.legendstyle = 'L'
                h_f4_maxTGC.title = 'f_{4}^{#gamma} = 0.0038, f_{4}^{Z} = 0.003'

                toPlot = [hStack, h_f4_maxTGC, h_f5_maxTGC, dataPts]

                legParams['entryheight'] = 0.02
                legParams['entrysep'] = 0.008
                legParams['textsize'] = 0.022
                legParams['leftmargin'] = 0.4

                cTGC = Canvas(1000,1000)
                cTGC.SetLogy()
                leg = makeLegend(cTGC, *toPlot, **legParams)

                (xaxis, yaxis), (xmin,xmax,ymin,ymax) = draw(toPlot, cTGC,
                                                             xtitle=xTitle,
                                                             ytitle=yTitle,
                                                             logy=True)
                # blinding box
                if binning4l['Mass'][-1] > 500.:
                    box = TBox(max(xmin,500.), ymin, min(binning4l['Mass'][-1], xmax), ymax)
                    box.SetFillColor(1)
                    box.SetFillStyle(3002)
                    box.Draw("same")
                    leg.SetFillStyle(1001)

                leg.Draw("same")

                style.setCMSStyle(cTGC, '', dataType='Preliminary', intLumi=lumi)
                cTGC.Print('{}/{}{}_aTGC.png'.format(outdir, chan, varName))


            except Exception as e:
                print "Adding aTGC signal failed with the following exception:"
                print e
                print "If you want that on your plots, you should fix this problem!"



binning2l = {
    'Mass' : [60, 60., 120.],
    'Pt' : [i * 25. for i in range(7)] + [200., 300.],
    'Eta' : [48,-6.,6.],
    'Phi' : [24, -3.15,3.15],
    }

binNormWidth2l = {
    'Mass' : 1.,
    'Pt' : 25.,
    'Eta' : 0.25,
    'Phi' : 0.25,
    }

if ana != 'smp':
    binning2l['Mass'] = [60, 0., 120.]
    binNormWidth2l['Mass'] = 2.

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

        hStack = stack.makeHist(var, selections2l[z], binning2l[varName],
                                postprocess=True,
                                perUnitWidth=binNormWidth2l[varName])
        dataPts = data.makeHist(var, selections2l[z], binning2l[varName],
                                poissonErrors=True,
                                perUnitWidth=binNormWidth2l[varName])

        if varName == 'Pt' and binning2l['Pt'][-1] > 200.:
            copy = dataPts.clone()
            x = copy.GetX()

            idx = 0
            for i in range(copy.GetN()):
                if x[i] >= 200.:
                    dataPts.RemovePoint(idx)
                    continue
                idx += 1


        # for ratio
        dataHist = data.makeHist(var, '', binning2l[varName],
                                 perUnitWidth=binNormWidth2l[varName])

        c = Canvas(1000,1000)

        legParams = {}
        if ana == 'full' and varName == 'Mass':
            legParams = legParamsLeft.copy()
        leg = makeLegend(c, hStack, dataPts, **legParams)

        xTitle = xTitles[varName]
        if '{obj}' in xTitle:
            xTitle = xTitle.format(obj=objNames[z])

        yTitle = 'Z bosons / {} {}'.format(binNormWidth2l[varName],
                                           units[varName])

        (xaxis, yaxis), (xmin,xmax,ymin,ymax) = draw([hStack, dataPts], c,
                                                     xtitle=xTitle,
                                                     ytitle=yTitle)

        # blinding box
        if varName == 'Pt' and binning2l['Pt'][-1] > 200.:
            box = TBox(max(xmin,200.), ymin, min(binning2l['Pt'][-1], xmax), ymax)
            box.SetFillColor(1)
            box.SetFillStyle(3002)
            box.Draw("same")
            leg.SetFillStyle(1001)

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

binNormWidth1l = {
    'Pt' : 10.,
    'Eta' : 0.25,
    'Phi' : 0.25,
    'Iso' : 0.05,
    'PVDXY' : 0.01,
    'PVDZ' : 0.02,
    'SIP3D' : 0.25,
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

        hStack = stack.makeHist(var, selections1l[lep], binning1l[varName],
                                postprocess=True,
                                perUnitWidth=binNormWidth1l[varName])
        dataPts = data.makeHist(var, selections1l[lep], binning1l[varName],
                                poissonErrors=True,
                                perUnitWidth=binNormWidth1l[varName])

        # for ratio
        dataHist = data.makeHist(var, selections1l[lep], binning1l[varName],
                                 poissonErrors=True,
                                 perUnitWidth=binNormWidth1l[varName])

        c = Canvas(1000,1000)

        leg = makeLegend(c, hStack, dataPts)

        xTitle = xTitles[varName]
        if '{obj}' in xTitle:
            xTitle = xTitle.format(obj=objNames[lep])

        yTitle = 'Leptons / {} {}'.format(binNormWidth1l[varName],
                                          units[varName])

        (xaxis, yaxis), (xmin,xmax,ymin,ymax) = draw([hStack, dataPts], c,
                                                     xtitle=xTitle,
                                                     ytitle=yTitle)
        leg.Draw("same")

        style.setCMSStyle(c, '', dataType='Preliminary', intLumi=lumi)
        c.Print('{}/{}{}.png'.format(outdir, lep, varName))
