
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
from Analysis import standardZZSamples, standardZZMC

from os import environ
from os import path as _path
from os import makedirs as _mkdir
from os.path import isdir as _isdir
from os.path import exists as _exists
from collections import OrderedDict
from math import sqrt



_objNames = {
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

_units = {
    'Pt' : 'GeV',
    'Eta' : '',
    'Phi' : '',
    'nJets' : '',
    'nJets_eta2p4' : '',
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

_xTitles = {
    'Mass' : 'm_{{{obj}}} \\, (\\text{{GeV}})',
    'Eta' : '\\eta_{{{obj}}}',
    'Phi' : '\\phi_{{{obj}}}',
    'Pt' : '{obj} \\, p_{{T}} \\, (\\text{{GeV}})',
    'nJets' : 'N_{\\text{jets}}',
    'nJets_eta2p4' : 'N_{\\text{jets}} \\left( \\left|\\eta\\right| < 2.4 \\right)',
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

for v,t in _xTitles.iteritems():
    if _units[v]:
        t += ' \\, [\\text{{{{{}}}}}]'.format(_units[v])

# some distributions need the legend moved
_legParamsLeft = {
    'leftmargin' : 0.05,
    'rightmargin' : 0.53,
    'textsize' : 0.029,
    }


_binning4l = {
    'Mass'  : [100.] + [200.+50.*i for i in range(5)] + [500.,600.,800., 1000., 1200.],
    'Pt'    : [0.,5.]+[25.+25.*i for i in range(3)] + [100., 150., 200., 300.], #[25.*i for i in range(4)] + [100., 150., 200., 300.],
    'Eta'   : [16, -5., 5.],
    'Phi'   : [12, -3.15, 3.15],
    'nvtx'  : [40, 0., 40.],
    'nJets' : [6, -0.5, 5.5],
    'nJets_eta2p4' : [6, -0.5, 5.5],
    'jet1Pt' : [0., 50., 100., 200., 300., 500.],
    'jet1Eta' : [0., 1.5, 3., 4.7],
    'jet2Pt' : [30., 100., 200., 500.],
    'jet2Eta' : [0., 1.5, 3., 4.7],
    'mjj' : [0., 100., 300., 800.],
    'deltaEtajj' : [6, 0.,6.],
    'deltaPhiZZ' : [0., 1.5] + [2.+.25*i for i in range(6)],
    'deltaRZZ' : [6, 0., 6.],
    }

_binNormWidth4l = {
    'Mass' : 50.,
    'Pt' : 25.,
    'Eta' : 1.,
    'Phi' : 1.,
    'nvtx' : 1.,
    'nJets' : False,
    'nJets_eta2p4' : False,
    'jet1Pt' : 50.,
    'jet2Pt' : 50.,
    'jet1Eta' : 1.,
    'jet2Eta' : 1.,
    'mjj' : 100.,
    'deltaEtajj' : 1.,
    'deltaPhiZZ' : 1.,
    'deltaRZZ' : 1.,
    }

_vars4l = {v:v for v in _binning4l}
_vars4l['jet1Pt'] = 'jetPt[0]'
_vars4l['jet2Pt'] = 'jetPt[1]'
_vars4l['jet1Eta'] = 'abs(jetEta[0])'
_vars4l['jet2Eta'] = 'abs(jetEta[1])'
_vars4l = {v:{c:_vars4l[v] for c in ['eeee','eemm','mmmm']} for v in _vars4l}
_vars4l['deltaPhiZZ'] = {
    'eeee' : '{}(e1_e2_Phi, e3_e4_Phi)'.format(deltaPhiString()),
    'eemm' : '{}(e1_e2_Phi, m1_m2_Phi)'.format(deltaPhiString()),
    'mmmm' : '{}(m1_m2_Phi, m3_m4_Phi)'.format(deltaPhiString()),
    }
_vars4l['deltaRZZ'] = {
    'eeee' : '{}(e1_e2_Eta, e1_e2_Phi, e3_e4_Eta, e3_e4_Phi)'.format(deltaRString()),
    'eemm' : '{}(e1_e2_Eta, e1_e2_Phi, m1_m2_Eta, m1_m2_Phi)'.format(deltaRString()),
    'mmmm' : '{}(m1_m2_Eta, m1_m2_Phi, m3_m4_Eta, m3_m4_Phi)'.format(deltaRString()),
    }

_selections4l = {v:'' for v in _vars4l}
_selections4l['jet1Pt'] = 'nJets >= 1'
_selections4l['jet2Pt'] = 'nJets >= 2'
_selections4l['jet1Eta'] = 'nJets >= 1'
_selections4l['jet2Eta'] = 'nJets >= 2'

_binning2l = {
    'Mass' : [60, 60., 120.],
    'Pt' : [i * 25. for i in range(7)] + [200., 300.],
    'Eta' : [48,-6.,6.],
    'Phi' : [24, -3.15,3.15],
    }

_binNormWidth2l = {
    'Mass' : 1.,
    'Pt' : 25.,
    'Eta' : 0.25,
    'Phi' : 0.25,
    }

ze1VarTemp = 'e1_e2_{var}'
ze2VarTemp = 'e3_e4_{var}'
zm1VarTemp = 'm1_m2_{var}'
zm2VarTemp = 'm3_m4_{var}'
_varTemplates2l = {
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

_selections2l = {z:'' for z in _varTemplates2l}
_selections2l['z1'] = {
    'eeee' : '',
    'mmmm' : '',
    'eemm' : ['abs(e1_e2_Mass - 91.1876) < abs(m1_m2_Mass - 91.1876)',
              'abs(e1_e2_Mass - 91.1876) > abs(m1_m2_Mass - 91.1876)']
    }
_selections2l['z2'] = {
    'eeee' : '',
    'mmmm' : '',
    'eemm' : ['abs(e1_e2_Mass - 91.1876) > abs(m1_m2_Mass - 91.1876)',
              'abs(e1_e2_Mass - 91.1876) < abs(m1_m2_Mass - 91.1876)']
    }

_binning1l = {
    'Pt' : [20, 0., 200.],
    'Eta' : [20, -2.5, 2.5],
    'Phi' : [24, -3.15, 3.15],
    'Iso' : [8, 0., .4],
    'PVDXY' : [20, -.1, .1],
    'PVDZ' : [20, -.2, .2],
    'SIP3D' : [20, 0., 5.],
    }

_binNormWidth1l = {
    'Pt' : 10.,
    'Eta' : 0.25,
    'Phi' : 0.25,
    'Iso' : 0.05,
    'PVDXY' : 0.01,
    'PVDZ' : 0.02,
    'SIP3D' : 0.25,
    }

ze1LepVarTemp = ['e1{var}','e2{var}']
ze2LepVarTemp = ['e3{var}','e4{var}']
zm1LepVarTemp = ['m1{var}','m2{var}']
zm2LepVarTemp = ['m3{var}','m4{var}']
_varTemplates1l = {
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

_selections1l = {l:'' for l in _varTemplates1l}
_selections1l['l1'] = {
    'eeee' : ['e1Pt > e3Pt', 'e3Pt > e1Pt'],
    'eemm' : ['e1Pt > m1Pt', 'm1Pt > e1Pt'],
    'mmmm' : ['m1Pt > m3Pt', 'm3Pt > m1Pt'],
    }


def main(inData, inMC, plotDir, ana, fakeRateFile, puWeightFile, lumi,
         eras='BCDEFGH', blind=False, amcatnlo=False, leptonSFFromHists=False):

    style = _Style()

    outdir = _path.join('/afs/cern.ch/user/n/nawoods/www/UWVVPlots', plotDir)
    if not _exists(outdir):
        _mkdir(outdir)
    elif not _isdir(outdir):
        raise IOError("There is already some non-directory object called {}.".format(outdir))

    if amcatnlo:
        outdir += '_amcatnlo'

    data, stack = standardZZSamples('zz', inData, inMC, ana, puWeightFile,
                                    fakeRateFile, lumi, amcatnlo=amcatnlo,
                                    higgs=(ana=='full'), eras=eras,
                                    scaleFactorsFromHists=leptonSFFromHists)

    if ana == 'smp':
        aTGCf4 = standardZZMC('zz', 'uwvvNtuples_mc_21feb2017_aTGC',
                              'ZZTo4L-aTGC-f4-fg0p0038-fz0p003',
                              ana, puWeightFile, lumi)
        aTGCf5 = standardZZMC('zz', 'uwvvNtuples_mc_21feb2017_aTGC',
                              'ZZTo4L-aTGC-f5-fg0p0038-fz0p003',
                              ana, puWeightFile, lumi)
        sherpa = standardZZMC('zz', 'uwvvNtuples_mc_21feb2017_aTGC',
                              'ZZTo4L-sherpa',
                              ana, puWeightFile, lumi)

    binning4l = _binning4l.copy()
    binNormWidth4l = _binNormWidth4l.copy()

    if ana == 'z4l':
        binning4l['Mass'] = [20, 80., 100.]
        binNormWidth4l['Mass'] = 1.
    elif ana == 'full':
        binning4l['Mass'] = [25.*i for i in range(17)] + [500.,600.,800.] #[80.,100.,120.,130.,150.,180.,200.,240.,300.,400.,1000]
        binNormWidth4l['Mass'] = 25. #10.

    for chan in ['zz', 'eeee', 'eemm', 'mmmm']:
        for varName, binning in binning4l.iteritems():
            print "Plotting {} {}".format(chan, varName)

            var = _vars4l[varName]
            if chan != 'zz':
                var = {chan:var[chan]}

            # blinding
            dataSelection = ''
            if blind and varName == 'Mass':
                dataSelection = 'Mass < 500.'
                if _selections4l[varName]:
                    dataSelection += ' && ' + _selections4l[varName]

            hStack = stack.makeHist(var, _selections4l[varName], binning,
                                    postprocess=True,
                                    perUnitWidth=binNormWidth4l[varName])
            dataPts = data.makeHist(var, dataSelection, binning,
                                    poissonErrors=True,
                                    perUnitWidth=binNormWidth4l[varName])
            toPlot = [hStack, dataPts]

            c = Canvas(1000,1000)

            legParams = {'textsize':0.029}
            if ana == 'z4l' and varName == 'Mass' or ana == 'smp' and varName == 'deltaRZZ':
                legParams = _legParamsLeft.copy()
            leg = makeLegend(c, *toPlot, **legParams)

            xTitle = _xTitles[varName]
            if 'obj' in xTitle:
                xTitle = xTitle.format(obj=_objNames[chan])

            yTitle = 'Events'
            if binNormWidth4l[varName]:
                yTitle += ' / {} {}'.format(makeNumberPretty(binNormWidth4l[varName], 2),
                                            _units[varName])

            # cure inexplicable crash with inexplicable fix
            if chan == 'eeee' and varName == 'deltaPhiZZ':
                cTemp = Canvas(1000,1000)
                cTemp2 = Canvas(1000,1000)
                c.cd()


            (xaxis, yaxis), (xmin,xmax,ymin,ymax) = draw(toPlot, c,
                                                         xtitle=xTitle,
                                                         ytitle=yTitle,
                                                         )#logx=(varName=='Mass' and ana=='full'))
            # blinding box
            if blind and varName == 'Mass' and binning4l['Mass'][-1] > 500.:
                box = TBox(max(xmin,500.), ymin, min(binning4l['Mass'][-1], xmax), ymax)
                box.SetFillColor(1)
                box.SetFillStyle(3002)
                box.Draw("same")
                leg.SetFillStyle(1001)

            leg.Draw("same")

            style.setCMSStyle(c, '', dataType='Preliminary', intLumi=lumi)
            c.Print('{}/{}{}.png'.format(outdir, chan, varName))


            if varName == 'Mass' and ana == 'smp' and chan == 'zz':

                hTGCf4 = aTGCf4.makeHist(var, _selections4l[varName], binning,
                                         perUnitWidth=binNormWidth4l[varName],
                                         mergeOverflow=True)
                hTGCf5 = aTGCf5.makeHist(var, _selections4l[varName], binning,
                                         perUnitWidth=binNormWidth4l[varName],
                                         mergeOverflow=True)
                hSherpa = sherpa.makeHist(var, _selections4l[varName], binning,
                                          perUnitWidth=binNormWidth4l[varName],
                                          mergeOverflow=True)

                for s in stack:
                    if s.name not in ['ZZTo4L', 'ZZTo4L-amcatnlo']:
                        hToAdd = s.makeHist(var, _selections4l[varName], binning,
                                            postprocess=True,
                                            perUnitWidth=binNormWidth4l[varName],
                                            mergeOverflow=True)
                        hTGCf4 += hToAdd
                        hTGCf5 += hToAdd
                        hSherpa += hToAdd

                hTGCf4.SetLineWidth(2*hTGCf4.GetLineWidth())
                hTGCf5.SetLineWidth(2*hTGCf5.GetLineWidth())

                hSherpa.fillstyle = 'hollow'
                hSherpa.linecolor = 'black'
                hSherpa.SetLineWidth(2*hSherpa.GetLineWidth())
                hSherpa.linestyle = 'dashed'

                toPlot = [hStack, hTGCf4, hTGCf5, hSherpa, dataPts]

                legParams['entryheight'] = 0.02
                legParams['entrysep'] = 0.008
                legParams['textsize'] = 0.02
                legParams['leftmargin'] = 0.4

                cTGC = Canvas(1000,1000)
                cTGC.SetLogy()
                leg = makeLegend(cTGC, *toPlot, **legParams)

                (xaxis, yaxis), (xmin,xmax,ymin,ymax) = draw(toPlot, cTGC,
                                                             xtitle=xTitle,
                                                             ytitle=yTitle,
                                                             logy=True)
                # blinding box
                if blind and binning4l['Mass'][-1] > 500.:
                    box = TBox(max(xmin,500.), ymin, min(binning4l['Mass'][-1], xmax), ymax)
                    box.SetFillColor(1)
                    box.SetFillStyle(3002)
                    box.Draw("same")
                    leg.SetFillStyle(1001)

                leg.Draw("same")

                style.setCMSStyle(cTGC, '', dataType='Preliminary', intLumi=lumi)
                cTGC.Print('{}/{}{}_aTGC.png'.format(outdir, chan, varName))




    binning2l = _binning2l.copy()
    binNormWidth2l = _binNormWidth2l.copy()

    if ana != 'smp':
        binning2l['Mass'] = [60, 0., 120.]
        binNormWidth2l['Mass'] = 2.


    for z in ['z', 'ze', 'zm', 'z1', 'z2']:
        for varName, binning in binning2l.iteritems():
            print "Plotting {} {}".format(z, varName)

            var = {c:[vt.format(var=varName) for vt in _varTemplates2l[z][c]] for c in _varTemplates2l[z]}

            hStack = stack.makeHist(var, _selections2l[z], binning2l[varName],
                                    postprocess=True,
                                    perUnitWidth=binNormWidth2l[varName])
            dataPts = data.makeHist(var, _selections2l[z], binning2l[varName],
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

            legParams = {'textsize':0.029}
            if ana == 'full' and varName == 'Mass':
                legParams = _legParamsLeft.copy()
            leg = makeLegend(c, hStack, dataPts, **legParams)

            xTitle = _xTitles[varName]
            if '{obj}' in xTitle:
                xTitle = xTitle.format(obj=_objNames[z])

            yTitle = 'Z bosons / {} {}'.format(binNormWidth2l[varName],
                                               _units[varName])

            (xaxis, yaxis), (xmin,xmax,ymin,ymax) = draw([hStack, dataPts], c,
                                                         xtitle=xTitle,
                                                         ytitle=yTitle)

            # blinding box
            if blind and varName == 'Pt' and binning2l['Pt'][-1] > 200.:
                box = TBox(max(xmin,200.), ymin, min(binning2l['Pt'][-1], xmax), ymax)
                box.SetFillColor(1)
                box.SetFillStyle(3002)
                box.Draw("same")
                leg.SetFillStyle(1001)

            leg.Draw("same")

            style.setCMSStyle(c, '', dataType='Preliminary', intLumi=lumi)
            c.Print('{}/{}{}.png'.format(outdir, z, varName))


    for lep in _varTemplates1l:
        for varName, binning in _binning1l.iteritems():
            print "Plotting {} {}".format(lep, varName)

            if varName == 'Iso':
                varStr = 'ZZIso'
            else:
                varStr = varName

            var = {c:[vt.format(var=varStr) for vt in _varTemplates1l[lep][c]] for c in _varTemplates1l[lep]}

            hStack = stack.makeHist(var, _selections1l[lep], _binning1l[varName],
                                    postprocess=True,
                                    perUnitWidth=_binNormWidth1l[varName])
            dataPts = data.makeHist(var, _selections1l[lep], _binning1l[varName],
                                    poissonErrors=True,
                                    perUnitWidth=_binNormWidth1l[varName])

            # for ratio
            dataHist = data.makeHist(var, _selections1l[lep], _binning1l[varName],
                                     poissonErrors=True,
                                     perUnitWidth=_binNormWidth1l[varName])

            c = Canvas(1000,1000)

            leg = makeLegend(c, hStack, dataPts, leftmargin=0.47)

            xTitle = _xTitles[varName]
            if '{obj}' in xTitle:
                xTitle = xTitle.format(obj=_objNames[lep])

            yTitle = 'Leptons / {} {}'.format(_binNormWidth1l[varName],
                                              _units[varName])

            (xaxis, yaxis), (xmin,xmax,ymin,ymax) = draw([hStack, dataPts], c,
                                                         xtitle=xTitle,
                                                         ytitle=yTitle)
            leg.Draw("same")

            style.setCMSStyle(c, '', dataType='Preliminary', intLumi=lumi)
            c.Print('{}/{}{}.png'.format(outdir, lep, varName))


if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Make reco 4l plots")
    parser.add_argument('--dataDir', type=str, nargs='?',
                        default='uwvvNtuples_data_20feb2017',
                        help='Directory where data ntuples live')
    parser.add_argument('--mcDir', type=str, nargs='?',
                        default='uwvvNtuples_mc_20feb2017',
                        help='Directory where MC ntuples live')
    parser.add_argument('--plotDir', type=str, nargs='?',
                        default='zzPlots',
                        help='Directory to put plots in, possibly relative '
                        'to ~/www/UWVVPlots')
    parser.add_argument('--fakeRateFile', type=str, nargs='?',
                        default='fakeRate_20feb2017',
                        help=('Name of fake rate file (assumed to be in usual '
                              'data directory unless full path is specified)'))
    parser.add_argument('--puWeightFile', type=str, nargs='?',
                        default='puWeight_69200_24jan2017.root',
                        help=('Name of pileup weight file (assumed to be in usual '
                              'data directory unless full path is specified)'))
    parser.add_argument('--lumi', type=float, nargs='?', default=35860.,
                        help='Integrated luminosity of sample (in pb^-1)')
    parser.add_argument('--amcatnlo', action='store_true',
                        help='Use MadGraph5_aMC@NLO as the primary MC instead '
                        'of Powheg.')
    parser.add_argument('--eras', type=str, nargs='?', default='BCDEFGH',
                        help='Data eras to use.')
    parser.add_argument('--analysis', '--ana', type=str, nargs='?',
                        default='smp',
                        help='Which set of cuts to use (full, smp, etc.).')
    parser.add_argument('--sfHists', action='store_true',
                        help='Get lepton scale factors from files instead of directly from the ntuples.')
    parser.add_argument('--blind', action='store_true',
                        help='Put blinding boxes on a few distributions.')

    args=parser.parse_args()

    main(args.dataDir, args.mcDir, args.plotDir, args.analysis,
         args.fakeRateFile, args.puWeightFile, args.lumi, args.eras,
         args.blind, args.amcatnlo, args.sfHists)
