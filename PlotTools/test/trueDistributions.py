import logging
from rootpy import log as rlog; rlog = rlog["/trueDistributions"]
# don't show most silly ROOT messages
logging.basicConfig(level=logging.WARNING)
rlog["/ROOT.TUnixSystem.SetDisplay"].setLevel(rlog.ERROR)

from rootpy import asrootpy
from rootpy.io import root_open
from rootpy.plotting import Canvas, Graph
from rootpy.plotting.utils import draw

from PlotTools import PlotStyle as _Style, pdfViaTex
from PlotTools import makeLegend, makeErrorBand
from Utilities import deltaRString, deltaPhiString, combineWeights
from Analysis.setupStandardSamples import genZZSamples

from os import makedirs as _mkdir
from os.path import join as _join
from os.path import isdir as _isdir
from os.path import exists as _exists
from math import sqrt

plotDir = '/afs/cern.ch/user/n/nawoods/www/trueDistributions'
inMC = 'uwvvNtuples_mc_23mar2017'

_style = _Style()

channels = ['eeee','eemm', 'mmmm']

# set up variables, selection, binnings etc.
# (jet-related variables and selections done later)
_variables = {
    'pt' : {c:'Pt' for c in channels},
    'mass' : {c:'Mass' for c in channels},
    'eta' : {c:'abs(Eta)' for c in channels},
    'zHigherPt' : {
        'eeee' : ['e1_e2_Pt', 'e3_e4_Pt'],
        'mmmm' : ['m1_m2_Pt', 'm3_m4_Pt'],
        'eemm' : ['e1_e2_Pt','m1_m2_Pt']
        },
    'zLowerPt' : {
        'eeee' : ['e1_e2_Pt', 'e3_e4_Pt'],
        'mmmm' : ['m1_m2_Pt', 'm3_m4_Pt'],
        'eemm' : ['e1_e2_Pt','m1_m2_Pt']
        },
    'zPt' : {
        'eeee' : ['e1_e2_Pt', 'e3_e4_Pt'],
        'mmmm' : ['m1_m2_Pt', 'm3_m4_Pt'],
        'eemm' : ['e1_e2_Pt','m1_m2_Pt']
        },
    'deltaPhiZZ' : {
        'eeee' : 'abs({}(e1_e2_Phi, e3_e4_Phi))'.format(deltaPhiString()),
        'eemm' : 'abs({}(e1_e2_Phi, m1_m2_Phi))'.format(deltaPhiString()),
        'mmmm' : 'abs({}(m1_m2_Phi, m3_m4_Phi))'.format(deltaPhiString()),
        },
    'deltaRZZ' : {
        'eeee' : '{}(e1_e2_Eta, e1_e2_Phi, e3_e4_Eta, e3_e4_Phi)'.format(deltaRString()),
        'eemm' : '{}(e1_e2_Eta, e1_e2_Phi, m1_m2_Eta, m1_m2_Phi)'.format(deltaRString()),
        'mmmm' : '{}(m1_m2_Eta, m1_m2_Phi, m3_m4_Eta, m3_m4_Phi)'.format(deltaRString()),
        },
    'l1Pt' : {
        'eeee' : 'max(e1Pt, e3Pt)',
        'eemm' : 'max(e1Pt, m1Pt)',
        'mmmm' : 'max(m1Pt, m3Pt)',
        },
    }

_binning = {
    'pt' : [25.*i for i in range(4)] + [100., 150., 200., 300.],
    'nJets' : [5,-0.5,4.5],
    'mass' : [100.] + [200.+50.*i for i in range(5)] + [500.,600.,800.],
    'eta' : [6,0.,6.],
    'zPt' : [i * 25. for i in range(7)] + [200., 300.],
    'zHigherPt' : [i * 25. for i in range(7)] + [200., 300.],
    'zLowerPt' : [i * 25. for i in range(7)] + [200., 300.],
    'deltaPhiZZ' : [0., 1.5] + [2.+.25*i for i in range(6)],
    'deltaRZZ' : [6, 0., 6.],
    'l1Pt' : [0.,15.,30.,40.,50.]+[60.+15.*i for i in range(9)]+[195.,225.],#[14,0.,210.],#[15, 0., 150.],
    }

_units = {
    'pt' : 'GeV',
    'mass' : 'GeV',
    'eta' : '',
    'zPt' : 'GeV',
    'zHigherPt' : 'GeV',
    'zLowerPt' : 'GeV',
    'deltaPhiZZ' : '',
    'deltaRZZ' : '',
    'l1Pt' : 'GeV',
    }

_prettyVars = {
    'pt' : 'p_T^{\\text{ZZ}}',
    'mass' : 'm_{\\text{ZZ}}',
    'eta' : '\\eta_{\\text{ZZ}}',
    'zPt' : 'p_T^{\\text{Z}}',
    'zHigherPt' : 'p_T^{\\text{Z}_{\\text{lead}}}',
    'zLowerPt' : 'p_T^{\\text{Z}_{\\text{sublead}}}',
    'deltaPhiZZ' : '\\Delta \\phi_{\\text{Z}_1,\\text{Z}_2}',
    'deltaRZZ' : '\\Delta \\text{R}_{\\text{Z}_1,\\text{Z}_2}',
    'l1Pt' : 'p_{T}^{\\ell_1}',
    }

_xTitle = {}
_yTitle = {}
_yTitleNoNorm = {}
_yTitleTemp = '{prefix} \\frac{{d\\sigma_{{\\text{{fid}}}}}}{{d{xvar}}} {units}'
for var, prettyVar in _prettyVars.iteritems():
    xt = prettyVar
    if _units[var]:
        xt += ' \\, \\text{{[{}]}}'.format(_units[var])
        yt = _yTitleTemp.format(xvar=prettyVar,
                                prefix='\\frac{1}{\\sigma_{\\text{fid}}}',
                                units='\\, \\left[ \\frac{{1}}{{\\text{{{unit}}}}} \\right]'.format(unit=_units[var]))
        ytnn = _yTitleTemp.format(xvar=prettyVar, prefix='',
                                  units='\\, \\left[ \\frac{{\\text{{fb}}}}{{\\text{{{unit}}}}} \\right]'.format(unit=_units[var]))
    else:
        yt = _yTitleTemp.format(prefix='\\frac{1}{\\sigma_{\\text{fid}}}',
                                xvar=prettyVar, units='')
        ytnn = _yTitleTemp.format(prefix='', xvar=prettyVar, units='\\left[ \\text{fb} \\right]')

    _xTitle[var] = xt
    _yTitle[var] = yt
    _yTitleNoNorm[var] = ytnn

_selections = {
    'pt' : {c:'' for c in channels},
    'mass' : {c:'' for c in channels},
    'eta' : {c:'' for c in channels},
    'zPt' : {c:'' for c in channels},
    'zHigherPt' : {
        'eeee' : ['e1_e2_Pt > e3_e4_Pt', 'e1_e2_Pt < e3_e4_Pt'],
        'mmmm' : ['m1_m2_Pt > m3_m4_Pt', 'm1_m2_Pt < m3_m4_Pt'],
        'eemm' : ['e1_e2_Pt > m1_m2_Pt', 'e1_e2_Pt < m1_m2_Pt'],
        },
     'zLowerPt' : {
        'eeee' : ['e1_e2_Pt < e3_e4_Pt', 'e1_e2_Pt > e3_e4_Pt'],
        'mmmm' : ['m1_m2_Pt < m3_m4_Pt', 'm1_m2_Pt > m3_m4_Pt'],
        'eemm' : ['e1_e2_Pt < m1_m2_Pt', 'e1_e2_Pt > m1_m2_Pt'],
        },
    'deltaPhiZZ' : {c:'' for c in channels},
    'deltaRZZ' : {c:'' for c in channels},
    'l1Pt' : {c:'' for c in channels},
    }

_trueSelections = {
    'eeee' : 'e1_e2_Mass > 60. && e3_e4_Mass > 60.',
    'eemm' : 'e1_e2_Mass > 60. && m1_m2_Mass > 60.',
    'mmmm' : 'm1_m2_Mass > 60. && m3_m4_Mass > 60.',
    }

for var, selections in _selections.iteritems():
    for chan, sel in selections.iteritems():
        if isinstance(sel,str):
            sel = combineWeights(sel, _trueSelections[chan], selections=True)
        else:
            sel = [combineWeights(s, _trueSelections[chan], selections=True) for s in sel]


_legDefaults = {
    'textsize' : .025,
    'leftmargin' : 0.377,
    }
_legParams = {v:_legDefaults.copy() for v in _variables}
_legParams['deltaRZZ'] = {
    'textsize' : .015,
    'leftmargin' : .03,
    'rightmargin' : .47,
    'entryheight' : .023,
    'entrysep' : .007,
    }
_legParams['deltaPhiZZ']['leftmargin'] = 0.05
_legParams['deltaPhiZZ']['rightmargin'] = 0.32
_legParams['eta'] = _legParams['deltaRZZ'].copy()
_legParams['eta']['leftmargin'] = .5
_legParams['eta']['rightmargin'] = .03
_legParams['eta']['topmargin'] = .05
_legParams['l1Pt']['topmargin'] = 0.06

_matrixNames = {
    'deltaRZZ' : 'dR.Z_0.25__NNLO_QCD',
    'deltaPhiZZ' : 'dphi.Z_0.25__NNLO_QCD',
    'mass' : 'm.ZZ_5.0__NNLO_QCD',
    'nJets' : 'n_jets__NNLO_QCD',
    'pt' : 'pT.ZZ_2.5__NNLO_QCD',
    }
_matrixXSecs = {
    '' : 17.5413,
    'up' : 17.1878,
    'dn' : 17.9555,
    }
_matrixPath = '/afs/cern.ch/user/k/kelong/www/ZZMatrixDistributions'


_logy = {
    'pt' : True,
    'nJets' : False,
    'mass' : True,
    'eta' : False,
    'zPt' : True,
    'zHigherPt' : True,
    'zLowerPt' : True,
    'deltaPhiZZ' : False,
    'deltaRZZ' : False,
    'l1Pt' : True,
    }


def _normalizeBins(h):
    binUnit = 1 # min(h.GetBinWidth(b) for b in range(1,len(h)+1))
    for ib in xrange(1,len(h)+1):
        w = h.GetBinWidth(ib)
        h.SetBinContent(ib, h.GetBinContent(ib) * binUnit / w)
        h.SetBinError(ib, h.GetBinError(ib) * binUnit / w)
        if h.GetBinError(ib) > h.GetBinContent(ib):
            h.SetBinError(ib, h.GetBinContent(ib))
    h.sumw2()

def _unnormalizeBins(h):
    binUnit = 1 # min(h.GetBinWidth(b) for b in range(1,len(h)+1))
    for ib in xrange(1,len(h)+1):
        w = h.GetBinWidth(ib)
        h.SetBinContent(ib, h.GetBinContent(ib) * w / binUnit)
        h.SetBinError(ib, h.GetBinError(ib) * w / binUnit)
        if h.GetBinError(ib) > h.GetBinContent(ib):
            h.SetBinError(ib, h.GetBinContent(ib))
    h.sumw2()

def _pdfError(sample, var, sel, binning):
    hTrueVariations = []
    for s in sample.values():
        if 'GluGluZZ' not in s.name and 'phantom' not in s.name:
            hTrueVariations.append(s.makeHist2(var, 'Iteration$', sel, binning,
                                               [100,0.,100.], 'pdfWeights/pdfWeights[0]', False))
    allTrueRMSes = [[Graph(h.ProjectionY('slice{}'.format(i), i+1,i+1)).GetRMS(2) for i in xrange(h.GetNbinsX())] for h in hTrueVariations]
    binTrueRMSes = [sum(rmses) for rmses in zip(*allTrueRMSes)]

    out = sample.values()[0].makeHist('0','0',binning)

    for i in xrange(out.GetNbinsX()):
        out[i+1].value = binTrueRMSes[i]

    return out

def _scaleErrors(sample, var, sel, binning, hNominal):
    variationIndices = [1,2,3,4,6,8]
    hTrues = [sample.makeHist(var, sel, binning,
                              {
                'ZZTo4L':'scaleWeights[{}]/scaleWeights[0]'.format(i),
                'ZZTo4L-amcatnlo':'scaleWeights[{}]/scaleWeights[0]'.format(i),
                'ZZJJTo4L_EWK':'scaleWeights[{}]/scaleWeights[0]'.format(i),
                },
                              perUnitWidth=False)
              for i in variationIndices]

    # save true-level uncertainty for later
    outUp = hNominal.empty_clone()
    outDn = hNominal.empty_clone()
    for bUp, bDn, bNom, variations in zip(outUp, outDn, hNominal,
                                          zip(*hTrues)):
        bUp.value = max(b.value for b in variations) - bNom.value
        bDn.value = abs(min(b.value for b in variations) - bNom.value)

    return outUp, outDn



if not (_exists(_join(plotDir,'texs')) and _exists(_join(plotDir,'pdfs'))):
    _mkdir(_join(plotDir,'texs'))
    _mkdir(_join(plotDir,'pdfs'))
elif not _isdir(_join(plotDir,'texs')):
    raise IOError("There is already some non-directory object called {}.".format(_join(plotDir,'texs')))
elif not _isdir(_join(plotDir,'pdfs')):
    raise IOError("There is already some non-directory object called {}.".format(_join(plotDir,'pdfs')))

ana = 'smp'


sPow = genZZSamples('zz', inMC, ana, 1000.)
sMG5 = genZZSamples('zz', inMC, ana, 1000., amcatnlo=True)

powName = r'\textbf{POWHEG+MCFM+Pythia8}'
mg5Name = r'\textbf{MG5\_aMC@NLO+MCFM+Pythia8}'
matName = r'\textbf{MATRIX}'


for varName in _variables:
    var = _variables[varName]
    sel = _selections[varName]
    bins = _binning[varName]

    hNominalPow = {
        c : sPow[c].makeHist(var[c], sel[c], bins, perUnitWidth=False)
        for c in channels
        }
    hNominalPow[''] = sPow.makeHist(var, sel, bins, perUnitWidth=False)
    hNominalMG5 = {
        c : sMG5[c].makeHist(var[c], sel[c], bins, perUnitWidth=False)
        for c in channels
        }
    hNominalMG5[''] = sMG5.makeHist(var, sel, bins, perUnitWidth=False)

    hPDFErrPow = {
        c : _pdfError(sPow[c], var[c], sel[c], bins)
        for c in channels
        }
    hPDFErrPow[''] = sum(hPDFErrPow.values())
    hPDFErrMG5 = {
        c : _pdfError(sMG5[c], var[c], sel[c], bins)
        for c in channels
        }
    hPDFErrMG5[''] = sum(hPDFErrMG5.values())

    hScaleUpPow = {}
    hScaleDnPow = {}
    hScaleUpMG5 = {}
    hScaleDnMG5 = {}

    for c in channels:
        hScaleUpPow[c], hScaleDnPow[c] = _scaleErrors(sPow[c], var[c], sel[c],
                                                      bins, hNominalPow[c])
        hScaleUpMG5[c], hScaleDnMG5[c] = _scaleErrors(sMG5[c], var[c], sel[c],
                                                      bins, hNominalMG5[c])

    hScaleUpPow[''] = sum(hScaleUpPow.values())
    hScaleDnPow[''] = sum(hScaleDnPow.values())
    hScaleUpMG5[''] = sum(hScaleUpMG5.values())
    hScaleDnMG5[''] = sum(hScaleDnMG5.values())

    if varName in _matrixNames:
        with root_open(_join(_matrixPath,_matrixNames[varName]+'.root')) as fMat:
            cMat = asrootpy(fMat.canvas)
            hMat = asrootpy(cMat.FindObject(_matrixNames[varName])).clone()
            # scaleDown is for the UPPER error, scaleUp is for the LOWER
            hUpMat = asrootpy(cMat.FindObject(_matrixNames[varName]+'__scaleDown')).clone()
            hDnMat = asrootpy(cMat.FindObject(_matrixNames[varName]+'__scaleUp')).clone()
            hMat.SetDirectory(0)
            hUpMat.SetDirectory(0)
            hDnMat.SetDirectory(0)

        # un-normalize the bins, rebin, renormalize
        _unnormalizeBins(hMat)
        hMat = hMat.rebinned([e for e in hNominalPow['']._edges(0)])
        _normalizeBins(hMat)
        _unnormalizeBins(hUpMat)
        hUpMat = hUpMat.rebinned([e for e in hNominalPow['']._edges(0)])
        _normalizeBins(hUpMat)
        _unnormalizeBins(hDnMat)
        hDnMat = hDnMat.rebinned([e for e in hNominalPow['']._edges(0)])
        _normalizeBins(hDnMat)

        hMat /= _matrixXSecs['']
        hUpMat /= _matrixXSecs['up']
        hDnMat /= _matrixXSecs['dn']

        hUpMat -= hMat
        hDnMat -= hMat

        hMat.title = r'MATRIX'
        hMat.color = 'forestgreen'
        hMat.drawstyle = 'hist'
        hMat.fillstyle = 'hollow'
        hMat.legendstyle = 'L'
        hMat.SetLineWidth(hMat.GetLineWidth()*2)

        errorBandMat = makeErrorBand(hMat, hUpMat, hDnMat)
        errorBandMat.fillstyle = 'solid'
        errorBandMat.SetFillColorAlpha(819, 0.5)

        legMat = hMat.clone()
        legMat.SetFillColorAlpha(819, 0.5)
        legMat.fillstyle = 'solid'
        legMat.legendstyle = 'LF'

    for c in hNominalPow:
        hPow = hNominalPow[c].clone()
        hPow.linecolor = '#000099'
        hPow.drawstyle = 'hist'
        hPow.fillstyle = 'hollow'
        hPow.legendstyle = 'L'
        hPow.title = '{}'.format(powName)
        hPow.SetLineWidth(hPow.GetLineWidth()*2)

        hMG5 = hNominalMG5[c].clone()
        hMG5.linecolor = 'red'
        hMG5.drawstyle = 'hist'
        hMG5.fillstyle = 'hollow'
        hMG5.legendstyle = 'L'
        hMG5.title = '{}'.format(mg5Name)
        hMG5.SetLineWidth(hMG5.GetLineWidth()*2)

        # uncertainty
        hUncUpPow = hPDFErrPow[c].empty_clone()
        hUncDnPow = hPDFErrPow[c].empty_clone()
        for bUp, bDn, bPDF, bScaleUp, bScaleDn in zip(hUncUpPow,
                                                      hUncDnPow,
                                                      hPDFErrPow[c],
                                                      hScaleUpPow[c],
                                                      hScaleDnPow[c]):
            bUp.value = sqrt(bPDF.value**2 + bScaleUp.value**2)
            bDn.value = sqrt(bPDF.value**2 + bScaleDn.value**2)

        intPow = hPow.Integral(0,hPow.GetNbinsX()+1)
        hPow /= intPow
        hUncUpPow /= (intPow + hUncUpPow.Integral(0,hUncUpPow.GetNbinsX()+1))
        hUncDnPow /= (intPow + hUncDnPow.Integral(0,hUncDnPow.GetNbinsX()+1))

        _normalizeBins(hPow)
        _normalizeBins(hUncUpPow)
        _normalizeBins(hUncDnPow)

        errorBandPow = makeErrorBand(hPow, hUncUpPow, hUncDnPow)
        errorBandPow.fillstyle = 'solid'
        errorBandPow.SetFillColorAlpha(600, 0.5)

        toPlot = [errorBandPow, hPow]
        legPow = hPow.clone()
        legPow.SetFillColorAlpha(600, 0.5)
        legPow.fillstyle = 'solid'
        legPow.legendstyle = 'LF'
        forLegend = [legPow]

        hUncUpMG5 = hPDFErrMG5[c].empty_clone()
        hUncDnMG5 = hPDFErrMG5[c].empty_clone()
        for bUp, bDn, bPDF, bScaleUp, bScaleDn in zip(hUncUpMG5,
                                                      hUncDnMG5,
                                                      hPDFErrMG5[c],
                                                      hScaleUpMG5[c],
                                                      hScaleDnMG5[c]):
            bUp.value = sqrt(bPDF.value**2 + bScaleUp.value**2)
            bDn.value = sqrt(bPDF.value**2 + bScaleDn.value**2)

        intMG5 = hMG5.Integral(0,hMG5.GetNbinsX()+1)
        hMG5 /= intMG5
        hUncUpMG5 /= (intMG5 + hUncUpMG5.Integral(0,hUncUpMG5.GetNbinsX()+1))
        hUncDnMG5 /= (intMG5 + hUncDnMG5.Integral(0,hUncDnMG5.GetNbinsX()+1))

        _normalizeBins(hMG5)
        _normalizeBins(hUncUpMG5)
        _normalizeBins(hUncDnMG5)

        errorBandMG5 = makeErrorBand(hMG5, hUncUpMG5, hUncDnMG5)
        errorBandMG5.fillstyle = 'solid'
        errorBandMG5.SetFillColorAlpha(628, 0.5)

        toPlot += [errorBandMG5, hMG5]
        legMG5 = hMG5.clone()
        legMG5.SetFillColorAlpha(628, 0.5)
        legMG5.fillstyle = 'solid'
        legMG5.legendstyle = 'LF'
        forLegend.append(legMG5)

        if varName in _matrixNames:
            toPlot += [errorBandMat, hMat]
            forLegend.append(legMat)


        drawOpts = {
            'xtitle' : _xTitle[varName],
            }
        drawOpts['ytitle'] = _yTitle[varName]
        drawOpts['logy'] = _logy[varName]
        drawOpts['yerror_in_padding'] = _logy[varName]
        drawOpts['ypadding'] = 0.04

        if varName == 'l1Pt':
            drawOpts['logy_crop_value'] = 1e-4

        can = Canvas(1000,1000)

        (xaxis, yaxis), (xmin,xmax,ymin,ymax) = draw(toPlot, can, **drawOpts)
        yaxis.SetTitleSize(0.7*yaxis.GetTitleSize())
        yaxis.SetTitleOffset(1.15*yaxis.GetTitleOffset())
        yaxis.SetLabelSize(0.72*yaxis.GetLabelSize())
        yaxis.SetMoreLogLabels()
        xaxis.SetLabelSize(0.82*xaxis.GetLabelSize())
        xaxis.SetTitleOffset(0.95*xaxis.GetTitleOffset())

        leg = makeLegend(can, *forLegend, **_legParams[varName])
        leg.Draw("same")

        _style.setCMSStyle(can, '', dataType='Preliminary Simulation',
                           intLumi=-1)
        outWords = ['true',varName,c]
        pdfViaTex(can, '_'.join(w for w in outWords if w),
                  _join(plotDir,'texs'), _join(plotDir, 'pdfs'))
