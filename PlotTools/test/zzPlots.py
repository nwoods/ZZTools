
import logging
from rootpy import log as rlog; rlog = rlog["/zzPlots"]
# don't show most silly ROOT messages
logging.basicConfig(level=logging.WARNING)
rlog["/ROOT.TUnixSystem.SetDisplay"].setLevel(rlog.ERROR)

from rootpy import asrootpy
from rootpy.io import root_open
from rootpy.plotting import Canvas, Legend, Graph
from rootpy.plotting.utils import draw, get_limits
from rootpy.ROOT import TBox, Double

from SampleTools import MCSample, DataSample, SampleGroup, SampleStack
from PlotTools import PlotStyle as _Style, pdfViaTex as _pdfViaTex
from PlotTools import makeLegend, addPadBelow, makeRatio, fixRatioAxes, \
    makeErrorBand
from Utilities import WeightStringMaker, deltaRString, deltaPhiString, \
    makeNumberPretty, parseChannels
from Analysis.setupStandardSamples import *
from Analysis.weightHelpers import baseMCWeight

from os import environ
from os import path as _path
from os import makedirs as _mkdir
from os.path import isdir as _isdir
from os.path import exists as _exists
from collections import OrderedDict
from math import sqrt



_objNames = {
    'zz' : '4\\ell',
    'eeee' : '4e',
    'eemm' : '2e2\\mu',
    'mmmm' : '4\\mu',
    'z' : '\\text{Z}',
    'z1' : '\\text{Z}_{1}',
    'z2' : '\\text{Z}_{2}',
    'ze' : '\\text{Z} \\rightarrow e^{+}e^{-}',
    'zm' : '\\text{Z} \\rightarrow \\mu^{+}\\mu^{-}',
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
    'Pt' : '{obj} \\, p_\\text{{T}} \\, (\\text{{GeV}})',
    'nJets' : 'N_{\\text{jets}}',
    'nJets_eta2p4' : 'N_{\\text{jets}} \\left( \\left|\\eta\\right| < 2.4 \\right)',
    'Iso' : 'R_{{Iso}} \\, ({obj})',
    'PVDXY' : '\\left| \\Delta_{{xy}} \\, ({obj}) \\right| \\, (\\text{{cm}})',
    'PVDZ' : '\\left| \\Delta_{{z}} \\, ({obj}) \\right| \\, (\\text{{cm}})',
    'nvtx' : 'N_{\\text{vtx}}',
    'SIP3D' : 'SIP_{{3D}} \\, ({obj})',
    'jet1Pt' : 'p_\\text{T}^\\text{j1} \\, (\\text{GeV})',
    'jet1Eta' : '\\eta_\\text{j1}',
    'jet2Pt' : 'p_\\text{T}^\\text{j2} \\, (\\text{GeV})',
    'jet2Eta' : '\\eta_\\text{j2}',
    'mjj' : 'm_\\text{jj} \\, (\\text{GeV})',
    'deltaEtajj' : '|\\Delta \\eta_{\\text{jj}}}|',
    'deltaPhiZZ' : '\\Delta \\phi (\\text{Z}_1, \\text{Z}_2)',
    'deltaRZZ' : '\\Delta \\text{R} (\\text{Z}_1, \\text{Z}_2)',
    }

for v,t in _xTitles.iteritems():
    if _units[v]:
        t += ' \\, (\\text{{{{{}}}}})'.format(_units[v])

# some distributions need the legend moved
_legParamsLeft = {
    'leftmargin' : 0.03,
    'rightmargin' : 0.48,
    'textsize' : 0.034,
    'topmargin' : 0.05,
    }


_binning4l = {
    'Mass'  : [100.] + [200.+50.*i for i in range(5)] + [500.,600.,800., 1000., 1200.],
    'Pt'    : [0.,5.]+[25.+25.*i for i in range(3)] + [100., 150., 200., 300.], #[25.*i for i in range(4)] + [100., 150., 200., 300.],
    'Eta'   : [16, -5., 5.],
    #'Phi'   : [12, -3.15, 3.15],
    'nvtx'  : [40, 0., 40.],
    'nJets' : [5, -0.5, 4.5],
    'nJets_eta2p4' : [6, -0.5, 5.5],
    'jet1Pt' : [30., 50., 100., 200., 300., 500.],
    'jet1Eta' : [0., 1.5, 2.4, 3.2, 4.7],
    'jet2Pt' : [30., 50., 100., 170., 300.],
    'jet2Eta' : [0., 1.5, 3., 4.7],
    'mjj' : [0., 200., 400., 600., 1000.],
    'deltaEtajj' : [4, 0.,4.7],
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
    'jet1Pt' : False, #50.,
    'jet2Pt' : 50.,
    'jet1Eta' : 1.,
    'jet2Eta' : 1.,
    'mjj' : 100.,
    'deltaEtajj' : 1.,
    'deltaPhiZZ' : 1.,
    'deltaRZZ' : 1.,
    }

_nDivisions4l = {
    'full' : {
        'Mass' : -408,
        },
    'z4l' : {
        'Mass' : -210,
        },
    'smp' : {
        },
    }

_nDivisions2l = {
    'full' : {
        'Mass' : -506,
        },
    'z4l' : {
        },
    'smp' : {
        },
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

for jSys in 'jer', 'jes':
    for shift in 'Up', 'Down':
        for nj in 1,2:
            _vars4l['jet{}Pt_{}{}'.format(nj,jSys,shift)] = {c:v.replace('jetPt',
                                                                         'jetPt_{}{}'.format(jSys,shift))
                                                             for c,v in _vars4l['jet{}Pt'.format(nj)].iteritems()}
            _selections4l['jet{}Pt_{}{}'.format(nj,jSys,shift)] = _selections4l['jet{}Pt'.format(nj)].replace('nJets',
                                                                                                              'nJets_{}{}'.format(jSys,shift))
            _vars4l['jet{}Eta_{}{}'.format(nj,jSys,shift)] = {c:v.replace('jetEta',
                                                                          'jetEta_{}{}'.format(jSys,shift))
                                                              for c,v in _vars4l['jet{}Eta'.format(nj)].iteritems()}
            _selections4l['jet{}Eta_{}{}'.format(nj,jSys,shift)] = _selections4l['jet{}Eta'.format(nj)].replace('nJets',
                                                                                                                'nJets_{}{}'.format(jSys,shift))

        _vars4l['mjj_{}{}'.format(jSys,shift)] = {c:v.replace('mjj',
                                                                 'mjj_{}{}'.format(jSys,shift))
                                                     for c,v in _vars4l['mjj'].iteritems()}
        _selections4l['mjj_{}{}'.format(jSys,shift)] = _selections4l['mjj'].replace('nJets',
                                                                                    'nJets_{}{}'.format(jSys,shift))
        _vars4l['deltaEtajj_{}{}'.format(jSys,shift)] = {c:v.replace('deltaEtajj',
                                                                     'deltaEtajj_{}{}'.format(jSys,shift))
                                                     for c,v in _vars4l['deltaEtajj'].iteritems()}
        _selections4l['deltaEtajj_{}{}'.format(jSys,shift)] = _selections4l['deltaEtajj'].replace('nJets',
                                                                                                  'nJets_{}{}'.format(jSys,shift))
        _vars4l['nJets_{}{}'.format(jSys,shift)] = {c:v.replace('nJets',
                                                                'nJets_{}{}'.format(jSys,shift))
                                                    for c,v in _vars4l['nJets'].iteritems()}
        _selections4l['nJets_{}{}'.format(jSys,shift)] = ''

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
    'PVDXY' : [50, 0., .5],
    'PVDZ' : [50, 0., 1.],
    'SIP3D' : [40, 0., 10.],
    }

_binNormWidth1l = {v:abs(b[2]-b[1])/b[0] for v,b in _binning1l.iteritems()}

_varTemps1l = {v:'{obj}'+v for v in _binning1l}
_varTemps1l['PVDXY'] = 'abs({obj}PVDXY)'
_varTemps1l['PVDZ'] = 'abs({obj}PVDZ)'
_varTemps1l['Iso'] = '{obj}ZZIso'


ze1Leps = ['e1','e2']
ze2Leps = ['e3','e4']
zm1Leps = ['m1','m2']
zm2Leps = ['m3','m4']

_vars1l = {
    v : {
        'l' : {
            'eeee' : [vt.format(obj=ob) for ob in ze1Leps + ze2Leps],
            'eemm' : [vt.format(obj=ob) for ob in ze1Leps + zm1Leps],
            'mmmm' : [vt.format(obj=ob) for ob in zm1Leps + zm2Leps],
            },
        'l1' : {
            'eeee' : [vt.format(obj=ze1Leps[0]), vt.format(obj=ze2Leps[0])],
            'eemm' : [vt.format(obj=ze1Leps[0]), vt.format(obj=zm1Leps[0])],
            'mmmm' : [vt.format(obj=zm1Leps[0]), vt.format(obj=zm2Leps[0])],
            },
        'e' : {
            'eeee' : [vt.format(obj=ob) for ob in ze1Leps + ze2Leps],
            'eemm' : [vt.format(obj=ob) for ob in ze1Leps],
            },
        'm' : {
            'mmmm' : [vt.format(obj=ob) for ob in zm1Leps + zm2Leps],
            'eemm' : [vt.format(obj=ob) for ob in zm1Leps],
            },
        } for v, vt in _varTemps1l.iteritems()
    }

_selections1l = {l:'' for l in _vars1l.values()[0]}
_selections1l['l1'] = {
    'eeee' : ['e1Pt > e3Pt', 'e3Pt > e1Pt'],
    'eemm' : ['e1Pt > m1Pt', 'm1Pt > e1Pt'],
    'mmmm' : ['m1Pt > m3Pt', 'm3Pt > m1Pt'],
    }


def _makeSystematics(varName, var, sel, binning,
                     sig, bkg, irr,
                     sigSyst, bkgSyst, irrSyst,
                     norm, puWeightFile, sfArgs={}):
    '''
    Get histograms of the systematic error sizes.

    Arguments:
    varName (str): Short name of the variable.
    var (dict): The variables to be used for each channel (channels are the
        keys)
    sel (dist): Selection to use, if any.
    binning (list of numbers): Binning to use.
    sig (SampleGroup): Signal samples keyed to channels. Should include qq, gg,
        and EWK samples if applicable.
    bkg (SampleGroup): Reducible background samples keyed to channels.
    irr (SampleGroup): Irreducible background samples keyed to channels.
    sigSyst (dict of SampleGroups): Signal samples with systematic shifts,
        keyed to the names of the shifts.
    bkgSyst (dict of SampleGroups): Reducible background samples with
        systematic shifts, keyed to the names of the shifts.
    irrSyst (dict of SampleGroups): Irreducible background samples with
        systematic shifts, keyed to the names of the shifts.
    norm (bool or float): Whether/how to normalize the bins to their width.
    puWeightFile (str): File with pileup weight scale factors to use.
    sfArgs (dict, optional): Files with efficiency scale factor histograms to
        use, if any are desired (otherwise, SFs will be taken from ntuples).

    Return:
    (hErrUp, hErrDown), histograms of the up- and down asymmetric systematics.
    '''
    nominalWeight = baseMCWeight('zz', puWeightFile, **sfArgs)

    hSigSyst = {}
    hBkgSyst = {}
    hIrrSyst = {}

    hSigNom = sig.makeHist(var, sel, binning, perUnitWidth=norm)
    hBkgNom = bkg.makeHist(var, sel, binning, postprocess=True,
                           perUnitWidth=norm)
    hIrrNom = irr.makeHist(var, sel, binning, perUnitWidth=norm)

    # PU weight
    hSigSyst['pu'] = {}
    hIrrSyst['pu'] = {}
    for sys in 'up','dn':
        wtStr = baseMCWeight('zz', puWeightFile, puSyst=sys, **sfArgs)

        sig.applyWeight(wtStr, True)
        irr.applyWeight(wtStr, True)
        hSigSyst['pu'][sys] = sig.makeHist(var, sel, binning,
                                           perUnitWidth=norm)
        hIrrSyst['pu'][sys] = irr.makeHist(var, sel, binning,
                                           perUnitWidth=norm)
        sig.applyWeight(nominalWeight, True)
        irr.applyWeight(nominalWeight, True)

    # lepton efficiency
    for lep in set(''.join(var.keys())):
        hSigSyst[lep+'Syst'] = {}
        hIrrSyst[lep+'Syst'] = {}
        for sys in ['up','dn']:
            wtArg = {lep+'Syst':sys}
            wtArg.update(sfArgs)
            wtStr = baseMCWeight('zz', puWeightFile,
                                 **wtArg)
            sig.applyWeight(wtStr, True)
            irr.applyWeight(wtStr, True)
            hSigSyst[lep+'Syst'][sys] = sig.makeHist(var, sel, binning,
                                                     perUnitWidth=norm)
            hIrrSyst[lep+'Syst'][sys] = irr.makeHist(var, sel, binning,
                                                     perUnitWidth=norm)

    sig.applyWeight(nominalWeight, True)
    irr.applyWeight(nominalWeight, True)

    # luminosity
    hSigSyst['lumi'] = {'up':hSigNom * 1.025,'dn':hSigNom * 0.975}
    hIrrSyst['lumi'] = {'up':hIrrNom * 1.025,'dn':hIrrNom * 0.975}

    relevantLeptons = set(''.join(var.keys()))

    # lepton fake rate
    for lep in relevantLeptons:
        relevantVar = {c:var[c] for c in var if lep in c}
        irrelevantVar = {c:var[c] for c in var if lep not in c}
        hBkgSyst[lep+'FR'] = {}
        for sys in ['up','dn']:
            hBkgSyst[lep+'FR'][sys] = bkgSyst[lep+sys].makeHist(
                relevantVar, sel, binning, perUnitWidth=norm,
                postprocess=True
                )
            if irrelevantVar:
                hBkgSyst[lep+'FR'][sys] += bkg.makeHist(
                    irrelevantVar, sel, binning, perUnitWidth=norm,
                    postprocess=True
                    )


    # lepton momentum
    eVar = {c:var[c] for c in var if 'e' in c}
    noEVar = {c:var[c] for c in var if 'e' not in c}
    if eVar:
        for sys in ['eScale', 'eRhoRes', 'ePhiRes']:
            hSigSyst[sys] = {}
            hIrrSyst[sys] = {}
            for shift in ['up','dn']:
                if sys == 'ePhiRes' and shift == 'dn':
                    continue
                sysStr = 'Up' if shift == 'up' else 'Dn'

                hSigSyst[sys][shift] = sigSyst[sys+sysStr].makeHist(
                    eVar, sel, binning, perUnitWidth=norm
                    )
                hIrrSyst[sys][shift] = irrSyst[sys+sysStr].makeHist(
                    eVar, sel, binning, perUnitWidth=norm
                    )
                if noEVar:
                    hSigSyst[sys][shift] += sig.makeHist(noEVar, sel, binning,
                                                         perUnitWidth=norm)
                    hIrrSyst[sys][shift] += irr.makeHist(noEVar, sel, binning,
                                                         perUnitWidth=norm)

        hSigSyst['ePhiRes']['dn'] = (hSigNom * 2) - hSigSyst['ePhiRes']['up']
        hIrrSyst['ePhiRes']['dn'] = (hIrrNom * 2) - hIrrSyst['ePhiRes']['up']


    muVar = {c:var[c] for c in var if 'm' in c}
    noMuVar = {c:var[c] for c in var if 'm' not in c}
    sys = 'mClosure'
    if muVar:
        hSigSyst[sys] = {}
        hIrrSyst[sys] = {}
        for shift in ['up','dn']:
            sysStr = 'Up' if shift == 'up' else 'Dn'

            hSigSyst[sys][shift] = sigSyst[sys+sysStr].makeHist(
                muVar, sel, binning, perUnitWidth=norm
                )
            hIrrSyst[sys][shift] = irrSyst[sys+sysStr].makeHist(
                muVar, sel, binning, perUnitWidth=norm
                )
            if noMuVar:
                hSigSyst[sys][shift] += sig.makeHist(noMuVar, sel, binning,
                                                     perUnitWidth=norm)
                hIrrSyst[sys][shift] += irr.makeHist(noMuVar, sel, binning,
                                                     perUnitWidth=norm)

    # PDF
    hSigSyst['pdf'] = {}
    hSigVariations = []
    hGG = None
    hEWK = None
    for c, v in var.iteritems():
        if isinstance(sel, dict) and c in sel:
            thisSel = sel[c]
        else:
            thisSel = sel
        for s in sig[c].values():
            if 'GluGluZZ' in s.name:
                hThisGG = s.makeHist(v, thisSel, binning, perUnitWidth=norm)
                hThisGG.ClearUnderflowAndOverflow()
                if hGG is None:
                    hGG = hThisGG
                else:
                    hGG += hThisGG
                if 'up' not in hSigSyst['pdf']:
                    hSigSyst['pdf']['up'] = hThisGG.clone() # * 1.18
                else:
                    hSigSyst['pdf']['up'] += hThisGG # * 1.18
                if 'dn' not in hSigSyst['pdf']:
                    hSigSyst['pdf']['dn'] = hThisGG.clone() # * 0.82
                else:
                    hSigSyst['pdf']['dn'] += hThisGG # * 0.82
            elif 'phantom' in s.name:
                hThisEWK = s.makeHist(v, thisSel, binning,
                                      perUnitWidth=norm)
                if hEWK is None:
                    hEWK = hThisEWK
                else:
                    hEWK += hThisEWK
                if 'up' not in hSigSyst['pdf']:
                    hSigSyst['pdf']['up'] = hThisEWK.clone()
                else:
                    hSigSyst['pdf']['up'] += hThisEWK
                if 'dn' not in hSigSyst['pdf']:
                    hSigSyst['pdf']['dn'] = hThisEWK.clone()
                else:
                    hSigSyst['pdf']['dn'] += hThisEWK
            else:
                hSigVariations.append(
                    s.makeHist2(
                        v, 'Iteration$', thisSel, binning, [100,0.,100.],
                        'pdfWeights / pdfWeights[0]', False
                        )
                    )
                hThisQQ = s.makeHist(v, thisSel, binning,
                                     perUnitWidth=norm)
                if 'up' not in hSigSyst['pdf']:
                    hSigSyst['pdf']['up'] = hThisQQ.clone()
                else:
                    hSigSyst['pdf']['up'] += hThisQQ
                if 'dn' not in hSigSyst['pdf']:
                    hSigSyst['pdf']['dn'] = hThisQQ.clone()
                else:
                    hSigSyst['pdf']['dn'] += hThisQQ



    # for each var bin in each sample, get the RMS across all the variations
    allSigRMSes = [
        [
            Graph(
                h.ProjectionY('slice{}'.format(i), i+1,i+1)
                ).GetRMS(2)
            for i in xrange(h.GetNbinsX())
            ] for h in hSigVariations
        ]

    # for each var bin, add variations for all samples
    sigBinRMSes = [sum(rmses) for rmses in zip(*allSigRMSes)]

    # apply variations
    for i in xrange(hSigNom.GetNbinsX()):
        hSigSyst['pdf']['up'][i+1].value += sigBinRMSes[i]
        hSigSyst['pdf']['dn'][i+1].value = max(0.,hSigSyst['pdf']['dn'][i+1].value - sigBinRMSes[i])


    # QCD scale uncertainties
    variationIndices = [1,2,3,4,6,8]
    hSigVars = []
    hIrrVars = []
    for ind in variationIndices:
        hSigVars.append(hGG.empty_clone())
        for c, v in var.iteritems():
            if isinstance(sel, dict) and c in sel:
                thisSel = sel[c]
            else:
                thisSel = sel
            for s in sig[c].values():
                if 'GluGluZZ' not in s.name and 'phantom' not in s.name:
                    hSigVars[-1] += s.makeHist(v, thisSel, binning,
                                               'scaleWeights[{}] / scaleWeights[0]'.format(ind),
                                               perUnitWidth=norm)
        hIrrVars.append(irr.makeHist(var, sel, binning,
                                     'scaleWeights[{}] / scaleWeights[0]'.format(ind),
                                     perUnitWidth=norm))

    hSigScaleUp = hSigNom.empty_clone()
    hSigScaleDn = hSigNom.empty_clone()
    for bUp, bDn, bVars in zip(hSigScaleUp, hSigScaleDn, zip(*hSigVars)):
        # print [b.value for b in bVars]
        bUp.value = max(b.value for b in bVars)
        bDn.value = min(b.value for b in bVars)
    hSigScaleUp += hGG# * 1.23
    hSigScaleDn += hGG# * 0.82
    if hEWK:
        hSigScaleUp += hEWK
        hSigScaleDn += hEWK

    hIrrScaleUp = hIrrNom.empty_clone()
    hIrrScaleDn = hIrrNom.empty_clone()
    for bUp, bDn, bVars in zip(hIrrScaleUp, hIrrScaleDn, zip(*hIrrVars)):
        bUp.value = max(b.value for b in bVars)
        bDn.value = min(b.value for b in bVars)

    hSigSyst['scale'] = {'up':hSigScaleUp,'dn':hSigScaleDn}
    hIrrSyst['scale'] = {'up':hIrrScaleUp,'dn':hIrrScaleDn}


    # alpha_s uncertainties
    alphaSIndices = [100,101]
    hSigVars = []
    for ind in alphaSIndices:
        hSigVars.append(hSigNom.empty_clone())
        for c, v in var.iteritems():
            if isinstance(sel, dict) and c in sel:
                thisSel = sel[c]
            else:
                thisSel = sel
            for s in sig[c].values():
                if 'GluGluZZ' not in s.name and 'phantom' not in s.name:
                    hSigVars[-1] += s.makeHist(v, thisSel, binning,
                                               'pdfWeights[{}]'.format(ind),
                                               perUnitWidth=norm)

    aSUnc = hSigVars[0] - hSigVars[-1]
    aSUnc /= 2.
    for b in aSUnc:
        b.value = abs(b.value)

    hSigSyst['alphaS'] = {
        'up' : hSigNom + aSUnc,
        'dn' : hSigNom - aSUnc,
        }


    # We don't have LHE information for MCFM, so just vary its
    # normalization
    hSigSyst['mcfm'] = {
        'up':hSigNom + hGG*0.18,
        'dn':hSigNom - hGG*0.15,
        }

    # JES/JER
    if ('jet' in varName.lower() or 'jj' in varName.lower()) and 'eta2p4' not in varName.lower():
        for sys in ['jer','jes']:
            hSigSyst[sys] = {}
            hIrrSyst[sys] = {}

            for shift in ['up','dn']:
                sysStr = 'Up' if shift == 'up' else 'Down'

                shiftedVarName = varName + '_' + sys + sysStr
                varShifted = {c:_vars4l[shiftedVarName][c] for c in var}
                selShifted = _selections4l[shiftedVarName]

                hSigSyst[sys][shift] = sig.makeHist(varShifted,
                                                    selShifted,
                                                    binning,
                                                    perUnitWidth=norm)
                hIrrSyst[sys][shift] = irr.makeHist(varShifted,
                                                    selShifted,
                                                    binning,
                                                    perUnitWidth=norm)


    # print "Signal:"
    # print "    Nominal: {}".format(hSigNom.Integral())
    # for sysName, sys in hSigSyst.iteritems():
    #     print "    "+sysName+":"
    #     for upperdown, h in sys.iteritems():
    #         print "        {}: {}".format(upperdown, h.Integral())
    # print "\nReducible Background:"
    # print "    Nominal: {}".format(hBkgNom.Integral())
    # for sysName, sys in hBkgSyst.iteritems():
    #     print "    "+sysName+":"
    #     for upperdown, h in sys.iteritems():
    #         print "        {}: {}".format(upperdown, h.Integral())
    # print "\nIrreducible Background:"
    # print "    Nominal: {}".format(hIrrNom.Integral())
    # for sysName, sys in hIrrSyst.iteritems():
    #     print "    "+sysName+":"
    #     for upperdown, h in sys.iteritems():
    #         print "        {}: {}".format(upperdown, h.Integral())

    # put it all together
    hUp = hSigNom.empty_clone()
    hDn = hSigNom.empty_clone()

    for sys in set(hSigSyst.keys() + hBkgSyst.keys() + hIrrSyst.keys()):
        thisUp = hSigSyst.get(sys, {}).get('up', hSigNom) - hSigNom + \
            hBkgSyst.get(sys, {}).get('up', hBkgNom) - hBkgNom + \
            hIrrSyst.get(sys, {}).get('up', hIrrNom) - hIrrNom
        thisDn = hSigSyst.get(sys, {}).get('dn', hSigNom) - hSigNom + \
            hBkgSyst.get(sys, {}).get('dn', hBkgNom) - hBkgNom + \
            hIrrSyst.get(sys, {}).get('dn', hIrrNom) - hIrrNom

        for bUp, bDn, b1, b2 in zip(hUp, hDn, thisUp, thisDn):
            bUp.value += max(b1.value,b2.value)**2
            bDn.value += min(b1.value,b2.value)**2


    for bUp, bDn in zip(hUp, hDn):
        bUp.value = sqrt(bUp.value)
        bDn.value = sqrt(bDn.value)

    return hUp, hDn


def main(inData, inMC, plotDir, ana, fakeRateFile, puWeightFile, lumi,
         eras='BCDEFGH', blind=False, amcatnlo=False, doSyst=True, logy=False,
         looseSIP=False, noSIP=False, sfRemake=False, paper=False):

    sfArgs = {}
    sipForBkg = 4.
    if noSIP:
        sfArgs['eSelSFFile'] = 'eleSelectionSF_HZZ_NWRemake_NoSIP'
        sfArgs['eSelSFFileGap'] = 'eleSelectionSFGap_HZZ_NWRemake_NoSIP'
        sfArgs['eRecoSFFile'] = 'eleRecoSF_HZZ_Moriond17'
        sfArgs['mSFFile'] = 'muSelectionAndRecoSF_HZZ_Moriond17_NoSIP'
        sipForBkg = -1
        if looseSIP:
            raise ValueError("You can use scale factors for loose SIP cut or "
                             "no SIP cut, but not both.")
    elif looseSIP:
        sfArgs['eSelSFFile'] = 'eleSelectionSF_HZZ_NWRemake_LooseSIP'
        sfArgs['eSelSFFileGap'] = 'eleSelectionSFGap_HZZ_NWRemake_LooseSIP'
        sfArgs['mSFFile'] = 'muSelectionAndRecoSF_HZZ_Moriond17_LooseSIP'
        sfArgs['eRecoSFFile'] = 'eleRecoSF_HZZ_Moriond17'
        sipForBkg = 10.
    elif sfRemake:
        sfArgs['eSelSFFile'] = 'eleSelectionSF_HZZ_NWRemake'
        sfArgs['eSelSFFileGap'] = 'eleSelectionSFGap_HZZ_NWRemake'
        sfArgs['eRecoSFFile'] = 'eleRecoSF_HZZ_Moriond17'
        sfArgs['mSFFile'] = 'muSelectionAndRecoSF_HZZ_Moriond17'

    typeToPrint = 'Preliminary'
    if paper:
        typeToPrint = ''

    style = _Style()

    outdir = _path.join('/afs/cern.ch/user/n/nawoods/www/UWVVPlots', plotDir)
    if not _exists(outdir):
        _mkdir(outdir)
    elif not _isdir(outdir):
        raise IOError("There is already some non-directory object called {}.".format(outdir))

    if paper:
        pdfDir = _path.join(outdir,'pdfs')
        if not _exists(pdfDir):
            _mkdir(pdfDir)
        elif not _isdir(pdfDir):
            raise IOError("There is already some non-directory object called {}.".format(pdfDir))
        texDir = _path.join(outdir,'texs')
        if not _exists(texDir):
            _mkdir(texDir)
        elif not _isdir(texDir):
            raise IOError("There is already some non-directory object called {}.".format(texDir))

    data, stack = standardZZSamples('zz', inData, inMC, ana, puWeightFile,
                                    fakeRateFile, lumi, amcatnlo=amcatnlo,
                                    higgs=(ana=='full'), eras=eras,
                                    skipEWK=(ana!='smp'), sipCut=sipForBkg,
                                    **sfArgs)

    objNames = _objNames.copy()
    if ana == 'smp':
        objNames['zz'] = r'\text{ZZ}'

    if doSyst:
        sig = zzStackSignalOnly('zz', inMC, ana, puWeightFile, lumi,
                                amcatnlo=amcatnlo, higgs=(ana=='full'),
                                asGroup=True, skipEWK=(ana!='smp'),
                                **sfArgs)
        irr = zzIrreducibleBkg('zz', inMC, ana, puWeightFile, lumi,
                               **sfArgs)

        sigSyst = {}
        irrSyst = {}
        for syst in ['eScaleUp', 'eScaleDn', 'eRhoResUp',
                     'eRhoResDn', 'ePhiResUp']:
            sigSyst[syst] = zzStackSignalOnly('eeee,eemm',
                                              inMC.replace('mc_',
                                                           'mc_{}_'.format(syst)),
                                              ana, puWeightFile, lumi,
                                              amcatnlo=amcatnlo, asGroup=True,
                                              higgs=(ana=='full'),
                                              skipEWK=(ana!='smp'), **sfArgs)
            irrSyst[syst] = zzIrreducibleBkg('eeee,eemm',
                                             inMC.replace('mc_',
                                                          'mc_{}_'.format(syst)),
                                             ana, puWeightFile, lumi, **sfArgs)

        for syst in ['mClosureUp','mClosureDn']:
            sigSyst[syst] = zzStackSignalOnly('eemm,mmmm',
                                              inMC.replace('mc_',
                                                           'mc_{}_'.format(syst)),
                                              ana, puWeightFile, lumi,
                                              amcatnlo=amcatnlo, asGroup=True,
                                              higgs=(ana=='full'),
                                              skipEWK=(ana!='smp'), **sfArgs)
            irrSyst[syst] = zzIrreducibleBkg('eemm,mmmm', inMC.replace('mc_','mc_{}_'.format(syst)),
                                             ana, puWeightFile, lumi, **sfArgs)


        bkg = standardZZBkg('zz', inData, inMC, ana, puWeightFile,
                            fakeRateFile, lumi, eras=eras, sipCut=sipForBkg,
                            **sfArgs)
        bkgSyst = {
            'eup' : standardZZBkg('zz', inData, inMC, ana, puWeightFile,
                                  fakeRateFile, lumi, eFakeRateSyst='up',
                                  eras=eras, sipCut=sipForBkg,
                                  **sfArgs),
            'edn' : standardZZBkg('zz', inData, inMC, ana, puWeightFile,
                                  fakeRateFile, lumi, eFakeRateSyst='dn',
                                  eras=eras, sipCut=sipForBkg,
                                  **sfArgs),
            'mup' : standardZZBkg('zz', inData, inMC, ana, puWeightFile,
                                  fakeRateFile, lumi, mFakeRateSyst='up',
                                  eras=eras, sipCut=sipForBkg,
                                  **sfArgs),
            'mdn' : standardZZBkg('zz', inData, inMC, ana, puWeightFile,
                                  fakeRateFile, lumi, mFakeRateSyst='dn',
                                  eras=eras, sipCut=sipForBkg,
                                  **sfArgs),
            }

    if ana == 'smp':
        aTGCDataset = 'uwvvNtuples_mc_21feb2017_aTGC'
        if looseSIP:
            aTGCDataset += '_LooseSIP'
        aTGCf4 = standardZZMC('zz', aTGCDataset,
                              'ZZTo4L-aTGC-f4-fg0p0019-fz0p0015',
                              ana, puWeightFile, lumi)
        aTGCf5 = standardZZMC('zz', aTGCDataset,
                              'ZZTo4L-aTGC-f5-fg0p0019-fz0p0015',
                              ana, puWeightFile, lumi)
        sherpa = standardZZMC('zz', aTGCDataset,
                              'ZZTo4L-sherpa',
                              ana, puWeightFile, lumi)

    binning4l = _binning4l.copy()
    binNormWidth4l = _binNormWidth4l.copy()
    units = _units.copy()

    if ana == 'z4l':
        binning4l['Mass'] = [20, 80., 100.]
        binNormWidth4l['Mass'] = 1.
    elif ana == 'full':
        binning4l['Mass'] = [25.*i for i in range(17)] + [500.,600.,800.]
        binNormWidth4l['Mass'] = 25.
    elif ana == 'smp':
        binning4l['Mass'] = [b / 1000. for b in binning4l['Mass']]
        binNormWidth4l['Mass'] /= 1000.

    for chan in ['zz', 'eeee', 'eemm', 'mmmm']:
        if paper and chan != 'zz':
            continue

        for varName, binning in binning4l.iteritems():
            if paper and varName != 'Mass':
                continue

            print "Plotting {} {}".format(chan, varName)

            var = _vars4l[varName]
            if varName == 'Mass' and ana == 'smp':
                var = {c:v+'/1000.' for c,v in var.iteritems()}

            if chan != 'zz':
                var = {chan:var[chan]}

            sel = _selections4l[varName]

            norm = binNormWidth4l[varName]
            normForYAxis = norm
            if varName == 'Mass' and ana == 'smp':
                normForYAxis = int(norm*1000)

            # blinding
            dataSelection = ''
            if blind and varName == 'Mass':
                dataSelection = 'Mass < 500.'
                if sel:
                    dataSelection += ' && ' + sel

            hStack = stack.makeHist(var, sel, binning,
                                    postprocess=True,
                                    perUnitWidth=norm)
            dataPts = data.makeHist(var, dataSelection, binning,
                                    poissonErrors=True,
                                    perUnitWidth=norm)
            if paper:
                dataPts.title = r'\textbf{{{}}}'.format(dataPts.title)
            toPlot = [hStack, dataPts]

            if doSyst:

                hSystUp, hSystDn = _makeSystematics(varName, var, sel, binning,
                                                    sig, bkg, irr,
                                                    sigSyst, bkgSyst, irrSyst,
                                                    norm, puWeightFile,
                                                    sfArgs)

                band = makeErrorBand(hStack, hSystUp, hSystDn)
                band.title = r'\textbf{Systematic unc.}'

                toPlot = [hStack, band, dataPts]

            c = Canvas(1000,1000)

            legParams = {
                'textsize':0.035,
                'leftmargin' : 0.45,
                }
            if ana == 'z4l' and varName == 'Mass' or ana == 'smp' and varName == 'deltaRZZ':
                legParams = _legParamsLeft.copy()
            leg = makeLegend(c, *toPlot, **legParams)

            xTitle = _xTitles[varName]
            if 'obj' in xTitle:
                xTitle = xTitle.format(obj=objNames[chan])
            yTitle = 'Events'
            if norm:
                yTitle += ' / {} {}'.format(makeNumberPretty(normForYAxis, 2),
                                            units[varName])

            if varName == 'Mass' and ana == 'smp':
                xTitle = xTitle.replace('GeV','TeV')
                if norm:
                    yTitle = yTitle.replace(makeNumberPretty(normForYAxis, 2),
                                            makeNumberPretty(normForYAxis/1000.,4)).replace('GeV','TeV')

            # cure inexplicable crash with inexplicable fix
            if chan == 'eeee' and varName == 'deltaPhiZZ':
                cTemp = Canvas(1000,1000)
                cTemp2 = Canvas(1000,1000)
                c.cd()

            drawOpts = {
                'xtitle' : xTitle,
                'ytitle' : yTitle,
                'logy' : logy,
                }
            # xmin, xmax, ymin, ymax = get_limits(toPlot, logy=logy)
            # print "rootpy min: ", ymin
            # print "{} minimum: {}".format(hStack[0].GetTitle(), min(b.value for b in hStack[0].bins()))
            # ymin = min(ymin, max(0.8*min(b.value for b in hStack[1].bins()),
            #                      1e-5 if logy else 0))
            # print "minimum chosen: {}".format(ymin)
            # drawOpts['ylimits'] = (ymin,ymax)

            if 'nJets' in varName:
                if logy:
                    c.SetLogy()
                frame = hStack.hists[0].empty_clone()
                frame.xaxis.title = xTitle
                frame.Draw() # for stack axis creation
                #toPlot = toPlot[1:]
                drawOpts['same'] = True
                for i in xrange(frame.GetXaxis().GetNbins()):
                    frame.GetXaxis().SetBinLabel(i+1, str(i))
                drawOpts['yaxis'] = frame.yaxis
            if varName in _nDivisions4l[ana]:
                drawOpts['xdivisions'] = _nDivisions4l[ana][varName]

            (xaxis, yaxis), (xmin,xmax,ymin,ymax) = draw(toPlot, c, **drawOpts)

            # blinding box
            if blind and varName == 'Mass' and binning4l['Mass'][-1] > 500.:
                box = TBox(max(xmin,500.), ymin, min(binning4l['Mass'][-1], xmax), ymax)
                box.SetFillColor(1)
                box.SetFillStyle(3002)
                box.Draw("same")
                leg.SetFillStyle(1001)

            leg.Draw("same")

            style.setCMSStyle(c, '', dataType=typeToPrint, intLumi=lumi,
                              forLatex=paper)
            if paper:
                _pdfViaTex(c, chan+varName, texDir, pdfDir)
            else:
                c.Print('{}/{}{}.png'.format(outdir, chan, varName))


            if varName == 'Mass' and ana == 'smp' and chan == 'zz':

                hTGCf4 = aTGCf4.makeHist(var, _selections4l[varName], binning,
                                         perUnitWidth=norm,
                                         mergeOverflow=True)
                hTGCf5 = aTGCf5.makeHist(var, _selections4l[varName], binning,
                                         perUnitWidth=norm,
                                         mergeOverflow=True)
                hSherpa = sherpa.makeHist(var, _selections4l[varName], binning,
                                          perUnitWidth=norm,
                                          mergeOverflow=True)

                for s in stack:
                    if s.name not in ['ZZTo4L', 'ZZTo4L-amcatnlo']:
                        hToAdd = s.makeHist(var, _selections4l[varName], binning,
                                            postprocess=True,
                                            perUnitWidth=norm,
                                            mergeOverflow=True)
                        hTGCf4 += hToAdd
                        hTGCf5 += hToAdd
                        hSherpa += hToAdd

                hTGCf4.SetLineWidth(2*hTGCf4.GetLineWidth())
                hTGCf5.SetLineWidth(2*hTGCf5.GetLineWidth())

                hSherpa.fillstyle = 'hollow'
                hSherpa.linecolor = 'black'
                hSherpa.SetLineWidth(2*hSherpa.GetLineWidth())
                #hSherpa.linestyle = 'dashed'

                hSherpa.legendstyle = 'L'
                hTGCf4.legendstyle = 'L'
                hTGCf5.legendstyle = 'L'

                toPlot = [hStack, hSherpa, hTGCf4, hTGCf5, dataPts]

                legParams['entryheight'] = 0.032
                legParams['entrysep'] = 0.008
                legParams['textsize'] = 0.031
                legParams['topmargin'] = 0.04
                legParams['rightmargin'] = 0.04

                cTGC = Canvas(1000,1000)
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

                style.setCMSStyle(cTGC, '', dataType=typeToPrint, intLumi=lumi,
                                  forLatex=paper)
                if paper:
                    _pdfViaTex(cTGC, chan+varName+'_aTGC', texDir, pdfDir)
                else:
                    cTGC.Print('{}/{}{}_aTGC.png'.format(outdir, chan, varName))




    binning2l = _binning2l.copy()
    binNormWidth2l = _binNormWidth2l.copy()

    if ana != 'smp':
        binning2l['Mass'] = [60, 0., 120.]
        binNormWidth2l['Mass'] = 2.


    for z in ['z', 'ze', 'zm', 'z1', 'z2']:
        for varName, binning in binning2l.iteritems():
            if paper:
                if varName != 'Mass':
                    continue
                if ana == 'z4l':
                    continue
                if ana == 'smp' and z != 'z':
                    continue
                if ana == 'full' and z not in ['z','z1']:
                    continue

            print "Plotting {} {}".format(z, varName)

            var = {c:[vt.format(var=varName) for vt in _varTemplates2l[z][c]] for c in _varTemplates2l[z]}

            hStack = stack.makeHist(var, _selections2l[z], binning,
                                    postprocess=True,
                                    perUnitWidth=binNormWidth2l[varName])
            dataPts = data.makeHist(var, _selections2l[z], binning,
                                    poissonErrors=True,
                                    perUnitWidth=binNormWidth2l[varName])

            if paper:
                dataPts.title = r'\textbf{{{}}}'.format(dataPts.title)

            if varName == 'Pt' and binning2l['Pt'][-1] > 200.:
                copy = dataPts.clone()
                x = copy.GetX()

                idx = 0
                for i in range(copy.GetN()):
                    if x[i] >= 200.:
                        dataPts.RemovePoint(idx)
                        continue
                    idx += 1

            toPlot = [hStack, dataPts]

            if doSyst:

                hSystUp, hSystDn = _makeSystematics(varName, var,
                                                    _selections2l[z], binning,
                                                    sig, bkg, irr,
                                                    sigSyst, bkgSyst, irrSyst,
                                                    binNormWidth2l[varName],
                                                    puWeightFile, sfArgs)

                band = makeErrorBand(hStack, hSystUp, hSystDn)
                band.title = r'\textbf{Systematic unc.}'

                toPlot = [hStack, band, dataPts]

            c = Canvas(1000,1000)

            legParams = {'textsize':0.03}
            if ana in ('full','smp') and varName == 'Mass':
                legParams = _legParamsLeft.copy()
            leg = makeLegend(c, *toPlot, **legParams)

            xTitle = _xTitles[varName]
            if '{obj}' in xTitle:
                xTitle = xTitle.format(obj=objNames[z])

            yTitle = 'Z bosons / {} {}'.format(makeNumberPretty(binNormWidth2l[varName],2),
                                               units[varName])

            drawArgs = {
                'xtitle' : xTitle,
                'ytitle' : yTitle,
                'logy' : logy,
                }
            if varName in _nDivisions2l[ana]:
                drawArgs['xdivisions'] = _nDivisions2l[ana][varName]

            (xaxis, yaxis), (xmin,xmax,ymin,ymax) = draw(toPlot, c, **drawArgs)

            # blinding box
            if blind and varName == 'Pt' and binning2l['Pt'][-1] > 200.:
                box = TBox(max(xmin,200.), ymin, min(binning2l['Pt'][-1], xmax), ymax)
                box.SetFillColor(1)
                box.SetFillStyle(3002)
                box.Draw("same")
                leg.SetFillStyle(1001)

            leg.Draw("same")

            style.setCMSStyle(c, '', dataType=typeToPrint, intLumi=lumi,
                              forLatex=paper)
            if paper:
                _pdfViaTex(c, z+varName, texDir, pdfDir)
            else:
                c.Print('{}/{}{}.png'.format(outdir, z, varName))


    for varName, binning in _binning1l.iteritems():
        if paper:
            continue
        for lep in _vars1l[varName]:
            print "Plotting {} {}".format(lep, varName)

            var = _vars1l[varName][lep]

            hStack = stack.makeHist(var, _selections1l[lep], _binning1l[varName],
                                    postprocess=True,
                                    perUnitWidth=_binNormWidth1l[varName])
            dataPts = data.makeHist(var, _selections1l[lep], _binning1l[varName],
                                    poissonErrors=True,
                                    perUnitWidth=_binNormWidth1l[varName])
            if paper:
                dataPts.title = r'\textbf{{{}}}'.format(dataPts.title)

            toPlot = [hStack, dataPts]

            if doSyst:

                hSystUp, hSystDn = _makeSystematics(varName, var,
                                                    _selections1l[lep],
                                                    binning, sig, bkg, irr,
                                                    sigSyst, bkgSyst, irrSyst,
                                                    _binNormWidth1l[varName],
                                                    puWeightFile, sfArgs)

                band = makeErrorBand(hStack, hSystUp, hSystDn)
                band.title = r'\textbf{Systematic unc.}'

                toPlot = [hStack, band, dataPts]

            c = Canvas(1000,1000)

            leg = makeLegend(c, *toPlot, leftmargin=0.47)

            xTitle = _xTitles[varName]
            if '{obj}' in xTitle:
                xTitle = xTitle.format(obj=objNames[lep])

            yTitle = 'Leptons / {} {}'.format(makeNumberPretty(_binNormWidth1l[varName],2),
                                              units[varName])

            (xaxis, yaxis), (xmin,xmax,ymin,ymax) = draw(toPlot, c,
                                                         xtitle=xTitle,
                                                         ytitle=yTitle,
                                                         logy=logy)
            leg.Draw("same")

            style.setCMSStyle(c, '', dataType=typeToPrint, intLumi=lumi,
                              forLatex=paper)
            if paper:
                _pdfViaTex(c, lep+varName, texDir, pdfDir)
            else:
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
    parser.add_argument('--blind', action='store_true',
                        help='Put blinding boxes on a few distributions.')
    parser.add_argument('--noSyst', action='store_true',
                        help="Don't make a hatched band for systematic errors.")
    parser.add_argument('--logy', '--logY', '--log', action='store_true',
                        help="Put the vertical axis on a log scale")
    parser.add_argument('--looseSIP', action='store_true',
                        help='Use scale factors for SIP<10 with no extra IP cuts.')
    parser.add_argument('--noSIP', action='store_true',
                        help='Use scale factors for no SIP cut and no extra IP cuts.')
    parser.add_argument('--sfRemake', action='store_true',
                        help='Use homebrewed scale factors for electrons.')
    parser.add_argument('--paper', action='store_true',
                        help='Print as journal-style PDFs, only the plots for the paper.')

    args=parser.parse_args()

    main(args.dataDir, args.mcDir, args.plotDir, args.analysis,
         args.fakeRateFile, args.puWeightFile, args.lumi, args.eras,
         args.blind, args.amcatnlo, not args.noSyst, args.logy, args.looseSIP,
         args.noSIP, args.sfRemake, args.paper)
