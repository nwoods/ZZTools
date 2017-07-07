'''

Little script to plot Senka's 1D limit curves as a function of mass cutoff
in the same style as the rest of the paper.

Nate Woods, U. Wisconsin

'''

import logging
from rootpy import log as rlog; rlog = rlog["/aTGCLimits1D"]
# don't show most silly ROOT messages
logging.basicConfig(level=logging.WARNING)
rlog["/ROOT.TUnixSystem.SetDisplay"].setLevel(rlog.ERROR)

from rootpy import asrootpy
from rootpy.io import root_open
from rootpy.ROOT import gStyle

from PlotTools import PlotStyle, pdfViaTex
from PlotTools import CMS_lumi
CMS_lumi.lumiTextSize = 1.2
CMS_lumi.cmsTextSize = 1.3
CMS_lumi.lumiTextOffset = 0.01
CMS_lumi.cmsTextOffset = 0.01

from os.path import join as pjoin, isdir, exists
from os import makedirs as mkdirp


indir = '/data/nawoods/aTGCLimits1D'
fTemplate = 'aTGCCutoff_canvas_f{}{}.root'#'plot_expObs_Cutoff_f{}{}_noPrelim_May10_fit_.root'
fNames = {vf+str(nf):pjoin(indir, fTemplate.format(nf,vf)) for nf in (4,5) for vf in 'gz'}

outdir = '/afs/cern.ch/user/n/nawoods/www/aTGCLimits1D_paper'
texdir = pjoin(outdir, 'texs')
pdfdir = pjoin(outdir, 'pdfs')

if not exists(texdir):
    mkdirp(texdir)
elif not isdir(texdir):
    raise IOError("There is already some non-directory object called {}.".format(texdir))
if not exists(pdfdir):
    mkdirp(pdfdir)
elif not isdir(pdfdir):
    raise IOError("There is already some non-directory object called {}.".format(pdfdir))

style = PlotStyle(False)
gStyle.SetOptFit(0)
gStyle.SetPadColor(0)
gStyle.SetOptTitle(0)

vName = {
    'g' : r'\\gamma',
    'z' : r'\\text{Z}',
    }

for fType, fName in fNames.iteritems():
    with root_open(fName) as f:
        c = asrootpy(f.c1_par_pol2bs_fit)#c1_par_expobs_fit)

        sub = {
            # totally redo y-axis title
            r'$f^{}_{}\\ 95\\% \\ \\text{{CL}}$'.format(vName[fType[0]],fType[1]):r'\} 95\\%CL',
            # make x-axis title bold
            r'm_{4\\ell} \\ \\text{cutoff (GeV)}' : r'm_\{4\\ell\}\\ \\text\{cut-off \(GeV\)\}',
            # center infinity
            r'\\!\\!\\infty' : r'\\infty',
            }

        style.setCMSStyle(c, '', True, '', intLumi=35860., forLatex=True)
        c.Draw()
        c.Print(pjoin(outdir,'limits1DVsCutoff_{}.png'.format(fType)))
        pdfViaTex(c, 'limits1DVsCutoff_{}'.format(fType), texdir, pdfdir, **sub)

