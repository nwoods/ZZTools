'''

Little script to plot Senka's 2D limit curves in the same style as the rest of
the paper.

Nate Woods, U. Wisconsin

'''

import logging
from rootpy import log as rlog; rlog = rlog["/aTGCLimits2D"]
# don't show most silly ROOT messages
logging.basicConfig(level=logging.WARNING)
rlog["/ROOT.TUnixSystem.SetDisplay"].setLevel(rlog.ERROR)

from rootpy import asrootpy
from rootpy.plotting import Canvas, Hist
from rootpy.io import root_open
from rootpy.ROOT import gStyle, gPad

from PlotTools import PlotStyle, pdfViaTex, makeLegend

from os.path import join as pjoin, isdir, exists
from os import makedirs as mkdirp
from collections import OrderedDict


indir = '/data/nawoods/aTGCLimits2D'
fTemplate = 'output_contours_f{}.root'
fNames = {str(nf):pjoin(indir, fTemplate.format(nf)) for nf in (4,5)}

outdir = '/afs/cern.ch/user/n/nawoods/www/aTGCLimits2D'
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

style = PlotStyle()
gStyle.SetGridColor(1)

# names in the root files
graphNames = OrderedDict([
    ('exp68'   , 'contour_68exp'),
    ('exp95'   , 'contour_95exp'),
    ('exp99'   , 'contour_99exp'),
    ('obs'     , 'contour_obs'),
    ('bestFit' , 'bestFit'),
    ])

titles = {
    'exp68'   : r'\textbf{Expected 68\% C.L.}',
    'exp95'   : r'\textbf{Expected 95\% C.L.}',
    'exp99'   : r'\textbf{Expected 99\% C.L.}',
    'obs'     : r'\textbf{Observed 95\% C.L.}',
    'bestFit' : r'\textbf{Best Fit}',
    }

colors = {
    'exp68'   : 'b',
    'exp95'   : 'g',
    'exp99'   : 'r',
    'obs'     : 'black',
    'bestFit' : 'black',
    }

linestyle = {
    'exp68'   : 'verylongdash',
    'exp95'   : 'verylongdash',
    'exp99'   : 'verylongdash',
    'obs'     : 'solid',
    'bestFit' : '',
    }

for nf, fName in fNames.iteritems():
    graphs = OrderedDict()
    with root_open(fName) as f:
        for name, nameInFile in graphNames.iteritems():
            g = asrootpy(getattr(f,nameInFile))
            g.title = titles[name]
            if name == 'bestFit':
                g.legendstyle = 'P'
                g.SetMarkerStyle(20)
                g.SetMarkerSize(2)
                g.drawstyle = 'P'
            else:
                g.drawstyle = 'L'
                g.legendstyle = 'L'
            g.color = colors[name]
            if linestyle[name]:
                g.linestyle = linestyle[name]
                g.SetLineWidth(2*g.GetLineWidth())
            graphs[name] = g

    c = Canvas(1000,1000)

    c.SetGrid(1,1)

    axlimits = (-0.0045,0.0045)

    # bug in this version of rootpy, so draw ourselves instead of using util
    frame = Hist(1,*axlimits)
    frame.Draw()
    frame.SetLineWidth(0)
    xaxis = frame.xaxis
    yaxis = frame.yaxis

    for g in graphs.values():
        g.Draw('SAME')

    xaxis.SetLimits(*axlimits)
    xaxis.SetRangeUser(*axlimits)
    xaxis.SetTitle(r'\boldsymbol{{f}}_\mathbf{{{}}}^\mathbf{{\gamma}}'.format(nf))
    yaxis.SetLimits(*axlimits)
    yaxis.SetRangeUser(*axlimits)
    yaxis.SetTitle(r'\boldsymbol{{f}}_\mathbf{{{}}}^\textbf{{Z}}'.format(nf))

    xaxis.SetNoExponent(True)
    yaxis.SetNoExponent(True)
    xaxis.CenterTitle()
    yaxis.CenterTitle()
    xaxis.SetNdivisions(405)
    yaxis.SetNdivisions(405)

    leg = makeLegend(c, *graphs.values(),
                     leftmargin=0.045,
                     rightmargin=0.455,
                     topmargin=0.63,
                     textsize=0.03)

    leg.SetFillStyle(1001)
    leg.Draw("same")

    c.Print(pjoin(outdir,'limits2D_f{}.png'.format(nf)))
    style.setCMSStyle(c, '', True, '', intLumi=35860., forLatex=True)

    pdfViaTex(c, 'limits2D_f{}'.format(nf), texdir, pdfdir)
