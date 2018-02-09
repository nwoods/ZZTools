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
from rootpy.plotting import Canvas, Hist, Graph
from rootpy.io import root_open
from rootpy.ROOT import gStyle, gPad

from PlotTools import PlotStyle, pdfViaTex, makeLegend

from os.path import join as pjoin, isdir, exists
from os import makedirs as mkdirp
from collections import OrderedDict


indir = '/data/nawoods/aTGCLimits2D'
fTemplate = 'output_contours_f{}.root'
fNames = {nf:pjoin(indir, fTemplate.format(nf)) for nf in (4,5)}

outdir = '/afs/cern.ch/user/n/nawoods/www/aTGCLimits2D_postFR'
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
    #('exp99'   , 'contour_99exp'),
    ('obs'     , 'contour_obs'),
    ('bestFit' , 'bestFit'),
    ])

titles = {
    'exp68'   : r'\textbf{Exp.\  68\% CL (2D)}',
    'exp95'   : r'\textbf{Exp.\  95\% CL (2D)}',
    'exp99'   : r'\textbf{Exp.\  99\% CL (2D)}',
    'obs'     : r'\textbf{Obs.\  95\% CL (2D)}',
    'bestFit' : r'\textbf{Best Fit (2D)}',
    }

colors = {
    'exp68'   : 'b',
    'exp95'   : 'magenta',
    'exp99'   : 'r',
    'obs'     : 'black',
    'bestFit' : 'black',
    }

linestyle = {
    'exp68'   : 'dashed',
    'exp95'   : 'verylongdash',
    'exp99'   : 'longdash',
    'obs'     : 'solid',
    'bestFit' : '',
    }

limits1D = {
    4 : {
        'z' : [-0.0012,0.0010],
        'g' : [-0.0012,0.0013],
        },
    5 : {
        'z' : [-0.0010,0.0013],
        'g' : [-0.0012,0.0013],
        },
    }

for nf, fName in fNames.iteritems():
    graphs_ = OrderedDict()

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
            graphs_[name] = g

    bars1D = Graph(1, type='asymm', title=r'\textbf{Obs.\ 95\% CL (1D)}',
                   color='r', drawstyle='ZL', legendstyle='L')
    bars1D.SetLineWidth(4 * bars1D.GetLineWidth())
    bars1D.SetPointError(0, *(abs(e) for e in limits1D[nf]['g']+limits1D[nf]['z']))
    bars1D.SetMarkerStyle(7)

    # have to reorder. Ugh.
    graphs = OrderedDict()
    for name, g in graphs_.iteritems():
        if name == 'bestFit':
            graphs['1D'] = bars1D
            continue
        graphs[name] = g
    graphs['bestFit'] = graphs_['bestFit']

    c = Canvas(1000,1000)

    c.SetGrid(1,1)

    axlimits = (-0.0028,0.0028) #(-0.0045,0.0045)

    # bug in this version of rootpy, so draw ourselves instead of using util
    frame = Hist(1,*axlimits)
    frame.Draw()
    frame.SetLineWidth(0)
    xaxis = frame.xaxis
    yaxis = frame.yaxis

    leg = makeLegend(c, *graphs.values(),
                     leftmargin=0.023,
                     rightmargin=0.452,
                     topmargin=0.63,
                     textsize=0.036)

    leg.SetFillStyle(1001)
    leg.Draw("same")

    for g in graphs.values():
        g.Draw('SAME')

    xaxis.SetLimits(*axlimits)
    xaxis.SetRangeUser(*axlimits)
    xaxis.SetTitle(r'f_{}^\gamma'.format(nf))
    yaxis.SetLimits(*axlimits)
    yaxis.SetRangeUser(*axlimits)
    yaxis.SetTitle(r'f_{}^\text{{Z}}'.format(nf)) #r'\mathbf{{f}}_\mathbf{{{}}}^\textbf{{Z}}'.format(nf))

    xaxis.SetNoExponent(True)
    yaxis.SetNoExponent(True)
    xaxis.CenterTitle()
    yaxis.CenterTitle()
    xaxis.SetNdivisions(405)
    yaxis.SetNdivisions(405)
    xaxis.SetLabelOffset(-0.1*xaxis.GetLabelOffset())
    yaxis.SetLabelOffset(-0.1*yaxis.GetLabelOffset())
    xaxis.SetLabelSize(1.1*xaxis.GetLabelSize())
    yaxis.SetLabelSize(1.1*yaxis.GetLabelSize())
    xaxis.SetTitleOffset(0.9*xaxis.GetTitleOffset())
    yaxis.SetTitleOffset(0.9*yaxis.GetTitleOffset())
    xaxis.SetTitleSize(1.1*xaxis.GetTitleSize())
    yaxis.SetTitleSize(1.1*yaxis.GetTitleSize())
    nDiv = xaxis.GetNdivisions()
    nDiv += nDiv%100
    xaxis.SetNdivisions(nDiv)
    yaxis.SetNdivisions(nDiv)

    c.Print(pjoin(outdir,'limits2D_f{}.png'.format(nf)))
    style.setCMSStyle(c, '', True, '', intLumi=35860., forLatex=True)

    pdfViaTex(c, 'limits2D_f{}'.format(nf), texdir, pdfdir)
