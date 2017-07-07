'''

Make a scatter plot of m_Z2 vs m_Z1, skipping enough events that it's readable.

Author: Nate Woods, U. Wisconsin

'''

import logging
from rootpy import log as rlog; rlog = rlog["/zMassScatter"]
# don't show most silly ROOT messages
logging.basicConfig(level=logging.WARNING)
rlog["/ROOT.TUnixSystem.SetDisplay"].setLevel(rlog.ERROR)

from rootpy.plotting import Graph, Canvas, Legend
from rootpy.plotting.utils import draw

from PlotTools import PlotStyle as _Style, pdfViaTex
from Utilities import Z_MASS
from Analysis import standardZZData

from os.path import exists, isdir, join
from os import makedirs as mkdir


inDir = 'uwvvNtuples_data_10mar2017_LooseSIPLooseVtx'
outDir = '/afs/cern.ch/user/n/nawoods/www/UWVVPlots/zMassScatter_paper'
lumi = 35860.

# plot one event for every thinFactor events found
thinFactor = 3

colors = {'eeee':'b','eemm':'r','mmmm':'forestgreen'}
markers = {'eeee':20,'eemm':21,'mmmm':22}
titles = {'eeee':r'\mathbf{4e}','eemm':r'\mathbf{2e2}\boldsymbol{\mu}','mmmm':r'\boldsymbol{4\mu}'}

texDir = join(outDir, 'texs')
pdfDir = join(outDir, 'pdfs')
if not exists(texDir):
    mkdir(texDir)
elif not isdir(texDir):
    raise IOError("There is already some non-directory object called {}.".format(texDir))
if not exists(pdfDir):
    mkdir(pdfDir)
elif not isdir(pdfDir):
    raise IOError("There is already some non-directory object called {}.".format(pdfDir))

def getMZ1_2e2m(row):
    if abs(getattr(row, 'e1_e2_Mass') - Z_MASS) < abs(getattr(row, 'm1_m2_Mass') - Z_MASS):
        return getattr(row, 'e1_e2_Mass')
    return getattr(row, 'm1_m2_Mass')
def getMZ2_2e2m(row):
    if abs(getattr(row, 'e1_e2_Mass') - Z_MASS) < abs(getattr(row, 'm1_m2_Mass') - Z_MASS):
        return getattr(row, 'm1_m2_Mass')
    return getattr(row, 'e1_e2_Mass')

getMZ1 = {
    'eeee' : lambda row: getattr(row, 'e1_e2_Mass'),
    'eemm' : getMZ1_2e2m,
    'mmmm' : lambda row: getattr(row, 'm1_m2_Mass'),
}
getMZ2 = {
    'eeee' : lambda row: getattr(row, 'e3_e4_Mass'),
    'eemm' : getMZ2_2e2m,
    'mmmm' : lambda row: getattr(row, 'm3_m4_Mass'),
}

style = _Style()

legParams = {
    'full' : {
        'leftmargin' : 0.1,
        'rightmargin' : 0.55,
        },
    'z4l' : {
        'leftmargin' : 0.6,
        },
    }

xlimits = {
    'z4l' : (40.,90.),
    'full' : (40.,120.),
    }
ylimits = {
    'z4l' : (0.,60.),
    'full' : (0.,120.),
    }


data = standardZZData('zz', inDir, 'full')

g = {'full':{},'z4l':{}}
for ch, sample in data.itersamples():
    nPtsMax = int(sample.getEntries())
    for ana in g:
        g[ana][ch] = Graph(nPtsMax, title=titles[ch])
        g[ana][ch].color = colors[ch]
        g[ana][ch].markerstyle = markers[ch]
        g[ana][ch].drawstyle = 'P'
        g[ana][ch].SetMarkerSize(g[ana][ch].GetMarkerSize()*1.5)
        if ch == 'mmmm':
            g[ana][ch].SetMarkerSize(g[ana][ch].GetMarkerSize()*1.18)

    iFull = 0
    iZ4l = 0
    for i, row in enumerate(sample.rows()):
        plotZ4l =  80. < row.Mass < 100.
        plotFull = i % thinFactor == 0
        if not (plotZ4l or plotFull):
            continue

        mZ1 = getMZ1[ch](row)
        mZ2 = getMZ2[ch](row)

        if plotZ4l:
            g['z4l'][ch].SetPoint(iZ4l, mZ1, mZ2)
            iZ4l += 1
        if plotFull:
            g['full'][ch].SetPoint(iFull, mZ1, mZ2)
            iFull += 1

    g['z4l'][ch].Set(iZ4l)
    g['full'][ch].Set(iFull)

for ana in ['full','z4l']:#g:
    # proper order
    graphList = [g[ana][ch] for ch in ['eemm','mmmm','eeee']]

    c = Canvas(1000,1000)
    (xaxis,yaxis), things = draw(graphList, c,
                                 xtitle=r'm_{\text{Z}_1} \, \left(\text{GeV}\right)',
                                 ytitle=r'm_{\text{Z}_2} \, \left(\text{GeV}\right)',
                                 xlimits=xlimits[ana], ylimits=ylimits[ana])
    yaxis.SetTitleSize(yaxis.GetTitleSize()*0.9)
    c.Update()
    leg = Legend(4, c, textsize=0.04,
                 #header=r'\quad\textbf{Data}',
                 entrysep=0.01,
                 entryheight=0.04, **legParams[ana])
    leg.AddEntry(None, r'\textbf{Data}', '')
    for gr in graphList:
        leg.AddEntry(gr)
    leg.Draw("same")
    style.setCMSStyle(c, "", True, "", 13, lumi, forLatex=True)

    pdfViaTex(c, 'mZ2VsmZ1_{}'.format(ana), texDir, pdfDir)
    #c.Print(join(outDir, 'mZ2VsmZ1_{}.png'.format(ana)))

