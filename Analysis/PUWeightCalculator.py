'''

Uses a text file with the MC pileup distribution and histograms of the data
pileup distribution and makes a file with a histogram of the ratio.

Data distributions can be made in a CMS environment with a command like
$ pileupCalc.py -i dataLumi_DoubleMuon_23aug2016.json --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 50 --numPileupBins 50 ~/ZZTools/data/pileup/puDistData_69200_23aug2016.root

The json file passed as the first argument above can be created from farmout 
output with a command like
$ /data/cms/farmout/jobReportSummary.py --json-out lumisIUsed.json `find /path/to/submit/directories/*.xml`

Nate Woods, U. Wisconsin

'''


from rootpy.io import root_open
from rootpy.plotting import Hist
from rootpy import asrootpy
from os import path, environ


outputFileName = 'puWeight_69200_23aug2016.root'
mcDistFileName = 'mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU.txt'
dataDistFileTemplate = 'puDistData_{}_23aug2016.root'
dataDistFileName = dataDistFileTemplate.format(69200)
dataDistFileNameDn = dataDistFileTemplate.format(66017)
dataDistFileNameUp = dataDistFileTemplate.format(72383)
nbins = 50

# make a histogram from the MC distribution text file
puMC = Hist(50,0.,50.)
with open(path.join(environ['zzt'], 'data', 'pileup', mcDistFileName)) as fMC:
    for i, n in enumerate(fMC):
        puMC[i+1].value = float(n)

# Get the data distributions from root files
with root_open(path.join(environ['zzt'], 'data', 
                         'pileup', dataDistFileName)) as fData:
    puData = asrootpy(fData.pileup).clone()
    puData.SetDirectory(0)
with root_open(path.join(environ['zzt'], 'data', 
                         'pileup', dataDistFileNameDn)) as fDataDn:
    puDataDn = asrootpy(fDataDn.pileup).clone()
    puDataDn.SetDirectory(0)
with root_open(path.join(environ['zzt'], 'data', 
                         'pileup', dataDistFileNameUp)) as fDataUp:
    puDataUp = asrootpy(fDataUp.pileup).clone()
    puDataUp.SetDirectory(0)

# normalize them
puData.Scale(1./puData.Integral())
puDataDn.Scale(1./puDataDn.Integral())
puDataUp.Scale(1./puDataUp.Integral())

with root_open(path.join(environ['zzt'], 'data', 'pileup', outputFileName),
               'recreate') as f:
    out = puData.clone(name='puScaleFactor')
    out /= puMC
    out.write()

    outDn = puDataDn.clone(name='puScaleFactor_down')
    outDn /= puMC
    outDn.write()

    outUp = puDataUp.clone(name='puScaleFactor_up')
    outUp /= puMC
    outUp.write()
                         

