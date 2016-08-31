sampleInfo = {}

sampleInfo['TTJets'] = {
    'xsec' : 670.3,
    'shortName' : 'TTJets',
    'prettyName' : 't \\={t}\\text{ + Jets}',
    'isSignal' : False,
    'format' : {
        'color' : 'limegreen',
        },
    }

sampleInfo['DYJets'] = {
    'xsec' : 6104.,
    'shortName' : 'DYJets',
    'prettyName' : '\\text{Z + Jets}',
    'isSignal' : False,
    'format' : {
        'color' : '#669966',
        },
    }

sampleInfo['GluGluZZTo4e'] = {
    'xsec' : 0.001586,
    'prettyName' : '\\text{gg}\\!\\!\\rightarrow \\!\\! \\text{ZZ}\\!\\!\\rightarrow \\!\\! 4e',
    'isSignal' : True,
    'kFactor' : '1.7',
    'format' : {
        'fillcolor' : '#4b78ff',
        'linecolor' : '#000099',
        },
    }

sampleInfo['GluGluZZTo4mu'] = {
    'xsec' : 0.001586,
    'prettyName' : '\\text{gg}\\!\\!\\rightarrow \\!\\! \\text{ZZ}\\!\\!\\rightarrow \\!\\! 4\\mu',
    'isSignal' : True,
    'kFactor' : '1.7',
    'format' : {
        'fillcolor' : '#4b78ff',
        'linecolor' : '#000099',
        },
    }

sampleInfo['GluGluZZTo4tau'] = {
    'xsec' : 0.001586,
    'prettyName' : '\\text{gg}\\!\\!\\rightarrow \\!\\! \\text{ZZ}\\!\\!\\rightarrow \\!\\! 4\\tau',
    'isSignal' : False,
    'kFactor' : '1.7',
    'format' : {
        'color' : 'darkgray',
        },
    }

sampleInfo['GluGluZZTo2e2mu'] = {
    'xsec' : 0.00319364,
    'prettyName' : '\\text{gg}\\!\\!\\rightarrow \\!\\! \\text{ZZ}\\!\\!\\rightarrow \\!\\! 2e2\\mu',
    'isSignal' : True,
    'kFactor' : '1.7',
    'format' : {
        'fillcolor' : '#4b78ff',
        'linecolor' : '#000099',
        },
    }

sampleInfo['GluGluZZTo2e2tau'] = {
    'xsec' : 0.00319364,
    'prettyName' : '\\text{gg}\\!\\!\\rightarrow \\!\\! \\text{ZZ}\\!\\!\\rightarrow \\!\\! 2e2\\tau',
    'isSignal' : False,
    'kFactor' : '1.7',
    'format' : {
        'color' : 'darkgray',
        },
    }

sampleInfo['GluGluZZTo2mu2tau'] = {
    'xsec' : 0.00319364,
    'prettyName' : '\\text{gg}\\!\\!\\rightarrow \\!\\! \\text{ZZ}\\!\\!\\rightarrow \\!\\! 2\\mu 2\\tau',
    'isSignal' : False,
    'kFactor' : '1.7',
    'format' : {
        'color' : 'darkgray',
        },
    }

sampleInfo["ZZTo4L"] = {
    'xsec' : 1.256,
    'prettyName' : 'q#bar{q} #rightarrow ZZ, Z#gamma*', 
    'isSignal' : True,
    'format' : {
        'fillcolor' : '#99ccff', 
        'linecolor' : '#000099',
        },
    'kFactor' : '1.1', 
    }

sampleInfo["WZTo3LNu"] = {
    'xsec' : 4.42965,
    'prettyName' : '\\text{WZ + Jets}',
    'isSignal' : False,
    'format' : {
        'color' : 'violet',
        },
    }



groupInfo = {}

groupInfo['GluGluZZ'] = {
    'isSignal' : True,
    'prettyName' : 'gg #rightarrow ZZ, Z#gamma*', 
    'format' : {
        'fillcolor' : '#4b78ff', 
        'linecolor' : '#000099',
        },
    }

groupInfo['ggZZ2l2tau'] = {
    'isSignal' : False,
    'prettyName' : 'gg \\rightarrow ZZ, Z\\gamma* \\rightarrow 2\\ell 2\\tau', 
    'format' : {
        'color' : 'darkgray', 
        },
    }

groupInfo['Z+X'] = {
    'isSignal' : False,
    'format' : {
        'fillcolor' : '#669966',
        'linecolor' : '#003300',
        },
    }

