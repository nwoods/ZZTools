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

sampleInfo['TTJets-DiLept'] = {
    'xsec' : 56.7,
    'shortName' : 'TTJets',
    'prettyName' : 't \\={t}\\text{ + Jets}',
    'isSignal' : False,
    'format' : {
        'color' : 'limegreen',
        },
    }

sampleInfo['TTTo2L2Nu'] = {
    'xsec' : 87.31,
    'shortName' : 'TTJets',
    'prettyName' : 't \\={t}\\text{ + Jets}',
    'isSignal' : False,
    'format' : {
        'color' : 'limegreen',
        },
    }

sampleInfo['TTZ'] = {
    'xsec' : 0.2529,
    'prettyName' : 't#bar{t}Z',
    'isSignal' : False,
    'format' : {
        'color' : 'pink',
        },
}

sampleInfo['DYJets'] = {
    'xsec' : 5765.4, #5943.2, #6104.,
    'shortName' : 'DYJets',
    'prettyName' : '\\text{Z + Jets}',
    'isSignal' : False,
    'format' : {
        'color' : '#669966',
        },
    }

sampleInfo['DYJets-madgraph'] = {
    'xsec' : 5765.4, #5943.2, #6104.,
    'shortName' : 'DYJets',
    'prettyName' : '\\text{Z + Jets}',
    'isSignal' : False,
    'format' : {
        'color' : '#669966',
        },
    }

sampleInfo['DYToLL-0J'] = {
    'xsec' : 4759.,
    'shortName' : 'DY0J',
    'prettyName' : 'Z + 0 jets',
    'isSignal' : False,
    'format' : {
        'color' : '#669966',
        },
    }

sampleInfo['DYToLL-1J'] = {
    'xsec' : 881.,
    'shortName' : 'DY1J',
    'prettyName' : 'Z + 1 jet',
    'isSignal' : False,
    'format' : {
        'color' : '#669944',
        },
    }

sampleInfo['DYToLL-2J'] = {
    'xsec' : 336.3,
    'shortName' : 'DY2J',
    'prettyName' : 'Z + 2 jets',
    'isSignal' : False,
    'format' : {
        'color' : '#669922',
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

sampleInfo["ZZTo4L-amcatnlo"] = {
    'xsec' : 1.212,
    'prettyName' : 'q#bar{q} #rightarrow ZZ, Z#gamma*',
    'isSignal' : True,
    'format' : {
        'fillcolor' : '#99ccff',
        'linecolor' : '#000099',
        },
    'kFactor' : '1.1',
    }

sampleInfo["ZZTo4L-sherpa"] = {
    'xsec' : 0.455985,
    'prettyName' : 'q#bar{q} #rightarrow ZZ (SHERPA)',
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

sampleInfo["WZTo3LNu-amcatnlo"] = {
    'xsec' : 4.71,
    'prettyName' : '\\text{WZ + Jets}',
    'isSignal' : False,
    'format' : {
        'color' : 'violet',
        },
    }

sampleInfo["ZZJJTo4L_EWK"] = {
    'xsec' : 0.0004404,
    'prettyName' : 'ZZ + 2 jets EWK',
    'isSignal' : True,
    'format' : {
        'fillcolor' : '#bbeeff',
        'linecolor' : '#000099',
        },
    }

sampleInfo["ZZJJTo4L_QCD"] = {
    'xsec' : 0.009335,
    'prettyName' : 'ZZ + 2 jets QCD',
    'isSignal' : True,
    'format' : {
        'fillcolor' : '#bbeeff',
        'linecolor' : '#000099',
        },
    }

sampleInfo['ggHZZ'] = {
    'xsec' : 0.01212,
    'prettyName' : 'H #rightarrow ZZ',
    'isSignal' : True,
    'format' : {
        'fillcolor' : '#ffb2b2',
        'linecolor' : '#cc0000',
        },
    }

sampleInfo['WWZ'] = {
    'xsec' : 0.1651,
    'prettyName' : 'WWZ',
    'isSignal' : False,
    'format' : {
        'color' : 'purple',
        },
    }

sampleInfo['WZZ'] = {
    'xsec' : 0.5565,
    'prettyName' : 'WZZ',
    'isSignal' : False,
    'format' : {
        'color' : 'indigo',
        },
    }

sampleInfo['WZZ'] = {
    'xsec' : 0.5565,
    'prettyName' : 'WZZ',
    'isSignal' : False,
    'format' : {
        'color' : 'indigo',
        },
    }

sampleInfo['ZZZ'] = {
    'xsec' : 0.01398,
    'prettyName' : 'ZZZ',
    'isSignal' : False,
    'format' : {
        'color' : 'darkorchid',
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

groupInfo['irreducible'] = {
    'isSignal' : False,
    'format' : {
        'fillcolor' : '#ffcc00',
        'linecolor' : 'black',
        },
    'prettyName' : 't#bar{t}Z, WWZ',
    }

groupInfo['DYToLL'] = {
    'shortName' : 'DYJets',
    'prettyName' : '\\text{Z + Jets}',
    'isSignal' : False,
    'format' : {
        'fillcolor' : '#669966',
        'linecolor' : '#003300',
        },
    }
