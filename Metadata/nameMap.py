

class _TwoWayMap(dict):
    '''
    A poor man's bidirectional map.
    '''
    def __setitem__(self, k, v):
        try:
            del self[k]
            del self[v]
        except KeyError:
            pass

        super(_TwoWayMap, self).__setitem__(k, v)
        super(_TwoWayMap, self).__setitem__(v, k)

    def __delitem__(self, k):
        super(_TwoWayMap, self).__delitem__(self[k])
        super(_TwoWayMap, self).__delitem__(k)

    def __len__(self):
        return super(_TwoWayMap, self).__len__() // 2

    def add(self, k, v):
        self.__setitem__(k,v)


nameMap = _TwoWayMap()

nameMap.add("DYJets", "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8")
nameMap.add("DYJets-madgraph", "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
nameMap.add("DYToLL-0J", "DYToLL_0J_13TeV-amcatnloFXFX-pythia8")
nameMap.add("DYToLL-1J", "DYToLL_1J_13TeV-amcatnloFXFX-pythia8")
nameMap.add("DYToLL-2J", "DYToLL_2J_13TeV-amcatnloFXFX-pythia8")
nameMap.add("TTJets", "TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8")
nameMap.add("TTJets-TuneCUETP8M1", "TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8")
nameMap.add("TTJets-DiLept", "TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8")
nameMap.add("TTTo2L2Nu", "TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8")
nameMap.add("WZTo3LNu", "WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8")
nameMap.add("WZTo3LNu-amcatnlo", "WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8")
nameMap.add("ZZTo4L-amcatnlo", "ZZTo4L_13TeV-amcatnloFXFX-pythia8")
nameMap.add("ZZTo4L-sherpa", "ZZTo4L_13TeV-sherpa")
nameMap.add("ZZTo4L", "ZZTo4L_13TeV_powheg_pythia8")
nameMap.add("ggHZZ", "GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8")
nameMap.add("ggHZZ-minloHJJ", "GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJJ_JHUgenV6_pythia8")
nameMap.add("ZZJJTo4L_EWK", "ZZJJTo4L_EWK_13TeV-madgraph-pythia8")
nameMap.add("ZZJJTo4L_QCD", "ZZJJTo4L_QCD_13TeV-madgraph-pythia8")
nameMap.add("WWZ", "WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8")
nameMap.add("WZZ", "WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8")
nameMap.add("ZZZ", "ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8")
nameMap.add("TTZ", "TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8")

for fs in ['4e', '4mu', '4tau', '2e2mu', '2e2tau', '2mu2tau']:
    nameMap.add("GluGluZZTo{}".format(fs),
                "GluGluToContinToZZTo{}_13TeV_MCFM701_pythia8".format(fs))
nameMap.add("GluGluZZTo4mu-defaultShower", "GluGluToContinToZZTo4mu_13TeV_DefaultShower_MCFM701_pythia8")
nameMap.add("GluGluZZTo2e2mu-defaultShower", "GluGluToContinZZTo2e2mu_13TeV_DefaultShower_MCFM701_pythia8")

for width in ['0', '0p1', '0p2', '0p3']:
    for mass in [750, 800, 1200, 2000, 3000, 4000]:
        nameMap.add("Grav2PB_width{}_M{}".format(width, mass),
                    "Graviton2PBToZZTo4L_width{}_M-{}_13TeV-JHUgenV6-pythia8".format(width, mass))


systNameMap = _TwoWayMap()

systNameMap.add("ESCALEUP", "eScaleUp")
systNameMap.add("ESCALEDN", "eScaleDn")
systNameMap.add("ERHORESUP", "eRhoResUp")
systNameMap.add("ERHORESDN", "eRhoResDn")
systNameMap.add("EPHIRESUP", "ePhiResUp")
systNameMap.add("MCLOSUREUP", "mClosureUp")
systNameMap.add("MCLOSUREDN", "mClosureDn")


