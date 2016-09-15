

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
nameMap.add("TTJets", "TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8")
nameMap.add("WZTo3LNu", "WZTo3LNu_TuneCUETP8M1_13TeV-powheg-pythia8")
nameMap.add("ZZTo4L-amcatnlo", "ZZTo4L_13TeV-amcatnloFXFX-pythia8")
nameMap.add("ZZTo4L", "ZZTo4L_13TeV_powheg_pythia8")
nameMap.add("ggHZZ", "GluGluHToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8")
nameMap.add("ggHZZ-minloHJJ", "GluGluHToZZTo4L_M125_13TeV_powheg2_minloHJJ_JHUgenV6_pythia8")

for fs in ['4e', '4mu', '4tau', '2e2mu', '2e2tau', '2mu2tau']:
    nameMap.add("GluGluZZTo{}".format(fs), 
                "GluGluToContinToZZTo{}_13TeV_MCFM701_pythia8".format(fs))

for width in ['0', '0p1', '0p2', '0p3']:
    for mass in [750, 800, 1200, 2000, 3000, 4000]:
        nameMap.add("Grav2PB_width{}_M{}".format(width, mass), 
                    "Graviton2PBToZZTo4L_width{}_M-{}_13TeV-JHUgenV6-pythia8".format(width, mass))
