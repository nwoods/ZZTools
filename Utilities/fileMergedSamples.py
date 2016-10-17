
from argparse import ArgumentParser as _Args
from glob import glob as _glob
from os.path import join as _join
import re as _re

from Utilities.moveAndRename import moveAndRename
from Metadata.nameMap import nameMap

if __name__ == '__main__':
    parser = _Args(description=("Move groups of files, renaming them in a "
                                "consistent way."))

    parser.add_argument('id', type=str, nargs=1,
                        help='ID of samples, e.g. "08SEP2016_0"')
    parser.add_argument('newID', type=str, nargs='?',
                        help='New ID for samples in ntuple directory, if different')

    args = parser.parse_args()


    oldID = args.id[0]
    if not args.newID:
        newID = oldID.lower()
    else:
        newID = args.newID

    paths = _glob('/data/nawoods/MERGE_*{}*mergeFilesJob'.format(oldID))

    pattern = 'MERGE_(?P<ntupleType>[A-Za-z0-9]+)_(?P<dataset>\w+)_{}_(?P<sample>\S+)-mergeFilesJob'.format(oldID)

    groups = {}
    for p in paths:
        m = _re.match(pattern, p.split('/')[-1])
        sampleInfo = m.groupdict()

        if sampleInfo['ntupleType'] not in groups:
            groups[sampleInfo['ntupleType']] = {}

        if 'DATA' in sampleInfo['dataset']:
            sampleName = 'Run2016{}_{}'.format(sampleInfo['dataset'][-1], sampleInfo['sample'])
        else:
            sampleName = nameMap[sampleInfo['sample']]

        groups[sampleInfo['ntupleType']][sampleName] = _join(p, '*', '*.root')

    for nt in groups:
        if 'singlez' in nt.lower():
            ntupleType = 'uwvvSingleZ'
        elif 'zplusl' in nt.lower():
            ntupleType = 'uwvvZPlusl'
        else:
            ntupleType = 'uwvvNtuples'

        dataGroups = {}
        mcGroups = {}
        for sample in groups[nt]:
            if 'Run2016' in sample:
                dataGroups[sample] = groups[nt][sample]
            else:
                mcGroups[sample] = groups[nt][sample]

        if len(dataGroups):
            moveAndRename(_join('/data/nawoods/ntuples', '{}_data_{}'.format(ntupleType, newID)),
                          **dataGroups)
        if len(mcGroups):
            moveAndRename(_join('/data/nawoods/ntuples', '{}_mc_{}'.format(ntupleType, newID)),
                          **mcGroups)

        if not (len(dataGroups) or len(mcGroups)):
            print "No samples found for group {}, nothing filed".format(nt)

