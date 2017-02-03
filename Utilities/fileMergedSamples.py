
from argparse import ArgumentParser as _Args
from glob import glob as _glob
from os.path import join as _join
from re import compile as _reComp

from Utilities.moveAndRename import moveAndRename as _mv
from Metadata.nameMap import nameMap as _nameMap, systNameMap as _sysNameMap

if __name__ == '__main__':
    parser = _Args(description=("Move groups of files, renaming them in a "
                                "consistent way."))

    parser.add_argument('id', type=str, nargs=1,
                        help='ID of samples, e.g. "08SEP2016_0"')
    parser.add_argument('newID', type=str, nargs='?',
                        help='New ID for samples in ntuple directory, if different')
    parser.add_argument('--dryRun', '--dry-run', '--dry_run',
                        action='store_true',
                        help='Print the files that would move, but don\'t actually move them.')
    parser.add_argument('--cp', '--copy', action='store_true',
                        help='Copy the files instead of moving them.')

    args = parser.parse_args()


    oldID = str.lstrip(args.id[0])
    if not args.newID:
        newID = oldID.lower()
    else:
        newID = args.newID

    paths = _glob('/data/nawoods/MERGE_*{}*mergeFilesJob'.format(oldID))

    #pattern = 'MERGE_(?P<ntupleType>[A-Za-z0-9]+)_(?P<dataset>\w+)_{}_(?P<sample>\S+?)(?P<ext>(_((backup)|(ext\\d)))?)?-mergeFilesJob'.format(oldID)
    dirPattern = _reComp('MERGE_(?P<ntupleType>[A-Za-z0-9]+)_(?P<inputType>\w+)_{}_(?P<sample>\S+?)-mergeFilesJob'.format(oldID))
    sysPattern = _reComp('MC_?(?P<sys>\\w*)')
    extPattern = _reComp('(?P<sample>\\S+?)_?(?P<ext>(backup|ext\\d)?(?=$))')

    groups = {}
    for p in paths:
        m = dirPattern.match(p.split('/')[-1])
        if m is None:
            raise IOError("I don't know how to interpret the path {}".format(p))
        sampleInfo = m.groupdict()

        if sampleInfo['ntupleType'] not in groups:
            groups[sampleInfo['ntupleType']] = {}

        if 'DATA' in sampleInfo['inputType']:
            inputType = 'data'
            sampleName = 'Run2016{}_{}'.format(sampleInfo['inputType'][-1], sampleInfo['sample'])
        else:
            mMC = sysPattern.match(sampleInfo['inputType'])
            if mMC is None:
                raise IOError("I don't know how to interpret the path {}".format(p))
            inputType = 'mc'
            sys = mMC.groupdict()['sys']
            if sys:
                if sys in _sysNameMap:
                    sys = _sysNameMap[sys]
                inputType += '_' + sys

            mSample = extPattern.match(sampleInfo['sample'])
            if mSample is None:
                raise IOError("I don't know how to interpret the path {}".format(p))
            try:
                sampleName = _nameMap[mSample.groupdict()['sample']]
            except KeyError:
                print extPattern
                print p.split('/')[-1]
                print sampleInfo
                raise

            ext = mSample.groupdict()['ext']
            if ext:
                sampleName += '_' + ext

        if inputType not in groups[sampleInfo['ntupleType']]:
            groups[sampleInfo['ntupleType']][inputType] = {}

        assert sampleName not in groups[sampleInfo['ntupleType']][inputType], \
            "Trying to move 2 {} {} samples both called {}.".format(sampleInfo['ntupleType'], inputType, sampleName)
        groups[sampleInfo['ntupleType']][inputType][sampleName] = _join(p, '*', '*.root')

    for nt in groups:
        if 'singlez' in nt.lower():
            ntupleType = 'uwvvSingleZ'
        elif 'zplusl' in nt.lower():
            ntupleType = 'uwvvZPlusl'
        else:
            ntupleType = 'uwvvNtuples'

        for it in groups[nt]:
            dirName = '_'.join([ntupleType, it, newID])

            if len(groups[nt][it]):
                _mv(_join('/data/nawoods/ntuples', dirName), copy=args.cp,
                    dryRun = args.dryRun, **groups[nt][it])


