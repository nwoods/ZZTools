from SampleTools import MCSample
from Utilities import zMassDist, mapObjects, parseChannels
from Analysis.weightHelpers import puWeight

from os.path import join

def getEventInfo(row, *args):
    return {
        'run' : row.run,
        'event' : row.evt,
        'lumi' : row.lumi,
        }


def getCandInfo(row, fPUWt, *objects):
    numbers = getEventInfo(row)
    numbers['m4l'] = row.Mass
    numbers['mZ1'] = getattr(row, '_'.join([objects[0], objects[1], 'Mass']))
    numbers['mZ2'] = getattr(row, '_'.join([objects[2], objects[3], 'Mass']))

    # eemm channel may have masses swapped
    if zMassDist(numbers['mZ1']) > zMassDist(numbers['mZ2']):
        numbers['mZ1'], numbers['mZ2'] = numbers['mZ2'], numbers['mZ1']

    jetPts = row.jetPt
    numbers['nJets'] = jetPts.size()
    if jetPts.size():
        numbers['jet1pt'] = jetPts.at(0)
        if jetPts.size() > 1:
            numbers['jet2pt'] = jetPts.at(1)
        else:
            numbers['jet2pt'] = -1.
    else:
        numbers['jet1pt'] = -1.
        numbers['jet2pt'] = -1.

    numbers['mjj'] = max(-1.,row.mjj)

    numbers['weight'] = fPUWt(row.nTruePU) * reduce(lambda x,y : x*y,
                                                    (getattr(row, obj+'EffScaleFactor') for obj in objects), 1.)

    return numbers


def getCandInfo3l(row, *objects):
    numbers = getEventInfo(row)
    numbers['m3l'] = row.Mass
    numbers['mZ'] = getattr(row,  '_'.join([objects[0], objects[1], 'Mass']))
    numbers['ptL3'] = getattr(row, objects[2]+'Pt')
    numbers['l3Tight'] = 1 if getattr(row, objects[2]+'ZZTightID') and getattr(row, objects[2]+'ZZIsoPass') else 0

    return numbers


def getAllInfo(channel, sample, fInfo):
    found = set()
    objects = mapObjects(channel)
    if channel == 'emm':
        objects = objects[1:]+objects[:1]

        numbers = fInfo(row, *objects)
        evtID = (numbers['run'],numbers['lumi'],numbers['event'])
        if evtID in found:
            continue
        found.add(evtID)
        yield numbers


if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Dump information about the 4l candidates in an ntuple to a text file, for synchronization.')
    parser.add_argument('input', type=str, nargs=1,
                        help='Glob for files holding ntuples, either absolute or relative to /data/nawoods/ntuples')
    parser.add_argument('output', type=str, nargs='?', default='candSync.txt',
                        help='Name of the text file to output.')
    parser.add_argument('channels', nargs='?', type=str, default='zz',
                        help='Comma separated (no spaces) list of channels, or keyword indicated multiple channels')
    parser.add_argument('--listOnly', action='store_true',
                        help='Print only run:lumi:event with no further info')
    parser.add_argument('--zPlusL', '--zPlusl', action='store_true',
                        help='Use the Z+l control region format instead of the 4l format')
    parser.add_argument('--puWeightFile', type=str, nargs='?',
                        default='puWeight_69200_24jan2017.root',
                        help=('Name of pileup weight file (assumed to be in '
                              'usual  data directory unless full path is '
                              'specified)'))


    args = parser.parse_args()

    inputPath = join('/data/nawoods/ntuples', args.input[0])

    if args.zPlusL and args.channels == 'zz':
        args.channels = 'zl'
    channels = parseChannels(args.channels)

    outStrings = []

    if args.listOnly:
        outTemp = '{run}:{lumi}:{event}:{channel}\n'
        infoGetter = getEventInfo
    elif args.zPlusL:
        outTemp = ('{run}:{lumi}:{event}:{channel}:{m3l:.2f}:{mZ:.2f}:{ptL3:.2f}:'
                   '{l3Tight}\n')
        infoGetter = getCandInfo3l

    else:
        outTemp = ('{run}:{lumi}:{event}:{channel}:{m4l:.2f}:{mZ1:.2f}:{mZ2:.2f}:'
                   '{nJets}:{jet1pt:.2f}:{jet2pt:.2f}:{mjj:.2f}:{weight:.4f}\n')
        puWtStr, puWtFun = puWeight(args.puWeightFile)
        infoGetter = lambda row, *objects: getCandInfo(row, puWtFun, *objects)


    for channel in channels:
        if channel == 'emm':
            channelForStr = 'mme' # for sync with Torino
        else:
            channelForStr = channel

        sample = MCSample('SyncSample', channel, inputPath)

        for numbers in getAllInfo(channel, sample, infoGetter):
            outStrings.append(outTemp.format(channel=channelForStr, **numbers))

    with open(args.output, 'w') as fout:
        for s in sorted(outStrings, key=lambda x: [int(y) for y in x.split(':')[:3]]):
            fout.write(s)


