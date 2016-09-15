'''

Simple script to take groups of files and move them all to the same place, giving them standardized names

Nate Woods, U. Wisconsin

'''


from argparse import ArgumentParser as _Args
from glob import glob as _glob
from shutil import move as _mv
from shutil import copy as _cp
from os import makedirs as _mkdir
from os.path import join as _join


def moveAndRename(outdir, extension='.root', copy=False, **samples):
    try:
        _mkdir(outdir)
    except OSError:
        pass
    
    if copy:
        move = _cp
    else:
        move = _mv
    
    for name, filePath in samples.iteritems():
        files = _glob(filePath)
    
        for i, f in enumerate(files):
            newName = _join(outdir, '{}_{}{}'.format(name, i, extension))
            move(f, newName)
    

if __name__ == '__main__':
    parser = _Args(description=("Move groups of files, renaming them in a "
                                "consistent way."))
    
    parser.add_argument('outdir', type=str, nargs=1,
                        help=('Directory to place files (will be created '
                              'if needed).'))
    parser.add_argument('--cp', '--copy', action='store_true',
                        help='Copy the files instead of moving them.')
    parser.add_argument('--extension', type=str, nargs='?', default='.root',
                        help=('Ending for output files (in case they are not '
                              'root files).'))
    parser.add_argument('-g', nargs=2, type=str, action='append',
                        help=('(once per group) '
                              'NAME GLOB pair for one file group. All files '
                              'matching GLOB will be moved and renamed '
                              'NAME_n.root, where n is a counter so that each file '
                              'has a unique name.'))
    
    args = parser.parse_args()
    
    
    if not args.g:
        raise IOError("No files specified with -g option")
    
    groups = dict(args.g)

    moveAndRename(args.outdir[0], args.extension, args.cp, **groups)
