
from os import path as _path, system as _unix, makedirs as _mkdirp
from shutil import move as _mv
from re import compile as _reComp


_texTemplate = '''
\\documentclass[tikz]{{standalone}}
\\usepackage{{standalone}}
\\usetikzlibrary{{patterns}}
\\usetikzlibrary{{plotmarks}}
\\usepackage{{amsmath}}
\\begin{{document}}
  \\input{{{fname}}}
\\end{{document}}
'''

def pdfViaTex(c, fname, texDir, pdfDir):
    '''
    Print a Canvas as a PDF, via a ROOT-generated .tex file.

    c (Canvas): Canvas to print.
    fname (str): Files will be called fname.tex and fname.pdf.
    texDir (str): Directory for tex files and pdflatex output. Will be created
        if necessary.
    pdfDir (str): Directory for final PDF. Will be created if necessary.
    '''
    if not _path.exists(texDir):
        _mkdirp(texDir)

    imgFile = _path.join(texDir, fname+'_img.tex')

    c.Print(imgFile)

    if not _path.exists(imgFile):
        raise IOError("Something went wrong trying to print {} to a tex file.".format(fname))

    # Remove unwanted boxes from around hatched and transparent fill areas
    imgFileFixed = imgFile.replace('.tex','_fixed.tex')
    drawPattern = _reComp(r'\\draw(?= \[((pattern=)|(.+fill opacity=)))')
    # make transparency actually work for hatched areas
    # there's probably a way to combine with the previous regex...
    opacityPattern = _reComp(r'(?<=\\path \[pattern=crosshatch, pattern color=c, )fill (?=opacity=[01])')
    with open(imgFile, 'r') as fIm:
        with open(imgFileFixed, 'w') as fImFix:
            for line in fIm:
                fImFix.write(opacityPattern.sub('',drawPattern.sub(r'\path',line)))

    texFile = _path.join(texDir, fname+'.tex')

    with open(texFile, 'w') as f:
        f.write(_texTemplate.format(fname=imgFileFixed))

    _unix('pdflatex -halt-on-error -output-directory {} {}'.format(texDir, texFile))

    pdfFile = texFile.replace('.tex','.pdf')

    if not _path.exists(pdfFile):
        raise IOError("Something went wrong trying to make {} from {}.".format(pdfFile, texFile))

    if not _path.exists(pdfDir):
        _mkdirp(pdfDir)

    newPDFFile = _path.join(pdfDir, fname+'.pdf')
    _mv(pdfFile, newPDFFile)


# for purists
pdfViaTeX = pdfViaTex
