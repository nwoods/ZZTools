
from os import path as _path, system as _unix, makedirs as _mkdirp
from shutil import move as _mv
from re import compile as _reComp


_texTemplate = '''
\\documentclass[tikz]{{standalone}}
\\usepackage{{standalone}}
\\usetikzlibrary{{patterns}}
\\usetikzlibrary{{plotmarks}}
\\usepackage{{amsmath}}
\\usepackage{{lmodern}}
\\usepackage{{sansmath}}
\\sansmath
\\usepackage[scaled]{{helvet}}
\\usepackage{{lmodern}}
\\usepackage{{sansmath}}
\\sansmath
\\usepackage[scaled]{{helvet}}
\\renewcommand\\familydefault{{\\sfdefault}}
\\DeclareSymbolFont{{numbers}}{{T1}}{{phv}}{{m}}{{n}}
\\DeclareMathSymbol{{0}}{{0}}{{numbers}}{{"30}}
\\DeclareMathSymbol{{1}}{{0}}{{numbers}}{{"31}}
\\DeclareMathSymbol{{2}}{{0}}{{numbers}}{{"32}}
\\DeclareMathSymbol{{3}}{{0}}{{numbers}}{{"33}}
\\DeclareMathSymbol{{4}}{{0}}{{numbers}}{{"34}}
\\DeclareMathSymbol{{5}}{{0}}{{numbers}}{{"35}}
\\DeclareMathSymbol{{6}}{{0}}{{numbers}}{{"36}}
\\DeclareMathSymbol{{7}}{{0}}{{numbers}}{{"37}}
\\DeclareMathSymbol{{8}}{{0}}{{numbers}}{{"38}}
\\DeclareMathSymbol{{9}}{{0}}{{numbers}}{{"39}}
\\usepackage[T1]{{fontenc}}
\\begin{{document}}
  \\input{{{fname}}}
\\end{{document}}
'''

def _doSub(s, (sub,exp)):
    '''
    Replace regex exp with str sub in str s.
    '''
    return exp.sub(sub,s)


def pdfViaTex(c, fname, texDir, pdfDir, **extraSubs):
    '''
    Print a Canvas as a PDF, via a ROOT-generated .tex file.

    c (Canvas): Canvas to print.
    fname (str): Files will be called fname.tex and fname.pdf.
    texDir (str): Directory for tex files and pdflatex output. Will be created
        if necessary.
    pdfDir (str): Directory for final PDF. Will be created if necessary.
    extraSubs(str keyed to str): Value is a regular expression that will be
        replaced with key anywhere it appears in the output tex file, via
        re.sub().
    '''
    if not _path.exists(texDir):
        _mkdirp(texDir)

    imgFile = _path.join(texDir, fname+'_img.tex')

    c.Print(imgFile)

    if not _path.exists(imgFile):
        raise IOError("Something went wrong trying to print {} to a tex file.".format(fname))

    subList = []

    # Remove unwanted boxes from around hatched and transparent fill areas
    imgFileFixed = imgFile.replace('.tex','_fixed.tex')
    subList.append(('\path',_reComp(r'\\draw(?= \[((pattern=)|(.+fill opacity=)))')))
    # make transparency actually work for hatched areas
    # there's probably a way to combine with the previous regex...
    subList.append((r'',_reComp(r'(?<=\\path \[pattern=crosshatch, pattern color=c, )fill (?=opacity=[01])')))
    # anything else that needs to change
    subList += [(k, _reComp(v)) for k,v in extraSubs.iteritems()]
    with open(imgFile, 'r') as fIm:
        with open(imgFileFixed, 'w') as fImFix:
            for line in fIm:
                fImFix.write(reduce(_doSub, subList, line))

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
