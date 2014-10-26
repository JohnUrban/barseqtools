from sequences import *
from classes import *

def run(parser, args):
    filename, fastx = filetype(args)
    data = BarSeq(fastx, filename, args.indexes, args.barcodes, args.outprefix)
    data.analyze()

