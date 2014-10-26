#
from sequences import *

def run(parser, args):
    filename, fastx = filetype(args)
    if filename and fastx:
        counts = yeast_barcode_counter(filename, fastx)
        print "Number unique 20mer indexes: %d" % len(counts)
        print "Number 20mers encountered: %d" % sum(counts.values())
        f = open(args.barcodes)
        indexes = set()
        ORFs = dict()
        for line in f:
            if line[0] != "#":
                line = line.strip().split()
                indexes.add(line[1])
                ORFs[line[0]] = line[1]
        print "Number of barcodes anticipated: %d" % len(indexes)
        represented = 0
        total = 0
        for e in indexes:
            if e in counts:
                represented += 1
                total += counts[e]
        print "Number of anticipated barcodes found: %d" % represented
        percent = 100*float(represented)/float(len(indexes))
        print "Percent of number of antipated barcodes: %d" % percent
        print "Number of 20mers that make up anticipated barcodes: %d" % total
        percent = 100*float(total)/sum(counts.values())
        print "Percent of 20mers that make up anticipated barcodes: %d" % percent
        print "Counts seen per ORF"
        for ORF,barcode in ORFs.iteritems():
            if counts[barcode] > 0:
                print ('\t').join([ORF, str(counts[barcode])])
        
