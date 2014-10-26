#TODO!!!!!
from sequences import *

def run(parser, args):
    filename, fastx = filetype(args)
    if filename and fastx:
        counts = index_barcode_counter(filename, fastx)
        print "Number unique 5mer indexes: %d" % len(counts)
        print "Number 5mers encountered: %d" % sum(counts.values())
        f = open(args.indexes, 'r')
        indexes = set()
        samples = dict()
        for line in f:
            if line[0] != "#":
                line = line.strip().split()
                indexes.add(line[1])
                samples[line[0]] = line[1]
        print "Number of indexes anticipated: %d" % len(indexes)
        represented = 0
        total = 0
        for e in indexes:
            if e in counts:
                represented += 1
                total += counts[e]
        print "Number of anticipated indexes found: %d" % represented
        percent = 100*float(represented)/float(len(indexes))
        print "Percent of number of antipated indexes: %d" % percent
        print "Number of 5mers that make up anticipated indexes: %d" % total
        percent = 100*float(total)/sum(counts.values())
        print "Percent of 5mers that make up anticipated indexes: %d" % percent
        f.close()
        print "Counts seen per sample"
        for sample in sorted(samples.keys()):
            index = samples[sample]
            print ('\t').join([sample, str(counts[index])])
        
