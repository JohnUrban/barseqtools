#
from sequences import *
from collections import defaultdict
from useful_functions import *



def run(parser, args):
    filename, fastx = filetype(args)
    start = args.start
    length = len(args.sequence)
    if start+length > args.min_read_length:
        print "Start position and sequence length must be made to not exceed shortest read"
        print "Default assumes shortest read is 50 bp"
        print "Right now, your analysis is used %d bp" % args.min_read_length
        print "Adjust with --min-read-length __"
        quit()
    
    if filename and fastx:
        counts = subseq_counter(filename, fastx, start, length)
        print "Number unique %d-mers seen at position %d: %d" % (length, args.start, len(counts))
        print "Number %d-mers encountered: %d" % (length, sum(counts.values()))
        represented = counts[args.sequence]
        print "Number of anticipated sequences found at given start position: %d" % represented
        total = float(sum(counts.values()))
        percent = 100*float(represented)/total
        print "Percent of %d-mers that make up anticipated sequence: %d" % (length, percent)
        print "Hamming distance information:"
        hd = defaultdict(int)
        for kmer in counts.keys():
            hd[hammingDist(kmer, args.sequence)] += counts[kmer]
        print ("\t").join(["HammingDist", "%d-merCount", "%d-merPercent"]) % (length, length)
        for ham, count in sorted(hd.iteritems()):
            print ("\t").join(["%d", "%d", "%f"]) % (ham, count, 100*float(count)/total)

        
