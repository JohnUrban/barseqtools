#
from sequences import *
from useful_functions import *
from collections import defaultdict
from copy import deepcopy

def run(parser, args):
    filename, fastx = filetype(args)
    if filename and fastx:
        counts = sample_index_counter(filename, fastx)
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
        print "Percent of number of anticipated indexes: %d" % percent
        print "Number of 5mers that make up anticipated indexes: %d" % total
        percent = 100*float(total)/sum(counts.values())
        print "Percent of 5mers that make up anticipated indexes: %d" % percent
        f.close()
        print "Hamming dist info: minimum hamming distance for a 5mer against all anticipated indexes"
        print "..and Unique assignment info: how many could be uniquely assigned with a hamming dist of 0, 1, 2 - when choosing the match that minimizes the hamming distance?"
        print "...in other words, even if allowing a distance of 1 causes formerly uniquely assignable kmers with hd=0 to match more kmers with hd=1, they are still uniquely assignable to the lowest hd=0 match"
        hd = defaultdict(int)
        hdunique = defaultdict(int)
        kmerhd = dict()
        for kmer in counts:
            kmerhd[kmer] = minHammingDist(kmer, indexes)
            minkmers, minhd = kmerhd[kmer]
            hd[minhd] += counts[kmer]
            if len(minkmers) == 1:
                hdunique[minhd] += counts[kmer]
        print ("\t").join(["minHammingDist", "5-merCount", "5-merPercent", "unquelyAssignableCount", "unquelyAssignablePercentKmersWithSameHamDist", "unquelyAssignablePercentTotalKmers", "cumulativeUnquelyAssignablePercentTotalKmers"])
        total = float(sum(counts.values()))
        cumUniTotal = 0
        for ham in sorted(hd.keys()):
            uniTotal = 100*float(hdunique[ham])/total
            cumUniTotal += uniTotal
            print ("\t").join(["%d", "%d", "%f", "%d", "%f", "%f", "%f"]) % (ham, hd[ham], 100*float(hd[ham])/total, hdunique[ham], 100*float(hdunique[ham])/hd[ham], uniTotal, cumUniTotal)
        ## counts with 0, 1, 2 per sample
        print "Uniquely assignable counts seen per sample when allowing hamming distances of:"
        print ("\t").join(["%s", "%d", "%d", "%d"]) % ("sample",0,1,2)
        totalWhenhd2 = 0
        for sample in sorted(samples.keys()):
            index = samples[sample]
            hd0 = 0
            hd1 = 0
            hd2 = 0
            for kmer in counts:
                if len(kmerhd[kmer][0]) == 1 and kmerhd[kmer][0][0] == index:
                    ## IF uniquely assignable to the current index
                    if kmerhd[kmer][1] == 0:
                        hd0 += counts[kmer]
                    elif kmerhd[kmer][1] == 1:
                        hd1 += counts[kmer]
                    elif kmerhd[kmer][1] == 2:
                        hd2 += counts[kmer]                             
            print ('\t').join(["%s", "%d", "%d", "%d", "%d"]) % (sample, hd0, hd0, hd0+hd1, hd0+hd1+hd2)
            totalWhenhd2 += hd0+hd1+hd2
##        print 100*totalWhenhd2/total - cumUniTotal
        

        
        
