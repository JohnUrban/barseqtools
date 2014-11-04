# Nov 3: TODO -- 
from sequences import *
from useful_functions import *
from collections import defaultdict



def run(parser, args):
    filename, fastx = filetype(args)
    if filename and fastx:
        counts = yeast_barcode_counter(filename, fastx)
        print "Number unique 20mer barcodes: %d" % len(counts)
        print "Number 20mers encountered: %d" % sum(counts.values())
        f = open(args.barcodes)
        barcodes = set()
        ORFs = dict()
        ORFlist = list() ## delete: ORFlist is just for diagnosing a problem where there are redundant ORF:barcode pairs in file and repreated ORFs with diff barcodes
        for line in f:
            if line[0] != "#":
                line = line.strip().split()
                barcodes.add(line[1])
                ORFs[line[0]] = line[1]
                if line[0] in ORFlist:
                    print line[0]
                ORFlist.append(line[0])
        print len(barcodes), len(ORFs.values()), len(set(ORFs.values())), len(ORFlist), len(set(ORFlist))   
        quit()
        print "Number of barcodes anticipated: %d" % len(barcodes)
        represented = 0 # Find how many of the anticipated barcodes were encountered via exact match
        total = 0       # Find the number of 20mers that account for anticipated barcodes
        for e in barcodes:
            if e in counts:
                represented += 1
                total += counts[e]
        print "Number of anticipated barcodes found: %d" % represented
        percent = 100*float(represented)/float(len(barcodes))
        print "Percent of number of antipated barcodes: %d" % percent
        print "Number of 20mers that make up anticipated barcodes: %d" % total
        percent = 100*float(total)/sum(counts.values())
        print "Percent of 20mers that make up anticipated barcodes: %d" % percent
        f.close()
        if args.hamming_analysis:
            print "Hamming dist info: minimum hamming distance for a 20mer against all anticipated barcodes"
            print "..and Unique assignment info: how many 20mers could be uniquely assigned with a hamming dist of 0, 1, 2 - when choosing the match that minimizes the hamming distance?"
            print "...in other words, even if allowing a distance of 1 causes formerly uniquely assignable kmers with hd=0 to match more kmers with hd=1, they are still uniquely assignable to the lowest hd=0 match"
            hd = defaultdict(int)
            hdunique = defaultdict(int)
            kmerhd = dict()
            for kmer in counts:
                kmerhd[kmer] = minHammingDist(kmer, barcodes)
                minkmers, minhd = kmerhd[kmer]
                hd[minhd] += counts[kmer]
                if len(minkmers) == 1:
                    hdunique[minhd] += counts[kmer]
    #
            print ("\t").join(["minHammingDist", "20-merCount", "20-merPercent", "unquelyAssignableCount", "unquelyAssignablePercentKmersWithSameHamDist", "unquelyAssignablePercentTotalKmers", "cumulativeUnquelyAssignablePercentTotalKmers"])
            total = float(sum(counts.values()))
            cumUniTotal = 0
            for ham in sorted(hd.keys()): ## hamming distances from 0-n
                uniTotal = 100*float(hdunique[ham])/total
                cumUniTotal += uniTotal
                print ("\t").join(["%d", "%d", "%f", "%d", "%f", "%f", "%f"]) % (ham, hd[ham], 100*float(hd[ham])/total, hdunique[ham], 100*float(hdunique[ham])/hd[ham], uniTotal, cumUniTotal)
            ## counts with 0, 1, 2 per sample
            print "Uniquely assignable counts seen per ORF when allowing hamming distances of:"
            print ("\t").join(["%s", "%d", "%d", "%d"]) % ("ORF",0,1,2)
            totalWhenhd0 = 0
            totalWhenhd1 = 0
            totalWhenhd2 = 0
            represented = defaultdict(int)
            for ORF in sorted(ORFs.keys()):
                barcode = ORFs[ORF]
                hd0 = 0
                hd1 = 0
                hd2 = 0
                for kmer in counts:
                    if len(kmerhd[kmer][0]) == 1 and kmerhd[kmer][0][0] == barcode: ## IF uniquely assignable to the current index
                        if kmerhd[kmer][1] == 0:
                            hd0 += counts[kmer]
                        elif kmerhd[kmer][1] == 1:
                            hd1 += counts[kmer]
                        elif kmerhd[kmer][1] == 2:
                            hd2 += counts[kmer]
                if hd0+hd1+hd2 > 0:
                    print ('\t').join(["%s", "%d", "%d", "%d"]) % (ORF, hd0, hd0+hd1, hd0+hd1+hd2)
                    if hd0 == 0:
                        if hd1 > 0:
                            represented[1] += 1
                        elif hd2 > 0:
                            represented[2] += 1
                    else:
                        represented[0] += 1
                            
                totalWhenhd0 += hd0
                totalWhenhd1 += hd0+hd1
                totalWhenhd2 += hd0+hd1+hd2
            print "Total 20mers accepted when hd=0:", totalWhenhd0  ## TODO/FigOut: Why does this not agree with other counts (2326 vs 2333)
            print "Total 20mers accepted when hd=1:", totalWhenhd1
            print "Total 20mers accepted when hd=2:", totalWhenhd2
            print "Number of anticipated ORFs represented when HD = 0, 1, 2:"
            print ("\t").join(["HD", "numberNewORFs", "totalNumberORFs", "PercentAnticipated"])
            totalORFs = represented[0]
            print ("\t").join(["%d", "%d", "%d", "%f"]) % (0, represented[0], totalORFs, 100.0*totalORFs/len(barcodes))
            totalORFs += represented[1]
            print ("\t").join(["%d", "%d", "%d", "%f"]) % (1, represented[1], totalORFs, 100.0*totalORFs/len(barcodes))
            totalORFs += represented[2]
            print ("\t").join(["%d", "%d", "%d", "%f"]) % (2, represented[2], totalORFs, 100.0*totalORFs/len(barcodes))

            print len(barcodes), len(ORFs.values())
            print set(barcodes) == set(ORFs.values())
            
        else:
            print "Exact match counts seen per ORF"
            for ORF,barcode in ORFs.iteritems():
                if counts[barcode] > 0:
                    print ('\t').join([ORF, str(counts[barcode])])

    print "Diganosing why there is disagreement between 2 diff analyses"
    print "Need to adjust uptags file to be non-redundant and for repeated gene names with diff barcodes to have slightly diff names"
    print "I could also have the code 'sense' these repeated names and slightly rename them - and sense redundant pairs to ignore..."

    
