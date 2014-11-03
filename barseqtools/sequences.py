from Bio import SeqIO
from collections import defaultdict

def filetype(args):
    if args.fasta:
        filename = args.fasta
        fastx = 'fasta'
    elif args.fastq:
        filename = args.fastq
        fastx = 'fastq'
    else:
        return None, None
    return filename, fastx

def subseq_counter(filename, fastx, start, length):
    '''
    filename = string, path to fq/fa file
    fastx = "fasta" or "fastq"
    start = start coordinate in all reads
    length = length of subsequence to look at
    '''
    end = start+length
    counts = defaultdict(int)
    for fx in SeqIO.parse(filename, fastx):
        sequence = str(fx.seq)
        counts[sequence[start:end]] += 1
    return counts

def sample_index_counter(filename, fastx):
    return subseq_counter(filename, fastx, 0,5)

def yeast_barcode_counter(filename, fastx):
    return subseq_counter(filename, fastx, 20,20)

def combo_counter(filename, fastx, s1, l1, s2, l2):
    '''
    filename = string, path to fq/fa file
    fastx = "fasta" or "fastq"
    start = start coordinate in all reads
    length = length of subsequence to look at
    '''
    end1 = s1+l1
    end2 = s2+l2
    counts = defaultdict(dict)
    for fx in SeqIO.parse(filename, fastx):
        sequence = str(fx.seq)
        try:
            counts[sequence[s1:end1]][sequence[s2:end2]] += 1
        except KeyError:
            counts[sequence[s1:end1]][sequence[s2:end2]] = 1
    return counts

def index_barcode_counter(filename, fastx):
    return combo_counter(filename, fastx, s1=0, l1=5, s2=20, l2=20)




### it seems silly to force this...qc check...
def triple_checker(filename, fastx, s1, l1, s2, l2, qcStart, qcSeq):
    '''
    filename = string, path to fq/fa file
    fastx = "fasta" or "fastq"
    start = start coordinate in all reads
    length = length of subsequence to look at
    s1,l1 = start and length of seq1 search
    s1,s2 - start and length of seq2 search
    qcStart = start position of qc sequence in reads
    qcSeq = qc sequence
    1. gets seq1
    2. checks quality of qc seq (hamming dist)
    3. seq2
    '''
    end1 = s1+l1
    end2 = s2+l2
    counts = defaultdict(dict)
    for fx in SeqIO.parse(filename, fastx):
        sequence = str(fx.seq)
        try:
            counts[sequence[s1:end1]][sequence[s2:end2]] += 1
        except KeyError:
            counts[sequence[s1:end1]][sequence[s2:end2]] = 1
    return counts
    




