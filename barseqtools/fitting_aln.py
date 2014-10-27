#
from useful_functions import *
from Bio import SeqIO
from sequences import *

class FitAln(object):
    def __init__(self, fastx, fastx_file, query_seq, with_read_names=True):
        self.fastx = fastx
        self.fastx_file = fastx_file
        self.query_seq = query_seq
        self.query_seq_len = float(len(query_seq))
        self.with_read_names = with_read_names
    def fitting_aln(self, read):
        return generalizedFittingAln(read, self.query_seq)[:2]
    def print_fitting_aln_stats_over_all_reads(self):
        for fx in SeqIO.parse(self.fastx_file, self.fastx):
            readname = (fx.name+"\t")*self.with_read_names
            score, start = self.fitting_aln(str(fx.seq))
            print readname+("\t").join([str(start), str(score), str(score/self.query_seq_len)])



def run(parser, args):
    filename, fastx = filetype(args)
    if args.random_sequence:
        args.sequence = random_seq(args.random_sequence)
    fit = FitAln(fastx, filename, args.sequence, args.with_read_names)
    fit.print_fitting_aln_stats_over_all_reads()
    
