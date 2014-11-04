#
from useful_functions import *
from Bio import SeqIO
from sequences import *
import sys

##TODO:
## make it so multiple_seqs can report all maximum fitting aln infos
## do this by storing all max fit aln info in  a dictionary in "best_fitting_aln"
## then in "print_best_..." loop over dict keys.
class FitAln(object):
    def __init__(self, fastx, fastx_file, query_seq=None, multiple_query_sequences=None, with_read_names=True, with_edit_distance=True, with_aln_seqs=True):
        self.fastx = fastx
        self.fastx_file = fastx_file
        self.query_seq = query_seq
        if query_seq:
            self.query_seq_len = float(len(query_seq))
        self.multiple_query_sequences_file = multiple_query_sequences
        if multiple_query_sequences:
            self.__get_multiple_query_sequences_info__()
            sys.stderr.write(str(self.multiple_query_sequences)+"\n")
            ##DELETE print
        self.with_read_names = with_read_names
        self.with_edit_distance = with_edit_distance
        self.with_aln_seqs = with_aln_seqs
    def __get_multiple_query_sequences_info__(self):
        self.multiple_query_sequences = {}
        for line in open(self.multiple_query_sequences_file, 'r'):
            seq = line.strip()
            self.multiple_query_sequences[seq] = float(len(seq))
    def fitting_aln(self, read):
        return generalizedFittingAln(read, self.query_seq)[:4]
    def print_fitting_aln_stats_over_all_reads(self):
        editdist = ''
        for fx in SeqIO.parse(self.fastx_file, self.fastx):
            readname = (fx.name+"\t")*self.with_read_names
            score, start, alnseq1, alnseq2 = self.fitting_aln(str(fx.seq))
            if self.with_edit_distance:
                editdist = "\t" + str(editDistance(alnseq1, alnseq2)) 
            print readname+("\t").join([str(start), str(score), str(score/self.query_seq_len)]) + editdist + ("\t"+alnseq1+"\t"+alnseq2)*self.with_aln_seqs
    def best_fitting_aln(self, read):
        maxscore = float("-inf")
        num_with_max_score = 1
        for seq in self.multiple_query_sequences:
            score, start, alnseq1, alnseq2 = generalizedFittingAln(read, seq)[:4]
            if score > maxscore:
                num_with_max_score = 1
                maxscore = score
                fit_aln_info = (seq, score, start, alnseq1, alnseq2)
            elif score == maxscore:
                num_with_max_score += 1
                if random.randint(0,1):
                    fit_aln_info = (seq, score, start, alnseq1, alnseq2)
        return (fit_aln_info[0], fit_aln_info[1], fit_aln_info[2], fit_aln_info[3],fit_aln_info[4], num_with_max_score)
    def print_best_fitting_aln_stats_over_all_reads(self):
        editdist = ''
        n = 0
        for fx in SeqIO.parse(self.fastx_file, self.fastx):
            n += 1
            sys.stderr.write(str(n)+"\n")
            readname = (fx.name+"\t")*self.with_read_names
            winner, score, start, alnseq1, alnseq2, num_with_max_score = self.best_fitting_aln(str(fx.seq))
            if self.with_edit_distance:
                editdist = "\t" + str(editDistance(alnseq1, alnseq2))
            seq_len = self.multiple_query_sequences[winner]
            print readname+("\t").join([str(num_with_max_score), winner, str(start), str(score), str(score/seq_len)]) + editdist + ("\t"+alnseq1+"\t"+alnseq2)*self.with_aln_seqs


def run(parser, args):
    filename, fastx = filetype(args)
    if args.random_sequence:
        args.sequence = random_seq(args.random_sequence)
    ## seq trans?
    if args.reverse_complement:
        args.sequence = reverseComplement(args.sequence)
    elif args.complement:
        args.sequence = complement(args.sequence)
    elif args.reverse_sequence:
        args.sequence = reverse_seq(args.sequence)
    if not args.multiple_sequences:
        fit = FitAln(fastx, filename, args.sequence,None, args.with_read_names, args.with_edit_distances, args.with_aln_seqs)
        fit.print_fitting_aln_stats_over_all_reads()
    else:
        fit = FitAln(fastx, filename, None, args.multiple_sequences, args.with_read_names, args.with_edit_distances, args.with_aln_seqs)
        fit.print_best_fitting_aln_stats_over_all_reads()
    
