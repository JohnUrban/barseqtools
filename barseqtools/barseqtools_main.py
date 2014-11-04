#!/usr/bin/env python

import os.path
import sys
import argparse

#logger
import logging
logger = logging.getLogger('barseqtools')

# barseqtools imports
import barseqtools.version

def run_subtool(parser, args):
    if args.command == 'count_samples':
        import count_sample_indexes as submodule
    if args.command == 'count_barcodes':
        import count_yeast_barcodes as submodule
    if args.command == 'count_subseq':
        import count_universal_primer_bits as submodule
    if args.command == 'pilot':
        import pilot as submodule
    if args.command == 'fitting_aln':
        import fitting_aln as submodule
    elif args.command == 'kmer':
        import kmer as submodule
    elif args.command == 'kmerplot':
        import kmerplot as submodule
    elif args.command == 'kmerdiff':
        import kmerdiff as submodule

    # run the chosen submodule.
    submodule.run(parser, args)

class ArgumentParserWithDefaults(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)
	self.add_argument("-q", "--quiet", help="Do not output warnings to stderr",
                        action="store_true",
                        dest="quiet")

def main():
    logging.basicConfig()

    #########################################
    # create the top-level parser
    #########################################
    parser = argparse.ArgumentParser(prog='barseqtools', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v", "--version", help="Installed barseqtools version",
                        action="version",
                        version="%(prog)s " + str(barseqtools.version.__version__))
    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command', parser_class=ArgumentParserWithDefaults)

    #########################################
    # create the individual tool parsers
    #########################################



    ##########
    # sample index counting
    ##########
    parser_count_samples = subparsers.add_parser('count_samples',
                                        help='''Count unique sample indexes in first 5 bp of reads. Appropriate for early pilot analyses on a subsample of all reads.
                                                Right now, exact matches are used to count samples.
                                                In future, it will also analyze non-exact matches for lowest hamming distance index < some cutoff.
                                                Will give exact match information on the following:
                                                number unique 5mers found, total number 5mers counted (should be equal to the
                                                number of reads analyzed), number of indexes anticipated (should be equal to the number of indexes
                                                in index file provided), number/percent of anticipated indexes found,
                                                number/percent of 5mers encountered that make up anticipated indexes, counts seen per sample.''')
    parser_count_samples_filetype = parser_count_samples.add_mutually_exclusive_group()
    parser_count_samples_filetype.add_argument('-fa', '--fasta',
                              dest='fasta',
                              default=None,
                              type=str,
                              help=('''Specify "--fasta file.fa" for analyzing a fasta file.'''))
    parser_count_samples_filetype.add_argument('-fq', '--fastq',
                              dest='fastq',
                              default=None,
                              type=str,
                              help=('''Specify "--fasta file.fq" for analyzing a fastq file.'''))
    parser_count_samples.add_argument('--indexes',
                             dest='indexes',
                             type=str,
                             help='''Path to tab-delim file with: sample_name, index, ...any....''',
                             default=None, required=True)
    parser_count_samples.set_defaults(func=run_subtool)


    ##########
    # yeast barcode counting
    ##########
    parser_count_barcodes = subparsers.add_parser('count_barcodes',
                                        help='NEW FEATURE -- NOT YET STABLE/FINISHED. Count unique 20mers in position 20-40 of reads.')
    parser_count_barcodes_filetype = parser_count_barcodes.add_mutually_exclusive_group()
    parser_count_barcodes_filetype.add_argument('-fa', '--fasta',
                              dest='fasta',
                              default=None,
                              type=str,
                              help=('''Specify "--fasta file.fa" for analyzing a fasta file.'''))
    parser_count_barcodes_filetype.add_argument('-fq', '--fastq',
                              dest='fastq',
                              default=None,
                              type=str,
                              help=('''Specify "--fasta file.fq" for analyzing a fastq file.'''))
    parser_count_barcodes.add_argument('--barcodes',
                             dest='barcodes',
                             type=str,
                             help='''Path to tab-delim file with: ORF_name, 20mer-barcode, ...any....''',
                             default=None, required=True)
    parser_count_barcodes.add_argument('-hd', '--hamming',
                             dest='hamming_analysis',
                             action="store_true",
                             help='''Instead of looking at exact match stats (which is hamming distance of 0), perform
                                    a more thorough hamming distance analysis looking at unique assignments if up to 1 or 2 mismatches
                                    is allowed (and unique barcode that minimizes hamDist is selected if one exists).''',
                             default=False)
    parser_count_barcodes.set_defaults(func=run_subtool)



    ##########
    # subsequence_counting ---> originally used for universal primer counting (positions 5-19) shared by all reads
    ##########
    parser_count_subseq = subparsers.add_parser('count_subseq',
                                        help='NEW FEATURE -- NOT YET STABLE/FINISHED. Default: Count shared 15mers in positions 5-20 of reads (in pyhton-ese coordinates). Alt: specify any sequence and starting position such that the sequence is not longer than the read when starting at that position.')
    parser_count_subseq_filetype = parser_count_subseq.add_mutually_exclusive_group(required=True)
    parser_count_subseq_filetype.add_argument('-fa', '--fasta',
                              dest='fasta',
                              default=None,
                              type=str,
                              help=('''Specify "--fasta file.fa" for analyzing a fasta file.'''))
    parser_count_subseq_filetype.add_argument('-fq', '--fastq',
                              dest='fastq',
                              default=None,
                              type=str,
                              help=('''Specify "--fasta file.fq" for analyzing a fastq file.'''))
    parser_count_subseq.add_argument('--sequence',
                             dest='sequence',
                             type=str,
                             help='''Sequence to match specified portion of reads against. Default: GTCCACGAGGTCTCT.
                                    The default is the universal 15 base sequence after the first 5 index bases in the bar-seq primers
                                    and should be shared by all reads. The start position in pythonese is 5, which is the default for --start''',
                             default="GTCCACGAGGTCTCT")
    parser_count_subseq.add_argument('--start',
                             dest='start',
                             type=int,
                             help='''Start position (python coordinates are 0-based) in read - for a 50 bp read, this can be 0-49.
                                    Default = 5 for the default analysis of matching the universal 15 bp in bar-seq primers shared by all reads.''',
                             default=5)
    parser_count_subseq.add_argument('--min-read-length',
                             dest='min_read_length',
                             type=int,
                             help='''Min read length expected. Takes an integer. Default 50 bp. ''',
                             default=50)
    parser_count_subseq.set_defaults(func=run_subtool)




    ##########
    # pilot analysis - to be run on the 20k reads I sampled
    ##########
    parser_pilot = subparsers.add_parser('pilot',
                                        help='NEW FEATURE -- NOT YET STABLE/FINISHED. Count unique sample indexes in first 5 bp of reads.')
    parser_pilot_filetype = parser_pilot.add_mutually_exclusive_group()
    parser_pilot_filetype.add_argument('-fa', '--fasta',
                              dest='fasta',
                              default=None,
                              type=str,
                              help=('''Specify "--fasta file.fa" for analyzing a fasta file.'''))
    parser_pilot_filetype.add_argument('-fq', '--fastq',
                              dest='fastq',
                              default=None,
                              type=str,
                              help=('''Specify "--fasta file.fq" for analyzing a fastq file.'''))
    parser_pilot.add_argument('--outprefix',
                              type=str,
                              required=True,
                              help="Provide an output prefix that will be appended to the beginning of output files.")
    parser_pilot.add_argument('--indexes',
                             dest='indexes',
                             type=str,
                             help='''Path to tab-delim file with: sample_name, index, ...any....''',
                             default=None, required=True)
    parser_pilot.add_argument('--barcodes',
                             dest='barcodes',
                             type=str,
                             help='''Path to tab-delim file with: ORF_name, 20mer-barcode, ...any....''',
                             default=None, required=True)
## TO ADD IN -- get info on proportion that had perfect match primer bits
##    parser_pilot.add_argument('--sequence',
##                             dest='sequence',
##                             type=str,
##                             help='''Path to tab-delim file with: ORF_name, 20mer-barcode, ...any....''',
##                             default="GTCCACGAGGTCTCT")
##    parser_pilot.add_argument('--start',
##                             dest='start',
##                             type=int,
##                             help='''Start position (python coordinates are 0-based) in read - for a 50 bp read, this can be 0-49.''',
##                             default=5)
##    parser_pilot.add_argument('--min-read-length',
##                             dest='min_read_length',
##                             type=int,
##                             help='''Min read length expected. Takes an integer. Default 50 bp. ''',
##                             default=50)
    parser_pilot.set_defaults(func=run_subtool)





    ##########
    # fitting alignment
    ##########
    parser_fitting_aln = subparsers.add_parser('fitting_aln',
                                        help='NEW FEATURE -- NOT YET STABLE/FINISHED. Performs fitting alignment of some string over all reads -- returns stats')
    parser_fitting_aln_filetype = parser_fitting_aln.add_mutually_exclusive_group(required=True)
    parser_fitting_aln_filetype.add_argument('-fa', '--fasta',
                              dest='fasta',
                              default=None,
                              type=str,
                              help=('''Specify "--fasta file.fa" for analyzing a fasta file.'''))
    parser_fitting_aln_filetype.add_argument('-fq', '--fastq',
                              dest='fastq',
                              default=None,
                              type=str,
                              help=('''Specify "--fasta file.fq" for analyzing a fastq file.'''))
    parser_fitting_aln.add_argument('--sequence',
                             dest='sequence',
                             type=str,
                             help='''Provide a sequence < min read length. Default is the KanRR fragment: CGTACGCTGCAGGTCG''',
                             default='CGTACGCTGCAGGTCG') ##shortest version of KanR fragment that maximized scores on some pilot reads
    parser_fitting_aln.add_argument('-m', '--multiple-sequences',
                             dest='multiple_sequences',
                             type=str, default=None,
                             help='''Provide a path to a file with 1 sequence per line.
                                    For each read in the fastx file, it will report the fitting alignment for the
                                    sequence in this file with the best fitting aln.
                                    Each time it encounters a score that ties the current max score, it exchanges the older fiting aln
                                    info for the new fitting aln info with a 50%% probability.
                                    This way there is a random assignment of the best barcode.
                                    Use --all-scores instead to get an output with all max scores and barcodes returned.''')
    parser_fitting_aln.add_argument('-w', '--with-read-names',
                             dest='with_read_names', action="store_true",
                             default=False,
                             help='''If set, will print "readname, startPosInRead, fitAlnScore, fitAlnScore/queryLen";
                                else just "startPosInRead,fitAlnScore, fitAlnScore/queryLen".
                                Start position is in pythonese (0-based).''')
    parser_fitting_aln.add_argument('-e', '--with-edit-distances',
                             dest='with_edit_distances', action="store_true",
                             default=False,
                             help='''If set, edit dist will be incl in output''')
    parser_fitting_aln.add_argument('-a', '--with-aln-seqs',
                             dest='with_aln_seqs', action="store_true",
                             default=False,
                             help='''If set, the aligned versions of sequences 1 (read) and 2 (provided) will be printed.''')
    parser_fitting_aln.add_argument('-r', '--random-sequence',
                             dest='random_sequence', type=int,
                             default=False,
                             help='''Provide integer for random sequence length. This option overrides --sequence.''')
    parser_fitting_aln_seqtransform = parser_fitting_aln.add_mutually_exclusive_group()
    parser_fitting_aln_seqtransform.add_argument("-c", "--complement", action="store_true", default=False,
                                                 help=''' Use complement of provided sequence -- right now only works on single seq.
                                                            e.g. AACC -> TTGG''')
    parser_fitting_aln_seqtransform.add_argument("-rc", "--reverse_complement", action="store_true", default=False,
                                                 help=''' Use reverse complement of provided sequence -- right now only works on single seq.
                                                            e.g. AACC -> GGTT''')
    parser_fitting_aln_seqtransform.add_argument("-rs", "--reverse_sequence", action="store_true", default=False,
                                                 help=''' Use reverse sequence of provided sequence -- right now only works on single seq.
                                                            e.g. AACC -> CCAA''')

    parser_fitting_aln.set_defaults(func=run_subtool)
## I ran some tests and doing the fitting aln is more sensitive for barcodes that match 100% but start at 18 or 19 instead of 20
    ## saw one perfect aln at pos 55...






    ##########
    # kmerCounting
    ##########
    parser_kmer = subparsers.add_parser('kmer',
                                        help='NEW FEATURE -- NOT YET STABLE/FINISHED. Count kmers in reads or reference.')
    parser_kmer.add_argument('files', metavar='FILES', nargs='+',
                             help='The input FAST5 files.')
    parser_kmer.add_argument('-k', '--kmersize',
                              dest='k',
                              default=5,
                              type=int,
                              help=('Kmer size. Default = 5. Sizes 1-7 work well with kmerplot on regular Mac OS. Up to 10 is possible. After that it might require too much memory for kmerplot on regular Mac OS.'))
    parser_kmer.add_argument('--fasta',
                              dest='fasta',
                              default=None,
                              type=str,
                              help=('''Specify "--fasta file.fa" for analyzing a fasta file instead of fast5dir/.
                                    While min and max length arguments remain meaningful for fasta files, the following arguments do not: start time, end time, high quality, type, single read per molecule.'''))
    parser_kmer.add_argument('--fastq',
                              dest='fastq',
                              default=None,
                              type=str,
                              help=('''Specify "--fasta file.fq" for analyzing a fastq file instead of fast5dir/.
                                    While min and max length arguments remain meaningful for fastq files, the following arguments do not: start time, end time, high quality, type, single read per molecule.'''))
    parser_kmer.add_argument('--rev-comp',
                              dest='rev_comp',
                              default=False,
                              action="store_true",
                              help='''Created to be used with --fasta and --fastq options.
                                    When creating kmer counts, it counts both the fwd and reverse complement kmers.
                                    For now, it does nothing when used with fast5 dirs (minION data files).''')

##    parser_kmer_output = parser_kmer.add_mutually_exclusive_group()
##    parser_kmer_output.add_argument('-t', '--table',
##                              dest='table',
##                              default=True,
##                              action='store_true',
##                              help=('''Output option: report tab-delimited table of kmer, count, and proportion of all kmers seen.
##                                    Default = True (to stdout). Use --saveas to specify file to save to.'''))
##    parser_kmer_output.add_argument('-p', '--plot',
##                              dest='plot',
##                              default=False,
##                              action='store_true',
##                              help=('''Output option: show or write out plot.
##                                    Default = False (to stdout). Use --saveas to specify file to save to.'''))
    parser_kmer.add_argument('--min-length',
                              dest='min_length',
                              default=0,
                              type=int,
                              help=('Minimum read length to be included in analysis.'))
    parser_kmer.add_argument('--max-length',
                              dest='max_length',
                              default=1000000000,
                              type=int,
                              help=('Maximum read length to be included in analysis.'))
    parser_kmer.add_argument('--start',
                              dest='start_time',
                              default=None,
                              type=int,
                              help='Only analyze reads from after start timestamp')
    parser_kmer.add_argument('--end',
                              dest='end_time',
                              default=None,
                              type=int,
                              help='Only analyze reads from before end timestamp')
    parser_kmer.add_argument('--high-quality',
                              dest='high_quality',
                              default=False,
                              action='store_true',
                              help='Only analyze reads with more complement events than template.')
    parser_kmer.add_argument('--saveas',
                             dest='saveas',
                             metavar='STRING',
                             help='''Save tab-delimited kmer + counts to file.''',
                             default=None)
    parser_kmer_readfilter = parser_kmer.add_mutually_exclusive_group()
    parser_kmer_readfilter.add_argument('--type',
                              dest='type',
                              metavar='STRING',
                              choices=['all', 'fwd', 'rev', '2D', 'fwd,rev'],
                              default='all',
                              help='Which type of reads should be analyzed? Def.=all, choices=[all, fwd, rev, 2D, fwd,rev]. Is mutually exclusive with --one-read-per-molecule.')
    parser_kmer_readfilter.add_argument('-1', '--one-read-per-molecule',
                              dest='single_read',
                              default=False,
                              action='store_true',
                              help='''Only analyze one read per molecule in priority order: 2D -> template -> complement.
                                            That is, if there is a 2D read use that.If not, then try to use template. etc.
                                            Is mutually exclusive with --type.''')
    parser_kmer.set_defaults(func=run_subtool)

    
    ##########
    # kmerplotting
    ##########
    parser_kmerplot = subparsers.add_parser('kmerplot',
                                        help='NEW FEATURE -- NOT YET STABLE/FINISHED. plot kmer counts in reads or reference.')
##    parser_kmerplot.add_argument('files', metavar='FILES', nargs='+',
##                             help='The input FAST5 files.')

    parser_kmerplot.add_argument('-t1', '--kmer-count-in-reads',
                             dest='table1',
                             type=str,
                             help='''Provide path to file with kmer count table from reads (or any kmer count table).
                                    This argument is required and when used alone, just generates a bar plot of kmer counts.''',
                             default=None)
    
    parser_kmerplot.add_argument('-t2', '--kmer-count-in-reference',
                             dest='table2',
                             type=str,
                             help='''Provide path to file with kmer count table from reference sequence (or any second kmer count table).
                                    This argument is not required and if used, results in a scatterplot of the 2 kmer count tables.''',
                             default=None)
    parser_kmerplot.add_argument('--matplotlib',
                             dest='mpl',
                             action='store_true',
                             help='''Temp option: plot in matplotlib''',
                             default=False)
    parser_kmerplot.add_argument('--ggplot2',
                             dest='gg',
                             action='store_true',
                             help='''Temp option: plot in ggplot2''',
                             default=False)
    parser_kmerplot.add_argument('--saveas',
                             dest='saveas',
                             metavar='STRING',
                             help='''Save to file. e.g. --saveas "filename.extension" where extension can be only pdf and jpg for now.''',
                             default=None)
    parser_kmerplot.set_defaults(func=run_subtool)
    

    ##########
    # kmer diff abundance
    ##########
    parser_kmerdiff = subparsers.add_parser('kmerdiff',
                                        help='NEW FEATURE -- NOT YET STABLE/FINISHED. plot kmer counts in reads or reference.')

    parser_kmerdiff.add_argument('-t1', '--kmer-count-in-reads',
                             dest='table1',
                             type=str,
                             help='''Provide path to file with kmer count table from reads (or any kmer count table).
                                    This argument is required and when used alone, just generates a bar plot of kmer counts.''',
                             default=None)
    
    parser_kmerdiff.add_argument('-t2', '--kmer-count-in-reference',
                             dest='table2',
                             type=str,
                             help='''Provide path to file with kmer count table from reference sequence (or any second kmer count table).
                                    This argument is not required and if used, results in a scatterplot of the 2 kmer count tables.''',
                             default=None)
    parser_kmerdiff.add_argument('--saveas',
                             dest='saveas',
                             metavar='STRING',
                             help='''Save to file. e.g. --saveas "filename.extension" where extension can be only pdf and jpg for now.''',
                             default=None)
    parser_kmerdiff.add_argument('-bcv', '--square-root-dispersion',
                             dest='bcv',
                             type=float,
                             help='''When there are no replicates in edgeR, dispersion must be determined by the user.
                                    The default is 0.2. Other values to try could be 0.01-0.4 (or any).
                                    p-values will be sensitive to choice of bcv. Fold change will not.''',
                             default=0.2)

    parser_kmerdiff.add_argument('--volcano',
                             dest='volcano',
                             type=str,
                             help='''If you want the analysis to generate a volcano plot,
                                    (log(fold change) vs. -log10(pvalue)), then use this flag
                                    and provide the name and extension of volcano plot file (e.g. volcano.jpg).''',
                             default=None)
    parser_kmerdiff.add_argument('--smear',
                             dest='smear',
                             type=str,
                             help='''If you want the analysis to generate a smear plot,
                                    (log(fold change) vs. log(CPM)), then use this flag
                                    and provide the name and extension of smear plot file (e.g. smear.jpg).''',
                             default=None)
    parser_kmerdiff.add_argument('--nt-content',
                             dest='nt_content',
                             type=str,
                             help='''If you want the analysis to generate a table analyzing the nucleotide content
                                    of kmers >= abs(fold change) and pval <= p, then use this flag with those values as in these
                                    examples: (a) --nt-content fc:2.5,p:0.001 (b) --nt-content fc:2,fdr:0.1)''',
                             default=None)
    
    parser_kmerdiff.set_defaults(func=run_subtool)


    

    #######################################################
    # parse the args and call the selected function
    #######################################################
    args = parser.parse_args()

    if args.quiet:
        logger.setLevel(logging.ERROR)

    try:
      args.func(parser, args)
    except IOError, e:
         if e.errno != 32:  # ignore SIGPIPE
             raise

if __name__ == "__main__":
    main()
