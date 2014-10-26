#Table entries:
# Protein condition gene rep_num_1_count rep_num_2_count rep_num_3_count 

## OR
# Each row is a specific gene, with counts for all samples

# algorithm
# 1. read in sample indexes to make sample:index dcts...
# 2. do that for barcodes:genes
# 3. Go through reads and
#       a. get index
#       b. ensure universal primer bit is accurate to move on
#       c. get barcode
#       d. add barcode count to that index 
# 4. translate index+barcode counts to pro:cond:rep_num counts
# 5. write out table that can be read in later to bypass steps 1-4.
# 6. For each protein:
#           (a) Do diff expr with the cond1_reps vs. cond2_reps
#               ->  since when given A and B the FE=B/A, A=cond_2, B=cond_1
#               -> this way we see after/before -- enriched after and depelted after

## eventually:
# 1. Try median normalization
#   Followed by Z-score? Would Z-score be any diff before/aft medianNorm, I dont think so
#   
# 2. Try quantile norm
# 3. Add in alignments and edit distance allow edit dist of 2 in barcode, 1 in sample index.

## check the distance between all indexes and all barcodes....
## ... e.g. hamming/edit
from sequences import *
from collections import defaultdict

class BarSeq(object):
    def __init__(self, fastx, fastx_fh, indexes_fh, barcodes_fh, output_prefix):
        self.fastx = fastx # 'fasta' or 'fastq'
        self.fastx_fh = fastx_fh #path to file
        self.index_fh = indexes_fh #path to file
        self.barcodes_fh = barcodes_fh #path to file
        self.output_prefix = output_prefix
        # lazy start
        self.index2sample = None
        self.sample2index = None
        self.barcode2gene = None
        self.gene2barcode = None
        self.index_barcode_counts = None
        self.sample_gene_counts = None
        
        
    def make_index_to_sample_dict(self):
        ''' Assumes index file has following 4-5 columns:
                protein_name, condition, rep_num, index, index_inside_primer_used(opt)'''
        self.index2sample, self.sample2index = make_index_to_sample_dict(self.index_fh)
                
    def make_barcode_to_gene_dict(self):
        ''' Assumes index file has following 2 columns:
                gene_name barcode'''
        self.barcode2gene, self.gene2barcode = make_barcode_to_gene_dict(self.barcodes_fh)

    def make_index_barcode_counts_dict(self):
        ''' makes a dict d s.t. s[sample][barcode] gives count'''
        self.index_barcode_counts = index_barcode_counter(self.fastx_fh, self.fastx)

    def write_index_barcode_counts_table(self, out_fh=None):
        if not self.index_barcode_counts:
            self.make_index_barcode_counts_dict()
        if not out_fh:
            out_fh = self.output_prefix + "-index_barcode_counts.tsv"
        write_index_barcode_counts_table(self.index_barcode_counts, out_fh)

    def equalize_index_barcode_dicts(self):
        if not self.index_barcode_counts:
            self.make_index_barcode_counts_dict()
        self.index_barcode_counts, self.barcode2gene, self.gene2barcode, self.index2sample, self.sample2index = equalize_index_barcode_dicts(self.index_barcode_counts, self.barcode2gene, self.gene2barcode, self.index2sample, self.sample2index)

    def indexbarcode_to_samplegene_dict(self):
        self.sample_gene_counts = indexbarcode_to_samplegene_dict(self.index_barcode_counts, self.index2sample, self.barcode2gene)

    def indexbarcode_to_samplebarcode_dict(self):
        self.sample_barcode_counts = indexbarcode_to_samplebarcode_dict(self.index_barcode_counts, self.index2sample, self.barcode2gene)


    def write_sample_gene_counts_table_1(self, out_fh=None):
        if not self.index_barcode_counts:
            self.make_index_barcode_counts_dict()
            self.index_barcode_counts = equalize_index_barcode_dict(self.index_barcode_counts)
        if not out_fh:
            out_fh = self.output_prefix + "-sample_gene_counts_1.tsv"
        write_sample_gene_counts_table_1(self.index_barcode_counts, self.index2sample, self.barcode2gene, out_fh)

    def write_sample_gene_counts_table_2(self, out_fh=None):
        if not self.index_barcode_counts:
            self.make_index_barcode_counts_dict()
            self.index_barcode_counts = equalize_index_barcode_dict(self.index_barcode_counts)
        if not out_fh:
            out_fh = self.output_prefix + "-sample_gene_counts_2.tsv"
        write_sample_gene_counts_table_2(self.sample_barcode_counts, self.sample2index, self.barcode2gene, out_fh)

    def analyze(self):
        ''' Assumes index file has following 4-5 columns:
                protein_name, condition, rep_num, index, index_inside_primer_used(opt)'''
        self.make_index_to_sample_dict()
        self.make_barcode_to_gene_dict()
        self.make_index_barcode_counts_dict()
        self.write_index_barcode_counts_table()
        self.equalize_index_barcode_dicts()
        self.indexbarcode_to_samplegene_dict() ## will prob trim this one
        self.indexbarcode_to_samplebarcode_dict()
        self.write_sample_gene_counts_table_1() ## will prob rm this one -- ensure equivalent to 2
        self.write_sample_gene_counts_table_2()

    def add_protein(protein):
        pass
    def add_condition(protein, condition):
        ## might be best calling them condition 0 and 1.
        pass
    def add_replicate(protein, condition, rep_num):
        pass
    def add_barcode_count(protein, condition, rep_num):
        pass
    def get_barcode_count_for_sample(index, barcode):
        pass
    def get_barcode_count_for_protein_condition(protein, condition, barcode):
        pass ## all reps
    def get_barcode_count_for_protein(protein, barcode):
        pass ## all reps from both conditions

    





def make_index_to_sample_dict(index_fh):
        ''' Assumes index file has following 4-5 columns:
                protein_name, condition, rep_num, index, index_inside_primer_used(opt)'''
        index2sample = dict()
        sample2index = dict()
        for line in open(index_fh):
            if line[0] != "#":
                line = line.strip().split()
                sample = (".").join([line[0],line[1],line[2]])
                index = line[3]
                index2sample[index] = sample
                sample2index[sample] = index
        return index2sample, sample2index

def make_barcode_to_gene_dict(barcodes_fh):
        ''' Assumes index file has following 2 columns:
                gene_name barcode'''
        barcode2gene = dict()
        gene2barcode = dict()
        for line in open(barcodes_fh):
            if line[0] != "#":
                line = line.strip().split()
                gene = line[0]
                barcode = line[1]
                barcode2gene[barcode] = gene
                gene2barcode[gene] = barcode
        return barcode2gene, gene2barcode


def write_index_barcode_counts_table(index_barcode_counts, out_fh):
    ## writes table:
    ## index, barcode, count, proportion of index
    out = open(out_fh, 'w')
    for index in sorted(index_barcode_counts.keys()): #sorted indexes
        index_total_counts = float(sum([index_barcode_counts[index][e] for e in index_barcode_counts[index].keys()]))
        for barcode in sorted(index_barcode_counts[index].keys()):
            count = index_barcode_counts[index][barcode]
            line = ("\t").join([index, barcode, str(count), str(count/index_total_counts)])
            out.write(line+"\n")
    out.close()


def equalize_index_barcode_dicts(index_barcode_counts, barcode2gene, gene2barcode, index2sample, sample2index):
    ## possible also add un-anticipated barcodes to barcode2gene: "un-anticipated":"un-assigned"
    ## --> does not work for gene2barcode though: "un-assigned" --> many un-anticipated barcodes
    all_barcodes = set()
    all_indexes = set()
    
    ## add all antincipated barcodes, add all anticipated indexes
    for barcode in barcode2gene.keys():
        all_barcodes.add(barcode)
    for index in index2sample.keys():
        all_indexes.add(index)
        
    ## search for and add all un-anticipated indexes and all un-anticipated barcodes 
    for index in index_barcode_counts.keys():
        all_indexes.add(index)
        for barcode in index_barcode_counts[index].keys():
            all_barcodes.add(barcode)
            
    ## ensure all anticipated indexes are part of index_barcode_counts
    ##      and make sure all un-anticipated indexes become part of index2sample
    for index in all_indexes:
        if index not in index_barcode_counts:
            index_barcode_counts[index]
        if index not in index2sample:
            index2sample[index] = "unassigned_index_" + index + ".unknown_condition.unknown_rep_num"
            sample = index2sample[index]
            sample2index[sample] = index
            
    ## ensure all samples have all barcodes (if not, add barcode with 0 count)
    ##      and ensure all un-anticipated barcodes become part of barcode2gene
    for index in index_barcode_counts.keys():
        for barcode in all_barcodes:
            try:
                index_barcode_counts[index][barcode] += 0
            except KeyError:
                index_barcode_counts[index][barcode] = 0
            #update barcode2gene
            if barcode not in barcode2gene:
                barcode2gene[barcode] = "unassigned_barcode_" + barcode
                gene = barcode2gene[barcode]
                gene2barcode[gene] = barcode
                
##    ## update barcode2gene
##    for barcode in all_barcodes:
##        if barcode not in barcode2gene:
##            barcode2gene[barcode] = "unassigned_barcode_" + barcode
    return index_barcode_counts, barcode2gene, gene2barcode, index2sample, sample2index



def indexbarcode_to_samplegene_dict(index_barcode_counts, index2sample, barcode2gene):
    sample_gene_counts = defaultdict(dict)
    for index in index_barcode_counts.keys():
        for barcode in index_barcode_counts[index].keys():
            sample_gene_counts[index2sample[index]][barcode2gene[barcode]] = index_barcode_counts[index][barcode]
    return sample_gene_counts


## may need to use this translation because of the overload "unassigned" gene.
def indexbarcode_to_samplebarcode_dict(index_barcode_counts, index2sample, barcode2gene):
    sample_barcode_counts = defaultdict(dict)
##    for index in index_barcode_counts.keys():
##        sample_barcode_counts[index2sample[index]] = index2sample[index]
    for index in index_barcode_counts.keys():
        for barcode in index_barcode_counts[index].keys():
            sample_barcode_counts[index2sample[index]][barcode] = index_barcode_counts[index][barcode]
    return sample_barcode_counts       

def proportion_string(numerator, denominator):
    try:
        prop = numerator/denominator
    except ZeroDivisionError:
        if numerator == 0:
            prop = 0
        else:
            prop = 'undefined'
    return str(prop)

def write_sample_gene_counts_table_1(index_barcode_counts, index2sample, barcode2gene, out_fh):
    ''' index2sample is a dict that with index:sample, sample protein_name.condition.rep_num
    barcode to gene is a dict with barcode:gene_name
    Tip: equalize all samples to have equal sets of barcodes in dict first.'''
    ## writes table:
    ## protein_name condition rep index barcode gene count proportion_of_sample_counts
    ## This is essentially a more complete version of the "index_barcode_counts_table"
    ## returns a new dict...
    out = open(out_fh, 'w')
    for index in sorted(index_barcode_counts.keys()): #sorted indexes
        index_total_counts = float(sum([index_barcode_counts[index][e] for e in index_barcode_counts[index].keys()]))
        protein_name, condition, rep_num = index2sample[index].split(".")
        for barcode in sorted(index_barcode_counts[index].keys()):
            gene = barcode2gene[barcode]
            count = index_barcode_counts[index][barcode]
            proportion = proportion_string(count, index_total_counts)
            line = ("\t").join([protein_name, condition, rep_num, index, barcode, gene, str(count), str(proportion)])
            out.write(line+"\n")
    out.close()


def write_sample_gene_counts_table_2(sample_barcode_counts, sample2index, barcode2gene, out_fh):
    # this has benefit of writing to file in alphabetical order
    # however, the other output could be sorted by piping in unix --> | sort -k1,1 -k2,2 -k3,3n
    ''' index2sample is a dict that with index:sample, sample protein_name.condition.rep_num
    barcode to gene is a dict with barcode:gene_name
    Tip: equalize all samples to have equal sets of barcodes in dict first.'''
    ## writes table:
    ## protein_name condition rep index barcode gene count proportion_of_sample_counts
    ## This is essentially a more complete version of the "index_barcode_counts_table"
    ## returns a new dict...
    out = open(out_fh, 'w')
    for sample in sorted(sample_barcode_counts.keys()): #sorted protein.condition.rep_num
        sample_total_counts = float(sum([sample_barcode_counts[sample][e] for e in sample_barcode_counts[sample].keys()]))
        protein_name, condition, rep_num = sample.split(".")
        index = sample2index[sample]
        for barcode in sorted(sample_barcode_counts[sample].keys()):
            gene = barcode2gene[barcode]
            count = sample_barcode_counts[sample][barcode]
            proportion = proportion_string(count, sample_total_counts)
            line = ("\t").join([protein_name, condition, rep_num, index, barcode, gene, str(count), str(proportion)])
            out.write(line+"\n")
    out.close()


