## possibly useful functions
import copy
from collections import defaultdict
import numpy
import random

def random_seq(seqLen):
    '''Depends on random module'''
    bases = 'ACGT'
    seq = ''
    for j in range(seqLen):
        index = random.randint(0,3)
        seq += bases[index]
    return seq


def aligns(read, sequence, d):
    '''Returns true if read aligns anywhere on sequence with <= with hamming distance d mismatches
    d is max hamming distance'''
    for i in range(len(sequence)-len(read)+1):
        mm = 0
        for j in range(len(read)):
            mm += read[j] != sequence[i+j]
            if mm > d:
                break
        if mm <= d:
            return True
    return False    



def hammingDistSameLen(kmer1, kmer2):
    """assumes len(string1) == len(string2)
        assumes strings are 'aligned'
        automatically makes all letters uppercase for comparison"""
    if len(kmer1) != len(kmer2):
        return "Error: len(string1) != len(string2)"
    kmer1 = kmer1.upper()
    kmer2 = kmer2.upper()
    hamm=0
    for i in range(len(kmer1)):
        hamm = hamm + (kmer1[i] != kmer2[i])
    return int(hamm)



def hammingDistDiffLen(kmer1, kmer2):
    ''' Takes in 2 DNA strings of different lengths
    Returns the minimum Hamming distance
    ##NOTE: the minimum Hamming distance is considered the distance between the k-mer/read and the sequence/longer string'''
    #initialize
    if len(kmer1) > len(kmer2):
        sequence = kmer1
        read = kmer2
    else:
        sequence = kmer2
        read = kmer1
    minHam = hammingDistSameLen(read, sequence[0:len(read)])
    minIndex = 0
    startPositions = [0]
    for i in range(1, len(sequence)-len(read)+1):
        hamdist = hammingDistSameLen(read, sequence[i:i+len(read)])
        if hamdist < minHam:
            minHam = hamdist
    return minHam

def hammingDist(kmer1, kmer2):
    '''Takes in any 2 kmers where k1 not nec. equal to k2
    returns the hamming distance accordingly'''
    if len(kmer1) == len(kmer2):
        hamdist = hammingDistSameLen(kmer1, kmer2)
    else:
        hamdist = hammingDistDiffLen(kmer1, kmer2)
    return hamdist
        

 
def hammingDistMult(kmer1, kmerList):
    '''Takes in kmer and a list of kmers where lengths can vary
    Returns sum of Hamming distances between kmer1 and all kmers in list
    It scores distances as follows:
        if len(kmer) = len(kmer1) it just does hamming dist
        if len(ker) > len(kmer1) or len(kmer) < len(kmer1) it returns minHammDist using the shorter to scan the longer'''
    sumHamm = 0
    kmer1Len = len(kmer1)
    for kmer in kmerList:
        sumHamm += hammingDist(kmer1, kmer)
    return sumHamm


def alignToMinHammDist(read, sequence):
    ''' Takes in 2 DNA strings: a read and a sequence where |read| <= |sequence|
    Returns the start position(s) that minimize(s) the hamming distance (i.e. minimizes mismatches)
    Returns the minimum Hamming distance
    ##NOTE: the minimum Hamming distance is considered the distance between the k-mer/read and the sequence/longer string'''
    #initialize
    minHam = hammingDistSameLen(read, sequence[0:len(read)])
    minIndex = 0
    startPositions = [0]
    for i in range(1, len(sequence)-len(read)+1):
        hamdist = hammingDistSameLen(read, sequence[i:i+len(read)])
        if hamdist < minHam:
            startPositions = [i]
            minHam = hamdist
        elif hamdist == minHam:
            startPositions += [i]
    return startPositions, minHam





def LCSinfo(seq1, seq2):
    '''Longest common subsequence between seq1 and seq2
    letters in subsequence not nec. consecutive,but are in order'''
    n = len(seq1)
    m = len(seq2)
    longestPathTo = numpy.zeros([n+1, m+1])
    down = 0 
    right = 1
    diagonal = 2
    directions = [down, right, diagonal]
    backtrack = numpy.zeros([n+1, m+1])
    # for i <- 0 to |seq1|
        # s_i,0 <- 0
    # for j <- 0 to |seq2|
        # s_0,j <- 0
    ## dont need to implement above b/c start out with 0s
    # for i <- 1 to |seq1|
    for i in range(1,n+1):
        # for j <-10 to |seq2|
        for j in range(1,m+1):
            # s_i,j <- max{s_i-1,j ;  s_i,j-1 ;  s_i-1,j-1 (+1 if seq1_i == seq2_v)
            if seq1[i-1] == seq2[j-1]: match = 1
            ### Need to do i-1,j-1 here b/c the sequences start:end 1:len+1 -- the 0 row and 0 col are coming into sequences
            else: match = 0
            possible = [longestPathTo[i-1,j], longestPathTo[i,j-1], longestPathTo[i-1,j-1] + match]
            #print possible
            indexOfMax = possible.index(max(possible))
            #print indexOfMax
            longestPathTo[i,j] = possible[indexOfMax]
            # backtrack_i,j = 'down' if s_i,j = s_i-1,j ;  right if s_i,j = s_i,j-1 ; diagonal if s_i,j = si-1,j-1 + 1 (if match)
            backtrack[i,j] = directions[indexOfMax] ## can actually just use index of max b/c indexOfMax IS the direction as defined above
    return n, m, longestPathTo, backtrack

def outputLCS(backtrack, seq1, i, j, LCS=""):
    down = 0
    right = 1
    backDirection = backtrack[i,j]
    while True:
        if i+1 == 0 or j+1 == 0: ##Then it will be all dashes one way or another -- no more in common
            return LCS
        elif backtrack[i,j] == down:
            backDirection = backtrack[i-1,j]
            i = i-1
        elif backDirection == right:
            backDirection = backtrack[i,j-1]
            j = j-1
        else: #it is diagonal
            LCS = seq1[i-1] + LCS
            ## need i-1 b/c seqs in mtrices start at indices 1,1 -- whereas index of seq starts at 0
            #print LCS
            backDirection = backtrack[i-1,j-1]
            i,j = i-1, j-1

def findLCS(seq1, seq2):
    n, m, longestPathTo, backtrack = LCSinfo(seq1, seq2)
    LCS = outputLCS(backtrack, seq1, n, m)
    return longestPathTo[n,m], LCS



def makeMatchMismatchDict(m=1, mm=-1):
    matchMismatchDict = defaultdict(dict)
    for b1 in "ACGT":
        for b2 in "ACGT":
            if b1 != b2:
                matchMismatchDict[b1][b2] = mm
            else:
                matchMismatchDict[b1][b2] = m
    return matchMismatchDict
                

class ScoreMatrix(object):
    '''Depends on "import copy"'''
    def __init__(self, matchMismatchDict, indelPenalty=None):
        self.scoreMatrix = copy.deepcopy(matchMismatchDict)
        if indelPenalty != None:
            self.indelPenalty = float(indelPenalty)
            self.scoreMatrix["-"] = {}
            for symbol in self.scoreMatrix.keys():
                self.scoreMatrix[symbol]["-"] = self.indelPenalty
                self.scoreMatrix["-"][symbol] = self.indelPenalty
            
    def getScore(self, rowname, colname):
        return self.scoreMatrix[rowname][colname]

    def __eq__(self, other):
        for row in other.scoreMatrix.keys():
            if row not in self.scoreMatrix.keys():
                return False
        for row in self.scoreMatrix.keys():
            if row not in other.scoreMatrix.keys():
                return False           
        for row in self.scoreMatrix.keys():
            for col in self.scoreMatrix.keys():
                try:
                    if self.scoreMatrix[row][col] != other.scoreMatrix[row][col]:
                        return False
                except KeyError:
                    return False
        return True
            
    def __str__(self):
        out = []
        cols = " "*5
        for symb in self.scoreMatrix.keys():
            cols += symb + " "*(5-len(symb))
        out += [cols]
        for row in self.scoreMatrix.keys():
            scores = row + ' '*(5-len(row))
            for col in self.scoreMatrix.keys():
                score = str(self.scoreMatrix[row][col])
                scores += score + ' '*(5-len(score))
            out += [scores]
        return '\n'.join(out)

## Example
Score = ScoreMatrix(makeMatchMismatchDict(m=0,mm=-1), indelPenalty=-1)

class longestPathObject(object):
    '''This makes a longestPathTo matrix more user friendly -- esp for printing'''
    def __init__(self, longestPathToMatrix, seq1, seq2):
        self.lpm = longestPathToMatrix
        self.seq1 = seq1
        self.seq2 = seq2
    def __str__(self):
        out = []
        cols = " "*10
        for symb in self.seq2:
            cols += symb + " "*(4)
        out += [cols]
        for i in range(len(self.seq1)+1):
            if i == 0:
                scores = ' '*(5)
            else:
                #print i
                scores = self.seq1[i-1] + ' '*4
            for score in self.lpm[i,:]:
                score = str(score)
                scores += score + ' '*(5-len(score))
            out += [scores]
        return '\n'.join(out)

class backtrackObject(longestPathObject):
    def __init__(self, backtrackMatrix, seq1, seq2, backtrack=True):
        longestPathObject.__init__(self, backtrackMatrix, seq1, seq2)
        if backtrack:
            self.directions = {0:"up", 1:"left", 2:"diag", 3:"free"}
        else: ## read forward directions
            self.directions = {0:"down", 1:"right", 2:"diag", 3:"free"}
    def __str__(self):
        out = []
        cols = " "*16
        for symb in self.seq2:
            cols += symb + " "*(7)
        out += [cols]
        for i in range(len(self.seq1)+1):
            if i == 0:
                scores = ' '*(8)
            else:
                scores = self.seq1[i-1] + ' '*7
            for score in self.lpm[i,:]:
                try:
                    score = self.directions[score]
                except KeyError:
                    score = "Src"
                scores += score + ' '*(8-len(score))
            out += [scores]
        return '\n\n'.join(out)





#####GLOBAL
def generalizedLCSinfo(seq1, seq2, Score):
    '''Longest common subsequence between seq1 and seq2
    letters in subsequence not nec. consecutive,but are in order
    "Score" is a ScoreMatrix Object'''
    n = len(seq1)
    m = len(seq2)
    longestPathTo = numpy.zeros([n+1, m+1])
    down = 0 
    right = 1
    diagonal = 2
    directions = [down, right, diagonal]
    backtrack = numpy.zeros([n+1, m+1])
    ## Need to set first row and first column to cumulative indels (and backtrack as all downs or all rights)
    for i in range(1, n+1):
        longestPathTo[i,0] = longestPathTo[i-1,0] + Score.getScore(seq1[i-1], "-")
    for j in range(1, m+1):
        longestPathTo[0,j] = longestPathTo[0,j-1] + Score.getScore("-",seq2[j-1])
    backtrack[0,1:] = 1
    backtrack[1:,0] = 0
    backtrack[0,0] = None
    ### iterate
    for i in range(1,n+1):
        # for j <-10 to |seq2|
        for j in range(1,m+1):
            # s_i,j <- max{s_i-1,j ;  s_i,j-1 ;  s_i-1,j-1 (+1 if seq1_i == seq2_v)
            ### Need to do i-1,j-1 here b/c the sequences start:end 1:len+1 -- the 0 row and 0 col are coming into sequences
            possible = [longestPathTo[i-1,j] + Score.getScore(seq1[i-1], "-"), longestPathTo[i,j-1] + Score.getScore("-",seq2[j-1]), longestPathTo[i-1,j-1] + Score.getScore(seq1[i-1], seq2[j-1])]
             #print possible
            indexOfMax = possible.index(max(possible))
            #print indexOfMax
            longestPathTo[i,j] = possible[indexOfMax]
            # backtrack_i,j = 'down' if s_i,j = s_i-1,j ;  right if s_i,j = s_i,j-1 ; diagonal if s_i,j = si-1,j-1 + 1 (if match)
            backtrack[i,j] = directions[indexOfMax] ## can actually just use index of max b/c indexOfMax IS the direction as defined above
    return n, m, longestPathTo, backtrack

def generalizedOutputLCS(backtrack, seq1, seq2, i, j, LCS=""):
    down = 0
    right = 1
    backDirection = backtrack[i,j]
    while i > 0 or j > 0:
##        if i+1 == 0 or j+1 == 0: ##Then it will be all dashes one way or another -- no more in common
##            return LCS, seq1, seq2
        if backtrack[i,j] == down:
            seq2 = seq2[:j] + "-" + seq2[j:]
            backDirection = backtrack[i-1,j]
            if i > 0:
                i = i-1
        elif backDirection == right:
            seq1 = seq1[:i] + "-" + seq1[i:]
            backDirection = backtrack[i,j-1]
            if j > 0:
                j = j-1
        else: #it is diagonal
            LCS = seq1[i-1] + LCS
            ## need i-1 b/c seqs in mtrices start at indices 1,1 -- whereas index of seq starts at 0
            #print LCS
            backDirection = backtrack[i-1,j-1]
            if i > 0:
                i = i-1
            if j > 0:
                j = j-1
    return LCS, seq1, seq2

def generalizedFindLCS(seq1, seq2, Score):
    ''' "Score" is a ScoreMatrix Object '''
    n, m, longestPathTo, backtrack = generalizedLCSinfo(seq1, seq2, Score)
    LCS, seq1, seq2 = generalizedOutputLCS(backtrack, seq1, seq2, n, m)
    return longestPathTo[n,m], LCS, seq1, seq2




### LOCAL
def generalizedLocalAlnInfo(seq1, seq2, Score):
    '''Longest common subsequence between seq1 and seq2
    letters in subsequence not nec. consecutive,but are in order
    "Score" is a ScoreMatrix Object'''
    n = len(seq1)
    m = len(seq2)
    longestPathTo = numpy.zeros([n+1, m+1])
    down = 0 
    right = 1
    diagonal = 2
    freeRide = 3
    backtrack = numpy.zeros([n+1, m+1])
    maxNode = (0,0)
    def maxScore(maxNode):
        return longestPathTo[maxNode[0], maxNode[1]]
    ## Need to set first row and first column to cumulative indels (and backtrack as all downs or all rights)
    directions = [down, freeRide]
    for i in range(1, n+1):
        possible = [longestPathTo[i-1,0] + Score.getScore(seq1[i-1], "-"), 0]
        indexOfMax = possible.index(max(possible))
        longestPathTo[i,0] = possible[indexOfMax]
        if longestPathTo[i,0] > maxScore(maxNode):
            maxNode = (i,0)
        backtrack[i,0] = directions[indexOfMax]

    directions = [right, freeRide]
    for j in range(1, m+1):
        possible = [longestPathTo[0,j-1] + Score.getScore("-",seq2[j-1]), 0]
        indexOfMax = possible.index(max(possible))
        longestPathTo[0,j] = possible[indexOfMax]
        if longestPathTo[i,0] > maxScore(maxNode):
            maxNode = (0,j)
        backtrack[0,j] = directions[indexOfMax]
    backtrack[0,0] = None
    ### iterate
    directions = [down, right, diagonal, freeRide]
    for i in range(1,n+1):
        # for j <-10 to |seq2|
        for j in range(1,m+1):
            # s_i,j <- max{s_i-1,j ;  s_i,j-1 ;  s_i-1,j-1 (+1 if seq1_i == seq2_v)
            ### Need to do i-1,j-1 here b/c the sequences start:end 1:len+1 -- the 0 row and 0 col are coming into sequences
            possible = [longestPathTo[i-1,j] + Score.getScore(seq1[i-1], "-"), longestPathTo[i,j-1] + Score.getScore("-",seq2[j-1]), longestPathTo[i-1,j-1] + Score.getScore(seq1[i-1], seq2[j-1]), 0]
             #print possible
            indexOfMax = possible.index(max(possible))
            ## There are 2 ways to get a 0 -- freeRide from src and an adjacent node that is 0 probably due to freeRide
            ## Check -- and make sure freeRide is recorded if it should be
            if possible[indexOfMax] == 0:
                indexOfMax = freeRide
            longestPathTo[i,j] = possible[indexOfMax]
            ## is it maxnode?
            if longestPathTo[i,j] > maxScore(maxNode):
                maxNode = (i,j)
            # backtrack_i,j = 'down' if s_i,j = s_i-1,j ;  right if s_i,j = s_i,j-1 ; diagonal if s_i,j = si-1,j-1 + 1 (if match)
            backtrack[i,j] = directions[indexOfMax] ## can actually just use index of max b/c indexOfMax IS the direction as defined above
    i,j = maxNode
    return i, j, longestPathTo, backtrack, maxScore(maxNode)



def generalizedOutputLocalAln(backtrack, seq1, seq2, i, j, maxLocalAln=""):
    down = 0
    right = 1
    diagonal = 2
    backDirection = backtrack[i,j]
    ## first - the sequences only go up to last part of local aln
    seq1 = seq1[:i]
    seq2 = seq2[:j]
    while i > 0 or j > 0:
        if backtrack[i,j] == down:
            seq2 = seq2[:j] + "-" + seq2[j:]
            backDirection = backtrack[i-1,j]
            if i > 0:
                i = i-1
        elif backDirection == right:
            seq1 = seq1[:i] + "-" + seq1[i:]
            backDirection = backtrack[i,j-1]
            if j > 0:
                j = j-1
        elif backDirection == diagonal:
            maxLocalAln = seq1[i-1] + maxLocalAln
            ## need i-1 b/c seqs in mtrices start at indices 1,1 -- whereas index of seq starts at 0
            backDirection = backtrack[i-1,j-1]
            if i > 0:
                i = i-1
            if j > 0:
                j = j-1
        else: # it is a freeRide
            return maxLocalAln, seq1[i:], seq2[j:], i, j
    return maxLocalAln, seq1, seq2, i, j

def generalizedFindLocalAln(seq1, seq2, Score):
    ''' "Score" is a ScoreMatrix Object '''
    i, j, longestPathTo, backtrack, maxScore = generalizedLocalAlnInfo(seq1, seq2, Score)
    aligned, seq1, seq2, i, j = generalizedOutputLocalAln(backtrack, seq1, seq2, i, j)
    return longestPathTo, backtrack, maxScore, aligned, seq1, seq2




''' FITTING ALIGNMENT PROBLEM '''

def generalizedFittingAln(seq1, seq2):  ## derived from gen Loc
    '''Longest common subsequence between seq1 and seq2
    letters in subsequence not nec. consecutive,but are in order
    "Score" is a ScoreMatrix Object'''
    ## seq1 should be > seq2:
    if len(seq1) < len(seq2):
        print "Seq1 should be > Seq2"
        return   
    n = len(seq1)
    m = len(seq2)
    longestPathTo = numpy.zeros([n+1, m+1])
    down = 0 
    right = 1
    diagonal = 2
    freeRide = 3
    backtrack = numpy.zeros([n+1, m+1])
    maxNode = (0,0)
    def maxScore(maxNode):
        return longestPathTo[maxNode[0], maxNode[1]]
    ## Need to keep first column as all 0s so its possible to travel there via free Ride
    ## Need to do 1st row as normal -- cumulative indels
    for j in range(1, m+1):
        #longestPathTo[0,j] = float("-inf")
        longestPathTo[0,j] = longestPathTo[0,j-1] - 1
    backtrack[0,1:] = 1 ## only way is right 
    backtrack[1:,0] = 3 ## free rides straight down 1st columnn
    backtrack[0,0] = None
    ### iterate
    directions = [down, right, diagonal]
    for i in range(1,n+1):
        # for j <-10 to |seq2| 
        for j in range(1,m+1):
            # s_i,j <- max{s_i-1,j ;  s_i,j-1 ;  s_i-1,j-1 (+1 if seq1_i == seq2_v)
            ### Need to do i-1,j-1 here b/c the sequences start:end 1:len+1 -- the 0 row and 0 col are coming into sequences
            ### Down choice is represented as -inf to avoid choosing and to keep indexes same
            if seq1[i-1] == seq2[j-1]:
                match = 1
            else:
                match = -1
            if j == m+1: ## can only come from side or diagonal -- not above since something above would have gone straight to sink node
                possible = [float("-inf"), longestPathTo[i,j-1] - 1, longestPathTo[i-1,j-1] + (seq1[i-1] == seq2[j-1])]
            else:
                possible = [longestPathTo[i-1,j] - 1, longestPathTo[i,j-1] - 1, longestPathTo[i-1,j-1] + match]
            indexOfMax = possible.index(max(possible))
            longestPathTo[i,j] = possible[indexOfMax]
            ## is it maxnode on final column?
            if longestPathTo[i,j] > maxScore(maxNode) and j == m:
                maxNode = (i,j)
            # backtrack_i,j = 'down' if s_i,j = s_i-1,j ;  right if s_i,j = s_i,j-1 ; diagonal if s_i,j = si-1,j-1 + 1 (if match)
            backtrack[i,j] = directions[indexOfMax] ## can actually just use index of max b/c indexOfMax IS the direction as defined above
    ## will want to return coordinates for max node
    i,j = maxNode
    LP = longestPathObject(longestPathTo, seq1, seq2)
    BT = backtrackObject(backtrack, seq1, seq2)
    CS, alnSeq1, alnSeq2, i, j = generalizedOutputLocalAln(backtrack, seq1, seq2, i, j)
    return int(maxScore(maxNode)), i, j, alnSeq1, alnSeq2, LP, BT
    ## Note: seq1>seq2 and i is where alignment starts in seq1 (j should be 0)
    ## Trying to align ALL of seq2 to substring of seq1
    ## i.e. all of a 20 bp barcode to a 100 bp read
    ## "i" should be 20 based on logic of how exp was designed....



## EDIT DISTANCE
def editDistance(alnSeq1, alnSeq2):
    ''' Edit distance is essentially hamming distance on alnSeqs, which can include "-" (gaps)
        Input: Globally aligned sequences 1 and 2 -- or localAln segments
        Output: edit distance'''
    editDist = 0
    for i in xrange(len(alnSeq1)):
        editDist += (alnSeq1[i] != alnSeq2[i])
    return editDist


## To find minimum changes do -1 for every indel, -1 for ever mismatch, and 0 for every match
## This way the best will always be things that required the least mismatches and indels

## To do this I will use LCSinfo -- originally used to find the longest common subsequence
##  It gave mismatches and indels 0 while giving matches 1
## Example
Score = ScoreMatrix(makeMatchMismatchDict(m=0,mm=-1), indelPenalty=-1)



def minEditDistance(seq1, seq2):
    '''Longest common subsequence between seq1 and seq2
    letters in subsequence not nec. consecutive,but are in order
    "Score" is a ScoreMatrix Object'''
    n = len(seq1)
    m = len(seq2)
    longestPathTo = numpy.zeros([n+1, m+1])
    down = 0 
    right = 1
    diagonal = 2
    directions = [down, right, diagonal]
    backtrack = numpy.zeros([n+1, m+1])
    ## Need to set first row and first column to cumulative indels (and backtrack as all downs or all rights)
    for i in range(1, n+1):
        longestPathTo[i,0] = longestPathTo[i-1,0] - 1
    for j in range(1, m+1):
        longestPathTo[0,j] = longestPathTo[0,j-1] - 1
    backtrack[0,1:] = 1
    backtrack[1:,0] = 0
    backtrack[0,0] = None
    ### iterate
    for i in range(1,n+1):
        # for j <-10 to |seq2|
        for j in range(1,m+1):
            # s_i,j <- max{s_i-1,j ;  s_i,j-1 ;  s_i-1,j-1 (+1 if seq1_i == seq2_v)
            ### Need to do i-1,j-1 here b/c the sequences start:end 1:len+1 -- the 0 row and 0 col are coming into sequences
            possible = [longestPathTo[i-1,j] + - 1, longestPathTo[i,j-1] - 1, longestPathTo[i-1,j-1] - (seq1[i-1] != seq2[j-1])] ## will subtract 1 if mismatch, 0 if match
             #print possible
            indexOfMax = possible.index(max(possible))
            #print indexOfMax
            longestPathTo[i,j] = possible[indexOfMax]
            # backtrack_i,j = 'down' if s_i,j = s_i-1,j ;  right if s_i,j = s_i,j-1 ; diagonal if s_i,j = si-1,j-1 + 1 (if match)
            backtrack[i,j] = directions[indexOfMax] ## can actually just use index of max b/c indexOfMax IS the direction as defined above
##    print backtrack
##    print seq1
##    print seq2
##    print n
##    print m
    LCS, alnSeq1, alnSeq2 = generalizedOutputLCS(backtrack, seq1, seq2, n, m, LCS="")
##    print LCS
##    print alnSeq1
##    print alnSeq2
    editDist = editDistance(alnSeq1, alnSeq2)
    return alnSeq1, alnSeq2, editDist
