#!/usr/bin/env python
import argparse as ap
import math
import sys
import traceback
from Bio import pairwise2
from Bio.Seq import Seq

__author__ = 'Sami Pietila (sampie@iki.fi)'
__version__ = '0.5'
__date__ = '25th Aug 2013'

def read_params(args):
    p = ap.ArgumentParser( description= 
                           
            "DESCRIPTION\n"
            "  merge_r1r2.py version "+__version__+" ("+__date__+")\n" 
            
            "  merge_r1r2 is a tool to merge paired reads produced by illumina sequencer.\n"
            "  R1 is the forward read and R2 is the reverse read to be merged.\n"
            "  AUTHOR: "+__author__+"\n\n"
            "EXAMPLE\n"
            "  python merge_r1r2.py -i read-I1.fastq -r1 forward-R1.fastq -r2 reverse-R2.fastq -o merged.fastq -oi I1-merged.fastq\n",
            formatter_class=ap.RawTextHelpFormatter )
    arg = p.add_argument

    arg( '-i1', metavar='I1', type=str, nargs='?', default=None, help= 
         "Index file.\n")

    arg( '-r1', metavar='R1', type=str, nargs='?', default=None, help= 
         "Forward reads file.\n")
        
    arg( '-r2', metavar='R2', type=str, nargs='?', default=None, help= 
         "Reverse reads file.\n")
    
    arg( '-o', metavar="Output", type=str, nargs='?', default=None, help = 
         "The output file for merged R1 and R2 reads\n")

    arg( '-oi', metavar="Output index", type=str, nargs='?', default=None, help = 
         "The output file indexes for merged R1 and R2 files.\n")
    
    arg( '-m', metavar="Min overlap", type=str, nargs='?', default=8, help = 
         "Minimum ovelap that must be found. Default 8 chars\n")

    arg( '--revcomp', metavar='True/False', type=bool, nargs='?', default=True, help= 
         "Reverse complement R2 reads. Apply this for R2\n" 
         "reads if they are in reverse complement form\n"
         "so that they can be compared with R1 reads. Default = true\n")
    
    arg( '--include', metavar='treatment', type=str, nargs='?', default="discard", help= 
         "Treatment for reads that can not be merged.\n"
         "discard:  Discard non-mergeable reads\n"
         "R1:       Use R1.\n"
         "R2:       Use R2. Note that --revcomp affects the written R2.\n"
         "hiqa:     Use either R1 or R2 depending on which has\n" 
         "          a higher QA score.\n")
    
    return vars(p.parse_args())

# Calculate sum of QA sequence
def sumQA(R):
    tot = 0
    for q in R["Qa"]:
        tot += int(q)
    return tot    

# Main loop to read input files and write output files
def processFiles(I1FName, R1FName, R2FName, ROFName, IOFName, minOverlap, revCompR2, include):
    
    I1Fh = open(I1FName)
    R1Fh = open(R1FName)
    R2Fh = open(R2FName)
    ROFh = open(ROFName, "w")
    IOFh = open(IOFName, "w")
    
    done = False
    
    reads = 0
    overlappingReads = 0
    overlapLengthAverage = 0

    fastqFields = ["Id", "Seq", "Plus", "Qa"]

    while done == False:
        I1 = {}
        for field in fastqFields:
            lval = I1Fh.readline().strip('\n')
            if (lval == None or lval == ""):
                done = True
                break
            I1[field] = lval

        R1 = {}
        for field in fastqFields:
            lval = R1Fh.readline().strip('\n')
            if (lval == None or lval == ""):
                done = True
                break
            R1[field] = lval
        R2 = {}
        for field in fastqFields:
            lval = R2Fh.readline().strip('\n')
            if (lval == None or lval == ""):
                done = True
                break
            R2[field] = lval
        
        if done == False:
            reads += 1
            overlapFound, overlapLength, r12id, r12Seq, r12Plus, r12Qa = processARead(I1, R1, R2, minOverlap, revCompR2)
            if overlapFound == True:
                overlappingReads += 1
                overlapLengthAverage += overlapLength
                ROFh.write(r12id+"\n")
                ROFh.write(r12Seq+"\n")
                ROFh.write(r12Plus+"\n")
                ROFh.write(r12Qa+"\n")
                IOFh.write(I1["Id"]+"\n")
                IOFh.write(I1["Seq"]+"\n")
                IOFh.write(I1["Plus"]+"\n")
                IOFh.write(I1["Qa"]+"\n")
            else:
                if include == "hiqa": 
                    if sumQA(R2["Qa"]) > sumQA(R1["Qa"]):
                        if revCompR2 == True:
                            r2 = str(Seq(R2["Seq"]).reverse_complement())
                            r2q = R2["Qa"][::-1]
                        else:
                            r2 = R2["Seq"]
                            r2q = R2["Qa"]
                            
                        IOFh.write(R2["Id"]+"\n")
                        IOFh.write(r2+"\n")
                        IOFh.write(R2["Plus"]+"\n")
                        IOFh.write(r2q+"\n")   
                    else:
                        IOFh.write(R1["Id"]+"\n")
                        IOFh.write(R1["Seq"]+"\n")
                        IOFh.write(R1["Plus"]+"\n")
                        IOFh.write(R1["Qa"]+"\n")   

                    IOFh.write(I1["Id"]+"\n")
                    IOFh.write(I1["Seq"]+"\n")
                    IOFh.write(I1["Plus"]+"\n")
                    IOFh.write(I1["Qa"]+"\n")

                elif include == "R1":
                    
                    IOFh.write(R1["Id"]+"\n")
                    IOFh.write(R1["Seq"]+"\n")
                    IOFh.write(R1["Plus"]+"\n")
                    IOFh.write(R1["Qa"]+"\n")   

                    IOFh.write(I1["Id"]+"\n")
                    IOFh.write(I1["Seq"]+"\n")
                    IOFh.write(I1["Plus"]+"\n")
                    IOFh.write(I1["Qa"]+"\n")

                elif include == "R2":
                   
                    if revCompR2 == True:
                        r2 = str(Seq(R2["Seq"]).reverse_complement())
                        r2q = R2["Qa"][::-1]
                    else:
                        r2 = R2["Seq"]
                        r2q = R2["Qa"]
                        
                    IOFh.write(R2["Id"]+"\n")
                    IOFh.write(r2+"\n")
                    IOFh.write(R2["Plus"]+"\n")
                    IOFh.write(r2q+"\n")   

                    IOFh.write(I1["Id"]+"\n")
                    IOFh.write(I1["Seq"]+"\n")
                    IOFh.write(I1["Plus"]+"\n")
                    IOFh.write(I1["Qa"]+"\n")
                        

    IOFh.close()        
    ROFh.close()
    I1Fh.close()
    R1Fh.close()
    R2Fh.close()
    
    if (overlappingReads > 0):
        overlapLengthAverage /= overlappingReads

    return reads, overlappingReads, overlapLengthAverage

# alignScore = {('A','A'):1, ('A', 'C'):0, ('A', 'G'):0, ('A', 'T'):0, ('A', 'N'):0, 
#               ('C','A'):0, ('C', 'C'):1, ('C', 'G'):0, ('C', 'T'):0, ('C', 'N'):0,
#               ('G','A'):0, ('G', 'C'):0, ('G', 'G'):1, ('G', 'T'):0, ('G', 'N'):0,
#               ('T','A'):0, ('T', 'C'):0, ('T', 'G'):0, ('T', 'T'):1, ('T', 'N'):0,
#               ('N', 'N'):0}
# 
# def align(r1olap, r2olap, r1QaOlap, r2QaOlap):
#     alignment = pairwise2.align.globalds(r1olap, r2olap, alignScore,-1,-.5)
#     score = alignment[0][2]
#     return score


# Merge using nucleotides with best quality score
def merge(r1olap, r2olap, r1QaOlap, r2QaOlap):
    assert len(r1olap) == len(r2olap) == len(r1QaOlap) == len(r2QaOlap)
    olap = ""
    olapQa = ""
    for i in range(len(r1olap)):
        R1NTScore = ord(r1QaOlap[i:i+1])
        R2NTScore = ord(r2QaOlap[i:i+1])
        if (R2NTScore > R1NTScore):
            olap += r2olap[i:i+1]
            olapQa += r2QaOlap[i:i+1]
        else:
            olap += r1olap[i:i+1]
            olapQa += r1QaOlap[i:i+1]
    return (olap, olapQa)
    
# Calculate score for overlap. Score grows with probabilities of the given NT pair being true for matching pairs and lowers
# with probabilities of the given NT pair being true for non-matching pairs. Effect of ordering is not removed.
def overlapScore(r1ol, r1pol, r2ol, r2pol):
    assert len(r1ol) == len(r1pol) == len(r2ol) == len(r2pol)
    score = 0
    n = len(r1ol)
    for i in range(n):
        r1nuc = r1ol[i]
        r2nuc = r2ol[i]
        r1p = r1pol[i]
        r2p = r2pol[i]
        # In case of N, QA is very low to make big difference
        if r1nuc == r2nuc:
            score += r1p * r2p
        else:
            score -= r1p * r2p
    return score
    
# Convert each Q to probability nucleotide being true    
def qualityP(QString):
    PString = []
    for q in QString:
        PString.append(1-math.pow(10,-((ord(q)-33) / 10)))
    return PString
    
# Match strings trying with different overlap values to find best overlap    
def findOverlap(r1, r2, r1q, r2q, minOverlap):
    assert len(r1) == len(r2) == len(r1q) == len(r1q)
    tlen = len(r1)
    bestScore = 0
    bestIndex = 0

    r1p = qualityP(r1q)
    r2p = qualityP(r2q)

    for i in range(minOverlap, tlen/2):
        r1ol = r1[tlen-i:]
        r2ol = r2[:i]
        r1pol = r1p[tlen-i:]
        r2pol = r2p[:i]
        
        score = overlapScore(r1ol, r1pol, r2ol, r2pol)
        if (score > bestScore):
            bestScore = score
            bestIndex = i
            
    return bestIndex  

# Process a read. If overlap is found return read details 
# so that it can be written to a file.
def processARead(I1, R1, R2, minOverlap, revCompR2):
    
    r1 = R1["Seq"]
    r1q = R1["Qa"]
    if revCompR2 == True:
        r2 = str(Seq(R2["Seq"]).reverse_complement())
        r2q = R2["Qa"][::-1]
    else:
        r2 = R2["Seq"]
        r2q = R2["Qa"]
    
    overlapLength = findOverlap(r1, r2, r1q, r2q, minOverlap)

    if overlapLength <= 0: # No overlap found
        return False, None, None, None, None, None
    
    r1olap = r1[-overlapLength:]
    r2olap = r2[:overlapLength]
    
    r1QaOlap = r1q[-overlapLength:]
    r2QaOlap = r2q[:overlapLength]
    
    overlap, overlapQa = merge(r1olap, r2olap, r1QaOlap, r2QaOlap)
    
    r12id = R1["Id"] # Just use R1 header
    r12Seq = r1[0:-overlapLength] + overlap + r2[overlapLength:]
    r12Plus = "+" # Just use + in third Line
    r12Qa = r1q[0:-overlapLength] + overlapQa + r2q[overlapLength:]
    return True, overlapLength, r12id, r12Seq, r12Plus, r12Qa


if __name__=="__main__":
    pars = read_params( sys.argv )

    if pars['i1'] is None:
        raise Exception("Please use -i1 to give index file.")
    
    if pars['r1'] is None:
        raise Exception("Please use -r1 to give forward reads file.")
    
    if pars['r2'] is None:
        raise Exception("Please use -r2 to give reverse reads file.")

    if pars['o'] is None:
        raise Exception("Please use -o to give name of merged read file to be written.")

    if pars['oi'] is None:
        raise Exception("Please use -oi to give name of index file for merged reads to be written.")
    
    if pars['include'] not in ["hiqa", "R1", "R2", "discard"]:
        raise Exception("Use hiqa, R1, R2 or discard with include parameter.")
        

    reads, overlappingReads, overlapLengthAverage = processFiles(pars['i1'], pars['r1'], pars['r2'], pars['o'], pars['oi'], pars['m'], pars['revcomp'], pars['include'])
    
    print "Total ", reads, " reads, from which ", overlappingReads, " are overlapping. Average overlap length is ", overlapLengthAverage 
    
try:
    None    

except Exception as exception:
    sys.stderr.write("Error: {0}\n".format(str(exception)))
    for frame in traceback.extract_tb(sys.exc_info()[2]):
        fname,lineno,fn,text = frame
        sys.stderr.write("Error in %s on line %d\n" % (fname, lineno))
    
    
    
    
    
    
    
    
    
    
    