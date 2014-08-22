#!/usr/bin/env python3

import signal
import argparse

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from BioExt.io._SamBamIO import _to_seqrecord
from BioExt.misc import gapful

import pysam

thresholdNum = .05

# We need this because alitnedPairs may end
# with something of the form (x, None)
def getEndingRefIndex(alignedPairs): 
    index = len(alignedPairs)-1
    while(alignedPairs[index][1] == None):
        index -= 1
    return index

# Different alignments align to different ending points in the reference
# Still, there really should be a much better way to do this!
# (instead of looping through all reads)
# This step adds significant computation time
# (we add 1 because of 0-based indexing)
def getRefLength(samfile):
   
    maxLength = 0    
    
    for read in samfile.fetch():
        index = len(read.aligned_pairs)-1
        while(read.aligned_pairs[index][1] == None):
            index -= 1
        if(read.aligned_pairs[index][1] > maxLength):
            maxLength = read.aligned_pairs[index][1]
    return maxLength + 1


def main(bam_file, out_handle, threshold, insertGroups, keepGaps):

    alignedStrings = [] #2d array holding alignedStringsRow's

    try:
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    except ValueError:
        pass


    try:
        # Index bam file in order to samfile.fetch
        pysam.index(bam_file)
        samfile = pysam.Samfile(bam_file, 'rb')
        print("Reading file: " + samfile.filename.decode('ascii'))

        maxLength = None # this will be updated to a matrix of 0's with a length of len(ref) + 1
        # might be an easy way to access this with pysam, I couldn't find it
        # We add 1 because the first column holds what matches to opening gaps in the reference

        for read in samfile.fetch():

            query = read.query.decode('ascii')
            alignedPairs = read.aligned_pairs
            alignedStringsRow = ['']
            
            if(maxLength == None):
                maxLength = [0] * (getRefLength(samfile)+1)

            prevIndex = 0
            nextIndex = 0

            while(alignedPairs[prevIndex][1] == None):
                prevIndex += 1
                alignedStringsRow[0] += query[prevIndex-1]

            nextIndex = prevIndex + 1

            if(alignedPairs[0][1] != None):
                for i in range(alignedPairs[0][1]):
                    alignedStringsRow.append('-')

            while( prevIndex <= getEndingRefIndex(alignedPairs) ): 

                if(alignedPairs[prevIndex][0] != None):
                    alignedStringsRow.append(query[alignedPairs[prevIndex][0]:alignedPairs[prevIndex][0]+1])

                else:
                    alignedStringsRow.append("-")

                if(prevIndex == getEndingRefIndex(alignedPairs)):
                    for j in range((len(alignedPairs)-1-getEndingRefIndex(alignedPairs))):
                        if(alignedPairs[prevIndex+j+1] == None):
                            alignedStringsRow[len(alignedStringsRow)-1] += '-'
                        else:
                            alignedStringsRow[len(alignedStringsRow)-1] += query[alignedPairs[prevIndex+j+1][0]:alignedPairs[prevIndex+j+1][0]+1]
                    break
                
                while alignedPairs[nextIndex][1] == None:
                    nextIndex += 1

                for j in range(nextIndex - prevIndex - 1):
                    alignedStringsRow[len(alignedStringsRow)-1] += query[alignedPairs[prevIndex+1+j][0]:alignedPairs[prevIndex+1+j][0]+1]

                prevIndex = nextIndex
                nextIndex += 1

            while(len(alignedStringsRow) < len(maxLength)):
                alignedStringsRow.append("-")

            for i in range(len(alignedStringsRow)):
                if(len(alignedStringsRow[i]) > maxLength[i]):
                    maxLength[i] = len(alignedStringsRow[i])
            
            alignedStrings.append(alignedStringsRow)
        
        for i in range(len(alignedStrings)):
            for j in range(len(alignedStrings[i])):
                if(j==0):
                    print(alignedStrings[i][j])
                alignedStrings[i][j] += '-' * (maxLength[j]-len(alignedStrings[i][j]))
                if(insertGroups):
                    while(len(alignedStrings[i][j]) % 3 != 1):
                        alignedStrings[i][j] += '-'
                    #at this point the string we are working with is of length 1 mod 3
                    #we are not done, any group of 3 with nucleotides cannot have gaps
                    #so we look at groups of 3 (starting at position 1) and replace gaps 
                    #with nucleotides
                    newstr = alignedStrings[i][j][0:1]
                    for k in range(int((len(alignedStrings[i][j])-1)/3)):
                        s = alignedStrings[i][j][3*k+1:3*k+4]
                        if("A" in s or "C" in s or "G" in s or "T" in s):
                            st = s.replace('-', 'N')
                        newstr += st
                    #print(newstr)
        
    	#We are ignoring any characters not A, C, G, T, -, or N
	#Since we need to deal with a special case if we have majority gap, we have to use
	#a longer method
	#There is still almost certainly a better way to do this than I implemented
        consensus = ""
        
        for currentCol in range(len(alignedStrings[0])): 
            #print(str(currentCol) + " " + str(len(alignedStrings[0][currentCol])))
            for i in range(len(alignedStrings[0][currentCol])): #position in string
                a_count = c_count = g_count = t_count = gap_count = n_count = 0
                for j in range(len(alignedStrings)): #row in array
                    if(alignedStrings[j][currentCol][i] == 'A'):
                        a_count += 1
                    elif(alignedStrings[j][currentCol][i] == 'C'):
                        c_count += 1
                    elif(alignedStrings[j][currentCol][i] == 'G'):
                        g_count += 1
                    elif(alignedStrings[j][currentCol][i] == 'T'):
                        t_count += 1
                    elif(alignedStrings[j][currentCol][i] == '-'):
                        gap_count += 1
                    elif(alignedStrings[j][currentCol][i] == 'N'):
                        n_count += 1

                if (a_count >= c_count and a_count >= g_count and a_count >= t_count and a_count >= gap_count and a_count >= n_count): 
                    consensus += "A"
                elif (c_count >= g_count and c_count >= t_count and c_count >= gap_count and c_count >= n_count): 
                    consensus += "C"
                elif (g_count >= t_count and g_count >= gap_count and g_count >= n_count):
                    consensus += "G"
                elif (t_count >= gap_count and t_count >= n_count):
                    consensus += "T"
                elif (gap_count >= n_count):
                    if(threshold):
                        if( (gap_count / (gap_count + a_count + c_count + g_count + t_count) ) < (1 - thresholdNum) ):
                            consensus += ['A', 'C', 'G', 'T', 'N'][[a_count, c_count, g_count, t_count, n_count].index(max([a_count, c_count, g_count, t_count, n_count]))]
                            print("threshold reached")
                        else:
                            consensus += '-'
                    elif(not threshold and i == 0):
                        consensus += "*"
                    else:
                        consensus += "-"
                else:
                    consensus += 'N'
                    print('N added' + str(a_count) +" "+str(c_count)+" "+str(g_count)+" "+str(t_count)+" "+str(n_count))
        for x in range(14):
            print(len(alignedStrings[0][x]))
        print(consensus[0:100])

        if(keepGaps):
            SeqIO.write( SeqRecord(Seq(consensus),
                               id="", name="consensus", description="consensus"),
                      out_handle, 'fasta')
        
        else:
            #code to strip consensus of gap
            gapless_consensus = ""
            for i in range(len(consensus)):
                if(consensus[i:i+1] != "-"):
                    gapless_consensus += consensus[i:i+1]
           
            SeqIO.write( SeqRecord(Seq(gapless_consensus),
                               id="", name="consensus", description="gap stripped consensus"),
                        out_handle, 'fasta')

        #code to write alignment to output
        '''for currentRow in range(len(alignedStrings)):
            output_seq = ""
            for currentCol in range(len(alignedStrings[0])):
                output_seq += alignedStrings[currentRow][currentCol]
            SeqIO.write( SeqRecord(Seq(output_seq), 
                                   id=str(currentRow), name="query "+str(currentRow), description = "query "+str(currentRow)),
                         out_handle, 'fasta')
            print(currentRow)'''
            

    finally:
        if samfile is not None:
            samfile.close()

    return 0


if __name__ == '__main__':
    import sys

    parser = argparse.ArgumentParser(
        description='get a multiple sequence alignment'
        )

    parser.add_argument(
        'input',
        metavar='INPUT',
        type=argparse.FileType('rb'),
        help='input BAM file'
        )
    parser.add_argument(
        'output',
        default=sys.stdout,
        metavar='OUTPUT',
        type=argparse.FileType('w'),
        help='output FASTA MSA file'
        )
    parser.add_argument(
        '-t', '--threshold',
        action='store_true',
        help='Use a threshold of 5 percent. In tests this does not keep frame, but it is included for historical reasons.'
        )
    parser.add_argument(
        '-ig', '--insertGroups', 
        action='store_true', 
        help='Insertions must be in groups of 3'
        )
    parser.add_argument(
        '-kg', '--keepGaps', 
        action='store_true',
        help='Print consensus with gaps. (Default is to gap-strip).'
        )

    args = None
    retcode = -1
    try:
        args = parser.parse_args()
        bam_file = args.input.name
        args.input.close()
        retcode = main(bam_file, args.output, args.threshold, args.insertGroups, args.keepGaps)
    finally:
        if args is not None:
            if args.output != sys.stdout:
                args.output.close()

    sys.exit(retcode)
