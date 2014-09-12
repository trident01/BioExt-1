#!/usr/bin/env python3

import sys

from Bio import SeqIO

from BioExt.args import (
    add_alphabet,
    add_reference,
    add_scorematrix
    )
from BioExt.io import BamIO
from BioExt.misc import compute_cigar, gapless
from BioExt.scorematrices import (
    DNAScoreMatrix,
    FrequenciesError,
    ProteinScoreMatrix
    )
from BioExt.uds import _align_par
import bealign
import msaConsensus


def toreference(string):
        try:
            with open(string) as handle:
                ref = next(SeqIO.parse(handle, 'fasta'))
            return ref
        except:
            msg = "'{0}' does not exist or is not a valid FASTA file".format(string)
            raise ArgumentTypeError(msg)



def main(
        input_handle,
        output_handle,
        oneAlignment,
        reference,
        expected_identity,
        alphabet,
        reverse_complement,
        score_matrix,
        discard_handle,
        do_sort,
        quiet,
        extendGapPenalty,
        codonMatrix,
        globalStartingPoint,
        threshold,
        insertGroups,
        keepGaps
        ):

    retcode1 = bealign.main(
                input_handle,
                output_handle+"FIRST.bam",
                reference,
                expected_identity,
                alphabet,
                reverse_complement,
                score_matrix,
                discard_handle,
                do_sort,
                quiet,
                codonMatrix,
                globalStartingPoint,
                extendGapPenalty
                )
    retcode2 = retcode3 = 0
    
    if(not oneAlignment):
        retcode2 = msaConsensus.main(output_handle+"FIRST.bam", output_handle+"Consensus.fasta", threshold, insertGroups, keepGaps)

        retcode3 = bealign.main(
                input_handle,
                output_handle+"FINAL.bam",
                toreference(output_handle+"Consensus.fasta"),
                expected_identity,
                alphabet,
                reverse_complement,
                score_matrix,
                discard_handle,
                do_sort,
                quiet,
                codonMatrix,
                globalStartingPoint,
                extendGapPenalty
                )

    return retcode1+retcode2+retcode3


if __name__ == '__main__':
    import argparse

    from os import remove
    from os.path import getsize

    def probability(string):
        try:
            p = float(string)
            if p < 0 or p > 1:
                raise ValueError()
            return p
        except ValueError:
            msg = "'{0}' is not a probability in [0, 1]".format(string)
            raise argparse.ArgumentTypeError(msg)

    parser = argparse.ArgumentParser(
        description=(
            'align sequences to HXB2, then '
            'to a msaConsensus using a codon alignment algorithm. Options are those of bealign and msaCons.'
            ' This outputs 3 files: a .bam which is aligned to the users reference choice, a consensus '
            'based on this alignment, and a .bam aligned to this consensus'
            )
        )

    parser.add_argument(
        'input',
        metavar='INPUT',
        type=argparse.FileType('r'),
        help='INPUT FASTA file with sequences to be pairwise aligned'
        )
    parser.add_argument(
        'output',
        metavar='OUTPUT',
        type=argparse.FileType('wb'),
        help='send files to this output name (Eg. if user enters asdf, alignFinal will output asdfFIRST.bam, asdfConsensus.fasta, and asdfFINAL.bam)'
        )
    parser.add_argument(
        '-oa', '--oneAlignment', 
        action = 'store_true',
        help='Only do one alignment to chosen reference'
        )
    add_reference(parser, '-r', '--reference')
    parser.add_argument(
        '-e', '--expected-identity',
        type=probability,
        default=None,
        help='discard sequences that are insufficiently identical to the reference'
        )
    add_alphabet(parser, '-a', '--alphabet')
    add_scorematrix(parser, '-m', '--score-matrix')
    parser.add_argument(
        '-D', '--discard',
        metavar='DISCARD',
        type=argparse.FileType('w'),
        help='discarded sequences are sent to DISCARD'
        )
    parser.add_argument(
        '-R', '--reverse-complement',
        action='store_true',
        help=(
            "also align the reverse complement of each query sequence, "
            "returning it if the alignment is superior"
            )
        )
    parser.add_argument(
        '-S', '--no-sort',
        dest='sort',
        action='store_false',
        help='do NOT sort the resulting BAM file [the default is to sort]'
        )
    parser.add_argument(
        '-q', '--quiet',
        action='store_true',
        help='do not print status update messages'
        )
    parser.add_argument(
        '-egp', '--extendGapPenalty',
        action='store',
        help='set the extend gap penalty to this percentage of the range of the scoring matrix [default=2.5]' 
        )
    parser.add_argument(
        '-cm', '--codonMatrix',
        action='store_true',
        help='Use an empirical codon matrix for scoring alignments. If selected, overrides any different scoring matrix selected.'
        )
    parser.add_argument(
        '-gsp', '--globalStartingPoint', 
        action='store_true',
        help='Sequences are penalized for not starting at the starting point of the reference (the first row and column of the scoring matrix are initialized with penalties). Sequences are not penalized for ending early, with the caveat that at least one sequence be used fully (the backtrack starts with the max score in the bottom row or rightmost column of the dynamic matrix). This option is therefore different from a global or local alignment.'
        )

    parser.add_argument(
        '-t', '--threshold',
        action='store_true',
        help='Use a threshold of 5 percent when considering whether to take a gap majority in making the consensus. This does not keep frame, but it is included for historical reasons.'
        )
    parser.add_argument(
        '-ig', '--insertGroups', 
        action='store_true', 
        help='Insertions in the consensus must be in groups of 3 (not working/fully implemented).'
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
        output_file = args.output.name
        args.output.close()
        if(args.extendGapPenalty == None):
            args.extendGapPenalty = 2.5
        retcode = main(
            args.input,
            output_file,
            args.oneAlignment,
            args.reference,
            args.expected_identity,
            args.alphabet,
            args.reverse_complement,
            args.score_matrix,
            args.discard,
            args.sort,
            args.quiet,
            args.extendGapPenalty,
            args.codonMatrix,
            args.globalStartingPoint, 
            args.threshold,
            args.insertGroups,
            args.keepGaps
        )
    finally:
        if args is not None:
            # close input file handle
            if args.input and args.input != sys.stdin:
                args.input.close()
            if args.discard and args.discard != sys.stdout:
                args.discard.close()
                try:
                    if getsize(args.discard.name) == 0:
                        remove(args.discard.name)
                except:
                    pass

    sys.exit(retcode)
