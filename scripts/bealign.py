#!/usr/bin/env python3

from __future__ import division, print_function

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


def main(
        input_handle,
        output_handle,
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
        ):

    try:
        score_matrix_ = score_matrix.load()
    except:
        raise RuntimeError('could not load the score matrix')

    if ((alphabet == 'dna' and not isinstance(score_matrix, DNAScoreMatrix)) and
            not isinstance(score_matrix, ProteinScoreMatrix)):
        raise ValueError(
            'DNA alphabet requires a DNA score matrix, '
            'while amino and codon alphabets require a protein score matrix'
            )

    do_codon = alphabet == 'codon'

    records = SeqIO.parse(input_handle, 'fasta')

    # grab the first, make it gapless once and for all
    if reference is None:
        reference = gapless(next(records))
        def allseqs(records):
            yield compute_cigar(reference, reference)
            for record in records:
                yield record
    else:
        def allseqs(records):
            for record in records:
                yield record

    if discard_handle:
        def discard(record):
            SeqIO.write([gapless(record.upper())], discard_handle, 'fasta')
    else:
        discard = None

    def output(records):
        BamIO.write(
            allseqs(records),
            output_handle,
            reference
            )

    retcode = -1
    try:
        _align_par(
            reference,
            records,
            score_matrix_,
            do_codon,
            reverse_complement,
            expected_identity,
            discard,
            output,
            codonMatrix,
            globalStartingPoint,
            extendGapPenalty,
            quiet
            )
        if do_sort:
            BamIO.sort(output_handle)
        retcode = 0
    except FrequenciesError:
        print(
            'supplied score-matrix does not imply a frequency distribution:',
            'please choose another matrix if you must filter on expected identity',
            file=sys.stderr
            )

    return retcode


if __name__ == '__main__':
    import argparse

    from os import remove
    from os.path import getsize

    from BioExt import __version__ as version

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
            'align sequences to a reference using '
            'a codon alignment algorithm and output to a BAM file'
            )
        )

    parser.add_argument(
        'input',
        metavar='INPUT',
        type=argparse.FileType('r'),
        help='INPUT FASTA file'
        )
    parser.add_argument(
        'output',
        metavar='OUTPUT',
        type=argparse.FileType('wb'),
        help='send BAM to OUTPUT'
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
        '-v', '--version',
        action='version',
        version='BioExt version {0}'.format(version),
        help='print version information and exit'
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
        '-egp', '--extendGapPenalty',
        action='store',
        help='set the extend gap penalty to this percentage of the range of the scoring matrix [default=2.5]' 
        )

    args = None
    retcode = -1
    try:
        args = parser.parse_args()
        output_file = args.output.name
        args.output.close()
        retcode = main(
            args.input,
            output_file,
            args.reference,
            args.expected_identity,
            args.alphabet,
            args.reverse_complement,
            args.score_matrix,
            args.discard,
            args.sort,
            args.quiet,
            args.codonMatrix,
            args.globalStartingPoint,
            args.extendGapPenalty
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
