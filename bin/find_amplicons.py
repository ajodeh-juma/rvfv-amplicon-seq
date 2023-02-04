#!/usr/bin/env python

'''
identify primer sequences given a consensus fasta file
'''
from __future__ import division

import re
import sys

import logging
import argparse
from itertools import groupby

try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO

log_level = logging.DEBUG
logging.basicConfig(level=log_level,
                    format='[%(asctime)s] - %(''levelname)s - [%(funcName)s] - %(message)s',
                    datefmt='%Y-%m-%d %I:%M:%S %p'
                    )


def parse_args():
    """command line options"""
    parser = argparse.ArgumentParser(
        prog="find_amplicons.py",
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    required_group = parser.add_argument_group('required arguments')
    parser.add_argument('--bed', metavar="<FILE>", type=str,
                        required=True,
                        dest="bed",
                        help="path to the bed file having primer sequences and positions"
                        )
    parser.add_argument('--fasta', metavar="<FILE>", type=str,
                        required=True,
                        dest="fasta",
                        help="path to the fasta file"
                        )
    return parser


def fasta_iterator(fasta_name):
    """
    Given a fasta file. yield tuples of header, sequence
    Author: Brent Pedersen
    https://www.biostars.org/p/710/
    """
    fa_fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fa_fh, lambda line: line[0] == ">"))

    for header in faiter:
        # drop the ">"
        header = next(header)[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq
    fa_fh.close()


def get_amplicons(fasta_file, bed_file):
    """

    :param fasta_file:
    :param bed_file:
    :return:
    """

    logging.info("fasta file: {}".format(fasta_file))
    logging.info("primer file: {}".format(bed_file))


    # get primers
    l_dict = dict()
    r_dict = dict()
    for line in open(bed_file):
        if line.startswith('name'):
            continue
        else:
            l_match = re.match(r'.*_LEFT$', line.strip().split()[0])
            r_match = re.match(r'.*_RIGHT$', line.strip().split()[0])
            if l_match is not None:
                l_dict[l_match.group()] = line.strip().split()[2]
            if r_match is not None:
                r_dict[r_match.group()] = line.strip().split()[2]

    fasta_dic = dict((x[0].split()[0], x[1]) for x in fasta_iterator(fasta_name=fasta_file))
    fasta_seq = ''.join(fasta_dic.values())

    # find primer locations on the + strand
    for primer_id, primer_seq in l_dict.items():
        l_start = [m.start() for m in re.finditer('(?={})'.format(primer_seq), fasta_seq)]
        l_end = [m.start() + len(primer_seq) for m in re.finditer('(?={})'.format(primer_seq), fasta_seq)]
        print("{}\t{}\t{}\t{}".format(primer_id, primer_seq, l_start, l_end))

    # find primer locations on the - strand
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'Y': 'Y', 'R': 'R', 'W': 'W', 'M': 'M', 'D': 'D',
                  'S': 'S', 'K': 'K', 'H': 'H', 'V': 'V'}
    complement_seq = ''.join([complement[base] for base in fasta_seq])
    rev_comp = ''.join([complement[base] for base in fasta_seq[::-1]])

    for primer_id, primer_seq in r_dict.items():
        rev_primer = ''.join([complement[base] for base in primer_seq[::-1]])
        r_end = [len(rev_comp) - m.start() for m in re.finditer('(?={})'.format(primer_seq), rev_comp)]
        r_start = [len(rev_comp) - (m.start() + len(primer_seq) )for m in re.finditer('(?={})'.format(primer_seq), rev_comp)]
        print("{}\t{}\t{}\t{}".format(primer_id, primer_seq, r_start, r_end))


def main():
    parser = parse_args()
    args = parser.parse_args()
    # open input file:
    if args.fasta is not None:
        fasta_file = args.fasta
    else:
        print("Please specify the FASTA file!\n")
        sys.exit(2)
    if args.bed is not None:
        bed_file = args.bed
    else:
        print("Please specify the BED file!\n")
        sys.exit(2)

    # run the function
    # rc = reverse_complement(fasta_file)

    get_amplicons(fasta_file=fasta_file, bed_file=bed_file)


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        logging.error(str(e))
        logging.debug('', exc_info=True)
        try:
            sys.exit(e.errno)
        except AttributeError:
            sys.exit(1)
