#!/usr/bin/env python

'''
create bed file with positions and sequence identifiers   “<chrom>:<start>-<end>”.
'''


import os
import sys
import errno
import logging
import argparse
import subprocess
from itertools import groupby


log_level = logging.DEBUG
logging.basicConfig(level=log_level,
                    format='[%(asctime)s] - %(''levelname)s - [%(funcName)s] - %(message)s',
                    datefmt='%Y-%m-%d %I:%M:%S %p'
                    )


def parse_args():
    """command line options"""
    parser = argparse.ArgumentParser(
        prog="get_fasta_seq.py",
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    required_group = parser.add_argument_group('required arguments')
    parser.add_argument('--prefix', metavar="<str>", type=str,
                        required=True,
                        dest="prefix",
                        help="prefix for the output fasta and bed files"
                        )
    parser.add_argument('--fasta', metavar="<FILE>", type=str,
                        required=True,
                        dest="fasta",
                        help="path to the fasta file"
                        )
    parser.add_argument('--start', metavar="<int>", type=int,
                        required=True,
                        dest="start",
                        help="start position to subset the sequence"
                        )
    parser.add_argument('--end', metavar="<int>", type=int,
                        required=True,
                        dest="end",
                        help="end position to subset the sequence"
                        )
    parser.add_argument('--nthreads', metavar="<int>", type=int,
                        required=False,
                        default=4,
                        dest="nthreads",
                        help="number of threads/cpus to use"
                        )
    parser.add_argument('--out_dir', metavar="<DIR>", type=str,
                        required=True,
                        dest="out_dir",
                        help="path to the output directory"
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


def extract_sequences(fasta_file, prefix, start, end, out_dir):
    """

    :param fasta_file:
    :param prefix:
    :param start:
    :param end: 
    :return:
    """

    logging.info("fasta file: {}".format(fasta_file))


    fasta_dic = dict((x[0].split()[0], x[1]) for x in fasta_iterator(fasta_name=fasta_file))
    fasta_seq = ''.join(fasta_dic.values())

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    bed_file = os.path.join(out_dir, prefix+str(start)+'-'+str(end))+'.bed'
    fasta_out = os.path.join(out_dir, prefix+str(start)+'-'+str(end))+'.fasta'
    with open(bed_file, 'w') as bed:
        for seq_id, seq in fasta_dic.items():
            bed.write(("{}\t{}\t{}".format(seq_id, str(start), str(end))))
            bed.write('\n')
    logging.info("bed file written to: {}".format(bed_file))

    # use bedtools getfasta to extract sequences from the FASTA file for each of the intervals defined in the bed file
    call = ["bedtools getfasta -fi {} -bed {} -fo {}".format(fasta_file, bed_file, fasta_out)]
    cmd = " ".join(call)
    try:
        subprocess.check_call("set -euo pipefail; " + cmd, shell=True, universal_newlines=True, executable="/bin/bash")
    except (subprocess.CalledProcessError, OSError) as error:
        rc = error.returncode
        if rc == 127:
            extra = "Are you sure this program is installed?"
        else:
            extra = " "
            print("Error occurred: shell exited with return code: {}\ncommand running: {}\n{}".format(error.returncode, cmd, extra), file=sys.stderr)
    logging.info("fasta file written to: {}".format(fasta_out))
    return bed_file, fasta_out

def generate_alignment(fasta_seqs, nthreads, out_dir):
    """
    :param fasta_seqs: 
    """
    split_name = os.path.splitext(os.path.basename(fasta_seqs))[0]
    align_fname = os.path.join(out_dir, split_name)+'.align.fasta'
    log_fname = os.path.join(out_dir, split_name)+'.align.log'
    
    call = ["mafft --reorder --anysymbol --nomemsave --adjustdirection --thread {} {} 1> {} 2> {}".format(nthreads, fasta_seqs, align_fname, log_fname)]
    cmd = " ".join(call)
    try:
        subprocess.check_call("set -euo pipefail; " + cmd, shell=True, universal_newlines=True, executable="/bin/bash")
    except (subprocess.CalledProcessError, OSError) as error:
        rc = error.returncode
        if rc == 127:
            extra = "Are you sure this program is installed?"
        else:
            extra = " "
            print("Error occurred: shell exited with return code: {}\ncommand running: {}\n{}".format(error.returncode, cmd, extra), file=sys.stderr)
    logging.info("alignment file written to: {}".format(align_fname))
    logging.info("log file written to: {}".format(log_fname))



def main():
    parser = parse_args()
    args = parser.parse_args()
    # open input file:
    if args.fasta is not None:
        fasta_file = args.fasta
    else:
        print("Please specify the FASTA file!\n")
        sys.exit(2)
    
    bed, fasta = extract_sequences(fasta_file=fasta_file, prefix=args.prefix, start=args.start, end=args.end, out_dir=args.out_dir)
    generate_alignment(fasta_seqs=fasta, nthreads=args.nthreads, out_dir=args.out_dir)


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
