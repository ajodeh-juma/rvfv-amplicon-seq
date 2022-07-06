#!/usr/bin/env python3

'''rename FASTA headers using a specified name (prefix) in masked consensus FASTA file'''

import os
import argparse
import fileinput
from textwrap import dedent


def parse_args():
    """command line options"""
    parser = argparse.ArgumentParser(
        prefix_chars='-',
        description=__doc__
    )

    required_group = parser.add_argument_group(dedent('''mandatory arguments'''))
    required_group.add_argument('--fasta', metavar='<FILE>', required=True,
                                dest='fasta',
                                help="path to the masked consensus FASTA file"
                                )

    required_group.add_argument('--prefix', metavar='str', required=True,
                                dest='prefix',
                                help="prefix of the new sample's header names"
                                )

    return parser


def rename_headers(fasta, prefix):
    """
    param masked_fasta
    param prefix
    return
    """

    bak = os.path.join(os.path.dirname(fasta), os.path.basename(fasta) + '.bak')
    if os.path.exists(bak):
        print("renamed masked consensus file exists")
        pass
    else:
        try:
            print('renaming sequence headers in: {}'.format(fasta))
            with fileinput.FileInput(fasta, inplace=True, backup='.bak') as fn:
                for line in fn:
                    if line.startswith('>'):
                        print('>' + str(prefix), end='\n')
                    else:
                        print(line, end='')
        except Exception as error:
            print(f"Error {error} occurred")
    return fasta


def main():
    parser = parse_args()
    args = parser.parse_args()

    rename_headers(fasta=args.fasta, prefix=args.prefix)


if __name__ == '__main__':
    main()
