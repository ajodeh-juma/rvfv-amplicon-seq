#!/usr/bin/env python3

"""
generate samplesheet.csv file

"""

import os
import re
import sys
import csv
from pathlib import Path
import argparse
from os import path
from sys import argv

def parse_args():
    """
    """
    parser = argparse.ArgumentParser(description=__doc__)
    required_group = parser.add_argument_group('input/output options')
    required_group.add_argument('--dataDir', type=Path, dest='data_dir',
                        help="path to the data directory having your reads"
                        )
    required_group.add_argument('--outfile', type=argparse.FileType('w'), dest='outfile', 
                        help='output sample sheet file'
                        )
    return parser

def samplesheet_from_dir(data_dir, outfile):
    """
    """
    ext = 'fastq fq'.split()
    samples_dict = dict()
    if not os.path.exists(data_dir):
        print('data directory {} does not exist'.format(data_dir))
    else:
        for dirname, dirs, filenames in os.walk(os.path.abspath(data_dir)):
            filenames.sort
            for filename in filenames:
                if ext[0] in filename or ext[1] in filename:
                    sample_id = filename.split('.fastq')[0].split('.fq')[0]
                    sample_id = sample_id.replace("_R1", "").replace("_r1", "").replace("_R2", "").replace("_r2", "")
                    sample_id = sample_id.replace(".R1", "_R1").replace(".R2", "").replace("_R1", "").replace("_R2","")
                    sample_id = re.sub(r'_S(\d+)_L\d{3}_\d{3}', '', sample_id)
                    sample_id = sample_id.replace(".", "_").replace(" ", "_")
                    sample_id = re.sub(r'\_.*', '', sample_id)

                    file_path = os.path.join(dirname, filename)
                    if not sample_id in samples_dict:
                        samples_dict[sample_id] = [file_path]
                    else:
                        samples_dict[sample_id].append(file_path)

    header = ['sample','fastq_1','fastq_2']
    try:
        if not os.path.exists(os.path.dirname(outfile.name)):
            os.makedirs(os.path.dirname(outfile.name))
        with open(outfile.name, 'w') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(header)
            for k, v in samples_dict.items():
                writer.writerow([k, sorted(v)[0], sorted(v)[1]])
                
    except IOError:
        print("I/O error")

def main():
    """
    """
    parser = parse_args()
    args = parser.parse_args()
    samplesheet_from_dir(data_dir=args.data_dir, outfile=args.outfile)

if __name__ == '__main__':
    main()

    
