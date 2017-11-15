#!/usr/bin/env python3

import argparse
from Bio import SeqIO

def main():
    # parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta_file')
    args = vars(parser.parse_args())

    # open fasta file
    handle = SeqIO.parse(args['fasta_file'], 'fasta')

    # print tsv to stdout
    print("record_name\tlength")
    for x in handle:
        print('{0}\t{1}'.format(x.id, len(x.seq)))

if __name__ == '__main__':
    main()
