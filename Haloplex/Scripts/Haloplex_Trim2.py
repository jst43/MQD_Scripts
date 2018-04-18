#!/usr/bin/python3


"""Despite troubleshooting with Agilent, trimming errors still
persist in the pipeline. As a workaround (read: hack) we remove 1bp from each
read in the trimmed fastq files provided as output from SureCallTRIMMER"""


import os
import sys
import argparse
import gzip
import glob


# FUNCTIONS


def get_fastq_types():
    filename = 'fastq_types.txt'
    fastq_types = {}
    with open(filename, 'r') as f:
        for line in f:
            key, val = line.split(':')
            val = val.strip()
            fastq_types[key] = val
    return fastq_types


def check_for_trimmed_files(fastq_types):
    if not os.path.exists('./trimmed_fastq/'):
        sys.exit('trimmed_fastq directory not found')
    # Read in 'samples.txt' and check for existence of trimmed fastqs
    with open('samples.txt', 'r') as f:
        for line in f:
            line = line.strip()
            path_1 = '{}_{}.trimmed.fastq.gz'.format(line, fastq_types['R1 type'])
            path_2 = '{}_{}.trimmed.fastq.gz'.format(line, fastq_types['R2 type'])
            if not os.path.isfile('./trimmed_fastq/' + path_1):
                sys.exit("Can't find {}".format(path_1))
            if not os.path.isfile('./trimmed_fastq/' + path_2):
                sys.exit("Can't find {}".format(path_2))
    print('Trimmed fastqs found for all samples')


def trim_fastq(fastq_file):
    print('Trimming {}'.format(fastq_file))
    new_fastq_name = fastq_file.split('.f')[0] + '2.fastq.gz'
    with gzip.open(new_fastq_name, 'wt') as outfile:
        with gzip.open(fastq_file, 'rt') as infile:
            for i, line in enumerate(infile):
                if i % 2 == 1:
                    line = line[1:-2] + '\n'
                outfile.write(line)
    print('Done. Finished writing output to {}'.format(new_fastq_name))


def trim_fastqs():
    fastq_names = glob.glob('./*trimmed.fastq.gz')
    for name in fastq_names:
        trim_fastq(name)


def main():
    parser = argparse.ArgumentParser(description='Checks for trimmed fastq files and provides further trimming')
    parser.add_argument('-f', help='Specifies directory to original fastq files', type=str)
    args = parser.parse_args()
    filepath = args.f
    os.chdir(filepath)
    fastq_types = get_fastq_types()
    check_for_trimmed_files(fastq_types)
    os.chdir('./trimmed_fastq/')
    trim_fastqs()


if __name__ == '__main__':
    main()
