#!/usr/bin/env python3

import os
import re
import subprocess
import pathlib


#############
# FUNCTIONS #
#############

def resolve_path(x):
    mypath = pathlib.Path(x).resolve()
    return str(mypath)


def get_branch():
    p = subprocess.Popen(['git', 'rev-parse', '--abbrev-ref', 'HEAD'],
                         stdout=subprocess.PIPE)
    return p.stdout.readline().decode().rstrip('\n')


def get_hash():
    p = subprocess.Popen(['git', 'rev-parse', 'HEAD'],
                         stdout=subprocess.PIPE)
    return p.stdout.readline().decode().rstrip('\n')


###########
# GLOBALS #
###########

read_dir = 'data/reads/'
meraculous_config_file = 'src/meraculous_config.txt'
read_set = ['norm', 'trim_decon']
k = ['99']


#########
# SETUP #
#########

# parse the GIT info
print('git branch: {0}'.format(get_branch()))
print('git hash: {0}'.format(get_hash()))

# get a list of fastq files
read_dir_files = list((dirpath, filenames)
                      for (dirpath, dirnames, filenames)
                      in os.walk(read_dir, followlinks=True))

all_fastq_files = []

for dirpath, filenames in read_dir_files:
    for filename in filenames:
        if 'fastq.gz' in filename:
            all_fastq_files.append(os.path.join(dirpath, filename))
