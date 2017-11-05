#!/usr/bin/env python3

import os
import re
import pathlib


#############
# FUNCTIONS #
#############

def resolve_path(x):
    mypath = pathlib.Path(x).resolve()
    return str(mypath)


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

# get a list of fastq files
read_path = resolve_path(read_dir)
read_dir_files = list((dirpath, filenames)
                      for (dirpath, dirnames, filenames)
                      in os.walk(read_path))

all_fastq_files = []

for dirpath, filenames in read_dir_files:
    for filename in filenames:
        if 'fastq.gz' in filename:
            all_fastq_files.append(os.path.join(dirpath, filename))
