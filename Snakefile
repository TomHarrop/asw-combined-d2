#!/usr/bin/env python3

import datetime
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


###########
# GLOBALS #
###########

read_dir = 'data/reads/'
meraculous_config_file = 'src/meraculous_config.txt'
read_set = ['norm', 'trim_decon']
k = ['99']

log_code = ('printf "date,branch,hash\n%s,%s,%s\n" '
            '"$(date)" '
            '"$(git rev-parse --abbrev-ref HEAD)" '
            '$(git rev-parse HEAD)" '
            '&>> {log.run} ; ')


#########
# SETUP #
#########

# parse the GIT info
print('Job run at {0}'.format(time_now()))
print('Branch: {0}'.format(get_branch()))
print('Hash: {0}'.format(get_hash()))

# get a list of fastq files
read_dir_files = list((dirpath, filenames)
                      for (dirpath, dirnames, filenames)
                      in os.walk(read_dir, followlinks=True))

all_fastq_files = []

for dirpath, filenames in read_dir_files:
    for filename in filenames:
        if 'fastq.gz' in filename:
            all_fastq_files.append(os.path.join(dirpath, filename))
