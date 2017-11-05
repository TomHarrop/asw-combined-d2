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

run_log = ('printf "date,branch,hash\\n%s,%s,%s\\n" '
            '"$(date)" '
            '"$(git rev-parse --abbrev-ref HEAD)" '
            '"$(git rev-parse HEAD)" '
            '&>> {log.run} ; ')


#########
# SETUP #
#########

# get a list of fastq files
read_dir_files = list((dirpath, filenames)
                      for (dirpath, dirnames, filenames)
                      in os.walk(read_dir, followlinks=True))

all_fastq_files = []

for dirpath, filenames in read_dir_files:
    for filename in filenames:
        if 'fastq.gz' in filename:
            all_fastq_files.append(os.path.join(dirpath, filename))

sorted(all_fastq_files)

# separate list of R1 and R2
r1_files = sorted(list(x for x in all_fastq_files if '_1.fastq.gz' in x))
r2_files = sorted(list(x for x in all_fastq_files if '_2.fastq.gz' in x))


#########
# RULES #
#########

rule target:
    input:
        'output/trim_decon/reads.fastq.gz'

rule trim_decon:
    input:
        fq = sorted(all_fastq_files)
    output:
        fq = 'output/trim_decon/reads.fastq.gz',
        repair_singles = temp('output/trim_decon/singles.fastq.gz'),
        decon_stats = 'output/trim_decon/decon_stats.txt',
        trim_stats = 'output/trim_decon/trim_stats.txt',
        trim_lhist = 'output/trim_decon/lhist.txt'
    log:
        run = 'logs/trim_decon.run',
        repair = 'logs/repair.log',
        decon = 'logs/decon.log',
        trim = 'logs/trim.log'
    threads:
        25
    shell:
        run_log +
        'echo \''
        'zcat {input.fq}'
        ' | '
        'bin/bbmap/repair.sh '
        'in=stdin.fq '
        'outs={output.repair_singles} '
        'out=stdout.fq '
        'repair '
        '2> {log.repair}'
        ' | '
        'bin/bbmap/bbduk.sh '
        'threads={threads} '
        'in=stdin.fq '
        'out=stdout.fastq '
        'ref=bin/bbmap/resources/phix174_ill.ref.fa.gz '
        'hdist=1 '
        'stats={output.decon_stats} '
        '2> {log.decon}'
        ' | '
        'bin/bbmap/bbduk.sh '
        'threads={threads} '
        'in=stdin.fastq '
        'out={output.fq} '
        'ref=bin/bbmap/resources/adapters.fa '
        'ktrim=r k=23 mink=10 hdist=1 tpe tbo '
        'forcetrimmod=5 '
        'stats={output.trim_stats} '
        'lhist={output.trim_lhist}'
        '2> {log.trim}'
        '\''

