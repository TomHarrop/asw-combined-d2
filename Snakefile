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


def reformat_input(wildcards):
    r1 = [x for x in all_fastq_files
          if wildcards.pe in x and '_1' in x]
    r2 = [x for x in all_fastq_files
          if wildcards.pe in x and '_2' in x]
    return {'r1': r1, 'r2': r2}


###########
# GLOBALS #
###########

read_dir = 'data/reads/'
meraculous_config_file = 'src/meraculous_config.txt'
read_set = ['norm', 'trim_decon']
kmer_lengths = ['99']

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


#########
# RULES #
#########

rule target:
    input:
        'output/trim_decon/reads.fastq.gz',
        expand('output/k_{kmer}/norm/normalised.fastq.gz',
               kmer=kmer_lengths)

rule reformat:
    input:
        unpack(reformat_input)
    output:
        fq = temp('output/reformat/{pe}.fastq')
    log:
        run = 'logs/reformat.run',
        reformat = 'logs/reformat.log'
    threads:
        1
    shell:
        run_log +
        'bin/bbmap/reformat.sh '
        'verifypaired verifyinterleaved '
        'in={input.r1} '
        'in2={input.r2} '
        'out={output.fq} '
        '2> {log.reformat}'

rule concatenate:
    input:
        expand('output/reformat/{pe}.fastq',
               pe=['pe100', 'pe150'])
    output:
        temp('output/reformat/reads.fastq')
    shell:
        'cat {input} > {output}'

rule trim_decon:
    input:
        fq = 'output/reformat/reads.fastq'
    output:
        fq = 'output/trim_decon/reads.fastq.gz',
        repair_singles = temp('output/trim_decon/singles.fastq.gz'),
        decon_stats = 'output/trim_decon/decon_stats.txt',
        trim_stats = 'output/trim_decon/trim_stats.txt',
        trim_lhist = 'output/trim_decon/lhist.txt'
    log:
        run = 'logs/trim_decon.run',
        repair = 'logs/repair.log',
        decon = 'logs/bbduk_decon.log',
        trim = 'logs/bbduk_trim.log'
    threads:
        16
    shell:
        run_log +
        'bin/bbmap/repair.sh '
        'in={input.fq} '
        'outs={output.repair_singles} '
        'out=stdout.fq '
        'fixinterleaving repair=f '
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
        'zl=9 '
        'ref=bin/bbmap/resources/adapters.fa '
        'ktrim=r k=23 mink=10 hdist=1 tpe tbo '
        'forcetrimmod=5 '
        'stats={output.trim_stats} '
        'lhist={output.trim_lhist} '
        '2> {log.trim}'

# normalise
rule norm:
    input:
        fq = 'output/trim_decon/reads.fastq.gz'
    output:
        fq = 'output/k_{kmer}/norm/normalised.fastq.gz',
        hist = 'output/k_{kmer}/norm/hist_before.txt',
        hist_out = 'output/k_{kmer}/norm/hist_after.txt',
        peaks = 'output/k_{kmer}/norm/peaks.txt',
    log:
        run = 'logs/norm.run',
        norm = 'logs/bbnorm_k{kmer}.log'
    threads:
        50
    shell:
        run_log +
        'bin/bbmap/bbnorm.sh '
        'in={input.fq} '
        'threads={threads} '
        'out={output.fq} '
        'zl=9 '
        'hist={output.hist} '
        'histout={output.hist_out} '
        'target=50 '
        'min=5 '
        'prefilter ecc '
        'k={wildcards.kmer} '
        'peaks={output.peaks} '
        '2> {log.norm}'
