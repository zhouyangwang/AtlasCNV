#!/usr/bin/python3 
# --------------------------------------------------------------------------------
# main atlas_cnv.pl v0, Sep 11, 2018. Ted Chiang
# Copyright 2016-2018, Baylor College of Medicine Human Genome Sequencing Center.
# All rights reserved.
# This is python3 version of atlas_cnv written by wangzy@haplox.com. All rights belong to Ted Chiang, original author
# --------------------------------------------------------------------------------

import os
import sys
import datetime
import argparse
import re
import pysam
from test import test
from depth import depth_combine
from rpkm import rpkm
from call import call

def parse_args():
    parser = argparse.ArgumentParser(description = 'atlas_cnv: batch samples based CNV analysis tool.\nversion=0.1.0')
    parser.add_argument('-s', '--sample', help = 'file contains path of BAM')
    parser.add_argument('-b', '--bed', help = 'bed file')
    parser.add_argument('--threshold_del', type= float, default = -0.6, help = 'log2 no more than threshold_del woule be seen as deletion')
    parser.add_argument('--threshold_dup', type= float, default = 0.4,  help = 'log2 no less than threshold_dup woule be seen as insertion')
    parser.add_argument('-q', '--sampleQC', type = float, default = 0.2, help = 'sample with quality value less than sampleQC will be filterd')
    parser.add_argument('--read_length',help = 'It\'s used when calculate rpkm, single read', type = int , default = 150)
    parser.add_argument('-o', '--outdir', help = 'directory where results and sub sample analysis stored')
    parser.add_argument('--test', action='store_true', help = 'Running with demo data')
    return parser.parse_args()

options = parse_args()
if options.test:
    test()
    exit()
else:
    if options.sample is None or options.bed is None:
        raise Exception('-s/--sample, -b/--bed shoule be provided')

if options.outdir is None:
    options.outdir = os.path.abspath(os.path.curdir)

date = datetime.datetime.now()
sys.stderr.write('Atlas_cnv starts to work at %s\n' % date)

#### file fromat check ###################
""" sex and autochromosome might be useful in the future version
    file: batch name, female, patient name, bam file path
"""
midpool, name_list = {}, []
with open(options.sample) as F1:
    sex_list = ('F', 'M', 'FEMALE', 'MALE')
    for line in F1:
        sep = line.rstrip().split('\t')
        if sep[0] not in midpool:
            midpool[sep[0]] = [sep[3]]
        else:
            midpool[sep[0]].append(sep[3])
        name_list.append(sep[2])

for batch in midpool.keys():
    batch_out = os.path.join(options.outdir,  batch)
    if os.system('mkdir -p ' + batch_out) != 0:
        raise Exception('Failed to create dir for %s: %s\n' %(batch, batch_out))

    sys.stderr.write('Use samtools to get depth of each region for samples: %s\n' %('\n'.join(midpool[batch])))
    sample_depth = depth_combine(midpool[batch], options.bed, os.path.join(batch_out, batch + '.depth.txt'), name_list)
    sample_rpkm = rpkm(sample_depth, options.read_length, os.path.join(batch_out, batch + '.rpkm.txt'))
    sys.stderr.write(f'\nBegin to call CNV, threshold for deletion: {options.threshold_del}, duplication: {options.threshold_dup},sampleQC: {options.sampleQC}\n')
    call(sample_rpkm, options.threshold_del, options.threshold_dup, options.sampleQC, os.path.join(batch_out, batch + '.rpkm.txt'))

sys.stderr.close()

