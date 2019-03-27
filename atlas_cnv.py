#!/usr/bin/python3 
# --------------------------------------------------------------------------------
# main atlas_cnv.pl v0, Sep 11, 2018. Ted Chiang
# Copyright 2016-2018, Baylor College of Medicine Human Genome Sequencing Center.
# All rights reserved.
# This is python3 version of atlas_cnv written by wangzy@haplox.com. All rights belong to Ted Chiang, original author
# --------------------------------------------------------------------------------

import os,sys
import datetime
import argparse
import re

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config')
    parser.add_argument('--panel', nargs = 1)
    parser.add_argument('--sample', nargs = 1)
    parser.add_argument('--threshold_del', type= float)
    parser.add_argument('--threshold_dup', type= float)
    return parser.parse_args()

options = parse_args()

LOG = open('atlas_cnv.log', 'w')
date = datetime.datetime.now()
print('\n' + date + '#---------------------------', LOG)
pwd = os.path.abspath(os.path.curdir)

# --------------------
# Validate the inputs.
# --------------------

path_where_this_script_lives = sys.argv[0]
path_where_this_script_lives = path_where_this_script_lives.replace('/atlas_cnv.pl', '')

if not os.path.isfile(options.config):
    LOG.write("No config provided: $options{'config'}. Looking for a config file in CWD or in atlas_cnv dir.\n")
    if os.path.isfile('config'):
        options.config = 'config'
        print('Found local config: %s/%s\n' %(pwd, 'config'))
        LOG.write('Found config: %s/%s\n' %(pwd, 'config'))

    elif os.path.isfile(path_where_this_script_lives + '/config'):
        options.config = path_where_this_script_lives + '/config'
        print("Found config in atlas_cnv code dir: %s\n" %(options.config))
        LOG.write('Found config: %s\n' %(options.config))

    else:
        print("No config found, please define one.. Exit.\n")
        LOG.write("No config found, exiting.\n")
        exit()
else:
    LOG.write("config: %s\n" %(options.config))

config = {}
with open(options.config) as file:
    for line in file:
        sep = line.rstrip().split('=')
        config[sep[0]] = sep[1]

# -------------------------------------------------------------------------------
# (1a). Process sample file.
# (1b). Process the panel design file.
# (2).  Create midpool dirs.
# (3).  Copy GATK *.DATA.sample_interval_summary files to which midpool dirs.
# (4).  convert_GATK_DoC.R (GATK_DoC --> RPKM).
# (5).  Create a R rpkm matrix for each midpool.
# (6).  Call CNVS using atlas_cnv.R on each R rpkm matrix.
# (7).  Copy .cnv and .pdf files to the sample dir after making CNV dir. (SPECIFIC TO IPIPE)
# -------------------------------------------------------------------------------
# (1a).

sample, midpool = {}, {}
with open(options.sample) as file:
    for line in file:
        sep = line.rstrip.split('\t')
        sample[sep[0]] = sep[1] + '\t' + sep[2]
        if sep[2] in midpool:
            midpool[sep[2]] += 1
        else:
            midpool[sep[2]] = 1

# (1b).
panel_targets, auto_genes, sex_genes = [], [], []
AUTO_GENES, SEX_GENES = {}, {}
with open(options.panel) as file:
    for line in file:
        if re.search(r'\d', line):
            continue
        exon_coor, gene_exon, callcnv, refseq = line.rstrip.split('\t')
        if callcnv == 'Y':
            panel_targets.append(exon_coor)

        if 'X:' in exon_coor or 'Y:' in exon_coor:
            SEX_GENES[gene_exon] = exon_coor
            sex_genes.append(gene_exon)
        else:
            AUTO_GENES[gene_exon] = exon_coor
            auto_genes.append(gene_exon)

# (2).
for i in midpool.keys():
    if os.system('mkdir ' + i) != 0:
        raise Exception('Failed to create midpool dir: %s %s Maybe already exists? Exit and die.\n' %(i, options.sample))

for sub_sample in sample.keys():
    gender, sub_midpool = sample[sub_sample].split('\t')
    this_sample = config['GATKDIR']
    try:
        id, fclbc = sub_sample.split('_')
        this_sample = this_sample.replace('[FCLBC]', fclbc)
    except:
        continue

    this_sample = this_sample.replace('[SAMPLE_FCLBC]', sub_sample)
    cmd = ' '.join(('cp', this_sample, midpool + '/'))
    if os.system(cmd) != 0:
        LOG.write('Failed to copy %s GATK-DoC file to %s\n' %(sub_sample, midpool))

    cmd = ' '.join((config['RSCRIPT'], config['ATLASCNV'] + '/convert_GATK_DoC.R', midpool + '/' + sample + '.sample_interval_summary',
                options.panel, gender, sub_midpool))
    if os.system(cmd) != 0:
        LOG.write("Failed to convert_GATK_DoC.R on %s.\n" % sub_sample)

# (5)
matrix_header = '\t'.join('RPKM' + panel_targets)
for sub_midpool in midpool.keys():
    OUT = open(sub_midpool + '/RPKM_matrix.' + sub_midpool, 'w')
    OUT.write(matrix_header + '\n')

    files = os.listdir(sub_midpool)
    for file in files:
        rpkm = re.search(r'(.*).rpkm.txt$', file)
        if rpkm:
            sample = rpkm.group(1)
            with open(sub_midpool + '/' + sample + '.rpkm.txt') as IN:
                sub_sample = sub_sample.replace('.sample_interval_summary', '')
                OUT.write(sub_sample)
                for line in IN:
                    if 'Gene_Exon' in line:
                        continue
                    OUT.write('\t%.0f' % (float(line.rstrip.split('\t')[2] )))
                OUT.write('\n')
    OUT.close()

    # (6)
    cmd_threshold = ''
    if options.threshold_dup is not None:
        cmd_threshold += ' --threshold_dup ' + str(options.threshold_dup)
    if options.threshold_del is not None:
        cmd_threshold += ' --threshold_del ' + str(options.threshold_del)

    cmd = ' 'join((config['RSCRIPT'], config['ATLASCNV'] + '/atlas_cnv.R', '--rpkm_file ', sub_midpool + '/' + 'RPKM_matrix.' + sub_midpool,
        '--panel_file ', options.panel, cmd_threshold, ">> atlas_cnv.log"))

    if os.system(cmd) != 0:
        LOG.write("Failed to run atlas_cnv.R on %s/RPKM_matrix.%s.\n" %(sub_midpool, sub_midpool))

# (7)
for sub_sample in sample.keys():
    gender, sub_midpool = sub_sample.split('\t')
    sample_dir = re.sub(r'.*_', 'Sample_', sub_sample)
    if os.system('mkdir ' +  sample_dir + '/CNV') != 0:
        LOG.write("Failed to create CNV subdir in %s. Maybe already exists?\n" %(sample_dir))

    if os.system('cp ' + sub_midpool + '/' + sub_sample + '.cnv* ' + sample_dir + '/CNV/') != 0:
        LOG.write("Failed to copy $midpool/$sample.cnv to $sample_dir/CNV/.\n")
    if os.system('cp ' + sub_midpool + '/' + 'sub_sample' + '.pdf ' + sample_dir + '/CNV/') != 0:
        LOG.write("Failed to copy $midpool/$sample.pdf to $sample_dir/CNV/.\n")






LOG.close()