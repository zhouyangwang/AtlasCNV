# --------------------------------------------------------------------------------
# R script that computes RPKM based on the GATK_DoC *.interval_summary file.
# and adds gene_exon, gender, and midpool
# Called by main atlas_cnv.pl v0, Sep 11, 2018. Ted Chiang
# Copyright 2016-2018, Baylor College of Medicine Human Genome Sequencing Center.
# All rights reserved.
# --------------------------------------------------------------------------------

# This is python3 version writeen by zhouyangwang@outlook.com of Atlas CNV.
# As LICENSE indicates,
# All rights belong to Ted Chiang, Baylor College of Medicine Human Genome Sequencing Center

import re
import sys
import pandas as pd

def rpkm(raw, read_length, out):
    rpkm = raw[raw.columns[:4]]    
    exon_lengths = raw['end'] - raw['start']
    # RPKM=Total exon reads/[Mapped reads(Millions)*Exon length(Kb)] 
    rpkm = rpkm.join(pd.DataFrame(map(lambda sample: raw[sample]/ read_length/ exon_lengths * 1000/ (sum(raw[sample])/read_length/1e6), raw.columns[4:]), index = raw.columns[4:]).T)
    rpkm.to_csv(out, sep = '\t', index = False)
    return rpkm

