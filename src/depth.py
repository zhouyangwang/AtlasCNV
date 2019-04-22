import pysam
import pandas as pd
import os

def depth_combine(samples, bed, batch_out, name_list):
    bed_df = pd.read_csv(bed, header = None, sep = '\t')
    try:
        df = pd.DataFrame({'chr': bed_df[0],
                           'start': bed_df[1], 
                            'end': bed_df[2],
                            'name': bed_df[3]})
    except:
        raise Exception('bed file %s should contain at least 4 columns composed of chr, start, end, the name of the segment' %(bed))

    df = df.join(pd.DataFrame(map(lambda sample: pysam_depth(sample, bed), samples), index = name_list).T)
    df.to_csv(path_or_buf = batch_out, sep = '\t', index = False)
    return df


def pysam_depth(bam, bed):
    "get number of total base in bed region"
    if not os.path.isfile(bam + '.bai'):
        raise Exception('index for BAM file %s isn\'t found' %(bam))

    cmd = [bed, bam]# ,'-Q', bytes(5)]
    try:
        raw = pysam.bedcov(*cmd, split_lines=False)
    except pysam.SamtoolsError as exc:
        raise ValueError("Failed processing %r coverages in %r regions. "
                         "PySAM error: %s" % (bam, bed, exc))

    return map(lambda x: int(x.split('\t')[-1]), raw.rstrip().split('\n'))
