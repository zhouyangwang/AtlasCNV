# --------------------------------------------------------------------------------
# atlas_cnv.R v0. called by main atlas_cnv.pl v0, Sep 11, 2018. Ted Chiang
# Copyright 2016-2018, Baylor College of Medicine Human Genome Sequencing Center.
# All rights reserved.
# --------------------------------------------------------------------------------

# This is python3 version writeen by zhouyangwang@outlook.com of Atlas CNV.
# As LICENSE indicates,
# All rights belong to Ted Chiang, Baylor College of Medicine Human Genome Sequencing Center

import re
import sys
import os
import pandas as pd
import numpy as np
import scipy.stats
import researchpy as rp
import statsmodels.api as sm
from statsmodels.formula.api import ols
  

def call(raw, threshold_del, threshold_dup, threshold_sampleQC, outfile, threshold_sample_anova = 0.05):
    # 95%, 97.5%, 99%, 99.95%, 99.99% => zscore: 1.645, 1.96, 2.576, 3.291, 4

    rpkm = raw[raw.columns[4:]].T
    rpkm.columns = list(map(lambda x, y, z, w: ':'.join((x, str(y), str(z), str(w))),
                          raw['chr'], raw['start'], raw['end'],  raw['name']))

    rpkm_wo_outlier = outliner_remove(rpkm)
    median_wo_rpkm = median_df(rpkm_wo_outlier)

    exon_length = len(rpkm.columns)
    sample_length = len(rpkm.index)
    log2_rpkm = pd.DataFrame(map(lambda target_coor, median : map(lambda x,y: np.log2(x/y) if y != 0 and x != 0 else np.NaN,
                                rpkm[target_coor], [median] * exon_length), rpkm.columns, median_wo_rpkm.median_rpkm), 
                            index= rpkm.columns, columns = rpkm.index).T
    log2_rpkm_finite = np.isfinite(log2_rpkm)

    log2_rpkm_wo_outlier = pd.DataFrame(map(lambda target_coor, median : map(lambda x,y: np.log2(x/y) if y != 0 and x != 0 else np.NaN,
                                rpkm_wo_outlier[target_coor], [median] * exon_length), rpkm.columns, median_wo_rpkm.median_rpkm), 
                            index= rpkm.columns, columns = rpkm.index).T
    log2_rpkm_wo_outlier_finite = np.isfinite(log2_rpkm_wo_outlier)

    #sys.stderr.write(log2_rpkm.head(), log2_rpkm_finite.head(), log2_rpkm_wo_outlier.head())

    ExonQC = pd.DataFrame(map(lambda exon: std(log2_rpkm[exon]), log2_rpkm.columns),
                          index = log2_rpkm.columns, columns= ['std'])
    ExonQC_wo_outliers = pd.DataFrame(map(lambda exon: std(log2_rpkm_wo_outlier[exon]), log2_rpkm_wo_outlier.columns),
                          index = log2_rpkm_wo_outlier.columns, columns = ['std'])

    threshold_del_soft = pd.DataFrame(map(lambda exon: -2.576 * std(log2_rpkm_wo_outlier[exon]) + np.nanmean(log2_rpkm_wo_outlier[exon]), log2_rpkm_wo_outlier.columns),
                          index = log2_rpkm_wo_outlier.columns, columns = ['del'])
    threshold_dup_soft = pd.DataFrame(map(lambda exon:  2.576 * std(log2_rpkm_wo_outlier[exon]) + np.nanmean(log2_rpkm_wo_outlier[exon]), log2_rpkm_wo_outlier.columns),
                          index = log2_rpkm_wo_outlier.columns, columns = ['dup'])

    rpkm_matrix_file =  outfile
    exon_qc_outfile = rpkm_matrix_file.replace('rpkm.txt', 'ExonQC')
    exon_qc_wo_outliers_outfile = rpkm_matrix_file.replace('rpkm.txt', 'ExonQC_wo_outliers')
    threshold_del_soft_outfile = rpkm_matrix_file.replace('rpkm.txt', 'Exon_threshold_del_soft_wo_outliers')
    threshold_dup_soft_outfile = rpkm_matrix_file.replace('rpkm.txt', 'Exon_threshold_dup_soft_wo_outliers')

    ExonQC.to_csv(exon_qc_outfile, sep='\t',  index=True)
    ExonQC_wo_outliers.to_csv(exon_qc_wo_outliers_outfile, sep = '\t', index = True)
    threshold_del_soft.to_csv(threshold_del_soft_outfile, sep = '\t', index = True)
    threshold_dup_soft.to_csv(threshold_dup_soft_outfile, sep = '\t', index = True)

    threshold_exonQC = np.nanmean(ExonQC_wo_outliers) + 3.291 * std(ExonQC_wo_outliers)
    sys.stderr.write(f"\nExonQC threshold based on 99% of sd distribution is: {threshold_exonQC}\n")
    sys.stderr.write(f"CNV Exon threshold del/dup: {threshold_del}, {threshold_dup}\n")

    SampleQC                = list(map(lambda sample: std(log2_rpkm.loc[sample]), log2_rpkm.index))
    Sample_count_wo_InfNaN  = list(map(lambda sample: sum(log2_rpkm_finite.loc[sample]), log2_rpkm.index))
    Exon_count_wo_InfNaN    = list(map(lambda sample: sum(log2_rpkm_finite[sample]), log2_rpkm.columns))

    log2_QC = log2_rpkm.join(pd.DataFrame({'SampleQC': SampleQC, 'Sample_count_wo_InfNaN' :Sample_count_wo_InfNaN},
                            index = log2_rpkm.index))

    ExonQC_wo_outliers = list(ExonQC_wo_outliers[ExonQC_wo_outliers.columns[0]].values)
    log2_QC = pd.concat([log2_QC, pd.DataFrame({'ExonQC': ExonQC_wo_outliers + [np.NaN, np.NaN],
                                               'Exon_count_wo_InfNaN': Exon_count_wo_InfNaN + [np.NaN, np.NaN]},
                                               index = log2_QC.columns).T])

    rrr = pd.DataFrame({'rpkm_mean':list(map(lambda sample: np.mean(rpkm.loc[sample]), rpkm.index)), 
                        'rpkm_stddev': list(map(lambda sample: std(rpkm.loc[sample]), rpkm.index))},
                        index = rpkm.index)
    rrr = round(rrr, 2)
    sss = pd.DataFrame({'rpkm_mean': list(map(lambda sample: np.mean(rpkm.loc[sample]), rpkm.index)), 
                        'rpkm_stddev': list(map(lambda sample: std(rpkm.loc[sample]), rpkm.index)),
                        'SampleQC': SampleQC},
                        index = rpkm.index)
    sss = round(sss, 2)

    sample_mean_stddev_outfile = rpkm_matrix_file.replace('rpkm.txt', 'Sample_RPKM-means-stddevs_log2-stddevs')
    sss.to_csv(sample_mean_stddev_outfile, sep = '\t', index = False)

    cscore_outfile = rpkm_matrix_file.replace('rpkm.txt', 'Cscore_outfile')
    pval_outfile   = rpkm_matrix_file.replace('rpkm.txt', 'Pval_matrix') 

    c_scores = pd.DataFrame(map(lambda targ_coor, std_value: map(lambda x,y : x/y if y !=0 else np.NaN, log2_rpkm[targ_coor], [std_value] * sample_length),
                                log2_rpkm.columns, ExonQC_wo_outliers), 
                            index = log2_rpkm.columns, columns = log2_rpkm.index).T
    c_scores = round(c_scores, 2)
    pvals = pd.DataFrame(map(lambda targ_coor: map(lambda x: 1 - scipy.stats.norm.cdf(abs(x)), c_scores[targ_coor]), 
                                c_scores.columns), index = c_scores.columns, columns = c_scores.index)
    pvals = round(pvals, 2)
    c_scores.to_csv(cscore_outfile, sep = '\t', index = True)
    pvals.to_csv(pval_outfile, sep = '\t', index = True)

    fff = log2_QC.copy()
    #### work as R do ###################
    fff_idx = log2_QC.loc['ExonQC'] > threshold_exonQC
    failed_exons_by_ExonQC = log2_QC.loc['ExonQC'][fff_idx]

    count, Exon_count_wo_InfNaN_values = 0, []
    for i in fff_idx.values:
        if i:
            Exon_count_wo_InfNaN_values.append(Exon_count_wo_InfNaN[count])
        count += 1
    failed_exons_by_ExonQC = pd.DataFrame({'Exon': list(log2_QC.columns[fff_idx]), 
                                           'ExonQC': list(log2_QC.loc['ExonQC'][fff_idx]),
                                           'Exon_count_wo_InfNaN': Exon_count_wo_InfNaN_values})

    failed_samples = list(fff.index[fff['SampleQC'] > 0.2])
    failed_samples_by_SampleQC = pd.DataFrame({'failed_samples': failed_samples,
                                                'sd_SampleQC': fff['SampleQC'][fff['SampleQC'] > threshold_sampleQC],
                                                'Sample_count_wo_InfNaN': fff['Sample_count_wo_InfNaN'][fff['SampleQC'] > threshold_sampleQC],
                                                'rpkm_mean': rrr['rpkm_mean'][rrr.index.isin(failed_samples)]})
    ############# index self filling ####################

    sample_rpkm_stats = pd.DataFrame({'sample_stats':  [min(rrr['rpkm_mean']), max(rrr['rpkm_mean']), np.nanmedian(rrr['rpkm_mean']), 
                                                        min(rrr['rpkm_stddev']), max(rrr['rpkm_stddev']), np.nanmedian(rrr['rpkm_stddev'])]},
                                                        index = ['min_rpkm', 'max_rpkm', 'median_rpkm', 'min_stddev', 'max_stddev', 'median_stddev'])
    aaa = rpkm.T
    stat = scipy.stats.f_oneway(*list(map(lambda sample: aaa[sample], aaa.columns)))
    anova = pd.DataFrame({'fstatistic': [stat.statistic], 'pvalue': [stat.pvalue]})
    bbb =  aaa.melt(var_name='sample')
    results = ols('value ~ sample', data=bbb).fit().summary()

    midpool_summary_results = rpkm_matrix_file.replace('rpkm.txt', 'atlas_cnv_summary')
    with open(midpool_summary_results, 'w') as out:
        print(*['threshold_exonQC', str(threshold_exonQC)], sep = '\t', file = out)

    failed_exons_by_ExonQC.to_csv(midpool_summary_results, mode = 'a', sep ='\t', index =False)
    failed_samples_by_SampleQC.to_csv(midpool_summary_results, mode = 'a', sep ='\t', index =False)
    sample_rpkm_stats.to_csv(midpool_summary_results, mode = 'a', sep ='\t', index =False)
    anova.to_csv(midpool_summary_results, mode = 'a', sep ='\t', index =False)

    try:
        coefficient = results.tables[1].data
        failed_samples_by_anova_pval_lt_5pct = pd.DataFrame({'Prob_gt_t': list(map(lambda x: float(x[4]) if float(x[4]) < threshold_sample_anova else np.NaN,coefficient[1:]))}, 
                                                            index = list(map(lambda x: x[0] if float(x[4]) < threshold_sample_anova else np.NaN,coefficient[1:])))
        failed_samples_by_anova_pval_lt_5pct.to_csv(midpool_summary_results, mode = 'a', sep ='\t', index =True)
    except:
        pass

    list(map(lambda sample_id: call_cnvs(sample_id, os.path.dirname(outfile), failed_samples, failed_samples_by_anova_pval_lt_5pct, fff, rpkm, 
                                         threshold_del_soft, threshold_dup_soft, threshold_del, threshold_dup, median_wo_rpkm, failed_exons_by_ExonQC),
            rpkm.index))


def outliner_remove(df):
    rpkm_wo_outlier = df.copy()
    for sample in df.columns:
        top  = np.mean(df[sample]) + 1.96 * np.std(df[sample]) 
        down = np.mean(df[sample]) - 1.96 * np.std(df[sample]) 
        for i in range(len(df)):
            if df[sample][i] < down or df[sample][i] > top:
                rpkm_wo_outlier[sample][i] = np.NaN
    return rpkm_wo_outlier


def median_df(df):
    median, sample_id, sample_idx = [], [], []
    for name, exon in df.iteritems():
        rpkm_value = exon.values
        try:
            median_value = max(rpkm_value[rpkm_value <= np.nanmedian(rpkm_value)])
        except:
            median.append(np.NaN)
            sample_idx.append(np.NaN)
            sample_id.append(np.NaN)
            continue
        idx = int(np.where(rpkm_value == median_value)[0][0])
        median.append(median_value)
        sample_idx.append(idx+1)
        sample_id.append(exon.keys()[idx]) 
    median_rpkm = pd.DataFrame({'median_rpkm': median,'sample_id': sample_id, 'sample_idx': sample_idx},
                                index = df.columns)
    return median_rpkm


def std(values):
    ### freedom in R is n -1 while in numpy is n:
    return np.sqrt(np.power(np.nanstd(values), 2) * len(values) / (len(values) - 1))    

def call_cnvs(sample_id, batch_out, failed_samples, failed_samples_by_anova_pval_lt_5pct, fff, rpkm,
            threshold_del_soft, threshold_dup_soft, threshold_del, threshold_dup, median_wo_rpkm, failed_exons_by_ExonQC):
    anova_sample_id = 'sample[T.' + sample_id + ']'

    outdir = os.path.join(batch_out, sample_id)
    if os.system('mkdir -p ' + outdir) != 0:
        raise Exception('Failed to create dir: %s\n' %(outdir))
    
    if sample_id in failed_samples and anova_sample_id not in failed_samples_by_anova_pval_lt_5pct.index:
        outfile = os.path.join(outdir, sample_id + '.cnv.FAILED_sampleQC')
    elif sample_id not in failed_samples and anova_sample_id in failed_samples_by_anova_pval_lt_5pct.index:
        outfile = os.path.join(outdir, sample_id + '.cnv.FAILED_sampleANOVA')
    elif sample_id in failed_samples and anova_sample_id in failed_samples_by_anova_pval_lt_5pct.index:
        outfile = os.path.join(outdir, sample_id + '.cnv.FAILED_sampleQC_and_sampleANOVA')
    else:
        outfile = os.path.join(outdir, sample_id + '.cnv')
    header = ('Gene_Exon', 'cnv', 'log2R', 'rpkm', 'median_rpkm', 'Exon_Status', 'E_StDev', 'c_Score')

    with open(outfile, 'w') as out:
        print(*header, file = out, sep = '\t')

    sys.stderr.write(f"Calling  CNV on: {sample_id}\n")
    for cnv_type in ('del', 'dup'):
        sample_cnv = None
        if cnv_type is 'del':
            cnv_index = list(map(lambda x,y: True if x <= y else False, fff.loc[sample_id][:-2], threshold_del_soft['del'].values))
            cnv_index = list(map(lambda x,y: all((x,y)), fff.loc[sample_id][:-2] <= threshold_del, cnv_index))
        else:
            cnv_index = list(map(lambda x,y: True if x >= y else False, fff.loc[sample_id][:-2], threshold_dup_soft['dup'].values))
            cnv_index = list(map(lambda x,y: all((x,y)), fff.loc[sample_id][:-2] >= threshold_del, cnv_index))
        if sum(cnv_index) == 0:
            sys.stderr.write(f"                  {sample_id} has no cnv dels.\n")
        else:
            sys.stderr.write(f"                  {sample_id} has {sum(cnv_index)} cnv dels.\n")
        cnv_index.append(False)
        cnv_index.append(False)
        #### extra column added: 'SampleQC', 'Sample_count_wo_InfNaN'

        sample_cnv = pd.DataFrame({ 'Gene_Exon': fff.columns[cnv_index],
                                    'cnv'      : ['del'] * (sum(cnv_index)),
                                    'log2R'    : fff.loc[sample_id][cnv_index].values,
                                    'rpkm'     : rpkm.loc[sample_id][cnv_index[:-2]].values,
                                    'median_rpkm' : median_wo_rpkm['median_rpkm'][median_wo_rpkm.index.isin(fff.columns[cnv_index])],
                                    'Exon_Status' : list(map(lambda x: 'Fail' if x else 'Pass', fff.columns[cnv_index].isin(failed_exons_by_ExonQC['Exon']))),
                                    'E_StDev'     : fff.loc['ExonQC'][cnv_index].values,
                                    'c_Score'     : list((fff.loc[sample_id][cnv_index].values / fff.loc['ExonQC'][cnv_index].values).round(2))})

        sample_cnv.to_csv(outfile, sep = '\t', index = False, mode = 'a')
