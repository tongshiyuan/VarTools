import os
import numpy as np
import pandas as pd
import scipy.stats as st


def Z_binomial(n1, n2, N1, N2):
    r = N2 / N1
    n = n1 / (n1 + n2) - 1 / (1 + r)
    d = (r / ((1 + r) * (1 + r) * (n1 + n2))) ** 0.5
    return n / d


def gene_rank(s1, s2, N1, N2):
    """
    s1 = case score group
    s2 = control score group
    N1 = total number of case
    N2 = total number pf control
    """
    n1, n2 = len(s1), len(s2)
    if s2 == [0]:
        n2 = 0
    Z1 = Z_binomial(n1, n2, N1, N2)
    p1 = st.norm.logsf(Z1)
    _, p2 = st.mannwhitneyu(s1, s2,
                            method='asymptotic',
                            alternative='greater')
    p = st.chi2.sf(-2 * (np.log(p2) + p1), 4)
    s = st.chi2.isf(p, 4)
    return p, s


def cc_preprocess(case_dir, control_dir, out_dir, snvdb, mode):
    out_case = out_dir + '/case_out'
    out_control = out_dir + '/control_out'
    os.makedirs(out_case)
    os.makedirs(out_control)
    tmp_rank = open(out_dir + '/rank_raw.txt', 'w')
    for num, file in enumerate(os.listdir(case_dir)):
        if num == 0:
            with open(case_dir + '/' + file) as f:
                header = f.readline()
                tmp_rank.write('SampleName\t' + header)
        if not file.startswith('.'):
            filter_score(case_dir + '/' + file, out_case, mode, snvdb, tmp_rank)
    for file in os.listdir(control_dir):
        if not file.startswith('.'):
            filter_score(control_dir + '/' + file, out_control, mode, snvdb)
    tmp_rank.close()

    return out_case, out_control


def div_cc(in_case_dir, in_control_dir, gene_col_name='GENE', score_col_name='SCORE', cutoff=0):
    for num, file in enumerate(os.listdir(in_case_dir)):
        if num == 0:
            df_case = pd.read_table(in_case_dir + '/' + file, low_memory=False)
            df_case = df_case[[gene_col_name, score_col_name]].copy()
            df_case.rename(columns={score_col_name: 'CASE_' + str(num + 1)}, inplace=True)
        elif not file.startswith('.'):
            _df = pd.read_table(in_case_dir + '/' + file, low_memory=False)
            _df = _df[[gene_col_name, score_col_name]].copy()
            _df.rename(columns={score_col_name: 'CASE_' + str(num + 1)}, inplace=True)
            df_case = pd.merge(df_case, _df, how='outer', on=[gene_col_name])
    N1 = df_case.shape[1] - 1

    for num, file in enumerate(os.listdir(in_control_dir)):
        if not file.startswith('.'):
            _df = pd.read_table(in_control_dir + '/' + file, low_memory=False)
            _df = _df[[gene_col_name, score_col_name]].copy()
            _df.rename(columns={score_col_name: 'CONTROL_' + str(num + 1)}, inplace=True)
            df_matrix = pd.merge(df_matrix, _df, how='left', on=[gene_col_name])
    N2 = df_matrix.shape[1] - 1 - N1
    df_matrix.fillna(0, inplace=True)

    gene_list, pvalue_list, score_list = [], [], []
    for tup in df_matrix.itertuples():
        gene_list.append(tup[1])
        s1 = list(tup[2:N1 + 2])
        s1 = [s1[i] for i in range(0, len(s1)) if s1[i] > cutoff]
        if not s1:
            s1 = [0]
        s2 = list(tup[N1 + 2:])
        s2 = [s2[i] for i in range(0, len(s2)) if s2[i] > cutoff]
        if not s2:
            s2 = [0]
        p, s = gene_rank(s1, s2, N1, N2)
        pvalue_list.append(p)
        score_list.append(s)
    df = pd.DataFrame(columns=['#Gene', 'Twopart_P_value', 'Score'])
    df['#Gene'] = gene_list
    df['Twopart_P_value'] = pvalue_list
    df['Score'] = score_list
    df.sort_values(by='Score', ascending=False, inplace=True)
    return df


# 过滤和打分模块
def filter_score(file, outdir, mode, ctlf='', case=False):
    score_conf = {
        'splice': 80,
        'nsSNV': 40,
        'predscore': 60,
        'stopgain': 90,
        'stoploss': 90,
        'frameshift deletion': 90,
        'frameshift insertion': 90,
        'frameshift block substitution': 90,
        'nonframeshift insertion': 40,
        'nonframeshift deletion': 40,
        'nonframshift block substitution': 40,
        'unknown': 20,
        'UTR_InterGenic': 0
    }
    gene_dict = {}
    df = pd.read_table(file, low_memory=False)
    df.replace('.', 0, inplace=True)
    df[['1000g2015aug_all', 'ExAC_EAS', 'gnomAD_exome_EAS']] = df[
        ['1000g2015aug_all', 'ExAC_EAS', 'gnomAD_exome_EAS']].apply(pd.to_numeric, errors='ignore')
    df['Chr'] = df['Chr'].astype(str)
    _tmp_chr = df['Chr']
    df['Chr'] = [_.lstrip('chr') for _ in _tmp_chr]
    df.drop_duplicates(inplace=True)
    if mode == 'AD':
        df = df[(df['1000g2015aug_all'] < 0.0001) & (df['ExAC_EAS'] < 0.0001) & (df['gnomAD_exome_EAS'] < 0.0001)]
    elif mode == 'AR':
        df = df[(df['1000g2015aug_all'] < 0.005) & (df['ExAC_EAS'] < 0.005) & (df['gnomAD_exome_EAS'] < 0.005)]
    # 对照过滤
    if ctlf:
        try:
            dfctl = pd.read_table(ctlf, low_memory=False, header=None)
            dfctl.rename(columns={0: 'Chr', 1: 'Start', 2: 'End', 3: 'Ref', 4: 'Alt'}, inplace=True)
        except:
            print('[ W: 无法读取 control 文件，请检查 ！]')
            dfctl = pd.DataFrame(columns=['Chr', 'Start', 'End', 'Ref', 'Alt'])
        dfcon = pd.merge(df, dfctl, how='inner', on=['Chr', 'Start', 'End', 'Ref', 'Alt'])
        dfcon = dfcon[['Chr', 'Start', 'End', 'Ref', 'Alt']]
        dfcon['Chr'] = dfcon['Chr'].astype(str)
        df = pd.concat([df, dfcon])
        df.drop_duplicates(subset=['Chr', 'Start', 'End', 'Ref', 'Alt'], keep=False, inplace=True)
    col_list = df.columns.to_list()
    _revel_idx = col_list.index('REVEL') + 1
    _efc_index = col_list.index('ExonicFunc.refGene') + 1
    _fc_index = col_list.index('Func.refGene') + 1
    _gene_index = col_list.index('Gene.refGene') + 1
    _format_index = col_list.index('FORMAT') + 1
    _splice_index = []
    for db in ['scsnv', 'SpliceAI']:
        _splice_index.append(col_list.index(db) + 1)
    outfile = open(outdir + '/' + file.split('/')[-1] + '.score', 'w')
    outfile.write('#CHROM\tPOS\tREF\tALT\tGENE\tSCORE\n')
    for tup in df.itertuples():
        dp_index = tup[_format_index].split(':').index('DP')
        ad_index = tup[_format_index].split(':').index('AD')
        dp = int(tup[_format_index + 1].split(':')[dp_index])
        if dp < 10 or int(tup[_format_index + 1].split(':')[ad_index].split(',')[1]) / dp <= 0.2:
            continue
        score = 0
        varlist = tup[1:3] + tup[4:6]
        var = '\t'.join([str(i) for i in varlist])
        if tup[_efc_index] == 'synonymous SNV':  # synonymous or splicing
            for index in _splice_index:
                if tup[index] not in ['.', '0', 0]:
                    score = score_conf['splice']  # splice score = 80, synonymous SNV = 0
                    break
        elif tup[_fc_index].find('UTR') != -1 or tup[_fc_index].find('intergenic') != -1:
            score = score_conf['UTR_InterGenic']
        elif tup[_efc_index] == 'nonsynonymous SNV':  # nsSNV = 40, predscore = 60
            score = score_conf['nsSNV'] + score_conf['predscore'] * float(tup[_revel_idx])
        elif tup[_efc_index] == 'stopgain':
            score = score_conf['stopgain']  # stopgain = 90
        elif tup[_efc_index] == 'stoploss':
            score = score_conf['stoploss']  # stoploss = 90
        elif tup[_efc_index] == 'frameshift deletion':
            score = score_conf['frameshift deletion']  # frameshift deletion = 90
        elif tup[_efc_index] == 'frameshift insertion':
            score = score_conf['frameshift insertion']  # frameshift insertion = 90
        elif tup[_efc_index] == 'frameshift block substitution':
            score = score_conf['frameshift block substitution']  # frameshift block substitution = 90
        elif tup[_efc_index] == 'nonframeshift insertion':
            score = score_conf['nonframeshift insertion']  # nonframeshift insertion = 40
        elif tup[_efc_index] == 'nonframeshift deletion':
            score = score_conf['nonframeshift deletion']  # nonframeshift deletion = 40
        elif tup[_efc_index] == 'nonframshift block substitution':
            score = score_conf['nonframshift block substitution']  # nonframshift block substitution = 60
        elif tup[_efc_index] == 'unknown':
            for index in _splice_index:
                if tup[index] not in ['.', '0', 0]:
                    score = score_conf['splice']
                    break
            else:
                score = score_conf['unknown']  # unknown = 20
        elif tup[_efc_index] not in ['.', '0', 0] and tup[_efc_index].find('splicing') != -1:
            score = score_conf['splice']
        else:
            for index in _splice_index:
                if tup[index] not in ['.', '0', 0]:
                    score = score_conf['splice']
        genes = tup[_gene_index]
        gt = tup[_format_index + 1].split(':')[0]
        if score == 0:
            continue
        _tmp_tup = list(tup)
        for gene in genes.split(';'):
            if case:
                _tmp_tup[_gene_index] = gene
                case.write('\t'.join([col_list[_format_index]] + [str(_) for _ in _tmp_tup[1:]]) + '\n')
            #     if gene not in gene_dict.keys() or (gene in gene_dict.keys() and gene_dict[gene]['score'] <= score):
            #         gene_dict[gene] = {'var': var,
            #                            'score': score}
            if gt in ['1/1']:
                if gene not in gene_dict.keys():
                    gene_dict[gene] = {'var': [var, var],
                                       'score': [score, score]}
                else:
                    if score >= max(gene_dict[gene]['score']):
                        gene_dict[gene] = {'var': [var, var],
                                           'score': [score, score]}
                    elif score > min(gene_dict[gene]['score']) and gene_dict[gene]['score'][0] >= \
                            gene_dict[gene]['score'][1]:
                        gene_dict[gene]['var'][1] = var
                        gene_dict[gene]['score'][1] = score
                    elif score > min(gene_dict[gene]['score']) and gene_dict[gene]['score'][0] < \
                            gene_dict[gene]['score'][1]:
                        gene_dict[gene]['var'][0] = var
                        gene_dict[gene]['score'][0] = score
            elif gt in ['0/1']:
                if gene not in gene_dict.keys():
                    gene_dict[gene] = {'var': [var],
                                       'score': [score]}
                else:
                    if len(gene_dict[gene]['var']) == 1:
                        gene_dict[gene]['var'] += [var]
                        gene_dict[gene]['score'] += [score]
                    else:
                        if score > min(gene_dict[gene]['score']) and gene_dict[gene]['score'][0] >= \
                                gene_dict[gene]['score'][1]:
                            gene_dict[gene]['var'][1] = var
                            gene_dict[gene]['score'][1] = score
                        elif score > min(gene_dict[gene]['score']) and gene_dict[gene]['score'][0] < \
                                gene_dict[gene]['score'][1]:
                            gene_dict[gene]['var'][0] = var
                            gene_dict[gene]['score'][0] = score
            # if gt in ['1/1']:
            #     outfile.write(
            #         var + '\t' + gene + '\t' + str(score) + '\n' + var + '\t' + gene + '\t' + str(score) + '\n')
            # elif gt in ['0/1']:
            #     outfile.write(var + '\t' + gene + '\t' + str(score) + '\n')
    if mode == 'AD':
        for k, v in gene_dict.items():
            if (len(v['score']) == 1 and v['score'][0] == 0) or (len(v['score']) == 2 and max(v['score']) == 0):
                pass
            elif len(v['score']) == 1 or v['score'][0] >= v['score'][1]:
                outfile.write(v['var'][0] + '\t' + k + '\t' + str(v['score'][0] / 100) + '\n')
            else:
                outfile.write(v['var'][1] + '\t' + k + '\t' + str(v['score'][1] / 100) + '\n')
    else:
        for k, v in gene_dict.items():
            if (len(v['score']) == 1 and v['score'][0] == 0) or (len(v['score']) == 2 and max(v['score']) == 0):
                pass
            elif len(v['score']) == 1:
                outfile.write(v['var'][0] + '\t' + k + '\t' + str(v['score'][0] / 100) + '\n')
            elif v['score'][0] >= v['score'][1]:
                outfile.write(v['var'][0] + '\t' + k + '\t' + str(v['score'][0] / 100) + '\n')
                outfile.write(v['var'][1] + '\t' + k + '\t' + str(v['score'][1] / 100) + '\n')
            else:
                outfile.write(v['var'][1] + '\t' + k + '\t' + str(v['score'][1] / 100) + '\n')
                outfile.write(v['var'][0] + '\t' + k + '\t' + str(v['score'][0] / 100) + '\n')
    outfile.close()


def variants_with_rank(df_gene_rank, raw_variants, variants_out):
    df_gene_rank['rank'] = df_gene_rank.index.to_list()
    df_gene_rank['rank'] = df_gene_rank['rank'] + 1
    df_gene_rank.rename(columns={'#Gene': 'Gene.refGene'}, inplace=True)
    del df_gene_rank['Twopart_P_value'], df_gene_rank['Score']
    # 输出带rank的结果
    df_var_raw = pd.read_table(raw_variants, low_memory=False)
    df_merge = pd.merge(df_gene_rank, df_var_raw, how='left', on=['Gene.refGene'])
    # df_merge = df_merge[df_merge['rank'] <= 10]
    df_merge.to_csv(variants_out, sep='\t', index=False)
