import re
import os
import sys

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


def burden(case_matrix, control_matrix, gene_col_name='Gene', cutoff=0):
    N1 = case_matrix.shape[1] - 1
    N2 = control_matrix.shape[1] - 1
    df_matrix = pd.merge(case_matrix, control_matrix, how='left', on=[gene_col_name])
    df_matrix.fillna(0, inplace=True)
    gene_list, pvalue_list, score_list = [], [], []
    for tup in df_matrix.itertuples():
        gene_list.append(tup[1])
        s1 = list(tup[2:N1 + 2])
        s1 = [eval(x) for x in ','.join([str(_) for _ in s1]).split(',')]
        s1 = [s1[i] for i in range(0, len(s1)) if s1[i] > cutoff]
        if not s1:
            s1 = [0]
        s2 = list(tup[N1 + 2:])
        s2 = [eval(x) for x in ','.join([str(_) for _ in s2]).split(',')]
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


def uniq_gene_df(df_gene_col, df_score_col, mode):
    genes = list(df_gene_col)
    scores = list(df_score_col)
    gene_score_dict = {}
    for i, j in zip(genes, scores):
        gene_score_dict[i] = gene_score_dict.get(i, '') + str(j) + ','
    genes, scores = [], []
    if mode == 'AD':
        for k, v in gene_score_dict.items():
            genes.append(k)
            _score = v.rstrip(',').split(',')
            scores.append(max([eval(_) for _ in _score]))
    else:
        for k, v in gene_score_dict.items():
            _score = v.rstrip(',').split(',')
            if len(_score) == 2:
                genes.append(k)
                scores.append(v.rstrip(','))
            elif len(_score) > 2:
                genes.append(k)
                _score.sort()
                m1 = _score[-1]
                m2 = _score[-2]
                scores.append(','.join([m1, m2]))
    return genes, scores


def get_matrix(file_dir, tag, config, out_dir='', mode='AD',
               gene_col_name='Gene', score_col_name='Score', snvdb='', pattern=True):
    if pattern:
        for num, file in enumerate(os.listdir(file_dir)):
            if not file.startswith('.'):
                filter_score(file_dir + '/' + file, out_dir, config, mode, snvdb)
            else:
                continue
            df = pd.read_table(out_dir + '/' + file.split('/')[-1] + '.score', low_memory=False)
            genes, scores = uniq_gene_df(df[gene_col_name], df[score_col_name], mode)
            if num == 0:
                df_matrix = pd.DataFrame(columns=[gene_col_name, tag + str(num + 1)])
                df_matrix[gene_col_name] = genes
                df_matrix[tag + str(num + 1)] = scores
            elif not file.startswith('.'):
                _df = pd.DataFrame(columns=[gene_col_name, tag + str(num + 1)])
                _df[gene_col_name] = genes
                _df[tag + str(num + 1)] = scores
                df_matrix = pd.merge(df_matrix, _df, how='outer', on=[gene_col_name])
    else:
        for num, file in enumerate(os.listdir(file_dir)):
            if not file.startswith('.'):
                if snvdb:
                    df = rm_snvdb(file_dir + '/' + file, snvdb)
                else:
                    df = pd.read_table(file_dir + '/' + file, low_memory=False)
            else:
                continue
            genes, scores = uniq_gene_df(df[gene_col_name], df[score_col_name], mode)
            if num == 0:
                df_matrix = pd.DataFrame(columns=[gene_col_name, tag + str(num + 1)])
                df_matrix[gene_col_name] = genes
                df_matrix[tag + str(num + 1)] = scores
            elif not file.startswith('.'):
                _df = pd.DataFrame(columns=[gene_col_name, tag + str(num + 1)])
                _df[gene_col_name] = genes
                _df[tag + str(num + 1)] = scores
                df_matrix = pd.merge(df_matrix, _df, how='outer', on=[gene_col_name])
    return df_matrix


def filter_score(file, outdir, config, mode, ctlf=''):
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
    # 读取文件
    df = pd.read_table(file)
    # AF过滤
    df.replace('.', 0, inplace=True)
    AF_list = [_.strip() for _ in config['cc_AF'].split(',')]
    try:
        df[AF_list] = df[AF_list].apply(pd.to_numeric, errors='ignore')
    except:
        sys.exit('[ E: something error when read af information! ]')
    df['Chr'] = df['Chr'].astype(str)
    df['Chr'] = [_.lstrip('chr') for _ in list(df['Chr'])]
    df.drop_duplicates(inplace=True)
    if mode == 'AD':
        af_th = [eval(_.strip()) for _ in config['cc_AF_AD'].split(',')]
    else:
        af_th = [eval(_.strip()) for _ in config['cc_AF_AR'].split(',')]
    for af, th in zip(AF_list, af_th):
        df = df[df[af] < th]
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
    # 打分
    col_list = df.columns.to_list()
    _revel_idx = col_list.index('REVEL') + 1
    _efc_index = col_list.index('ExonicFunc.refGene') + 1
    _fc_index = col_list.index('Func.refGene') + 1
    _gene_index = col_list.index('Gene.refGene') + 1
    _format_index = col_list.index('FORMAT') + 1
    _splice_index = []
    splice_list = [_.strip() for _ in config['cc_splice'].split(',')]
    for db in splice_list:
        _splice_index.append(col_list.index(db) + 1)
    outfile = open(outdir + '/' + file.split('/')[-1] + '.score', 'w')
    outfile.write('Gene\tScore\n')
    for tup in df.itertuples():
        dp_index = tup[_format_index].split(':').index('DP')
        ad_index = tup[_format_index].split(':').index('AD')
        dp = int(tup[_format_index + 1].split(':')[dp_index])
        if dp < 10 or int(tup[_format_index + 1].split(':')[ad_index].split(',')[1]) / dp <= 0.2:
            continue
        score = 0
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
        for gene in genes.split(';'):
            if gt in ['1/1']:
                outfile.write(gene + '\t' + str(score / 100) + '\n')
                outfile.write(gene + '\t' + str(score / 100) + '\n')
            elif gt in ['0/1']:
                outfile.write(gene + '\t' + str(score / 100) + '\n')
    outfile.close()


def build_snvdb(file_path, out_dir, out_file_name, script_path, file_type='vcf', rate=0.5):
    total_var_dict = {}
    # 获得输入目录中的vcf文件
    file_names = []
    if file_type in ['vcf', 'VCF', 'Vcf']:
        for file in os.listdir(file_path):
            if re.findall(r'(\.vcf\.gz|\.vcf)$', file):
                file_names.append(file)
    elif file_type in ['avinput']:
        for file in os.listdir(file_path):
            if not file.startswith('.'):
                file_names.append(file)
    if len(file_names) < 1:
        sys.exit('[ E: lack input vcf file ！]')
    # 获得记录及出现次数
    if file_type in ['vcf', 'VCF', 'Vcf']:
        num = 0
        for file in file_names:
            file_name = file_path + '/' + file
            tmp_out = out_dir + '/' + file
            tran_file_cmd = 'perl %s/bin/annovar/convert2annovar.pl -format vcf4 %s -outfile %s -allsample' % (
                script_path, file_name, tmp_out)
            result = os.system(tran_file_cmd)
            if result:
                print(
                    '[ E: Something wrong with change format with input vcf file < %s > ！remove from program !]' % file)
            else:
                tmp_file_names = []
                for _file in os.listdir(out_dir):
                    if re.findall(r'(\.avinput)$', _file):
                        tmp_file_names.append(_file)
                num += len(tmp_file_names) * 2
                for _ in tmp_file_names:
                    with open(out_dir + '/' + _) as f:
                        for line in f:
                            records = line.strip().split('\t')
                            key = '\t'.join(records[:5]).lstrip('chr')
                            if records[5].startswith('het'):
                                total_var_dict[key] = total_var_dict.get(key, 0) + 1
                            elif records[5].startswith('hom'):
                                total_var_dict[key] = total_var_dict.get(key, 0) + 2
            clear_tmp_cmd = 'rm -rf %s/* ' % out_dir
            os.system(clear_tmp_cmd)
        out_file = open('./' + out_file_name, 'a+')
        for k, v in total_var_dict.items():
            if v / num >= rate:
                out_file.write(k + '\n')
        out_file.close()
        clear_tmp_cmd = 'rm -rf %s ' % out_dir
        result = os.system(clear_tmp_cmd)
        if result:
            sys.exit('[ E: Something wrong with clear temp file ！]')
        else:
            print('[ S: False positive database build successfully ！]')
    elif file_type in ['avinput']:
        num = len(file_names)
        for file in file_names:
            file_name = file_path + '/' + file
            with open(file_name) as f:
                for line in f:
                    key = '\t'.join(line.strip().split('\t')[:5]).lstrip('chr')
                    total_var_dict[key] = total_var_dict.get(key, 0) + 1
        out_file = open('./' + out_file_name, 'a+')
        for k, v in total_var_dict.items():
            if v / num >= rate:
                out_file.write(k + '\n')
        out_file.close()
        print('[ S: False positive database build successfully ！]')


def rm_snvdb(file, snvdb):
    df = pd.read_table(file, low_memory=False)
    df_col = df.columns.to_list()
    new_col = ['Chr', 'Start', 'End', 'Ref', 'Alt'] + df_col[5:]
    df.columns = new_col
    df['Chr'] = df['Chr'].astype(str)
    df['Chr'] = [_.lstrip('chr') for _ in list(df['Chr'])]
    try:
        dfctl = pd.read_table(snvdb, low_memory=False, header=None)
        dfctl.rename(columns={0: 'Chr', 1: 'Start', 2: 'End', 3: 'Ref', 4: 'Alt'}, inplace=True)
        dfctl['Chr'] = dfctl['Chr'].astype(str)
        dfctl['Chr'] = [_.lstrip('chr') for _ in list(dfctl['Chr'])]
    except:
        print('[ W: 无法读取 control 文件，请检查 ！]')
        dfctl = pd.DataFrame(columns=['Chr', 'Start', 'End', 'Ref', 'Alt'])
    dfcon = pd.merge(df, dfctl, how='inner', on=['Chr', 'Start', 'End', 'Ref', 'Alt'])
    dfcon = dfcon[['Chr', 'Start', 'End', 'Ref', 'Alt']]
    dfcon['Chr'] = dfcon['Chr'].astype(str)
    df = pd.concat([df, dfcon])
    df.drop_duplicates(subset=['Chr', 'Start', 'End', 'Ref', 'Alt'], keep=False, inplace=True)
    return df
