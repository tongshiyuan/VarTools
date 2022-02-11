def bioinfo_scores(file, sample, gene_range,
                   freq_index, var_type_index, var_func_index, var_gene_index, pred_index, splicing_index,
                   report_clin_index, report_xy_index,
                   type_score_matrix, omim_index, format_index,
                   pred_score=60,
                   clin_score='-90,-80,80,90', xy_score='80,100', report_rate='1:1',
                   freq_score='90,80,60,40', freq_threshold='0.01', freq_site='0.0001,0.001,0.005',
                   freq_null_score=100, bio_rate='1:1:0.2', dp=10, vaf=0.2):
    with open(file) as f:
        header = f.readline()
        db_split = header.rstrip().split('\t')
        try:
            sample_index = db_split.index(sample)
        except:
            print('[ E: Can not find sample < %s > ！]' % sample)
            return 0
        freq_site_list = freq_site.split(',')
        freq_score_list = freq_score.split(',')
        scores_dict = {}
        # freq_score_dict = {}
        # pred_score_dict = {}
        # report_score_dict = {}
        _t = sum([eval(i) for i in bio_rate.split(':')])
        _r1 = eval(bio_rate.split(':')[0]) / _t
        _r2 = eval(bio_rate.split(':')[1]) / _t
        _r3 = eval(bio_rate.split(':')[2]) / _t
        for i in f:
            colunoms = i.rstrip().split('\t')
            try:
                if colunoms[sample_index].split(':')[0] in ['0/0', './.']:
                    pass
                else:
                    DP_index = colunoms[format_index].split(':').index('DP')
                    AD_index = colunoms[format_index].split(':').index('AD')
                    tmp_dp = eval(colunoms[sample_index].split(':')[DP_index])
                    tmp_ap = eval(colunoms[sample_index].split(':')[AD_index].split(',')[-1])
            except:
                print('[ W: Can not find qc info in this variation <%s> ! ]' % i)
                DP_index = AD_index = 0
            if colunoms[sample_index].split(':')[0] in ['0/0', './.']:
                pass
            elif DP_index != AD_index and (tmp_dp < dp or tmp_ap / tmp_dp < vaf):
                pass
            else:
                records = i.rstrip().split('\t')
                key = '\t'.join(records[:5])
                annotation = '\t'.join(records[5:])
                scores_dict[key] = {'pred_score': 0,
                                    'freq_score': 0,
                                    'report_score': 0,
                                    'bio_score': 0,
                                    'all_score': 0,
                                    'genetic_score': False,
                                    'anno': annotation,
                                    'gene': records[var_gene_index]
                                    }
                # 类型打分
                if records[sample_index].split(':')[0] in ['1/1']:
                    scores_dict[key]['genetic_score'] = True
                # 类型、预测打分
                if records[var_type_index] == 'synonymous SNV':
                    # scores_dict[key]['pred_score'] = 0
                    for index in splicing_index:
                        if records[index] in ['.', '0']:
                            continue
                            # scores_dict[key]['pred_score'] = 0
                            # pred_score_dict[key] = type_score_matrix[10]
                        else:
                            # _splicing_score = max([eval(i) for i in records[splicing_index].split(';') if i != '.'])
                            # scores_dict[key]['pred_score'] = splicing_score * _splicing_score
                            scores_dict[key]['pred_score'] = type_score_matrix[11]
                            break
                    else:
                        scores_dict[key]['pred_score'] = 0
                    # pred_score_dict[key] = type_score_matrix[0]
                # elif gene_range == 'omim' and not bool(
                #         [gene for gene in records[var_gene_index].split(',') if gene in gene_list]):
                #     scores_dict[key]['pred_score'] = 0
                elif gene_range == 'omim' and records[omim_index] in ['.', '0']:
                    scores_dict[key]['pred_score'] = 0
                elif records[var_func_index].find('UTR') != -1 or records[var_func_index].find('intergenic') != -1:
                    scores_dict[key]['pred_score'] = 0
                elif records[var_type_index] == 'nonsynonymous SNV':
                    if records[pred_index] == '.':
                        scores_dict[key]['pred_score'] = type_score_matrix[1]
                        # pred_score_dict[key] = type_score_matrix[1]
                    else:
                        scores_dict[key]['pred_score'] = type_score_matrix[1] + pred_score * eval(records[pred_index])
                        # pred_score_dict[key] = type_score_matrix[1] + eval(pred_score) * eval(records[pred_index])
                elif records[var_type_index] == 'stopgain':
                    scores_dict[key]['pred_score'] = type_score_matrix[2]
                    # pred_score_dict[key] = type_score_matrix[2]
                elif records[var_type_index] == 'frameshift deletion':
                    scores_dict[key]['pred_score'] = type_score_matrix[3]
                    # pred_score_dict[key] = type_score_matrix[3]
                elif records[var_type_index] == 'frameshift insertion':
                    scores_dict[key]['pred_score'] = type_score_matrix[4]
                    # pred_score_dict[key] = type_score_matrix[4]
                elif records[var_type_index] == 'frameshift block substitution':
                    scores_dict[key]['pred_score'] = type_score_matrix[5]
                    # pred_score_dict[key] = type_score_matrix[5]
                elif records[var_type_index] == 'stoploss':
                    scores_dict[key]['pred_score'] = type_score_matrix[6]
                    # pred_score_dict[key] = type_score_matrix[6]
                elif records[var_type_index] == 'nonframeshift insertion':
                    scores_dict[key]['pred_score'] = type_score_matrix[7]
                    # pred_score_dict[key] = type_score_matrix[7]
                elif records[var_type_index] == 'nonframeshift deletion':
                    scores_dict[key]['pred_score'] = type_score_matrix[8]
                    # pred_score_dict[key] = type_score_matrix[8]
                elif records[var_type_index] == 'nonframshift block substitution':
                    scores_dict[key]['pred_score'] = type_score_matrix[9]
                    # pred_score_dict[key] = type_score_matrix[9]
                elif records[var_type_index] == 'unknown':
                    # pred_score_dict[key] = type_score_matrix[10]
                    # if records[splicing_index] == '.':
                    #     scores_dict[key]['pred_score'] = type_score_matrix[10]
                    #     # pred_score_dict[key] = type_score_matrix[10]
                    # else:
                    #     _splicing_score = max([eval(i) for i in records[splicing_index].split(';') if i != '.'])
                    #     scores_dict[key]['pred_score'] = type_score_matrix[10] + splicing_score * _splicing_score
                    # pred_score_dict[key] = type_score_matrix[10] + eval(splicing_score) * eval(records[splicing_index])

                    for index in splicing_index:
                        if records[index] in ['.', '0']:
                            continue
                        else:
                            scores_dict[key]['pred_score'] = type_score_matrix[11]
                            break
                    else:
                        scores_dict[key]['pred_score'] = type_score_matrix[10]

                elif records[var_func_index].find('splicing') != -1:
                    scores_dict[key]['pred_score'] = type_score_matrix[11]
                    # if records[splicing_index] in ['.', '0']:
                    #     scores_dict[key]['pred_score'] = type_score_matrix[11]
                    #     # pred_score_dict[key] = type_score_matrix[11]
                    # else:
                    #     _splicing_score = max([eval(i) for i in records[splicing_index].split(';') if i != '.'])
                    #     scores_dict[key]['pred_score'] = type_score_matrix[11] + (
                    #             100 - type_score_matrix[11]) * _splicing_score
                    # pred_score_dict[key] = type_score_matrix[11] + (100 - type_score_matrix[11]) * eval(records[splicing_index])
                    # pred_score_dict[key] = type_score_matrix[11]
                else:
                    # if records[splicing_index] != '.':
                    #     _splicing_score = max([eval(i) for i in records[splicing_index].split(';') if i != '.'])
                    #     scores_dict[key]['pred_score'] = splicing_score * _splicing_score
                    # pred_score_dict[key] = eval(splicing_score) * eval(records[splicing_index])
                    for index in splicing_index:
                        if records[index] in ['.', '0']:
                            continue
                        else:
                            scores_dict[key]['pred_score'] = type_score_matrix[11]

                # 频率打分
                af_max = 0.00000000001
                for af_index in freq_index:
                    if records[af_index] != '.' and eval(records[af_index]) >= af_max:
                        af_max = eval(records[af_index])
                if af_max == 0.00000000001:
                    scores_dict[key]['freq_score'] = freq_null_score
                    # freq_score_dict[key] = eval(freq_null_score)
                elif af_max >= freq_threshold:
                    pass
                else:
                    for site_index in range(1, len(freq_site_list) + 1):
                        if af_max >= eval(freq_site_list[-site_index]):
                            scores_dict[key]['freq_score'] = eval(freq_score_list[-site_index])
                            # freq_score_dict[key] = freq_score_list[-site_index]
                            break
                    else:
                        scores_dict[key]['freq_score'] = eval(freq_score_list[0])
                        # freq_score_dict[key] = freq_score_list[0]
                # 报道打分
                _clin_score = [eval(num) for num in clin_score.split(',')]
                _xy_score = [eval(num) for num in xy_score.split(',')]
                if records[report_clin_index].find('Likely_benign') != -1:
                    clin_tmp = _clin_score[1]
                elif records[report_clin_index].find('Benign') != -1:
                    clin_tmp = _clin_score[0]
                elif records[report_clin_index].find('Pathogenic') != -1:
                    clin_tmp = _clin_score[3]
                elif records[report_clin_index].find('Likely_pathogenic') != -1:
                    clin_tmp = _clin_score[2]
                elif records[report_clin_index].find('Uncertain_significance') != -1:
                    clin_tmp = 0.0001
                else:
                    clin_tmp = 0
                if records[report_xy_index] not in ['.', '0']:
                    if records[report_xy_index].split(';')[1] == 'DM':
                        xy_tmp = _xy_score[1]
                    elif records[report_xy_index].split(';')[1] == 'DM?':
                        xy_tmp = _xy_score[0]
                else:
                    xy_tmp = 0
                scores_dict[key]['report_score'] = clin_tmp * (
                        eval(report_rate.split(':')[0]) / (eval(report_rate.split(':')[0]) + eval(
                    report_rate.split(':')[1]))) + xy_tmp * (eval(report_rate.split(':')[1]) / (eval(
                    report_rate.split(':')[0]) + eval(report_rate.split(':')[1])))
                # report_score_dict[key] = clin_tmp * (eval(report_rate.split(':')[0]) /
                # eval(report_rate.split(':')[0]) + eval(report_rate.split(':')[1])) +
                # xy_tmp * (eval(report_rate.split(':')[1]) / eval(report_rate.split(':')[0]) +
                # eval(report_rate.split(':')[1]))
                scores_dict[key]['bio_score'] = scores_dict[key]['freq_score'] * _r1 + scores_dict[key][
                    'pred_score'] * _r2 + scores_dict[key]['report_score'] * _r3

    del_var = []
    for k in scores_dict.keys():
        if scores_dict[k]['freq_score'] == 0 or (
                scores_dict[k]['pred_score'] == 0 and scores_dict[k]['report_score'] <= 0):
            del_var.append(k)
    for var in del_var:
        del scores_dict[var]

    # total_score_dict = {}
    # _t = sum([eval(i) for i in bio_rate.split(':')])
    # _r1 = eval(bio_rate.split(':')[0]) / _t
    # _r2 = eval(bio_rate.split(':')[1]) / _t
    # _r3 = eval(bio_rate.split(':')[2]) / _t
    # for k in scores_dict.keys():
    #     scores_dict[k]['bio_score'] = scores_dict[k]['freq_score'] * _r1 + \
    #                           scores_dict[k]['pred_score'] * _r2 + \
    #                           scores_dict[k]['report_score'] * _r3

    return scores_dict


# 表型分析打分
def pheno_scores(phenotype):
    pre_phe_list = phenotype.split(';')
    fin_phe_list = []
    pheno_score_dict = {}
    # 形式：
    # 0000238
    # Hydrocephalus
    # 0000238：Hydrocephalus[神经系统异常]
    # Hydrocephalus[神经系统异常]
    for i in pre_phe_list:
        if i.endswith(']'):
            _tmp = i.split('[')[0]
            if _tmp.find('：') != -1:
                fin_phe_list += ('HP:' + _tmp).split('：')
            else:
                fin_phe_list.append(_tmp)
        else:
            fin_phe_list.append(i)
    pheno_str = ';'.join(fin_phe_list)
    print(pheno_str)
    pheno_tmp_out = Tmp_dir + '/pheno_tmp'
    pheno_cmd = 'perl ./resource/src/phenolyzer-master/disease_annotation.pl "%s" -p -ph -logistic -out %s' % (
        pheno_str, pheno_tmp_out)
    result = os.system(pheno_cmd)
    if result:
        msg = '[ E: Something wrong with analysis phenotype ！]'
        print(msg)
        return msg
    else:
        with open(pheno_tmp_out + '.final_gene_list') as f:
            _ = f.readline()
            for _ in f:
                pheno_score_dict[_.split('\t')[1]] = eval(_.split('\t')[3]) * 100
        return pheno_score_dict


# 综合打分
def total_scores(in_file, file_out_name, bio_score, ph_score, total_rate=total_rate, extra_rate=extra_rate,
                 genetic_rate=genetic_rate):
    file_out = Out_dir + '/' + file_out_name + '.txt'
    outf = open(file_out, 'w')
    # 输出文件头
    header_file = Tmp_dir + '/header.txt'
    if os.path.isfile(header_file):
        os.remove(header_file)
    get_header_cmd = 'head -n 1 %s > %s' % (in_file, header_file)
    os.system(get_header_cmd)
    f = open(header_file)
    line = f.readline()
    outf.write('#Gene_score\ttotal_score\tBioinfo_score\tphenotype_score\t' + line)
    # 表型打分
    if not ph_score:
        for i in bio_score.keys():
            bio_score[i]['pheno_score'] = 0
            bio_score[i]['all_score'] = bio_score[i]['bio_score']
    else:
        _t = sum([eval(i) for i in total_rate.split(':')])
        _r1 = eval(total_rate.split(':')[0]) / _t
        _r2 = eval(total_rate.split(':')[1]) / _t
        pheno_max = 0
        for i in bio_score.keys():
            _tmp = 0
            for gene in bio_score[i]['gene'].split(','):
                if ph_score.get(gene, 0) >= _tmp:
                    _tmp = ph_score.get(gene, 0)
                bio_score[i]['pheno_score'] = _tmp
            if _tmp >= pheno_max:
                pheno_max = _tmp
        for i in bio_score.keys():
            bio_score[i]['pheno_score'] = bio_score[i].get('pheno_score', 0) / pheno_max * 100
            bio_score[i]['all_score'] = bio_score[i]['bio_score'] * _r1 + bio_score[i].get('pheno_score', 0) * _r2
    # 基因打分
    gene_score = {}
    extra_score = {}
    for i in bio_score.keys():
        gene_name = bio_score[i]['gene']
        gene_score[gene_name] = gene_score.get(gene_name, []) + [bio_score[i]['all_score']]
    for k in gene_score.keys():
        if len(gene_score[k]) == 1:
            extra_score[k] = gene_score[k][0]
        elif len(gene_score[k]) == 2:
            tmp_max = max(gene_score[k])
            gene_score[k].remove(tmp_max)
            tmp_max2 = gene_score[k][0]
            extra_score[k] = tmp_max + tmp_max2 * extra_rate
        else:
            tmp_max = max(gene_score[k])
            gene_score[k].remove(tmp_max)
            tmp_max2 = max(gene_score[k])
            extra_score[k] = tmp_max + tmp_max2 * extra_rate
    # max_gene = 0
    for i in bio_score.keys():
        gene_name = bio_score[i]['gene']
        if bio_score[i]['genetic_score']:
            bio_score[i]['gene_score'] = extra_score[gene_name] * (1 + genetic_rate)
        else:
            bio_score[i]['gene_score'] = extra_score[gene_name]
        # if bio_score[i]['gene_score'] >= max_gene:
        #     max_gene = bio_score[i]['gene_score']
    # 输出
    for k in sorted(bio_score.items(), key=lambda j: (j[1]['gene_score'], j[1]['all_score']), reverse=True):
        Total_score = str(round(bio_score[k[0]]['all_score'], 4))
        Pheno_score = str(round(bio_score[k[0]]['pheno_score'], 4))
        Bio_score = str(round(bio_score[k[0]]['bio_score'], 4))
        anno = bio_score[k[0]]['anno']
        Gene_score = str(round(bio_score[k[0]]['gene_score'], 4))
        outf.write(
            Gene_score + '\t' + Total_score + '\t' + Bio_score + '\t' + Pheno_score + '\t' + k[0] + '\t' + anno + '\n')
    outf.close()

    clear_tmp_cmd = 'rm -rf %s/* ' % Tmp_dir
    result = os.system(clear_tmp_cmd)
    if result:
        msg = '[ E: Something wrong with clear temp file ！]'
        print(msg)
        return msg
    else:
        msg = '运行完成 ！'
        print(msg)
        return msg
