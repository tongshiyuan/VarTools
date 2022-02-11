import sys
import pandas as pd
from script.common import execute_system


def frequency_filter(infile, out_file, fredb_list, threshold_list, retain_line):
    th_list = [eval(i) for i in threshold_list.strip(',').split(',')]
    fredb_list = fredb_list.strip(',').split(',')
    if len(th_list) == 1:
        th_list = th_list * len(fredb_list)
        print(f'[ Msg: Start frequency filter, threshold: < {th_list[0]} > .]')
    else:
        if len(fredb_list) != len(th_list):
            return
        else:
            print('[ Msg: Start frequency filter, with muliple threshold.]')
    freq_index = []
    out_ = open(out_file, 'w')
    with open(infile) as f:
        line = f.readline()
        out_.write(line)
        col_list = line.strip().split('\t')
        for _db in fredb_list:
            try:
                _idx = col_list.index(_db)
                freq_index.append(_idx)
                col_list[_idx] = _db + '_indexed_Oxo_'
            except:
                print('[ Warn: Can not find database <%s> in frequency filter ! ]' % _db)
        if not freq_index:
            sys.exit('[ Error: Can not find allele frequency database in frequency filter ! ]')
        db_col, db_retain = [], []
        if retain_line:
            # retail_col1:retain_tag1|retain_tag2;retail_col2:retain_tag1|retain_tag2
            for db_ln in retain_line.strip().rstrip(';').split(';'):
                try:
                    col = db_ln.strip().split(':')[0]
                    col_idx = line.strip().split('\t').index(col)
                    db_col.append(col_idx)
                    db_retain.append(db_ln.strip().split(':')[1].strip().split('|'))
                except:
                    print('[ Warn: Something wrong with parsing retain expressions <%s> ! ]' % col)
            if not db_col:
                print('[ Warn: Program can not identify retain database in frequency filter ! ]')
        for line in f:
            records = line.strip().split('\t')
            if db_col:
                for idx, _ in enumerate(db_col):
                    if records[_] in db_retain[idx]:
                        out_.write(line)
                        break
                else:
                    for af_index, threshold in zip(freq_index, th_list):
                        if records[af_index] != '.' and eval(records[af_index]) >= threshold:
                            break
                    else:
                        out_.write(line)
            else:
                for af_index, threshold in zip(freq_index, th_list):
                    if records[af_index] != '.' and eval(records[af_index]) >= threshold:
                        break
                else:
                    out_.write(line)
    out_.close()
    print('[ Msg: Frequency filter done !]')


def exonic_filter(infile, out_file, gene_db, retain_line):
    print('[ Msg: Start exonic filter. ]')
    keys, fun_index = [], []
    out_ = open(out_file, 'w')
    with open(infile) as f:
        line = f.readline()
        out_.write(line)
        for _db in gene_db.strip().strip(',').split(','):
            try:
                fun_index.append(line.strip().split('\t').index('Func.' + _db))
            except:
                print('[ Warn: Can not find database <Func.%s> in exonic filter ! ]' % _db)
        if not fun_index:
            sys.exit('[ Error: Can not find gene function database in exonic filter ! ]')
        db_col, db_retain = [], []
        if retain_line:
            for db_ln in retain_line.strip().rstrip(';').split(';'):
                try:
                    col = db_ln.strip().split(':')[0]
                    col_idx = line.strip().split('\t').index(col)
                    db_col.append(col_idx)
                    db_retain.append(db_ln.strip().split(':')[1].strip().split('|'))
                except:
                    print('[ Warn: Something wrong with parsing retain expressions <%s> ! ]' % col)
            if not db_col:
                print('[ Warn: Program can not identify retain database in exonic filter ! ]')
        for line in f:
            records = line.strip().split('\t')
            for i in fun_index:
                if records[i].find('exonic') != -1 or records[i].find('splicing') != -1 or records[i] == 'intron':
                    out_.write(line)
                    break
            else:
                if db_col:
                    for idx, _ in enumerate(db_col):
                        if records[_] in db_retain[idx]:
                            out_.write(line)
                            break
    out_.close()
    print('[ Msg: Exonic filter done !]')


def parent_filter(file, outfile, probend_sample_name='', father_sample_name='', mother_sample_name=''):
    with open(file) as f:
        for line in f:
            if line.startswith('##'):
                pass
            else:
                header = line
                break
    num_left = len(header.strip().split('\tFORMAT\t')[0].split('\t'))
    header_list = header.strip().split('\tFORMAT\t')[1].split('\t')
    if not (father_sample_name and mother_sample_name):
        for i, j in enumerate(header_list):
            if 'F' in j or 'Fa' in j or 'Father' in j or 'f' in j or 'fa' in j or 'father' in j:
                fa_idx = i + num_left + 1
                break
        for i, j in enumerate(header_list):
            if 'M' in j or 'Mo' in j or 'Mother' in j or 'm' in j or 'mo' in j or 'mother' in j:
                mo_idx = i + num_left + 1
                break
    else:
        fa_idx = num_left + header_list.index(father_sample_name) + 1
        mo_idx = num_left + header_list.index(mother_sample_name) + 1
    if probend_sample_name:
        p_idx = num_left + header_list.index(probend_sample_name) + 1
    else:
        p_idx = num_left + 1
    of = open(outfile, 'w')
    with open(file) as f:
        head = f.readline()
        while head.startswith('##'):
            head = f.readline()
        for line in f:
            rec = line.strip().split('\t')
            pinfo = rec[p_idx]
            fainfo = rec[fa_idx]
            moinfo = rec[mo_idx]
            # 过滤条件
            l1 = pinfo.split(':')[0] in ['0/0', '0|0']
            l2 = fainfo.split(':')[0] in ['1/1', '1|1']
            l3 = moinfo.split(':')[0] in ['1/1', '1|1']
            if l1 or l2 or l3:
                pass
            else:
                of.write(line)
    of.close()
    print('[ Msg: Homozygote filtered done! ]')


def inheritance_pattern_filter():
    pass


def synonymous_filter(infile, out_file, db_list, clindb_dict, con_dict):
    print('[ Msg: Start synonymous filter. ]')
    # 变异类型
    type_index = []
    # 临床
    clin_index = []
    # 保守性
    con_index = []
    clin_keys = []
    con_keys = []
    out_ = open(out_file, 'w')
    with open(infile) as f:
        line = f.readline()
        out_.write(line)
        for _db in db_list:
            try:
                type_index.append(line.strip().split('\t').index(_db))
            except:
                print('[ W: Can not find database <%s> in synonymous filter ! ]' % _db)
        if not type_index:
            sys.exit('[ E: Can not find variant function database in synonymous filter ! ]')
        if clindb_dict:
            for _db in clindb_dict.keys():
                try:
                    clin_index.append(line.strip().split('\t').index(_db))
                    clin_keys.append(_db)
                except:
                    print('[ W: Can not find database <%s> in synonymous filter ! ]' % _db)
        else:
            print('[ W: Program have not reported database in synonymous filter ! ]')
        if con_dict:
            for _db in con_dict.keys():
                try:
                    con_index.append(line.strip().split('\t').index(_db))
                    con_keys.append(_db)
                except:
                    print('[ W: Can not find database <%s> in synonymous filter ! ]' % _db)
        else:
            print('[ W: Program have not regulation database in synonymous filter ! ]')
        for line in f:
            if line:
                records = line.strip().split('\t')
                for i in type_index:
                    if records[i] != 'synonymous SNV':
                        out_.write(line)
                        break
                else:
                    if clin_index:
                        for index, db in enumerate(clin_index):
                            if records[db] in clindb_dict[clin_keys[index]]:
                                out_.write(line)
                                break
                        else:
                            if con_index:
                                for ind, i in enumerate(con_index):
                                    if records[i] != '.' and eval(records[i]) > con_dict[con_keys[ind]]:
                                        out_.write(line)
                                        break
    out_.close()
    print('[ Msg: Synonymous filter done !]')


def genetic_filter(infile, out_file, trio_list, gander):
    print('[ Msg: Start genetic model filter. ]')
    start = time.perf_counter()
    out_ = open(out_file, 'w')
    with open(infile) as f:
        line = f.readline()
        out_.write(line)
        try:
            F = line.strip().split('\t').index(trio_list[2])
            M = line.strip().split('\t').index(trio_list[1])
            P = line.strip().split('\t').index(trio_list[0])
        except:
            sys.exit('[ E: Something wrong of sample index in genetic model filter ! ]')
        if gander in ['F', 'f', 'female', 'Female']:
            for line in f:
                if line:
                    records = line.strip().split('\t')
                    if records[P].startswith('1/1') and (records[F].startswith('1/1') or records[M].startswith('1/1')):
                        pass
                    elif records[P].startswith('0/1') and (
                            records[F].startswith('1/1') or records[M].startswith('1/1')):
                        pass
                    else:
                        out_.write(line)
        elif gander in ['M', 'm', 'male', 'Male']:
            for line in f:
                records = line.strip().split('\t')
                if records[0] == 'X':
                    if records[P].startswith('1/1') and (records[F].startswith('1/1') or records[M].startswith('1/1')):
                        pass
                    elif records[P].startswith('0/1') and (
                            records[F].startswith('1/1') or records[F].startswith('0/1')):
                        pass
                    else:
                        out_.write(line)
                elif records[0] == 'Y':
                    if records[P].startswith('1/1') and records[F].startswith('1/1'):
                        pass
                    else:
                        out_.write(line)
                else:
                    if records[P].startswith('1/1') and (records[F].startswith('1/1') or records[M].startswith('1/1')):
                        pass
                    elif records[P].startswith('0/1') and (
                            records[F].startswith('1/1') or records[M].startswith('1/1')):
                        pass
                    else:
                        out_.write(line)
    out_.close()
    end = time.perf_counter()
    print('[ Msg: Genetic model filter done ! use time : < %d > s ]' % (end - start))


def pheno_type():
    pass


def snv_indel_score():
    pass


def clinical_filter(infile, outdir, gender, reference):
    # infile: output/W001/family/gatk/W001.HC.filter.vcf.gz
    # outdir: output/W001
    tmp_dir = infile.split('/family/')[0]
    sampleName = get_sort_sampleName(os.listdir(tmp_dir))
    # W001
    familyName = outdir.split('/')[-1]
    # output/W001/filter
    outputdir = outdir + '/filter'
    os.makedirs(outputdir)
    # 对结果从多等位位点转换为单等位位点
    normCmd = 'bcftools norm -Ov -m-any -f ' + reference + ' ' + infile + ' > ' + \
              outputdir + '/' + familyName + '.norm.vcf'
    notice = 'something wrong with norm'
    execute_system(normCmd, notice)
    print(familyName, 'norm variations done!')
    # 提取gt信息
    gtCmd = 'grep -v "##" ' + outputdir + '/' + familyName + '.norm.vcf | cut -f1-9 --complement > ' + \
            outputdir + '/' + familyName + '.gt.txt'
    notice = 'something wrong with get gt information'
    execute_system(gtCmd, notice)
    print(familyName, 'get gt information done!')
    # 频率与基因注释
    aviCmd = 'perl ./bin/convert2annovar.pl -includeinfo -format vcf4old ' + \
             outputdir + '/' + familyName + '.norm.vcf > ' + \
             outputdir + '/' + familyName + '.avinput'
    notice = 'something wrong with format transition (vcf to avinput)'
    execute_system(aviCmd, notice)
    print(familyName, 'format transition done!')
    annoCmd = 'perl ./bin/table_annovar_splicing.pl ' + outputdir + '/' + familyName + '.avinput ' + annovardb + \
              ' --buildver hg19 -out ' + outputdir + '/' + familyName + \
              ' -remove -protocol refGene,EAS.sites.2015_08,ALL.sites.2015_08,kaviar_20150923,hrcr1,cg69,' + \
              'gnomad_genome,exac03,exac03nonpsych,esp6500siv2_all,cg46,omim201806,clinvar_20190305,hgmd ' + \
              '-operation g,f,f,f,f,f,f,f,f,f,f,r,f,f -nastring . --thread 12 > /dev/null 2>&1'
    notice = 'something wrong with frequency annotation'
    execute_system(annoCmd, notice)
    print(familyName, 'frequency annotation done!')
    # 合并结果
    pasteCmd = 'paste ' + outputdir + '/' + familyName + '.hg19_multianno.txt ' + \
               outputdir + '/' + familyName + '.gt.txt > ' + \
               outputdir + '/' + familyName + '.anno.txt'
    notice = 'something wrong with paste result'
    execute_system(pasteCmd, notice)
    print(familyName, 'paste result done!')
    # 搜素已经报道的可能致病位点
    searchCmd = 'grep ' + disease + ' ' + outputdir + '/' + familyName + '.anno.txt > ' + \
                outputdir + '/' + familyName + '.' + disease + '.txt'
    # notice = 'something wrong with search disease variations'
    # execute_system(searchCmd, notice)
    os.system(searchCmd)
    print(familyName, 'search disease variations done!')
    # 筛选
    freqCmd = './bin/annotools.py freq -i ' + outputdir + '/' + familyName + '.anno.txt -o ' + \
              outputdir + '/' + familyName + '.anno -t ' + threshold
    notice = 'something wrong with frequency filter'
    execute_system(freqCmd, notice)
    print(familyName, 'frequency filterdone!')
    exonCmd = './bin/annotools.py exon -i ' + outputdir + '/' + familyName + '.anno_freqFilter.txt -o ' + \
              outputdir + '/' + familyName + '.anno_freqFilter'
    notice = 'something wrong with coding/noncoding split'
    execute_system(exonCmd, notice)
    print(familyName, 'coding/noncoding split done!')
    trioCmd = './bin/annotools.py trio -i ' + outputdir + '/' + familyName + '.anno_freqFilter_coding.txt -o ' + \
              outputdir + '/' + familyName + '.anno_freqFilter_coding -g ' + gender + \
              ' -f ' + sampleName[1] + ' -m ' + sampleName[2] + ' -c ' + sampleName[0]
    notice = 'something wrong with trio analysis'
    execute_system(trioCmd, notice)
    print(familyName, 'trio analysis done!')
    # 备份结果，去除header
    Cmd = 'cp ' + outputdir + '/' + familyName + '.anno_freqFilter_coding_trio.txt ' + outputdir + '/' + familyName + \
          ".filted.aviput && sed -i '1d' " + outputdir + '/' + familyName + '.filted.aviput'
    notice = 'something wrong with input file prepare in annotation'
    execute_system(Cmd, notice)
    print(familyName, 'input file prepare done!')
    # 第二次注释
    annoCmd = 'perl ./bin/table_annovar_splicing.pl ' + outputdir + '/' + familyName + '.filted.aviput ' + annovardb + \
              ' --buildver hg19 -remove -out ' + outputdir + '/' + familyName + \
              '.filted -protocol refGene,cytoBand,omim201806,gwasCatalog,CCRS,phastConsElements100way,' + \
              'tfbsConsSites,rmsk,pseudogene,' + \
              'snp151,clinvar_20190305,hgmd,gerp++gt2,spidex,dbscsnv11,revel,mcap13,dbnsfp35a,' + \
              'ReVe,PrimateAI,eigen,gwava,avsnp150,regsnpintron,scape,scap3cd,scap3cr,scap3i,scap5cd,scap5cr,' + \
              'scap5i,scap5s,intervar_20180118,clinpred ' + \
              '-operation g,r,r,r,r,r,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f ' + \
              '-nastring . --thread 12 > /dev/null 2>&1'
    notice = 'something wrong with annotation'
    execute_system(annoCmd, notice)
    print(familyName, 'annotation done!')
    # 截取gt以及遗传模式
    sampleNum = len(sampleName)
    with open(outputdir + '/' + familyName + '.anno_freqFilter_coding_trio.txt', 'r') as f:
        header = f.readline()
        rowNum = len(header.strip().split('\t'))
    gtNum = rowNum - sampleNum - 1
    gtCmd = 'cut -f1-' + str(
        gtNum) + ' --complement ' + outputdir + '/' + familyName + '.anno_freqFilter_coding_trio.txt > ' + \
            outputdir + '/' + familyName + '.filted.gt.txt'
    notice = 'something wrong with cut genotyping'
    execute_system(gtCmd, notice)
    print(familyName, 'cut genotyping done!')
    # paste 结果
    pasteCmd = 'paste ' + outputdir + '/' + familyName + '.filted.hg19_multianno.txt ' + \
               outputdir + '/' + familyName + '.filted.gt.txt > ' + \
               outputdir + '/' + familyName + '.filted.merge.txt'
    notice = 'something wrong with paste result'
    execute_system(pasteCmd, notice)
    print(familyName, 'paste result done!')
    # 对照过滤
    if sampleNum > 3:
        ctlCmd = './bin/annotools.py control -i ' + outputdir + '/' + familyName + '.filted.merge.txt -o ' + \
                 outputdir + '/' + familyName + '.filted.merge'
        notice = 'something wrong with control analysis'
        execute_system(ctlCmd, notice)
        print(familyName, 'control analysis done!')
    else:
        mvCmd = 'mv ' + outputdir + '/' + familyName + '.filted.merge.txt ' + \
                outputdir + '/' + familyName + '.filted.merge_ctl.txt'
        notice = 'something wrong with change name'
        execute_system(mvCmd, notice)
        print(familyName, 'change name done!')
    # 良性过滤
    benignCmd = './bin/annotools.py pathogenic -i ' + outputdir + '/' + familyName + '.filted.merge_ctl.txt -o ' + \
                outputdir + '/' + familyName
    notice = 'something wrong with benign filter'
    execute_system(benignCmd, notice)
    print(familyName, 'benign filter done!')
    # 非保守区的同义过滤
    synCmd = './bin/annotools.py filter -i ' + outputdir + '/' + familyName + '_harmful.txt -o ' + \
             outputdir + '/' + familyName + '.harmful'
    notice = 'something wrong with synonymous SNV filter'
    execute_system(synCmd, notice)
    print(familyName, 'synonymous SNV filter done!')
    # 可变剪切
    spliceCmd = './bin/annotools.py splice -i ' + outputdir + '/' + familyName + '.harmful_myfilter.txt -o ' + \
                outputdir + '/' + familyName + '.harmful.syn -c ' + sampleName[0] + ' -g ' + gender
    notice = 'something wrong with splicing filter'
    execute_system(spliceCmd, notice)
    print(familyName, 'splicing filter done!')
    print('')
    print(familyName, 'basic filter done!')
    # 删除过程文件
    rmCmd = 'rm -f ' + outputdir + '/' + familyName + '.norm.vcf ' + \
            outputdir + '/' + familyName + '.hg19_multianno.txt ' + \
            outputdir + '/' + familyName + '.avinput ' + \
            outputdir + '/' + familyName + '*.gt.txt ' + \
            outputdir + '/' + familyName + '*.invalid_input'
    notice = 'something wrong with delete process files'
    execute_system(rmCmd, notice)
    print(familyName, 'delete process files done !')

    return outputdir + '/' + familyName + '.harmful.syn.spliced_comhet.txt'


def tran_format(file, outfile, bp=100):
    output = open(outfile + '.region.bed', "w")
    with open(file, 'r') as f:
        whole = f.readline()
        if whole.rstrip().split('\t')[0].startswith("Chr"):
            whole = f.readline()
        while whole:
            line = whole.rstrip().split('\t')
            output.write(
                line[0].split('hr')[-1] + '\t' + str(eval(line[1]) - bp) + '\t' + str(eval(line[1]) + bp) + '\n')
            whole = f.readline()
    output.close()
    print("bed file format done!")


def merge(df, distance=0):
    # 合并需求：在miao相同的情况下
    # 1.有交集合并: chr1's max > chr2's min
    #            如果 chr1's max < chr2's max :
    #                二者合并为chr1's min , chr2's max
    #            否则：
    #                chr1's min , chr1's max
    # 2.两集合间距离小于x合并:
    #       0 < chr2's start - chr1's end < x
    # 合并后，删除后一行，合并至前一行，掩码实现。
    # 或者 ：chr2's start - chr1's end < x
    i = 0
    length = len(df.index)
    while (i < length - 1):
        if df.loc[i, 'chr'] == df.loc[i + 1, 'chr']:
            if df.loc[i + 1, 'start'] - df.loc[i, 'end'] <= distance:
                if df.loc[i, 'end'] < df.loc[i + 1, 'end']:  # 下一项的end大于此项，使用下一项的end
                    df.loc[i, 'end'] = df.loc[i + 1, 'end']
                df = df.drop([i + 1], axis=0)  # 删除合并的行
                length -= 1  # df长度变短
                df.index = range(length)  # 调整index
                i -= 1

        i += 1
    return df


def sort_file(outfile, distance=1000):
    # path = 'tong.csv'  # 文件路径
    path = outfile + '.region.bed'
    columns = ['chr', 'start', 'end']  # 列名
    df = pd.read_csv(path, header=None, sep='\t')  # 读取文件
    df.columns = columns
    df = df.astype(str)
    for i in range(0, len(df.index)):
        # print(df.loc[i,'chr'])
        # print(type(df.loc[i,'chr']))
        if (df.loc[i, 'chr'] == 'X'):  # 调整x为100，用于排序
            df.loc[i, 'chr'] = 23
        elif (df.loc[i, 'chr'] == 'Y'):  # 调整y为101，用于排序
            df.loc[i, 'chr'] = 24
        elif (df.loc[i, 'chr'].find('_') != -1):
            if df.loc[i, 'chr'] == '1_gl000191_random':
                df.loc[i, 'chr'] = 30
            elif df.loc[i, 'chr'] == '1_gl000192_random':
                df.loc[i, 'chr'] = 31
            elif df.loc[i, 'chr'] == '4_ctg9_hap1':
                df.loc[i, 'chr'] = 32
            elif df.loc[i, 'chr'] == '4_gl000193_random':
                df.loc[i, 'chr'] = 33
            elif df.loc[i, 'chr'] == '4_gl000194_random':
                df.loc[i, 'chr'] = 34
            elif df.loc[i, 'chr'] == '6_apd_hap1':
                df.loc[i, 'chr'] = 35
            elif df.loc[i, 'chr'] == '6_cox_hap2':
                df.loc[i, 'chr'] = 36
            elif df.loc[i, 'chr'] == '6_dbb_hap3':
                df.loc[i, 'chr'] = 37
            elif df.loc[i, 'chr'] == '6_mann_hap4':
                df.loc[i, 'chr'] = 38
            elif df.loc[i, 'chr'] == '6_mcf_hap5':
                df.loc[i, 'chr'] = 39
            elif df.loc[i, 'chr'] == '6_qbl_hap6':
                df.loc[i, 'chr'] = 40
            elif df.loc[i, 'chr'] == '6_ssto_hap7':
                df.loc[i, 'chr'] = 41
            elif df.loc[i, 'chr'] == '7_gl000195_random':
                df.loc[i, 'chr'] = 42
            elif df.loc[i, 'chr'] == '8_gl000196_random':
                df.loc[i, 'chr'] = 43
            elif df.loc[i, 'chr'] == '8_gl000197_random':
                df.loc[i, 'chr'] = 44
            elif df.loc[i, 'chr'] == '9_gl000198_random':
                df.loc[i, 'chr'] = 45
            elif df.loc[i, 'chr'] == '9_gl000199_random':
                df.loc[i, 'chr'] = 46
            elif df.loc[i, 'chr'] == '9_gl000200_random':
                df.loc[i, 'chr'] = 47
            elif df.loc[i, 'chr'] == '9_gl000201_random':
                df.loc[i, 'chr'] = 48
            elif df.loc[i, 'chr'] == '11_gl000202_random':
                df.loc[i, 'chr'] = 49
            elif df.loc[i, 'chr'] == '17_ctg5_hap1':
                df.loc[i, 'chr'] = 50
            elif df.loc[i, 'chr'] == '17_gl000203_random':
                df.loc[i, 'chr'] = 51
            elif df.loc[i, 'chr'] == '17_gl000204_random':
                df.loc[i, 'chr'] = 52
            elif df.loc[i, 'chr'] == '17_gl000205_random':
                df.loc[i, 'chr'] = 53
            elif df.loc[i, 'chr'] == '17_gl000206_random':
                df.loc[i, 'chr'] = 54
            elif df.loc[i, 'chr'] == '18_gl000207_random':
                df.loc[i, 'chr'] = 55
            elif df.loc[i, 'chr'] == '19_gl000208_random':
                df.loc[i, 'chr'] = 56
            elif df.loc[i, 'chr'] == '19_gl000209_random':
                df.loc[i, 'chr'] = 57
            elif df.loc[i, 'chr'] == '21_gl000210_random':
                df.loc[i, 'chr'] = 58
            elif df.loc[i, 'chr'] == 'Un_gl000211':
                df.loc[i, 'chr'] = 59
            elif df.loc[i, 'chr'] == 'Un_gl000212':
                df.loc[i, 'chr'] = 60
            elif df.loc[i, 'chr'] == 'Un_gl000213':
                df.loc[i, 'chr'] = 61
            elif df.loc[i, 'chr'] == 'Un_gl000214':
                df.loc[i, 'chr'] = 62
            elif df.loc[i, 'chr'] == 'Un_gl000215':
                df.loc[i, 'chr'] = 63
            elif df.loc[i, 'chr'] == 'Un_gl000216':
                df.loc[i, 'chr'] = 64
            elif df.loc[i, 'chr'] == 'Un_gl000217':
                df.loc[i, 'chr'] = 65
            elif df.loc[i, 'chr'] == 'Un_gl000218':
                df.loc[i, 'chr'] = 66
            elif df.loc[i, 'chr'] == 'Un_gl000219':
                df.loc[i, 'chr'] = 67
            elif df.loc[i, 'chr'] == 'Un_gl000220':
                df.loc[i, 'chr'] = 68
            elif df.loc[i, 'chr'] == 'Un_gl000221':
                df.loc[i, 'chr'] = 69
            elif df.loc[i, 'chr'] == 'Un_gl000222':
                df.loc[i, 'chr'] = 70
            elif df.loc[i, 'chr'] == 'Un_gl000223':
                df.loc[i, 'chr'] = 71
            elif df.loc[i, 'chr'] == 'Un_gl000224':
                df.loc[i, 'chr'] = 72
            elif df.loc[i, 'chr'] == 'Un_gl000225':
                df.loc[i, 'chr'] = 73
            elif df.loc[i, 'chr'] == 'Un_gl000226':
                df.loc[i, 'chr'] = 74
            elif df.loc[i, 'chr'] == 'Un_gl000227':
                df.loc[i, 'chr'] = 75
            elif df.loc[i, 'chr'] == 'Un_gl000228':
                df.loc[i, 'chr'] = 76
            elif df.loc[i, 'chr'] == 'Un_gl000229':
                df.loc[i, 'chr'] = 77
            elif df.loc[i, 'chr'] == 'Un_gl000230':
                df.loc[i, 'chr'] = 78
            elif df.loc[i, 'chr'] == 'Un_gl000231':
                df.loc[i, 'chr'] = 79
            elif df.loc[i, 'chr'] == 'Un_gl000232':
                df.loc[i, 'chr'] = 80
            elif df.loc[i, 'chr'] == 'Un_gl000233':
                df.loc[i, 'chr'] = 81
            elif df.loc[i, 'chr'] == 'Un_gl000234':
                df.loc[i, 'chr'] = 82
            elif df.loc[i, 'chr'] == 'Un_gl000235':
                df.loc[i, 'chr'] = 83
            elif df.loc[i, 'chr'] == 'Un_gl000236':
                df.loc[i, 'chr'] = 84
            elif df.loc[i, 'chr'] == 'Un_gl000237':
                df.loc[i, 'chr'] = 85
            elif df.loc[i, 'chr'] == 'Un_gl000238':
                df.loc[i, 'chr'] = 86
            elif df.loc[i, 'chr'] == 'Un_gl000239':
                df.loc[i, 'chr'] = 87
            elif df.loc[i, 'chr'] == 'Un_gl000240':
                df.loc[i, 'chr'] = 88
            elif df.loc[i, 'chr'] == 'Un_gl000241':
                df.loc[i, 'chr'] = 89
            elif df.loc[i, 'chr'] == 'Un_gl000242':
                df.loc[i, 'chr'] = 90
            elif df.loc[i, 'chr'] == 'Un_gl000243':
                df.loc[i, 'chr'] = 91
            elif df.loc[i, 'chr'] == 'Un_gl000244':
                df.loc[i, 'chr'] = 92
            elif df.loc[i, 'chr'] == 'Un_gl000245':
                df.loc[i, 'chr'] = 93
            elif df.loc[i, 'chr'] == 'Un_gl000246':
                df.loc[i, 'chr'] = 94
            elif df.loc[i, 'chr'] == 'Un_gl000247':
                df.loc[i, 'chr'] = 95
            elif df.loc[i, 'chr'] == 'Un_gl000248':
                df.loc[i, 'chr'] = 96
            elif df.loc[i, 'chr'] == 'Un_gl000249':
                df.loc[i, 'chr'] = 97
    df = df.astype(int)  # data强制转'换为int型
    # 按列chr,start,end进行排序
    df = df.sort_values(by=columns, ascending=(1, 1, 1))
    # 重排index，由小到大
    df.index = range(0, len(df.index))
    df = merge(df, distance)
    for i in range(0, len(df.index)):
        if (df.loc[i, 'chr'] in range(1, 23)):  # 调整100为x，存回
            df.loc[i, 'chr'] = 'chr' + str(df.loc[i, 'chr'])
        if (df.loc[i, 'chr'] == 23):  # 调整100为x，存回
            df.loc[i, 'chr'] = 'chrX'
        elif (df.loc[i, 'chr'] == 24):  # 调整101为y，存回
            df.loc[i, 'chr'] = 'chrY'
        elif df.loc[i, 'chr'] == 30:
            df.loc[i, 'chr'] = 'chr1_gl000191_random'
        elif df.loc[i, 'chr'] == 31:
            df.loc[i, 'chr'] = 'chr1_gl000192_random'
        elif df.loc[i, 'chr'] == 32:
            df.loc[i, 'chr'] = 'chr4_ctg9_hap1'
        elif df.loc[i, 'chr'] == 33:
            df.loc[i, 'chr'] = 'chr4_gl000193_random'
        elif df.loc[i, 'chr'] == 34:
            df.loc[i, 'chr'] = 'chr4_gl000194_random'
        elif df.loc[i, 'chr'] == 35:
            df.loc[i, 'chr'] = 'chr6_apd_hap1'
        elif df.loc[i, 'chr'] == 36:
            df.loc[i, 'chr'] = 'chr6_cox_hap2'
        elif df.loc[i, 'chr'] == 37:
            df.loc[i, 'chr'] = 'chr6_dbb_hap3'
        elif df.loc[i, 'chr'] == 38:
            df.loc[i, 'chr'] = 'chr6_mann_hap4'
        elif df.loc[i, 'chr'] == 39:
            df.loc[i, 'chr'] = 'chr6_mcf_hap5'
        elif df.loc[i, 'chr'] == 40:
            df.loc[i, 'chr'] = 'chr6_qbl_hap6'
        elif df.loc[i, 'chr'] == 41:
            df.loc[i, 'chr'] = 'chr6_ssto_hap7'
        elif df.loc[i, 'chr'] == 42:
            df.loc[i, 'chr'] = 'chr7_gl000195_random'
        elif df.loc[i, 'chr'] == 43:
            df.loc[i, 'chr'] = 'chr8_gl000196_random'
        elif df.loc[i, 'chr'] == 44:
            df.loc[i, 'chr'] = 'chr8_gl000197_random'
        elif df.loc[i, 'chr'] == 45:
            df.loc[i, 'chr'] = 'chr9_gl000198_random'
        elif df.loc[i, 'chr'] == 46:
            df.loc[i, 'chr'] = 'chr9_gl000199_random'
        elif df.loc[i, 'chr'] == 47:
            df.loc[i, 'chr'] = 'chr9_gl000200_random'
        elif df.loc[i, 'chr'] == 48:
            df.loc[i, 'chr'] = 'chr9_gl000201_random'
        elif df.loc[i, 'chr'] == 49:
            df.loc[i, 'chr'] = 'chr11_gl000202_random'
        elif df.loc[i, 'chr'] == 50:
            df.loc[i, 'chr'] = 'chr17_ctg5_hap1'
        elif df.loc[i, 'chr'] == 51:
            df.loc[i, 'chr'] = 'chr17_gl000203_random'
        elif df.loc[i, 'chr'] == 52:
            df.loc[i, 'chr'] = 'chr17_gl000204_random'
        elif df.loc[i, 'chr'] == 53:
            df.loc[i, 'chr'] = 'chr17_gl000205_random'
        elif df.loc[i, 'chr'] == 54:
            df.loc[i, 'chr'] = 'chr17_gl000206_random'
        elif df.loc[i, 'chr'] == 55:
            df.loc[i, 'chr'] = 'chr18_gl000207_random'
        elif df.loc[i, 'chr'] == 56:
            df.loc[i, 'chr'] = 'chr19_gl000208_random'
        elif df.loc[i, 'chr'] == 57:
            df.loc[i, 'chr'] = 'chr19_gl000209_random'
        elif df.loc[i, 'chr'] == 58:
            df.loc[i, 'chr'] = 'chr21_gl000210_random'
        elif df.loc[i, 'chr'] == 59:
            df.loc[i, 'chr'] = 'chrUn_gl000211'
        elif df.loc[i, 'chr'] == 60:
            df.loc[i, 'chr'] = 'chrUn_gl000212'
        elif df.loc[i, 'chr'] == 61:
            df.loc[i, 'chr'] = 'chrUn_gl000213'
        elif df.loc[i, 'chr'] == 62:
            df.loc[i, 'chr'] = 'chrUn_gl000214'
        elif df.loc[i, 'chr'] == 63:
            df.loc[i, 'chr'] = 'chrUn_gl000215'
        elif df.loc[i, 'chr'] == 64:
            df.loc[i, 'chr'] = 'chrUn_gl000216'
        elif df.loc[i, 'chr'] == 65:
            df.loc[i, 'chr'] = 'chrUn_gl000217'
        elif df.loc[i, 'chr'] == 66:
            df.loc[i, 'chr'] = 'chrUn_gl000218'
        elif df.loc[i, 'chr'] == 67:
            df.loc[i, 'chr'] = 'chrUn_gl000219'
        elif df.loc[i, 'chr'] == 68:
            df.loc[i, 'chr'] = 'chrUn_gl000220'
        elif df.loc[i, 'chr'] == 69:
            df.loc[i, 'chr'] = 'chrUn_gl000221'
        elif df.loc[i, 'chr'] == 70:
            df.loc[i, 'chr'] = 'chrUn_gl000222'
        elif df.loc[i, 'chr'] == 71:
            df.loc[i, 'chr'] = 'chrUn_gl000223'
        elif df.loc[i, 'chr'] == 72:
            df.loc[i, 'chr'] = 'chrUn_gl000224'
        elif df.loc[i, 'chr'] == 73:
            df.loc[i, 'chr'] = 'chrUn_gl000225'
        elif df.loc[i, 'chr'] == 74:
            df.loc[i, 'chr'] = 'chrUn_gl000226'
        elif df.loc[i, 'chr'] == 75:
            df.loc[i, 'chr'] = 'chrUn_gl000227'
        elif df.loc[i, 'chr'] == 76:
            df.loc[i, 'chr'] = 'chrUn_gl000228'
        elif df.loc[i, 'chr'] == 77:
            df.loc[i, 'chr'] = 'chrUn_gl000229'
        elif df.loc[i, 'chr'] == 78:
            df.loc[i, 'chr'] = 'chrUn_gl000230'
        elif df.loc[i, 'chr'] == 79:
            df.loc[i, 'chr'] = 'chrUn_gl000231'
        elif df.loc[i, 'chr'] == 80:
            df.loc[i, 'chr'] = 'chrUn_gl000232'
        elif df.loc[i, 'chr'] == 81:
            df.loc[i, 'chr'] = 'chrUn_gl000233'
        elif df.loc[i, 'chr'] == 82:
            df.loc[i, 'chr'] = 'chrUn_gl000234'
        elif df.loc[i, 'chr'] == 83:
            df.loc[i, 'chr'] = 'chrUn_gl000235'
        elif df.loc[i, 'chr'] == 84:
            df.loc[i, 'chr'] = 'chrUn_gl000236'
        elif df.loc[i, 'chr'] == 85:
            df.loc[i, 'chr'] = 'chrUn_gl000237'
        elif df.loc[i, 'chr'] == 86:
            df.loc[i, 'chr'] = 'chrUn_gl000238'
        elif df.loc[i, 'chr'] == 87:
            df.loc[i, 'chr'] = 'chrUn_gl000239'
        elif df.loc[i, 'chr'] == 88:
            df.loc[i, 'chr'] = 'chrUn_gl000240'
        elif df.loc[i, 'chr'] == 89:
            df.loc[i, 'chr'] = 'chrUn_gl000241'
        elif df.loc[i, 'chr'] == 90:
            df.loc[i, 'chr'] = 'chrUn_gl000242'
        elif df.loc[i, 'chr'] == 91:
            df.loc[i, 'chr'] = 'chrUn_gl000243'
        elif df.loc[i, 'chr'] == 92:
            df.loc[i, 'chr'] = 'chrUn_gl000244'
        elif df.loc[i, 'chr'] == 93:
            df.loc[i, 'chr'] = 'chrUn_gl000245'
        elif df.loc[i, 'chr'] == 94:
            df.loc[i, 'chr'] = 'chrUn_gl000246'
        elif df.loc[i, 'chr'] == 95:
            df.loc[i, 'chr'] = 'chrUn_gl000247'
        elif df.loc[i, 'chr'] == 96:
            df.loc[i, 'chr'] = 'chrUn_gl000248'
        elif df.loc[i, 'chr'] == 97:
            df.loc[i, 'chr'] = 'chrUn_gl000249'
    df.to_csv(outfile + '.region.sorted.bed', sep="\t", header=None, index=None)
    print("bed file sort done!")


def cut_multi_task(infile, sample, path, bam, outdir, familyName, ref):
    cutCmd = 'bedtools intersect -a ' + infile + '/' + sample + path + '/' + sample + '.' + ref + \
             bam + ' -b ' + outdir + '/' + familyName + '.region.sorted.bed > ' + \
             outdir + '/' + sample + '.target_reads.bam && samtools index ' + \
             outdir + '/' + sample + '.target_reads.bam'
    notice = 'something wrong with cut bam or bulid index'
    execute_system(cutCmd, notice)
    print(sample, 'cut bam and bulid index done!')


def cut_region(varfile, infile, calltools, mapping, reference):
    # W001
    familyName = infile.split('/')[-1]
    sampleName = get_sort_sampleName(os.listdir(infile))
    ref = reference.split('/')[-1].split('.')[0]
    # output/W001/SIRegion
    outdir = infile + '/SIRegion'
    os.makedirs(outdir)
    tran_format(varfile, outdir + '/' + familyName)
    sort_file(outdir + '/' + familyName)
    if calltools == 'gatk':
        path = '/callSI/gatk'
        bam = '.sorted.merged.markdup.BQSR.bam'
    else:
        path = '/mapping/' + mapping
        bam = '.sorted.merged.markdup.bam'
    process = Pool(8)
    for sample in sampleName:
        process.apply_async(cut_multi_task, args=(infile, sample, path, bam, outdir, familyName, ref,))
    print('Waiting for all sample cut region done...')
    process.close()
    process.join()
    print(familyName + ' all sample cut bam and bulid index done!')
