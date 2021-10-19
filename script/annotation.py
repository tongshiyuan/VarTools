import os
import gzip
from script.common import execute_system, get_row_num
from script.filter import frequency_filter


# 模块1：对vcf文件的格式转换
# 模块2：针对snvs和indels的频率注释
# 模块3：针对snvs和indels的有害性noncoding位点注释
# 模块4：针对snvs和indels的所有的注释
# 模块5：针对SVs的注释
def short_variants_convert_format(vcf, prefix, out_dir, script_path):
    avinput = '%s/%s_annovar.avinput' % (out_dir, prefix)
    info_tmp_file = '%s/%s_tmp.info' % (out_dir, prefix)
    info_file = '%s/%s_sample.info' % (out_dir, prefix)
    head_file = '%s/%s_head.txt' % (out_dir, prefix)
    # 转换格式的时候保留vcf信息
    cv_cmd = 'perl %s/bin/annovar/convert2annovar.pl -format vcf4 %s -outfile %s -includeinfo -allsample -withfreq' % (
        script_path, vcf, avinput)
    execute_system(cv_cmd, '[ Msg: Transfer format done ! ]', '[ Error: Something wrong with transfer format ! ]')
    # 截取前五列以后的信息
    get_info_cmd = 'cut -f 9- %s > %s' % (avinput, info_tmp_file)
    execute_system(get_info_cmd, '[ Msg: Get genotype done ! ]', '[ Error: Something wrong with get genotype ! ]')
    # 截取头信息
    head = open(head_file, 'w')
    if vcf.endswith('vcf.gz'):
        with gzip.open(vcf, 'rb') as f:
            for i in f:
                if i.decode().startswith('#CHROM'):
                    head.write(i.split('\tALT\t')[1])
                    break
    else:
        with open(vcf) as f:
            for i in f:
                if i.startswith('#CHROM'):
                    head.write(i.split('\tALT\t')[1])
                    break
    head.close()
    os.system('cat %s %s > %s' % (head_file, info_tmp_file, info_file))
    return avinput, info_file


def short_variants_anno():
    pass


def short_variants_step_anno():
    pass


def anno_frequency(avinput, gt, out_dir, db_list, type_list, anno_dir, thread, ver, script_path):
    # 注释
    fileName = os.path.basename(avinput).split('.')[0]
    file_out = out_dir.rstrip('/') + '/' + fileName
    annoCmd = 'perl %s/bin/table_annovar.pl %s %s --buildver %s -out %s -remove -protocol %s -operation %s -nastring ' \
              '. --thread %d > /dev/null 2>&1' % (
                  script_path, avinput, anno_dir, ver, file_out, db_list, type_list, thread)
    execute_system(annoCmd, '[ Msg: <%s> frequency annotation done ! ]' % fileName,
                   '[ E: Something wrong with <%s> frequency annotation ! ]' % fileName)
    # 合并
    tmp_result = file_out + '.hg19_multianno.txt'
    result = file_out + '.freq_anno'
    paste_cmd = 'paste %s %s > %s' % (tmp_result, gt, result)
    execute_system(paste_cmd, '[ Msg: <%s> frequency annotation and genotype paste done ! ]' % fileName,
                   '[ E: Something wrong with <%s> paste frequency annotation and genotype ! ]' % fileName)
    # 获得列数
    header = file_out + '.hg19_multianno.header'
    num = get_row_num(tmp_result, header)

    return result, num


def anno_gene(infile, out_dir, num, geneList, geneType, anno_dir, thread, ver, scriptPath):
    out_dir = out_dir.rstrip('/') + '/'
    file_name = os.path.basename(infile).split('.freq_anno')[0]
    # get genotyping
    gt_file = out_dir + file_name + '.gene_gt'
    get_genotype = 'cut -f %d- %s > %s' % (num, infile, gt_file)
    execute_system(get_genotype, '[ Msg: Get genotype done ! ]',
                   '[ E: Something wrong with get genotype ! ]')
    # anno
    file_out = out_dir + file_name + '_gene_tmp'
    anno_cmd = 'perl %s/bin/table_annovar_splicing.pl %s %s --buildver %s -out %s -remove -protocol %s -operation %s ' \
               '-nastring . --thread %d > /dev/null 2>&1' % (
                   scriptPath, infile, anno_dir, ver, file_out, geneList, geneType, thread)
    execute_system(anno_cmd, '[ Msg: Gene annotation done ! ]',
                   '[ E: Something wrong with Gene annotation ! ]')
    # 合并
    tmp_result = file_out + '.hg19_multianno.txt'
    result = out_dir + file_name + '.gene_anno'
    paste_cmd = 'paste %s %s > %s' % (tmp_result, gt_file, result)
    execute_system(paste_cmd, '[ Msg: Gene annotation and genotype paste done ! ]',
                   '[ E: Something wrong with paste gene annotation and genotype ! ]')
    # 获得列数
    header = file_out + '.hg19_multianno.header'
    num = get_row_num(tmp_result, header)

    return result, num


def anno_all(infile, out_dir, db_list, type_list, anno_dir, anno_thread=12, ver='hg19'):
    num = 6
    out_dir = out_dir.rstrip('/') + '/'
    file_name = os.path.basename(infile).split('.gene_anno')[0]
    # get genotyping
    gt_file = out_dir + file_name + '.all_gt'
    get_genotype = 'cut -f %d- %s > %s' % (num, infile, gt_file)
    execute_system(get_genotype, '[ Msg: <%s> get genotype done ! ]' % file_name,
                   '[ E: Something wrong with <%s> get genotype ! ]' % file_name)
    # anno
    file_out = out_dir + file_name + '_all_tmp'
    anno_cmd = 'perl ./bin/table_annovar_splicing.pl %s %s --buildver %s -out %s -remove -protocol %s -operation %s ' \
               '-nastring . --thread %d > /dev/null 2>&1' % (
                   infile, anno_dir, ver, file_out, db_list, type_list, anno_thread)
    execute_system(anno_cmd, '[ Msg: <%s> annotation done ! ]' % file_name,
                   '[ E: Something wrong with <%s> annotation ! ]' % file_name)
    # 合并
    tmp_result = file_out + '.hg19_multianno.txt'
    result = out_dir + file_name + '.all_anno'
    paste_cmd = 'paste %s %s > %s' % (tmp_result, gt_file, result)
    execute_system(paste_cmd, '[ Msg: <%s> annotation and genotype paste done ! ]' % file_name,
                   '[ E: Something wrong with <%s> paste annotation and genotype ! ]' % file_name)

    return result


def trio_short_variants_filter(vcf, prefix, out_dir, script_path,
                               AF_db, AF_type, AF_list, af_threshold,
                               geneDb, geneType, clin_list,
                               anno_dir, ref_version, thread):
    # vcf -> avinput
    avinput, info = short_variants_convert_format(vcf, prefix, out_dir, script_path)
    af_anno, num = anno_frequency(avinput, info, out_dir, AF_db, AF_type, anno_dir, thread, ref_version, script_path)
    af_filted = out_dir + '/AF_filted.txt'
    frequency_filter(af_anno, af_filted, AF_list, af_threshold)
    geneAnno, num = anno_gene(af_filted, out_dir, num, geneDb, geneType, anno_dir, thread, ref_version, script_path)


def single_short_variants_filter(vcf):
    pass


def trio_structure_variants_filter(vcf):
    pass


def single_structure_variants_filter(vcf):
    pass


def anno_sv():
    pass
