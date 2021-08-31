import os
import gzip
from script.common import execute_system, get_row_num
from script.filter import frequency_filter


# 模块1：对vcf文件的格式转换
# 模块2：针对snvs和indels的频率注释
# 模块3：针对snvs和indels的有害性noncoding位点注释
# 模块4：针对snvs和indels的所有的注释
# 模块5：针对SVs的注释

def trio_short_variants_filter(vcf,
                               outDir,
                               scriptPath,
                               AFDb,
                               AFType,
                               geneDb,
                               geneType,
                               AFlist,
                               clinList,
                               afThreshold,
                               annoDir,
                               ref_version,
                               thread):
    # vcf -> avinput
    avinput, gt = short_variants_transfer_format(vcf, outDir, scriptPath)
    afAnno, num = anno_frequency(avinput, gt, outDir, AFDb, AFType, annoDir, thread, ref_version, scriptPath)
    afFilted = outDir + '/AF_filted.txt'
    frequency_filter(afAnno, afFilted, AFlist, afThreshold)
    geneAnno, num = anno_gene(afFilted, outDir, num, geneDb, geneType, annoDir, thread, ref_version, scriptPath)



def single_short_variants_filter(vcf):
    pass


def trio_structure_variants_filter(vcf):
    pass


def single_structure_variants_filter(vcf):
    pass


def short_variants_transfer_format(vcf, outDir, scriptPath):
    """paste input1.anno.vcf genotype > input.anno.vcf"""
    # out_dir = out_dir.rstrip('/') + '/'
    # file_name = os.path.basename(vcf_file).split('.vcf')[0]
    avinput = outDir + '/cohort.avinput'
    gt_file_tmp = outDir + '/cohort.gt_tmp'
    gt_file = outDir + '/cohort.gt'
    # 转换格式的时候保留vcf信息
    transCmd = 'perl %s/bin/annovar/convert2annovar.pl -format vcf4 %s -outfile %s -includeinfo -allsample -withfreq' \
               % (scriptPath, vcf, avinput)
    execute_system(transCmd, '[ Msg: Transfer format done ! ]', '[ E: Something wrong with transfer format ! ]')
    # 截取前五列以后的信息
    get_genotype = 'cut -f 14- %s > %s' % (avinput, gt_file_tmp)
    execute_system(get_genotype, '[ Msg: Get genotype done ! ]', '[ E: Something wrong with get genotype ! ]')
    # 截取头信息
    head = open(outDir + '/cohort.head.txt', 'w')
    with gzip.open(vcf, 'rb') as f:
        for i in f:
            if i.startswith('#CHROM'):
                head.write(i.split('\tALT\t')[1])
                break
    head.close()
    os.system('cat %s %s > %s' % (outDir + '/cohort.head.txt', gt_file_tmp, gt_file))

    return avinput, gt_file


def anno_frequency(avinput, gt, out_dir, db_list, type_list, anno_dir, thread, ver, scriptPath):
    # 注释
    fileName = os.path.basename(avinput).split('.')[0]
    file_out = out_dir.rstrip('/') + '/' + fileName
    annoCmd = 'perl %s/bin/table_annovar.pl %s %s --buildver %s -out %s -remove -protocol %s -operation %s -nastring ' \
              '. --thread %d > /dev/null 2>&1' % (
                  scriptPath, avinput, anno_dir, ver, file_out, db_list, type_list, thread)
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


def anno_sv():
    pass
