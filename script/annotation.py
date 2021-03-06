import os
import gzip
from script.common import execute_system, get_row_num
from script.filter import frequency_filter, exonic_filter


def short_variants_convert_format(vcf, prefix, out_dir, script_path):
    avinput = f'{out_dir}/{prefix}_annovar.avinput'
    info_tmp_file = f'{out_dir}/{prefix}_tmp.info'
    info_file = f'{out_dir}/{prefix}_sample.info'
    head_file = f'{out_dir}/{prefix}_head.txt'
    # 转换格式的时候保留vcf信息
    cv_cmd = 'perl %s/bin/annovar/convert2annovar.pl -format vcf4 %s -outfile %s -includeinfo -allsample -withfreq' % (
        script_path, vcf, avinput)
    execute_system(cv_cmd, '[ Msg: Transfer format done ! ]', '[ Error: Something wrong with transfer format ! ]')
    # 截取前9列以后的信息(vcf)
    get_info_cmd = 'cut -f 9- %s > %s' % (avinput, info_tmp_file)
    execute_system(get_info_cmd, '[ Msg: Get genotype done ! ]', '[ Error: Something wrong with get genotype ! ]')
    # 截取头信息
    head = open(head_file, 'w')
    if vcf.endswith('vcf.gz'):
        with gzip.open(vcf, 'rb') as f:
            for i in f:
                if i.decode().startswith('#CHROM'):
                    head.write(i.decode())  # .split('\tALT\t')[1])
                    break
    else:
        with open(vcf) as f:
            for i in f:
                if i.startswith('#CHROM'):
                    head.write(i)  # .split('\tALT\t')[1])
                    break
    head.close()
    os.system('cat %s %s > %s' % (head_file, info_tmp_file, info_file))
    return avinput, info_file


def anno_frequency(avinput, info, out_dir, prefix, db_list, type_list, anno_dir, ver, script_path, thread):
    # 注释
    file_out = out_dir + '/' + prefix
    anno_cmd = 'perl %s/bin/annovar/table_annovar.pl %s %s --buildver %s -out %s -remove ' \
               '-protocol %s -operation %s -nastring . --thread %d > /dev/null 2>&1' % (
                   script_path, avinput, anno_dir, ver, file_out, db_list, type_list, thread)
    execute_system(anno_cmd, '[ Msg: <%s> frequency annotation done ! ]' % prefix,
                   '[ Error: Something wrong with <%s> frequency annotation ! ]' % prefix)
    # 合并
    tmp_result = f'{file_out}.{ver}_multianno.txt'
    result = file_out + '.freq_anno'
    paste_cmd = f'paste {tmp_result} {info} > {result}'
    execute_system(paste_cmd, '[ Msg: <%s> frequency annotation and genotype paste done ! ]' % prefix,
                   '[ Error: Something wrong with <%s> paste frequency annotation and genotype ! ]' % prefix)
    # 获得列数
    header = f'{file_out}.{ver}_multianno.header'
    num = get_row_num(tmp_result, header)

    return result, num + 1


def anno_db(infile, info_file, out_dir, prefix, postfix, db_list, type_list, anno_dir, ver, script_path, thread,
            splice_distance, cn=True, first=True):
    file_out = out_dir + '/' + prefix
    # get genotyping
    # info_file = out_dir + '/' + prefix + '.info'
    # get_info = 'cut -f %d- %s > %s' % (num, infile, info_file)
    # execute_system(get_info, '[ Msg: Get genotype done ! ]', '[ Error: Something wrong with get genotype ! ]')
    # anno
    if first:
        cut_cmd = f"sed -i '1d' {infile}"
        execute_system(cut_cmd)
    anno_cmd = 'perl %s/bin/annovar/table_annovar.pl %s %s --buildver %s -out %s -remove -protocol %s -operation %s ' \
               '-nastring . --thread %d --intronhgvs %d > /dev/null 2>&1' % (
                   script_path, infile, anno_dir, ver, file_out, db_list, type_list, thread, splice_distance)
    execute_system(anno_cmd, '[ Msg: Gene annotation done ! ]', '[ Error: Something wrong with Gene annotation ! ]')
    # 合并
    tmp_result = '%s.%s_multianno.txt' % (file_out, ver)
    result = file_out + postfix
    paste_cmd = 'paste %s %s > %s' % (tmp_result, info_file, result)
    execute_system(paste_cmd, '[ Msg: Gene annotation and genotype paste done ! ]',
                   '[ Error: Something wrong with paste gene annotation and genotype ! ]')
    # 获得列数
    if cn:
        header = '%s.%s_multianno.header' % (file_out, ver)
        num = get_row_num(tmp_result, header)
        return result, num + 1
    else:
        return result


def db_format(gene_db='', region_db='', af_db='', filter_db='', dd_db=''):
    db_list = []
    ty_list = []
    if gene_db:
        for i in gene_db.strip().strip(',').split(','):
            db_list.append(i)
            ty_list.append('g')
    if dd_db:
        for i in dd_db.strip().strip(',').split(','):
            db_list.append(i)
            ty_list.append('f')
    if region_db:
        for i in region_db.strip().strip(',').split(','):
            db_list.append(i)
            ty_list.append('r')
    if af_db:
        for i in af_db.strip().strip(',').split(','):
            db_list.append(i)
            ty_list.append('f')
    if filter_db:
        for i in filter_db.strip().strip(',').split(','):
            db_list.append(i)
            ty_list.append('f')

    return ','.join(db_list), ','.join(ty_list)


def short_variants_filter(vcf, prefix, out_dir, gene_db, region_db, af_db, filter_db, dd_db, splice_distance, af_list,
                          af_th_list, retain_line, anno_dir, ref_version, script_path, thread):
    # vcf -> avinput
    avinput, info = short_variants_convert_format(vcf, prefix, out_dir, script_path)
    # af anno & filter
    db_for_af, ty_for_af = db_format(af_db=af_db, dd_db=dd_db)
    af_anno, num = anno_frequency(avinput, info, out_dir, prefix, db_for_af, ty_for_af, anno_dir, ref_version,
                                  script_path, thread)
    af_filted = f'{out_dir}/{prefix}_AF_filted.txt'
    af_list = af_list
    frequency_filter(af_anno, af_filted, af_list, af_th_list, retain_line)
    # exonic filter
    db_for_gene, ty_for_gene = db_format(gene_db=gene_db, dd_db=dd_db)
    af_info = f'{out_dir}/{prefix}_AF_info.txt'
    get_info = 'cut -f %d- %s > %s' % (num, af_filted, af_info)
    execute_system(get_info, '[ Msg: Get genotype done ! ]', '[ Error: Something wrong with get genotype ! ]')
    gene_anno, num = anno_db(af_filted, af_info, out_dir, prefix, '.gene_anno.txt', db_for_gene, ty_for_gene, anno_dir,
                             ref_version, script_path, thread, splice_distance=splice_distance)
    exonic_file = '%s/%s_exonic.txt' % (out_dir, prefix)
    exonic_filter(gene_anno, exonic_file, gene_db, retain_line)
    # annotation
    all_db, all_ty = db_format(gene_db=gene_db, region_db=region_db, dd_db=dd_db, af_db=af_db, filter_db=filter_db)
    ex_info = f'{out_dir}/{prefix}_exonic_info.txt'
    get_info = 'cut -f %d- %s > %s' % (num, exonic_file, ex_info)
    execute_system(get_info, '[ Msg: Get genotype done ! ]', '[ Error: Something wrong with get genotype ! ]')
    anno_file = anno_db(exonic_file, ex_info, out_dir, prefix, '.complete_anno.txt', all_db, all_ty, anno_dir,
                        ref_version, script_path, thread, splice_distance=splice_distance, cn=False)
    return anno_file


def anno_all_short_variants(vcf, prefix, out_dir,
                            gene_db, region_db, af_db, filter_db, dd_db,
                            splice_distance, anno_dir, ref_version, script_path, thread,
                            af_list, af_th, retain_line
                            ):
    avinput, info = short_variants_convert_format(vcf, prefix, out_dir, script_path)
    all_db, all_ty = db_format(gene_db=gene_db, region_db=region_db, dd_db=dd_db, af_db=af_db, filter_db=filter_db)
    _info = f'{out_dir}/{prefix}_info.txt'
    get_info = 'cut -f %d- %s > %s' % (9, avinput, _info)
    execute_system(get_info, '[ Msg: Get genotype done ! ]', '[ Error: Something wrong with get genotype ! ]')
    anno_file = anno_db(avinput, _info, out_dir, prefix, '.complete_anno.txt', all_db, all_ty, anno_dir, ref_version,
                        script_path, thread, splice_distance=splice_distance, cn=False, first=False)
    af_filted = f'{out_dir}/{prefix}_AF_filted.txt'
    af_list = af_list.strip(',').split(',')
    frequency_filter(anno_file, af_filted, af_list, af_th, retain_line)
    exonic_file = f'{out_dir}/{prefix}_exonic.txt'
    exonic_filter(af_filted, exonic_file, gene_db, retain_line)
    return exonic_file


def single_short_variants_ranking(vcf):
    pass


def trio_short_variants_ranking(vcf):
    pass


def single_structure_variants_filter(vcf):
    pass


def trio_structure_variants_filter(vcf):
    pass


def anno_sv():
    pass


def short_variants_anno():
    pass


def short_variants_step_anno():
    pass
