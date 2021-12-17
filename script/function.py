import os
import sys
import configparser
from script.fastq_qc import fastp_qc
from script.bam_qc import bam_stats
from script.mapping import bam_deal
from script.calling import gatk_pre, gatk, gatk_hard_filter, vardict, deep_variant_single, bcftools, strelka
from script.common import check_software, affinity, execute_system, gender_determined
from script.case_control import get_matrix, burden


def read_config(script_path, config_file):
    if config_file:
        if os.path.isfile(config_file):
            pass
        else:
            sys.exit('[ Error: Can not open <%s>.]' % config_file)
    else:
        config_file = script_path + '/lib/config.ini'
    config = configparser.ConfigParser()
    config.read(config_file)
    conf_dict = {}
    conf_dict['reference'] = config.get('global', 'reference')
    conf_dict['adapterR1'] = config.get('fastqc', 'adapter_r1')
    conf_dict['adapterR2'] = config.get('fastqc', 'adapter_r2')
    conf_dict['fastp_cmd'] = config.get('fastqc', 'fastp_cmd')
    conf_dict['gatk_bundle'] = config.get('database', 'GATK_bundle')
    conf_dict['mapping'] = config.get('parameter', 'mapping')
    conf_dict['platform'] = config.get('map', 'platform')
    conf_dict['caller'] = config.get('call', 'short_var')
    conf_dict['af_db'] = config.get('anno', 'af_db')
    conf_dict['gene_db'] = config.get('anno', 'gene_db')
    conf_dict['anno_dir'] = config.get('anno', 'anno_dir')
    conf_dict['ref_version'] = config.get('anno', 'ref_version')
    conf_dict['AFTh'] = eval(config.get('filter', 'af_threshold'))
    # case-control
    conf_dict['cc_AF'] = config.get('cc_filter', 'AF_list')
    conf_dict['cc_AF_AD'] = config.get('cc_filter', 'AF_th_AD')
    conf_dict['cc_AF_AR'] = config.get('cc_filter', 'AF_th_AR')
    conf_dict['cc_splice'] = config.get('cc_filter', 'splice_list')
    return conf_dict


def check_tmp(tmp_dir, out_dir):
    if tmp_dir:
        if not os.path.isdir(tmp_dir):
            os.makedirs(tmp_dir)
    else:
        tmp_dir = out_dir + '/tmp'
        os.makedirs(tmp_dir)
    tmp_dir = os.path.abspath(tmp_dir)
    return tmp_dir


def f2v(in_dir, out_dir, bed, prefix,
        vcf, fqc, bq,
        fmd, rm_dup, frmd,
        thread, script_path, config_file, tmp_dir, keep_tmp, gender_rate):
    rst = check_software('samtools')
    if rst:
        sys.exit('[ Error: Can not open <samtools>.]')
    # 读取参数
    config = read_config(script_path, config_file)
    # 输出目录
    in_dir = os.path.abspath(in_dir)
    out_dir = os.path.abspath(out_dir)
    # prefix
    if prefix:
        sample_name = prefix
    else:
        sample_name = os.path.split(os.path.abspath(in_dir))[-1]
    # temp directory
    tmp_dir = check_tmp(tmp_dir, out_dir)
    # QC
    print('[ Msg: QC running ... ]')
    fq_rep_dir = out_dir + '/reports/fastqQC'
    os.makedirs(fq_rep_dir)
    # how to mark duplication
    if frmd:
        fp_rmd = True
        sbb, rmd = False, False
    else:
        fp_rmd = False
        sbb = True if fmd else False
        rmd = True if rm_dup else False
    fq1, fq2 = fastp_qc(in_dir, tmp_dir, fq_rep_dir, config['adapterR1'], config['adapterR2'],
                        thread, script_path, sample_name, fp_rmd, fqc, config['fastp_cmd'])
    # mapping
    print('[ Msg: Waiting for mapping done ... ]')
    rst_out_dir = out_dir + '/results'
    os.makedirs(rst_out_dir)
    bam_rep_dir = out_dir + '/reports/bamQC'
    os.makedirs(bam_rep_dir)
    bam_file = bam_deal(fq1, fq2, rst_out_dir, bam_rep_dir, tmp_dir,
                        sample_name, config['reference'],
                        config['mapping'],
                        thread, script_path, bq, config['ref_version'],
                        bed, sbb, rmd, fp_rmd, config['platform'], keep_tmp, gender_rate)
    print('[ Msg: All sample mapping done ! ]')
    if not keep_tmp:
        rm_cmd = 'rm -rf %s %s' % (fq1, fq2)
        execute_system(rm_cmd, '[Msg: Delete <%s> process fastq done ... ]' % sample_name,
                       '[ Error: Something wrong with delete <%s> process fastq ! ]' % sample_name)
    print('[ Msg: Waiting for gatk calling done ... ]')
    # short variants calling
    if config['caller'] == 'gatk':
        gvcf = gatk_pre(bam_file, rst_out_dir, tmp_dir,
                        config['reference'], sample_name, config['gatk_bundle'], script_path, bed, keep_tmp)
        if vcf:
            vcf_rep_dir = out_dir + '/reports/vcfQC'
            os.makedirs(vcf_rep_dir)
            vcf_file = gatk_hard_filter(gvcf, rst_out_dir, vcf_rep_dir, tmp_dir, sample_name, config['reference'],
                                        script_path, bed, keep_tmp)
    elif config['caller'] == 'vardict':
        vcf_rep_dir = out_dir + '/reports/vcfQC'
        os.makedirs(vcf_rep_dir)
        vcf_file = vardict(bam_file, rst_out_dir, config['reference'], bed, vcf_rep_dir, tmp_dir, script_path,
                           sample_name, thread, filter_freq=0.01)
    elif config['caller'] == 'strelka2':
        vcf_rep_dir = out_dir + '/reports/vcfQC'
        os.makedirs(vcf_rep_dir)
        vcf_file = strelka(bam_file, rst_out_dir, vcf_rep_dir, config['reference'], script_path, thread, bed, tmp_dir,
                           sample_name)
    elif config['caller'] == 'deepvariants':
        vcf_rep_dir = out_dir + '/reports/vcfQC'
        os.makedirs(vcf_rep_dir)
        vcf_file = deep_variant_single(bam_file, rst_out_dir, config['reference'], bed, vcf_rep_dir, script_path,
                                       sample_name, thread, version="1.2.0")
    elif config['caller'] == 'bcftools':
        vcf_rep_dir = out_dir + '/reports/vcfQC'
        os.makedirs(vcf_rep_dir)
        vcf_file = bcftools(bam_file, rst_out_dir, tmp_dir, vcf_rep_dir, config['reference'], sample_name, thread,
                            script_path, bed)
    else:
        sys.exit('[Msg: Can not identify caller <%s>. ]' % config['caller'])


def bamQC(bam, bed, out_dir, tmp_dir, script_path, thread, bq, keep_tmp, gender_rate):
    tmp_dir = check_tmp(tmp_dir, out_dir)
    bam_stats(bam, out_dir, thread, tmp_dir, script_path, bq, '', bed, keep_tmp, gender_rate)


def trio_gt(p_gvcf, f_gvcf, m_gvcf, s_gvcfs, out_dir, script_path, config_file, keep_tmp, tmp_dir, bed, prefix):
    config = read_config(script_path, config_file)
    # 输入目录
    if s_gvcfs:
        gatk_gvcf = [p_gvcf, f_gvcf, m_gvcf] + s_gvcfs
    else:
        gatk_gvcf = [p_gvcf, f_gvcf, m_gvcf]
    # 输出目录
    out_dir = os.path.abspath(out_dir)
    tmp_dir = check_tmp(tmp_dir, out_dir)
    report_dir = out_dir + '/vcfQC'
    aff_dir = out_dir + '/affinity'
    os.makedirs(report_dir)
    os.makedirs(aff_dir)
    if prefix:
        prefix += '.'
    else:
        prefix = ''
    vcf = gatk(gatk_gvcf, out_dir, report_dir,
               config['reference'], config['gatk_bundle'], script_path, prefix, bed,
               tmp_dir, keep_tmp)
    affinity(vcf, aff_dir, script_path)


def single_gt(gvcf, out_dir, script_path, bed, tmp_dir, keep_tmp, prefix, config_file, caller='gatk_hard', noqc=False,
              noflt=False):
    flt = not noflt
    if not tmp_dir:
        tmp_dir = '%s/tmp_dir' % out_dir
    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir)
    config = read_config(script_path, config_file)
    out_dir = os.path.abspath(out_dir)  # + '/result'
    report_dir = out_dir + '/vcfQC'
    os.makedirs(report_dir)
    if caller == 'gatk_hard':
        vcf = gatk_hard_filter(gvcf, out_dir, report_dir, tmp_dir, prefix, config['reference'], script_path, bed,
                               keep_tmp)
    elif caller == 'vqsr':
        vcf = gatk([gvcf], out_dir, report_dir, config['reference'], config['gatk_bundle'], script_path, prefix,
                   bed, tmp_dir, keep_tmp, noqc, flt)
    else:
        sys.exit('[Error: Can not identify caller <%s>. ]' % caller)


def burden_test(case, control, case_matrix, control_matrix,
                out_dir, fp, mode, cutoff, method, gene, score, script_path, config_file):
    config = read_config(script_path, config_file)
    # path
    out_dir = os.path.abspath(out_dir)
    out_case = out_dir + '/case_out'
    out_control = out_dir + '/control_out'
    os.makedirs(out_case)
    os.makedirs(out_control)

    if method:
        # case
        if case and not case_matrix:
            case_matrix = get_matrix(case, 'CASE_', config, out_case, mode, snvdb=fp)
            case_matrix.to_csv(out_case + '/case_matrix.txt', index=False, sep='\t')
        # control
        if control and not control_matrix:
            control_matrix = get_matrix(control, 'CONTROL_', config, out_control, mode, snvdb=fp)
            control_matrix.to_csv(out_control + '/control_matrix.txt', index=False, sep='\t')
    else:
        # case
        if case and not case_matrix:
            case_matrix = get_matrix(case, 'CASE_', config, mode=mode, gene_col_name=gene, score_col_name=score,
                                     snvdb=fp, pattern=False)
            case_matrix.to_csv(out_case + '/case_matrix.txt', index=False, sep='\t')
        # control
        if control and not control_matrix:
            control_matrix = get_matrix(control, 'CONTROL_', config, mode=mode, gene_col_name=gene,
                                        score_col_name=score, snvdb=fp, pattern=False)
            control_matrix.to_csv(out_control + '/control_matrix.txt', index=False, sep='\t')
    rank_df = burden(case_matrix, control_matrix, cutoff=cutoff)
    rank_df.to_csv(out_dir + '/gene_rank.txt', index=False, sep='\t')


def gender_identify(bam, bed, out_dir, thread, script_path, rate):
    if os.path.isdir(out_dir):
        rm_tmp = False
    else:
        os.makedirs(out_dir)
        rm_tmp = True
    bed_cmd = 'grep \'[XY]\' %s > %s' % (bed, out_dir + '/xy.bed')
    bed = out_dir + '/xy.bed'
    os.system(bed_cmd)
    x, y = 0, 0
    with open(bed) as f:
        for i in f:
            if i.startswith('X') or i.startswith('chrX'):
                x += 1
            elif i.startswith('Y') or i.startswith('chrY'):
                y += 1
    if not x or not y:
        if rm_tmp:
            os.system('rm -rf %s' % out_dir)
        else:
            os.system('rm -rf %s' % bed)
        sys.exit('[ Error: XY chromosome coverage incomplete.]')
    depth_file_prefix = out_dir + '/' + bam.split('/')[-1].rstrip('.bam')
    depth_cmd = script_path + '/bin/mosdepth -n -t %d -b %s -T 1,10,20,30 %s %s' % (
        thread, bed, depth_file_prefix, bam)
    execute_system(depth_cmd, '[ Msg: bam stat depth done ! ]', '[ Error: Something wrong with bam stat depth ! ]')
    gender = gender_determined(depth_file_prefix + '.thresholds.bed.gz', rate)
    if rm_tmp:
        os.system('rm -rf %s' % out_dir)
    else:
        os.system('rm -rf %s* %s' % (depth_file_prefix, bed))


def trio_analysis():
    pass
    # annoDir = outDir + '/annotation'
    # os.makedirs(annoDir)
    # trio_short_variants_filter(vcf, annoDir, scriptPath,
    #                            config['af_db'], config['af_type'], config['gene_db'], config['gene_type'],
    #                            config['AFList'], config['ClinList'],
    #                            config['AFTh'],
    #                            config['anno_dir'], config['ref_version'], thread)
    # print('[ Msg: All sample gatk calling done ! ]')


def variants_call(bam, out_dir, caller, bed, prefix, thread, tmp_dir, keep_tmp, config_file, script_path, noqc, noflt):
    flt = not noflt
    config = read_config(script_path, config_file)
    out_dir = os.path.abspath(out_dir)
    if not tmp_dir:
        tmp_dir = '%s/tmp_dir' % out_dir
    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir)
    vcf_rep_dir = out_dir + '/reports/vcfQC'
    os.makedirs(vcf_rep_dir)
    if caller == 'gatk':
        gvcf = gatk_pre(bam, out_dir, tmp_dir,
                        config['reference'], prefix, config['gatk_bundle'], script_path, bed, keep_tmp)
        vcf_file = gatk([gvcf], out_dir, vcf_rep_dir, config['reference'], config['gatk_bundle'], script_path, prefix,
                        bed, tmp_dir, keep_tmp, noqc, flt)
    elif caller == 'gatk_hard':
        gvcf = gatk_pre(bam, out_dir, tmp_dir,
                        config['reference'], prefix, config['gatk_bundle'], script_path, bed, keep_tmp)
        vcf_file = gatk_hard_filter(gvcf, out_dir, vcf_rep_dir, tmp_dir, prefix, config['reference'], script_path, bed,
                                    keep_tmp, noqc, flt)
    elif caller == 'vardict':
        vcf_file = vardict(bam, out_dir, config['reference'], bed, vcf_rep_dir, tmp_dir, script_path, prefix, thread,
                           0.01, noqc, flt)
    elif caller == 'strelka2':
        vcf_file = strelka(bam, out_dir, vcf_rep_dir, config['reference'], script_path, thread, bed, tmp_dir, prefix,
                           noqc, flt)
    elif caller == 'deepvariant':
        vcf_file = deep_variant_single(bam, out_dir, config['reference'], bed, vcf_rep_dir, script_path, prefix, thread,
                                       noqc, flt, version="1.2.0")
    elif caller == 'bcftools':
        vcf_file = bcftools(bam, out_dir, tmp_dir, vcf_rep_dir, config['reference'], prefix, thread, script_path, bed,
                            noqc, flt)
    else:
        sys.exit('[Error: Can not identify caller <%s>. ]' % caller)
