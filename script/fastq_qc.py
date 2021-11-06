import os
from script.common import execute_system, fastq_prework


def fastp_qc(in_dir, out_dir, report_dir, adapter1, adapter2, thread, script_path, prefix, rmd, fqc, fp_cmd_add):
    fq1_list, fq2_list = fastq_prework(in_dir)
    if len(fq1_list) > 1:
        fq1_file, fq2_file = '', ''
        for fq1, fq2 in zip(fq1_list, fq2_list):
            fq1_file += '%s/%s ' % (in_dir, fq1)
            fq2_file += '%s/%s ' % (in_dir, fq2)
        fq1 = out_dir + '/merge_R1.fq.gz'
        fq2 = out_dir + '/merge_R2.fq.gz'
        merge_cmd = 'cat %s > %s && cat %s > %s ' % (fq1_file, fq1, fq2_file, fq2)
        execute_system(merge_cmd, '[ Msg: Fastq merge done !]', '[ Error: something wrong with merge fastq ！]')
    else:
        fq1 = in_dir + '/' + fq1_list[0]
        fq2 = in_dir + '/' + fq2_list[0]
    print('[ Msg: Waiting for fastq qc ... ]')
    out_fq1 = out_dir + '/merge_clean_R1.fq.gz'
    out_fq2 = out_dir + '/merge_clean_R2.fq.gz'

    html = report_dir + '/' + prefix + '.html'
    json_file = report_dir + '/' + prefix + '.json'
    # cmd
    fp_cmd = script_path + '/bin/fastp -i %s -I %s -o %s -O %s -h %s -j %s -w %d ' % (
        fq1, fq2, out_fq1, out_fq2, html, json_file, thread)
    fp_cmd += fp_cmd_add
    if rmd:
        fp_cmd += ' --dedup '
    if adapter1 and adapter2:
        fp_cmd += ' --adapter_sequence %s --adapter_sequence_r2 %s ' % (adapter1, adapter2)

    execute_system(fp_cmd, '[ Msg: Fastq quality control done with fastp !]',
                   '[ Error: Run fastq quality control fail with fastp ！]')
    # 使用fastqc进行质控
    if fqc:
        report_dir_raw_fq = report_dir + '/fastqc_raw_fastq'
        report_dir_clean_fq = report_dir + '/fastqc_clean_fastq'
        os.makedirs(report_dir_raw_fq)
        os.makedirs(report_dir_clean_fq)
        raw_fastqc_cmd = script_path + '/bin/FastQC/fastqc -t 2 -o %s %s %s' % (report_dir_raw_fq, fq1, fq2)
        execute_system(raw_fastqc_cmd, '[ Msg: %s Raw fastq QC done ! ]' % prefix,
                       '[ Error: Something wrong with raw data fastqc in QC ! ]')
        clean_fastqc_cmd = script_path + '/bin/FastQC/fastqc -t 2 -o %s %s %s' % (report_dir_clean_fq, out_fq1, out_fq2)
        execute_system(clean_fastqc_cmd, '[ Msg: %s Clean fastq QC done ! ]' % prefix,
                       '[ Error: Something wrong with clean data fastqc in QC ! ]')
    return out_fq1, out_fq2
