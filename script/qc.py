import os
from multiprocessing import Pool
from script.common import execute_system, fastq_prework


# 输入是包含一对（或多对）fastq文件的目录，文件名须以1/2.fq.gz或1/2.fastq.gz或1/2.fq或1/2.fastq结尾，同对文件前缀一致
# 输出为clean以后的文件目录（新建，已存在则报错退出），及其报告（单独的目录）
def fastp_qc(indir,
             outdir,
             report_dir,
             adapter1,
             adapter2,
             fastp_thread,
             script_path,
             max_process):
    indir, max_process, sample_name, state = fastq_prework(indir, max_process)
    if max_process:
        process = Pool(max_process)
        for file in sample_name:
            process.apply_async(single_fastq_qc, args=(
                indir, outdir, file, report_dir, adapter1, adapter2, fastp_thread, state, script_path,))
        print('[ Msg: Waiting for all sample qc done ... ]')
        process.close()
        process.join()
    else:
        for file in sample_name:
            single_fastq_qc(indir, outdir, file, report_dir, adapter1, adapter2, fastp_thread, state, script_path)
    print('[ Msg: All sample qc done ! ]')


def single_fastq_qc(indir, outdir, file, report_dir, adapter1, adapter2, fastp_thread, name_state, scriptPath):
    # 使用fastp进行质控
    # W2018005_NZTD180700064_H5MYLDSXX_L3_1.fq.gz
    # NA24695_CTTGTA_L002_R1_014.fastq.gz
    fastq1 = indir + '/' + file
    if name_state == 'novo':
        fastq2 = indir + '/' + file[:file.rfind('1')] + '2' + file[file.rfind('1') + 1:]
        out_fastq1 = outdir + '/' + file
        out_fastq2 = outdir + '/' + file[:file.rfind('1')] + '2' + file[file.rfind('1') + 1:]
        html = report_dir + '/' + file[:file.rfind('1')] + '.html'
        json_file = report_dir + '/' + file[:file.rfind('1')] + '.json'
    else:
        fastq2 = indir + '/' + file[:file.rfind('R1')] + 'R2' + file[file.rfind('R1') + 2:]
        out_fastq1 = outdir + '/' + file
        out_fastq2 = outdir + '/' + file[:file.rfind('R1')] + 'R2' + file[file.rfind('R1') + 2:]
        html = report_dir + '/' + file[:file.find('.')] + '.html'
        json_file = report_dir + '/' + file[:file.find('.')] + '.json'

    if not adapter1 or not adapter2:
        command = scriptPath + '/bin/fastp -i %s -I %s -o %s -O %s -h %s -j %s -c -q 20 -u 50 -n 15 -5 20 -3 20 -w %d' \
                  % (fastq1, fastq2, out_fastq1, out_fastq2, html, json_file, fastp_thread)
    else:
        command = scriptPath + '/bin/fastp -i %s -I %s -o %s -O %s -h %s -j %s --adapter_sequence %s ' \
                               '--adapter_sequence_r2 %s -c -q 20 -u 50 -n 15 -5 20 -3 20 -w %d' \
                  % (fastq1, fastq2, out_fastq1, out_fastq2, html, json_file, adapter1, adapter2, fastp_thread)
    execute_system(command, '[ Msg: Fastq quality control done with fastp !]',
                   '[ E: Run fastq quality control fail with fastp ！]')
    # 使用fastqc进行质控
    report_dir_raw_fq = report_dir + '/raw_fastq'
    report_dir_clean_fq = report_dir + '/clean_fastq'
    os.makedirs(report_dir_raw_fq)
    os.makedirs(report_dir_clean_fq)
    raw_fastqc_cmd = scriptPath + '/bin/FastQC/fastqc -o %s %s %s' % (report_dir_raw_fq, fastq1, fastq2)
    execute_system(raw_fastqc_cmd, '[ Msg: %s Raw fastq QC done ! ]' % file,
                   '[ E: Something wrong with raw data fastqc in QC ! ]')
    clean_fastqc_cmd = scriptPath + '/bin/FastQC/fastqc -o %s %s %s' % (report_dir_clean_fq, out_fastq1, out_fastq2)
    execute_system(clean_fastqc_cmd, '[ Msg: %s Clean fastq QC done ! ]' % file,
                   '[ E: Something wrong with clean data fastqc in QC ! ]')
