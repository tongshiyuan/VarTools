import os
import pandas as pd
from script.common import check_software, gender_determined


def bam_stats(duped_bam, report_dir, thread, tmp_dir, script_path, bq, ver, bed, keep_tmp, gender_rate):
    if bq:
        rst = check_software('qualimap')
        if rst:
            print('[ Warn: Can not find qualimap. Skipping the bam quality check with qualimap ! ]')
        else:
            bamqc_cmd = 'qualimap bamqc --java-mem-size=20G -bam %s -c -nt %d -outdir %s -outformat PDF:HTML' % (
                duped_bam, thread, report_dir)
            rst = os.system(bamqc_cmd)
            if rst:
                print('[ Error: fail to bam QC with qualimap! ]')
            else:
                print('[ Msg: qualimap bam QC done ! ]')
    if not bed and ver == 'hg38':
        bed_file = script_path + '/lib/Hg38.genome.bed'
    elif not bed and ver == 'hg19':
        bed_file = script_path + '/lib/Hg19.genome.bed'
    else:
        bed_file = bed
    bam_qc(duped_bam, report_dir, tmp_dir, script_path, thread, bed_file, keep_tmp, gender_rate)


def bam_qc(bam, report_dir, tmp_dir, script_path, thread, bed_file, keep_tmp, gender_rate):
    base_cov, _start, end = 0, 0, 0
    _chr = ['chr1', '1']
    # 目标区大小
    new_target_file = open(tmp_dir + '/tmp_target.txt', 'w')
    with open(bed_file) as f:
        for line in f:
            record = line.strip().split('\t')
            if eval(record[1]) > eval(record[2]):
                print('[ Error: The positions are not in chromosomal order (%s:%s comes after %s) ]' % (
                    record[0], record[2], record[1]))
                new_target_file.close()
                return 0
            new_target_file.write('\t'.join([record[0], str(eval(record[1]) + 1), record[2]]) + '\n')
            if record[0] in _chr and eval(record[1]) <= end:
                if _start > eval(record[1]):
                    print('[ Error: The positions are not in chromosomal order (%s:%s comes before %s) ]'
                          % (record[0], _start, record[1]))
                    new_target_file.close()
                    return 0
                base_cov = base_cov - (end - eval(record[1]))
            base_cov += eval(record[2]) - eval(record[1])
            # 更新
            _chr, _start, end = record[0], eval(record[1]), eval(record[2])
    new_target_file.close()
    # 比对率、平均测序深度
    new_target_file = tmp_dir + '/tmp_target.txt'
    single_bam_qc(bam, bed_file, base_cov, tmp_dir, report_dir, new_target_file, script_path, thread, keep_tmp,
                  gender_rate)
    if not keep_tmp:
        os.system('rm -rf %s' % new_target_file)
    print('[ Msg: Bam QC done！]')


def single_bam_qc(bam, bed_file, base_cov, tmp_dir, report_dir, new_target_file, script_path, thread, keep_tmp,
                  gender_rate):
    result_dict = {
        'target_dep_num': 0,
        'dep_1x': 0,
        'dep_10x': 0,
        'dep_20x': 0,
        'dep_30x': 0,
        'mapping_rate': 0,
        'average_coverage': 0,
        'gender': ''
    }
    # flagstat
    flag_cmd = 'samtools flagstat -@ %d %s > %s' % (thread, bam, report_dir + '/flagstat.txt')
    rst = os.system(flag_cmd)
    if rst:
        print('[ Error: Something wrong with bam flagstat ! ]')
    else:
        print('[ Msg: bam flagstat done ! ]')
        with open(report_dir + '/flagstat.txt') as f:
            for line in f:
                if line.find('mapped (') != -1:
                    mapping_rate = line.split('(')[-1].split(':')[0].strip()
                    result_dict['mapping_rate'] = mapping_rate
                    break
    # stats
    stats_cmd = 'samtools stats -d -@ %d -t %s %s > %s' % (thread, new_target_file, bam, report_dir + '/stats.txt')
    rst = os.system(stats_cmd)
    if rst:
        print('[ Error: Something wrong with bam stats ! ]')
    else:
        print('[ Msg: bam stats done ! ]')
        with open(report_dir + '/stats.txt') as f:
            for i in f:
                if i.startswith('SN') and i.find('bases mapped:') != -1:
                    result_dict['target_dep'] = eval(i.split(':')[-1].split('#')[0].strip())
    tmp_dir = tmp_dir + '/mosdepth/'
    os.makedirs(tmp_dir)
    # coverage
    depth_cmd = script_path + '/bin/mosdepth -n -t %d -b %s -T 1,10,20,30 %s %s' % (
        thread, bed_file, tmp_dir + '/' + bam.split('/')[-1].rstrip('.bam'), bam)
    rst = os.system(depth_cmd)
    if rst:
        print('[ Error: Something wrong with bam stat depth ! ]')
    else:
        print('[ Msg: bam stat depth done ! ]')
        gender = gender_determined(tmp_dir + '/' + bam.split('/')[-1].rstrip('.bam') + '.thresholds.bed.gz',
                                   gender_rate)
        result_dict['gender'] = gender
        data = pd.read_table(tmp_dir + '/' + bam.split('/')[-1].rstrip('.bam') + '.thresholds.bed.gz',
                             low_memory=False,
                             compression='gzip')
        sum_end_start = sum(data['end'] - data['start'])
        sum_1X = sum(data['1X'])
        sum_10X = sum(data['10X'])
        sum_20X = sum(data['20X'])
        sum_30X = sum(data['30X'])
        p1 = sum_1X / sum_end_start
        p10 = sum_10X / sum_end_start
        p20 = sum_20X / sum_end_start
        p30 = sum_30X / sum_end_start
        result_dict['dep_1x'] = p1
        result_dict['dep_10x'] = p10
        result_dict['dep_20x'] = p20
        result_dict['dep_30x'] = p30
        df = pd.read_table(tmp_dir + '/' + bam.split('/')[-1].rstrip('.bam') + '.mosdepth.summary.txt')
        df = df[df['chrom'].isin(
            ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
             'chr14', 'chr15', 'chr16', 'chr17', 'chr88', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY'])]
        result_dict['target_dep_num'] = sum(df['bases'])
        result_dict['average_coverage'] = sum(df['bases']) / base_cov
    # 区域总base 太费时了
    # bed_cov_list = []
    # bedcov_cmd = 'samtools bedcov %s %s > %s' % (bedfile, bam, report_dir + '/bedcov.txt')
    # rst = os.system(bedcov_cmd)
    # if rst:
    #     print('[ Error: Something wrong with bam bedcov! ]')
    # else:
    #     print('[ Msg: bam bedcov done ! ]')
    #     with open(report_dir + '/bedcov.txt') as f:
    #         for i in f:
    #             bed_cov_list.append(eval(i.split('\t')[-1]))
    #         result_dict['target_dep_num'] = sum(bed_cov_list)
    # total_dep_num 太耗时
    # stats_all_cmd = 'samtools stats -d -@ %d %s > %s' % (thread, bam, tmp_dir + '/stats_arst.tmp')
    # rst = os.system(stats_all_cmd)
    # if rst:
    #     print('[ Error: Something wrong with bam stats all ! ]')
    # else:
    #     print('[ Msg: bam stats all done ! ]')
    #     with open(tmp_dir + 'stats_arst.tmp') as f:
    #         for i in f:
    #             if i.startswith('SN') and i.find('bases mapped:') > 0:
    #                 result_dict['total_dep_num'] = eval(i.split(':')[-1].split('#')[0].strip())
    fo = open(report_dir + '/Bam_QC.txt', 'w')
    for k, v in result_dict.items():
        fo.write(k + '\t' + str(v) + '\n')
    fo.close()
    if not keep_tmp:
        os.system('rm -rf %s' % tmp_dir)


def ref_n_counts(reference, script_path):
    # 计算基因组每条染色体N的个数, 输出是统计覆盖率要用的，存在于当前目录下的lib目录中
    N_count = {}
    _chr = ''
    with open(reference) as f:
        for line in f:
            if line.startswith('>'):
                _chr = line.lstrip('>').rstrip().split('\t')[0].split(' ')[0]
            else:
                N_count[_chr] = line.count('N') + N_count.get(_chr, 0)
    n_file = open(script_path + '/lib/' + reference + '.N.txt', 'w')
    for k, v in N_count.items():
        n_file.write(k + '\t' + str(v) + '\n')
    n_file.close()
    print('[ Msg: N number file created done !]')


def n_region(fasta):
    from Bio import SeqIO
    bed_list = []
    # import the SeqIO module from Biopython
    with open(fasta) as fasta_handle:
        for record in SeqIO.parse(fasta_handle, "fasta"):
            start_pos, counter, gap, gap_length = 0, 0, False, 0
            for char in record.seq:
                if char in ['N', 'n']:
                    if gap_length == 0:
                        start_pos = counter
                        gap_length = 1
                        gap = True
                    else:
                        gap_length += 1
                else:
                    if gap:
                        bed_list.append(record.id + "\t" + str(start_pos) + "\t" + str(start_pos + gap_length))
                        gap_length = 0
                        gap = False
                counter += 1
    return bed_list
