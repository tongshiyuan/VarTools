import os
import sys
import platform
import re
import gzip


def execute_system(cmd, done_notice='', fail_notice=''):
    # 函数用于执行命令行，输入命令行语句以及警告信息，一旦运行失败，输出命令以及警告信息，然后退出程序
    result = os.system(cmd)
    if result:
        print('\n[ Error cmd: < %s > ]\n' % cmd)
        sys.exit(fail_notice)
    else:
        print(done_notice)


def check_software(soft):
    # 检查软件是否存在，不存在就退出
    cmd = 'command -v %s >/dev/null 2>&1 || { echo >&2 "[ E: VarTools require %s install it and add to system path. ' \
          'Aborting. ] "; exit 1; }' % (soft, soft)
    # type foo >/dev/null 2>&1 || { echo >&2 "I require foo but it's not installed.  Aborting."; exit 1; }
    # hash foo 2>/dev/null || { echo >&2 "I require foo but it's not installed.  Aborting."; exit 1; }
    # execute_system(cmd)
    result = os.system(cmd)
    return result


def check_system():
    # 查看系统python
    if platform.python_version() < '3.0.0':
        sys.exit('[ Error: Please use it by python3.x ! ]')
    # 确保这些一定会使用到的软件存在
    software = ['samtools', 'plot-bamstats', 'qualimap', 'bcftools', 'R', 'bgzip']
    for _ in software:
        check_software(_)


def paired_fq(fq1, kw1, kw2):
    if fq1.count(kw1) == 1:
        fq2 = fq1.replace(kw1, kw2)
        sfn = fq1.replace(kw1, '_')
        sfn = re.sub(r'(\.fq\.gz|\.fastq\.gz|\.fq|\.fastq|_fq\.gz|_fq|_fastq|_fastq\.gz)$', '', sfn)
    else:
        rec = fq1.split(kw1)
        fq2 = kw1.join(rec[:-1]) + kw2 + rec[-1]
        sfn = kw1.join(rec[:-1]) + '_' + rec[-1]
        sfn = re.sub(r'(\.fq\.gz|\.fastq\.gz|\.fq|\.fastq|_fq\.gz|_fq|_fastq|_fastq\.gz)$', '', sfn)
    return fq2, sfn


def fastq_prework(in_dir):
    tmp_file_name, fq1_list, fq2_list = [], [], []
    for file in os.listdir(in_dir):
        if re.findall(r'(\.fq\.gz|\.fastq\.gz|\.fq|\.fastq|_fq\.gz|_fq|_fastq|_fastq\.gz)$', file):
            tmp_file_name.append(file)
    if len(tmp_file_name) < 2 or len(tmp_file_name) % 2 == 1:
        sys.exit('[ Error: Can not identify paired fastq data. Please check! ]')
    # check name: xxx_1.fq.gz, xxx.1.fq.gz, xxx_1_xxx.fq.gz, xxx_R1_xxx.fq.gz, xxx.R1.xxx.fq.gz
    # W2018005_NZTD180700064_H5MYLDSXX_L3_1.fq.gz
    # NA24695_CTTGTA_L002_R1_014.fastq.gz
    for fq in tmp_file_name:
        if fq.find('.R1.') != -1:
            fq1_list.append(fq)
            fq2, sfn = paired_fq(fq, '.R1.', '.R2.')
            if fq2 in tmp_file_name:
                fq2_list.append(fq2)
            else:
                sys.exit('[ Error: Can not find paired file of <%s>]' % fq)
        elif fq.find('_R1_') != -1:
            fq1_list.append(fq)
            fq2, sfn = paired_fq(fq, '_R1_', '_R2_')
            if fq2 in tmp_file_name:
                fq2_list.append(fq2)
            else:
                sys.exit('[ Error: Can not find paired file of <%s>]' % fq)
        elif fq.find('_R1.') != -1:
            fq1_list.append(fq)
            fq2, sfn = paired_fq(fq, '_R1.', '_R2.')
            if fq2 in tmp_file_name:
                fq2_list.append(fq2)
            else:
                sys.exit('[ Error: Can not find paired file of <%s>]' % fq)
        elif fq.find('.1.') != -1:
            fq1_list.append(fq)
            fq2, sfn = paired_fq(fq, '.1.', '.2.')
            if fq2 in tmp_file_name:
                fq2_list.append(fq2)
            else:
                sys.exit('[ Error: Can not find paired file of <%s>]' % fq)
        elif fq.find('_1_') != -1:
            fq1_list.append(fq)
            fq2, sfn = paired_fq(fq, '_1_', '_2_')
            if fq2 in tmp_file_name:
                fq2_list.append(fq2)
            else:
                sys.exit('[ Error: Can not find paired file of <%s>]' % fq)
        elif fq.find('_1.') != -1:
            fq1_list.append(fq)
            fq2, sfn = paired_fq(fq, '_1.', '_2.')
            if fq2 in tmp_file_name:
                fq2_list.append(fq2)
            else:
                sys.exit('[ Error: Can not find paired file of <%s>]' % fq)
        # else:
        #     sys.exit('[ Error: Can not identify fastq file structure of <%s>]' % fq)
    if len(tmp_file_name) != len(fq1_list) + len(fq2_list):
        sys.exit('[ Error: Some fastq file structure can not be identified. ]')
    if len(fq1_list) != len(fq2_list):
        sys.exit('[ Error: Something wrong with identify fastq data. ]')
    return fq1_list, fq2_list


def get_raw_info(fq1):
    if fq1.endswith('.gz'):
        with gzip.open(fq1, 'rb') as f:
            header = f.readline().decode()
            _lb = header.split(':')[2]
            _id = header.split(':')[3]
    else:
        with open(fq1) as f:
            header = f.readline()
            _lb = header.split(':')[2]
            _id = header.split(':')[3]
    return _lb, _id


def get_row_num(tmp_result, header):
    num_cmd = 'head -n 1 %s > %s ' % (tmp_result, header)
    os.system(num_cmd)
    with open(header) as f:
        num = len(f.readline().strip().split('\t'))
    return num + 1


def affinity(vcf, aff_dir, script_path):
    plink_cmd = script_path + '/bin/plink --double-id --vcf %s --make-bed --out %s --allow-extra-chr' \
                % (vcf, aff_dir + '/aff')
    rst = os.system(plink_cmd)
    if rst:
        print('[ Warn: fail to make bed！]')
        return 0
    else:
        print('[ Msg: make bed done! ]')
        king_cmd = script_path + '/bin/king -b %s --kinship --prefix %s ' % (aff_dir + '/aff.bed', aff_dir + '/aff')
        rst = os.system(king_cmd)
        if rst:
            print('[ Warn: fail to run kingship ！]')
        else:
            print('[ Msg: calculate affinity done! ]')
            rm = open(aff_dir + '/readme.txt', 'w')
            rm.write('说明：当 Kinship > 0.354，可能为同一样本或者孪生兄弟姐妹；\n')
            rm.write('     当 Kinship 在[0.177, 0.354]内，为一级亲缘；\n')
            rm.write('     当 Kinship 在[0.0884, 0.177]内，为二级亲缘；\n')
            rm.write('     当 Kinship 在[0.0442, 0.0884]内，为三级亲缘；\n')
            rm.write('note: ENGLISH note please read KING\'s document, https://www.kingrelatedness.com\n')
            rm.close()


def gender_determined(depth_file, rate=20):
    x_region, x_bases, y_region, y_bases = 0, 0, 0, 0
    if depth_file.endswith('gz'):
        with gzip.open(depth_file, 'rb') as f:
            for line in f:
                de_line = line.decode().split('\t')
                if de_line[0].startswith('X') or de_line[0].startswith('chrX'):
                    x_region += eval(de_line[2]) - eval(de_line[1])
                    x_bases += eval(de_line[4])
                elif de_line[0].startswith('Y') or de_line[0].startswith('chrY'):
                    y_region += eval(de_line[2]) - eval(de_line[1])
                    y_bases += eval(de_line[4])
    else:
        with open(depth_file) as f:
            for line in f:
                de_line = line.split('\t')
                if de_line[0].startswith('X') or de_line[0].startswith('chrX'):
                    x_region += eval(de_line[2]) - eval(de_line[1])
                    x_bases += eval(de_line[4])
                elif de_line[0].startswith('Y') or de_line[0].startswith('chrY'):
                    y_region += eval(de_line[2]) - eval(de_line[1])
                    y_bases += eval(de_line[4])
    if not x_region or not y_region:
        print('[ Error: XY chromosome coverage incomplete.]')
        return ''
    else:
        if (x_bases / x_region) / (y_bases / y_region) > rate:
            print('[ Msg: This sample may FEMALE.]')
            return 'female'
        else:
            print('[ Msg: This sample may MALE.]')
            return 'male'
