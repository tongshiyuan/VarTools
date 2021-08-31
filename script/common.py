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
        sys.exit('[ E: Please use it by python3.x ! ]')
    # 确保这些一定会使用到的软件存在
    software = ['samtools', 'plot-bamstats', 'qualimap', 'bcftools', 'R', 'bgzip']
    for _ in software:
        check_software(_)


def fastq_prework(indir, max_process):
    # indir = os.path.abspath(indir)
    # outdir = os.path.abspath(outdir)
    # report_dir = os.path.abspath(report_dir)
    sample_name = []
    # for file in os.listdir(indir):
    #     if re.findall(r'(1\.fq\.gz|1\.fastq\.gz|1\.fq|1\.fastq)$', file):
    #         sample_name.append(file)
    # if max_process >= len(sample_name):
    #     max_process = len(sample_name)
    tmp_file_name = []
    for file in os.listdir(indir):
        if re.findall(r'(\.fq\.gz|\.fastq\.gz|\.fq|\.fastq)$', file):
            tmp_file_name.append(file)
    # check name
    _ex = tmp_file_name[0]
    name_info = _ex.split('_')
    # for info in name_info:
    if name_info[-1][0] in ['1', '2'] and name_info[-2].startswith('L'):
        for file in tmp_file_name:
            if re.findall(r'(1\.fq\.gz|1\.fastq\.gz|1\.fq|1\.fastq)$', file):
                sample_name.append(file)
        name_state = 'novo'
    elif name_info[-2] in ['R1', 'R2'] and name_info[-3].startswith('L'):
        for file in tmp_file_name:
            if file.split('_')[-2] == 'R1':
                sample_name.append(file)
        name_state = 'giab'
    else:
        sys.exit('[ E: Can not analysis input raw data name <%s>]' % _ex)
    if max_process >= len(sample_name):
        max_process = len(sample_name)
    return indir, max_process, sample_name, name_state


# def get_raw_info(name_in, state):
#     # W2018005_NZTD180700064_H5MYLDSXX_L3_1.fq.gz
#     # NA24695_CTTGTA_L002_R2_014.fastq.gz
#     info = name_in.split('_')
#     if state == 'novo':
#         return info[0], info[0], info[3]
#     elif state == 'giab':
#         return info[0], info[0], info[2]
#     else:
#         return info[0], info[0], info[0]

def get_raw_info(fq1, state):
    # W2018005_NZTD180700064_H5MYLDSXX_L3_1.fq.gz
    # NA24695_CTTGTA_L002_R2_014.fastq.gz
    name_in = fq1.split('/')[-1]
    info = name_in.split('_')
    if state is 'novo':
        return info[0], info[2], info[3]
    else:
        with gzip.open(fq1, 'rb') as f:
            header = f.readline().decode()
            _lb = header.split(':')[2]
            _id = header.split(':')[3]
        return info[0], _lb, _id


def get_row_num(tmp_result, header):
    num_cmd = 'head -n 1 %s > %s ' % (tmp_result, header)
    os.system(num_cmd)
    with open(header) as f:
        num = len(f.readline().strip().split('\t'))
    return num + 1


def affinity(vcf, affDir, scriptPath):
    plinkCmd = scriptPath + '/bin/plink --double-id --vcf %s --make-bed --out %s --allow-extra-chr' \
               % (vcf, affDir + '/aff')
    rst = os.system(plinkCmd)
    if rst:
        print('[ W: fail to make bed！]')
        return 0
    else:
        print('[ Msg: make bed done! ]')
        kingCmd = scriptPath + '/bin/king -b %s --kinship --prefix %s ' % (affDir + '/aff.bed', affDir + '/aff')
        rst = os.system(kingCmd)
        if rst:
            print('[ W: fail to run kingship ！]')
        else:
            print('[ Msg: calculate affinity done! ]')
            rm = open(affDir + '/readme.txt', 'w')
            rm.write('说明：当 Kinship > 0.354，可能为同一样本或者孪生兄弟姐妹；\n')
            rm.write('     当 Kinship 在[0.177, 0.354]内，为一级亲缘；\n')
            rm.write('     当 Kinship 在[0.0884, 0.177]内，为二级亲缘；\n')
            rm.write('     当 Kinship 在[0.0442, 0.0884]内，为三级亲缘；\n')
            rm.write('note: ENGLISH note please read KING\'s document, https://www.kingrelatedness.com\n')
            rm.close()


def gender_determined(vcf):
    pass