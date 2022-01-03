#!/usr/bin/python3
# -*- coding:utf-8 -*-

##################################################
# Author: Shi-Yuan Tong
# Email: tongshiyuan@foxmail.com
# Created Time: 2021-5-14
# log:
# 21.05.31: 根据新的框架修改，为了拓展性、适配性更好
# 21.08.31: 完善框架，修改了一些已知错误
# 21.09.03: 修改传入参数，添加GRIPT版case-control (未根据源码而是根据论文复现，会存在与原工具的结果差异)
# 21.09.08: GRIPT方法的case-control测试完毕，开始复现TRAPD版的burden
# 21.10.08: 重构 bam qc
# 21.10.14: 引入高速模式，调整传参模式, 重构部分代码
# 21.10.29: 集成多款calling
##################################################
import os
import sys
import time
import argparse
from script.function import f2v, trio_gt, single_gt, burden_test, bamQC, gender_identify, variants_call, anno_variants
from script.case_control import build_snvdb


def format_time(seconds):
    if seconds < 60:
        tim = '%d s.' % seconds
    elif seconds < 3600:
        m = seconds // 60
        s = seconds % 60
        tim = '%d min %d s.' % (m, s)
    elif seconds < 86400:
        h = seconds // 3600
        m = (seconds % 3600) // 60
        s = (seconds % 3600) % 60
        tim = '%d h %d min %d s.' % (h, m, s)
    else:
        d = seconds // 86400
        h = (seconds % 86400) // 3600
        m = (seconds % 86400) % 3600 // 60
        s = (seconds % 86400) % 3600 % 60
        tim = '%d day %d h %d min %d s.' % (d, h, m, s)
    return tim


def check_bed(bed):
    if bed and not os.path.isfile(bed):
        sys.exit('[ Error: Can not find bed file!]')


def f2v_args(args):
    argsd = {}
    if not args.in_dir:
        sys.exit('[ Error: Parameter is incomplete ! ]')
    else:
        argsd['in_dir'] = args.in_dir
    argsd['out_dir'] = args.out_dir
    argsd['bed'] = args.bed
    check_bed(argsd['bed'])
    argsd['prefix'] = args.prefix
    argsd['vcf'] = args.vcf
    argsd['fastqc'] = args.fastqc
    argsd['qualimap'] = args.qualimap
    argsd['fast_mark_dup'] = args.fast_mark_dup
    argsd['rm_dup'] = args.rm_dup
    argsd['fast_rm_dup'] = args.fast_rm_dup
    argsd['thread'] = args.thread
    argsd['tmp_dir'] = args.tmp_dir
    argsd['keep_tmp'] = args.keep_tmp
    argsd['config'] = args.config
    argsd['gender_rate'] = args.gender_rate
    return argsd


def bqc_args(args):
    argsd = {}
    if not args.bam:
        sys.exit('[ Error: Parameter is incomplete ! ]')
    else:
        argsd['bam'] = args.bam
    argsd['out_dir'] = args.out_dir
    argsd['bed'] = args.bed
    check_bed(argsd['bed'])
    argsd['qualimap'] = args.qualimap
    argsd['thread'] = args.thread
    argsd['tmp_dir'] = args.tmp_dir
    argsd['keep_tmp'] = args.keep_tmp
    argsd['gender_rate'] = args.gender_rate
    return argsd


def tGT_args(args):
    argsd = {}
    if not args.proband or not args.father or not args.mother:
        sys.exit('[ Error: Trio sample incomplete ! ]')
    else:
        argsd['p_gvcf'] = os.path.realpath(args.proband)
        argsd['f_gvcf'] = os.path.realpath(args.father)
        argsd['m_gvcf'] = os.path.realpath(args.mother)
        if args.sibling:
            siblings = args.sibling
            argsd['s_gvcfs'] = [os.path.realpath(_.strip(' ')) for _ in siblings.split(',')]
            print('[ Msg: find %d sibling. ]' % len(argsd['s_gvcfs']))
    argsd['out_dir'] = args.out_dir
    argsd['bed'] = args.bed
    check_bed(argsd['bed'])
    argsd['prefix'] = args.prefix
    argsd['tmp_dir'] = args.tmp_dir
    argsd['keep_tmp'] = args.keep_tmp
    argsd['config'] = args.config
    return argsd


def sGT_args(args):
    argsd = {}
    if not args.gvcf:
        sys.exit('[ Error: Sample incomplete ! ]')
    else:
        argsd['gvcf'] = os.path.realpath(args.gvcf)
    argsd['out_dir'] = args.out_dir
    argsd['bed'] = args.bed
    check_bed(argsd['bed'])
    argsd['tmp_dir'] = args.tmp_dir
    argsd['keep_tmp'] = args.keep_tmp
    argsd['prefix'] = args.prefix
    argsd['config'] = args.config
    argsd['caller'] = args.caller
    argsd['novcfqc'] = args.noqc
    argsd['noflt'] = args.noflt
    return argsd


def fp_args(args):
    argsd = {}
    if not args.in_dir:
        sys.exit('[ Error: Parameter is incomplete ! ]')
    else:
        argsd['in_dir'] = args.in_dir
    argsd['outfile'] = args.outfile
    argsd['file_type'] = args.file_type
    argsd['overlap_rate'] = args.overlap_rate
    argsd['tmp_dir'] = args.tmp_dir
    return argsd


def cc_args(args):
    argsd = {}
    if (args.case or args.case_matrix) and (args.control or args.control_matrix):
        argsd['case'] = args.case
        argsd['case_matrix'] = args.case_matrix
        argsd['control'] = args.control
        argsd['control_matrix'] = args.control_matrix
    else:
        sys.exit('[ Error: Parameter is incomplete ! ]')
    argsd['out_dir'] = args.out_dir
    argsd['cutoff'] = args.cutoff
    argsd['mode'] = args.mode
    argsd['fp'] = args.fp
    argsd['gene'] = args.gene
    argsd['score'] = args.score
    if args.gene and args.score:
        print('[ Msg: Use the metrics that user set to calculate case-control. ]')
        argsd['cc_default'] = False
    elif not args.gene and not args.score:
        argsd['cc_default'] = True
        print('[ Msg: Use default method to calculate case-control. ]')
    else:
        print('[ Warn: Parameter is incomplete ! and use default method to calculate case-control. ]')
        argsd['cc_default'] = True
    argsd['config'] = args.config
    return argsd


def gd_args(args):
    argsd = {}
    if not args.bam:
        sys.exit('[ Error: Parameter is incomplete ! ]')
    else:
        argsd['bam'] = args.bam
    argsd['bed'] = args.bed
    check_bed(argsd['bed'])
    argsd['thread'] = args.thread
    argsd['tmp_dir'] = args.tmp_dir
    argsd['rate'] = args.rate
    return argsd


def call_args(args):
    argsd = {}
    if not args.bam:
        sys.exit('[ Error: Sample incomplete ! ]')
    else:
        argsd['bam'] = os.path.realpath(args.bam)
    argsd['out_dir'] = args.out_dir
    argsd['bed'] = args.bed
    check_bed(argsd['bed'])
    argsd['thread'] = args.thread
    argsd['tmp_dir'] = args.tmp_dir
    argsd['keep_tmp'] = args.keep_tmp
    argsd['prefix'] = args.prefix
    argsd['config'] = args.config
    argsd['caller'] = args.caller
    argsd['novcfqc'] = args.noqc
    argsd['noflt'] = args.noflt
    return argsd


def anno_args(args):
    argsd = {}
    if not args.invcf:
        sys.exit('[ Error: Args incomplete ! ]')
    else:
        argsd['vcf'] = os.path.realpath(args.invcf)
    argsd['out_dir'] = args.out_dir
    argsd['thread'] = args.thread
    argsd['keep_tmp'] = args.keep_tmp
    if args.prefix:
        argsd['prefix'] = args.prefix
    else:
        argsd['prefix'] = os.path.basename(argsd['vcf']).split('.vcf')[0]
    argsd['config'] = args.config
    argsd['mode'] = args.mode
    return argsd


def ana_args():
    description = '=' * 77 + '\nVarTools 0.1.0 20211011\nWorkflow of WGS/WES analysis.\n' + '=' * 77
    print(description)
    func_description = '''
Usage:   VarTools.py <command> [options]
function of VarTools:
    (1) f2v: analysis from fastq to gvcf. 
    (2) tGT: from gvcf created by GATK to vcf in trio mode.
    (3) sGT: from gvcf created by GATK to vcf in single mode.
    (4) fp: create false positive database from vcf files or avinput files.
    (5) cc: case-control analysis with GRIPT.
    (6) bqc: bam quality check.
    (7) gd: gender identify.
    (8) call: call variants from bam.
    ----------
    (9) tSV: call trio SV with clinSV (only for WGS).
    (10) sSV: call single SV with clinSV (only for WGS).
    (11) tA: trio analysis.
    (12) sA: single case analysis.
    (13) bt: Burden testing with TRAPD.
    -----------
    (14) anno: annotation for small variants.
    To get help on a particular command, call it with -h/--help.
    '''
    function = {
        'f2v': 'analysis from fastq to gvcf.',
        'tGT': 'from gvcf created by GATK to vcf in trio mode.',
        'sGT': 'from gvcf created by GATK to vcf in single mode.',
        'bqc': 'bam quality check.',
        'fp': 'create false positive database from vcf files or avinput files.',
        'cc': 'case-control analysis with GRIPT.',
        'bt': 'Burden testing with TRAPD.',
        'tSV': 'call trio SV with clinSV (only for WGS).',
        'sSV': 'call single SV with clinSV (only for WGS).',
        'sA': 'single case analysis.',
        'tA': 'trio analysis.',
        'gd': 'identify gender from bam coverage.',
        'call': 'call variants from bam.',
        'anno': 'annotation for small variants.'
    }

    if len(sys.argv) == 1 or sys.argv[1] in ['--help', 'help', '-h']:
        sys.exit(func_description)
    elif sys.argv[1] == 'f2v':
        parser = argparse.ArgumentParser(prog='VarTool.py', usage='%(prog)s f2v [options] -i INDIR')
        parser.description = function['f2v']
        parser.add_argument('f2v')
        parser.add_argument('-i', '--in_dir', required=True, help='directory of single sample raw data.')
        parser.add_argument('-o', '--out_dir', default='./', help='output directory of result, [./].')
        parser.add_argument('-b', '--bed', default=False, help='regions of interest.')
        parser.add_argument('-p', '--prefix', default=False,
                            help='prefix of output file, if not, will use input directory name.')
        parser.add_argument('--vcf', action="store_true", help='Whether to generate vcf file.')
        parser.add_argument('--fastqc', action="store_true",
                            help='in default, fastp will give reports, '
                                 'set it if you want fastqc to check fastq with raw/clean data.')
        parser.add_argument('--qualimap', action="store_true",
                            help='bamQC with qualimap. but maybe it is slowly with large bam.')
        parser.add_argument('--fast_mark_dup', action='store_true',
                            help='use sambamba to mark duplication, but when bam is very big, may get something wrong.')
        parser.add_argument('--rm_dup', action='store_true',
                            help='remove duplication rather than mark.')
        parser.add_argument('--fast_rm_dup', action='store_true',
                            help='use fastp to remove duplication, and will skip mark duplication in follow-up steps. '
                                 'if --fast_rm_dup option is enabled, '
                                 'then --fast_mark_dup and --rm_dup options are ignored.')
        parser.add_argument('--tmp_dir', default=False,
                            help='temp directory, if not, it will create in the result directory.')
        parser.add_argument('--keep_tmp', action='store_true', help='keep temp directory.')
        parser.add_argument('-t', '--thread', default=1, type=int, help='thread of component softwares, [1].')
        parser.add_argument('--config', default=False, help='you can change config in \'lib\' or set by your need.')
        parser.add_argument('-r', '--gender_rate', default=20, type=float,
                            help='coverage rate of X/Y for calculate gender [20].')
        args = parser.parse_args()
        args_dict = f2v_args(args)
    elif sys.argv[1] == 'bqc':
        parser = argparse.ArgumentParser(prog='VarTool.py', usage='%(prog)s bqc [options] -b in.bam --bed bed')
        parser.description = function['bqc']
        parser.add_argument('bqc')
        parser.add_argument('-b', '--bam', required=True, help='bam file for qc. (after sort and index).')
        parser.add_argument('--bed', required=True,
                            help='regions of interest. If WGS file, can use bed in \'lib\' or set by your self.')
        parser.add_argument('-o', '--out_dir', default='./result', help='output directory of result, [./result].')
        parser.add_argument('--qualimap', action="store_true",
                            help='bamQC with qualimap. but maybe it is slowly with large bam.')
        parser.add_argument('-t', '--thread', default=1, type=int, help='thread of component softwares, [1].')
        parser.add_argument('--tmp_dir', default=False,
                            help='temp directory, if not, it will create in the result directory.')
        parser.add_argument('--keep_tmp', action='store_true', help='keep temp directory.')
        parser.add_argument('-r', '--gender_rate', default=20, type=float,
                            help='coverage rate of X/Y for calculate gender [20].')
        args = parser.parse_args()
        args_dict = bqc_args(args)
    elif sys.argv[1] == 'tGT':
        parser = argparse.ArgumentParser(prog='VarTool.py',
                                         usage='%(prog)s tGT [options] -p g.vcf.gz -f g.vcf.gz -m g.vcf.gz')
        parser.description = function['tGT']
        parser.add_argument('tGT')
        parser.add_argument('-p', '--proband', required=True, help='g.vcf of proband.')
        parser.add_argument('-f', '--father', required=True, help='g.vcf of father.')
        parser.add_argument('-m', '--mother', required=True, help='g.vcf of mother.')
        parser.add_argument('-s', '--sibling', help='g.vcf of siblings, more than one use \',\' to split.')
        parser.add_argument('-o', '--out_dir', default='./result', help='output directory of result, [./result].')
        parser.add_argument('-b', '--bed', default=False, help='regions of interest.')
        parser.add_argument('--prefix', default=False,
                            help='prefix of output file.[].')
        parser.add_argument('--tmp_dir', default=False,
                            help='temp directory, if not, it will create in the result directory.')
        parser.add_argument('--keep_tmp', action='store_true', help='keep temp directory.')
        parser.add_argument('--config', default=False, help='you can change config in \'lib\' or set by your need.')
        args = parser.parse_args()
        args_dict = tGT_args(args)
    elif sys.argv[1] == 'sGT':
        parser = argparse.ArgumentParser(prog='VarTool.py', usage='%(prog)s sGT [options] -g g.vcf.gz')
        parser.description = function['sGT']
        parser.add_argument('sGT')
        parser.add_argument('-g', '--gvcf', required=True, help='directory of proband raw data.')
        parser.add_argument('-o', '--out_dir', default='./result', help='output directory of result, [./result].')
        parser.add_argument('-b', '--bed', default=False, help='regions of interest.')
        parser.add_argument('-c', '--caller', default='gatk_hard', help='caller gatk_hard/vqsr, [gatk_hard].')
        parser.add_argument('-p', '--prefix', default='', help='prefix of output file.[].')
        parser.add_argument('--tmp_dir', default=False,
                            help='temp directory, if not, it will create in the result directory.')
        parser.add_argument('--keep_tmp', action='store_true', help='keep temp directory.')
        parser.add_argument('--config', default=False, help='you can change config in \'lib\' or set by your need.')
        parser.add_argument('--noqc', action='store_true', help='do not vcf quality check.')
        parser.add_argument('--noflt', action='store_true', help='do not filter raw vcf by base line.')
        args = parser.parse_args()
        args_dict = sGT_args(args)
    elif sys.argv[1] == 'fp':
        parser = argparse.ArgumentParser(prog='VarTool.py', usage='%(prog)s fp [options] -i INDIR')
        parser.description = function['fp']
        parser.add_argument('fp')
        parser.add_argument('-i', '--in_dir', required=True, help='directory of single sample raw data.')
        parser.add_argument('-o', '--outfile', default='fp.txt',
                            help='out file of result, if exist, will add result to end, [fp.txt].')
        parser.add_argument('--file_type', default='vcf',
                            help='file type , vcf/region. region:chr,start,end,ref,alt, 1-based, split by tab. [vcf].')
        parser.add_argument('--overlap_rate', default=0.5, type=float,
                            help='overlap rate of variants to build false positive database, [0.5].')
        parser.add_argument('--tmp_dir', default='./.tmp_dir_for_fp', help='temp directory [./.tmp_dir_for_fp].')
        args = parser.parse_args()
        args_dict = fp_args(args)
    elif sys.argv[1] == 'cc':
        parser = argparse.ArgumentParser(prog='VarTool.py',
                                         usage='%(prog)s cc [options] -c/-m case.dir/case.matrix '
                                               '-C/-M control.dir/control.matrix')
        parser.description = function['cc']
        parser.add_argument('cc')
        parser.add_argument('-c', '--case', help='input directory of case.')
        parser.add_argument('-C', '--control', help='input directory of control.')
        parser.add_argument('-m', '--case_matrix', help='matrix of case create by this program.')
        parser.add_argument('-M', '--control_matrix', help='matrix of control create by this program.')
        parser.add_argument('-o', '--out_dir', default='./result', help='output directory of result, [./result].')
        parser.add_argument('-t', '--cutoff', default=0, type=float, help='variant score cutoff, [0].')
        parser.add_argument('--mode', default='AD', help='mode of disease, AD/AR, [AD].')
        parser.add_argument('--fp', help='the false positive database that filted.')
        parser.add_argument('--gene', help='the columns name of gene in file.')
        parser.add_argument('--score', help='the columns name of metrics score in file.')
        parser.add_argument('--config', default=False, help='you can change config in \'lib\' or set by your need.')
        args = parser.parse_args()
        args_dict = cc_args(args)
    elif sys.argv[1] == 'gd':
        parser = argparse.ArgumentParser(prog='VarTool.py', usage='%(prog)s gd [options] -b bam -d bed')
        parser.description = function['gd']
        parser.add_argument('gd')
        parser.add_argument('-b', '--bam', required=True, help='bam file for gender identify. (after sort and index).')
        parser.add_argument('-d', '--bed', required=True,
                            help='regions of interest. If WGS file, can use bed in \'lib\' or set by your self.')
        parser.add_argument('--tmp_dir', default='./.tmp_dir_for_gd', help='temp directory [./.tmp_dir_for_gd].')
        parser.add_argument('-t', '--thread', default=1, type=int, help='thread of component softwares, [1].')
        parser.add_argument('-r', '--rate', default=20, type=float, help='coverage rate of X/Y [20].')
        args = parser.parse_args()
        args_dict = gd_args(args)
    elif sys.argv[1] == 'call':
        parser = argparse.ArgumentParser(prog='VarTool.py', usage='%(prog)s call [options] -i bam')
        parser.description = function['call']
        parser.add_argument('call')
        parser.add_argument('-i', '--bam', required=True, help='bam file for calling, must be sorted and index.')
        parser.add_argument('-o', '--out_dir', default='./result', help='output directory of result, [./result].')
        parser.add_argument('-b', '--bed', default=False, help='regions of interest.')
        parser.add_argument('-c', '--caller', default='gatk',
                            help='caller, now support gatk(_hard)/deepvariant/bcftools/vardict/strelka2 [gatk].')
        parser.add_argument('-p', '--prefix', default='result', help='prefix of output file.[result].')
        parser.add_argument('--tmp_dir', default=False,
                            help='temp directory, if not, it will create in the result directory.')
        parser.add_argument('--keep_tmp', action='store_true', help='keep temp directory.')
        parser.add_argument('--config', default=False, help='you can change config in \'lib\' or set by your need.')
        parser.add_argument('--noqc', action='store_true', help='do not vcf quality check.')
        parser.add_argument('-t', '--thread', default=1, type=int, help='thread of component softwares, [1].')
        parser.add_argument('--noflt', action='store_true', help='do not filter raw vcf by base line.')
        args = parser.parse_args()
        args_dict = call_args(args)
    elif sys.argv[1] == 'anno':
        parser = argparse.ArgumentParser(prog='VarTool.py', usage='%(prog)s anno [options] -i VCF')
        parser.description = function['anno']
        parser.add_argument('-i', '--invcf', required=True, help='input vcf file.')
        parser.add_argument('-o', '--out_dir', default='./', help='output directory of result, [./].')
        parser.add_argument('-p', '--prefix', default=False,
                            help='prefix of output file, if not, will use input vcf name.')
        parser.add_argument('--keep_tmp', action='store_true', help='keep temp directory.')
        parser.add_argument('--config', default=False, help='you can change config in \'lib\' or set by your need.')
        parser.add_argument('-t', '--thread', default=1, type=int, help='thread of component softwares, [1].')
        parser.add_argument('-m', '--mode', default='FA', choices=['FA', 'TA'],
                            help='two mode TA/FA, TA: annotate all info, FA: filter with annotation,[FA].')
        args = parser.parse_args()
        args_dict = anno_args(args)
    else:
        sys.exit('[ Error: Can not identify the function of <%s>]' % sys.argv[1])

    # 脚本所在路径
    args_dict['fun'] = sys.argv[1]
    args_dict['script_path'] = os.path.split(os.path.realpath(__file__))[0]
    return args_dict


def main():
    # 设置参数
    args = ana_args()
    if args['fun'] == 'f2v':
        f2v(args['in_dir'], args['out_dir'], args['bed'], args['prefix'],
            args['vcf'], args['fastqc'], args['qualimap'],
            args['fast_mark_dup'], args['rm_dup'], args['fast_rm_dup'],
            args['thread'], args['script_path'], args['config'], args['tmp_dir'], args['keep_tmp'], args['gender_rate'])

    elif args['fun'] == 'bqc':
        bamQC(args['bam'], args['bed'], args['out_dir'], args['tmp_dir'],
              args['script_path'], args['thread'], args['qualimap'], args['keep_tmp'], args['gender_rate'])

    elif args['fun'] == 'tGT':
        trio_gt(args['p_gvcf'], args['f_gvcf'], args['m_gvcf'], args['s_gvcfs'],
                args['out_dir'], args['script_path'], args['config'],
                args['keep_tmp'], args['tmp_dir'], args['bed'], args['prefix'])

    elif args['fun'] == 'sGT':
        single_gt(args['gvcf'], args['out_dir'], args['script_path'], args['bed'], args['tmp_dir'], args['keep_tmp'],
                  args['prefix'], args['config'], args['caller'], args['novcfqc'], args['noflt'])

    elif args['fun'] == 'fp':
        rm_tmp = False
        if not os.path.isdir(args['tmp_dir']):
            rm_tmp = True
            os.makedirs(args['tmp_dir'])
        tmp_dir = os.path.abspath(args['tmp_dir'])
        build_snvdb(args['in_dir'], args['outfile'], tmp_dir,
                    args['script_path'], args['file_type'], args['overlap_rate'], rm_tmp)

    elif args['fun'] == 'cc':
        burden_test(args['case'], args['control'], args['case_matrix'], args['control_matrix'],
                    args['out_dir'], args['fp'], args['mode'], args['cutoff'],
                    args['cc_default'], args['gene'], args['score'], args['script_path'], args['config'])

    elif args['fun'] == 'gd':
        gender_identify(args['bam'], args['bed'], args['tmp_dir'], args['thread'], args['script_path'], args['rate'])
    elif args['fun'] == 'tA':
        pass
    elif args['fun'] == 'sA':
        pass
    elif args['fun'] == 'call':
        variants_call(args['bam'], args['out_dir'], args['caller'], args['bed'], args['prefix'], args['thread'],
                      args['tmp_dir'], args['keep_tmp'], args['config'], args['script_path'],
                      args['novcfqc'], args['noflt'])
    elif args['fun'] == 'anno':
        anno_variants(args['vcf'], args['prefix'], args['out_dir'],
                      args['script_path'], args['config'], args['thread'], args['mode'])


if __name__ == '__main__':
    start_time = time.perf_counter()
    main()
    end_words = '=' * 77 + '\nThanks for using VarTools! \nYou can report bugs to tongshiyuan@foxmail.com\n' + '=' * 77
    end_time = time.perf_counter()
    tim = format_time(end_time - start_time)
    print('[ Msg: Use time : %s ]' % tim)
    print(end_words)
