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
##################################################
import os
import sys
import time
import argparse
from script.function import f2v, trio_gt, single_gt, burden_test
from script.case_control import build_snvdb


def ana_args():
    args_dict = {}
    description = '=' * 77 + '\nVarTools 0.1.0 20210513\nWorkflow of WGS/WES trio analysis.\n' + '=' * 77
    func = '''
    function of VarTools:
    (1) f2v: analysis from fastq to gvcf. 
    (2) tGT: from gvcf created by GATK to vcf in trio mode.
    (3) sGT: from gvcf created by GATK to vcf in single mode.
    (4) tSV: call trio SV with clinSV (only for WGS).
    (5) sSV: call single SV with clinSV (only for WGS).
    (6) tA: trio analysis.
    (7) sA: single case analysis.
    (8) fp: create false positive database from vcf files or avinput files.
    (9) cc: case-control analysis with GRIPT.
    (10) bt: Burden testing with TRAPD.
    
    '''
    print(description)
    parser = argparse.ArgumentParser()
    parser.description = description
    parser.add_argument('function', help=func)
    # common
    parser.add_argument('-o', '--output', help='output directory of result, [./]')
    parser.add_argument('-t', '--thread', help='thread of component softwares, [6]')
    # f2v
    parser.add_argument('-i', '--indir', help='directory of single sample raw data')
    parser.add_argument('--vcf', action="store_true", help='Whether to generate vcf file.')
    # tGT && sGT
    parser.add_argument('-p', '--proband', help='directory of proband raw data')
    parser.add_argument('-f', '--father', help='directory of father raw data')
    parser.add_argument('-m', '--mother', help='directory of mother raw data')
    parser.add_argument('-s', '--sibling', help='directory of siblings raw data, more one use \',\' split')
    # trio
    # parser.add_argument('-v', '--vcf', help='vcf file')
    # parser.add_argument('-p', '--proband', help='sample name of proband in vcf')
    # parser.add_argument('-f', '--father', help='directory of father raw data')
    # parser.add_argument('-m', '--mother', help='directory of mother raw data')
    # parser.add_argument('-s', '--sibling', help='directory of siblings raw data, more one use \',\' split')
    # case-control
    parser.add_argument('--case', help='input directory of case')
    parser.add_argument('--control', help='input directory of control')
    parser.add_argument('--case_matrix', help='matrix of case create by this program')
    parser.add_argument('--control_matrix', help='matrix of control create by this program')
    parser.add_argument('--cutoff', default=0, type=float, help='variant score cutoff. [0]')
    parser.add_argument('--mode', default='AD', help='mode of disease, AD/AR, [AD]')
    parser.add_argument('--snvdb', help='the false positive database that filted')
    parser.add_argument('--gene', help='the columns name of gene in file')
    parser.add_argument('--score', help='the columns name of metrics score in file')
    # build false positive data base
    parser.add_argument('--file_type', default='vcf', help='file type of input files, vcf/avinput. [vcf]')
    parser.add_argument('--overlap_rate', default=0.5, type=float,
                        help='overlap rate of variants to build false positive database, [0.5]')
    args = parser.parse_args()
    # thread
    if not args.thread:
        args_dict['thread'] = 8
    else:
        try:
            thread = int(args.thread)
            args_dict['thread'] = thread
        except:
            sys.exit('[ E: Please input integer in thread ! ]')
    # 脚本所在路径
    args_dict['scriptPath'] = os.path.split(os.path.realpath(__file__))[0]
    # 输出目录，会重新创建一个
    if not args.output:
        args_dict['outPath'] = './'
    else:
        args_dict['outPath'] = args.output
    # 输入 & fun
    args_dict['fun'] = args.function
    if args_dict['fun'] == 'f2v':
        if not args.indir:
            sys.exit('[ E: Parameter is incomplete ! ]')
        else:
            args_dict['inDir'] = args.indir
        if args.vcf:
            args_dict['vcf'] = True
        else:
            args_dict['vcf'] = False
    elif args_dict['fun'] == 'tGT':
        if args.proband and args.father and args.mother:
            args_dict['pgVCF'] = os.path.realpath(args.proband)
            args_dict['fgVCF'] = os.path.realpath(args.father)
            args_dict['mgVCF'] = os.path.realpath(args.mother)
            if args.sibling:
                siblings = args.sibling
                args_dict['sgVCF'] = [os.path.realpath(_.strip(' ')) for _ in siblings.split(',')]
                print('[ M: find %d sibling. ]' % len(args_dict['sgVCF']))
        else:
            sys.exit('[ E: trio sample incomplete ! ]')
    elif args_dict['fun'] == 'sGT':
        args_dict['pgVCF'] = os.path.realpath(args.proband)
    elif args_dict['fun'] == 'cc':
        if (args.case or args.case_matrix) and (args.control or args.control_matrix):
            args_dict['case'] = args.case
            args_dict['case_matrix'] = args.case_matrix
            args_dict['control'] = args.control
            args_dict['control_matrix'] = args.control_matrix
        else:
            sys.exit('[ E: Parameter is incomplete ! ]')
        args_dict['cutoff'] = args.cutoff
        args_dict['mode'] = args.mode
        args_dict['snvdb'] = args.snvdb
        args_dict['gene'] = args.gene
        args_dict['score'] = args.score
        if args.gene and args.score:
            print('[ M: Use the metrics that user set to calculate case-control. ]')
            args_dict['cc_default'] = False
        elif not args.gene and not args.score:
            args_dict['cc_default'] = True
            print('[ M: Use default method to calculate case-control. ]')
        else:
            print('[ W: Parameter is incomplete ! and use default method to calculate case-control. ]')
            args_dict['cc_default'] = True
    elif args_dict['fun'] == 'fp':
        args_dict['file_type'] = args.file_type
        args_dict['overlap_rate'] = args.overlap_rate
    else:
        sys.exit('[ Err: can not identify the function of <%s>]' % args_dict['fun'])

    return args_dict


def main():
    start_time = time.perf_counter()
    end_words = '=' * 77 + '\nThanks for using VarTools! \nReport bugs to tongshiyuan@foxmail.com\n' + '=' * 77
    # 设置参数
    args = ana_args()
    if args['fun'] == 'f2v':
        f2v(args['inDir'], args['outPath'], args['thread'], args['scriptPath'], args['vcf'])
    elif args['fun'] == 'tGT':
        trio_gt(args['pgVCF'], args['fgVCF'], args['mgVCF'], args['sgVCF'], args['outPath'], args['scriptPath'])
    elif args['fun'] == 'sGT':
        single_gt(args['pgVCF'], args['outPath'], args['scriptPath'])
    elif args['fun'] == 'tA':
        pass
    elif args['fun'] == 'sA':
        pass
    elif args['fun'] == 'cc':
        burden_test(args['case'], args['control'], args['case_matrix'], args['control_matrix'],
                    args['outPath'], args['snvdb'], args['mode'], args['cutoff'],
                    args['cc_default'], args['gene'], args['score'], args['scriptPath'])
    elif args['fun'] == 'fp':
        build_snvdb(args['inDir'], args['outPath'], args['snvdb'],
                    args['scriptPath'], args['file_type'], args['overlap_rate'],
                    )
    end_time = time.perf_counter()
    print('[ Msg: Use time : <%d> s]' % (end_time - start_time))
    print(end_words)


if __name__ == '__main__':
    main()
