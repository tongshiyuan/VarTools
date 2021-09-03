#!/usr/bin/python3
# -*- coding:utf-8 -*-

##################################################
# Author: Shi-Yuan Tong
# Email: tongshiyuan@foxmail.com
# Created Time: 2021-5-14
# log:
# 21.05.31: 根据新的框架修改，为了拓展性、适配性更好
# 21.08.31: 完善框架，修改了一些已知错误
##################################################
import os
import sys
import time
import argparse
from script.function import f2v, trio_gt, single_gt


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
    (8) cc: case-control analysis with 2 ways.
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
    parser.add_argument('-v', '--vcf', help='Whether to generate vcf file, True/False, [False]')
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
    parser.add_argument('-i', '--', help='')
    parser.add_argument('-o', '--', help='')

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
        if args.vcf == 'True':
            args_dict['vcf'] = True
        elif args.vcf == 'False' or not args.vcf:
            args_dict['vcf'] = False
        else:
            print('[ W: Can not identify vcf parameter, and use "False". ]')
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
        pass
    else:
        sys.exit('[ Err: can not identify the function of %s]' % args['fun'])
    end_time = time.perf_counter()
    print('[ Msg: Use time : <%d> s]' % (end_time - start_time))
    print(end_words)


if __name__ == '__main__':
    main()
