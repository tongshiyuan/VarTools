import os
import sys
import configparser
from script.qc import fastp_qc
from script.mapping import align_deal
from script.calling import gatk_pre, gatk, gatk_hard
from script.common import check_software, affinity
from script.annotation import trio_short_variants_filter
from script.case_control import cc_preprocess, div_cc


def readConfig(configFile):
    config = configparser.ConfigParser()
    config.read(configFile)
    confDict = {}
    confDict['reference'] = config.get('global', 'reference')
    confDict['adapterR1'] = config.get('fastqc', 'adapter_r1')
    confDict['adapterR2'] = config.get('fastqc', 'adapter_r2')
    confDict['gatkBundle'] = config.get('database', 'GATK_bundle')
    confDict['mapping'] = config.get('parameter', 'mapping')
    confDict['qcProcess'] = int(config.get('fastqc', 'qc_process'))
    confDict['alnProcess'] = int(config.get('map', 'map_process'))
    confDict['af_db'] = config.get('anno', 'af_db')
    confDict['af_type'] = config.get('anno', 'af_type')
    confDict['gene_db'] = config.get('anno', 'gene_db')
    confDict['gene_type'] = config.get('anno', 'gene_type')
    confDict['anno_dir'] = config.get('anno', 'anno_dir')
    confDict['ref_version'] = config.get('anno', 'ref_version')
    _ = config.get('filter', 'AFList')
    confDict['AFList'] = [i.strip() for i in _.split(',')]
    confDict['AFTh'] = eval(config.get('filter', 'AFThreshold'))
    _ = config.get('filter', 'ClinList')
    confDict['ClinList'] = [i.strip() for i in _.split(',')]
    return confDict


def f2v(inDir, outDir, thread, scriptPath, vcf):
    rst = check_software('samtools')
    if rst:
        sys.exit()
    # 读取参数
    config = readConfig(scriptPath + '/lib/config.ini')
    # 输出目录
    inDir = os.path.abspath(inDir)
    outDir = os.path.abspath(outDir)
    sampleName = os.path.split(os.path.abspath(inDir))[-1]
    # QC
    print('[ Msg: QC running ... ]')
    fqOutDir = outDir + '/' + sampleName + '/tmp'
    os.makedirs(fqOutDir)
    fqReportDir = outDir + '/' + sampleName + '/reports/fastqQC'
    os.makedirs(fqReportDir)
    fastp_qc(inDir, fqOutDir, fqReportDir, config['adapterR1'], config['adapterR2'], thread, scriptPath,
             config['qcProcess'])
    # mapping
    print('[ Msg: Waiting for all sample mapping done ... ]')
    bamOutdir = outDir + '/' + sampleName + '/result'
    os.makedirs(bamOutdir)
    bam_report_dir = outDir + '/' + sampleName + '/reports/bamQC'
    os.makedirs(bam_report_dir)
    bam_tmp_dir = outDir + '/' + sampleName + '/tmp'
    # os.makedirs(bam_tmp_dir)
    bamfile = align_deal(fqOutDir, bamOutdir, sampleName, config['reference'], bam_report_dir, bam_tmp_dir,
                         config['mapping'], config['alnProcess'], thread, scriptPath)
    print('[ Msg: All sample mapping done ! ]')
    print('[ Msg: Waiting for gatk calling done ... ]')
    # short variants calling
    gvcf_outdir = outDir + '/' + sampleName + '/result'
    # os.makedirs(gvcf_outdir)
    gvcf_tmp_dir = outDir + '/' + sampleName + '/tmp'
    # os.makedirs(gvcf_tmp_dir)
    gatk_pre(bamfile, gvcf_outdir, gvcf_tmp_dir, config['reference'], sampleName, config['gatkBundle'], scriptPath, vcf)


def trio_gt(p_gvcf, f_gvcf, m_gvcf, s_gvcfs, outDir, scriptPath):
    config = readConfig(scriptPath + '/lib/config.ini')
    # 输入目录
    if s_gvcfs:
        # s_gvcfs = s_gvcfs.split(',')
        # for i, v in enumerate(s_gvcfs):
        #     s_gvcfs[i] = os.path.realpath(v)
        gatk_gvcf = [p_gvcf, f_gvcf, m_gvcf] + s_gvcfs
    else:
        gatk_gvcf = [p_gvcf, f_gvcf, m_gvcf]
    # 输出目录
    outDir = os.path.abspath(outDir)  # + '/result'
    report_dir = outDir + '/variantQC'
    aff_dir = outDir + '/affinity'
    vcfDir = outDir + '/vcf'
    os.makedirs(vcfDir)
    os.makedirs(report_dir)
    os.makedirs(aff_dir)
    vcf = gatk(gatk_gvcf, vcfDir, report_dir, config['reference'], config['gatkBundle'], scriptPath)
    affinity(vcf, aff_dir, scriptPath)


def single_gt(gvcf, outDir, scriptPath):
    config = readConfig(scriptPath + '/lib/config.ini')
    outDir = os.path.abspath(outDir)  # + '/result'
    report_dir = outDir + '/variantQC'
    vcfDir = outDir + '/vcf'
    os.makedirs(vcfDir)
    os.makedirs(report_dir)
    gatk_hard(gvcf, vcfDir, report_dir, config['reference'], scriptPath)


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


def burden_test(case, control, out_dir, snvdb, mode, cutoff, method, gene, score):
    tmp_dir = out_dir + '/tmp'
    os.makedirs(tmp_dir)
    if method:
        out_case, out_control = cc_preprocess(case, control, tmp_dir, snvdb, mode)
        df = div_cc(out_case, out_control, cutoff=cutoff)
    else:
        df = div_cc(case, control, gene, score, cutoff)
    df.to_csv(out_dir + '/result.txt', sep='\t', index=False)
