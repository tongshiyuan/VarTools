import os
import sys
from script.common import execute_system


# 输入是bam文件，
# 输出是对应模式的vcf
# 所有单个软件输入bam文件，输出vcf文件
def gatk_pre(bam, outdir, tmpDir, reference, sampleName, gatk_bundle_dir, scriptPath, vcf=False):
    # outdir = os.path.abspath(outdir)
    # gatk_bundle_dir = os.path.abspath(gatk_bundle_dir)
    tmp_dir = os.path.abspath(tmpDir) + '/'
    _tmp_file1 = tmp_dir + sampleName + '.recal_data.table'
    _tmp_file2 = tmp_dir + sampleName + '.BQSR.bam'
    g_vcf = outdir + '/' + sampleName + '.gatk.g.vcf.gz'
    # call snvs/indels
    # BQSR
    # BaseRecalibrator
    if os.path.basename(reference).find('hg19') != -1:
        _1kg_indel = '--known-sites ' + gatk_bundle_dir + '/1000G_phase1.indels.hg19.sites.vcf '
        _mill = '--known-sites ' + gatk_bundle_dir + '/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf '
        _dbsnp = '--known-sites ' + gatk_bundle_dir + '/dbsnp_138.hg19.vcf '
        _1kg_snv = '--known-sites ' + gatk_bundle_dir + '/1000G_phase1.snps.high_confidence.hg19.sites.vcf'
        gatk_db = _1kg_indel + _mill + _dbsnp + _1kg_snv
    elif os.path.basename(reference).find('v37') != -1:
        _1kg_indel = '--known-sites ' + gatk_bundle_dir + '/1000G_phase1.indels.b37.vcf '
        _mill = '--known-sites ' + gatk_bundle_dir + '/Mills_and_1000G_gold_standard.indels.b37.vcf '
        _dbsnp = '--known-sites ' + gatk_bundle_dir + '/dbsnp_138.b37.vcf '
        _1kg_snv = '--known-sites ' + gatk_bundle_dir + '/1000G_phase1.snps.high_confidence.b37.vcf'
        gatk_db = _1kg_indel + _mill + _dbsnp + _1kg_snv
    elif os.path.basename(reference).find('38') != -1:
        _mill = '--known-sites ' + gatk_bundle_dir + '/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz '
        _dbsnp = '--known-sites ' + gatk_bundle_dir + '/dbsnp_146.hg38.vcf.gz '
        _1kg_vcf = '--known-sites ' + gatk_bundle_dir + '/1000G_phase1.snps.high_confidence.hg38.vcf.gz'
        gatk_db = _mill + _dbsnp + _1kg_vcf
    else:
        sys.exit('[ E: Can not find bundle file with <%s>]' % os.path.basename(reference))
    BR_cmd = scriptPath + '/bin/gatk/gatk BaseRecalibrator -R %s -I %s %s -O %s' % (
        reference, bam, gatk_db, _tmp_file1)
    execute_system(BR_cmd, '[ Msg: <%s> baseRecalibrator done ! ]' % sampleName,
                   '[ E: Something wrong with <%s> BaseRecalibrator ! ]' % sampleName)
    # ApplyBQSRCmd
    apply_BQSR_cmd = scriptPath + '/bin/gatk/gatk ApplyBQSR --bqsr-recal-file %s -R %s -I %s -O %s' % (
        _tmp_file1, reference, bam, _tmp_file2)
    execute_system(apply_BQSR_cmd, '[ Msg: <%s> ApplyBQSR done ! ]' % sampleName,
                   '[ E: Something wrong with <%s> ApplyBQSR ! ]' % sampleName)
    # index
    index_cmd = 'samtools index %s' % _tmp_file2
    execute_system(index_cmd, '[ Msg: <%s> BQSR bam index done ! ]' % sampleName,
                   '[ E: Something wrong with <%s> BQSR bam index ! ]' % sampleName)
    # HaplotypeCaller
    HC_cmd = scriptPath + '/bin/gatk/gatk HaplotypeCaller --emit-ref-confidence GVCF -R %s -I %s -O %s' % (
        reference, bam, g_vcf)
    execute_system(HC_cmd, '[ Msg: <%s> HaplotypeCaller done ! ]' % sampleName,
                   '[ E: Something wrong with <%s> HaplotypeCaller ! ]' % sampleName)
    if vcf:
        raw_vcf = outdir + '/' + sampleName + '.gatk.raw.vcf.gz'
        snp_vcf = outdir + '/' + sampleName + '.gatk.snp.vcf.gz'
        indel_vcf = outdir + '/' + sampleName + '.gatk.indel.vcf.gz'
        filter_snp = outdir + '/' + sampleName + '.gatk.snp.filter.vcf.gz'
        filter_indel = outdir + '/' + sampleName + '.gatk.indel.filter.vcf.gz'
        final_vcf = outdir + '/' + sampleName + '.gatk.final.vcf.gz'
        genotype_cmd = scriptPath + '/bin/gatk/gatk GenotypeGVCFs -R %s -V %s -O %s' % (
            reference, g_vcf, raw_vcf)
        execute_system(genotype_cmd, '[ Msg: Genotype GVCFs done ! ]',
                       '[ E: Something wrong with genotype GVCFs ! ]')
        # select snp
        sv_cmd = scriptPath + '/bin/gatk/gatk SelectVariants -select-type SNP -V %s -O %s' % (
            raw_vcf, snp_vcf)
        execute_system(sv_cmd, '[ Msg: Select SNP done ! ]',
                       '[ E: Something wrong with select SNP ! ]')
        # select indel
        sv_cmd = scriptPath + '/bin/gatk/gatk SelectVariants -select-type INDEL -V %s -O %s' % (
            raw_vcf, indel_vcf)
        execute_system(sv_cmd, '[ Msg: Select InDel done ! ]',
                       '[ E: Something wrong with select InDel ! ]')
        # snp filter
        snp_filter_cmd = scriptPath + '/bin/gatk/gatk VariantFiltration -V %s -O %s ' % (snp_vcf, filter_snp) + \
                         '-filter "QD < 2.0" --filter-name "QD2" ' + \
                         '-filter "QUAL < 30.0" --filter-name "QUAL30" ' + \
                         '-filter "SOR > 3.0" --filter-name "SOR3" ' + \
                         '-filter "FS > 60.0" --filter-name "FS60" ' + \
                         '-filter "MQ < 40.0" --filter-name "MQ40" ' + \
                         '-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" ' + \
                         '-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"'
        execute_system(snp_filter_cmd, '[ Msg: Filter snp done ! ]',
                       '[ E: Something wrong with filter snp ! ]')
        # indel filter
        indel_filter_cmd = scriptPath + '/bin/gatk/gatk VariantFiltration -V %s -O %s ' % (indel_vcf, filter_indel) + \
                           '-filter "QD < 2.0" --filter-name "QD2" ' + \
                           '-filter "QUAL < 30.0" --filter-name "QUAL30" ' + \
                           '-filter "FS > 200.0" --filter-name "FS200" ' + \
                           '-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"'
        execute_system(indel_filter_cmd, '[ Msg: Filter indel done ! ]',
                       '[ E: Something wrong with filter indel ! ]')
        # merge vcf
        merge_cmd = scriptPath + '/bin/gatk/gatk MergeVcfs -I %s -I %s -O %s' % (filter_snp, filter_indel, final_vcf)
        execute_system(merge_cmd, '[ Msg: merge variants done ! ]',
                       '[ E: Something wrong with merge variants ! ]')
        # 删除中间文件
        rm_cmd = 'rm -f %s* %s* %s* %s* ' % (snp_vcf, filter_snp, indel_vcf, filter_indel)
        execute_system(rm_cmd, '[ Msg: Delete process file done in gatk! ]',
                       '[ E: Something wrong with delete process file in gatk ! ]')
        # vcf stats
        report_dir = outdir + '/variantQC'
        os.makedirs(report_dir)
        vcf_stats(final_vcf, report_dir, reference, scriptPath)


def gatk(gvcf_list, outdir, report_dir, reference, gatk_bundle_dir, scriptPath):
    # outdir = os.path.abspath(outdir)
    # report_dir = os.path.abspath(report_dir)
    # gatk_bundle_dir = os.path.abspath(gatk_bundle_dir)
    # reference = os.path.abspath(reference)
    merged_gvcf = outdir + '/cohort.g.vcf.gz'
    cohortVCF = outdir + '/cohort.vcf.gz'
    excesshet = outdir + '/cohort_excesshet.vcf.gz'
    sitesonly = outdir + '/cohort_sitesonly.vcf.gz'
    indels_recal = outdir + '/cohort_indels.recal'
    indels_tranches = outdir + '/cohort_indels.tranches'
    indels_vcf = outdir + '/indel.recalibrated.vcf.gz'
    snp_recal = outdir + '/cohort_snps.recal'
    snp_tranches = outdir + '/cohort_snps.tranches'
    final_vcf = outdir + '/cohort.filter.vcf.gz'
    # mergeGVCFs
    sample_gvcf = ''
    for i in gvcf_list:
        sample_gvcf += '-V %s ' % i
    merge_cmd = scriptPath + '/bin/gatk/gatk CombineGVCFs -R %s %s -O %s' % (reference, sample_gvcf, merged_gvcf)
    execute_system(merge_cmd, '[ Msg: GATK combine gvcfs done ! ]',
                   '[ E: Something wrong with combine gvcfs in GATK ! ]')
    # GenotypeGVCFs
    genotype_cmd = scriptPath + '/bin/gatk/gatk GenotypeGVCFs -R %s -V %s -O %s' % (reference, merged_gvcf, cohortVCF)
    execute_system(genotype_cmd, '[ Msg: Genotype GVCFs done ! ]',
                   '[ E: Something wrong with genotype GVCFs ! ]')
    # VariantFiltration
    filtration_cmd = scriptPath + '/bin/gatk/gatk VariantFiltration -V %s --filter-expression "ExcessHet > 54.69" ' \
                                  '--filter-name ExcessHet -O %s' % (cohortVCF, excesshet)
    execute_system(filtration_cmd, '[ Msg: Cohort VCF variantFiltration done ! ]',
                   '[ E: Something wrong with cohort VCF variantFiltration ! ]')
    # MakeSitesOnly
    sites_only_cmd = scriptPath + '/bin/gatk/gatk MakeSitesOnlyVcf -I %s -O %s' % (excesshet, sitesonly)
    execute_system(sites_only_cmd, '[ Msg: Make sites only vcf done ! ]',
                   '[ E: Something wrong with make sites only vcf ! ]')
    # VariantRecalibrator
    if os.path.basename(reference).find('hg19') != -1:
        _mill = '-resource:mills,known=false,training=true,truth=true,prior=12 ' \
                '%s/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf ' % gatk_bundle_dir
        _dbsnp = '-resource:dbsnp,known=true,training=false,truth=false,prior=2 ' \
                 '%s/dbsnp_138.hg19.vcf' % gatk_bundle_dir
        _dbsnp2 = '%s/dbsnp_138.hg19.vcf' % gatk_bundle_dir
        _hapmap = '%s/hapmap_3.3.hg19.sites.vcf' % gatk_bundle_dir
        _omni = '%s/1000G_omni2.5.hg19.sites.vcf' % gatk_bundle_dir
        _1kg = '%s/1000G_phase1.snps.high_confidence.hg19.sites.vcf' % gatk_bundle_dir
        indel_db = _mill + _dbsnp
    elif os.path.basename(reference).find('v37') != -1:
        _mill = '-resource:mills,known=false,training=true,truth=true,prior=12 ' \
                '%s/Mills_and_1000G_gold_standard.indels.b37.vcf' % gatk_bundle_dir
        _dbsnp = '-resource:dbsnp,known=true,training=false,truth=false,prior=2 ' \
                 '%s/dbsnp_138.b37.vcf' % gatk_bundle_dir
        _dbsnp2 = '%s/dbsnp_138.b37.vcf' % gatk_bundle_dir
        _hapmap = '%s/hapmap_3.3.b37.vcf' % gatk_bundle_dir
        _omni = '%s/1000G_omni2.5.b37.vcf' % gatk_bundle_dir
        _1kg = '%s/1000G_phase1.snps.high_confidence.b37.vcf' % gatk_bundle_dir
        indel_db = _mill + _dbsnp
    elif os.path.basename(reference).find('38') != -1:
        _mill = '-resource:mills,known=false,training=true,truth=true,prior=12 ' \
                '%s/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz ' % gatk_bundle_dir
        _dbsnp = '-resource:dbsnp,known=true,training=false,truth=false,prior=2 ' \
                 '%s/dbsnp_146.hg38.vcf.gz' % gatk_bundle_dir
        _axiom = '-resource:axiomPoly,known=false,training=true,truth=false,prior=10 ' \
                 '%s/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz ' % gatk_bundle_dir
        _dbsnp2 = '%s/dbsnp_146.hg38.vcf.gz' % gatk_bundle_dir
        _hapmap = '%s/hapmap_3.3.hg38.vcf.gz' % gatk_bundle_dir
        _omni = '%s/1000G_omni2.5.hg38.vcf.gz' % gatk_bundle_dir
        _1kg = '%s/1000G_phase1.snps.high_confidence.hg38.vcf.gz' % gatk_bundle_dir
        indel_db = _mill + _axiom + _dbsnp
    else:
        sys.exit('[ E: Can not find bundle file with <%s>]' % os.path.basename(reference))
    recalibratorCmd = scriptPath + '/bin/gatk/gatk VariantRecalibrator -V %s --trust-all-polymorphic ' \
                                   '-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 ' \
                                   '-tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 ' \
                                   '-tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 ' \
                                   '-an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP -mode INDEL ' \
                                   '--max-gaussians 4 %s -O %s --tranches-file %s' % (
                          sitesonly, indel_db, indels_recal, indels_tranches)
    execute_system(recalibratorCmd, '[ Msg: Indel recalibrator done ! ]',
                   '[ E: Something wrong with indel recalibrator ! ]')
    apply_vqsr_cmd = scriptPath + '/bin/gatk/gatk ApplyVQSR -V %s --recal-file %s --tranches-file %s ' \
                                  '--truth-sensitivity-filter-level 99.7 --create-output-variant-index true ' \
                                  '-mode INDEL -O %s' % (excesshet, indels_recal, indels_tranches, indels_vcf)
    execute_system(apply_vqsr_cmd, '[ Msg: Indel apply VQSR done ! ]',
                   '[ E: Something wrong with indel apply VQSR ! ]')
    recalibrator_cmd = scriptPath + '/bin/gatk/gatk VariantRecalibrator -V %s --trust-all-polymorphic ' \
                                    '-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 ' \
                                    '-tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 ' \
                                    '-tranche 97.0 -tranche 90.0 ' \
                                    '-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP -mode SNP ' \
                                    '--max-gaussians 6 ' \
                                    '-resource:hapmap,known=false,training=true,truth=true,prior=15 %s ' \
                                    '-resource:omni,known=false,training=true,truth=true,prior=12 %s ' \
                                    '-resource:1000G,known=false,training=true,truth=false,prior=10 %s ' \
                                    '-resource:dbsnp,known=true,training=false,truth=false,prior=7 %s ' \
                                    '-O %s --tranches-file %s' % (
                           sitesonly, _hapmap, _omni, _1kg, _dbsnp2, snp_recal, snp_tranches)
    execute_system(recalibrator_cmd, '[ Msg: Snp recalibrator done ! ]',
                   '[ E: Something wrong with snp recalibrator ! ]')
    filter_cmd = scriptPath + '/bin/gatk/gatk ApplyVQSR -V %s --recal-file %s --tranches-file %s ' \
                              '--truth-sensitivity-filter-level 99.7 ' \
                              '--create-output-variant-index true -mode SNP -O %s' % (
                     indels_vcf, snp_recal, snp_tranches, final_vcf)
    execute_system(filter_cmd, '[ Msg: snps recalibrated done ! ]', '[ E: Something wrong with snps recalibrated ! ]')
    # 建立索引
    # index_cmd = 'bcftools index -t %s' % final_vcf
    # execute_system(index_cmd, '[ Msg: Build vcf index done in gatk ! ]',
    #                '[ E: Something wrong with build vcf index file in gatk ! ]')
    # 删除无用中间文件
    rm_cmd = 'rm -f %s %s %s %s %s %s %s %s' % (
        sitesonly, excesshet, merged_gvcf, indels_tranches, indels_recal, indels_vcf, snp_recal,
        snp_tranches)
    execute_system(rm_cmd, '[ Msg: Delete process file done in gatk! ]',
                   '[ E: Something wrong with delete process file in gatk ! ]')
    # vcf stats
    vcf_stats(final_vcf, report_dir, reference, scriptPath)
    return final_vcf


def gatk_hard(gvcf, outdir, report_dir, reference, scriptPath):
    raw_vcf = outdir + '/raw.vcf.gz'
    snp_vcf = outdir + '/snvs.vcf.gz'
    snp_filted_vcf = outdir + '/snvs_filtered.vcf.gz'
    indelVcf = outdir + '/indels.vcf.gz'
    indel_filted_vcf = outdir + '/indels_filtered.vcf.gz'
    vcf = outdir + '/filtered.vcf.gz'
    genotype_cmd = scriptPath + '/bin/gatk/gatk GenotypeGVCFs -R %s -V %s -O %s' % (reference, gvcf, raw_vcf)
    execute_system(genotype_cmd, '[ Msg: Genotype GVCF done ! ]', '[ E: Something wrong with genotype GVCF ! ]')
    select_snp_cmd = scriptPath + '/bin/gatk/gatk SelectVariants -V %s -select-type SNP -O %s' % (raw_vcf, snp_vcf)
    execute_system(select_snp_cmd, '[ Msg: select SNV done ! ]', '[ E: Something wrong with select SNV ! ]')
    select_indel_cmd = scriptPath + '/bin/gatk/gatk SelectVariants -V %s -select-type INDEL -O %s' % (raw_vcf, indelVcf)
    execute_system(select_indel_cmd, '[ Msg: select INDEL done ! ]', '[ E: Something wrong with select INDEL ! ]')
    snv_filter_cmd = scriptPath + '/bin/gatk/gatk VariantFiltration ' \
                                  '-filter "QD < 2.0" --filter-name "QD2" ' \
                                  '-filter "QUAL < 30.0" --filter-name "QUAL30" ' \
                                  '-filter "SOR > 3.0" --filter-name "SOR3" ' \
                                  '-filter "FS > 60.0" --filter-name "FS60" ' \
                                  '-filter "MQ < 40.0" --filter-name "MQ40" ' \
                                  '-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" ' \
                                  '-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" ' \
                                  '-V %s -O %s' % (snp_vcf, snp_filted_vcf)
    execute_system(snv_filter_cmd, '[ Msg: SNV filtered done ! ]', '[ E: Something wrong with filter SNV ! ]')
    indel_filter_cmd = scriptPath + '/bin/gatk/gatk VariantFiltration ' \
                                    '-filter "QD < 2.0" --filter-name "QD2"  ' \
                                    '-filter "QUAL < 30.0" --filter-name "QUAL30" ' \
                                    '-filter "FS > 200.0" --filter-name "FS200" ' \
                                    '-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" ' \
                                    '-V %s -O %s' % (indelVcf, indel_filted_vcf)
    execute_system(indel_filter_cmd, '[ Msg: INDEL filtered done ! ]', '[ E: Something wrong with filter INDEL ! ]')
    merge_cmd = scriptPath + '/bin/gatk/gatk MergeVcfs -I %s -I %s -O %s' % (snp_filted_vcf, indel_filted_vcf, vcf)
    execute_system(merge_cmd, '[ Msg: merge filtered vcf done ! ]', '[ E: Something wrong with merge filtered vcf ! ]')

    rm_cmd = 'rm -f %s* %s* %s* %s* ' % (snp_vcf, snp_filted_vcf, indelVcf, indel_filted_vcf)
    execute_system(rm_cmd, '[ Msg: Delete process file done in gatk! ]',
                   '[ E: Something wrong with delete process file in gatk ! ]')
    # vcf stats
    vcf_stats(vcf, report_dir, reference, scriptPath)


def bcftools(infile, outdir, tmp_dir, report_dir, reference, sample_name, call_thread=12, zip_thread=12):
    outdir = outdir.rstrip('/') + '/'
    report_dir = report_dir.rstrip('/') + '/'
    tmp_dir = tmp_dir + sample_name + '/bcftools/'
    os.makedirs(tmp_dir)
    _tmp_vcf1 = tmp_dir + sample_name + '.bcftools.raw.vcf'
    _tmp_vcf2 = tmp_dir + sample_name + '.bcftools.flt.vcf'
    _ziped_vcf = outdir + sample_name + '.bcftools.vcf.gz'
    # call snvs/indels
    call_cmd = 'bcftools mpileup --threads {} -q 20 -Q 20 -Ou -f {} {} | bcftools call --threads {} -mv -Ov | ' \
               'bcftools filter --threads {} -s FILTER -g 10 -G 10 -i "%QUAL>20 && DP>6 && MQ>=40 && ' \
               '(DP4[2]+DP4[3])>4" > {}'.format(call_thread, reference, infile, call_thread, call_thread, _tmp_vcf1)
    execute_system(call_cmd, '[ Msg: <%s> call snvs/indels done by bcftools ! ]' % sample_name,
                   '[ E: Something wrong with <%s> call snvs/indels by bcftools ! ]' % sample_name)
    # filter variations
    filt_cmd = r'''awk -F "\t" '{if($1~/#/){print}else if($7~/PASS/){print}}' %s > %s''' % (_tmp_vcf1, _tmp_vcf2)
    execute_system(filt_cmd, '[ Msg: <%s> filter snvs/indels done in bcftools ! ]' % sample_name,
                   '[ E: Something wrong with filter <%s> raw variations in bcftools ! ]' % sample_name)
    # 合成前准备，包括合成列表文件，压缩，索引
    # 压缩
    zip_cmd = 'bgzip -c -f -@ %d %s > %s' % (zip_thread, _tmp_vcf2, _ziped_vcf)
    execute_system(zip_cmd, '[ Msg: Bgzip <%s> vcf done in bcfrools ! ]' % sample_name,
                   '[ E: Something wrong with bgzip <%s> vcf in bcftools ! ]' % sample_name)
    # 建立索引
    index_cmd = 'bcftools index -t %s' % _ziped_vcf
    execute_system(index_cmd, '[ Msg: Build <%s> vcf file index done in bcftools ! ]' % sample_name,
                   '[ E: Something wrong with build <%s> vcf index file in bcftools ! ]' % sample_name)
    # 删除过程文件
    rm_cmd = 'rm -f %s %s' % (_tmp_vcf1, _tmp_vcf2)
    execute_system(rm_cmd, '[ Msg: Delete <%s> process file done in bcftools ! ]' % sample_name,
                   '[ E: Something wrong with delete <%s> process file in bcftools ! ]' % sample_name)
    # vcf stats
    vcf_stats(_ziped_vcf, report_dir, sample_name, reference)
    return _ziped_vcf


def vardict(infile, outdir, tmp_dir, report_dir, reference, sample_name, filter_freq=0.1, call_thread=12,
            zip_thread=12):
    outdir = outdir.rstrip('/') + '/'
    report_dir = report_dir.rstrip('/') + '/'
    tmp_dir = tmp_dir + sample_name + '/vardict/'
    os.makedirs(tmp_dir)
    _tmp_vcf1 = tmp_dir + sample_name + '.vardict.raw.vcf'
    _tmp_vcf2 = tmp_dir + sample_name + '.vardict.flt.vcf'
    _ziped_vcf = outdir + sample_name + '.vardict.vcf.gz'
    # call snvs/indels
    call_cmd = './bin/VarDict-1.7.0/bin/VarDict -G %s -b %s -f %.2f -N %s -z -c 1 -S 2 -E 3 -g 4 -th %d ' \
               './lib/VarDict_assembly19_fromBroad_5k_150bpOL_seg.bed | teststrandbias.R | var2vcf_valid.pl -N %s -E ' \
               '-f %.2f > %s' % (reference, infile, filter_freq, sample_name, call_thread, sample_name, filter_freq,
                                 _tmp_vcf1)
    execute_system(call_cmd, '[ Msg: <%s> call snvs/indels done by VarDict ! ]' % sample_name,
                   '[ E: Something wrong with <%s> call snvs/indels by VarDict ! ]' % sample_name)
    # filter variations
    filt_cmd = r'''perl -e 'map{chomp; if($_=~/^#/){print $_."\n";} elsif($_ =~/<dup-/){} else''' \
               r'''{$_=~/(.*)TYPE=(\w+);/;if(($2 eq "SNV") || ($2 eq "Insertion") || ($2 eq "Deletion"))''' \
               r'''{ print $_."\n"}}}`cat %s`' | java -jar ./bin/snpEff/SnpSift.jar filter "(QUAL >= 20) ''' \
               r'''& (DP > 6) & (VD > 4) & (MQ >= 40) & ((FILTER='PASS')|(FILTER='InDelLikely'))" > %s''' % (
                   _tmp_vcf1, _tmp_vcf2)
    execute_system(filt_cmd, '[ Msg: <%s> filter snvs/indels done in VarDict ! ]' % sample_name,
                   '[ E: Something wrong with filter <%s> raw variations in VarDict ! ]' % sample_name)
    # 合成前准备，包括合成列表文件，压缩，索引
    # 压缩
    zip_cmd = 'bgzip -c -f -@ %d %s > %s' % (zip_thread, _tmp_vcf2, _ziped_vcf)
    execute_system(zip_cmd, '[ Msg: Bgzip <%s> vcf done in VarDict ! ]' % sample_name,
                   '[ E: Something wrong with bgzip <%s> vcf in VarDict ! ]' % sample_name)
    # 建立索引
    index_cmd = 'bcftools index -t %s' % _ziped_vcf
    execute_system(index_cmd, '[ Msg: Build <%s> vcf file index done in VarDict ! ]' % sample_name,
                   '[ E: Something wrong with build <%s> vcf index file in VarDict ! ]' % sample_name)
    # 删除过程文件
    rm_cmd = 'rm -f %s %s' % (_tmp_vcf1, _tmp_vcf2)
    execute_system(rm_cmd, '[ Msg: Delete <%s> process file done in VarDict ! ]' % sample_name,
                   '[ E: Something wrong with delete <%s> process file in VarDict ! ]' % sample_name)
    # vcf stats
    vcf_stats(_ziped_vcf, report_dir, sample_name, reference)
    return _ziped_vcf


def deepVariant_pre(infile, outdir, reference, sample_name, call_thread=12):
    indir = os.path.dirname(infile)
    outdir = os.path.abspath(outdir)
    ref_dir = os.path.dirname(reference)
    ref_name = os.path.basename(reference)
    bam_name = os.path.basename(infile)
    g_vcf = outdir + '/' + sample_name + '.deep.g.vcf.gz'
    _vcf = outdir + '/' + sample_name + '.deep.vcf.gz'
    call_cmd = 'docker run -v "%s":"/input" -v "%s":"/output" -v "%s":"/ref" google/deepvariant:"0.9.0" ' \
               '/opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=/ref/%s --reads=/input/%s ' \
               '--output_vcf=/output/%s --output_gvcf=/output/%s --num_shards=%d' % (
                   indir, outdir, ref_dir, ref_name, bam_name, _vcf, g_vcf, call_thread)
    execute_system(call_cmd, '[ Msg: <%s> deepvariant calling done ! ]' % sample_name,
                   '[ E: Something wrong with <%s> deepvariant calling ! ]' % sample_name)


def deepVariant(gvcf_list, outdir, report_dir, reference):
    outdir = os.path.abspath(outdir)
    merged_vcf = outdir + '/cohort.merged.vcf.gz'
    filter_vcf = outdir + '/cohort.filter.vcf.gz'
    colist = gvcf_list[0].split('/')
    for gvcf in gvcf_list[1:]:
        for num, i in enumerate(gvcf.split('/')):
            if len(colist) > num and i != colist[num]:
                colist = colist[:num]
    indir = '/'.join(colist)
    gvcf_str = ''
    for i in gvcf_list:
        gvcf_str += '/data%s ' % i[len(indir):]
    merge_cmd = 'docker run -v "%s":"/data" -v "%s":"/out" quay.io/mlin/glnexus:v1.2.2 ' \
                '/usr/local/bin/glnexus_cli --config DeepVariantWGS ' \
                '%s| bcftools view - | bgzip -c > /out/%s' % (indir, outdir, gvcf_str, merged_vcf)
    execute_system(merge_cmd, '[ Msg: DeepVariant combine gvcfs done ! ]',
                   '[ E: Something wrong with combine gvcfs in deepVariant ! ]')
    # filter
    filter_cmd = 'java -jar ./bin/snpEff/SnpSift.jar filter "(QUAL >= 25) & (GEN[0].DP > 3) & (GEN[0].DP < 300) & ' \
                 '(GEN[0].GQ > 20) & (GEN[0].VAF > 0.2) & (FILTER=\'PASS\')" %s > %s' % (merged_vcf, filter_vcf)
    execute_system(filter_cmd, '[ Msg: Filter snvs/indels done in DeepVariant ! ]',
                   '[ E: Something wrong with filter raw variations in DeepVariant ! ]')
    # vcf stats
    # vcf_stats(merged_vcf, report_dir, reference)
    return merged_vcf


def vcf_stats(vcf, report_dir, reference, scriptPath):
    outFile = report_dir + '/variantQC.html'
    vcfQCCmd = 'java -jar %s/bin/DISCVRSeq.jar VariantQC -O %s -R %s -V %s' % (scriptPath, outFile, reference, vcf)
    rst = os.system(vcfQCCmd)
    if rst:
        print('[ E: fail to vcf QC ! ]')
    else:
        print('[ Msg: vcf QC done ! ]')


def local_db_():
    cmd = 'gatk GenomicsDBImport ' \
          '--genomicsdb-workspace-path my_database ' \
          '--batch-size 50' \
          '-L bed' \
          '-V .g.vcf.gz' \
          '--tmp-dir=/path/to/large/tmp' \
          '--reader-threads 5'
    cmd = 'gatk GenomicsDBImport' \
          '-V data/gvcfs/mother.g.vcf.gz' \
          '-V data/gvcfs/father.g.vcf.gz' \
          '-V data/gvcfs/son.g.vcf.gz' \
          '--genomicsdb-update-workspace-path my_database' \
          '--tmp-dir=/path/to/large/tmp'
    cmd = 'gatk GenotypeGVCFs' \
          '-R Homo_sapiens_assembly38.fasta' \
          '-V gendb://my_database' \
          '-O output.vcf.gz' \
          '--tmp-dir=/path/to/large/tmp'


def clinSV():
    pass
