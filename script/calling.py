import os
import sys
from script.common import execute_system


# 输入是bam文件，
# 输出是对应模式的vcf
# 所有单个软件输入bam文件，输出vcf文件
def gatk_pre(bam, out_dir, tmp_dir, reference, prefix, gatk_bundle_dir, script_path, bed, keep_tmp):
    table = '%s/%s.recal_data.table' % (out_dir, prefix)
    bqsr_bam = '%s/%s.BQSR.bam' % (tmp_dir, prefix)
    g_vcf = '%s/%s.gatk.g.vcf.gz' % (out_dir, prefix)
    # BQSR
    # BaseRecalibrator
    if os.path.basename(reference).find('hg19') != -1:
        db_list = ['1000G_phase1.indels.hg19.sites.vcf.gz',
                   'Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz',
                   'dbsnp_138.hg19.vcf.gz',
                   '1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz']
    elif os.path.basename(reference).find('v37') != -1:
        db_list = ['1000G_phase1.indels.b37.vcf.gz',
                   'Mills_and_1000G_gold_standard.indels.b37.vcf.gz',
                   'dbsnp_138.b37.vcf.gz',
                   '1000G_phase1.snps.high_confidence.b37.vcf.gz']
    elif os.path.basename(reference).find('38') != -1:
        db_list = ['Mills_and_1000G_gold_standard.indels.hg38.vcf.gz',
                   'dbsnp_146.hg38.vcf.gz',
                   '1000G_phase1.snps.high_confidence.hg38.vcf.gz']
    else:
        sys.exit('[ Error: Can not find bundle file with <%s>]' % os.path.basename(reference))
    gatk_db = ''
    for db in db_list:
        if os.path.isfile(gatk_bundle_dir + '/' + db):
            gatk_db += '--known-sites %s/%s ' % (gatk_bundle_dir, db)
        else:
            print(' [ Warn: Can not find %s .]' % db)
    if not gatk_db:
        sys.exit('[ Error: Can not find any database for gatk.]')

    br_cmd = script_path + '/bin/gatk/gatk BaseRecalibrator --tmp-dir %s -R %s -I %s %s -O %s ' % (
        tmp_dir, reference, bam, gatk_db, table)
    if bed:
        br_cmd += '-L %s' % bed
    execute_system(br_cmd, '[ Msg: <%s> baseRecalibrator done ! ]' % prefix,
                   '[ Error: Something wrong with <%s> BaseRecalibrator ! ]' % prefix)
    # ApplyBQSRCmd
    bqsr_cmd = script_path + '/bin/gatk/gatk ApplyBQSR --tmp-dir %s --bqsr-recal-file %s -R %s -I %s -O %s ' % (
        tmp_dir, table, reference, bam, bqsr_bam)
    if bed:
        bqsr_cmd += '-L %s' % bed
    execute_system(bqsr_cmd, '[ Msg: <%s> ApplyBQSR done ! ]' % prefix,
                   '[ Error: Something wrong with <%s> ApplyBQSR ! ]' % prefix)
    # index
    index_cmd = 'samtools index %s' % bqsr_bam
    execute_system(index_cmd, '[ Msg: <%s> BQSR bam index done ! ]' % prefix,
                   '[ Error: Something wrong with <%s> BQSR bam index ! ]' % prefix)
    # HaplotypeCaller
    hc_cmd = script_path + '/bin/gatk/gatk HaplotypeCaller --tmp-dir %s -ERC GVCF -R %s -I %s -O %s ' % (
        tmp_dir, reference, bam, g_vcf)
    if bed:
        hc_cmd += '-L %s' % bed
    execute_system(hc_cmd, '[ Msg: <%s> HaplotypeCaller done ! ]' % prefix,
                   '[ Error: Something wrong with <%s> HaplotypeCaller ! ]' % prefix)
    if not keep_tmp:
        rm_cmd = 'rm -f %s* ' % bqsr_bam
        execute_system(rm_cmd, '[ Msg: Delete process file done in gatk! ]',
                       '[ Error: Something wrong with delete process file in gatk ! ]')
    return g_vcf


def gatk(gvcf_list, out_dir, report_dir, reference, gatk_bundle_dir, script_path, prefix, bed, tmp_dir, keep_tmp):
    merged_gvcf = '%s/%scohort.g.vcf.gz' % (out_dir, prefix)
    cohort_vcf = '%s/%scohort.raw.vcf.gz' % (out_dir, prefix)
    excesshet = '%s/%scohort_excesshet.vcf.gz' % (tmp_dir, prefix)
    sitesonly = '%s/%scohort_sitesonly.vcf.gz' % (tmp_dir, prefix)
    indels_recal = '%s/%cohort_indels.recal' % (tmp_dir, prefix)
    indels_tranches = '%s/%scohort_indels.tranches' % (tmp_dir, prefix)
    indels_vcf = '%s/%sindel.recalibrated.vcf.gz' % (tmp_dir, prefix)
    snp_recal = '%s/%cohort_snps.recal' % (tmp_dir, prefix)
    snp_tranches = '%s/%scohort_snps.tranches' % (tmp_dir, prefix)
    final_vcf = '%s/%scohort.filter.vcf.gz' % (out_dir, prefix)
    # mergeGVCFs
    sample_gvcf = ''
    for i in gvcf_list:
        sample_gvcf += '-V %s ' % i
    merge_cmd = script_path + '/bin/gatk/gatk CombineGVCFs --tmp-dir %s -R %s %s -O %s ' % (
        tmp_dir, reference, sample_gvcf, merged_gvcf)
    if bed:
        merge_cmd += '-L %s' % bed
    execute_system(merge_cmd, '[ Msg: GATK combine gvcfs done ! ]',
                   '[ Error: Something wrong with combine gvcfs in GATK ! ]')
    # GenotypeGVCFs
    gt_cmd = script_path + '/bin/gatk/gatk GenotypeGVCFs --tmp-dir %s -R %s -V %s -O %s' % (
        tmp_dir, reference, merged_gvcf, cohort_vcf)
    execute_system(gt_cmd, '[ Msg: Genotype GVCFs done ! ]',
                   '[ Error: Something wrong with genotype GVCFs ! ]')
    # VariantFiltration
    filtration_cmd = script_path + '/bin/gatk/gatk VariantFiltration ' \
                                   '--tmp-dir %s -V %s --filter-expression "ExcessHet > 54.69" ' \
                                   '--filter-name ExcessHet -O %s' % (tmp_dir, cohort_vcf, excesshet)
    execute_system(filtration_cmd, '[ Msg: Cohort VCF variantFiltration done ! ]',
                   '[ Error: Something wrong with cohort VCF variantFiltration ! ]')
    # MakeSitesOnly
    sites_only_cmd = script_path + '/bin/gatk/gatk MakeSitesOnlyVcf --tmp-dir %s -I %s -O %s' % (
        tmp_dir, excesshet, sitesonly)
    execute_system(sites_only_cmd, '[ Msg: Make sites only vcf done ! ]',
                   '[ Error: Something wrong with make sites only vcf ! ]')
    # VariantRecalibrator
    if os.path.basename(reference).find('hg19') != -1:
        _mill = '-resource:mills,known=false,training=true,truth=true,prior=12 ' \
                '%s/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz ' % gatk_bundle_dir
        _dbsnp = '-resource:dbsnp,known=true,training=false,truth=false,prior=2 ' \
                 '%s/dbsnp_138.hg19.vcf.gz' % gatk_bundle_dir
        _dbsnp2 = '%s/dbsnp_138.hg19.vcf.gz' % gatk_bundle_dir
        _hapmap = '%s/hapmap_3.3.hg19.sites.vcf.gz' % gatk_bundle_dir
        _omni = '%s/1000G_omni2.5.hg19.sites.vcf.gz' % gatk_bundle_dir
        _1kg = '%s/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz' % gatk_bundle_dir
        indel_db = _mill + _dbsnp
    elif os.path.basename(reference).find('v37') != -1:
        _mill = '-resource:mills,known=false,training=true,truth=true,prior=12 ' \
                '%s/Mills_and_1000G_gold_standard.indels.b37.vcf.gz' % gatk_bundle_dir
        _dbsnp = '-resource:dbsnp,known=true,training=false,truth=false,prior=2 ' \
                 '%s/dbsnp_138.b37.vcf.gz' % gatk_bundle_dir
        _dbsnp2 = '%s/dbsnp_138.b37.vcf.gz' % gatk_bundle_dir
        _hapmap = '%s/hapmap_3.3.b37.vcf.gz' % gatk_bundle_dir
        _omni = '%s/1000G_omni2.5.b37.vcf.gz' % gatk_bundle_dir
        _1kg = '%s/1000G_phase1.snps.high_confidence.b37.vcf.gz' % gatk_bundle_dir
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
        sys.exit('[ Error: Can not find bundle file with <%s>]' % os.path.basename(reference))
    recalibrator_cmd = script_path + '/bin/gatk/gatk VariantRecalibrator -V %s --trust-all-polymorphic ' \
                                     '-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 ' \
                                     '-tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 ' \
                                     '-tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 ' \
                                     '-an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP -mode INDEL ' \
                                     '--max-gaussians 4 %s -O %s --tranches-file %s' % (
                           sitesonly, indel_db, indels_recal, indels_tranches)
    execute_system(recalibrator_cmd, '[ Msg: Indel recalibrator done ! ]',
                   '[ Error: Something wrong with indel recalibrator ! ]')
    apply_vqsr_cmd = script_path + '/bin/gatk/gatk ApplyVQSR -V %s --recal-file %s --tranches-file %s ' \
                                   '--truth-sensitivity-filter-level 99.7 --create-output-variant-index true ' \
                                   '-mode INDEL -O %s' % (excesshet, indels_recal, indels_tranches, indels_vcf)
    execute_system(apply_vqsr_cmd, '[ Msg: Indel apply VQSR done ! ]',
                   '[ Error: Something wrong with indel apply VQSR ! ]')
    recalibrator_cmd = script_path + '/bin/gatk/gatk VariantRecalibrator -V %s --trust-all-polymorphic ' \
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
                   '[ Error: Something wrong with snp recalibrator ! ]')
    filter_cmd = script_path + '/bin/gatk/gatk ApplyVQSR -V %s --recal-file %s --tranches-file %s ' \
                               '--truth-sensitivity-filter-level 99.7 ' \
                               '--create-output-variant-index true -mode SNP -O %s' % (
                     indels_vcf, snp_recal, snp_tranches, final_vcf)
    execute_system(filter_cmd, '[ Msg: snps recalibrated done ! ]',
                   '[ Error: Something wrong with snps recalibrated ! ]')
    # 建立索引
    # index_cmd = 'bcftools index -t %s' % final_vcf
    # execute_system(index_cmd, '[ Msg: Build vcf index done in gatk ! ]',
    #                '[ E: Something wrong with build vcf index file in gatk ! ]')
    # 删除无用中间文件
    if not keep_tmp:
        rm_cmd = 'rm -f %s %s %s %s %s %s %s %s' % (
            sitesonly, excesshet, merged_gvcf, indels_tranches, indels_recal, indels_vcf, snp_recal,
            snp_tranches)
        execute_system(rm_cmd, '[ Msg: Delete process file done in gatk! ]',
                       '[ Error: Something wrong with delete process file in gatk ! ]')
    # vcf stats
    vcf_stats(final_vcf, report_dir, reference, script_path)
    return final_vcf


def gatk_hard_filter(g_vcf, out_dir, report_dir, tmp_dir, prefix, reference, script_path, bed, keep_tmp):
    raw_vcf = '%s/%s.gatk.raw.vcf.gz' % (out_dir, prefix)
    snv_vcf = '%s/%s.gatk.snvs.vcf.gz' % (tmp_dir, prefix)
    snv_filted_vcf = '%s/%s.gatk.snvs.filtered.vcf.gz' % (tmp_dir, prefix)
    indel_vcf = '%s/%s.gatk.indels.vcf.gz' % (tmp_dir, prefix)
    indel_filted_vcf = '%s/%s.gatk.indels.filtered.vcf.gz' % (tmp_dir, prefix)
    vcf = '%s/%s.gatk.final.vcf.gz' % (out_dir, prefix)
    # GT
    gt_cmd = script_path + '/bin/gatk/gatk GenotypeGVCFs -R %s -V %s -O %s ' % (reference, g_vcf, raw_vcf)
    if bed:
        gt_cmd += '-L %s' % bed
    execute_system(gt_cmd, '[ Msg: Genotype GVCF done ! ]', '[ Error: Something wrong with genotype GVCF ! ]')
    # select
    select_snv_cmd = script_path + '/bin/gatk/gatk SelectVariants -V %s -select-type SNP -O %s' % (raw_vcf, snv_vcf)
    execute_system(select_snv_cmd, '[ Msg: select SNV done ! ]', '[ Error: Something wrong with select SNV ! ]')
    select_indel_cmd = script_path + '/bin/gatk/gatk SelectVariants -V %s -select-type INDEL -O %s' % (
        raw_vcf, indel_vcf)
    execute_system(select_indel_cmd, '[ Msg: select INDEL done ! ]', '[ Error: Something wrong with select INDEL ! ]')
    # filter
    snv_filter_cmd = script_path + '/bin/gatk/gatk VariantFiltration ' \
                                   '-filter "QD < 2.0" --filter-name "QD2" ' \
                                   '-filter "QUAL < 30.0" --filter-name "QUAL30" ' \
                                   '-filter "SOR > 3.0" --filter-name "SOR3" ' \
                                   '-filter "FS > 60.0" --filter-name "FS60" ' \
                                   '-filter "MQ < 40.0" --filter-name "MQ40" ' \
                                   '-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" ' \
                                   '-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" ' \
                                   '-V %s -O %s' % (snv_vcf, snv_filted_vcf)
    execute_system(snv_filter_cmd, '[ Msg: SNV filtered done ! ]', '[ Error: Something wrong with filter SNV ! ]')
    indel_filter_cmd = script_path + '/bin/gatk/gatk VariantFiltration ' \
                                     '-filter "QD < 2.0" --filter-name "QD2"  ' \
                                     '-filter "QUAL < 30.0" --filter-name "QUAL30" ' \
                                     '-filter "FS > 200.0" --filter-name "FS200" ' \
                                     '-filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" ' \
                                     '-V %s -O %s' % (indel_vcf, indel_filted_vcf)
    execute_system(indel_filter_cmd, '[ Msg: INDEL filtered done ! ]', '[ Error: Something wrong with filter INDEL ! ]')
    merge_cmd = script_path + '/bin/gatk/gatk MergeVcfs -I %s -I %s -O %s' % (snv_filted_vcf, indel_filted_vcf, vcf)
    execute_system(merge_cmd, '[ Msg: merge filtered vcf done ! ]',
                   '[ Error: Something wrong with merge filtered vcf ! ]')
    if not keep_tmp:
        rm_cmd = 'rm -f %s* %s* %s* %s* ' % (snv_vcf, snv_filted_vcf, indel_vcf, indel_filted_vcf)
        execute_system(rm_cmd, '[ Msg: Delete process file done in gatk! ]',
                       '[ Error: Something wrong with delete process file in gatk ! ]')
    # vcf stats
    vcf_stats(vcf, report_dir, reference, script_path)
    return vcf


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


def vcf_stats(vcf, report_dir, reference, script_path):
    out_file = report_dir + '/variantQC.html'
    vcfQC_cmd = 'java -jar %s/bin/DISCVRSeq.jar VariantQC -O %s -R %s -V %s' % (script_path, out_file, reference, vcf)
    rst = os.system(vcfQC_cmd)
    if rst:
        print('[ Error: fail to vcf QC ! ]')
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
