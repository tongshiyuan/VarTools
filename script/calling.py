import os
import sys
from script.common import execute_system, check_software


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


def gatk(gvcf_list, out_dir, report_dir, reference, gatk_bundle_dir, script_path, prefix, bed, tmp_dir, keep_tmp,
         nqc=False, flt=True):
    merged_gvcf = '%s/%scohort.g.vcf.gz' % (out_dir, prefix)
    cohort_vcf = '%s/%scohort.raw.vcf.gz' % (out_dir, prefix)
    excesshet = '%s/%scohort_excesshet.vcf.gz' % (tmp_dir, prefix)
    sitesonly = '%s/%scohort_sitesonly.vcf.gz' % (tmp_dir, prefix)
    indels_recal = '%s/%scohort_indels.recal' % (tmp_dir, prefix)
    indels_tranches = '%s/%scohort_indels.tranches' % (tmp_dir, prefix)
    indels_vcf = '%s/%sindel.recalibrated.vcf.gz' % (tmp_dir, prefix)
    snp_recal = '%s/%scohort_snps.recal' % (tmp_dir, prefix)
    snp_tranches = '%s/%scohort_snps.tranches' % (tmp_dir, prefix)
    final_vcf = '%s/%scohort.filter.vcf.gz' % (out_dir, prefix)
    # mergeGVCFs
    sample_gvcf = ''
    if len(gvcf_list) == 1:
        merged_gvcf = gvcf_list[0]
    else:
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
    if flt:
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
                                         '-an QD -an MQRankSum -an ReadPosRankSum ' \
                                         '-an FS -an MQ -an SOR -an DP -mode SNP ' \
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
    else:
        final_vcf = cohort_vcf
    # vcf stats
    if not nqc:
        vcf_stats(final_vcf, report_dir, reference, script_path)
    return final_vcf


def gatk_hard_filter(g_vcf, out_dir, report_dir, tmp_dir, prefix, reference, script_path, bed, keep_tmp, nqc=False,
                     flt=True):
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
    if flt:
        # select
        select_snv_cmd = script_path + '/bin/gatk/gatk SelectVariants -V %s -select-type SNP -O %s' % (raw_vcf, snv_vcf)
        execute_system(select_snv_cmd, '[ Msg: select SNV done ! ]', '[ Error: Something wrong with select SNV ! ]')
        select_indel_cmd = script_path + '/bin/gatk/gatk SelectVariants -V %s -select-type INDEL -O %s' % (
            raw_vcf, indel_vcf)
        execute_system(select_indel_cmd, '[ Msg: select INDEL done ! ]',
                       '[ Error: Something wrong with select INDEL ! ]')
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
        execute_system(indel_filter_cmd, '[ Msg: INDEL filtered done ! ]',
                       '[ Error: Something wrong with filter INDEL ! ]')
        merge_cmd = script_path + '/bin/gatk/gatk MergeVcfs -I %s -I %s -O %s' % (snv_filted_vcf, indel_filted_vcf, vcf)
        execute_system(merge_cmd, '[ Msg: merge filtered vcf done ! ]',
                       '[ Error: Something wrong with merge filtered vcf ! ]')
        if not keep_tmp:
            rm_cmd = 'rm -f %s* %s* %s* %s* ' % (snv_vcf, snv_filted_vcf, indel_vcf, indel_filted_vcf)
            execute_system(rm_cmd, '[ Msg: Delete process file done in gatk! ]',
                           '[ Error: Something wrong with delete process file in gatk ! ]')
    else:
        vcf = raw_vcf
    # vcf stats
    if not nqc:
        vcf_stats(vcf, report_dir, reference, script_path)
    return vcf


def strelka(bam, out_dir, report_dir, reference, script_path, thread, bed, tmp_dir, nqc=False, flt=True):
    rst1 = check_software('bgzip')
    rst2 = check_software('tabix')
    if rst1 or rst2:
        sys.exit('[ Error: Can not open <bgzip> or <tabix>.]')
    out_dir += '/strelka_workplace'
    if bed:
        compress_cmd = 'cp %s %s && bgzip %s/%s && tabix -b 2 -e 3 -p bed %s/%s.gz' % (
            bed, tmp_dir, tmp_dir, bed, tmp_dir, bed)
        execute_system(compress_cmd, '[ Msg: Make bed for strelka ! ]',
                       '[ Error: Something wrong with Make bed for strelka ! ]')
        bed_cmd = '--callRegions %s/%s.gz' % (tmp_dir, bed)
    else:
        bed_cmd = ''
    cf_cmd = '%s/bin/strelka2/bin/configureStrelkaGermlineWorkflow.py --bam %s --referenceFasta %s --runDir %s %s' % (
        script_path, bam, reference, out_dir, bed_cmd)
    execute_system(cf_cmd, '[ Msg: Configuration for strelka done ! ]',
                   '[ Error: Something wrong with strelka configuration ! ]')
    call_cmd = '%s/runWorkflow.py -m local -j %d' % (out_dir, thread)
    execute_system(call_cmd, '[ Msg: Call variants with strelka done ! ]',
                   '[ Error: Something wrong with strelka call variants ! ]')
    final_vcf = out_dir + '/results/variants/variants.vcf.gz'
    if flt:
        pass
    if not nqc:
        vcf_stats(final_vcf, report_dir, reference, script_path)
    return final_vcf


def vardict(bam, out_dir, reference, bed, report_dir, tmp_dir, script_path, prefix, thread, filter_freq, nqc=False,
            flt=True):
    rst1 = check_software('bgzip')
    rst2 = check_software('tabix')
    if rst1 or rst2:
        sys.exit('[ Error: Can not open <bgzip> or <tabix>.]')
    if not bed:
        if reference.find('38') != -1:
            bed = script_path + '/lib/VarDict_assembly38_fromBroad_5k_150bpOL_seg.bed'
        else:
            bed = script_path + '/lib/VarDict_assembly19_fromBroad_5k_150bpOL_seg.bed'
    raw_vcf = '%s/%s.vardict.raw.vcf' % (tmp_dir, prefix)
    sort_var_vcf = '%s/%s.vardict.short_var.vcf' % (tmp_dir, prefix)
    final_vcf = '%s/%s.vardict.flt.vcf.gz' % (out_dir, prefix)
    # call snvs/indels
    call_cmd = '%s/bin/vardict/bin/VarDict -G %s -b %s -f %.2f -N %s -c 1 -S 2 -E 3 -g 4 -th %d %s | Rscript ' \
               '%s/bin/vardict/bin/teststrandbias.R | %s/bin/vardict/bin/var2vcf_valid.pl -N %s -E -f %.2f > %s' % (
                   script_path, reference, bam, filter_freq, prefix, thread, bed, script_path, script_path, prefix,
                   filter_freq, raw_vcf)
    execute_system(call_cmd, '[ Msg: <%s> call snvs/indels done by VarDict ! ]' % prefix,
                   '[ Error: Something wrong with <%s> call snvs/indels by VarDict ! ]' % prefix)
    filter_cmd = r'''perl -e 'map{chomp; if($_=~/^#/){print $_."\n";} elsif($_ =~/<dup-/){} else''' \
                 r'''{$_=~/(.*)TYPE=(\w+);/;if(($2 eq "SNV") || ($2 eq "Insertion") || ($2 eq "Deletion"))''' \
                 r'''{ print $_."\n"}}}`cat %s`' > %s''' % (raw_vcf, sort_var_vcf)
    execute_system(filter_cmd, '[ Msg: <%s> filter snvs/indels done in VarDict ! ]' % prefix,
                   '[ Error : Something wrong with filter <%s> raw variations in VarDict ! ]' % prefix)
    if flt:
        # filter variations
        # filter_cmd = r'''perl -e 'map{chomp; if($_=~/^#/){print $_."\n";} elsif($_ =~/<dup-/){} else''' \
        #              r'''{$_=~/(.*)TYPE=(\w+);/;if(($2 eq "SNV") || ($2 eq "Insertion") || ($2 eq "Deletion"))''' \
        #              r'''{ print $_."\n"}}}`cat %s`' | java -jar %s/bin/snpEff/SnpSift.jar filter "(QUAL >= 20) ''' \
        #              r'''& (DP > 6) & (VD > 4) & (MQ >= 40) & ((FILTER='PASS')|(FILTER='InDelLikely'))" | ''' \
        #              r'''bgzip -c -f -@ %d > %s''' % (raw_vcf, script_path, thread, final_vcf)
        filter_cmd = r'''java -jar %s/bin/snpEff/SnpSift.jar filter "(QUAL >= 20) ''' \
                     r'''& (DP > 6) & (VD > 4) & (MQ >= 40) & ((FILTER='PASS')|(FILTER='InDelLikely'))" %s | ''' \
                     r'''bgzip -c -f -@ %d > %s''' % (script_path, sort_var_vcf, thread, final_vcf)
        execute_system(filter_cmd, '[ Msg: <%s> filter snvs/indels done in VarDict ! ]' % prefix,
                       '[ Error : Something wrong with filter <%s> raw variations in VarDict ! ]' % prefix)
    else:
        final_vcf = sort_var_vcf
    # 建立索引
    index_cmd = 'tabix %s' % final_vcf
    execute_system(index_cmd, '[ Msg: Build <%s> vcf file index done in VarDict ! ]' % prefix,
                   '[ Error: Something wrong with build <%s> vcf index file in VarDict ! ]' % prefix)
    # vcf stats
    if not nqc:
        vcf_stats(final_vcf, report_dir, reference, script_path)
    return final_vcf


def bcftools(bam, out_dir, tmp_dir, report_dir, reference, prefix, thread, script_path, bed, nqc=False, flt=True):
    rst = check_software('bcftools')
    rst1 = check_software('bgzip')
    rst2 = check_software('tabix')
    if rst or rst1 or rst2:
        sys.exit('[ Error: Can not open <bgzip> or <tabix> or <bcftools>.]')
    if bed:
        bed_cmd = '-T %s' % bed
    else:
        bed_cmd = ''
    raw_vcf = '%s/%s.bcftools.raw.vcf.gz' % (tmp_dir, prefix)
    final_vcf = '%s/%s.bcftools.flt.vcf.gz' % (out_dir, prefix)
    # call snvs/indels
    call_cmd = 'bcftools mpileup --threads {} {} -Ou -f {} {} | bcftools call --threads {} -mv -Oz -o {}'.format(
        thread, bed_cmd, reference, bam, thread, raw_vcf)
    execute_system(call_cmd, '[ Msg: <%s> call snvs/indels done by bcftools ! ]' % prefix,
                   '[ Error: Something wrong with <%s> call snvs/indels by bcftools ! ]' % prefix)
    if flt:
        # filter
        filter_cmd = 'bcftools filter -s FILTER -g 10 -G 10 -i "%QUAL>20 && DP>6 && MQ>=40 && (DP4[2]+DP4[3])>4" ' + \
                     '--threads %d -Ov %s | awk -F"\t" \'{if($1~/#/){print}else if($7~/PASS/){print}}\' | ' \
                     'bgzip -c -f -@ %d > %s' % (thread, raw_vcf, thread, final_vcf)
        execute_system(filter_cmd, '[ Msg: <%s> filter snvs/indels done in bcftools ! ]' % prefix,
                       '[ Error: Something wrong with filter <%s> raw variations in bcftools ! ]' % prefix)
    else:
        final_vcf = raw_vcf
    # 建立索引
    index_cmd = 'tabix %s' % final_vcf
    execute_system(index_cmd, '[ Msg: Build <%s> vcf file index done in bcftools ! ]' % prefix,
                   '[ Error: Something wrong with build <%s> vcf index file in bcftools ! ]' % prefix)
    # vcf stats
    if not nqc:
        vcf_stats(final_vcf, report_dir, reference, script_path)
    return final_vcf


def deep_variant_single(bam, out_dir, reference, bed, report_dir, script_path, prefix, thread, nqc=False, flt=True,
                        version="1.2.0"):
    rst1 = check_software('bgzip')
    rst2 = check_software('tabix')
    if rst1 or rst2:
        sys.exit('[ Error: Can not open <bgzip> or <tabix>.]')
    in_dir = os.path.dirname(bam)
    out_dir = os.path.abspath(out_dir)
    ref_dir = os.path.dirname(reference)
    ref_name = os.path.basename(reference)
    bam_name = os.path.basename(bam)
    raw_vcf = prefix + '.deep.raw.vcf.gz'
    g_vcf = prefix + '.deep.g.vcf.gz'
    os.makedirs(out_dir + '/intermediate_results_dir')
    final_vcf = '%s/%s.deep.flt.vcf.gz' % (out_dir, prefix)
    if bed:
        mt = 'WES'
        cp_bed_cmd = 'cp %s %s/%s' % (bed, in_dir, bed)
        execute_system(cp_bed_cmd, '[Msg: Copy bed for deep variants done! ]',
                       '[Error: Something wrong with copy bed for deepvariants]')
        region_cmd = '--regions /input/%s' % bed
    else:
        mt = 'WGS'
        region_cmd = ''
    rst = check_software('docker')
    if rst:
        print('[ Warn: Can not open <docker>. Try to use Singularity.]')
        rst = check_software('singularity')
        if rst:
            sys.exit('[ Error: Can not open <Singularity>.]')
        else:
            if not os.path.isfile('%s/bin/deepvariant_latest.sif' % script_path):
                pull_cmd = 'singularity pull docker:google/deepvariant && ' \
                           'mv deepvariant_latest.sif %s/bin/' % script_path
                execute_system(pull_cmd, '[ Msg: pull deepvariant image done ! ]' % prefix,
                               '[ Error: Can not pull deepvariant image ! ]' % prefix)
            # singularity call
            call_cmd = 'singularity exec --bind %s:/input,%s:/output,%s:/ref %s/bin/deepvariant_latest.sif ' \
                       '/opt/deepvariant/bin/run_deepvariant -model_type %s ' \
                       '--ref /ref/%s ' \
                       '--reads /input/%s ' \
                       '--output_vcf /output/%s ' \
                       '--output_gvcf /output/%s ' \
                       '--intermediate_results_dir /output/intermediate_results_dir ' \
                       '--num_shards %d %s' % (
                           in_dir, out_dir, ref_dir, script_path, mt, ref_name, bam_name, raw_vcf, g_vcf, thread,
                           region_cmd)
    else:
        # docker call
        call_cmd = 'docker run -v "%s":"/input" -v "%s":"/output" -v "%s":"/ref" google/deepvariant:"%s" ' \
                   '/opt/deepvariant/bin/run_deepvariant --model_type %s --ref /ref/%s --reads /input/%s ' \
                   '--output_vcf /output/%s --output_gvcf /output/%s --num_shards %d ' \
                   '--intermediate_results_dir /output/intermediate_results_dir %s' % (
                       in_dir, out_dir, ref_dir, version, mt, ref_name, bam_name, raw_vcf, g_vcf, thread, region_cmd)
    execute_system(call_cmd, '[ Msg: <%s> deepvariant calling done ! ]' % prefix,
                   '[ Error: Something wrong with <%s> deepvariant calling ! ]' % prefix)
    if flt:
        # filter
        filter_cmd = 'java -jar %s/bin/snpEff/SnpSift.jar filter "(QUAL >= 20) & (GEN[0].DP > 6) & ' \
                     '(GEN[0].GQ > 20) & (GEN[0].VAF > 0.2) & (FILTER=\'PASS\')" %s/%s | bgzip -c -f -@ %d > %s' % (
                         script_path, out_dir, raw_vcf, thread, final_vcf)
        execute_system(filter_cmd, '[ Msg: Filter snvs/indels done in DeepVariant ! ]',
                       '[ Error: Something wrong with filter raw variations in DeepVariant ! ]')
    else:
        final_vcf = raw_vcf
    # 建立索引
    index_cmd = 'tabix %s' % final_vcf
    execute_system(index_cmd, '[ Msg: Build <%s> vcf file index done in VarDict ! ]' % prefix,
                   '[ Error: Something wrong with build <%s> vcf index file in VarDict ! ]' % prefix)
    # vcf stats
    if not nqc:
        vcf_stats(final_vcf, report_dir, reference, script_path)
    return final_vcf


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
