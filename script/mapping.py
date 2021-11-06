import os
from script.bam_qc import bam_stats
from script.common import execute_system, get_raw_info


def bam_deal(fq1, fq2, out_dir, report_dir, tmp_dir,
             sample_name, reference,
             software, thread, script_path,
             bq, ver, bed, fmd, rmd, fp_rmd, platform, keep_tmp, gender_rate):
    sorted_bam = align(fq1, fq2, tmp_dir, sample_name, reference, software, thread, script_path, platform, tmp_dir)
    # mark duplicate
    duped_bam = '%s/%s.sorted.merge.markdup.bam' % (out_dir, sample_name)
    if fp_rmd:
        # index
        index_cmd = 'mv %s %s && samtools index %s' % (sorted_bam, duped_bam, duped_bam)
        execute_system(index_cmd, '[ Msg: <%s> bam index done ! ]' % sample_name,
                       '[ Error: Something wrong with index <%s> bam ! ]' % sample_name)
    else:
        if fmd:
            if rmd:
                dup_cmd = script_path + '/bin/sambamba markdup -r -t %d --tmpdir %s %s %s' % (
                    thread, tmp_dir, sorted_bam, duped_bam)
            else:
                dup_cmd = script_path + '/bin/sambamba markdup -t %d --tmpdir %s %s %s' % (
                    thread, tmp_dir, sorted_bam, duped_bam)
            execute_system(dup_cmd, '[ Msg: <%s> mark duplication done with sambamba ! ]' % sample_name,
                           '[ Error: Something wrong with <%s> mark duplication with sambamba ! ]' % sample_name)
        else:
            dup_metrics = '%s/%s.markdup_metrics.txt' % (out_dir, sample_name)
            dup_cmd = script_path + '/bin/gatk/gatk MarkDuplicates ' \
                                    '-INPUT %s -OUTPUT %s -METRICS_FILE %s --TMP_DIR %s' % (
                          sorted_bam, duped_bam, dup_metrics, tmp_dir)
            execute_system(dup_cmd, '[ Msg: <%s> mark duplication done ! ]' % sample_name,
                           '[ Error: Something wrong with <%s> mark duplication ! ]' % sample_name)
            # index 使用sambamba会自动产生索引文件
            index_cmd = 'samtools index %s' % duped_bam
            execute_system(index_cmd, '[ Msg: <%s> bam index done ! ]' % sample_name,
                           '[ Error: Something wrong with index <%s> marked dup bam file ! ]' % sample_name)
    bam_stats(duped_bam, report_dir, thread, tmp_dir, script_path, bq, ver, bed, keep_tmp, gender_rate)
    return duped_bam


def align(fq1, fq2, out_dir, prefix, reference, software, thread, script_path, platform, tmp_dir):
    _lb, _id = get_raw_info(fq1)
    # 比对、转换、排序
    tmp_prefix = tmp_dir + '/' + prefix
    if software == 'bwa':
        out_bam = bwa_mem2(fq1, fq2, out_dir, prefix, reference, thread, _id, _lb, script_path, platform, tmp_prefix)
    elif software == 'bowtie2':
        out_bam = bowtie2(fq1, fq2, out_dir, prefix, reference, _id, _lb, script_path, platform, thread, tmp_prefix)
    elif software == 'gg':  # 由于暂止不支持hg38，暂时不进行维护使用
        out_bam = graph_genome(fq1, fq2, out_dir, prefix, reference, _id, _lb, thread)
    else:
        exit('[ Error: Do not support <%s> to mapping ! ]' % software)
        out_bam = ''
    return out_bam


def bwa_mem2(fq1, fq2, out_dir, prefix, reference, thread, _id, _lb, script_path, platform, tmp_prefix):
    # 检查索引文件
    out_bam = '%s/%s.sorted.bam' % (out_dir, prefix)
    suffix_list = ['.0123', '.amb', '.ann', '.bwt.2bit.64', '.bwt.8bit.32', '.pac']
    for suffix in suffix_list:
        if not os.path.exists(reference + suffix):
            print('[ Msg: Do not find bwa-mem2 index file of reference <%s> .]' % os.path.basename(reference))
            print('[ Msg: begin to build bwa-mem2 index file of reference ... ]')
            index_cmd = script_path + '/bin/bwa-mem2/bwa-mem2 index %s' % reference
            execute_system(index_cmd, '[ Msg: Build bwa-mem2 index file done ! ]',
                           '[ Error: Fail to build bwa-mem2 index file of reference ! ]')
            break
    # 开始比对
    map_cmd = script_path + r'/bin/bwa-mem2/bwa-mem2 mem -t %d -M -R ' \
                            r'"@RG\tID:%s\tPL:%s\tLB:%s\tSM:%s" %s %s %s | ' \
                            r'samtools sort -@ %d -T %s -m 4G - > %s' % (
                  thread, _id, platform, _lb, prefix, reference, fq1, fq2, thread, tmp_prefix, out_bam)
    execute_system(map_cmd, '[ Msg: <%s> mapping and sort done ! ]' % prefix,
                   '[ Error: <%s> fail to mapping with bwa-mem2 or sort bam file ! ]' % prefix)
    return out_bam


def bowtie2(fq1, fq2, out_dir, prefix, reference, _id, _lb, script_path, platform, thread, tmp_prefix):
    # 检查索引文件
    out_bam = '%s/%s.sorted.bam' % (out_dir, prefix)
    suffix_list = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']
    for suffix in suffix_list:
        if not os.path.exists(reference + suffix):
            print('[ Msg: Do not find bowtie2 index file of reference <%s> .]' % os.path.basename(reference))
            print('[ Msg: begin to build bowtie2 index file of reference ... ]')
            index_cmd = script_path + '/bin/bowtie2/bowtie2-build %s %s' % (reference, reference)
            execute_system(index_cmd, '[ Msg: Build bowtie2 index file done ! ]',
                           '[ Error: Fail to build bowtie2 index file of reference ! ]')
            break
    # 开始比对
    # grep -v "XS:" %s > %s && samtools view -q 1 -Shb %s
    map_cmd = script_path + '/bin/bowtie2/bowtie2 -p %d -x %s --no-unal --rg-id %s --rg PL:%s --rg LB:%s --rg SM:%s ' \
                            '-1 %s -2 %s | samtools sort -@ %d -T %s -m 4G - > %s' \
              % (thread, reference, _id, platform, _lb, prefix, fq1, fq2, thread, tmp_prefix, out_bam)
    execute_system(map_cmd, '[ Msg: <%s> mapping and sort done ! ]' % prefix,
                   '[ Error: <%s> fail to mapping with bowtie2 or sort bam file ! ]' % prefix)
    return out_bam


def graph_genome(fq1, fq2, out_dir, reference, prefix, _id, _lb, thread, gg_version='0.9.1'):
    out_bam = '%s/%s.sorted.bam' % (out_dir, prefix)
    # 挂载目录
    input_dir = os.path.dirname(os.path.abspath(fq1))
    reference_dir = os.path.dirname(os.path.abspath(reference))
    reference_basename = os.path.basename(reference)
    _fq1 = os.path.basename(fq1)
    _fq2 = os.path.basename(fq2)
    _bam = os.path.basename(out_bam)
    _tmp_bam = _bam + '.tmp0'
    map_cmd = 'docker run -v "%s":"/input" -v "%s":"/ref" -v "%s":"/out" gral-bpa:"%s" /usr/local/bin/aligner ' \
              '--vcf /ref/SBG.Graph.B37.V6.rc6.vcf.gz --reference /ref/%s -q /input/%s -Q /input/%s -o /out/%s ' \
              '--read_group_sample "%s" --read_group_library "%s" --read_group_id "%s" --threads %d && ' \
              'samtools sort -@ %d -m 4G -O bam -o %s %s' % (
                  input_dir, reference_dir, out_dir, gg_version, reference_basename, _fq1, _fq2, _tmp_bam, prefix, _lb,
                  _id, thread, thread, out_bam, out_dir + '/' + _tmp_bam)
    execute_system(map_cmd, '[ Msg: <%s> mapping and sort done ! ]' % prefix,
                   '[ Error: <%s> fail to mapping with graph genome or sort bam file ! ]' % prefix)
    # 删除中间文件
    rm_cmd = 'rm -f %s' % _tmp_bam
    execute_system(rm_cmd, '[Msg: Delete <%s> process file done after mapping and sort !]' % prefix,
                   '[ Error: Fail delete <%s> process file after mapping and sort ! ]' % prefix)
    return out_bam
