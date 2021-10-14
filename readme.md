# VarTools

## requests

- python 3.x
- numpy
- pandas
- scipy

## funcion

1. f2v: analysis from fastq to gvcf.
2. tGT: from gvcf created by GATK to vcf in trio mode.
3. sGT: from gvcf created by GATK to vcf in single mode.
4. cc: case-control analysis with 2 ways.
5. fp: build false positive database for filter.
6. bqc: bam quality check.
7. gd: gender identify.

### f2v

~~~(shell)
python3 VarTools.py f2v -i in_dir -o out_dir -b bed -p prefix --vcf --fastqc --qualimap --keep_tmp -t 6
~~~

### trio_genotype

~~~(shell)
python3 VarTools tGT -p .gvcf -f .gvcf -m .gvcf -s .gvcf,.gvcf -o ./result/outName
~~~

### single_genotype

~~~(shell)
python3 VarTools sGT -p .gvcf -o ./result/outName
~~~

### case-control

~~~(shell)
python3 VarTools cc --case case_dir --control control_dir --mode AD -o result
~~~

### build false positive databse

~~~(shell)
python3 VarTools fp -i indir -o ./result --snvdb false_positive.txt --overlap_rate 0.3 --file_type vcf
~~~

### bam quality check

~~~(shell)
python3 VarTools bqc -b in.bam --bed bed
~~~

### gender identify

~~~(shell)
python VarTools gd -b in.bam -d bed
~~~

### software

- [GATK](https://github.com/broadinstitute/gatk/releases/download/4.2.0.0/gatk-4.2.0.0.zip)
- [fastp](http://opengene.org/fastp/fastp)
- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip)
- [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2)
- [samtools](https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2)
- [qualimap](https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.2.1.zip)
- [mosdepth](https://github.com/brentp/mosdepth/releases/download/v0.3.2/mosdepth)
- [sambamba](https://github.com/biod/sambamba/releases/download/v0.8.1/sambamba-0.8.1-linux-amd64-static.gz)
- [plink1.9](https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20210606.zip)
- [kings](https://www.kingrelatedness.com/Linux-king.tar.gz)
- [VariantQC](https://github.com/BimberLab/DISCVRSeq/releases/download/1.3.2/DISCVRSeq-1.3.2.jar)
- [annovar](http://www.openbioinformatics.org/annovar)
- [clinSV](https://github.com/KCCG/ClinSV)

### database

- [GATK](https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references)
- [Reference gap](https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/gap.txt.gz)

### some files in `lib`

- hg19.genome.bed

~~~(shell)
# Hg38 from GATK
cut -f 1,2 Hg38.fasta.fai | head -n 24 | sort -k1,1  > Hg38.genome
# Hg19 from GATK
cut -f 1,2 Hg19.fasta.fai | sed -n '2,25p' | sort -k1,1  > Hg19.genome
zcat gap.txt.gz | awk -F"\t" '{if ($2~/^.{4,5}$/) print $2"\t"$3"\t"$4 }' | bedtools sort -i stdin > gap.bed
bedtools complement -i gap.bed -g Hg38.genome > Hg38.genome.bed
~~~