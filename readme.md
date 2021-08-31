# VarTools

## funcion
1. f2v: analysis from fastq to gvcf. 
2. gt: from gvcf created by GATK to vcf. 
3. trio: trio analysis

### f2v
~~~(shell)
python3 VarTools f2v -i indir  -o ./result
~~~

### trio_genotype
~~~(shell)
python3 VarTools tGT -p .gvcf -f .gvcf -m .gvcf -s .gvcf,.gvcf -o ./result/outName
~~~

### single_genotype
~~~(shell)
python3 VarTools sGT -p .gvcf -o ./result/outName
~~~

### 软件列表
- [GATK](https://github.com/broadinstitute/gatk/releases/download/4.2.0.0/gatk-4.2.0.0.zip)
- [fastp](http://opengene.org/fastp/fastp)
- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip)
- [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2)
- [qualimap](https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v2.2.1.zip)
- [plink2](https://s3.amazonaws.com/plink2-assets/alpha2/plink2_linux_x86_64.zip)
- [kings](https://www.kingrelatedness.com/Linux-king.tar.gz)
- [VariantQC](https://github.com/BimberLab/DISCVRSeq/releases/download/1.3.2/DISCVRSeq-1.3.2.jar)
- [annovar](http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz)
- [clinSV](https://github.com/KCCG/ClinSV)

### 数据库列表
- [GATK](https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references)