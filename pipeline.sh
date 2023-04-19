#!/bin/bash

# Make directory tree structure realative to pipeline.sh directory

mkdir ngs_pipeline_assessment
cd ngs_pipeline_assessment
mkdir dnaseq
cd dnaseq
mkdir data logs meta results
cd data
mkdir untrimmed_fastq
mkdir trimmed_fastq

wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

cd untrimmed_fastq

wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz

fastqc *.fastq.qz

mkdir ../../results/fastqc_untrimmed_reads
mv *fastqc* ../../results/fastqc_untrimmed_reads/

cd ../../results/fastqc_untrimmed_reads

for zip in *.zip; do unzip $zip; done

cat */summary.txt > ~/ngs_pipeline_assessment/dnaseq/logs/fastqc_summaries.txt

mv NGS0001.R1.fastq.qz NGS0001.R1.fastq.gz
mv NGS0001.R2.fastq.qz NGS0001.R2.fastq.gz

cd ../../../../

trimmomatic PE -threads 4 \
 -phred33 \
 ngs_pipeline_assessment/dnaseq/data/untrimmed_fastq/NGS0001.R1.fastq.gz \
 ngs_pipeline_assessment/dnaseq/data/untrimmed_fastq/NGS0001.R2.fastq.gz \
 -baseout ngs_pipeline_assessment/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R \
 ILLUMINACLIP:anaconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 \
 TRAILING:25 MINLEN:50

cd ngs_pipeline_assessment/dnaseq/data/trimmed_fastq

fastqc NGS0001_trimmed_R_1P
fastqc NGS0001_trimmed_R_2P

mkdir ../../results/fastqc_trimmed_reads
mv *fastqc* ../../results/fastqc_trimmed_reads/

for zip in *.zip; do unzip $zip; done

cat */summary.txt > ~/ngs_pipeline_assessment/dnaseq/logs/trimmed_fastqc_summaries.txt

mkdir -p ~/ngs_pipeline_assessment/dnaseq/data/reference
mv ~/ngs_pipeline_assessment/dnaseq/data/hg19.fa.gz ~/ngs_pipeline_assessment/dnaseq/data/reference/

bwa index ~/ngs_pipeline_assessment/dnaseq/data/reference/hg19.fa.gz
bwa mem -t 4 -v 1 -R '@RG\tID:hiSeq2500\tPL:ILLUMINA\tSM:S_0\tLB:Nextera\tPU:111' \
 ~/ngs_pipeline_assessment/dnaseq/data/reference/hg19.fa.gz \
 ~/ngs_pipeline_assessment/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R_1P \
 ~/ngs_pipeline_assessment/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R_2P > ~/ngs_pipeline_assessment/dnaseq/data/aligned_data/NGS0001_aligned.sam

java -jar ~/anaconda3/share/picard-2.18.29-0/picard.jar SamFormatConverter INPUT=NGS0001_aligned.sam OUTPUT=NGS0001_aligned.bam

java -jar ~/anaconda3/share/picard-2.18.29-0/picard.jar SortSam INPUT=NGS0001_aligned.bam OUTPUT=NGS0001_aligned_sorted.bam SORT_ORDER=coordinate

java -jar ~/anaconda3/share/picard-2.18.29-0/picard.jar BuildBamIndex INPUT=NGS0001_aligned_sorted.bam

java -jar ~/anaconda3/share/picard-2.18.29-0/picard.jar MarkDuplicates I=NGS0001_aligned_sorted.bam O=NGS0001_aligned_dupmarked.bam  M=marked_dup_metrics.txt

java -jar ~/anaconda3/share/picard-2.18.29-0/picard.jar BuildBamIndex INPUT=NGS0001_aligned_dupmarked.bam

samtools view -F 1796 -q 20 -o NGS0001_aligned_filtered_dupmarked.bam NGS0001_aligned_dupmarked.bam

samtools index NGS0001_aligned_filtered_dupmarked.bam

samtools flagstat NGS0001_aligned_dupmarked.bam

samtools idxstats ~/ngs_pipeline_assessment/dnaseq/data/aligned_data/NGS0001_aligned_dupmarked.bam > alignment_stats.txt

bedtools coverage -a annotation.bed -b NGS0001_aligned_dupmarked.bam > depth_coverage.txt

java -jar ~/anaconda3/share/picard-2.18.29-0/picard.jar CollectInsertSizeMetrics \
 I=NGS0001_aligned_dupmarked.bam \
 O=insert_size_metrics.txt \
 H=insert_size_histogram.pdf

zcat ~/ngs_piepline_assessment/dnaseq/data/reference/hg19.fa.gz > ~/ngs_course/dnaseq/data/reference/hg19.fa 

samtools faidx ~/ngs_piepline_assessment /dnaseq/data/reference/hg19.fa

freebayes --bam ~/ngs_piepline_assessment /dnaseq/data/aligned_data/ NGS0001_aligned_dupmarked.bam --fasta-reference ~/ngs_piepline_assessment /dnaseq/data/reference/hg19.fa --vcf ~/ngs_piepline_assessment /dnaseq/results/ NGS0001_aligned_dupmarked.vcf

bgzip ~/ngs_piepline_assessment /dnaseq/results/ NGS0001_aligned_dupmarked.vcf

tabix -p vcf ~/ngs_piepline_assessment /dnaseq/results/ NGS0001_aligned_dupmarked.vcf.gz

vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \
 ~/ngs_pipeline_assessment/dnaseq/results/ NGS0001_aligned_dupmarked.vcf.gz > ~/ngs_pipeline_assessment/dnaseq/results/ NGS0001_aligned_dupmarked_filtered.vcf

bedtools intersect -header -wa -a ~/ngs_pipeline_assessment/dnaseq/results/ NGS0001_aligned_dupmarked_filtered.vcf -b ../annotation.bed \
 ~/ngs_pipeline_assessment/dnaseq/results/ NGS0001_aligned_dupmarked_filtered _annotation.vcf

bgzip ~/ngs_pipeline_assessment/dnaseq/results/NGS0001_aligned_dupmarked_filtered _annotation.vcf

tabix -p vcf ~/ngs_pipeline_assessment/dnaseq/results/NGS0001_aligned_dupmarked_filtered _annotation.vcf.gz

tar -zxvf annovar.latest.tar.gz
cd annovar

./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/

./convert2annovar.pl \
 -format vcf4 ~/ngs_pipeline_assessment/dnaseq/results/NGS0001_aligned_dupmarked_filtered _annotation.vcf.gz > \
 ~/ngs_pipeline_assessment/dnaseq/results/NGS0001_aligned_dupmarked_filtered_annotation.avinput

./table_annovar.pl ~/ngs_pipeline_assessment/dnaseq/results/NGS0001_aligned_dupmarked_filtered_annotation.avinput humandb/ -buildver hg19 \
 -out ~/ngs_pipeline_assessment/dnaseq/results/ NGS0001_aligned_dupmarked_filtered _annotation -remove \
 -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout
