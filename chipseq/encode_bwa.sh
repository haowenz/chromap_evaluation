#!/bin/bash
#
#SBATCH --job-name=chip-seq_bwa
#SBATCH --output=slurm_chip-seq_bwa.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --nodelist node05
#SBATCH --open-mode=append

module load samtools/1.10

date
/usr/bin/time -v bwa mem -t 8 GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta VCaP/CTCF/r1/ENCFF950MHM.fastq.gz VCaP/CTCF/r1/ENCFF858GIB.fastq.gz | samtools view -S -h -b -F 1804 -f 2 -q 30 /dev/stdin | samtools sort -n --threads 7 -o bwa1_q30.bam
date

/usr/bin/time -v samtools fixmate --threads 7 -r bwa1_q30.bam bwa1_q30_fixmate.bam
/usr/bin/time -v samtools view -h -b -F 1804 -f 2 bwa1_q30_fixmate.bam | samtools sort --threads 7 -o bwa1_q30_1.bam
/usr/bin/time -v java -jar ~/Softwares/picard.jar MarkDuplicates INPUT=bwa1_q30_1.bam OUTPUT=/dev/stdout METRICS_FILE=bwa1_marked_dup_metrics.txt VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false QUIET=true COMPRESSION_LEVEL=0 | samtools view -h -b -F 1804 -f 2 /dev/stdin > bwa1_q30_processed.bam 

date

