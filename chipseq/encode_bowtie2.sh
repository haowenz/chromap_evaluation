#!/bin/bash
#
#SBATCH --job-name=chip-seq_bowtie2
#SBATCH --output=slurm_chip-seq_bowtie2.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --open-mode=append

module load samtools/1.10
date
/usr/bin/time -v bowtie2 -X2000 --threads 8 -x GRCh38_no_alt_analysis_set_GCA_000001405.15 -1 VCaP/CTCF/r1/ENCFF950MHM.fastq.gz -2 VCaP/CTCF/r1/ENCFF858GIB.fastq.gz | samtools view -S -h -b -F 1804 -f 2 -q 30 /dev/stdin | samtools sort --threads 7 -n -o bowtie21_q30.bam

date
/usr/bin/time -v samtools fixmate --threads 7 -r bowtie21_q30.bam bowtie21_q30_fixmate.bam
/usr/bin/time -v samtools view -h -b -F 1804 -f 2 bowtie21_q30_fixmate.bam | samtools sort --threads 7 -o bowtie21_q30_1.bam
/usr/bin/time -v java -jar ~/Softwares/picard.jar MarkDuplicates INPUT=bowtie21_q30_1.bam OUTPUT=/dev/stdout METRICS_FILE=bowtie21_marked_dup_metrics.txt VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false QUIET=true COMPRESSION_LEVEL=0 |  samtools view -h -b -F 1804 -f 2 /dev/stdin > bowtie21_q30_processed.bam 
date

