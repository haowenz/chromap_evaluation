#!/bin/sh

module load bedtools

awk '{print $1"\t"$2"\t"$3"\t"1}' bwa_filtered_peaks.narrowPeak > bwa_peak_preunion.bed
awk '{print $1"\t"$2"\t"$3"\t"2}' bowtie2_filtered_peaks.narrowPeak > bowtie2_peak_preunion.bed
awk '{print $1"\t"$2"\t"$3"\t"4}' chromap_filtered_peaks.narrowPeak > chromap_peak_preunion.bed

cat *_preunion.bed > preunion.bed
bedtools sort -i preunion.bed > preunion.sorted.bed
bedtools merge -o distinct -c 4 -i preunion.sorted.bed > union.bed
cat union.bed| cut -f4 | sort | uniq -c
