read_length=${1}

ref=GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
ref_prefix=GRCh38_no_alt_analysis_set_GCA_000001405.15
num_frags=1000000
read1_output=${read_length}_simulated1.fq
read2_output=${read_length}_simulated2.fq
read_output=${read_length}_simulated.fq
sam_output=${read_length}_simulated.sam

paftools="k8 minimap2/misc/paftools.js"
SEQTK="seqtk"

mason_simulator -ir ${ref} -n ${num_frags} -o ${read1_output} -or ${read2_output} -oa ${sam_output} --illumina-read-length ${read_length} --illumina-prob-mismatch-scale 0.25

${paftools} mason2fq ${sam_output} > ${read_output}
${SEQTK} seq -1 ${read_output} > ${read1_output}
${SEQTK} seq -2 ${read_output} > ${read2_output}

bwa=bwa
bwa_output=${read_length}_simulated_bwa.sam
bwa_eval_output=${read_length}_simulated_bwa.eval
${bwa} index ${ref} 
${bwa} mem -t 32 ${ref} ${read1_output} ${read2_output} > ${bwa_output} 
${paftools} mapeval ${bwa_output} > ${bwa_eval_output} 

bowtie2=bowtie2
bowtie2_build=bowtie2-build
bowtie2_output=${read_length}_simulated_bowtie2.sam
bowtie2_eval_output=${read_length}_simulated_bowtie2.eval
${bowtie2_build} --threads 32 ${ref} ${ref_prefix}
${bowtie2} -X2000 --threads 32 -x ${ref_prefix} -1 ${read1_output} -2 ${read2_output} > ${bowtie2_output}
${paftools} mapeval ${bowtie2_output} > ${bowtie2_eval_output} 

minimap2=minimap2
minimap2_output=${read_length}_simulated_minimap2.sam
minimap2_eval_output=${read_length}_simulated_minimap2.eval
${minimap2} -ax sr -t 32 ${ref} ${read1_output} ${read2_output} > ${minimap2_output}
${paftools} mapeval ${minimap2_output} > ${minimap2_eval_output} 

chromap=chromap
chromap_index=GRCh38_no_alt_analysis_set_GCA_000001405.15_k17_w7.chromap # For 50bp reads
#chromap_index=GRCh38_no_alt_analysis_set_GCA_000001405.15_k23_w11.chromap # For 150bp reads and longer
chromap_output=${read_length}_simulated_chromap.paf
chromap_eval_output=${read_length}_simulated_chromap.eval
/usr/bin/time -v ${chromap} -t 1 --PAF -q 0 -r ${ref} -x ${chromap_index} -1 ${read1_output} -2 ${read2_output} -o ${chromap_output}
/usr/bin/time -v ${paftools} mapeval ${chromap_output} > ${chromap_eval_output} 
