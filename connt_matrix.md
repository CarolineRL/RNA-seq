# count_matrix
```
#!/bin/bash
#WangHongke 2020-10-14
## Quantify gene expression
path="/BioII/lulab_b/wanghongke/lulab_data/THCA_Huaxi/FTC"
gff3="/BioII/lulab_b/shared/genomes/hg38/gff3/miRBase.gff3"
count_reads="/BioII/lulab_b/wanghongke/exSeek/bin/count_reads.py"
dataset="FTC"
mkdir -p ${path}/output/${dataset}/count_matrix
mkdir -p ${path}/temp
input="${path}/output/${dataset}/mapping"
output="${path}/output/${dataset}/count_matrix"
temp="${path}/temp"

# Construct miRNA expression matrix
target="miRNA"
files=$(ls ${input})
for file in ${files}
do
sample=$(echo ${file} | cut -d "_" -f 1)
echo "Start count ${target} for ${file} in `date`"
#samtools sort --threads 8 -T ${temp}/${file} -o ${input}/${file}/bam-deduped/miRNA_sort_by_coordinate.bam ${input}/${file}/bam-deduped/miRNA.bam
samtools index -@ 8 ${input}/${file}/bam-deduped/${target}.bam
${count_reads} count_mature_mirna -a ${gff3} -i ${input}/${file}/bam-deduped/${target}.bam -o ${output}/${sample}.${target}.count
cut -f 1 ${output}/${sample}.${target}.count | sed -e "1i all_transcript_id" > ${output}/${target}.id.txt
cut -f 2 ${output}/${sample}.${target}.count | sed -e "1i \\${sample}"> ${output}/${sample}.${target}.txt
done
paste ${output}/${target}.id.txt  ${output}/*${target}.txt | sort -t $'\t' -k 1 -u > ${output}/${target}_count_matrix.txt
echo "Finish ${target}_count_matrix in `date`"
# Construct circRNA expression matrix
target="circRNA"
files=$(ls ${input})
for file in ${files}
do
sample=$(echo ${file} | cut -d "_" -f 1)
echo "Start count ${target} for ${file} in `date`"
samtools index -@ 8 ${input}/${file}/bam-deduped/${target}.bam
${count_reads} count_circrna -s forward -i ${input}/${file}/bam-deduped/${target}.bam -o ${output}/${sample}.${target}.count
cut -f 1 ${output}/${sample}.${target}.count > ${output}/${target}.id.txt
cut -f 2 ${output}/${sample}.${target}.count > ${output}/${sample}.${target}.txt
done
paste ${output}/${target}.id.txt  ${output}/*${target}.txt | sort -t $'\t' -k 1 -u > ${output}/${target}_count_matrix.txt
echo "Finish ${target}_count_matrix in `date`"
# Construct transcript expression matrix
for file in ${files}
do
sample=$(echo ${file} | cut -d "_" -f 1)
for target in piRNA lncRNA mRNA tRNA srpRNA snRNA snoRNA tucpRNA
do
samtools index -@ 8 ${input}/${file}/bam-deduped/${target}.bam
${count_reads} count_transcript -s forward -i ${input}/${file}/bam-deduped/${target}.bam -o ${output}/${sample}.${target}.count
#sed -e "1i \\${sample}"
cut -f 1 ${output}/${sample}.${target}.count > ${output}/${target}.transcript_id.txt
cut -f 2 ${output}/${sample}.${target}.count > ${output}/${sample}.transcript_${target}.txt
done
cat ${output}/*.transcript_id.txt >${output}/transcript.id.txt
cat ${output}/${sample}.transcript_*.txt >${output}/${sample}.transcript.count.txt
done
paste ${output}/transcript.id.txt ${output}/*transcript.count.txt >${output}/transcript.matrix.txt
cat ${output}/miRNA_count_matrix.txt ${output}/circRNA_count_matrix.txt ${output}/transcript.matrix.txt>${output}/transcript.count_matrix.txt
echo "Finish transcript_count_matrix in `date`"
for file in `ls ${output} | grep -v 'count_matrix.txt'`
do
rm ${output}/${file}
done
rm -r ${temp}
```
