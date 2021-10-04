# diff_exp
```
#!/bin/bash
#WangHongke 2021-08-15 Differential expression analysis
path="/BioII/lulab_b/wanghongke/lulab_data/THCA_Huaxi/FTC"
dataset="FTC"
mkdir -p ${path}/output/${dataset}/diff_exp/16vs14
input="${path}/output/${dataset}/count_matrix"
output="${path}/output/${dataset}/diff_exp/16vs14"
positive_group="FTC"
negative_group="FTA"
# miRNA diff_exp
python ${path}/scripts/filter.py -i ${input}/miRNA_count_matrix.txt --method by_value --threshold 0 --proportion 0.2 -o ${output}/filtered_miRNA_count_matrix.txt
/usr/bin/Rscript ${path}/scripts/RPM.R -i ${output}/filtered_miRNA_count_matrix.txt -o ${output}/filtered_miRNA_RPM.txt
/usr/bin/Rscript ${path}/scripts/differential_expression.R -i ${output}/filtered_miRNA_count_matrix.txt -p ${path}/sample_class/${positive_group}.txt -n ${path}/sample_class/${negative_group}.txt -m edger_glmlrt --norm-method TMM -o ${output}/${positive_group}vs${negative_group}.edgeR.glmlrt.txt
awk -F "\t" 'NR==1{print}' ${output}/${positive_group}vs${negative_group}.edgeR.glmlrt.txt >${output}/${positive_group}vs${negative_group}_up_down.txt
awk -F "\t" '{if($6 < 0.05 && ($2 >= 1 || $2 <= -1)) print $0}' ${output}/${positive_group}vs${negative_group}.edgeR.glmlrt.txt >>${output}/${positive_group}vs${negative_group}_up_down.txt
```
