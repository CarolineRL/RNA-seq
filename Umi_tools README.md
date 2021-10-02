# RNA-seq
  用来处理带有分子标签(Unique Molecular Identifiers，UMIs)和细胞条形码(cell barcode)的数据。目前来说，比较主流的单细胞测序数据都包含UMIs和barcode，下面就来演示这款工具是如何处理单细胞数据的。
## 功能概况
  ![image](https://user-images.githubusercontent.com/91002411/135707809-df1928e1-e530-48dd-b9df-4adbe795ca34.png)
### 软件安装方式
```
#第一种
conda install -c bioconda -c conda-forge umi_tools
#第二种
pip install umi_tools
#第三种
unzip 1.0.0.zip
cd UMI-tools-1.0.0
python setup.py install --user
```

### 具体使用步骤
1. 提取cell barcode白名单
  whitelist 命令会从原始数据种提取去可能的cell barcode。通常情况下，10X的barcode长度为16nt，umi长度为10nt；Drop-seq的barcode长度为12nt，umi长度为8nt。示例代码如下：
```
umi_tools whitelist --stdin hgmm_100_R1.fastq.gz \  
                    --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \ 
                    --set-cell-number=100 \    
                    --log2stderr > whitelist.txt
```

  常用参数解释
  --stdin：指定输入文件或者标准输入。
  --plot-prefix：指定QC统计图的前缀，可用于判断细胞数阈值设置是否合理。
  --bc-pattern：指定cell barcode和umi的位置，默认情况下，他们位于序列的5'端，否则可用--3prime参数改变。cell barcode用同等数量的字符"C"表示；umi用同等数量的字符"N"表示。
  --set-cell-number：设置检测到的cell barcode数量，若事先知道数量可设定，否则省略软件会自动判断(结合QC图判读阈值是否合理，若不合理可人为设定阈值重新分析)。
  --expect-cells=200：设置检测到的cell barcode数据上限，该值是根据捕获效率预估得到(一般10X捕获效率不低于10%)，结合QC图判断阈值是否合理，若不合理可人为设定阈值重新分析。
  --stdout/-S：指定输出barcode结果到文件，log信息到还是标准输出。
  --stdlog/-L：指定log信息到文件，barcode结果还是到标准输出。
  --log2stderr：指定log信息到标准错误输出，barcode结果还是到标准输出。
  -v 0：关闭log信息的输出。

  结果文件解释：
  whitelist生成的结果文件包含四列：1、可接受的cell barcode；2、与可接受的barcode汉明距离相差1的barcode；3、第一列barcode的umi数；4、第3列barcode的umi数

  结果文件格式如下：
```
AAAGATGAGAAACGAG AAAAATGAGAAACGAG,AAACATGAGAAACGAG,... 53122 4,6,...
AAAGCAAGTACCTACA AAAACAAGTACCTACA,AAACCAAGTACCTACA,... 36255 2,3,...
AACACGTCAGCGTAAG AAAACGTCAGCGTAAG,AACAAGTCAGCGTAAG,... 53133 4,11,...
```

2. 提取barcode和过滤reads
extract命令会从fastq文件中提取包含可接受barcode的reads，默认情况下extract命令会忽略umi的reads质量情况而不做处理。示例代码如下：
```
umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                  --stdin hgmm_100_R1.fastq.gz \
                  --stdout hgmm_100_R1_extracted.fastq.gz \
                  --read2-in hgmm_100_R2.fastq.gz \
                  --read2-out=hgmm_100_R2_extracted.fastq.gz \
                  --filter-cell-barcode \
                  --whitelist=whitelist.txt
```

常用参数解释：
--bc-pattern：指定cell barcode和umi的位置，同whilelist。
--stdin：指定输入文件或者标准输入，同whilelist。
--stdout：指定read1的输出文件。
--read2-in：指定read2的输入文件，即基因表达的原始reads文件。
--read2-out：指定read2的输出文件。
--read2-stdout：设置read2的输出到标准输出，同时不会生成提取后read1文件。
--filter-cell-barcode：指定只提取包含可接受barcoded的reads。
--error-correct-cell：指定提取包含与可接受barcode汉明距离相差1的barcode的reads。
--whitelist：barcode白名单文件。
--quality-filter-threshold=[FILTER]：设定通过阈值来过滤umi的read，低于阈值的umi会被丢弃。
--quality-filter-mask=[FILTER]：设置一个阈值来把umi reads中低于阈值的碱基替换为"N"。

3. 比对
使用STAR软件来将reads比对到参考基因组。示例代码如下：
```
STAR --runThreadN 4 \
       --genomeDir hg38_noalt_junc85_99.dir \
       --readFilesIn hgmm_100_R2_extracted.fastq.gz \
       --readFilesCommand zcat \
       --outFilterMultimapNmax 1 \
       --outSAMtype BAM SortedByCoordinate
```

4. 比对到基因
使用软件featureCounts确定每一条read比对到基因的位置，会生成一个新的包含每条read回帖到基因位置的bam文件，该软件来自于Subread软件包，下载时应选择版本大于1.5.3的，subread链接：https://sourceforge.net/projects/subread/files/subread-2.0.1/。
示例代码如下：
```
featureCounts -a geneset.gtf -o gene_assigned -R BAM Aligned.sortedByCoord.out.bam -T 4
```
featureCounts生成的bam文件没有排序，需要用samtools排序并建索引方便后续使用。示例代码如下：
```
samtools sort Aligned.sortedByCoord.out.bam.featureCounts.bam -o assigned_sorted.bam
samtools index assigned_sorted.bam
```

常用参数解释：
--per-gene：指定给每个基因统计umi数。
--gene-tag=XT：指定基因标记。
--assigned-status-tag=XS：指定基因标记 。
--per-cell：指定给每个细胞统计umi数。
--skip-tags-regex：指定跳过的标记，默认值为^[__|Unassigned]。
--wide-format-cell-counts：设置输出结果的格式为宽矩阵，即行为基因，列为细胞。

默认结果格式如下：
```
$ zcat counts.tsv.gz | head
gene cell count
ENSG00000000003 AAAGATGAGAAACGAG 3
ENSG00000000003 AACTCTTGTTCTGAAC 4
ENSG00000000003 ACACCGGGTACGACCC 2
ENSG00000000003 ACACTGAGTCGGGTCT 5
ENSG00000000003 ACTATCTCAAGGTGTG 2
```

## 最后
官网有英文教程，官网链接：https://umi-tools.readthedocs.io/en/latest/index.html

