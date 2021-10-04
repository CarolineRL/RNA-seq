# RNA-seq Trim_galore
conda install -c grst trim_galore

## 优劣势分析
### 优势
1. 安装和使用都非常简单；
2. 代码较短
3. 参数更直观，不用去死记硬背
4. 默认下paired和unpaired，运行速度更快
### 劣势
1. 可调参数较少

   除此之外，trim_galore还有一个特色就是可以对RRBS（Reduced Representation Bisulfite Sequencing， DNA甲基化测序）构建的文本库进行相应的质控，比如“–rrbs” 参数对于 RRBS测序中引入的碱基C进行处理，会去除3’端的2个碱基。具体参照[说明书](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md)。

## Trim_galore使用说明
### Trim_galore问价相关参数
Trim_galore处理序列的流程可以分为以下几个步骤：
   1. quality trimming质量分数相关处理
   2. adapter trimming街头处理
   3. removing Short Sequences过短reads处理
   4. specialised Trimming 特殊处理

输入参数是最好根据上述流程输入各自对应的参数，常用命令：
```
trim_galore -q 25 --phred33 --fastqc --length 36 --stringency 3 --paired -o ../cleanData *${sample}*.gz 
```

参数解释

-q/–quality    
对RRBS样本来说，quanlity trimming结束后才会进行adapter trimming；在一般样本中，两者同时进行。此处默认质量分数为20。

–phred33   -phred64 
默认编码为phred33编码，具体用33还是64取决于测序平台。

-fastqc
在处理结束后自动运行FastQC

-illumina/nextera/small_rna
跟trimmomatic一样，根据adapter comtent报告要去除的街头序列类型，如果不选择该类参数，trim_galore会根据前一百万个剪辑自动判断街头序列；演奏注意的是，当选择-small_rna时，-length参数会自动设置为18bp。

-stringency严格度
街头序列最小配对碱基数：简单来说就是最多能允许末端残留多少个街头序列的碱基，默认值为极端值1；该参数与trimmomatic中ILLUMINACLIP<minAdapterLength>含义相同。 
-length
保留reads的最短长度。设置为0时代表不选择该选项。默认为20bp。对双端测序来说，一旦有一段低于length值，另一端必须选择是否保留。

  ### 其他参数
  –clip_R1
切掉Read 1 5端的N个碱基

–clip_R2
切掉Read 2 5端的N个碱基

–three_prime_clip_R1
在接头序列/质量差碱基切除后，切掉read 1 3端的N个碱基

–three_prime_clip_R2
在接头序列/质量差碱基切除后，切掉read 2 3端的N个碱基

以上选项主要是reads两端出现一些Bias的碱基（非接头序列或质量分数相关）时可以用到。

–basename <PREFERRED_NAME>

–basename类似trimmomatic中的baseout，对双端测序会自动输出 PREFERRED_NAME_val_1.fq(.gz) and PREFERRED_NAME_val_2.fq(.gz) 。

-j/–cores INT

运行程序的核心数。Python2环境下只能为1。如果环境中已安装python3和pigz（parallel gzip）gzip，可以选择多核处理（但要注意gzip会严重拖累多核处理的速度）。一般选择–cores = 4就可以了，这里的core数虽然是4，但其实代表的是每个过程使用的核心数。
–cores 4 代表: 4 个用于读取 + 4 个用于写入 + 4 个用于Cutadapt+ 2 个用于Cutadapt的其他任务+ 1 个用于trim_galore本身 = 总共15个。

–paired
表明是双端测序文件

–retain_unpaired
在双端测序某一端处理后过短时，默认为舍弃该段的两端reads, 添加此参数则保留reads，长一点的那一端文件会自动保存为.unpaired_1.fq 或者 .unpaired_2.fq
   
   
   
   
   ## trim.sh(建议的方法)
   ```
   vi trim.sh
   ```
   
   写入以下内容
   ```
   dir=/BioII/lulab_b/qiaoruolin/FTC2/output/FTC/trimmed/
cat config |while read id
do
      arr=${id}
      fq1=${arr[0]}
      fq2=${arr[1]}
      nohup trim_galore -q 30 --phred33 --length 15 --trim-n -o $dir $fq1 $fq2 &
done
   ```
   
