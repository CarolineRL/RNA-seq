# RNA-seq
```
#dataset = "GSE131483"
dataset = "FTC_new"
sample_ids = open ("sample_id/{}.txt".format(dataset)).read().strip().split("\n")
sequences = ["spikein_small","univec","rRNA","miRNA","lncRNA","mRNA","piRNA","snoRNA","snRNA","srpRNA","tRNA","tucpRNA","Y_RNA","genome","circRNA"]




rule all:
    input:
        qc1 = expand("output/{dataset}/qc1/{sample_id}_fastqc.html",dataset=[dataset],sample_id = sample_ids),
        bams = expand("output/{dataset}/mapping/{sample_id}/bam-deduped/{sequence}.bam", sample_id=sample_ids, dataset= dataset, sequence =  sequences)

rule extractUMI:
    input:
        fastq = "data/{dataset}/{sample_id}.fastq.gz"
    output:
        fastq = "output/{dataset}/extract-UMI/{sample_id}.fastq.gz" 
    log:
        "output/{dataset}/log/extract-UMI/{sample_id}.txt" 
    shell:
        """
        umi_tools extract --extract-method=regex  -p "^.+(?P<discard_1>AACTGTAGGCACCATCAAT)(?P<umi_1>.{{12}})(?P<discard_2>.+)" -I {input.fastq} --stdout={output.fastq} -L {log}
        """


rule trim:
    input:
        fastq = "output/{dataset}/extract-UMI/{sample_id}.fastq.gz"
    output:
        fastq = "output/{dataset}/trimmed/{sample_id}.fastq.gz"
    log:
        "output/{dataset}/log/trim/{sample_id}.txt"
    shell:
        """
        trim_galore --quality 30 --length 15 --trim-n --phred33 -o output/{wildcards.dataset}/trimmed --basename {wildcards.sample_id} {input.fastq} > {log} 2>&1
        mv output/{wildcards.dataset}/trimmed/{wildcards.sample_id}_trimmed.fq.gz {output.fastq}
        """


rule qc1:
    input:
        fastq = "output/{dataset}/trimmed/{sample_id}.fastq.gz"
    output:
        report = "output/{dataset}/qc1/{sample_id}_fastqc.html"
    params:
        outdir = "output/{dataset}/qc1"
    shell:
        """ls
        fastqc -o {params.outdir} {input.fastq}
        """

rule mapping:
    input:
        fastq = "output/{dataset}/trimmed/{sample_id}.fastq.gz"
    output:
        bams = [ "output/{dataset}/mapping/{sample_id}/bam/" + sequence + ".bam" for sequence in sequences ]
    threads:
        4
    shell:
        """
        scripts/mapping.py --fastq {input.fastq} --bam-dir output/{wildcards.dataset}/mapping/{wildcards.sample_id}/bam \
        --log-dir output/{wildcards.dataset}/log/mapping/{wildcards.sample_id} --threads {threads} \
        --unmapped-dir output/{wildcards.dataset}/mapping/{wildcards.sample_id}/unmapped
        """ 
        
        
rule map_spikein:
    input:
        fastq = "output/{dataset}/trimmed/{sample_id}.fastq.gz"
        index= '/index/spikein_small.1.bt2'
    output:
        unmapped='output/{dataset}/mapping/{sample_id}/unmapped/spikein.fa.gz',
        bams=[ 'output/{dataset}/mapping/{sample_id}/bam/spikein.bam' + sequence + ".bam" for sequence in sequences ]
    log:
        "output/{dataset}/log/mapping/{sample_id}/spikein.txt"
    params:
        index='/index/spikein_small'
    threads: 
        4
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --norc --sensitive --no-unal \
            --un-gz {output.unmapped} -x {params.index} - -S - > {log}\
        | samtools view -b -o {output.bam}
        '''



rule dedup:
    input:
        bam = "output/{dataset}/mapping/{sample_id}/bam/{sequence}.bam"
    output:
        bam = "output/{dataset}/mapping/{sample_id}/bam-deduped/{sequence}.bam"
    log:
        "output/{dataset}/log/dedup/{sample_id}/{sequence}.log"
    shell:
        """
        umi_tools dedup -I {input.bam} -S {output.bam} > {log}
        """
```
