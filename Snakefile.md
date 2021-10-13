# Snakefile
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
fastqc -o output/qc1 {dataset}/{sample_id}.fastq.gz

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
        index=genome_dir + '/index/spikein_small.1.bt2'
    output:
        unmapped='output/{dataset}/mapping/{sample_id}/unmapped/spikein.fa.gz',
        bams=[ 'output/{dataset}/mapping/{sample_id}/bam/spikein.bam' + sequence + ".bam" for sequence in sequences ]
    log:
        "output/{dataset}/log/mapping/{sample_id}/spikein.txt"
    params:
        index=genome_dir + '/index/spikein_small'
    threads: 
        4
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --norc --sensitive --no-unal \
            --un-gz {output.unmapped} -x {params.index} - -S - > {log}\
        | samtools view -b -o {output.bam}
        '''

rule map_univec:
    input:
        reads='output/{dataset}/mapping/{sample_id}/unmapped/spikein.fa.gz',
        index=genome_dir + '/index/univec.1.bt2'
    output:
        unmapped='output/{dataset}/mapping/{sample_id}/unmapped/univec.fa.gz',
        bams=[ 'output/{dataset}/mapping/{sample_id}/bam/univec.bam' + sequence + ".bam" for sequence in sequences ]
    log:
        "output/{dataset}/log/mapping/{sample_id}/univec.txt"
    params:
        index=genome_dir + '/index/univec'
    threads: 
        4
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --norc --sensitive --no-unal \
            --un-gz {output.unmapped} -x {params.index} - -S - > {log}\
        | samtools view -b -o {output.bam}
        '''

rule map_rRNA:
    input:
        reads='output/{dataset}/mapping/{sample_id}/unmapped/univec.fa.gz',
        index=genome_dir + '/index/rRNA.1.bt2'
    output:
        unmapped='output/{dataset}/mapping/{sample_id}/unmapped/rRNA.fa.gz',
        bam=[ 'output/{dataset}/mapping/{sample_id}/bam/rRNA.bam' + sequence + ".bam" for sequence in sequences ]
    log:
        "output/{dataset}/log/mapping/{sample_id}/rRNA.txt"
    params:
        index=genome_dir + '/index/rRNA'
    threads: 
        4
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --norc --sensitive --no-unal \
            --un-gz {output.unmapped} -x {params.index} - -S - > {log}\
        | samtools view -b -o {output.bam}
        '''

rule map_lncRNA:
    input:
        reads='output/{dataset}/mapping/{sample_id}/unmapped/rRNA.fa.gz',
        index=genome_dir + '/index/bowtie2/lncRNA.1.bt2'
    output:
        unmapped='output/{dataset}/mapping/{sample_id}/unmapped/lncRNA.fa.gz',
        bam=[ 'output/{dataset}/mapping/{sample_id}/bam/lncRNA.bam' + sequence + ".bam" for sequence in sequences ]
    log:
        "output/{dataset}/log/mapping/{sample_id}/lncRNA.txt"
    params:
        index=genome_dir + '/index/lncRNA'
    threads: 
        4
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --norc --sensitive --no-unal \
            --un-gz {output.unmapped} -x {params.index} - -S - > {log}\
        | samtools view -b -o {output.bam}
        '''

rule map_miRNA:
    input:
        reads='output/{dataset}/mapping/{sample_id}/unmapped/lncRNA.fa.gz',
        index=genome_dir + '/index/miRNA.1.bt2'
    output:
        unmapped='output/{dataset}/mapping/{sample_id}/unmapped/miRNA.fa.gz',
        bam=[ 'output/{dataset}/mapping/{sample_id}/bam/miRNA.bam' + sequence + ".bam" for sequence in sequences ]
    log:
        "output/{dataset}/log/mapping/{sample_id}/miRNA.txt"
    params:
        index=genome_dir + '/indexmiRNA'
    threads: 
        4
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --norc --sensitive --no-unal \
            --un-gz {output.unmapped} -x {params.index} - -S - > {log}\
        | samtools view -b -o {output.bam}
        '''

rule map_mRNA:
    input:
        reads='output/{dataset}/mapping/{sample_id}/unmapped/miRNA.fa.gz',
        index=genome_dir + '/index//mRNA.1.bt2'
    output:
        unmapped='output/{dataset}/mapping/{sample_id}/unmapped/mRNA.fa.gz',
        bam=[ 'output/{dataset}/mapping/{sample_id}/bam/mRNA.bam' + sequence + ".bam" for sequence in sequences ]
    log:
        "output/{dataset}/log/mapping/{sample_id}/mRNA.txt"
    params:
        index=genome_dir + '/index/mRNA'
    threads: 
        4
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --norc --sensitive --no-unal \
            --un-gz {output.unmapped} -x {params.index} - -S - > {log}\
        | samtools view -b -o {output.bam}
        '''

rule map_piRNA:
    input:
        reads='output/{dataset}/mapping/{sample_id}/unmapped/mRNA.fa.gz',
        index=genome_dir + '/index/piRNA.1.bt2'
    output:
        unmapped='output/{dataset}/mapping/{sample_id}/unmapped/piRNA.fa.gz',
        bam=[ 'output/{dataset}/mapping/{sample_id}/bam/piRNA.bam' + sequence + ".bam" for sequence in sequences ]
    log:
        "output/{dataset}/log/mapping/{sample_id}/piRNA.txt"
    params:
        index=genome_dir + '/index/piRNA'
    threads: 
        4
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --norc --sensitive --no-unal \
            --un-gz {output.unmapped} -x {params.index} - -S - > {log}\
        | samtools view -b -o {output.bam}
        '''

rule map_snoRNA:
    input:
        reads='output/{dataset}/mapping/{sample_id}/unmapped/piRNA.fa.gz',
        index=genome_dir + '/index/snoRNA.1.bt2'
    output:
        unmapped='output/{dataset}/mapping/{sample_id}/unmapped/snoRNA.fa.gz',
        bam=[ 'output/{dataset}/mapping/{sample_id}/bam/snoRNA.bam' + sequence + ".bam" for sequence in sequences ]
    log:
        "output/{dataset}/log/mapping/{sample_id}/snoRNA.txt"
    params:
        index=genome_dir + '/index/snoRNA'
    threads: 
        4
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --norc --sensitive --no-unal \
            --un-gz {output.unmapped} -x {params.index} - -S - > {log}\
        | samtools view -b -o {output.bam}
        '''

rule map_snRNA:
    input:
        reads='output/{dataset}/mapping/{sample_id}/unmapped/snoRNA.fa.gz',
        index=genome_dir + '/index/snRNA.1.bt2'
    output:
        unmapped='output/{dataset}/mapping/{sample_id}/unmapped/snRNA.fa.gz',
        bam=[ 'output/{dataset}/mapping/{sample_id}/bam/snRNA.bam' + sequence + ".bam" for sequence in sequences ]
    log:
        "output/{dataset}/log/mapping/{sample_id}/snRNA.txt"
    params:
        index=genome_dir + '/index/snRNA'
    threads: 
        4
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --norc --sensitive --no-unal \
            --un-gz {output.unmapped} -x {params.index} - -S - > {log}\
        | samtools view -b -o {output.bam}
        '''

rule map_srpRNA:
    input:
        reads='output/{dataset}/mapping/{sample_id}/unmapped/snRNA.fa.gz',
        index=genome_dir + '/index/srpRNA.1.bt2'
    output:
        unmapped='output/{dataset}/mapping/{sample_id}/unmapped/srpRNA.fa.gz',
        bam=['output/{dataset}/mapping/{sample_id}/bam/srpRNA.bam' + sequence + ".bam" for sequence in sequences ]
    log:
        "output/{dataset}/log/mapping/{sample_id}/srpRNA.txt"
    params:
        index=genome_dir + '/index/srpRNA'
    threads: 
        4
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --norc --sensitive --no-unal \
            --un-gz {output.unmapped} -x {params.index} - -S - > {log}\
        | samtools view -b -o {output.bam}
        '''

rule map_tRNA:
    input:
        reads='output/{dataset}/mapping/{sample_id}/unmapped/srpRNA.fa.gz',
        index=genome_dir + '/index/tRNA.1.bt2'
    output:
        unmapped='output/{dataset}/mapping/{sample_id}/unmapped/tRNA.fa.gz',
        bam=[ 'output/{dataset}/mapping/{sample_id}/bam/tRNA.bam' + sequence + ".bam" for sequence in sequences ]
    log:
        "output/{dataset}/log/mapping/{sample_id}/tRNA.txt"
    params:
        index=genome_dir + '/index/tRNA'
    threads: 
        4
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --norc --sensitive --no-unal \
            --un-gz {output.unmapped} -x {params.index} - -S - > {log}\
        | samtools view -b -o {output.bam}
        '''

rule map_tucpRNA:
    input:
        reads='output/{dataset}/mapping/{sample_id}/unmapped/tRNA.fa.gz',
        index=genome_dir + '/index/tucpRNA.1.bt2'
    output:
        unmapped='output/{dataset}/mapping/{sample_id}/unmapped/tucpRNA.fa.gz',
        bam=[ 'output/{dataset}/mapping/{sample_id}/bam/tucpRNA.bam' + sequence + ".bam" for sequence in sequences ]
    log:
        "output/{dataset}/log/mapping/{sample_id}/tucpRNA.txt"
    params:
        index=genome_dir + '/index/tucpRNA'
    threads: 
        4
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --norc --sensitive --no-unal \
            --un-gz {output.unmapped} -x {params.index} - -S - > {log}\
        | samtools view -b -o {output.bam}
        '''

rule map_Y_RNA:
    input:
        reads='output/{dataset}/mapping/{sample_id}/unmapped/tucpRNA.fa.gz',
        index=genome_dir + '/index/Y_RNA.1.bt2'
    output:
        unmapped='output/{dataset}/mapping/{sample_id}/unmapped/Y_RNA.fa.gz',
        bam=[ 'output/{dataset}/mapping/{sample_id}/bam/Y_RNA.bam' + sequence + ".bam" for sequence in sequences ]
    log:
        "output/{dataset}/log/mapping/{sample_id}/Y_RNA.txt"
    params:
        index=genome_dir + '/index/Y_RNA'
    threads: 
        config['threads_mapping']
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --norc --sensitive --no-unal \
            --un-gz {output.unmapped} -x {params.index} - -S - \
        | samtools view -b -o {output.bam}
        '''

rule map_circRNA:
    input:
        reads='output/{dataset}/mapping/{sample_id}/unmapped/Y_RNA.fa.gz',
        index=genome_dir + '/index/circRNA.1.bt2'
    output:
        unmapped='output/{dataset}/mapping/{sample_id}/unmapped/circRNA.fa.gz',
        bam=[ 'output/{dataset}/mapping/{sample_id}/bam/circRNA.bam' + sequence + ".bam" for sequence in sequences ]
    log:
        "output/{dataset}/log/mapping/{sample_id}/circ_RNA.txt"
    params:
        index=genome_dir + '/index/circRNA'
    threads: 
        4
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --norc --sensitive --no-unal \
            --un-gz {output.unmapped} -x {params.index} - -S - > {log}\
        | samtools view -b -o {output.bam}
        '''

rule map_other:
    input:
        reads='output/{dataset}/mapping/{sample_id}/unmapped/circRNA.fa.gz',
        index=genome_dir + '/index/genome.1.bt2'
    output:
        unmapped='output/{dataset}/mapping/{sample_id}/unmapped/other.fa.gz',
        bam=[ 'output/{dataset}/mapping/{sample_id}/bam/other.bam' + sequence + ".bam" for sequence in sequences ]
    log:
        "output/{dataset}/log/mapping/{sample_id}/other.txt"
    params:
        index=genome_dir + '/index/genome'
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --sensitive --no-unal \
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
