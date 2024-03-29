# cfRNA-SEEK
```
rule map_spikein:
    input:
        reads='{output_dir}/unmapped/{sample_id}/clean.fa.gz',
        index=genome_dir + '/index/bowtie2/spikein_small.1.bt2'
    output:
        unmapped='{output_dir}/unmapped/{sample_id}/spikein.fa.gz',
        bam='{output_dir}/tbam/{sample_id}/spikein.bam'
    params:
        index=genome_dir + '/index/bowtie2/spikein_small'
    threads: 
        config['threads_mapping']
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --norc --sensitive --no-unal \
            --un-gz {output.unmapped} -x {params.index} - -S - \
        | samtools view -b -o {output.bam}
        '''

rule map_univec:
    input:
        reads='{output_dir}/unmapped/{sample_id}/spikein.fa.gz',
        index=genome_dir + '/rsem_index/bowtie2/univec.1.bt2'
    output:
        unmapped='{output_dir}/unmapped/{sample_id}/univec.fa.gz',
        bam='{output_dir}/tbam/{sample_id}/univec.bam'
    params:
        index=genome_dir + '/rsem_index/bowtie2/univec'
    threads: 
        config['threads_mapping']
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --norc --sensitive --no-unal \
            --un-gz {output.unmapped} -x {params.index} - -S - \
        | samtools view -b -o {output.bam}
        '''

rule map_rRNA:
    input:
        reads='{output_dir}/unmapped/{sample_id}/univec.fa.gz',
        index=genome_dir + '/index/bowtie2/rRNA.1.bt2'
    output:
        unmapped='{output_dir}/unmapped/{sample_id}/rRNA.fa.gz',
        bam='{output_dir}/tbam/{sample_id}/rRNA.bam'
    params:
        index=genome_dir + '/index/bowtie2/rRNA'
    threads: 
        config['threads_mapping']
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --norc --sensitive --no-unal \
            --un-gz {output.unmapped} -x {params.index} - -S - \
        | samtools view -b -o {output.bam}
        '''

rule map_lncRNA:
    input:
        reads='{output_dir}/unmapped/{sample_id}/rRNA.fa.gz',
        index=genome_dir + '/rsem_index/bowtie2/lncRNA.1.bt2'
    output:
        unmapped='{output_dir}/unmapped/{sample_id}/lncRNA.fa.gz',
        bam='{output_dir}/tbam/{sample_id}/lncRNA.bam'
    params:
        index=genome_dir + '/rsem_index/bowtie2/lncRNA'
    threads: 
        config['threads_mapping']
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --norc --sensitive --no-unal \
            --un-gz {output.unmapped} -x {params.index} - -S - \
        | samtools view -b -o {output.bam}
        '''

rule map_miRNA:
    input:
        reads='{output_dir}/unmapped/{sample_id}/lncRNA.fa.gz',
        index=genome_dir + '/index/bowtie2/miRNA.1.bt2'
    output:
        unmapped='{output_dir}/unmapped/{sample_id}/miRNA.fa.gz',
        bam='{output_dir}/tbam/{sample_id}/miRNA.bam'
    params:
        index=genome_dir + '/index/bowtie2/miRNA'
    threads: 
        config['threads_mapping']
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --norc --sensitive --no-unal \
            --un-gz {output.unmapped} -x {params.index} - -S - \
        | samtools view -b -o {output.bam}
        '''

rule map_mRNA:
    input:
        reads='{output_dir}/unmapped/{sample_id}/miRNA.fa.gz',
        index=genome_dir + '/rsem_index/bowtie2/mRNA.1.bt2'
    output:
        unmapped='{output_dir}/unmapped/{sample_id}/mRNA.fa.gz',
        bam='{output_dir}/tbam/{sample_id}/mRNA.bam'
    params:
        index=genome_dir + '/rsem_index/bowtie2/mRNA'
    threads: 
        config['threads_mapping']
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --norc --sensitive --no-unal \
            --un-gz {output.unmapped} -x {params.index} - -S - \
        | samtools view -b -o {output.bam}
        '''

rule map_piRNA:
    input:
        reads='{output_dir}/unmapped/{sample_id}/mRNA.fa.gz',
        index=genome_dir + '/rsem_index/bowtie2/piRNA.1.bt2'
    output:
        unmapped='{output_dir}/unmapped/{sample_id}/piRNA.fa.gz',
        bam='{output_dir}/tbam/{sample_id}/piRNA.bam'
    params:
        index=genome_dir + '/rsem_index/bowtie2/piRNA'
    threads: 
        config['threads_mapping']
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --norc --sensitive --no-unal \
            --un-gz {output.unmapped} -x {params.index} - -S - \
        | samtools view -b -o {output.bam}
        '''

rule map_snoRNA:
    input:
        reads='{output_dir}/unmapped/{sample_id}/piRNA.fa.gz',
        index=genome_dir + '/rsem_index/bowtie2/snoRNA.1.bt2'
    output:
        unmapped='{output_dir}/unmapped/{sample_id}/snoRNA.fa.gz',
        bam='{output_dir}/tbam/{sample_id}/snoRNA.bam'
    params:
        index=genome_dir + '/rsem_index/bowtie2/snoRNA'
    threads: 
        config['threads_mapping']
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --norc --sensitive --no-unal \
            --un-gz {output.unmapped} -x {params.index} - -S - \
        | samtools view -b -o {output.bam}
        '''

rule map_snRNA:
    input:
        reads='{output_dir}/unmapped/{sample_id}/snoRNA.fa.gz',
        index=genome_dir + '/rsem_index/bowtie2/snRNA.1.bt2'
    output:
        unmapped='{output_dir}/unmapped/{sample_id}/snRNA.fa.gz',
        bam='{output_dir}/tbam/{sample_id}/snRNA.bam'
    params:
        index=genome_dir + '/rsem_index/bowtie2/snRNA'
    threads: 
        config['threads_mapping']
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --norc --sensitive --no-unal \
            --un-gz {output.unmapped} -x {params.index} - -S - \
        | samtools view -b -o {output.bam}
        '''

rule map_srpRNA:
    input:
        reads='{output_dir}/unmapped/{sample_id}/snRNA.fa.gz',
        index=genome_dir + '/rsem_index/bowtie2/srpRNA.1.bt2'
    output:
        unmapped='{output_dir}/unmapped/{sample_id}/srpRNA.fa.gz',
        bam='{output_dir}/tbam/{sample_id}/srpRNA.bam'
    params:
        index=genome_dir + '/rsem_index/bowtie2/srpRNA'
    threads: 
        config['threads_mapping']
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --norc --sensitive --no-unal \
            --un-gz {output.unmapped} -x {params.index} - -S - \
        | samtools view -b -o {output.bam}
        '''

rule map_tRNA:
    input:
        reads='{output_dir}/unmapped/{sample_id}/srpRNA.fa.gz',
        index=genome_dir + '/rsem_index/bowtie2/tRNA.1.bt2'
    output:
        unmapped='{output_dir}/unmapped/{sample_id}/tRNA.fa.gz',
        bam='{output_dir}/tbam/{sample_id}/tRNA.bam'
    params:
        index=genome_dir + '/rsem_index/bowtie2/tRNA'
    threads: 
        config['threads_mapping']
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --norc --sensitive --no-unal \
            --un-gz {output.unmapped} -x {params.index} - -S - \
        | samtools view -b -o {output.bam}
        '''

rule map_tucpRNA:
    input:
        reads='{output_dir}/unmapped/{sample_id}/tRNA.fa.gz',
        index=genome_dir + '/rsem_index/bowtie2/tucpRNA.1.bt2'
    output:
        unmapped='{output_dir}/unmapped/{sample_id}/tucpRNA.fa.gz',
        bam='{output_dir}/tbam/{sample_id}/tucpRNA.bam'
    params:
        index=genome_dir + '/rsem_index/bowtie2/tucpRNA'
    threads: 
        config['threads_mapping']
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --norc --sensitive --no-unal \
            --un-gz {output.unmapped} -x {params.index} - -S - \
        | samtools view -b -o {output.bam}
        '''

rule map_Y_RNA:
    input:
        reads='{output_dir}/unmapped/{sample_id}/tucpRNA.fa.gz',
        index=genome_dir + '/rsem_index/bowtie2/Y_RNA.1.bt2'
    output:
        unmapped='{output_dir}/unmapped/{sample_id}/Y_RNA.fa.gz',
        bam='{output_dir}/tbam/{sample_id}/Y_RNA.bam'
    params:
        index=genome_dir + '/rsem_index/bowtie2/Y_RNA'
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
        reads='{output_dir}/unmapped/{sample_id}/Y_RNA.fa.gz',
        index=genome_dir + '/index/bowtie2/circRNA.1.bt2'
    output:
        unmapped='{output_dir}/unmapped/{sample_id}/circRNA.fa.gz',
        bam='{output_dir}/tbam/{sample_id}/circRNA.bam'
    params:
        index=genome_dir + '/index/bowtie2/circRNA'
    threads: 
        config['threads_mapping']
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --norc --sensitive --no-unal \
            --un-gz {output.unmapped} -x {params.index} - -S - \
        | samtools view -b -o {output.bam}
        '''

rule map_other:
    input:
        reads='{output_dir}/unmapped/{sample_id}/circRNA.fa.gz',
        index=genome_dir + '/genome_index/bowtie2/genome.1.bt2'
    output:
        unmapped='{output_dir}/unmapped/{sample_id}/other.fa.gz',
        bam='{output_dir}/gbam/{sample_id}/other.bam'
    params:
        index=genome_dir + '/genome_index/bowtie2/genome'
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --sensitive --no-unal \
            --un-gz {output.unmapped} -x {params.index} - -S - \
        | samtools view -b -o {output.bam}
        '''
```
