rule umitools_extract_SE:
    input:
        r1=fastq_dir+'/{experiment}/{run}.fastq.gz'
    output:
        r1=temp(tmp_dir+"/{experiment}/umitools/{run}.fastq")
    params:
        tmp=tmp_dir,
        config=lambda wildcards: get_params(wildcards.run, 'umitools_extract')
    log:
        logs_dir+'/{experiment}/umitools/{run}.extract.log'
    threads: 1
    conda:
        '../envs/sm_umitools.yaml'
    shell:
        '''
        umi_tools extract \
        {params.config} \
        --stdin={input.r1} \
        --log={log} \
        --temp-dir={params.tmp} \
        --stdout={output.r1}
        '''

rule umitools_extract_PE:
    input:
        r1=fastq_dir+'/{experiment}/{run}_1.fastq.gz',
        r2=fastq_dir+'/{experiment}/{run}_2.fastq.gz'
    output:
        r1=temp(tmp_dir+"/{experiment}/umitools/{run}_1.fastq"),
        r2=temp(tmp_dir+"/{experiment}/umitools/{run}_2.fastq")
    params:
        tmp=tmp_dir,
        config=lambda wildcards: get_params(wildcards.run, 'umitools_extract')
    log:
        logs_dir+'/{experiment}/umitools/{run}.extract.log'
    threads: 1
    conda:
        '../envs/sm_umitools.yaml'
    shell:
        '''
        umi_tools extract \
        {params.config} \
        --stdin={input.r1} \
        --read2-in={input.r2} \
        --log={log} \
        --stdout={output.r1} \
        --read2-out={output.r2} \
        --temp-dir={params.tmp}
        '''

rule umitools_dedup:
    input:
        bam=results_dir+'/{experiment}/star/{species}_{build}_{release}/{run}.Aligned.sortedByCoord.out.bam',
        idx=results_dir+'/{experiment}/star/{species}_{build}_{release}/{run}.Aligned.sortedByCoord.out.bam.bai'
    output:
        bam=results_dir+'/{experiment}/star/{species}_{build}_{release}/{run}.Aligned.sortedByCoord.out.dedup.bam'
    params:
        config=lambda wildcards: get_params(wildcards.run, 'umitools_dedup')
    log:
        logs_dir+'/{experiment}/umitools/{species}_{build}_{release}/{run}.dedup.log'
    threads: 1
    conda:
        '../envs/sm_umitools.yaml'
    shell:
        '''
        umi_tools dedup \
        -I {input.bam} \
        -L {log} \
        -S {output.bam} \
        {params.config}
        '''
rule index_dedup:
    input:
        bam=results_dir+'/{experiment}/star/{species}_{build}_{release}/{run}.Aligned.sortedByCoord.out.dedup.bam'
    output:
        idx=results_dir+'/{experiment}/star/{species}_{build}_{release}/{run}.Aligned.sortedByCoord.out.dedup.bam.bai'
    conda:
        '../envs/sm_preprocess.yaml'
    shell:
        '''
        samtools index {input.bam}
        '''