ruleorder: fastqc_PE > fastqc_SE


rule fastqc_SE:
    wildcard_constraints:
        run=constraints_SE
    input:
        r1=fastq_dir+'/{experiment}/{run}.fastq.gz'
    output:
        html=logs_dir+'/{experiment}/fastqc/{run}_fastqc.html',
        zip=logs_dir+'/{experiment}/fastqc/{run}_fastqc.zip'
    params:
        outDir=logs_dir+'/{experiment}/fastqc'
    threads: 6
    conda:
        '../envs/sm_preprocess.yaml'
    shell:
        '''
        fastqc \
        --quiet \
        --outdir  {params.outDir} \
        --threads {threads} \
        --adapters workflow/resources/fastqc/adapter_list.txt \
        --contaminants workflow/resources/fastqc/contaminant_list.txt \
        {input.r1}
        '''

rule fastqc_PE:
    wildcard_constraints:
        run=constraints_PE
    input:
        r1=fastq_dir+'/{experiment}/{run}_1.fastq.gz'
    output:
        html=logs_dir+'/{experiment}/fastqc/{run}_1_fastqc.html',
        zip=logs_dir+'/{experiment}/fastqc/{run}_1_fastqc.zip'
    params:
        outDir=logs_dir+'/{experiment}/fastqc'
    threads: 6
    conda:
        '../envs/sm_preprocess.yaml'
    shell:
        '''
        fastqc \
        --quiet \
        --outdir  {params.outDir} \
        --threads {threads} \
        --adapters workflow/resources/fastqc/adapter_list.txt \
        --contaminants workflow/resources/fastqc/contaminant_list.txt \
        {input.r1}
        '''
