rule flexbar_SE:
    wildcard_constraints:
        run=constraints_SE
    input:
        tmp_dir+"/{experiment}/umitools/{run}.fastq"
    output:
        temp(tmp_dir+"/{experiment}/flexbar/{run}.fastq")
    params:
        config=lambda wildcards: get_params(wildcards.run, 'flexbar')
    log:
        logs_dir+'/{experiment}/flexbar/{run}.log'
    params:
        config=config['defaults']['flexbar']
    threads: 3
    conda:
        '../envs/sm_flexbar.yaml'
    shell:
        '''
        flexbar \
        --threads {threads} \
        --reads  {input} \
        --output-reads  {output} \
        --output-log {log} \
        {params.config}
        '''

rule flexbar_PE:
    wildcard_constraints:
        run=constraints_PE
    input:
        r1=tmp_dir+"/{experiment}/umitools/{run}_1.fastq",
        r2=tmp_dir+"/{experiment}/umitools/{run}_2.fastq"
    output:
        r1=temp(tmp_dir+"/{experiment}/flexbar/{run}_1.fastq"),
        r2=temp(tmp_dir+"/{experiment}/flexbar/{run}_2.fastq")
    params:
        config=lambda wildcards: get_params(wildcards.run, 'flexbar')
    log:
        logs_dir+'/{experiment}/flexbar/{run}.log'
    threads: 3
    conda:
        '../envs/sm_flexbar.yaml'
    shell:
        '''
        flexbar \
        --threads {threads} \
        --reads  {input.r1} \
        --reads2 {input.r2} \
        --output-reads  {output.r1} \
        --output-reads2 {output.r2} \
        --output-log {log} \
        {params.config}
        '''