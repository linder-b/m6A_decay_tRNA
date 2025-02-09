rule kallisto_index:
    input:
        genome_dir+"/{species}_{build}_{release}_{spikein}/transcripts.fa"
    output:
        genome_dir+"/{species}_{build}_{release}_{spikein}/transcriptome_kallisto.idx"
    log:
        "logs/kallisto/{species}_{build}_{release}_{spikein}/kallisto_index.log"
    conda: 'envs/sm_kallisto.yaml'
    threads: 12
    shell:
        '''
        kallisto index \
        --index={output} \
        --kmer-size=5 \
        {input} \
        > {log}
        '''

rule kallisto_quant_SE:
    input:
        fq1=temp_dir+"/flexbar/{sample}.fastq",
        index = genome_dir+"/{species}_{build}_{release}_{spikein}/transcriptome_kallisto.idx"
    output:
        directory('results/kallisto/{species}_{build}_{release}_{spikein}/quant_results_{sample}')
    log:
        "logs/kallisto/{species}_{build}_{release}_{spikein}/kallisto_quant_{sample}.log"
    conda: 'envs/sm_kallisto.yaml'
    threads: 24
    shell:
        '''
        kallisto quant \
        --threads={threads} \
        --index={input.index} \
        --single \
        -l 180 \
        -s 20 \
        --output-dir={output} \
        {input.fq1} \
        > {log}
        '''