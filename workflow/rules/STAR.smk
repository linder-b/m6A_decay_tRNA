rule STAR_index:
    input:
        fasta=genome_dir+'/{species}_{build}_{release}_{spikein}/genome.fasta',
        gtf=genome_dir+'/{species}_{build}_{release}_{spikein}/annotation.gtf'
    output:
        directory(genome_dir+'/{species}_{build}_{release}_{spikein}/STAR_INDEX')
    threads: 24
    conda: '../envs/sm_preprocess.yaml'
    shell:
        """
        mkdir -p {output};\
        STAR\
        --runThreadN {threads}\
        --runMode genomeGenerate\
        --genomeDir {output}\
        --genomeFastaFiles {input.fasta}\
        --sjdbGTFfile {input.gtf}
        """
        
rule STAR_SE:
    wildcard_constraints:
        run=constraints_SE
    input:
        idx=genome_dir+'/{species}_{build}_{release}_{spikein}/STAR_INDEX',
        r1=tmp_dir+'/{experiment}/flexbar/{run}.fastq'
    output:
        results_dir+'/{experiment}/star/{species}_{build}_{release}_{spikein}/{run}.Aligned.sortedByCoord.out.bam'
    log:
        logs_dir+'/{experiment}/star/{species}_{build}_{release}_{spikein}/{run}.Log.final.out'
    params:
        prefix=results_dir+'/{experiment}/star/{species}_{build}_{release}_{spikein}/{run}.', ## trailing dot required
        tmpbase=tmp_dir+'/{experiment}/star/{species}_{build}_{release}_{spikein}',
        tmpdir=tmp_dir+'/{experiment}/star/{species}_{build}_{release}_{spikein}/{run}',
        logfile=results_dir+'/{experiment}/star/{species}_{build}_{release}_{spikein}/{run}.Log.final.out',
        config=lambda wildcards: get_params(wildcards.run, 'STAR')
    threads: 12
    conda:
        '../envs/sm_preprocess.yaml'
    shell:
        '''
        mkdir -p {params.tmpbase} &&
        STAR \
        --runThreadN {threads} \
        --genomeDir {input.idx} \
        --readFilesIn {input.r1} \
        --outFileNamePrefix {params.prefix} \
        --outTmpDir {params.tmpdir} \
        {params.config} &&
        mv {params.logfile} {log}
        '''

rule STAR_PE:
    wildcard_constraints:
        run=constraints_PE
    input:
        idx=genome_dir+'/{species}_{build}_{release}_{spikein}/STAR_INDEX',
        r1=tmp_dir+"/{experiment}/flexbar/{run}_1.fastq",
        r2=tmp_dir+"/{experiment}/flexbar/{run}_2.fastq"
    output:
        results_dir+'/{experiment}/star/{species}_{build}_{release}_{spikein}/{run}.Aligned.sortedByCoord.out.bam'
    log:
        logs_dir+'/{experiment}/star/{species}_{build}_{release}_{spikein}/{run}.Log.final.out'
    params:
        prefix=results_dir+'/{experiment}/star/{species}_{build}_{release}_{spikein}/{run}.', ## trailing dot required
        tmpbase=tmp_dir+'/star/{species}_{build}_{release}_{spikein}',
        tmpdir=tmp_dir+'/star/{species}_{build}_{release}_{spikein}/{run}',
        logfile=results_dir+'/{experiment}/star/{species}_{build}_{release}_{spikein}/{run}.Log.final.out',
        config=lambda wildcards: get_params(wildcards.run, 'STAR')
    threads: 12
    conda:
        '../envs/sm_preprocess.yaml'
    shell:
        '''
        mkdir -p {params.tmpbase} &&
        STAR \
        --runThreadN {threads} \
        --genomeDir {input.idx} \
        --readFilesIn {input.r1} {input.r2}\
        --outFileNamePrefix {params.prefix} \
        --outTmpDir {params.tmpdir} \
        {params.config} &&
        mv {params.logfile} {log}
        '''

rule index_bam:
    input:
        bam=results_dir+'/{experiment}/star/{species}_{build}_{release}_{spikein}/{run}.Aligned.sortedByCoord.out.bam'
    output:
        idx=results_dir+'/{experiment}/star/{species}_{build}_{release}_{spikein}/{run}.Aligned.sortedByCoord.out.bam.bai'
    conda:
        '../envs/sm_preprocess.yaml'
    shell:
        '''
        samtools index {input.bam}
        '''