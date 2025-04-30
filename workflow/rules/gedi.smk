
rule gedi_genome:
    input:
        fasta=genome_dir+"/{species}_{build}_{release}/genome.fasta",
        gtf=genome_dir+"/{species}_{build}_{release}/annotation.gtf"
    output:
        genome_dir+"/{species}_{build}_{release}/gedi/{species}_{build}_{release}.oml"
    params:
        name="{species}_{build}_{release}"
    threads: 12
    conda:
        '../envs/sm_gedi.yaml'
    shell:
        '''
        export PATH="$PATH:workflow/resources/bin/GRAND-SLAM_2.0.5f" &&
        gedi \
        -e IndexGenome \
        -ignoreMulti \
        -nokallisto \
        -nobowtie \
        -nostar \
        -s {input.fasta} \
        -a {input.gtf} \
        -n {params.name} \
        -o {output}
        '''

rule gedi_bamlist:
    input:
        expand_runs(results_dir+'/{experiment}/star/{species}_{build}_{release}/{run}.Aligned.sortedByCoord.out.bam', runs_slamseq)
    output:
        bamlist=results_dir+'/{experiment}/gedi/{species}_{build}_{release}/slam_input.bamlist'
    run:
        f=open(str(list({output.bamlist})[0]),'w')
        for b in list({input})[0]:
            print(b+'\n')
            f.write(b+'\n')
        f.close()

rule gedi_bam2cit:
    input:
        bamlist=results_dir+'/{experiment}/gedi/{species}_{build}_{release}/slam_input.bamlist'
    output:
        cit=results_dir+'/{experiment}/gedi/{species}_{build}_{release}/slam_input.bamlist.cit',
    params:
        tmpdir=tmp_dir
    threads: 24
    conda:
        '../envs/sm_gedi.yaml'
    shell:
        '''
        export PATH="$PATH:workflow/resources/bin/GRAND-SLAM_2.0.5f" &&
        bamlist2cit \
        -t {params.tmpdir} \
        {input.bamlist}
        '''

rule gedi_grandslam:
    wildcard_constraints:
        run='|'.join(samples.query(f" type == 'slamseq' ")['run'].tolist()+["dummy_constraint"])
    input:
        cit=expand_runs(results_dir+'/{experiment}/gedi/{species}_{build}_{release}/slam_input.bamlist.cit', runs_slamseq),
        bamindex=expand_runs(results_dir+'/{experiment}/star/{species}_{build}_{release}/{run}.Aligned.sortedByCoord.out.bam.bai', runs_slamseq),
        oml=genome_dir+"/{species}_{build}_{release}/gedi/{species}_{build}_{release}.oml"
    output:
        tsv=results_dir+'/{experiment}/gedi/{species}_{build}_{release}/slam.tsv'
    params:
        prefix=results_dir+'/{experiment}/gedi/{species}_{build}_{release}/slam',
        config= lambda wildcards: config['experiments'][wildcards.experiment]['gedi_grandslam']
    threads: 24
    conda:
        '../envs/sm_gedi.yaml'
    shell:
        '''
        export PATH="$PATH:workflow/resources/bin/GRAND-SLAM_2.0.5f"
        gedi \
        -e Slam \
        -nthreads {threads} \
        -genomic {input.oml} \
        -reads {input.cit} \
        -prefix {params.prefix} \
        {params.config}
        '''
