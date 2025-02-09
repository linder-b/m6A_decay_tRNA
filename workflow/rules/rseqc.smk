rule rseqc_gtf2bed:
    input:
        gtf=genome_dir+"/{species}_{build}_{release}/annotation.gtf"
    output:
        bed=genome_dir+"/{species}_{build}_{release}/annotation.protein_coding.bed"
    conda:
        '../envs/sm_rseqc.yaml'
    shell:
        """
        cat {input.gtf} | grep "protein_coding" | bedparse gtf2bed - > {output.bed}
        """
        
rule rseqc_infer_experiment:
    input:
        bam=results_dir+'/{experiment}/star/{species}_{build}_{release}_{spikein}/{run}.Aligned.sortedByCoord.out.bam',
        bed=genome_dir+"/{species}_{build}_{release}/annotation.protein_coding.bed"
    output:
        log=logs_dir+'/{experiment}/rseqc/{species}_{build}_{release}_{spikein}/{run}.infer_experiment.txt'
    conda:
        '../envs/sm_rseqc.yaml'
    shell:
        """
        infer_experiment.py -r {input.bed} -i {input.bam} -s 1000000 > {output.log}
        """

rule rseqc_mismatch_profile:
    input:
        bam=results_dir+'/{experiment}/star/{species}_{build}_{release}_{spikein}/{run}.Aligned.sortedByCoord.out.bam',
        bed=genome_dir+"/{species}_{build}_{release}/annotation.protein_coding.bed"
    output:
        log=logs_dir+'/{experiment}/rseqc/{species}_{build}_{release}_{spikein}/{run}.mismatch_profile.r'
    params:
        prefix=logs_dir+'/{experiment}/rseqc/{species}_{build}_{release}_{spikein}/{run}'
    conda:
        '../envs/sm_rseqc.yaml'
    shell:
        """
        mismatch_profile.py -l 43 -i {input.bam} -o {params.prefix}
        """
        