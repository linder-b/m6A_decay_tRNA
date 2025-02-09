rule htseq_count:
    input:
        gtf=genome_dir+"/{species}_{build}_{release}_{spikein}/annotation.gtf",
        bam=expand('results/star/{species}_{build}_{release}_{spikein}/{sample}.Aligned.sortedByCoord.out.bam', sample = runs, species = species, build = build, release = release, spikein = spikein)
    output:
        'results/htseq_count/{species}_{build}_{release}_{spikein}/counts.tsv'
    threads: 24
    conda:
        'envs/sm_htseq.yaml'
    shell:
        '''
                htseq-count \
                --nprocesses {threads} \
                --format bam \
                --order pos \
                --additional-attr gene_name \
                --mode intersection-strict \
                --counts_output {output} \
                {input.bam} {input.gtf}
        '''