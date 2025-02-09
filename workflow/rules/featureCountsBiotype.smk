rule featureCountsBiotype_SE:
    ## use featureCounts to summarize biotype counts for MultiQC. from: https://github.com/ewels/MultiQC/issues/1078#issuecomment-566802321
    wildcard_constraints:
       run=constraints_SE
    input:
        bam=results_dir+'/{experiment}/star/{species}_{build}_{release}_{spikein}/{run}.Aligned.sortedByCoord.out.bam',
        gtf=genome_dir+"/{species}_{build}_{release}_{spikein}/annotation.gtf"
    output:
        counts=logs_dir+'/{experiment}/featureCountsBiotype/{species}_{build}_{release}_{spikein}/{run}.featureCountsBiotype.txt',
        summary=logs_dir+'/{experiment}/featureCountsBiotype/{species}_{build}_{release}_{spikein}/{run}.featureCountsBiotype.txt.summary',
        mqc=logs_dir+'/{experiment}/featureCountsBiotype/{species}_{build}_{release}_{spikein}/{run}_biotype_counts_mqc.txt',
        gs_mqc=logs_dir+'/{experiment}/featureCountsBiotype/{species}_{build}_{release}_{spikein}/{run}_biotype_counts_gs_mqc.tsv'
    threads: 6
    conda:
        '../envs/sm_subread.yaml'
    shell:
        '''
        samtools view -h -s 42.1 {input.bam} | \
        featureCounts \
        -T {threads} \
        --tmpDir {tmp_dir} \
        -a {input.gtf} \
        -g gene_biotype \
        -M \
        --fraction \
        -s 2 \
        -o {output.counts} &&
        sed -i "2s/STDIN/{wildcards.run}/g" {output.counts} &&
        cut -f 1,7 {output.counts} \
            | tail -n +3 \
            | cat workflow/resources/featureCountsBiotype/biotypes_header.txt - \
            >> {output.mqc} && 
        workflow/scripts/featureCountsBiotype/mqc_features_stat.py \
            {output.mqc} \
            -s {wildcards.run} \
            -f protein_coding spike-in \
            -o {output.gs_mqc} &&
        sed -i "1s/STDIN/{wildcards.run}/g" {output.summary}
        '''

rule featureCountsBiotype_PE:
    ## use featureCounts to summarize biotype counts for MultiQC. from: https://github.com/ewels/MultiQC/issues/1078#issuecomment-566802321
    wildcard_constraints:
       run=constraints_PE
    input:
        bam=results_dir+'/{experiment}/star/{species}_{build}_{release}_{spikein}/{run}.Aligned.sortedByCoord.out.bam',
        gtf=genome_dir+"/{species}_{build}_{release}_{spikein}/annotation.gtf"
    output:
        counts=logs_dir+'/{experiment}/featureCountsBiotype/{species}_{build}_{release}_{spikein}/{run}.featureCountsBiotype.txt',
        summary=logs_dir+'/{experiment}/featureCountsBiotype/{species}_{build}_{release}_{spikein}/{run}.featureCountsBiotype.txt.summary',
        mqc=logs_dir+'/{experiment}/featureCountsBiotype/{species}_{build}_{release}_{spikein}/{run}_biotype_counts_mqc.txt',
        gs_mqc=logs_dir+'/{experiment}/featureCountsBiotype/{species}_{build}_{release}_{spikein}/{run}_biotype_counts_gs_mqc.tsv'
    threads: 6
    conda:
        '../envs/sm_subread.yaml'
    shell:
        '''
        samtools view -h -s 42.1 {input.bam} | \
        featureCounts \
        -T {threads} \
        --tmpDir {tmp_dir} \
        -a {input.gtf} \
        -g gene_biotype \
        -M \
        --fraction \
        -s 2 \
        -o {output.counts} \
        -p &&
        sed -i "2s/STDIN/{wildcards.run}/g" {output.counts} &&
        cut -f 1,7 {output.counts} \
        | tail -n +3 \
        | cat workflow/resources/featureCountsBiotype/biotypes_header.txt - \
        >> {output.mqc} && 
        workflow/scripts/featureCountsBiotype/mqc_features_stat.py \
        {output.mqc} \
        -s {wildcards.run} \
        -f protein_coding spike-in \
        -o {output.gs_mqc} &&
        sed -i "1s/STDIN/{wildcards.run}/g" {output.summary}
        '''