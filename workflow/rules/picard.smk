
# --- generate metadata --- #

rule picard_dict:
    input:
        fasta=genome_dir+"/{species}_{build}_{release}/genome.fasta"
    output:
        dict=genome_dir+"/{species}_{build}_{release}/genome.dict"
    threads:1
    params:
        tmp_directory=tmp_dir
    conda: '../envs/sm_picard.yaml'
    shell:
        """
        picard \
        CreateSequenceDictionary \
        -REFERENCE {input} \
        -OUTPUT {output}
        """


rule dropseq_refFlat:
    """
    Generate refFlat file. Does not include spike-ins.
    """
    input:
        annotation=genome_dir+"/{species}_{build}_{release}/annotation.gtf",
        reference_dict=genome_dir+"/{species}_{build}_{release}/genome.dict"
    params:
        memory=config['local']['memory'],
        tmp_directory=tmp_dir
    output:
        flat=genome_dir+"/{species}_{build}_{release}/annotation.refFlat"
    conda: '../envs/sm_dropseq_tools.yaml'
    shell:
        """
        export _JAVA_OPTIONS=-Djava.io.tmpdir={params.tmp_directory} &&
        ConvertToRefFlat -m {params.memory} \
        ANNOTATIONS_FILE={input.annotation} \
        OUTPUT={output.flat} \
        SEQUENCE_DICTIONARY={input.reference_dict}
        """

rule curate_annotation:
    input:
        # biotypes=config['META']['gtf_biotypes'],
        #biotypes=['rRNA'],
        annotation=genome_dir+"/{species}_{build}_{release}/annotation.gtf"
    output:
        temp(genome_dir+"/{species}_{build}_{release}/curated_annotation.gtf")
    params:
        patterns='|'.join('rRNA')
    shell:
        """cat {input.annotation} | grep -E "{params.patterns}" > {output}"""

rule dropseq_reduce_gtf:
    input:
        reference_dict=genome_dir+"/{species}_{build}_{release}/genome.dict",
        annotation=genome_dir+"/{species}_{build}_{release}/curated_annotation.gtf"
    params:
        memory=config['local']['memory'],
        tmp_directory=tmp_dir
    output:
        genome_dir+"/{species}_{build}_{release}/curated_reduced_annotation.gtf"
    conda: '../envs/sm_dropseq_tools.yaml'
    shell:
        """export _JAVA_OPTIONS=-Djava.io.tmpdir={params.tmp_directory} && ReduceGtf -m {params.memory}\
        GTF={input.annotation}\
        OUTPUT={output}\
        SEQUENCE_DICTIONARY={input.reference_dict}\
        IGNORE_FUNC_TYPE='null'\
        ENHANCE_GTF='false'"""

rule dropseq_create_intervals:
    input:
        annotation_reduced=genome_dir+"/{species}_{build}_{release}/curated_reduced_annotation.gtf",
        reference_dict=genome_dir+"/{species}_{build}_{release}/genome.dict"
    params:
        memory=config['local']['memory'],
        reference_directory=genome_dir+"/{species}_{build}_{release}",
        tmp_directory=tmp_dir,
        prefix="annotation"
    output:
        intervals=genome_dir+"/{species}_{build}_{release}/annotation.rRNA.intervals"
    conda: '../envs/sm_dropseq_tools.yaml'
    shell:
        """export _JAVA_OPTIONS=-Djava.io.tmpdir={params.tmp_directory} && CreateIntervalsFiles -m {params.memory}\
        REDUCED_GTF={input.annotation_reduced}\
        SEQUENCE_DICTIONARY={input.reference_dict}\
        O={params.reference_directory}\
        PREFIX={params.prefix}
        """

## collect metrics

rule picard_CollectRnaSeqMetrics:
    input:
        bam=results_dir+'/{experiment}/star/{species}_{build}_{release}/{run}.Aligned.sortedByCoord.out.bam',
        flat=genome_dir+"/{species}_{build}_{release}/annotation.refFlat",
        rRNA_intervals=genome_dir+"/{species}_{build}_{release}/annotation.rRNA.intervals"
    params:     
        tmp_directory=tmp_dir,
        memory=config['local']['memory']
    output:
        metrics=logs_dir+'/{experiment}/picard/{species}_{build}_{release}/{run}.rna_metrics.txt'
    conda: '../envs/sm_picard.yaml'
    shell:
        """
        picard \
        CollectRnaSeqMetrics \
        -INPUT {input.bam} \
        -OUTPUT {output.metrics} \
        -STRAND SECOND_READ_TRANSCRIPTION_STRAND \
        -REF_FLAT {input.flat} \
        -RIBOSOMAL_INTERVALS {input.rRNA_intervals}
        """
        
rule picard_CollectInsertSizeMetrics:
    wildcard_constraints:
        sample='|'.join(runs_PE)
    input:
        bam=results_dir+'/{experiment}/star/{species}_{build}_{release}/{run}.Aligned.sortedByCoord.out.bam',
        idx=results_dir+'/{experiment}/star/{species}_{build}_{release}/{run}.Aligned.sortedByCoord.out.bam.bai'
    output:
        txt=logs_dir+'/{experiment}/picard/{species}_{build}_{release}/{run}.isize.txt',
        pdf=logs_dir+'/{experiment}/picard/{species}_{build}_{release}/{run}.isize.pdf'
    conda: '../envs/sm_picard.yaml'
    shell:
        #"""
        #picard \
        #CollectInsertSizeMetrics \
        #-INPUT {input.bam} \
        #-OUTPUT {output.txt} \
        #-HISTOGRAM_FILE {output.pdf}
        #"""
        """
        picard CollectInsertSizeMetrics \
        INPUT={input.bam} \
        OUTPUT={output.txt} \
        HISTOGRAM_FILE={output.pdf}
        """
