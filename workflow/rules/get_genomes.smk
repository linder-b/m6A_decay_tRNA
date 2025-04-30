rule download_genome:
    """\
    Download genome FASTA file.
    """
    params:
        url=lambda wildcards: config['genomes'][wildcards.species]['fasta_url']
    output:
        genome_dir+"/{species}_{build}_{release}/genome.fasta"
    shell:
        """\
        wget -O - {params.url} | gunzip -d -c > {output}
        """

rule download_annotation:
    """\
    Download GTF file & fix biotype field from "gene_type" to "gene_biotype".
    """
    params:
        url=lambda wildcards: config['genomes'][wildcards.species]['gtf_url']
    output:
        genome_dir+"/{species}_{build}_{release}/annotation.gtf"
    shell:
        """\
        wget -O - {params.url} | gunzip -d -c > {output} &&
        sed -i 's/gene_type/gene_biotype/g' {output} &&
        sed -i 's/transcript_type/transcript_biotype/g' {output}
        """

rule download_transcripts:
    """\
    Download transcript FASTA file.
    """
    params:
        url=lambda wildcards: config['genomes'][wildcards.species]['tx_fasta_url']
    output:
        genome_dir+"/{species}_{build}_{release}/transcripts.fasta"
    shell:
        """\
        wget -O - {params.url} | gunzip -d -c > {output}
        """



