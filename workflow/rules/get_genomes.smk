rule download_genome:
    """\
    Download genome FASTA file.
    """
    params:
        url=lambda wildcards: config['genomes'][wildcards.species]['fasta_url']
    output:
        genome_dir+"/{species}_{build}_{release}_nospike/genome.fasta"
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
        genome_dir+"/{species}_{build}_{release}_nospike/annotation.gtf"
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
        genome_dir+"/{species}_{build}_{release}_nospike/transcripts.fasta"
    shell:
        """\
        wget -O - {params.url} | gunzip -d -c > {output}
        """


# empty strings considered falsy (from: https://stackoverflow.com/a/9573259)
if any(list(set(expand_runs('{spikein}', runs = runs)) - set(['nospike']))):
    rule download_spikein:
        """\
        Download spike-in FASTA and GTF files and add a "gene_biotype" termed "spike-in".
        """
        params:
            url=lambda wildcards: config['genomes']['spikein'][wildcards.spikein.lstrip('_')]['url'] # we need to remove leading wildcard separator
        output:
            fasta=genome_dir+"/{spikein}/{spikein}.fasta",
            gtf=genome_dir+"/{spikein}/{spikein}.gtf",
        shell:
            """
            wget -O tmp.zip {params.url} &&
            unzip tmp.zip &&
            mv ERCC92.fa {output.fasta} &&
            mv ERCC92.gtf {output.gtf} &&
            awk '{{print $0,"gene_biotype \\"spike-in\\";"}}' {output.gtf} > {output.gtf}.tmp &&
            mv {output.gtf}.tmp {output.gtf} &&
            rm tmp.zip
            """

    rule insert_spikein:
        """\
        Concatenate genome and spike-in files.
        """
        wildcard_constraints:
            spikein='|'.join(samples.query(f" spikein != 'nospike' ")['spikein'].tolist()+["dummy_constraint"])
        input:
            fasta=genome_dir+"/{species}_{build}_{release}/genome.fasta",
            gtf=  genome_dir+"/{species}_{build}_{release}/annotation.gtf",
            spikein_fasta=genome_dir+"/{spikein}/{spikein}.fasta",
            spikein_gtf=genome_dir+"/{spikein}/{spikein}.gtf"
        output:
            fasta=genome_dir+"/{species}_{build}_{release}_{spikein}/genome.fasta",
            gtf=genome_dir+"/{species}_{build}_{release}_{spikein}/annotation.gtf"
        shell:
            '''
            cat {input.fasta} {input.spikein_fasta} > {output.fasta}
            cat {input.gtf} {input.spikein_gtf} > {output.gtf}
            '''


