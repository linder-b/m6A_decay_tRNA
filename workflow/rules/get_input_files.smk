# --- download from ENA --- #
def get_ena_path(run):
    """\
    Get the path for FTP downloads from ENA using an SRR accession.
    """
    path=f"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{run[0:6]}/00{run[-1]}/{run}/"
    #path=f"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{run[0:6]}/{run}/"
    return(path)

    
rule download_ena_ftp_SE:
    """\
    Download single end fastq file from ENA.
    """
    wildcard_constraints:
        run='|'.join(samples.query(f" source == 'ENA' and layout == 'SINGLE' ")['run'].tolist()+["dummy_constraint"])
    params:
        lambda wildcards: get_ena_path(wildcards.run)+samples.loc[wildcards.run, 'read1']
    output:
        fastq_dir+"/{experiment}/{run}.fastq.gz"
    shell:
        """\
        wget -O - {params} > {output}
        """

rule download_ena_ftp_PE:
    """\
    Download paired end fastq files from ENA.
    """
    wildcard_constraints:
        run='|'.join(samples.query(f" source == 'ENA' and layout == 'PAIRED' ")['run'].tolist()+["dummy_constraint"])
    params:
        r1=lambda wildcards: get_ena_path(wildcards.run)+samples.loc[wildcards.run, 'read1'],
        r2=lambda wildcards: get_ena_path(wildcards.run)+samples.loc[wildcards.run, 'read2']
    output:
        r1=fastq_dir+"/{experiment}/{run}_1.fastq.gz",
        r2=fastq_dir+"/{experiment}/{run}_2.fastq.gz"
    shell:
        """\
        wget -O - {params.r1} > {output.r1} && \
        wget -O - {params.r2} > {output.r2}
        """

rule download_sra_SE:
    """\
    Download single end fastq files from SRA.
    """
    wildcard_constraints:
        run='|'.join(samples.query(f" source == 'SRA' and layout == 'SINGLE' ")['run'].tolist()+["dummy_constraint"])
    params:
        tmpdir=tmp_dir,
        outdir=fastq_dir+"/{experiment}/",
        r1=fastq_dir+"/{experiment}/{run}.fastq"
    output:
        r1=temp(fastq_dir+"/{experiment}/{run}.fastq.gz")
    conda: '../envs/sm_preprocess.yaml'
    shell:
        """\
        fasterq-dump {wildcards.run} -O {params.outdir} -t {params.tmpdir} && \
        gzip {params.r1}
        """

rule download_sra_PE:
    """\
    Download paired end fastq files from SRA.
    """
    wildcard_constraints:
        run='|'.join(samples.query(f" source == 'SRA' and layout == 'PAIRED' ")['run'].tolist()+["dummy_constraint"])
    params:
        tmpdir=tmp_dir,
        outdir=fastq_dir+"/{experiment}/",
        r1=fastq_dir+"/{experiment}/{run}_1.fastq",
        r2=fastq_dir+"/{experiment}/{run}_2.fastq"
    output:
        r1=temp(fastq_dir+"/{experiment}/{run}_1.fastq.gz"),
        r2=temp(fastq_dir+"/{experiment}/{run}_2.fastq.gz")
    conda: '../envs/sm_preprocess.yaml'
    shell:
        """\
        fasterq-dump {wildcards.run} -O {params.outdir} -t {params.tmpdir} && \
        gzip {params.r1} && \
        gzip {params.r2}
        """
