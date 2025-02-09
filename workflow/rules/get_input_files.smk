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

# --- custom downlaods --- #

rule get_fastq_mcm5s2U_slamseq_2020_08_11:
    """\
    A rule that gets the fastq files for the mcm5s2U_slamseq_2020-08-11 experiment.
    At the moment, this ist just a link to a local path.
    """
    wildcard_constraints:
        run='|'.join(samples.query(f" experiment == 'mcm5s2U_slam'")['run'].tolist()+["dummy_constraint"])
    input:
        r1=lambda wildcards: "/g/steinmetz/linder/projects/m6A/data/raw/leidel_mcm5s2U_slamseq_2020-08-11/fastq/"+samples.loc[wildcards.run, 'read1'],
        r2=lambda wildcards: "/g/steinmetz/linder/projects/m6A/data/raw/leidel_mcm5s2U_slamseq_2020-08-11/fastq/"+samples.loc[wildcards.run, 'read2']
    output:
        r1=fastq_dir+"/{experiment}/{run}_1.fastq.gz",
        r2=fastq_dir+"/{experiment}/{run}_2.fastq.gz"
    shell:
        """\
        ln -s {input.r1} {output.r1} && \
        ln -s {input.r2} {output.r2}
        """

rule get_fastq_mcm5s2U_riboseq_2020_06_09:
    """\
    A rule that gets the fastq files for the mcm5s2U_riboseq_2020-06-09 experiment.
    At the moment, this ist just a link to a local path.
    """
    wildcard_constraints:
        run='|'.join(samples.query(f" experiment == 'mcm5s2U_riboseq'")['run'].tolist()+["dummy_constraint"])
    input:
        r1=lambda wildcards: "/g/steinmetz/linder/projects/m6A/data/raw/leidel_mcm5s2U_riboseq_2020-06-09/fastq/"+samples.loc[wildcards.run, 'read1']
    output:
        r1=fastq_dir+"/{experiment}/{run}.fastq.gz"
    shell:
        """\
        ln -s {input.r1} {output.r1}
        """

rule get_fastq_mcm5s2U_rnaseq_2020_04_02:
    """\
    A rule that gets the fastq files for the mcm5s2U_rnaseq_2020-04-02 experiment.
    At the moment, this ist just a link to a local path.
    """
    wildcard_constraints:
        run='|'.join(samples.query(f" experiment == 'mcm5s2U_rnaseq'")['run'].tolist()+["dummy_constraint"])
    input:
        r1=lambda wildcards: "/g/steinmetz/linder/projects/m6A/data/raw/leidel_mcm5s2U_rnaseq_2020-04-02/fastq/"+samples.loc[wildcards.run, 'read1']
    output:
        r1=fastq_dir+"/{experiment}/{run}.fastq.gz"
    shell:
        """\
        ln -s {input.r1} {output.r1}
        """

rule get_fastq_mcm5s2U_quantseq_2021_07_12:
    """\
    A rule that gets the fastq files for the mcm5s2U_quantseq_2021-07-12 experiment.
    At the moment, this ist just a link to a local path.
    """
    wildcard_constraints:
        run='|'.join(samples.query(f" experiment == 'mcm5s2U_quantseq'")['run'].tolist()+["dummy_constraint"])
    input:
        r1=lambda wildcards: "/g/steinmetz/linder/projects/m6A/data/raw/leidel_mcm5s2U_quantseq_2021-07-12/fastq/"+samples.loc[wildcards.run, 'read1']
    output:
        r1=fastq_dir+"/{experiment}/{run}.fastq.gz"
    shell:
        """\
        ln -s {input.r1} {output.r1}
        """
rule get_fastq_mcm5s2U_quantseq_2021_09_08:
    """\
    A rule that gets the fastq files for the mcm5s2U_quantseq_2021-09-08 experiment.
    At the moment, this ist just a link to a local path.
    """
    wildcard_constraints:
        run='|'.join(samples.query(f" experiment == 'mcm5s2U_quantseq'")['run'].tolist()+["dummy_constraint"])
    input:
        r1=lambda wildcards: "/g/steinmetz/linder/projects/m6A/data/raw/leidel_mcm5s2U_quantseq_2021-09-08/fastq/"+samples.loc[wildcards.run, 'read1']
    output:
        r1=fastq_dir+"/{experiment}/{run}.fastq.gz"
    shell:
        """\
        ln -s {input.r1} {output.r1}
        """

rule get_fastq_mcm5s2U_STM_riboseq:
    """\
    A rule that gets the fastq files for the mcm5s2U_quantseq_2021-09-08 experiment.
    At the moment, this ist just a link to a local path.
    """
    wildcard_constraints:
        run='|'.join(samples.query(f" experiment == 'mcm5s2U_STM_riboseq'")['run'].tolist()+["dummy_constraint"])
    input:
        r1=lambda wildcards: "/g/steinmetz/linder/projects/m6A/data/raw/leidel_mcm5s2U_STM2457_riboseq_2023-08-17/fastq/"+samples.loc[wildcards.run, 'read1']
    output:
        r1=fastq_dir+"/{experiment}/{run}.fastq.gz"
    shell:
        """\
        ln -s {input.r1} {output.r1}
        """

rule get_mcm5s2U_actinomycinD_quantseq:
    """\
    A rule that gets the fastq files for the mcm5s2U_actinomycinD_quantseq_2023-12-20 experiment.
    At the moment, this ist just a link to a local path.
    """
    wildcard_constraints:
        run='|'.join(samples.query(f" experiment == 'mcm5s2U_actD'")['run'].tolist()+["dummy_constraint"])
    input:
        r1=lambda wildcards: "/g/steinmetz/linder/projects/m6A/data/raw/leidel_mcm5s2U_actinomycinD_quantseq_2023-12-20/fastq/"+samples.loc[wildcards.run, 'read1']
    output:
        r1=fastq_dir+"/{experiment}/{run}.fastq.gz"
    shell:
        """\
        ln -s {input.r1} {output.r1}
        """