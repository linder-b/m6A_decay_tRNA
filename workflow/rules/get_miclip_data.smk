

# --- B116 data --- #

rule download_B116:
    """\
    Download the tar archive for the lane(s) of the B116 miCLIP experiment.
    """
    params:
        outdir=tmp_dir,
        pwd="--password=PASSWORD"
    output:
        temp(tmp_dir+"/Sample_B116.tar")
    shell:
        '''\
        wget \
        -O {output} \
        --no-check-certificate \
        --content-disposition \
        --continue \
        --user=bal2008@med.cornell.edu \
        {params.pwd} \
        "https://abc.med.cornell.edu/pubshare/data/2/62986/Project_EC-BL-3508/Sample_B116/*.fastq.gz"
        '''

rule untar_B116:
    """\
    Extract the tar archive for the lane(s) of the B116 miCLIP experiment.
    """
    input:
        tmp_dir+"/Sample_B116.tar"
    output:
        r1=tmp_dir+"/Sample_B116/B116_S1_L006_R1_001.fastq.gz",
        r2=tmp_dir+"/Sample_B116/B116_S1_L006_R2_001.fastq.gz"
    params:
        outdir=tmp_dir
    shell:
        '''\
        tar -C {params.outdir} -xvf {input}
        '''

rule demux_B116:
    """\
    Demultiplex the B116 fastq file using the barcode fasta file provided in workflow/resources.
    Then remove target prefix (_barcode_).
    """
    wildcard_constraints:
        run='|'.join(samples.query(f" experiment == 'B116' ")['run'].tolist()+["dummy_constraint"])
    input:
        r1=tmp_dir+"/Sample_B116/B116_S1_L006_R1_001.fastq.gz",
        r2=tmp_dir+"/Sample_B116/B116_S1_L006_R2_001.fastq.gz",
        barcodes="workflow/resources/flexbar_demux/B116_barcodes.fasta",
        adapter_r1="workflow/resources/flexbar_demux/adapters_r1.fasta",
        adapter_r2="workflow/resources/flexbar_demux/adapters_r2.fasta"
    output:
        r1=expand(tmp_dir+"/Sample_B116/{run}_1.fastq", run = samples[samples['experiment']=='B116']['run'].values),
        r2=expand(tmp_dir+"/Sample_B116/{run}_2.fastq", run = samples[samples['experiment']=='B116']['run'].values)
    params:
        target=tmp_dir+"/Sample_B116/"
    log: logs_dir+"/B116/B116_demux.log"
    conda:
        '../envs/sm_flexbar.yaml'
    threads: 24
    shell:
        '''\
        flexbar \
        --threads {threads} \
        --reads {input.r1} \
        --reads2 {input.r2} \
        --target {params.target} \
        --output-log {log} \
        \
        --barcodes {input.barcodes} \
        --barcode-trim-end LTAIL \
        --barcode-error-rate 0 \
        --barcode-keep \
        \
        --adapters {input.adapter_r1} \
        --adapters2 {input.adapter_r2} \
        --adapter-pair-overlap ON \
        --adapter-min-poverlap 20 \
        --adapter-add-barcode \
        --adapter-trim-end RIGHT \
        --adapter-error-rate 0.1 \
        --adapter-min-overlap 1 \
        \
        --min-read-length 21 \
        --single-reads-paired &&
        rename "_barcode_" "" {params.target}_barcode_*.fastq
        '''

rule umi_extract_B116:
    wildcard_constraints:
        run='|'.join(samples.query(f" source != 'ENA' and experiment == 'B116'")['run'].tolist()+["dummy_constraint"])
    input:
        r1=tmp_dir+"/Sample_B116/{run}_1.fastq",
        r2=tmp_dir+"/Sample_B116/{run}_2.fastq"
    output:
        r1=fastq_dir+'/{experiment}/{run}_1.fastq.gz',
        r2=fastq_dir+'/{experiment}/{run}_2.fastq.gz'
    params:
        tmp=tmp_dir
    log:
        logs_dir+'/{experiment}/umitools/{run}.extract.log'
    threads: 1
    conda:
        '../envs/sm_umitools.yaml'
    shell:
        '''
        umi_tools extract \
        --bc-pattern='^(?P<umi_1>.{{3}})(?P<discard_1>.{{4}})(?P<umi_2>.{{2}})' \
        --extract-method=regex \
        --stdin={input.r1} \
        --stdout={output.r1} \
        --read2-in={input.r2} \
        --read2-out={output.r2} \
        --log={log} \
        --temp-dir={params.tmp}
        '''

# --- DRB data --- #

rule download_DRB:
    """\
    Download the tar archive for the lane(s) of the DRB miCLIP experiment.
    """
    params:
        outdir=tmp_dir,
        pwd="--password=PASSWORD"
    output:
        temp(tmp_dir+"/Sample_DRB_miCLIP.tar")
    shell:
        '''\
        wget \
        -O {output} \
        --no-check-certificate \
        --content-disposition \
        --continue \
        --user=bal2008@med.cornell.edu \
        {params.pwd} \
        "https://abc.med.cornell.edu/pubshare/data/2/44038/Project_EC-AG-3017/Sample_DRB_miCLIP/*.fastq.gz"
        '''

rule untar_DRB:
    """\
    Extract the tar archive for the lane(s) of the DRB miCLIP experiment. This archive has split fastq files, so we concatenate them.
    """
    input:
        tmp_dir+"/Sample_DRB_miCLIP.tar"
    output:
        r1=tmp_dir+"/Sample_DRB_miCLIP/DRB_miCLIP_R1.fastq.gz",
        r2=tmp_dir+"/Sample_DRB_miCLIP/DRB_miCLIP_R2.fastq.gz"
    params:
        outdir=tmp_dir
    shell:
        '''\
        tar -C {params.outdir} -xvf {input} && 
        cat {params.outdir}/Sample_DRB_miCLIP/*_R1_*.fastq.gz > {params.outdir}/Sample_DRB_miCLIP/DRB_miCLIP_R1.fastq.gz && 
        rm  {params.outdir}/Sample_DRB_miCLIP/*_R1_*.fastq.gz &&
        cat {params.outdir}/Sample_DRB_miCLIP/*_R2_*.fastq.gz > {params.outdir}/Sample_DRB_miCLIP/DRB_miCLIP_R2.fastq.gz &&
        rm  {params.outdir}/Sample_DRB_miCLIP/*_R2_*.fastq.gz
        '''

rule demux_DRB:
    """\
    Demultiplex the DRB fastq file using the barcode fasta file provided in workflow/resources.
    Then remove target prefix ()
    """
    wildcard_constraints:
        run='|'.join(samples.query(f" experiment == 'DRB' and layout == 'PAIRED' ")['run'].tolist()+["dummy_constraint"])
    input:
        r1=tmp_dir+"/Sample_DRB_miCLIP/DRB_miCLIP_R1.fastq.gz",
        r2=tmp_dir+"/Sample_DRB_miCLIP/DRB_miCLIP_R2.fastq.gz",
        barcodes="workflow/resources/flexbar_demux/DRB_barcodes.fasta",
        adapter_r1="workflow/resources/flexbar_demux/adapters_r1.fasta",
        adapter_r2="workflow/resources/flexbar_demux/adapters_r2.fasta"
    output:
        r1=temp(expand(tmp_dir+"/Sample_DRB_miCLIP/{run}_1.fastq", run = samples[samples['experiment']=='DRB']['run'].values)),
        r2=temp(expand(tmp_dir+"/Sample_DRB_miCLIP/{run}_2.fastq", run = samples[samples['experiment']=='DRB']['run'].values))
    params:
        target=tmp_dir+"/Sample_DRB_miCLIP/"
    log:
        logs_dir+"/DRB/DRB_demux.log"
    conda:
        '../envs/sm_flexbar.yaml'
    threads: 24
    shell:
        '''\
        flexbar \
        --threads {threads} \
        --reads {input.r1} \
        --reads2 {input.r2} \
        --target {params.target} \
        --output-log {log} \
        \
        --barcodes {input.barcodes} \
        --barcode-trim-end LTAIL \
        --barcode-error-rate 0 \
        --barcode-keep \
        \
        --adapters {input.adapter_r1} \
        --adapters2 {input.adapter_r2} \
        --adapter-pair-overlap ON \
        --adapter-min-poverlap 20 \
        --adapter-add-barcode \
        --adapter-trim-end RIGHT \
        --adapter-error-rate 0.1 \
        --adapter-min-overlap 1 \
        \
        --min-read-length 21 \
        --single-reads-paired \
        &&
        rename "_barcode_" "" {params.target}_barcode_*.fastq
        '''

rule umi_extract_DRB:
    wildcard_constraints:
        run='|'.join(samples.query(f" source != 'ENA' and experiment == 'DRB'")['run'].tolist()+["dummy_constraint"])
    input:
        r1=tmp_dir+"/Sample_DRB_miCLIP/{run}_1.fastq",
        r2=tmp_dir+"/Sample_DRB_miCLIP/{run}_2.fastq"
    output:
        r1=fastq_dir+'/{experiment}/{run}_1.fastq.gz',
        r2=fastq_dir+'/{experiment}/{run}_2.fastq.gz'
    params:
        tmp=tmp_dir
    log:
        logs_dir+'/{experiment}/umitools/{run}.extract.log'
    threads: 1
    conda:
        '../envs/sm_umitools.yaml'
    shell:
        '''
        umi_tools extract \
        --bc-pattern='^(?P<umi_1>.{{3}})(?P<discard_1>.{{4}})(?P<umi_2>.{{2}})' \
        --extract-method=regex \
        --stdin={input.r1} \
        --stdout={output.r1} \
        --read2-in={input.r2} \
        --read2-out={output.r2} \
        --log={log} \
        --temp-dir={params.tmp}
        '''


# --- Nature Methods data --- #

rule download_NatMeth_PE:
    """\
    """
    wildcard_constraints:
        run='|'.join(samples.query(f" experiment == 'NatMeth' and layout == 'PAIRED' ")['run'].tolist()+["dummy_constraint"])
    params:
        url1=lambda wildcards: get_ena_path(wildcards.run)+samples.loc[wildcards.run, 'read1'],
        url2=lambda wildcards: get_ena_path(wildcards.run)+samples.loc[wildcards.run, 'read2']
    output:
        r1=tmp_dir+"/NatMeth/download/{run}_1.fastq.gz",
        r2=tmp_dir+"/NatMeth/download/{run}_2.fastq.gz"
    shell:
        '''\
        wget -O - {params.url1} > {output.r1} &&
        wget -O - {params.url2} > {output.r2}
        '''

rule demux_NatMeth_PE:
    """\
    """
    wildcard_constraints:
        run='|'.join(samples.query(f" experiment == 'NatMeth' and layout == 'PAIRED' ")['run'].tolist()+["dummy_constraint"])
    input:
        r1=tmp_dir+"/NatMeth/download/{run}_1.fastq.gz",
        r2=tmp_dir+"/NatMeth/download/{run}_2.fastq.gz",
        barcodes="workflow/resources/flexbar_demux/NatMeth_barcodes.fasta",
        adapter_r1="workflow/resources/flexbar_demux/adapters_r1.fasta",
        adapter_r2="workflow/resources/flexbar_demux/adapters_r2.fasta"
    output:
        r1=tmp_dir+"/NatMeth/demux_{run}/{run}_1.fastq",
        r2=tmp_dir+"/NatMeth/demux_{run}/{run}_2.fastq"
    params:
        target=tmp_dir+"/NatMeth/demux_{run}/"
    log: logs_dir+"/NatMeth_demux_{run}.log"
    conda:
        '../envs/sm_flexbar.yaml'
    threads: 24
    shell:
        '''\
        flexbar \
        --threads {threads} \
        --reads {input.r1} \
        --reads2 {input.r2} \
        --target {params.target} \
        --output-log {log} \
        \
        --barcodes {input.barcodes} \
        --barcode-trim-end LTAIL \
        --barcode-error-rate 0 \
        --barcode-keep \
        \
        --adapters {input.adapter_r1} \
        --adapters2 {input.adapter_r2} \
        --adapter-pair-overlap ON \
        --adapter-min-poverlap 20 \
        --adapter-add-barcode \
        --adapter-trim-end RIGHT \
        --adapter-error-rate 0.1 \
        --adapter-min-overlap 1 \
        \
        --min-read-length 21 \
        --single-reads-paired \
        &&
        rename "_barcode_" "" {params.target}_barcode_*.fastq
        '''

rule umi_extract_NatMeth_PE:
    wildcard_constraints:
        run='|'.join(samples.query(f" source != 'ENA' and experiment == 'NatMeth' and layout == 'PAIRED' ")['run'].tolist()+["dummy_constraint"])
    input:
        r1=tmp_dir+"/NatMeth/demux_{run}/{run}_1.fastq",
        r2=tmp_dir+"/NatMeth/demux_{run}/{run}_2.fastq"
    output:
        r1=fastq_dir+'/NatMeth/{run}_1.fastq.gz',
        r2=fastq_dir+'/NatMeth/{run}_2.fastq.gz'
    params:
        tmp=tmp_dir
    log:
        logs_dir+'/NatMeth/umitools/{run}.extract.log'
    threads: 1
    conda:
        '../envs/sm_umitools.yaml'
    shell:
        '''
        umi_tools extract \
        --bc-pattern='^(?P<umi_1>.{{3}})(?P<discard_1>.{{4}})(?P<umi_2>.{{2}})' \
        --extract-method=regex \
        --stdin={input.r1} \
        --stdout={output.r1} \
        --read2-in={input.r2} \
        --read2-out={output.r2} \
        --log={log} \
        --temp-dir={params.tmp}
        '''

rule download_NatMeth_SE:
    """\
    """
    wildcard_constraints:
        run='|'.join(samples.query(f" experiment == 'NatMeth' and layout == 'SINGLE' ")['run'].tolist()+["dummy_constraint"])
    params:
        url1=lambda wildcards: get_ena_path(wildcards.run)+samples.loc[wildcards.run, 'read1']
    output:
        r1=tmp_dir+"/NatMeth/download/{run}.fastq.gz"
    shell:
        '''\
        wget -O - {params.url1} > {output.r1}
        '''

rule demux_NatMeth_SE:
    """\
    """
    wildcard_constraints:
        run='|'.join(samples.query(f" experiment == 'NatMeth' and layout == 'SINGLE' ")['run'].tolist()+["dummy_constraint"])
    input:
        r1=tmp_dir+"/NatMeth/download/{run}.fastq.gz",
        barcodes="workflow/resources/flexbar_demux/NatMeth_barcodes.fasta",
        adapter_r1="workflow/resources/flexbar_demux/adapters_r1.fasta"
    output:
        r1=tmp_dir+"/NatMeth/demux_{run}/{run}.fastq"
    params:
        target=tmp_dir+"/NatMeth/demux_{run}/"
    log:
        logs_dir+"/NatMeth/NatMeth_demux_{run}.log"
    conda:
        '../envs/sm_flexbar.yaml'
    threads: 12
    shell:
        '''\
        flexbar \
        --threads {threads} \
        --reads {input.r1} \
        --target {params.target} \
        --output-log {log} \
        \
        --barcodes {input.barcodes} \
        --barcode-trim-end LTAIL \
        --barcode-error-rate 0 \
        --barcode-keep \
        \
        --adapters {input.adapter_r1} \
        --adapter-add-barcode \
        --adapter-trim-end RIGHT \
        --adapter-error-rate 0.1 \
        --adapter-min-overlap 1 \
        \
        --min-read-length 21 \
        &&
        rename "_barcode_" "" {params.target}_barcode_*.fastq 
        '''

rule umi_extract_NatMeth_SE:
    wildcard_constraints:
        run='|'.join(samples.query(f" source != 'ENA' and experiment == 'NatMeth' and layout == 'SINGLE' ")['run'].tolist()+["dummy_constraint"])
    input:
        r1=tmp_dir+"/NatMeth/demux_{run}/{run}.fastq"
    output:
        r1=fastq_dir+'/NatMeth/{run}.fastq.gz'
    params:
        tmp=tmp_dir
    log:
        logs_dir+'/NatMeth/umitools/{run}.extract.log'
    threads: 1
    conda:
        '../envs/sm_umitools.yaml'
    shell:
        '''
        umi_tools extract \
        --bc-pattern='^(?P<umi_1>.{{3}})(?P<discard_1>.{{4}})(?P<umi_2>.{{2}})' \
        --extract-method=regex \
        --stdin={input.r1} \
        --stdout={output.r1} \
        --log={log} \
        --temp-dir={params.tmp}
        '''

