ruleorder: fastqscreen_PE > fastqscreen_SE

rule fastqscreen_get_genomes:
    input:
    output:
        conf=genome_dir+"/FastQ_Screen_Genomes/fastq_screen.conf"
    threads: 6
    conda:
        '../envs/sm_fastqscreen.yaml'
    shell:
        '''
        fastq_screen --get_genomes --outdir {genome_dir}
        '''

rule fastqscreen_add_mycoplasma:
    input:
        conf=genome_dir+"/FastQ_Screen_Genomes/fastq_screen.conf",
    output:
        conf=genome_dir+"/FastQ_Screen_Genomes/fastq_screen_myco.conf",
        bt2=genome_dir+"/FastQ_Screen_Genomes/Mycoplasma/Mycoplasma.1.bt2"
    threads: 6
    conda:
        '../envs/sm_fastqscreen.yaml'
    shell:
        '''
        mkdir -p {genome_dir}/FastQ_Screen_Genomes/Mycoplasma &&
        wget \
        --directory-prefix="{genome_dir}/FastQ_Screen_Genomes/Mycoplasma" \
        -r \
        --include-directories="genomes/archive/old_refseq/Bacteria/Mycoplasma_*" \
        --no-parent \
        -nH \
        --cut-dir=5 \
        -A "*fna" \
        ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Bacteria/ &&
        cat {genome_dir}/FastQ_Screen_Genomes/Mycoplasma/*.fna > {genome_dir}/FastQ_Screen_Genomes/Mycoplasma/mycoplasma.fa &&
        rm {genome_dir}/FastQ_Screen_Genomes/Mycoplasma/*.fna &&
        bowtie2-build {genome_dir}/FastQ_Screen_Genomes/Mycoplasma/mycoplasma.fa {genome_dir}/FastQ_Screen_Genomes/Mycoplasma/Mycoplasma &&
        cat {input.conf} workflow/resources/fastqscreen/mycoplasma_conf.txt > {output.conf}
        '''
rule fastqscreen_SE:
    wildcard_constraints:
        run=constraints_SE
    input:
        conf=genome_dir+"/FastQ_Screen_Genomes/fastq_screen_myco.conf",
        r1=tmp_dir+'/{experiment}/flexbar/{run}.fastq'
    output:
        txt=logs_dir+"/{experiment}/fastqscreen/{run}_screen.txt"
    params:
        outdir = logs_dir+"/{experiment}/fastqscreen/"
    threads: 6
    conda:
        '../envs/sm_fastqscreen.yaml'
    shell:
        '''
        fastq_screen \
        --force \
        --threads {threads} \
        --outdir {params.outdir} \
        --conf  {input.conf} \
        {input.r1} 
        '''

rule fastqscreen_PE:
    wildcard_constraints:
        run=constraints_PE
    input:
        conf=genome_dir+"/FastQ_Screen_Genomes/fastq_screen_myco.conf",
        r1=tmp_dir+"/{experiment}/flexbar/{run}_1.fastq"
    output:
        txt=logs_dir+"/{experiment}/fastqscreen/{run}_1_screen.txt"
    params:
        outdir = logs_dir+"/{experiment}/fastqscreen/"
    threads: 6
    conda:
        '../envs/sm_fastqscreen.yaml'
    shell:
        '''
        fastq_screen \
        --force \
        --threads {threads} \
        --outdir {params.outdir} \
        --conf  {input.conf} \
        {input.r1} 
        '''