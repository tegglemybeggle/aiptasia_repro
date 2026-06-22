configfile:"config/samples.yaml"
configfile:"config/config.yaml"

big_data_dir = config["big_data_dir"]

with open("config/samples.txt", "r") as f:
    samples = [line.strip() for line in f]

final_figures = [
    "figures/2A.png",
    "figures/2B_1.pdf",
    "figures/2B_2.pdf",
    "figures/2C_1.pdf",
    "figures/2C_2.pdf",
    "figures/2D.png",
    "figures/2E.pdf"
]


rule aall:
    input:
        final_figures


rule get_archive:
    output:
        directory(big_data_dir + "sras/{run_id}")
    shell:
        "prefetch {wildcards.run_id} && mv {wildcards.run_id} {big_data_dir}sras"

rule fasterq_dump:
    input:
        big_data_dir + "sras/{run_id}"
    output:
        temp("{run_id}.fastq")
    resources:
        io=1
    shell:
        "cp -r {big_data_dir}sras/{wildcards.run_id} . && fasterq-dump {wildcards.run_id} && rm -r {wildcards.run_id}"

rule trim:
    input:
        "{run_id}.fastq"
    output:
        temp("{run_id}_trimmed.fq"),
        temp("{run_id}.fastq_trimming_report.txt")
    threads: 4
    resources:
        io=1
    shell:
        "trim_galore --cores {threads} {input}"

rule STAR:
    input:
        "{run_id}_trimmed.fq"
    output:
        protected(big_data_dir + "stars/{run_id}/{run_id}_Aligned.sortedByCoord.out.bam")
    threads: 8
    shell:
        '''
        pwd
        ls -lh
        mkdir {wildcards.run_id}
        STAR \
        --runThreadN {threads} \
        --genomeDir {config[genome_dir]} \
        --readFilesIn {input} \
        --outFileNamePrefix ./{wildcards.run_id}/{wildcards.run_id}_ \
        --outSAMtype BAM SortedByCoordinate
        mv {wildcards.run_id} {big_data_dir}stars
        '''

def srr_to_count_input(wildcards):
    run_id = config["samples"][wildcards.sample]
    return f"{big_data_dir}stars/" + run_id + "/" + run_id + "_Aligned.sortedByCoord.out.bam"

rule htseq_count:
    input:
        srr_to_count_input
    output:
        big_data_dir + "counts/{sample}.txt"
    resources:
        io=1
    shell:
        '''
        cp {input} .
        htseq-count \
        -f bam \
        -t exon \
        -i gene_id \
        $(basename {input}) {config[genome_dir]}genomic.gtf > {wildcards.sample}.txt
        chmod +w $(basename {input})
        rm $(basename {input}) 
        mv {wildcards.sample}.txt {big_data_dir}counts
        '''


rule build_dds:
    input:
        expand(big_data_dir + "counts/{sample}.txt", sample=samples)
    output:
        "data/dds.rds"
    shell:
        '''
        Rscript scripts/intermediate/build_dds.R
        '''

rule build_rlog:
    input:
        "data/dds.rds"
    output:
        "data/rlog.rds"
    shell:
        '''
        Rscript scripts/intermediate/build_rlog.R
        ''' 

rule get_early_genes:
    input:
        "data/dds.rds"
    output:
        "data/early_genes.txt"
    shell:
        '''
        Rscript scripts/intermediate/get_early_genes.R
        '''

rule get_clusters:
    input:
        "data/early_genes.txt",
        "data/dds.rds",
        "data/rlog.rds"
    output:
        "data/cluster1.txt",
        "data/cluster2.txt"
    shell:
        '''
        Rscript scripts/intermediate/get_clusters.R
        '''


rule get_annotations:
    input:
        "genome/genomic.gtf"
    output:
        "data/gene_info.rds"
    shell:
        '''
        Rscript scripts/intermediate/get_annotations.R
        '''


rule build_all_results:
    input:
        "data/dds.rds"
    output:
        "data/all_results.rds"
    shell:
        '''
        Rscript scripts/intermediate/build_all_results.R
        '''


rule meme:
    input:
        "data/fixed_top_promoters.fa"
    output:
        directory("data/meme_out")
    shell:
        '''
        meme \
        -nmotifs 4 \
        -dna \
        -maxw 25 \
        -oc /home/liamt/cleves/new/amid/data/meme_out \
        data/fixed_top_promoters.fa
        '''


rule fig_2A:
    input:
        "data/cluster1.txt",
        "data/cluster2.txt",
        "data/dds.rds",
        "data/rlog.rds"
    output:
        "figures/2A.png"
    shell:
        '''
        Rscript scripts/final/2A.R
        '''



rule fig_2B:
    input:
        "data/dds.rds",
        "data/gene_info.rds"
    output:
        "figures/2B_1.pdf",
        "figures/2B_2.pdf"
    shell:
        '''
        Rscript scripts/final/2B.R
        '''


rule fig_2C:
    input:
        "data/dds.rds",
        "genome/go_annotation.gaf"
    output:
        "figures/2C_1.pdf",
        "figures/2C_2.pdf",
        "gene_2_go.map"
    shell:
        '''
        Rscript scripts/final/2C.R
        '''


rule fig_2DE:
    input:
        "data/cluster1.txt",
        "data/cluster2.txt",
        "data/dds.rds",
        "data/rlog.rds",
        "data/gene_info.rds"
    output:
        "figures/2D.png",
        "figures/2E.pdf"
    shell:
        '''
        Rscript scripts/final/2D,E.R
        '''


rule fig_2F:
    input:
        "data/cluster1.txt",
        "data/cluster2.txt",
        "data/dds.rds",
        "data/rlog.rds",
        "data/all_results.rds"
    output:
        "figures/2F.png"
    shell:
        '''
        Rscript scripts/final/2F.R
        '''


rule fig_2G:
    input:
        "data/cluster1.txt",
        "data/cluster2.txt",
        "data/all_results.rds"
    output:
        "figures/2G.pdf"
    shell:
        '''
        Rscript scripts/final/2G.R
        '''


'''rule MEME:'''

