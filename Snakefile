configfile:"config/samples.yaml"
configfile:"config/config.yaml"

big_data_dir = config["big_data_dir"]
tomtom_db_dir = config["tomtom_db_dir"]

with open("config/samples.txt", "r") as f:
    samples = [line.strip() for line in f]

final_figures = [
    "figures/2A.png",
    "figures/2B_1.pdf",
    "figures/2B_2.pdf",
    "figures/2C_1.pdf",
    "figures/2C_2.pdf",
    "figures/2D.png",
    "figures/2E.pdf",
    "figures/2F.png",
    "figures/2G.pdf",
    "figures/meme.html",
    "figures/tomtom.html",
    "data/fimo_out",
    "figures/3DE",
    "figures/3FG",
    "figures/4A.pdf"
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

rule get_diff_genes:
    input:
        "data/all_results.rds"
    output:
        "data/sym_up_genes.txt",
        "data/sym_down_genes.txt"
    shell:
        '''
        Rscript scripts/intermediate/get_diff_genes.R
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


rule get_motif_counts:
    input:
        "data/fimo_out"
    output:
        "data/motif_counts.rds"
    shell:
        '''
        Rscript scripts/intermediate/get_motif_counts.R
        '''


rule meme:
    input:
        "data/fixed_top_promoters.fa"
    output:
        directory("data/meme_out"),
        "figures/meme.html"
    shell:
        '''
        meme \
        -nmotifs 4 \
        -dna \
        -maxw 25 \
        -oc data/meme_out \
        data/fixed_top_promoters.fa
        cp data/meme_out/meme.html figures
        '''


rule tomtom:
    input:
        "data/meme_out"
    output:
        directory("data/tomtom_out"),
        "figures/tomtom.html"
    shell:
        '''
        tomtom \
        -oc data/tomtom_out \
        data/meme_out/meme.txt \
        /home/liamt/meme/db/motif_databases/JASPAR/JASPAR_CORE_2016.meme
        cp data/tomtom_out/tomtom.html figures
        '''


rule fimo:
    input:
        "data/meme_out",
        "data/fixed_all_promoters.fa"
    output:
        directory("data/fimo_out"),
    shell:
        '''
        fimo \
        -oc data/fimo_out \
        --no-pgc \
        data/meme_out/meme.html \
        data/fixed_all_promoters.fa
        '''


rule get_top100:
    input:
        "data/dds.rds"
    output:
        "data/top100_genes.txt"
    shell:
        '''
        Rscript scripts/intermediate/get_top100.R
        '''


rule build_bed:
    input:
        "genome/genomic.gtf"
    output:
        "data/genes.bed"
    shell:
        r'''
        awk -F'\t' '
        $3=="gene" {{
            gene = $9
            sub(/.*gene_id "/, "", gene)
            sub(/".*/, "", gene)
            print $1 "\t" ($4-1) "\t" $5 "\t" gene "\t.\t" $7
        }}
        ' genome/genomic.gtf > data/genes.bed
        '''


rule get_top_promoters:
    input:
        "data/genes.bed",
        "genome/GCF_001417965.1_Aiptasia_genome_1.1_genomic.fna",
        "data/top100_genes.txt"
    output:
        temp("data/top100_promoters.fa")
    shell:
        '''
        ./scripts/intermediate/get_top_promoters
        '''


rule get_all_promoters:
    input:
        "data/genes.bed",
        "genome/GCF_001417965.1_Aiptasia_genome_1.1_genomic.fna"
    output:
        temp("data/all_promoters.fa")
    shell:
        '''
        ./scripts/intermediate/get_all_promoters
        '''


rule fix_top_promoters:
    input:
        "data/top100_promoters.fa"
    output:
        "data/fixed_top_promoters.fa"
    shell:
        '''
        python scripts/intermediate/fix_promoters.py data/top100_promoters.fa data/fixed_top_promoters.fa
        '''


rule fix_all_promoters:
    input:
        "data/all_promoters.fa"
    output:
        "data/fixed_all_promoters.fa"
    shell:
        '''
        python scripts/intermediate/fix_promoters.py data/all_promoters.fa data/fixed_all_promoters.fa
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
        "data/gene_info.rds",
        "data/cluster1.txt",
        "data/cluster2.txt"
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
        "genome/go_annotation.gaf",
        "data/cluster1.txt",
        "data/cluster2.txt"
    output:
        "figures/2C_1.pdf",
        "figures/2C_2.pdf",
        "data/gene2go.map"
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
        "data/gene_info.rds",
        "figures/2A.png"
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
        "data/all_results.rds",
        "figures/2A.png",
        "figures/2D.png"
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


rule fig_3DE:
    input:
        "data/motif_counts.rds"
    output:
        directory("figures/3DE")
    shell:
        '''
        Rscript scripts/final/3D,E.R
        '''


rule fig_3FG:
    input:
        "data/motif_counts.rds"
    output:
        directory("figures/3FG")
    shell:
        '''
        Rscript scripts/final/3F,G.R
        '''

rule fig_4A:
    input:
        "data/all_results.rds",
        "data/sym_up_genes.txt",
        "data/sym_down_genes.txt"
    output:
        "figures/4A.pdf"
    shell:
        '''
        Rscript scripts/final/4A.R
        '''
