import glob
import os

# Define input path and output path
input_path = "/home/shared/raw_data/Ecoli_antibiotics_GSE220559"

# Get the list of sample names
files_R1 = glob.glob(os.path.join(input_path, "*_1.fastq.gz"))
samples = sorted(set(os.path.basename(f).replace("_1.fastq.gz", "") for f in files_R1))

print(samples)

genome_index_dir = "/home/mif/references/NC_000913/genome"
genome_gff = "/home/mif/references/NC_000913/NC_000913.gff3"

# Define the all rule: The final target files
rule all:
    input:
        #expand("results/trim_galore/{sample}_1_val_1.fq.gz", sample=samples),
        #expand("results/trim_galore/{sample}_2_val_2.fq.gz", sample=samples), 
        #expand("results/hisat/{sample}.bam", sample=samples), 
        "results/multiqc_report.html",
        "results/counts.txt"


# Define the trimGalore rule
rule trimGalore:
    input:
        fq1 = lambda wildcards: f"{input_path}/{wildcards.sample}_1.fastq.gz",
        fq2 = lambda wildcards: f"{input_path}/{wildcards.sample}_2.fastq.gz"
    output:
        fq1 = "results/trim_galore/{sample}_1_val_1.fq.gz",
        fq2 = "results/trim_galore/{sample}_2_val_2.fq.gz"
    conda: "envs/preprocess_rnaseq.yaml"
    threads: 4
    params:
        quality = 25,
        length = 20
    shell:
        """
        trim_galore --paired \
        -o results/trim_galore \
        -q {params.quality} \
        -j {threads} \
        --length {params.length} \
        {input.fq1} {input.fq2}
        """

rule hisat2:
    input:
        fq1 = "results/trim_galore/{sample}_1_val_1.fq.gz",
        fq2 = "results/trim_galore/{sample}_2_val_2.fq.gz"
    output:
        sam = "results/hisat/{sample}.sam",
        bam = "results/hisat/{sample}.bam",
        bai = "results/hisat/{sample}.bam.bai",
        stat = "results/hisat/{sample}_hisat.stat"
    conda: "envs/preprocess_rnaseq.yaml"
    threads: 8
    shell:
        """
        hisat2 -p {threads} -x {genome_index_dir} --no-spliced-alignment --no-unal --no-mixed --no-discordant -1 {input.fq1} -2 {input.fq2} -S {output.sam} --new-summary --summary-file {output.stat}
        samtools view -@ {threads} -bS {output.sam} | samtools sort -@ {threads} -o {output.bam} 
        samtools index {output.bam}
        """

rule featureCounts:
    input:
        bam = expand("results/hisat/{sample}.bam", sample=samples)
    output:
        counts = "results/counts.txt",
        stats = "results/counts.txt.summary"
    conda: "envs/preprocess_rnaseq.yaml"
    threads: 8
    shell:
        """
        featureCounts -a {genome_gff} -F GTF -t CDS -g Parent -o {output.counts} -s 2 -T {threads} -p  --countReadPairs -B  {input.bam}
        """

############
#  QC part #
############

rule fastqc:
    input:
        fq1 = lambda wildcards: f"{input_path}/{wildcards.sample}_1.fastq.gz",
        fq2 = lambda wildcards: f"{input_path}/{wildcards.sample}_2.fastq.gz",  
        fq3 = "results/trim_galore/{sample}_1_val_1.fq.gz",
        fq4 = "results/trim_galore/{sample}_2_val_2.fq.gz"
    output:
        html1 = "results/fastqc/{sample}_1_fastqc.html",
        zip1  = "results/fastqc/{sample}_1_fastqc.zip",
        html2 = "results/fastqc/{sample}_2_fastqc.html",
        zip2  = "results/fastqc/{sample}_2_fastqc.zip",
        html3 = "results/fastqc/{sample}_1_val_1_fastqc.html",
        zip3  = "results/fastqc/{sample}_1_val_1_fastqc.zip",
        html4 = "results/fastqc/{sample}_2_val_2_fastqc.html",
        zip4  = "results/fastqc/{sample}_2_val_2_fastqc.zip"
    conda: "envs/preprocess_rnaseq.yaml"
    threads: 4
    shell:
        """
            fastqc -t {threads} --outdir results/fastqc {input.fq1} {input.fq2} {input.fq3} {input.fq4}
        """ 

rule multiqc:
    input:
        expand("results/fastqc/{sample}_1_fastqc.zip", sample=samples),
        expand("results/fastqc/{sample}_1_val_1_fastqc.zip", sample=samples),
        expand("results/fastqc/{sample}_2_fastqc.zip", sample=samples),
        expand("results/fastqc/{sample}_2_val_2_fastqc.zip", sample=samples),
        "results/counts.txt.summary"
       
    output:
        "results/multiqc_report.html",
        directory("results/multiqc_data")
    log:
        "logs/multiqc.log",
    shell:
        """
        multiqc --outdir results/ {input}
        """
        