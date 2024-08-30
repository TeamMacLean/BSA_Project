import os

# Define the ref "solanum_americanum" directory and other directories. Using solanum_americanum as my ref_directory
sol = "solanum_americanum"
reference_genome = os.path.join(sol, "ref_genome")
fastqc_dir = os.path.join(sol, "fastqc_dir")
trim_dir = os.path.join(sol, "trim_dir")
align_dir = os.path.join(sol, "align_dir")
bam_dir = os.path.join(sol, "bam_dir")
index_dir = os.path.join(sol, "index_dir")
vcf_dir = os.path.join(sol, "vcf_dir")
homosnps_dir = os.path.join(sol, "homosnps_dir")
ref_genome = os.path.join(reference_genome, "GCA_03040_.fna")

# List of all samples
samples = [
    'a2plus_L2',
    'a2plus_L2',
    'a2minus_L2',
    'a2minus_L2'   
]

# List all reads
reads = ['1', '2']

rule all:
    input:
        expand(os.path.join(fastqc_dir, '{sample}_{read}_fastqc.html'), sample=samples, read=reads),
        expand(os.path.join(trim_dir, '{sample}_{read}.trimmed.fq.gz'), sample=samples, read=reads),
        expand(os.path.join(align_dir, '{sample}_{read}.sam'), sample=samples, read=reads),
        expand(os.path.join(bam_dir, '{sample}_{read}.sorted.bam'), sample=samples, read=reads),
        expand(os.path.join(index_dir, '{sample}_{read}.sorted.bam.bai'), sample=samples, read=reads),
        expand(os.path.join(vcf_dir, '{sample}_{read}.raw.vcf'), sample=samples, read=reads),
        expand(os.path.join(homosnps_dir, '{sample}_{read}.homozygous_snps.vcf'), sample=samples, read=reads),
        expand(os.path.join(reference_genome, 'GCA_03040_H.fna.bwt'))


rule fastqc:
    input:
        os.path.join(sol, '{sample}_{read}.fq.gz')
    output:
        os.path.join(fastqc_dir, '{sample}_{read}_fastqc.html')
    shell:
        """
        fastqc {input} --outdir={fastqc_dir}
        """

rule trim:
    input:
        os.path.join(sol, '{sample}_{read}.fq.gz')
    output:
        os.path.join(trim_dir, '{sample}_{read}.trimmed.fq.gz')
    shell:
        """
        trimmomatic SE {input} {output} SLIDINGWINDOW:4:20 MINLEN:36
        """

rule index_ref:
    input:
        ref=ref_genome
    output:
        ref_index=os.path.join(reference_genome, 'GCA_03040_H.fna.bwt')
    shell:
        """
        bwa index {input.ref}
        """

rule align:
    input:
        fq=os.path.join(sol, '{sample}_{read}.fq.gz'),
        ref=ref_genome,
        index=os.path.join(reference_genome, 'GCA_03040_H.fna.bwt')
    output:
        sam=os.path.join(align_dir, '{sample}_{read}.sam')
    shell:
        """
        bwa mem {input.ref} {input.fq} > {output.sam}
        """

rule sam_bam:
    input:
        sam=os.path.join(align_dir, '{sample}_{read}.sam')
    output:
        bam=os.path.join(bam_dir, '{sample}_{read}.sorted.bam')
    shell:
        """
        samtools view -b {input.sam} | samtools sort -o {output.bam}
        """

rule index_bam:
    input:
        bam=os.path.join(bam_dir, '{sample}_{read}.sorted.bam')
    output:
        bam_index=os.path.join(index_dir, '{sample}_{read}.sorted.bam.bai')
    shell:
        """
        samtools index {input.bam} {output.bam_index}
        """

rule call_variants:
    input:
        bam=os.path.join(bam_dir, '{sample}_{read}.sorted.bam'),
        ref=ref_genome,
        bam_index=os.path.join(index_dir, '{sample}_{read}.sorted.bam.bai')
    output:
        vcf=os.path.join(vcf_dir, '{sample}_{read}.raw.vcf')
    shell:
        """
        bcftools mpileup -f {input.ref} {input.bam} | bcftools call -mv -Ov -o {output.vcf}
        """

rule homozygous_snps:
    input:
        vcf=os.path.join(vcf_dir, '{sample}_{read}.raw.vcf')
    output:
        vcf=os.path.join(homosnps_dir, '{sample}_{read}.homozygous_snps.vcf')
    shell:
        """
        "bcftools view -i 'INFO/DP > 10' {input.vcf} | bcftools filter -e 'GT="0/1"' -Ov -o {output.vcf}"
        """

        
