import os

configfile: "config.yaml"

with open(config["samples"]) as f:
    SAMPLES = [line.strip() for line in f if line.strip()]

output_dir = config["output_dir"].rstrip("/")

rule all:
    input:
        expand(
            os.path.join(output_dir, "{sample}/Variants_Octopus_v0.5.2/matrix/output/ID/{sample}.ID83.all"),
            sample=SAMPLES
        )

rule call_mutation:
    threads: 4
    resources:
        mem_mb=16000
    input:
        bam = lambda wc: os.path.join(config["bam_dir"], f"{wc.sample}.bam"),
        ref = config["ref"]
    output:
        vcf = os.path.join(output_dir, "{sample}/Variants_Octopus_v0.5.2/{sample}.octopus.vcf")
    params:
        outdir = lambda wc: os.path.join(output_dir, wc.sample, "Variants_Octopus_v0.5.2")
    log:
        os.path.join(output_dir, "logs/{sample}_mutation_calling.log")
    conda:
        "envs/octopus.yaml"
    shell:
        """
        mkdir -p {params.outdir}

        # Remove .tmp file if it exists from a previous failed run
        [ -f {output.vcf}.tmp ] && rm -f {output.vcf}.tmp

        octopus -R {input.ref} -I {input.bam} --output {output.vcf}.tmp > {log} 2>&1 && \
        mv {output.vcf}.tmp {output.vcf}
        """

rule convert_vcf_to_avinput:
    threads: 4
    resources:
        mem_mb=16000
    input:
        vcf = os.path.join(output_dir, "{sample}/Variants_Octopus_v0.5.2/{sample}.octopus.vcf")
    output:
        avinput = os.path.join(output_dir, "{sample}/Variants_Octopus_v0.5.2/{sample}.octopus.vcf.avinput")
    log:
        os.path.join(output_dir, "logs/{sample}_convert2avinput.log")
    shell:
        """
        {config[convert2anno]} -format vcf4 {input.vcf} -includeinfo > {output.avinput} 2> {log}
        """

rule filter_1000g:
    threads: 4
    resources:
        mem_mb=16000
    input:
        avinput = os.path.join(output_dir, "{sample}/Variants_Octopus_v0.5.2/{sample}.octopus.vcf.avinput")
    output:
        filtered = os.path.join(output_dir, "{sample}/Variants_Octopus_v0.5.2/{sample}.1.hg19_ALL.sites.2015_08_filtered")
    log:
        os.path.join(output_dir, "logs/{sample}_filter_1000g.log")
    shell:
        """
        set -euo pipefail
        outdir=$(dirname {output.filtered})
        outprefix=$outdir/{wildcards.sample}.1

        {config[annovar]} -filter -dbtype 1000g2015aug_all -buildver hg19 \
            -out $outprefix {input.avinput} {config[annovar_db]} \
            -maf 0.001 -thread 5 > {log} 2>&1
        """

rule filter_esp6500:
    threads: 4
    resources:
        mem_mb=16000
    input:
        filtered_prev = rules.filter_1000g.output.filtered
    output:
        filtered = os.path.join(output_dir, "{sample}/Variants_Octopus_v0.5.2/{sample}.2.hg19_esp6500siv2_all_filtered")
    log:
        os.path.join(output_dir, "logs/{sample}_filter_esp6500.log")
    shell:
        """
        set -euo pipefail
        outdir=$(dirname {output.filtered})
        outprefix=$outdir/{wildcards.sample}.2

        {config[annovar]} -filter -dbtype esp6500siv2_all -buildver hg19 \
            -out $outprefix {input.filtered_prev} {config[annovar_db]} \
            -score_threshold 0.001 -thread 5 > {log} 2>&1
        """

rule filter_exac03:
    threads: 4
    resources:
        mem_mb=16000
    input:
        filtered_prev = rules.filter_esp6500.output.filtered
    output:
        filtered = os.path.join(output_dir, "{sample}/Variants_Octopus_v0.5.2/{sample}.3.hg19_exac03nontcga_filtered")
    log:
        os.path.join(output_dir, "logs/{sample}_filter_exac03.log")
    shell:
        """
        set -euo pipefail
        outdir=$(dirname {output.filtered})
        outprefix=$outdir/{wildcards.sample}.3

        {config[annovar]} -filter -dbtype exac03nontcga -buildver hg19 \
            -out $outprefix {input.filtered_prev} {config[annovar_db]} \
            -score_threshold 0.001 -thread 5 > {log} 2>&1
        """

rule filter_gnomad:
    threads: 4
    resources:
        mem_mb=16000
    input:
        filtered_prev = rules.filter_exac03.output.filtered
    output:
        filtered = os.path.join(output_dir, "{sample}/Variants_Octopus_v0.5.2/{sample}.4.hg19_gnomad211_genome_filtered")
    log:
        os.path.join(output_dir, "logs/{sample}_filter_gnomad.log")
    shell:
        """
        set -euo pipefail
        outdir=$(dirname {output.filtered})
        outprefix=$outdir/{wildcards.sample}.4

        {config[annovar]} -filter -dbtype gnomad211_genome -buildver hg19 \
            -out $outprefix {input.filtered_prev} {config[annovar_db]} \
            -score_threshold 0.005 -thread 5 > {log} 2>&1
        """

rule filter_kaviar:
    threads: 4
    resources:
        mem_mb=16000
    input:
        filtered_prev = rules.filter_gnomad.output.filtered
    output:
        filtered     = os.path.join(output_dir, "{sample}/Variants_Octopus_v0.5.2/{sample}.5.hg19_kaviar_20150923_filtered"),
        filtered_vcf = os.path.join(output_dir, "{sample}/Variants_Octopus_v0.5.2/{sample}.octopus.filtered.vcf")
    log:
        os.path.join(output_dir, "logs/{sample}_filter_kaviar.log")
    shell:
        """
        set -euo pipefail
        outdir=$(dirname {output.filtered})
        outprefix=$outdir/{wildcards.sample}.5

        # Run Annovar filter
        {config[annovar]} -filter -dbtype kaviar_20150923 -buildver hg19 \
            -out $outprefix {input.filtered_prev} {config[annovar_db]} \
            -score_threshold 0.001 -thread 5 > {log} 2>&1

        # Cut columns 6-16 to produce the filtered VCF
        cut -f 6-16 {output.filtered} > {output.filtered_vcf}
        """

rule matrix_generation:
    threads: 4
    resources:
        mem_mb=16000
    input:
        filtered_vcf   = os.path.join(output_dir, "{sample}/Variants_Octopus_v0.5.2/{sample}.octopus.filtered.vcf"),
        matrix_generator = config["matrix_generator"],
        reorder_sbs      = config["reorder_sbs"],
        channels_csv     = config["channels_csv"],
    output:
        # UPDATED: final expected output is the ID83.all file
        id83 = os.path.join(output_dir, "{sample}/Variants_Octopus_v0.5.2/matrix/output/ID/{sample}.ID83.all")
    log:
        os.path.join(output_dir, "logs/{sample}_matrix_generation.log")
    conda:
        "envs/matrix.yaml"
    shell:
        r"""
        set -euo pipefail

        VARIANT_DIR=$(dirname {input.filtered_vcf})
        MATRIX_DIR=$VARIANT_DIR/matrix

        rm -rf $MATRIX_DIR
        mkdir -p $MATRIX_DIR

        cp {input.filtered_vcf} $MATRIX_DIR/{wildcards.sample}.octopus.filtered.vcf

        # Generate matrices (this should create $MATRIX_DIR/output/ID/{wildcards.sample}.ID83.all)
        python {input.matrix_generator} {wildcards.sample} $MATRIX_DIR

        # Optional post-processing (kept, but not required for the new final output)
        rm -rf $MATRIX_DIR/logs $MATRIX_DIR/input
        rm -f $MATRIX_DIR/{wildcards.sample}.octopus.filtered.vcf

        # Ensure the declared output exists
        test -s {output.id83}
        """

