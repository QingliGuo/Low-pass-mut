## ----- Qingli Guo --------
## ----- 4 Sep 2025 --------
## ----- TEST Workflow: mut calling -> filtering -> sigature generation ------
## ----- Low pass ------

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
    input:
        bam = lambda wc: os.path.join(config["bam_dir"], f"{wc.sample}.bam"),
        ref = config["ref"]
    output:
        vcf = os.path.join(output_dir, "{sample}/Variants_Octopus_v0.5.2/{sample}.octopus.vcf")
    params:
        outdir = lambda wc: os.path.join(output_dir, wc.sample, "Variants_Octopus_v0.5.2")
    log:
        "logs/{sample}_mutation_calling.log"
    conda:
        "envs/octopus.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        [ -f {output.vcf}.tmp ] && rm -f {output.vcf}.tmp
        octopus -R {input.ref} -I {input.bam} --output {output.vcf}.tmp > {log} 2>&1 && \\
        mv {output.vcf}.tmp {output.vcf}
        """

rule filter_mutation:
    threads: 4
    input:
        vcf = os.path.join(output_dir, "{sample}/Variants_Octopus_v0.5.2/{sample}.octopus.vcf")
    output:
        filtered_vcf = os.path.join(output_dir, "{sample}/Variants_Octopus_v0.5.2/{sample}.octopus.filtered.vcf")
    log:
        "logs/{sample}_filter_mutation.log"
    shell:
        """
        set -euo pipefail

        # --- Define temporary file paths as shell variables ---
        OUTDIR=$(dirname {output.filtered_vcf})
        
        AVINPUT=$OUTDIR/{wildcards.sample}.octopus.vcf.avinput
        F1=$OUTDIR/{wildcards.sample}.1.hg38_ALL.sites.2015_08_filtered
        F2=$OUTDIR/{wildcards.sample}.2.hg38_esp6500siv2_all_filtered
        F3=$OUTDIR/{wildcards.sample}.3.hg38_exac03nontcga_filtered
        F4=$OUTDIR/{wildcards.sample}.4.hg38_gnomad211_genome_filtered
        F5=$OUTDIR/{wildcards.sample}.5.hg38_kaviar_20150923_filtered

        merge_chunks() {{
            local final="$1"
            if compgen -G "${{final}}."* > /dev/null; then
                printf '%s\\n' "${{final}}."* | sort -V | xargs cat -- > "${{final}}"
                rm -f "${{final}}."*
            fi
        }}

        # --- Step 0: Convert VCF to AVINPUT ---
        {config[convert2anno]} -format vcf4 {input.vcf} -includeinfo > $AVINPUT 2> {log}

        # --- Step 1: 1000g ---
        {config[annovar]} -filter -dbtype 1000g2015aug_all -buildver hg38 \\
            -out $OUTDIR/{wildcards.sample}.1 $AVINPUT {config[annovar_db]} \\
            -maf 0.001 -thread {threads} >> {log} 2>&1
        merge_chunks $F1

        # --- Step 2: esp6500 ---
        {config[annovar]} -filter -dbtype esp6500siv2_all -buildver hg38 \\
            -out $OUTDIR/{wildcards.sample}.2 $F1 {config[annovar_db]} \\
            -score_threshold 0.001 -thread {threads} >> {log} 2>&1
        merge_chunks $F2

        # --- Step 3: exac03 ---
        {config[annovar]} -filter -dbtype exac03nontcga -buildver hg38 \\
            -out $OUTDIR/{wildcards.sample}.3 $F2 {config[annovar_db]} \\
            -score_threshold 0.001 -thread {threads} >> {log} 2>&1
        merge_chunks $F3

        # --- Step 4: gnomad ---
        {config[annovar]} -filter -dbtype gnomad211_genome -buildver hg38 \\
            -out $OUTDIR/{wildcards.sample}.4 $F3 {config[annovar_db]} \\
            -score_threshold 0.005 -thread {threads} >> {log} 2>&1
        merge_chunks $F4

        # --- Step 5: kaviar ---
        {config[annovar]} -filter -dbtype kaviar_20150923 -buildver hg38 \\
            -out $OUTDIR/{wildcards.sample}.5 $F4 {config[annovar_db]} \\
            -score_threshold 0.001 -thread {threads} >> {log} 2>&1
        merge_chunks $F5
        
        # --- Final Step: Create VCF and fix chromosome names for hg38 ---
        # Convert chromosome names from hg19 (1,2,3,X,Y,MT) to hg38 (chr1,chr2,chr3,chrX,chrY,chrM)
        cut -f 6-16 $F5 | awk 'BEGIN {{OFS="\\t"}} /^#/ {{print; next}} {{
            if ($1 == "MT") $1="chrM"
            else if ($1 !~ /^chr/) $1="chr"$1
            print
        }}' > {output.filtered_vcf}

        # --- Clean up ---
        rm -f $OUTDIR/{wildcards.sample}.*_dropped $OUTDIR/{wildcards.sample}.*.log $OUTDIR/{wildcards.sample}.*.invalid_input $OUTDIR/{wildcards.sample}.*_filtered $OUTDIR/{wildcards.sample}*avinput
        """

rule matrix_generation:
    input:
        filtered_vcf   = os.path.join(output_dir, "{sample}/Variants_Octopus_v0.5.2/{sample}.octopus.filtered.vcf"),
        matrix_generator = config["matrix_generator"]
    output:
        id83 = os.path.join(output_dir, "{sample}/Variants_Octopus_v0.5.2/matrix/output/ID/{sample}.ID83.all")
    log:
        "logs/{sample}_matrix_generation.log"
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

        python {input.matrix_generator} {wildcards.sample} $MATRIX_DIR >> {log} 2>&1

        ## clean up
        rm -rf $MATRIX_DIR/logs $MATRIX_DIR/input
        rm -f $MATRIX_DIR/{wildcards.sample}.octopus.filtered.vcf
        rm -rf $MATRIX_DIR/output/vcf_files

        test -s {output.id83}
        """
