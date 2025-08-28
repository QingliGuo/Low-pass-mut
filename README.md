# Variant Calling & Filtering Pipeline for Low-Pass WGS Data

This Snakemake pipeline performs **variant calling using Octopus**, followed by **multi-step filtering with ANNOVAR** against several population databases, and finally generates mutation matrices for downstream analysis.  
It is designed for reproducibility and tested on SLURM-based HPC clusters.

---

## Repository Contents

- `Snakefile`: Main pipeline logic.
- `config.yaml`: Configuration file specifying samples, references, and tool/database paths.
- `submit_snakemake.sh`: SLURM job submission script for running the pipeline on a cluster.
- `scripts/`: Contains helper scripts for variant annotation and matrix generation.
- `sample_names.txt`: (not included) A list of sample names used in the pipeline.

---

## Dependencies

### 1. Snakemake environment
Create the Snakemake environment with mamba:
```bash
mamba env create -f envs/snakemake.yaml
```

### 2. ANNOVAR setup
This pipeline requires the following ANNOVAR databases (GRCh37/hg19):

- 1000g2015aug
- esp6500siv2_all
- exac03nontcga
- gnomad211_genome
- kaviar_20150923

**To install:**
1. Download ANNOVAR from [http://annovar.openbioinformatics.org/en/latest/user-guide/download/](http://annovar.openbioinformatics.org/en/latest/user-guide/download/)
2. Use these commands to download the databases:
```bash
perl annotate_variation.pl -downdb -webfrom annovar 1000g2015aug humandb/
perl annotate_variation.pl -downdb -webfrom annovar esp6500siv2_all humandb/
perl annotate_variation.pl -downdb -webfrom annovar exac03nontcga humandb/
perl annotate_variation.pl -downdb -webfrom annovar gnomad211_genome humandb/
perl annotate_variation.pl -downdb -webfrom annovar kaviar_20150923 humandb/
```

---

## Setup

1. **Clone the repository**
   ```bash
   git clone https://github.com/QingliGuo/Low-pass-mut.git
   cd Low-pass-mut
   ```

2. **Edit the configuration file**  
   Update `config.yaml` to specify:
   - `samples`: path to sample list (`sample_names.txt`)
   - `bam_dir`: directory containing BAM files
   - `ref`: reference genome path
   - `output_dir`: output directory
   - `annovar_db`: path to ANNOVAR database directory
   - Paths to helper scripts

3. **Prepare the sample list**  
   Create `sample_names.txt` containing sample names (one per line, without extensions).

---

## Usage

Submit the pipeline to SLURM:
```bash
sbatch submit_snakemake.sh
```

The submission script runs Snakemake with the following features:
- Uses `--use-conda --conda-frontend mamba` for reproducible environments
- Submits jobs to SLURM with defined resource requirements
- Resumes incomplete jobs and prints shell commands for debugging
- Keeps going on errors and shows failed logs

---

## Output

Results are written to the `output_dir` specified in `config.yaml`. Outputs include:
- Annotated variant files (from ANNOVAR)
- Mutation matrices
- Log files

---

## Clean up (optional)

To remove intermediate files and force reruns:
```bash
snakemake --delete-all-output
```

---

## Contact

By: **Qingli Guo**  
Email: qingli.guo@icr.ac.uk | qingliguo@outlook.com

---

## License

This work is licensed under the **Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)** license.  
See: [https://creativecommons.org/licenses/by-nc/4.0/](https://creativecommons.org/licenses/by-nc/4.0/)
