# Variant Calling & Filtering Pipeline for Low-Pass WGS Data

This Snakemake pipeline performs **variant calling using [Octopus](https://luntergroup.github.io/octopus/)** from BAM files, followed by **multi-step filtering with [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/)** against several population databases, and finally generates mutation matrices for downstream analysis, such as our MSI prediction tool **[MILO](https://github.com/QingliGuo/MILO)**. This pipeline is designed for low-pass DNA sequencing data (GRCh37/hg19) and tested on SLURM-based HPC clusters.

---

## Folder structure

Clone or download this git repository and ensure your folder structure is as follows:

```bash
.
├── Snakefile                 # The main workflow definition
├── config.yaml               # Configuration file for paths and parameters
├── submit_snakemake.sh       # SLURM submission script
├── data/                     # A directory for input data, including BAM files and the sample list.
│   ├── sample1.bam           # Toy bam files (md5sum of ref GRCh37-lite.fa: 89cdf2f85ef8559a3834a82d6bdb1b2d) 
│   ├── sample1.bam.bai
│   ├── sample2.bam
│   ├── sample2.bam.bai
│   └── sample_names.txt      # List of sample names, one per line
├── envs/                     # Conda environment for matrix generation
│   ├── matrix.yaml           # Conda environment for matrix generation
│   ├── milo.yaml             # Conda environment for MILO
│   └── octopus.yaml          # Conda environment for Octopus
└── logs/                     # Directory for log files
```
---

## Dependencies

### 1. Snakemake environment
Create the Snakemake environment with mamba:
```bash
mamba env create -f envs/snakemake.yaml --channel-priority flexible
```

### 2. ANNOVAR setup
This pipeline requires the following ANNOVAR databases (GRCh37/hg19):

- 1000g2015aug
- esp6500siv2_all
- exac03nontcga
- gnomad211_genome
- kaviar_20150923

**To install:**
1. Download ANNOVAR from [its official website](http://annovar.openbioinformatics.org/en/latest/user-guide/download/)
2. Use the `annotate_variation.pl` script to download the databases into your `humandb` directory:
```bash
perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar 1000g2015aug humandb/
perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar esp6500siv2_all humandb/
perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar exac03nontcga humandb/
perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar gnomad211_genome humandb/
perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar kaviar_20150923 humandb/
```

### 3. Matrix environment
To run the [SigProfilerMatrixGenerator](https://github.com/AlexandrovLab/SigProfilerMatrixGenerator), you need to install the reference genome in Python within the matrix environment.

```bash
from SigProfilerMatrixGenerator import install as genInstall
genInstall.install('GRCh37', rsync=False, bash=True)
```
---

## Usage

1. **Edit the configuration file**  
   Update `config.yaml` to specify correct paths for your system:
   - `samples`: path to sample list file, containing sample names (one per line, without extensions).
   - `bam_dir`: directory containing your input BAM files.
   - `ref`: Path to the reference genome FASTA file.
   - `output_dir`: output directory
   - `annovar_db`: path to ANNOVAR database directory
   - Paths to helper scripts

2. **Update job submission script**  
   Update `submit_snakemake.sh` to specify correct settings for your system.
   Submit the pipeline to SLURM:
```bash
cd <path-to-your-workdir>
sbatch submit_snakemake.sh
```
---

## Output

Results are written to the `output_dir` specified in `config.yaml`. Outputs include:
- Variant Calling files (from [Octopus](https://luntergroup.github.io/octopus/))
- Filetered variant files (from [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/))
- Mutation matrices (from [SigProfilerMatrixGenerator](https://github.com/AlexandrovLab/SigProfilerMatrixGenerator))
 
---

## License

This work is licensed under the **Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)** license.  
See: [https://creativecommons.org/licenses/by-nc/4.0/](https://creativecommons.org/licenses/by-nc/4.0/)

---

## Citation

If you use this pipeline in your research, please cite:

Guo, Q. *et al.*  *Long deletion signatures in repetitive genomic regions track somatic evolution and enable sensitive detection of microsatellite instability.* *bioRxiv* 2024.10.03.616572 (2024) doi [10.1101/2024.10.03.616572](https://doi.org/10.1101/2024.10.03.616572)

[Link to Preprint](https://doi.org/10.1101/2024.10.03.616572)

---

## Contact

By: **Qingli Guo**  
Email: qingliguo@outlook.com

