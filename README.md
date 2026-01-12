# Bulk RNA-seq & ChIP-seq Local Pipeline
*README prepared by Janghyun Choi*

This repository provides a unified, fully local execution pipeline for bulk RNA-seq and ChIP-seq analysis, designed to run reproducibly on **Linux**, **macOS**, and **Windows systems via WSL2**. The pipeline does not rely on any cloud services or job schedulers and executes all steps sequentially using modular shell scripts. For high-throughput or SLURM-based batch execution, a separate HPC (High-Performance Computing)-oriented pipeline is maintained in a different repository.
Two main entry points are provided:

- `pipeline_rna.sh` for RNA-seq  
- `pipeline_chip.sh` for ChIP-seq  

Each pipeline orchestrates a series of modular scripts located in the `bin/` directory.

![Python](https://img.shields.io/badge/Python-3.10-blue?logo=python) ![micromamba](https://img.shields.io/badge/micromamba-env-green?logo=anaconda) ![conda-forge](https://img.shields.io/badge/channel-conda--forge-orange?logo=conda-forge) ![status](https://img.shields.io/badge/status-stable-green)

---

## Directory structure

```
.
├── environment.yml
├── pipeline_rna.sh
├── pipeline_chip.sh
├── bin/
│   ├── trimmomatic_rna.sh
│   ├── hisat2_rna.sh
│   ├── filterQC_rna.sh
│   ├── htseq_rna.sh
│   ├── generate_report.sh
│   ├── trimmomaticChip.sh
│   ├── bowtie2Chip.sh
│   ├── filterQCChip.sh
│   ├── convertingBEDChip.sh
│   ├── macs2Chip.sh
│   └── generate_reportChip.sh
├── refGenome/
│   ├── hisat2_index/
│   │   └── os_IRGSP/
│   │       └── genome.*
│   ├── bowtie2_index/
│   │   └── os_IRGSP/
│   │       └── genome.*
│   └── gtf/
│       └── os_IRGSP/
│           └── *.gtf
└── programs/
    ├── Trimmomatic-0.39/
    ├── qualimap_v2.3/
    ├── ucsc_tools/
    └── picard/
```

---

## Environment setup (micromamba)

This pipeline is designed to run in a micromamba-based isolated environment to guarantee reproducibility across servers, WSL2, and local machines.
### 1) Install micromamba
Download the latest static micromamba binary and place it in your personal `~/bin` directory:

```
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
| tar -xvj bin/micromamba

mkdir -p ~/bin
mv bin/micromamba ~/bin
```

Register it in your shell:

```
echo 'export PATH="$HOME/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

### 2) Initialize micromamba
Create a micromamba root and enable shell integration:
```
micromamba shell init -s bash -p ~/micromamba
source ~/.bashrc

```
This enables micromamba activate to behave like conda activate.

### 3) Create the analysis environment
Create the pipeline environment from the provided YAML file:
`micromamba create -n bulkSeq -f enviornment.yml`

Activate it:
`micromamba activate bulkSeq`

This environment contains all required tools (HISAT2, Bowtie2, Samtools, MACS2, HTSeq, FastQC, etc.).

### 4) Configure channels
Set channel priority for stable bioinformatics dependency resolution:

```
micromamba config append channels conda-forge
micromamba config append channels bioconda
micromamba config set channel_priority strict
```

You can verify:
`micromamba config list`

This ensures all packages are resolved in the following order:
    1. conda-forge
    2. bioconda
with strict priority, preventing mixed or incompatible builds.

### Final environment state
After setup:
- micromamba root: `~/micromamba`
- active environment: `bulkSeq`
- executable path: `~/bin/micromamba`
- channels: `conda-forge > bioconda (strict)`

With this environment active, both `pipeline_rna.sh` and `pipeline_chip.sh` run in a fully reproducible and portable execution context.

### Programs (third-party tools)
This pipeline relies on a small set of external tools provided under the local `programs/` directory.

- Trimmomatic (v0.39) — provided via Google Drive [[download, 130KB](https://drive.google.com/file/d/1lYFu4Z5v5yPycNKzApQJrcoSes4oukky/view?usp=sharing)], MD5 (Trimmomatic-0.39.zip) = 271ed9dca91132eee0c960e0ae487bcd
- Qualimap (v2.3) — provided via Google Drive [[download, 27.1MB](https://drive.google.com/file/d/1N9E_6ajaxS_9mV9UYmxBrNZ_4BmyptS5/view?usp=sharing)], MD5 (qualimap_v2.3.zip) = dfe659cb950b7902ee5c57aeab35eaff
- UCSC utilities — provided via Google Drive [[download, 1.39GB](https://drive.google.com/file/d/13VYoe7LlNZRODmk9NwAR1wRmJyHCIc6x/view?usp=sharing)], MD5 (ucsc_tools.tar.gz) = 687761e66aa9bd682efdd853350f2353
- Picard (v2.26.10) — bundled in this repository as `programs/picard/picard.jar`

After downloading the archives from the provided links, extract them into the `programs/` directory so that the folder names match those referenced by the pipeline scripts.
If required, make the binaries executable:

```
chmod +x programs/qualimap_v2.3/qualimap
chmod +x programs/ucsc_tools/*
```

### Optional: Dedicated JupyterLab environment
JupyterLab should **never be installed inside the analysis environment** (`bulkSeq`).
It must run in a separate, clean environment to avoid Python and library conflicts.
Create a dedicated JupyterLab environment:

`micromamba create -n jlab -c conda-forge python=3.10 jupyterlab ipykernel`

Activate it:

`micromamba activate jlab`

Launch JupyterLab:

`jupyter lab --no-browser`

### Why a separate JupyterLab environment is required
The bulkSeq environment contains many low-level bioinformatics libraries (HTSlib, pysam, numpy, samtools bindings, etc.).
Installing JupyterLab into this environment often breaks Python ABI compatibility and leads to errors such as:
- missing `_sysconfigdata_*`
- broken `ipykernel`
- mismatched `glibc`, `libstdc++`, or `openssl`

By running JupyterLab in its own environment (`jlab`) and connecting kernels from `bulkSeq`, you get:
- a stable Jupyter UI
- a reproducible bioinformatics environment
- no dependency conflicts
This is the recommended production setup.

---

## Reference genome system

The pipeline is species-agnostic and uses a **folder-based genome selection system**.
The rice reference genome (IRGSP, `os_IRGSP`) is provided separately via Google Drive and is not included in this repository.
- HISAT2 (`tar.gz`): [[download, 519.8MB](https://drive.google.com/file/d/1ymDnZXkzls3CNEPQ-YCFLNtCY_dXLlb6/view?usp=sharing)], MD5 (hisat2_rice.tar.gz) = 34af26246c04e5e5bc49d964fb16a980
- Bowtie2 (`tar.gz`): [[download, 460.2MB](https://drive.google.com/file/d/1Blxsuo9spnK00ufgiXBdF_DpO1c99R-B/view?usp=sharing)], MD5 (bowtie2_rice.tar.gz) = 4d2e9d9ce2373ce69a166e2eb50d4771
- rice(MSU)_gtf (`gtf`): [[download, 67MB](https://drive.google.com/file/d/1iamoHElymZHQ4qzM6QbNqCXwtWgZv5bA/view?usp=sharing)], MD5 (all.gtf) = 0bd3e7f6ea7097fd8e9933c9aba36ca4

All reference data are stored under `refGenome/` in three subdirectories:

- `hisat2_index/`
- `bowtie2_index/`
- `gtf/`

Each of these contains one or more **species folders**, for example:

```
refGenome/hisat2_index/os_IRGSP/
refGenome/bowtie2_index/os_IRGSP/
refGenome/gtf/os_IRGSP/
```

The folder name (`os_IRGSP`) is the **genome key** used by the pipeline.  
When the user inputs `os_IRGSP`, the scripts automatically resolve:

- HISAT2 index → `refGenome/hisat2_index/os_IRGSP/genome*`
- Bowtie2 index → `refGenome/bowtie2_index/os_IRGSP/genome*`
- GTF → `refGenome/gtf/os_IRGSP/*.gtf`

To add support for another species, simply create matching folders (refer to `Adding new species` section):

```
refGenome/hisat2_index/<newGenome>/
refGenome/bowtie2_index/<newGenome>/
refGenome/gtf/<newGenome>/
```

and place the appropriate index and GTF files inside.

---

## RNA-seq workflow

Order of execution:

1. Trimming and FastQC  
2. HISAT2 alignment  
3. BAM filtering, sorting, duplicate removal  
4. (Optional) HTSeq counting  
5. Report generation  

Run:
```
./pipeline_rna.sh
```

The script will ask for:

- Forward FASTQ
- Reverse FASTQ
- Genome key (e.g. `os_IRGSP`)
- Experiment name
- GTF file
- Thread count
- Whether to run HTSeq

Outputs are written to `<experiment_name>/`.

---

## ChIP-seq workflow

Order of execution:

1. Trimming and FastQC  
2. Bowtie2 alignment  
3. BAM filtering and duplicate removal  
4. BAM → BED conversion  
5. (Optional) MACS2 peak calling  
6. Report generation  

Run:
```
./pipeline_chip.sh
```

The script will ask for:

- Forward FASTQ
- Reverse FASTQ
- Genome key (e.g. `os_IRGSP`)
- Experiment name
- Assay type (TF or histone)
- Thread count
- Whether to run MACS2

---

## Adding new species

To support a new organism:

1. Build HISAT2 index and place in:
   ```
   refGenome/hisat2_index/<genome_key>/
   ```

2. Build Bowtie2 index and place in:
   ```
   refGenome/bowtie2_index/<genome_key>/
   ```

3. Place GTF in:
   ```
   refGenome/gtf/<genome_key>/
   ```

4. Use `<genome_key>` when prompted by the pipeline.

No code modification is required.

---

## Contact

**Main Developer and Code Architect:**  
- Janghyun Choi, Ph.D. — implemented the main pipeline modules and designed major functional components.

**Project Supervisor:**
- Janghyun Choi, Ph.D. — provided conceptual guidance and overall project supervision.

**Contributors:**  
- Hyeokchan Kwon — contributed to code implementation, module development, and pipeline testing.  
- Soobin Hwang – supported data processing, logic refinement, and validation.

Pipeline maintained by Janghyun Choi  
jchoi@inha.ac.kr
