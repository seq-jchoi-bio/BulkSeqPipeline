# Changelog

## 2026.02 - 2026-06-09

### Added
- Added shared resource utilities in `bin/resource_utils.sh`.
- Added automatic CPU and memory safety checks for RNA-seq and ChIP-seq workflows.
- Added `~/` input path expansion for FASTQ and ChIP input-control paths.
- Added automatic `samtools sort -m` calculation based on detected system memory.

### Changed
- Updated `pipeline_rna.sh` and `pipeline_chip.sh` to report adjusted resource settings when user-requested threads exceed safe local workstation limits.
- Unified CPU handling across Trimmomatic, FastQC, HISAT2, Bowtie2, filtering/QC, and ChIP BAM-to-BED conversion scripts.
- Removed the fixed RNA Trimmomatic thread override (`TCORE=10`) so the user-requested value is handled through the shared safety policy.
- Updated modular `bin/` scripts so direct step-level execution uses the same CPU and memory policy as the main pipelines.

### Fixed
- Fixed interactive input paths beginning with `~/` not resolving to the user's home directory inside the pipeline scripts.
- Reduced the risk of `samtools sort` out-of-memory failures by setting explicit per-thread memory limits.

## 2026.01

- Initial public local RNA-seq and ChIP-seq pipeline release.
