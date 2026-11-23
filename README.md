# lr_pipeline

Generic long-read (PacBio HiFi / ONT) isoform **discovery → QC → quantification**
Snakemake pipeline, reproducing the methodology of Shimada et al. 2024
([sciadv.adj5279](https://www.science.org/doi/10.1126/sciadv.adj5279)).

## Flow
```
[optional download] -> reads (HiFi | ONT) -> seqkit stats (QC gate 1)
  -> [ONT: pychopper, optionally in parallel chunks] -> length filter (100-15000)
  -> minimap2 -ax splice --MD --secondary=no -> MAPQ>=10 -> sorted BAM (kept)
  -> identity gate 2 (mean de:f) -> IsoQuant discovery+quant -> novel GTF
  -> SQANTI3 QC -> SQANTI3 filter (filtering.json) -> filtered GTF
  -> IsoQuant re-quant vs filtered GTF (one job per sample, merged) -> FINAL counts/TPM
  -> APP isoform table [+ optional SUPPA2 / DESeq2]
```
Reference: GRCh38 primary assembly + GENCODE v41. IsoQuant discovery runs once per platform, with
all of its BAMs as one experiment (IsoQuant >= 3.7 would otherwise treat each BAM as a sample).

## Setup — two conda envs
```bash
conda env create -f workflow/envs/lr_pipeline.yaml          # main env (name: lr_pipeline)
snakemake --profile profiles/local clone_sqanti3            # clone SQANTI3 v5.5.4 (sqanti3.repo_ref) -> resources/SQANTI3
conda env create -f resources/SQANTI3/SQANTI3.conda_env.yml # SQANTI3's own env (name: sqanti3)
```
Rules reference their env **by name**: most use `config["conda_env"]` (`lr_pipeline`); the two
SQANTI3 rules use `config["conda_env_sqanti"]` (`sqanti3`). With `--use-conda` the pre-built
envs are reused, never re-created. SQANTI3 is not a package — it's cloned and run by script path.

## Configure
- `config/samples.tsv` — one row per sample: `sample, group, platform, oriented, fastq`. To download
  reads instead, set `fastq` to `<outdir>/raw/<sample>.fastq.gz` and add `url` (+ optional `md5`) or,
  with `download_source: sra_aws`, `run_accession`.
- `config/config.yaml` — reference URLs, QC thresholds, per-platform params, downstream toggles.

## Run
```bash
snakemake --profile profiles/local -n     # dry run
snakemake --profile profiles/local        # local
snakemake --profile profiles/slurm        # cluster (set account/partition first)
```

## Test data (SG-NEx, matched PacBio + ONT)
Matched A549 pair from the [SG-NEx](https://github.com/GoekeLab/sg-nex-data) open S3 bucket
(public, no credentials): PacBio HiFi cDNA (~5.9 GB) + ONT PCR-cDNA (~2.3 GB).
```bash
bash workflow/scripts/download_sgnex_test.sh          # -> resources/sgnex_test/
snakemake --profile profiles/local --config samplesheet=config/samples_sgnex.tsv -n
```
Samplesheet: `config/samples_sgnex.tsv`.

## Notes
- Input is **full-length reads FASTQ** (HiFi for PacBio; basecalled cDNA for ONT).
  Raw PacBio subreads need `ccs` first (out of scope).
- Some tool CLIs/output filenames (IsoQuant, SQANTI3) vary by version — the rules
  normalize outputs to canonical paths; verify the `cp` lines after the first real run.
- `filtering.json` is vendored from `mihshimada/Brain_Iso-Seq` (MIT), without its `predicted_NMD`
  conditions: QC runs with `--skipORF`, so that column is empty and SQANTI3 5.5 rejects every
  isoform tested on it.
- On shared/parallel filesystems, point `isoquant_local_dir` at node-local disk (default `/tmp`):
  IsoQuant's gene database and per-read temp files are very slow there.
