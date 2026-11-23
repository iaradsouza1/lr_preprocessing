import pandas as pd

# --- paths -----------------------------------------------------------------
OUT = config["outdir"].rstrip("/")
REF = config["reference"]
SQ = config["sqanti3"]

# --- samplesheet -----------------------------------------------------------
samples = (
    pd.read_csv(config["samplesheet"], sep="\t", comment="#", dtype=str)
    .dropna(subset=["sample"])
    .set_index("sample", drop=False)
)
samples.index.name = "sample_id"

# platforms present in the cohort (IsoQuant runs once per platform group)
PLATFORMS = sorted(samples["platform"].unique())


wildcard_constraints:
    sample="|".join(samples.index),
    platform="|".join(PLATFORMS) if PLATFORMS else "pacbio|ont",


# --- per-sample / per-platform helpers -------------------------------------
def maybe_temp(path, key):
    # temp() when config[key] is true: large intermediates are deleted once consumed
    return temp(path) if config.get(key, False) else path


def platform_of(sample):
    return samples.loc[sample, "platform"]


def samples_for_platform(p):
    return list(samples.index[samples["platform"] == p])


def raw_fastq(wildcards):
    return samples.loc[wildcards.sample, "fastq"]


def is_ont_unoriented(sample):
    row = samples.loc[sample]
    return row["platform"] == "ont" and str(row.get("oriented", "yes")).lower() == "no"


def mm2_preset(p):
    return config["platforms"][p]["minimap2_preset"]


def mm2_extra(p):
    return config["platforms"][p]["minimap2_extra"]


def isoquant_dtype(p):
    return config["platforms"][p]["isoquant_data_type"]


def identity_threshold(p):
    return config["qc"]["identity_thresholds"][p]


# BAMs feeding an IsoQuant platform-group run
def platform_bams(wildcards):
    return [
        f"{OUT}/align/{s}.sorted.bam" for s in samples_for_platform(wildcards.platform)
    ]


def platform_bam_labels(wildcards):
    return samples_for_platform(wildcards.platform)


# --- final targets ---------------------------------------------------------
def all_targets():
    t = []
    for s in samples.index:
        t.append(f"{OUT}/qc/seqkit/{s}.stats.tsv")
        t.append(f"{OUT}/qc/nanoplot/{s}/NanoPlot-report.html")
        t.append(f"{OUT}/qc/identity/{s}.identity.txt")
    for p in PLATFORMS:
        t.append(f"{OUT}/sqanti/{p}/filtered.gtf")
        t.append(f"{OUT}/quant/{p}/transcript_grouped_counts.tsv")
    if config["downstream"].get("app_isoform_table"):
        t.append(f"{OUT}/downstream/app_isoform_table.tsv")
    return t
