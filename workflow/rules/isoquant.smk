# IsoQuant: (1) isoform discovery + quantification per platform group,
#           (2) re-quantification against the SQANTI3-filtered annotation.
#
# IsoQuant takes ONE --data_type, so each platform is a separate invocation.
# Recent IsoQuant (>= 3.7) treats every file given with --bam as a separate sample (separate
# models and outputs). To keep ONE joint discovery over all samples, the BAMs go in as a single
# experiment through a YAML file, with `--read_group file_name` giving one count column per BAM
# (column names = sample labels).
# config `isoquant_reference_flags` (e.g. --complete_genedb for GENCODE) is added to discovery only.
# NOTE: verify IsoQuant's internal output filenames against your installed version
#       (the copy step below normalizes them to canonical paths).


localrules:
    isoquant_yaml,
    isoquant_requant,


rule isoquant_yaml:
    input:
        bams=platform_bams,
    output:
        yaml=f"{OUT}/isoquant/{{platform}}.input.yaml",
    params:
        labels=platform_bam_labels,
    run:
        files = ",\n      ".join(f'"{b}"' for b in input.bams)
        labels = ",\n      ".join(f'"{l}"' for l in params.labels)
        with open(output.yaml, "w") as fh:
            fh.write(
                '[\n  data format: "bam",\n  {\n    name: "EXPERIMENT",\n'
                f"    long read files: [\n      {files}\n    ],\n"
                f"    labels: [\n      {labels}\n    ]\n  }}\n]\n"
            )


# IsoQuant reads its gene database (SQLite) from many worker processes; on a parallel filesystem
# (BeeGFS in our runs) that is extremely slow (chr22 of a 100k-read test: 6 min with a node-local DB, still
# running after 60 min on BeeGFS). So the DB is converted on / copied to node-local disk
# (config `isoquant_local_dir`, default /tmp) and a copy is cached next to the outputs.
LOCAL_DB = r"""
        LOCAL={params.localdir}/isoquant_${{SLURM_JOB_ID:-$$}}_$RANDOM
        mkdir -p $LOCAL
        trap 'rm -rf "$LOCAL"' EXIT
        DBNAME=$(basename {input.gtf} .gtf).db
        CACHE={params.dbcache}/$DBNAME
        # a cached DB is stale once its GTF is regenerated (e.g. a re-run SQANTI3 filter)
        if [ -s "$CACHE" ] && [ "$CACHE" -nt {input.gtf} ]; then cp "$CACHE" $LOCAL/; GENEDB=$LOCAL/$DBNAME; else GENEDB={input.gtf}; fi
"""
SAVE_DB = r"""
        if [ -s "$LOCAL/$DBNAME" ] && ! {{ [ -s "$CACHE" ] && [ "$CACHE" -nt {input.gtf} ]; }}; then mkdir -p {params.dbcache}; cp "$LOCAL/$DBNAME" "$CACHE.$$" && mv -f "$CACHE.$$" "$CACHE"; fi
"""


# copies the first file matching one of the globs (space-separated; IsoQuant output names vary
# between versions, e.g. 3.13 writes *_grouped_file_name_counts.tsv)
# (ls exits non-zero when a fallback glob does not match; `|| true` keeps pipefail from aborting)
COPY_FIRST = r"""copy_first() {{ f=$(ls -d $1 2>/dev/null | head -1 || true); [ -n "$f" ] || {{ echo "no file matching $1" >&2; exit 1; }}; cp "$f" "$2"; }}"""


rule isoquant_discovery:
    input:
        yaml=f"{OUT}/isoquant/{{platform}}.input.yaml",
        bams=platform_bams,
        bais=lambda w: [b + ".bai" for b in platform_bams(w)],
        fa=REF["genome_fasta"],
        fai=REF["genome_fasta"] + ".fai",
        gtf=REF["gtf"],
    output:
        gtf=f"{OUT}/isoquant/{{platform}}/models.gtf",
        counts=f"{OUT}/isoquant/{{platform}}/discovery_grouped_counts.tsv",
        # per-sample counts on the reference (GENCODE) transcripts from the same run
        ref_counts=f"{OUT}/isoquant/{{platform}}/reference_transcript_grouped_counts.tsv",
    params:
        dtype=lambda w: isoquant_dtype(w.platform),
        extra=lambda w: platform_extra(w.platform, "isoquant_extra"),
        genedb=config.get("isoquant_reference_flags", ""),
        outdir=lambda w: f"{OUT}/isoquant/{w.platform}",
        localdir=config.get("isoquant_local_dir", "/tmp"),
        dbcache=f"{OUT}/isoquant/genedb",
    threads: config["threads"]["isoquant"]
    log:
        f"{OUT}/logs/isoquant/discovery_{{platform}}.log",
    conda:
        config["conda_env"]
    shell:
        COPY_FIRST
        + r"""
        set -euo pipefail
        """
        + LOCAL_DB
        + r"""
        IQ=$(command -v isoquant.py || command -v isoquant)
        sed 's/EXPERIMENT/isoquant/' {input.yaml} > {params.outdir}.discovery.yaml
        $IQ --data_type {params.dtype} {params.extra} \
            --reference {input.fa} --genedb $GENEDB {params.genedb} --genedb_output $LOCAL \
            --yaml {params.outdir}.discovery.yaml --read_group file_name \
            --threads {threads} --force \
            --output {params.outdir} > {log} 2>&1
        """
        + SAVE_DB
        + r"""
        copy_first "{params.outdir}/isoquant/isoquant.transcript_models.gtf" {output.gtf}
        copy_first "{params.outdir}/isoquant/isoquant.discovered_transcript_grouped*counts.tsv {params.outdir}/isoquant/isoquant.transcript_model_grouped_counts.tsv" {output.counts}
        copy_first "{params.outdir}/isoquant/isoquant.transcript_grouped*counts.tsv" {output.ref_counts}
        """


# Re-quantification runs one IsoQuant job PER SAMPLE and merges the tables: the joint run is bound by
# its slowest chromosome (~2 days on large ONT data); a per-sample job reads only its own BAM.
# With --no_model_construction reads are assigned the same way, but before writing counts IsoQuant
# zeroes transcripts without a "confirming" read in ANY sample of the run. The per-sample jobs run
# IsoQuant through scripts/isoquant_requant_unfiltered.py, which writes counts without that step and
# lists each sample's confirmed transcripts; merge_requant.py drops transcripts confirmed in no
# sample. On the test data this reproduced the joint requant exactly (counts and TPM).
# Small-chromosome subsampling is switched off: its rate depends on coverage, so per-sample runs
# would subsample chrM differently in each sample (the joint discovery run kept 1 read in 14).
rule isoquant_requant_sample:
    input:
        bam=f"{OUT}/align/{{sample}}.sorted.bam",
        bai=f"{OUT}/align/{{sample}}.sorted.bam.bai",
        fa=REF["genome_fasta"],
        fai=REF["genome_fasta"] + ".fai",
        gtf=f"{OUT}/sqanti/{{platform}}/filtered.gtf",
    output:
        counts=f"{OUT}/quant/{{platform}}/samples/{{sample}}.transcript_counts_unfiltered.tsv",
        confirmed=f"{OUT}/quant/{{platform}}/samples/{{sample}}.confirmed_transcripts.txt",
    params:
        dtype=lambda w: isoquant_dtype(w.platform),
        extra=lambda w: platform_extra(w.platform, "isoquant_extra"),
        outdir=lambda w: f"{OUT}/quant/{w.platform}/samples/{w.sample}",
        localdir=config.get("isoquant_local_dir", "/tmp"),
        dbcache=f"{OUT}/quant/genedb",
        script="workflow/scripts/isoquant_requant_unfiltered.py",
    threads: config["threads"].get("isoquant_requant", 16)
    resources:
        # chrM is not subsampled: give a retry more memory in case its pile-up needs it
        mem_mb=lambda w, attempt: 64000 * attempt,
    log:
        f"{OUT}/logs/isoquant/requant_{{platform}}_{{sample}}.log",
    conda:
        config["conda_env"]
    shell:
        COPY_FIRST
        + r"""
        set -euo pipefail
        """
        + LOCAL_DB
        + r"""
        rm -rf {params.outdir} {params.outdir}.confirmed.d
        mkdir -p {params.outdir} $LOCAL/confirmed
        printf '[\n  data format: "bam",\n  {{\n    name: "requant",\n    long read files: [\n      "%s"\n    ],\n    labels: [\n      "%s"\n    ]\n  }}\n]\n' \
            {input.bam} {wildcards.sample} > {params.outdir}.requant.yaml
        # IsoQuant writes (and keeps flushing) per-read temp files from every worker; with all
        # samples running, ~150 workers doing that on a parallel filesystem (BeeGFS) crawled (7% of a sample in 6 h), so
        # the run writes to node-local disk (~300 B per read; on tmpfs this counts toward mem_mb).
        # --high_memory keeps each region's alignments in memory; without it IsoQuant re-opens the
        # BAM for every region (pysam multiple_iterators), re-reading the whole .bai from the shared filesystem each
        # time, which was the other bottleneck. Counts are unchanged (checked on the test data).
        IQ_CONFIRMED_DIR=$LOCAL/confirmed python {params.script} \
            --data_type {params.dtype} {params.extra} \
            --reference {input.fa} --genedb $GENEDB --genedb_output $LOCAL \
            --yaml {params.outdir}.requant.yaml --read_group file_name \
            --no_model_construction --high_memory --max_coverage_small_chr 1000000000 \
            --threads {threads} --force \
            --output $LOCAL/out > {log} 2>&1
        """
        + SAVE_DB
        + r"""
        n=$(ls $LOCAL/confirmed | wc -l)
        [ "$n" -gt 0 ] || {{ echo "no confirmed-transcript lists written (IsoQuant patch not applied?)" >&2; exit 1; }}
        cat $LOCAL/confirmed/*.txt | sort -u > {output.confirmed}
        # keep the per-read assignments (read-level plots) and the IsoQuant log
        cp $LOCAL/out/requant/requant.read_assignments.tsv.gz $LOCAL/out/isoquant.log {params.outdir}/ || true
        copy_first "$LOCAL/out/requant/requant.transcript_grouped*counts.tsv" {output.counts}
        """


def _requant_inputs(kind):
    ext = {"counts": "transcript_counts_unfiltered.tsv", "confirmed": "confirmed_transcripts.txt"}[kind]
    return lambda w: [
        f"{OUT}/quant/{w.platform}/samples/{s}.{ext}" for s in samples_for_platform(w.platform)
    ]


rule isoquant_requant:
    # FINAL counts against the cleaned (SQANTI3-filtered) annotation, one column per sample
    input:
        counts=_requant_inputs("counts"),
        confirmed=_requant_inputs("confirmed"),
    output:
        counts=f"{OUT}/quant/{{platform}}/transcript_grouped_counts.tsv",
        tpm=f"{OUT}/quant/{{platform}}/transcript_grouped_tpm.tsv",
    params:
        script="workflow/scripts/merge_requant.py",
    log:
        f"{OUT}/logs/isoquant/requant_{{platform}}_merge.log",
    conda:
        config["conda_env"]
    shell:
        "python {params.script} {output.counts} {output.tpm}"
        " --counts {input.counts} --confirmed {input.confirmed} > {log} 2>&1"
