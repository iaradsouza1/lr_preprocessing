# IsoQuant: (1) isoform discovery + quantification per platform group,
#           (2) re-quantification against the SQANTI3-filtered annotation.
#
# IsoQuant takes ONE --data_type, so each platform is a separate invocation.
# Reads are grouped back to samples with `--read_group file_name` (column per BAM).
# NOTE: verify IsoQuant's internal output filenames against your installed version
#       (the `cp` lines below normalize them to canonical paths).


rule isoquant_discovery:
    input:
        bams=platform_bams,
        bais=lambda w: [b + ".bai" for b in platform_bams(w)],
        fa=REF["genome_fasta"],
        fai=REF["genome_fasta"] + ".fai",
        gtf=REF["gtf"],
    output:
        gtf=f"{OUT}/isoquant/{{platform}}/models.gtf",
        counts=f"{OUT}/isoquant/{{platform}}/discovery_grouped_counts.tsv",
    params:
        dtype=lambda w: isoquant_dtype(w.platform),
        outdir=lambda w: f"{OUT}/isoquant/{w.platform}",
    threads: config["threads"]["isoquant"]
    log:
        f"{OUT}/logs/isoquant/discovery_{{platform}}.log",
    conda:
        config["conda_env"]
    shell:
        r"""
        set -euo pipefail
        isoquant.py --data_type {params.dtype} \
            --reference {input.fa} --genedb {input.gtf} \
            --bam {input.bams} --read_group file_name \
            --threads {threads} --prefix isoquant \
            --output {params.outdir} > {log} 2>&1
        cp {params.outdir}/isoquant/isoquant.transcript_models.gtf {output.gtf}
        cp {params.outdir}/isoquant/isoquant.transcript_model_grouped_counts.tsv {output.counts}
        """


rule isoquant_requant:
    # Quantify against the cleaned (SQANTI3-filtered) annotation -> FINAL counts.
    input:
        bams=platform_bams,
        bais=lambda w: [b + ".bai" for b in platform_bams(w)],
        fa=REF["genome_fasta"],
        fai=REF["genome_fasta"] + ".fai",
        gtf=f"{OUT}/sqanti/{{platform}}/filtered.gtf",
    output:
        counts=f"{OUT}/quant/{{platform}}/transcript_grouped_counts.tsv",
        tpm=f"{OUT}/quant/{{platform}}/transcript_grouped_tpm.tsv",
    params:
        dtype=lambda w: isoquant_dtype(w.platform),
        outdir=lambda w: f"{OUT}/quant/{w.platform}",
    threads: config["threads"]["isoquant"]
    log:
        f"{OUT}/logs/isoquant/requant_{{platform}}.log",
    conda:
        config["conda_env"]
    shell:
        r"""
        set -euo pipefail
        isoquant.py --data_type {params.dtype} \
            --reference {input.fa} --genedb {input.gtf} \
            --bam {input.bams} --read_group file_name \
            --no_model_construction \
            --threads {threads} --prefix requant \
            --output {params.outdir} > {log} 2>&1
        cp {params.outdir}/requant/requant.transcript_grouped_counts.tsv {output.counts}
        cp {params.outdir}/requant/requant.transcript_grouped_tpm.tsv {output.tpm}
        """
