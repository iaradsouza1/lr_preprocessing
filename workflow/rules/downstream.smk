# Optional downstream analyses (config-gated). Kept lightweight / extensible.


rule app_isoform_table:
    # Gene-of-interest isoform table (default APP) from final counts + SQANTI3 class.
    input:
        counts=expand(
            f"{OUT}/quant/{{platform}}/transcript_grouped_counts.tsv", platform=PLATFORMS
        ),
        classification=expand(
            f"{OUT}/sqanti/{{platform}}/qc/{{platform}}_classification.txt",
            platform=PLATFORMS,
        ),
        gtf=REF["gtf"],
    output:
        tsv=f"{OUT}/downstream/app_isoform_table.tsv",
    params:
        gene=config["downstream"].get("gene_of_interest", "APP"),
    log:
        f"{OUT}/logs/downstream/app_isoform_table.log",
    conda:
        config["conda_env"]
    shell:
        "python workflow/scripts/app_isoform_table.py "
        "--gene {params.gene} "
        "--counts {input.counts} "
        "--classification {input.classification} "
        "--gtf {input.gtf} "
        "--out {output.tsv} 2> {log}"


# --- SUPPA2 / DESeq2 stubs (enable in config; implement when needed) --------
# SUPPA2: generateEvents on the filtered GTF -> psiPerIsoform on the TPM table.
# DESeq2: differential isoform expression across `group` using the counts matrix.
