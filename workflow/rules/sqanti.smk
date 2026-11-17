# SQANTI3 QC: classify IsoQuant transcript models vs GENCODE v41 + CAGE + polyA.
# NOTE: SQANTI3 CLI differs slightly across versions — verify flags against the
#       cloned resources/SQANTI3 (v6.x) and adjust the output filenames if needed.


rule sqanti3_qc:
    input:
        cloned=f"{SQ['repo_dir']}/.cloned",
        gtf=f"{OUT}/isoquant/{{platform}}/models.gtf",
        ref_gtf=REF["gtf"],
        fa=REF["genome_fasta"],
        fai=REF["genome_fasta"] + ".fai",
    output:
        classification=f"{OUT}/sqanti/{{platform}}/qc/{{platform}}_classification.txt",
        corrected_gtf=f"{OUT}/sqanti/{{platform}}/qc/{{platform}}_corrected.gtf",
    params:
        script=lambda w: f"{SQ['repo_dir']}/sqanti3_qc.py",
        cage=lambda w: f"{SQ['repo_dir']}/{SQ['cage_peak']}",
        polya=lambda w: f"{SQ['repo_dir']}/{SQ['polya_motifs']}",
        outdir=lambda w: f"{OUT}/sqanti/{w.platform}/qc",
        prefix=lambda w: w.platform,
    threads: config["threads"]["sqanti"]
    log:
        f"{OUT}/logs/sqanti/qc_{{platform}}.log",
    conda:
        config["conda_env_sqanti"]
    shell:
        r"""
        set -euo pipefail
        python {params.script} {input.gtf} {input.ref_gtf} {input.fa} \
            --CAGE_peak {params.cage} --polyA_motif_list {params.polya} \
            --cpus {threads} --report skip --skipORF \
            -d {params.outdir} -o {params.prefix} > {log} 2>&1
        """
