# SQANTI3 QC: classify IsoQuant transcript models vs GENCODE v41 + CAGE + polyA.
# NOTE: SQANTI3 CLI differs slightly across versions — verify flags against the
#       cloned resources/SQANTI3 and adjust the output filenames if needed.
#       The command below is the SQANTI3 5.3-5.5 syntax (tested with v5.5.4, config sqanti3.repo_ref).


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
        # SQANTI3's env pip-installs pyarrow, which loads the first libstdc++ it finds (its own
        # folder, then the system's); the system one is too old for pandas, so put the env first
        export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${{LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}}"
        python {params.script} --isoforms {input.gtf} --refGTF {input.ref_gtf} --refFasta {input.fa} \
            --CAGE_peak {params.cage} --polyA_motif_list {params.polya} \
            -t {threads} --report skip --skipORF \
            -d {params.outdir} -o {params.prefix} > {log} 2>&1
        """
