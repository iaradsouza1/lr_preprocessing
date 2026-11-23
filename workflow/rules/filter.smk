# SQANTI3 rules filter driven by the paper's filtering.json (vendored).
# Produces the artifact-filtered annotation that IsoQuant re-quantifies against.
# NOTE: verify sqanti3_filter.py flags / output name against your SQANTI3 version.


rule sqanti3_filter:
    input:
        cloned=f"{SQ['repo_dir']}/.cloned",
        classification=f"{OUT}/sqanti/{{platform}}/qc/{{platform}}_classification.txt",
        corrected_gtf=f"{OUT}/sqanti/{{platform}}/qc/{{platform}}_corrected.gtf",
        json=SQ["filter_json"],
    output:
        gtf=f"{OUT}/sqanti/{{platform}}/filtered.gtf",
    params:
        script=lambda w: f"{SQ['repo_dir']}/sqanti3_filter.py",
        outdir=lambda w: f"{OUT}/sqanti/{w.platform}/filter",
        prefix=lambda w: w.platform,
    log:
        f"{OUT}/logs/sqanti/filter_{{platform}}.log",
    conda:
        config["conda_env_sqanti"]
    shell:
        r"""
        set -euo pipefail
        # SQANTI3's env pip-installs pyarrow, which loads the first libstdc++ it finds (its own
        # folder, then the system's); the system one is too old for pandas, so put the env first
        export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${{LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}}"
        # filtering.json must not test predicted_NMD: QC runs with --skipORF, so the column is
        # empty and SQANTI3 (5.5.4) fails any rule on an empty value, even one that expects "NA";
        # such a rule dropped every ISM/NIC/NNIC transcript
        python {params.script} rules \
            --sqanti_class {input.classification} \
            --filter_gtf {input.corrected_gtf} \
            --json_filter {input.json} --skip_report \
            -d {params.outdir} -o {params.prefix} > {log} 2>&1
        # normalize SQANTI3's filtered-GTF output name to a canonical path
        cp {params.outdir}/{params.prefix}.filtered.gtf {output.gtf}
        """
