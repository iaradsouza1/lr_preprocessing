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
        python {params.script} rules \
            --sqanti_class {input.classification} \
            --gtf {input.corrected_gtf} \
            --json_filter {input.json} \
            -d {params.outdir} -o {params.prefix} > {log} 2>&1
        # normalize SQANTI3's filtered-GTF output name to a canonical path
        cp {params.outdir}/{params.prefix}.filtered.gtf {output.gtf}
        """
