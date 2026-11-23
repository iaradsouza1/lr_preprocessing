# Read preprocessing: optional ONT orientation (pychopper), then length filter (100-15000 bp).


rule prepare_reads:
    input:
        fq=raw_fastq,
    output:
        fq=f"{OUT}/reads/{{sample}}.filtered.fastq.gz",
    params:
        do_pychopper=lambda w: is_ont_unoriented(w.sample),
        minlen=config["length_filter"]["min"],
        maxlen=config["length_filter"]["max"],
        tmp=lambda w: f"{OUT}/reads/{w.sample}.pychopper.fastq",
        # cDNA kit of the library (pychopper's default is PCS109); "" keeps the default
        kit=lambda w: f"-k {config['pychopper_kit']}" if config.get("pychopper_kit") else "",
        # e.g. "-m edlib": much faster than the default phmm on large runs
        extra=config.get("pychopper_extra", ""),
    threads: 8
    log:
        f"{OUT}/logs/preprocess/{{sample}}.log",
    conda:
        config["conda_env"]
    shell:
        r"""
        set -euo pipefail
        if [ "{params.do_pychopper}" = "True" ]; then
            pychopper -t {threads} {params.kit} {params.extra} {input.fq} {params.tmp} > {log} 2>&1
            seqkit seq -j {threads} -m {params.minlen} -M {params.maxlen} {params.tmp} \
                | gzip -c > {output.fq} 2>> {log}
            rm -f {params.tmp}
        else
            seqkit seq -j {threads} -m {params.minlen} -M {params.maxlen} {input.fq} \
                | gzip -c > {output.fq} 2> {log}
        fi
        """
