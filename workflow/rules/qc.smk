# QC gate 1 (read-level: seqkit stats + NanoPlot) and gate 2 (alignment identity).


rule seqkit_stats:
    input:
        fq=raw_fastq,
    output:
        tsv=f"{OUT}/qc/seqkit/{{sample}}.stats.tsv",
    log:
        f"{OUT}/logs/qc/seqkit_{{sample}}.log",
    conda:
        config["conda_env"]
    shell:
        "seqkit stats -a -T {input.fq} > {output.tsv} 2> {log}"


rule nanoplot:
    input:
        fq=raw_fastq,
    output:
        html=f"{OUT}/qc/nanoplot/{{sample}}/NanoPlot-report.html",
    params:
        outdir=lambda w: f"{OUT}/qc/nanoplot/{w.sample}",
    threads: 4
    log:
        f"{OUT}/logs/qc/nanoplot_{{sample}}.log",
    conda:
        config["conda_env"]
    shell:
        "NanoPlot -t {threads} --fastq {input.fq} -o {params.outdir} 2> {log}"


rule identity_check:
    # QC gate 2: mean minimap2 de:f identity vs the platform threshold.
    input:
        bam=f"{OUT}/align/{{sample}}.sorted.bam",
        bai=f"{OUT}/align/{{sample}}.sorted.bam.bai",
    output:
        txt=f"{OUT}/qc/identity/{{sample}}.identity.txt",
    params:
        thr=lambda w: identity_threshold(platform_of(w.sample)),
        gate=config["qc"]["identity_gate"],
    log:
        f"{OUT}/logs/qc/identity_{{sample}}.log",
    conda:
        config["conda_env"]
    shell:
        "python workflow/scripts/identity_check.py "
        "--bam {input.bam} --threshold {params.thr} --gate {params.gate} "
        "--out {output.txt} 2> {log}"
