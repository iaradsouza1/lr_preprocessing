# Splice-aware alignment -> MAPQ>=10 -> sorted, indexed BAM (kept whole-genome).


rule minimap2_align:
    input:
        fq=f"{OUT}/reads/{{sample}}.filtered.fastq.gz",
        fa=REF["genome_fasta"],
    output:
        bam=f"{OUT}/align/{{sample}}.sorted.bam",
        bai=f"{OUT}/align/{{sample}}.sorted.bam.bai",
    params:
        preset=lambda w: mm2_preset(platform_of(w.sample)),
        extra=lambda w: mm2_extra(platform_of(w.sample)),
    threads: config["threads"]["align"]
    log:
        f"{OUT}/logs/align/{{sample}}.log",
    conda:
        config["conda_env"]
    shell:
        r"""
        set -euo pipefail
        minimap2 -t {threads} -ax {params.preset} {params.extra} --MD --secondary=no \
            {input.fa} {input.fq} 2> {log} \
          | samtools view -@ {threads} -q 10 -b - \
          | samtools sort -@ {threads} -o {output.bam} - 2>> {log}
        samtools index {output.bam}
        """
