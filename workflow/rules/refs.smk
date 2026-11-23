# Reference preparation: genome, GENCODE v41 GTF, minimap2 index, SQANTI3 clone.


rule download_genome:
    output:
        fa=REF["genome_fasta"],
    params:
        url=REF["genome_url"],
    log:
        f"{OUT}/logs/refs/download_genome.log",
    conda:
        config["conda_env"]
    shell:
        "curl -L {params.url} 2> {log} | gunzip -c > {output.fa}"


rule genome_faidx:
    input:
        fa=REF["genome_fasta"],
    output:
        fai=REF["genome_fasta"] + ".fai",
    conda:
        config["conda_env"]
    shell:
        "samtools faidx {input.fa}"


rule download_gtf:
    output:
        gtf=REF["gtf"],
    params:
        url=REF["gtf_url"],
    log:
        f"{OUT}/logs/refs/download_gtf.log",
    conda:
        config["conda_env"]
    shell:
        "curl -L {params.url} 2> {log} | gunzip -c > {output.gtf}"


rule minimap2_index:
    input:
        fa=REF["genome_fasta"],
    output:
        mmi=REF["genome_fasta"] + ".splice.mmi",
    threads: 8
    log:
        f"{OUT}/logs/refs/minimap2_index.log",
    conda:
        config["conda_env"]
    shell:
        "minimap2 -t {threads} -x splice -d {output.mmi} {input.fa} 2> {log}"


rule clone_sqanti3:
    output:
        marker=f"{SQ['repo_dir']}/.cloned",
    params:
        url=SQ["repo_url"],
        dest=SQ["repo_dir"],
        ref=lambda w: f"--branch {SQ['repo_ref']}" if SQ.get("repo_ref") else "",
    log:
        f"{OUT}/logs/refs/clone_sqanti3.log",
    shell:
        "( [ -d {params.dest}/.git ] || git clone --depth 1 {params.ref} {params.url} {params.dest} ) > {log} 2>&1 && touch {output.marker}"
