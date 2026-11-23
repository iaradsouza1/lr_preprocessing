# Optional: fetch reads instead of reading local FASTQs. The samplesheet `fastq` path of such
# samples must be f"{outdir}/raw/{sample}.fastq.gz".
#
# config `download_source`:
#   "url"     (default) samplesheet column `url` (e.g. ENA fastq_ftp over https); optional column
#             `md5` is checked after a full download (else gzip -t).
#   "sra_aws" samplesheet column `run_accession`: the SRA object from NCBI's AWS open-data bucket
#             (often much faster than ENA), checked with vdb-validate and converted with
#             fasterq-dump (env config `conda_env_download`, needs sra-tools, pigz and aria2).
#             `--qual-defline '+'` keeps the '+' lines bare: fasterq-dump repeats the read name there
#             by default, and pychopper counts records by '\n+\n' (0 records -> it loads the
#             whole file to tune its cutoff and runs out of memory).
#             config `download_connections` (default 8): parallel byte-range connections per file.
# config `subsample_reads: N` (> 0) keeps only the first N reads (quick end-to-end tests; streams,
#   so use it with "url").
# config `temp_raw_reads: true` deletes each raw file once QC and preprocessing have read it.
# Jobs use the `download` resource, so `--resources download=2` caps parallel transfers.


rule download_fastq:
    output:
        fq=maybe_temp(f"{OUT}/raw/{{sample}}.fastq.gz", "temp_raw_reads"),
    params:
        source=config.get("download_source", "url"),
        url=lambda w: samples.loc[w.sample, "url"] if "url" in samples.columns else "",
        md5=lambda w: samples.loc[w.sample, "md5"] if "md5" in samples.columns else "",
        acc=lambda w: samples.loc[w.sample, "run_accession"] if "run_accession" in samples.columns else "",
        n=config.get("subsample_reads", 0),
        sra=lambda w: f"{OUT}/raw/{w.sample}.sra",
        conn=config.get("download_connections", 8),
    threads: 8
    resources:
        download=1,
    retries: 3
    log:
        f"{OUT}/logs/download/{{sample}}.log",
    conda:
        config.get("conda_env_download", config["conda_env"])
    shell:
        r"""
        if [ "{params.source}" = "sra_aws" ]; then
            # segmented download (resumable: -c continues a partial file after a retry)
            aria2c -q -c -x {params.conn} -s {params.conn} -k 64M --file-allocation=none \
                --max-tries=20 --retry-wait=30 -d $(dirname {params.sra}) -o $(basename {params.sra}) \
                https://sra-pub-run-odp.s3.amazonaws.com/sra/{params.acc}/{params.acc} > {log} 2>&1
            vdb-validate {params.sra} >> {log} 2>&1
            fasterq-dump --stdout --split-spot --skip-technical --qual-defline '+' --threads {threads} \
                --temp ${{TMPDIR:-/tmp}} {params.sra} 2>> {log} \
              | pigz -p {threads} -c > {output.fq}.part
            gzip -t {output.fq}.part 2>> {log}
            mv {output.fq}.part {output.fq}
            rm -f {params.sra}
        elif [ {params.n} -gt 0 ]; then
            # the first N reads only; head closes the pipe early, so pipefail is off here
            set +o pipefail
            curl -sL --retry 5 {params.url} | gunzip -c 2> /dev/null | head -n $((4 * {params.n})) | gzip -c > {output.fq}
            set -o pipefail
            gzip -t {output.fq}
        else
            wget -q -c --tries=20 --timeout=120 -O {output.fq}.part {params.url} > {log} 2>&1
            if [ -n "{params.md5}" ]; then
                echo "{params.md5}  {output.fq}.part" | md5sum -c - >> {log} 2>&1
            else
                gzip -t {output.fq}.part 2>> {log}
            fi
            mv {output.fq}.part {output.fq}
        fi
        """
