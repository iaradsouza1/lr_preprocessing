# Read preprocessing: optional ONT orientation (pychopper), then length filter (100-15000 bp).
#
# pychopper is bound by its single main process (parsing/writing reads), not by its worker threads:
# with -t 32 on a 120M-read sample the main process sat at 100% and the workers at ~10% (~5k reads/s,
# ~6-10 h per sample). config `pychopper_chunks: N` (default 1) splits the input into N parts
# (seqkit split2, round robin) and runs N pychopper processes side by side, each with
# threads/N - 1 workers; the cutoff is autotuned per part. Count statistics are summed over the
# parts; the per-part tuning rows (Parameter, AutotuneSample) are kept with the part number.


rule prepare_reads:
    input:
        fq=raw_fastq,
    output:
        fq=maybe_temp(f"{OUT}/reads/{{sample}}.filtered.fastq.gz", "temp_filtered_reads"),
    params:
        do_pychopper=lambda w: is_ont_unoriented(w.sample),
        minlen=config["length_filter"]["min"],
        maxlen=config["length_filter"]["max"],
        chunks=int(config.get("pychopper_chunks", 1)),
        chunkdir=lambda w: f"{OUT}/reads/{w.sample}.chunks",
        # cDNA kit of the library (pychopper's default is PCS109); "" keeps the default
        kit=lambda w: f"-k {config['pychopper_kit']}" if config.get("pychopper_kit") else "",
        stats=lambda w: f"{OUT}/qc/pychopper/{w.sample}.stats.tsv",
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
            mkdir -p $(dirname {params.stats})
            # pychopper counts records by '\n+\n'; a '+' line repeating the read name (fasterq-dump's
            # default) makes it count 0 reads and load the whole file into memory
            plus=$(set +o pipefail; gzip -cdf {input.fq} 2>/dev/null | sed -n '3{{p;q}}')
            [ "$plus" = "+" ] || {{ echo "{input.fq}: FASTQ '+' lines must be bare for pychopper (found: ${{plus:0:40}})" > {log}; exit 1; }}
            rm -rf {params.chunkdir}
            mkdir -p {params.chunkdir}
            if [ {params.chunks} -gt 1 ]; then
                seqkit split2 -p {params.chunks} -j {threads} --compress-level 1 \
                    -O {params.chunkdir} {input.fq} > {log} 2>&1
            else
                ln -s {input.fq} {params.chunkdir}/{wildcards.sample}.part_001.fastq.gz
            fi
            tper=$(( {threads} / {params.chunks} - 1 )); [ $tper -ge 1 ] || tper=1
            pids=()
            for c in {params.chunkdir}/*.part_*.f*q.gz; do
                b=${{c%.f*q.gz}}
                pychopper -t $tper {params.kit} {params.extra} -S $b.stats.tsv $c $b.pychopper.fastq > $b.log 2>&1 &
                pids+=($!)
            done
            st=0; for p in "${{pids[@]}}"; do wait $p || st=1; done
            for l in {params.chunkdir}/*.part_*.log; do echo "== $(basename $l)"; tr '\r' '\n' < $l | grep -v '%|' || true; done >> {log}
            [ $st -eq 0 ] || {{ echo "a pychopper part failed" >> {log}; exit 1; }}
            awk -F'\t' -v OFS='\t' '
                FNR == 1 {{ part++; next }}
                $1 == "Parameter" || $1 == "AutotuneSample" {{ keep[++nk] = $1 OFS $2 "_part" part OFS $3; next }}
                {{ k = $1 OFS $2; if (!(k in sum)) order[++n] = k; sum[k] += $3 }}
                END {{ print "Category", "Name", "Value"; for (i = 1; i <= n; i++) print order[i], sum[order[i]];
                       for (i = 1; i <= nk; i++) print keep[i] }}' {params.chunkdir}/*.part_*.stats.tsv > {params.stats}
            cat {params.chunkdir}/*.part_*.pychopper.fastq \
                | seqkit seq -j {threads} -m {params.minlen} -M {params.maxlen} -o {output.fq} 2>> {log}
            rm -rf {params.chunkdir}
        else
            seqkit seq -j {threads} -m {params.minlen} -M {params.maxlen} -o {output.fq} {input.fq} 2> {log}
        fi
        """
