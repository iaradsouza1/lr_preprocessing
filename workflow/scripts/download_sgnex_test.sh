#!/usr/bin/env bash
# Download the matched PacBio + ONT A549 pair (full files) from the SG-NEx open
# S3 bucket (public, no credentials) for pipeline testing, using the AWS CLI
# (--no-sign-request, parallel multipart -> much faster than plain HTTPS).
#
# Requires `aws` on PATH (awscli is in the lr_pipeline conda env):
#   conda activate lr_pipeline        # or: conda run -n lr_pipeline bash <this script>
#
# Usage:
#   bash download_sgnex_test.sh [DEST_DIR]
#
#   DEST_DIR   where to write the fastqs (default: resources/sgnex_test). Created
#              if missing. `aws s3 cp` resumes/overwrites cleanly on re-run.
#
# Files (full):
#   PacBio HiFi cDNA : ~5.9 GB  -> SGNex_A549_PacBio_rep7.fastq.gz
#   ONT PCR-cDNA     : ~2.3 GB  -> SGNex_A549_ONTcDNA_rep1.fastq.gz
set -euo pipefail

DEST="${1:-resources/sgnex_test}"
BUCKET="s3://sg-nex-data"

PB_KEY="data/sequencing_data_pacbio/fastq/SGNex_A549_PacBio-SMRTcell_replicate7_run1/SGNex_A549_PacBio-SMRTcell_replicate7_run1.fastq.gz"
ONT_KEY="data/sequencing_data_ont/fastq/SGNex_A549_cDNA_replicate1_run2/SGNex_A549_cDNA_replicate1_run2.fastq.gz"

command -v aws >/dev/null 2>&1 || {
  echo "ERROR: 'aws' not found. Activate the env first:  conda activate lr_pipeline" >&2
  exit 1
}

mkdir -p "$DEST"
echo "[download_sgnex_test] destination: $DEST"

fetch() {  # s3key dest
  local key="$1" out="$2"
  if [ -s "$out" ]; then
    echo "  exists, skipping: $out"
  else
    echo "  downloading: $out"
    aws s3 cp --no-sign-request "$BUCKET/$key" "$out"
  fi
}

fetch "$PB_KEY"  "$DEST/SGNex_A549_PacBio_rep7.fastq.gz"
fetch "$ONT_KEY" "$DEST/SGNex_A549_ONTcDNA_rep1.fastq.gz"

echo "[download_sgnex_test] done:"
ls -lh "$DEST"
echo
echo "Verify integrity with:  gzip -t <file>.fastq.gz"
