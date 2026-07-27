#!/bin/bash
# Stage a Logan-vs-ours head-to-head: 3 paired-end samples, identical downstream.
set -euo pipefail
BASE=~/Projects/logan-test
F=~/Projects/samples_middle_fill/fastq
C=~/Projects/samples_middle_fill/results/assembly/contigs
mkdir -p $BASE/contigs/logan $BASE/contigs/ours $BASE/logs

# acc <TAB> sample (sample = our MEGAHIT contig id)
cat > $BASE/samples.tsv <<EOF
ERR2726406	CCMD15562448ST_11_0
ERR2726443	CCMD38158721ST_11_0
SRR2912786	CA24
EOF

# inputs.tsv: one row per (source, acc, sample); array jobs index into this
: > $BASE/inputs.tsv
while read A M; do
  echo -e "logan\t$A\t$M" >> $BASE/inputs.tsv
  echo -e "ours\t$A\t$M"  >> $BASE/inputs.tsv
done < $BASE/samples.tsv

echo "=== download + decompress Logan contigs (login node has S3) ==="
while read A M; do
  if [ ! -s $BASE/contigs/logan/$A.fa ]; then
    curl -sS --connect-timeout 30 -o $BASE/contigs/logan/$A.fa.zst \
      https://s3.amazonaws.com/logan-pub/c/$A/$A.contigs.fa.zst
    zstd -dq -f $BASE/contigs/logan/$A.fa.zst -o $BASE/contigs/logan/$A.fa
    rm -f $BASE/contigs/logan/$A.fa.zst
  fi
  # our contigs decompressed under the SAME acc name so paths are uniform
  [ -s $BASE/contigs/ours/$A.fa ] || zcat $C/MEGAHIT_${M}_contigs.fa.gz > $BASE/contigs/ours/$A.fa
  lg=$(grep -c '^>' $BASE/contigs/logan/$A.fa); ou=$(grep -c '^>' $BASE/contigs/ours/$A.fa)
  printf "  %-12s logan=%d contigs  ours=%d contigs  reads=%s\n" $A $lg $ou \
    "$([ -f $F/${A}_1.fastq.gz ] && echo ok || echo MISSING)"
done < $BASE/samples.tsv

echo "=== locate CheckM2 diamond DB ==="
find ~ -maxdepth 6 -name "*.dmnd" 2>/dev/null | grep -iE "checkm2|uniref100" | head
echo "=== confirm DB dirs ==="
ls -d ~/databases/genomad-db-v1.7/genomad_db ~/databases/checkv-db-v1.5 2>&1
echo "STAGED at $BASE"