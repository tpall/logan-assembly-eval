#!/bin/bash
#SBATCH --job-name=logan-vir
#SBATCH --partition=main
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --array=1-6
#SBATCH -o logs/vir-%A_%a.out
#SBATCH -e logs/vir-%A_%a.out
# Head-to-head viral recovery: geNomad end-to-end + CheckV on Logan vs our contigs.
# Same DBs, same params; only the contig source differs.
set -uo pipefail
module load singularity/3.8.5
BASE=$HOME/Projects/logan-test
GENOMAD_IMG=$HOME/.singularity_cache/quay.io-biocontainers-genomad-1.11.0--pyhdfd78af_0.img
CHECKV_IMG=$HOME/.singularity_cache/quay.io-biocontainers-checkv-1.0.3--pyhdfd78af_0.img
GENOMAD_DB=$HOME/databases/genomad-db-latest/genomad_db  # v1.9 (v1.7 marker IDs crash geNomad 1.11)
CHECKV_DB=$HOME/databases/checkv-db-v1.5
CPUS=$SLURM_CPUS_PER_TASK
SB="singularity exec -B /path/to/hostfs"

read SRC ACC SAMPLE < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $BASE/inputs.tsv)
CONTIGS=$BASE/contigs/$SRC/$ACC.fa
OUT=$BASE/virome/${SRC}_${ACC}; rm -rf $OUT; mkdir -p $OUT
echo "[$(date)] geNomad start src=$SRC acc=$ACC sample=$SAMPLE contigs=$(grep -c '^>' $CONTIGS)"

$SB $GENOMAD_IMG genomad end-to-end --cleanup --splits 4 --threads $CPUS \
    $CONTIGS $OUT/genomad $GENOMAD_DB
echo "[$(date)] geNomad rc=$?"

VFNA=$OUT/genomad/${ACC}_summary/${ACC}_virus.fna
if [ -s "$VFNA" ]; then
  echo "[$(date)] CheckV on $(grep -c '^>' $VFNA) viral contigs"
  $SB $CHECKV_IMG checkv end_to_end $VFNA $OUT/checkv -d $CHECKV_DB -t $CPUS
  echo "[$(date)] checkv rc=$?"
else
  echo "[$(date)] ERROR: geNomad produced no virus.fna ($VFNA) — skipping CheckV"
fi
echo "[$(date)] DONE src=$SRC acc=$ACC"