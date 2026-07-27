#!/bin/bash
#SBATCH --job-name=logan-gtdb
#SBATCH --partition=main
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --time=12:00:00
#SBATCH --array=1-6
#SBATCH -o logs/gtdb-%A_%a.out
#SBATCH -e logs/gtdb-%A_%a.out
# Taxonomic concordance: GTDB-Tk r226 classify_wf on the MAGs from each binning
# task. Same reference, same params; do Logan and our contigs recover the same
# genomes with the same classifications?
set -uo pipefail
module load singularity/3.8.5
BASE=$HOME/Projects/logan-test
GTDBTK=$HOME/simgs/ecogenomic-gtdbtk-2.4.1.img
GTDBDB=$HOME/databases/release226
CPUS=$SLURM_CPUS_PER_TASK
SB="singularity exec -B /path/to/hostfs --env GTDBTK_DATA_PATH=$GTDBDB"

read SRC ACC SAMPLE < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $BASE/inputs.tsv)
BINS=$BASE/binning/${SRC}_${ACC}/bins
OUT=$BASE/gtdbtk/${SRC}_${ACC}; mkdir -p $OUT/tmp
NB=$(ls $BINS/*.fa 2>/dev/null | wc -l)
echo "[$(date)] gtdbtk src=$SRC acc=$ACC sample=$SAMPLE bins=$NB"
[ "$NB" -eq 0 ] && { echo "no bins produced, nothing to classify"; exit 0; }

$SB $GTDBTK gtdbtk classify_wf --genome_dir $BINS -x fa --out_dir $OUT \
    --cpus $CPUS --pplacer_cpus 1 --skip_ani_screen --tmpdir $OUT/tmp
echo "[$(date)] gtdbtk rc=$?"
rm -rf $OUT/tmp
echo "[$(date)] DONE src=$SRC acc=$ACC"