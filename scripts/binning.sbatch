#!/bin/bash
#SBATCH --job-name=logan-bin
#SBATCH --partition=main
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --array=1-6
#SBATCH -o logs/bin-%A_%a.out
#SBATCH -e logs/bin-%A_%a.out
# Head-to-head MAG recovery: minimap2 short-read map -> jgi depth -> MetaBAT2 ->
# CheckM2, on Logan vs our contigs. Same reads, same params; contig source is the
# only variable. (minimap2 -ax sr replicates coverm's internal mapper; both
# minimap2 and samtools live in the coverm container.)
set -uo pipefail
module load singularity/3.8.5
BASE=$HOME/Projects/logan-test
F=$HOME/Projects/samples_middle_fill/fastq
COVERM=$HOME/.singularity_cache/quay.io-biocontainers-coverm-0.7.0--hcb7b614_4.img
METABAT=$HOME/simgs/depot.galaxyproject.org-singularity-metabat2-2.17--hd498684_0.img
CHECKM2=$HOME/simgs/community.wave.seqera.io-library-checkm2_diamond_pip_binette-d2744f87fb84aab7.img
CHECKM2DB=$HOME/databases/CheckM2_database/uniref100.KO.1.dmnd
CPUS=$SLURM_CPUS_PER_TASK
SB="singularity exec -B /path/to/hostfs"

read SRC ACC SAMPLE < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $BASE/inputs.tsv)
CONTIGS=$BASE/contigs/$SRC/$ACC.fa
OUT=$BASE/binning/${SRC}_${ACC}; rm -rf $OUT; mkdir -p $OUT/bins
R1=$F/${ACC}_1.fastq.gz; R2=$F/${ACC}_2.fastq.gz
BAM=$OUT/map.bam

echo "[$(date)] map minimap2 src=$SRC acc=$ACC contigs=$(grep -c '^>' $CONTIGS)"
$SB $COVERM bash -c "set -o pipefail; minimap2 -ax sr -t $CPUS $CONTIGS $R1 $R2 2>$OUT/minimap2.log \
    | samtools sort -@ $CPUS -m 1G -o $BAM - && samtools index $BAM"
echo "[$(date)] map rc=$? bam=$(du -h $BAM 2>/dev/null | cut -f1)"

echo "[$(date)] jgi depth"
$SB $METABAT jgi_summarize_bam_contig_depths --outputDepth $OUT/depth.tsv $BAM
echo "[$(date)] depth lines=$(wc -l < $OUT/depth.tsv 2>/dev/null)"

echo "[$(date)] metabat2"
$SB $METABAT metabat2 -i $CONTIGS -a $OUT/depth.tsv -o $OUT/bins/bin -t $CPUS -v
NB=$(ls $OUT/bins/*.fa 2>/dev/null | wc -l)
echo "[$(date)] metabat2 produced $NB bins"

if [ "$NB" -gt 0 ]; then
  echo "[$(date)] checkm2 predict"
  $SB $CHECKM2 checkm2 predict --threads $CPUS --input $OUT/bins -x fa \
      --output-directory $OUT/checkm2 --database_path $CHECKM2DB --force
  echo "[$(date)] checkm2 rc=$?"
fi
rm -f $BAM $BAM.bai   # drop the large BAM once depth is computed
echo "[$(date)] DONE src=$SRC acc=$ACC bins=$NB"