#!/bin/bash
#SBATCH --job-name=logan-sanity
#SBATCH --partition=main
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH -o logs/sanity-%j.out
#SBATCH -e logs/sanity-%j.out
set -uo pipefail
module load singularity/3.8.5
H=$HOME
SB="singularity exec -B /path/to/hostfs"
echo "== genomad =="; $SB $H/.singularity_cache/quay.io-biocontainers-genomad-1.11.0--pyhdfd78af_0.img genomad --version 2>&1 | head -2
echo "== checkv =="; $SB $H/.singularity_cache/quay.io-biocontainers-checkv-1.0.3--pyhdfd78af_0.img checkv --help 2>&1 | head -2
echo "== coverm =="; $SB $H/.singularity_cache/quay.io-biocontainers-coverm-0.7.0--hcb7b614_4.img coverm --version 2>&1 | head -2
echo "== metabat2 =="; $SB $H/simgs/depot.galaxyproject.org-singularity-metabat2-2.17--hd498684_0.img metabat2 --help 2>&1 | head -2
echo "== checkm2 =="; $SB $H/simgs/community.wave.seqera.io-library-checkm2_diamond_pip_binette-d2744f87fb84aab7.img checkm2 --help 2>&1 | head -2
echo "== DBs readable =="
ls -d $H/databases/genomad-db-v1.7/genomad_db $H/databases/checkv-db-v1.5 $H/databases/CheckM2_database/uniref100.KO.1.dmnd 2>&1
echo "== inputs.tsv =="; cat $H/Projects/logan-test/inputs.tsv
echo "ALL-SANITY-DONE"