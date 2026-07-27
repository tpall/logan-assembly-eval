#!/bin/bash
#SBATCH --job-name=logan-san2
#SBATCH --partition=main
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH -o logs/sanity2-%j.out
#SBATCH -e logs/sanity2-%j.out
set -uo pipefail
module load singularity/3.8.5
H=$HOME
SB="singularity exec -B /path/to/hostfs"
probe(){ for t in "${@:2}"; do echo -n "  $t -> "; $SB "$1" bash -lc "command -v $t || echo MISSING"; done; }
echo "== bbmap img =="; probe $H/Projects/extraves/images/bbmap_39.10.sif bbmap.sh pileup.sh reformat.sh samtools java
echo "== bowtie2 img =="; probe $H/Projects/saja-aastased/singularity/depot.galaxyproject.org-singularity-bowtie2-2.4.2--py38h1c8e9b9_1.img bowtie2 bowtie2-build samtools
echo "== samtools img =="; probe $H/simgs/depot.galaxyproject.org-singularity-samtools-1.22.1--h96c455f_0.img samtools minimap2
echo "== metabat2 img =="; probe $H/simgs/depot.galaxyproject.org-singularity-metabat2-2.17--hd498684_0.img metabat2 jgi_summarize_bam_contig_depths samtools
echo "== coverm img: is minimap2 actually inside? =="; probe $H/.singularity_cache/quay.io-biocontainers-coverm-0.7.0--hcb7b614_4.img minimap2 samtools coverm
echo "ALL-SANITY2-DONE"