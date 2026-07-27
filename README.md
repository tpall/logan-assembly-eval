# Logan pre-assembled contigs vs. de-novo assembly — a downstream evaluation

**Question:** can the [Logan](https://github.com/IndexThePlanet/Logan) database of
pre-assembled contigs (one contig set per SRA accession, hosted on public S3)
replace the per-sample de-novo assembly step for publicly available SRA
metagenomes? Logan would let us skip assembly entirely (download contigs
directly), so the question is whether the *downstream* results hold up.

**Short answer:** it depends on the downstream.
- **Coverage of the SRA is excellent** (~97% of our cohort accessions present).
- **Contigs are far more fragmented** than de-novo (MEGAHIT) assemblies — Logan
  keeps a long tail of short contigs, so its N50 is ~7x lower.
- **Genome-resolved binning suffers badly:** ~20x fewer high-quality MAGs.
- **Gene/contig-level and viral discovery are much more forgiving** (see viral
  results below), because those live on the longer contigs Logan does retain.

This repo is a self-contained record of that test: scripts, per-sample result
tables, and a reproducible recipe.

---

## What is Logan

Logan (Chikhi et al., IndexThePlanet) is a bulk assembly of ~all of the SRA.
Per accession it publishes **unitigs** (`/u/`), **contigs** (`/c/`), and
**protein clusters**, on a no-auth S3 bucket:

```bash
# contigs for one accession
aws s3 cp s3://logan-pub/c/<ACC>/<ACC>.contigs.fa.zst . --no-sign-request
# or
wget https://s3.amazonaws.com/logan-pub/c/<ACC>/<ACC>.contigs.fa.zst
zstd -d <ACC>.contigs.fa.zst
```

Version tested: **Logan v1.2** (December 2025 SRA freeze, ~38 M accessions).
S3 is reachable directly from our HPC login node — no local round-trip needed.

---

## Test design

A controlled head-to-head where the **only** variable is the contig source.
For each of 3 deep, paired-end samples we run the *identical* downstream on:

- **`logan`** — Logan `/c/` contigs, and
- **`ours`** — our own MEGAHIT assembly of the same reads.

| sample (SRA) | our sample id | study | reads |
|---|---|---|---|
| ERR2726406 | CCMD15562448ST_11_0 | curatedMetagenomicData | 2.1 GB |
| ERR2726443 | CCMD38158721ST_11_0 | curatedMetagenomicData | 1.6 GB |
| SRR2912786 | CA24 | (public SRA) | 2.7 GB |

Downstream, run per (source × sample):

- **Viral recovery** — geNomad `end-to-end` → CheckV `end_to_end`.
- **MAG recovery** — minimap2 (`-ax sr`) → `jgi_summarize_bam_contig_depths`
  → MetaBAT2 → CheckM2.
- **Taxonomy** — GTDB-Tk `classify_wf` (r226) on the resulting MAGs.

**Tools / DBs (HPC, SLURM + singularity):** geNomad 1.11.0 + **DB v1.9**
(`genomad-db-latest`; note v1.7 is *incompatible* with geNomad 1.11 — see
Gotchas), CheckV 1.0.3 + checkv-db-v1.5, minimap2/CoverM 0.7.0 container,
MetaBAT2 2.17, CheckM2 1.1.0 + uniref100.KO, GTDB-Tk 2.4.1 + release226.

---

## Results

### 1. Coverage of the SRA
29 / 30 sampled cohort accessions were present in Logan (~97%). A few are
missing (404), e.g. `SRR2565770`, so Logan cannot be the *sole* source.

### 2. Contig-level fragmentation
Logan assembles **more total sequence** but at **much lower contiguity** — it
retains a long tail of short contigs. Long-contig content (>=5 kb, >=10 kb) is
comparable; the difference is almost all sub-5 kb.

Example (SRR2565980):

| metric | Logan | ours (MEGAHIT) |
|---|---|---|
| contigs | 182,027 | 24,187 |
| total | 87.0 Mb | 55.5 Mb |
| **N50** | **1,596 bp** | **11,899 bp** |
| >=5 kb | 2,603 | 1,991 |
| >=10 kb | 1,051 | 1,103 |

Same pattern on the 3 test samples (Logan / ours contig counts):
ERR2726406 472,998 / 78,455 · ERR2726443 450,597 / 64,935 ·
SRR2912786 348,394 / 63,767.

### 3. MAG recovery (CheckM2)  — the fragmentation cost lands here

| sample | src | bins | HQ (>=90/<5) | MQ+ (>=50/<10) | median completeness |
|---|---|---|---|---|---|
| ERR2726406 | logan | 32 | 0 | 17 | 53% |
| ERR2726406 | ours | 45 | 7 | 25 | 69% |
| ERR2726443 | logan | 34 | 0 | 8 | 27% |
| ERR2726443 | ours | 51 | 9 | 24 | 66% |
| SRR2912786 | logan | 32 | 1 | 12 | 34% |
| SRR2912786 | ours | 48 | 5 | 17 | 38% |
| **total** | **logan** | **98** | **1** | **37** | |
| **total** | **ours** | **144** | **21** | **66** | |

**~20x fewer high-quality MAGs from Logan (1 vs 21).** Contamination is
uniformly low (<1%) in both — Logan does not make chimeric bins, it makes
*incomplete* ones. That is the signature of low contiguity: MetaBAT2 cannot
stitch a genome from short contigs.

### 4. Taxonomic concordance (GTDB-Tk r226, medium+ MAGs)

Where Logan and ours overlap, classification is identical — Logan does not
recover *different* organisms, it recovers a **near-strict subset**. Per sample,
genera among medium+ MAGs:

| sample | logan genera | ours genera | shared | ours-only | logan-only |
|---|---|---|---|---|---|
| ERR2726406 | 15 | 22 | 15 | 7 | 0 |
| ERR2726443 | 7 | 21 | 7 | 14 | 0 |
| SRR2912786 | 11 | 15 | 8 | 7 | 3 |

Species over the 3 samples: **Logan 37, ours 65, shared 32, ours-only 33,
logan-only 5.** In 2 of 3 samples every Logan genus is a subset of ours (zero
logan-only).

And what Logan drops is **not random — it loses core gut genera**:
*Bacteroides, Faecalibacterium, Roseburia, Alistipes, Parabacteroides,
Gemmiger, Phascolarctobacterium* ... (ERR2726443 collapses from 21 genera to 7).
So a MAG-based community-composition read-out from Logan would be systematically
depleted and biased, not merely noisier — the missing genomes fell below the
medium+ quality bar because fragmentation left their bins incomplete.

### 5. Viral recovery (geNomad + CheckV)

Totals across the 3 samples:

| src | viral contigs | >=5 kb | >=10 kb | viral Mb | CheckV Complete/High/Med/Low/NA |
|---|---|---|---|---|---|
| logan | **3,673** | 501 | 219 | 10.0 | 2 / 10 / 32 / 2568 / 1061 |
| ours | 2,271 | **618** | **303** | **13.4** | **4 / 45 / 49** / 1666 / 507 |

More nuanced than "viral is fine". Logan flags **more** viral candidates
(3,673 vs 2,271) — that is its short-contig tail, and it is genuinely useful for
**detection / screening** (presence of a phage in a library). But for a
**high-quality viral genome catalogue** our assembly wins on every axis that
matters: more long viral contigs (>=10 kb: 303 vs 219), more viral sequence
(13.4 vs 10.0 Mb), and **~4x more Complete + High-quality genomes (49 vs 12)**.
Logan's viral set is dominated by Low-quality / Not-determined fragments.

So even viral discovery pays a fragmentation tax at the **complete-genome** end —
just far less catastrophic than binning, and it is exactly the tax that
graph-aware recovery (previous section) is designed to pay down.

---

## Verdict

Logan is *fragmented contigs over an intact graph*. Whether it can replace
assembly depends entirely on whether your downstream reads contigs **literally**
or reads the **graph / targets a gene** — and it never costs you assembly compute
and covers ~97% of the SRA.

- **Genome-resolved binning / MAGs — do not use Logan.** ~20x fewer high-quality
  MAGs (1 vs 21); the recovered taxa are a *subset* of a de-novo assembly and it
  systematically loses core gut genera. Disqualifying for MAG-based work.
- **High-quality viral genome catalogue — Logan is ~4x weaker** at Complete/High
  genomes (12 vs 49), though it flags *more* candidates overall. Fine for viral
  **detection / screening**, weak for representative genomes.
- **Targeted gene / marker recovery — Logan is strong, and this is the point.**
  The graph is intact, so graph-aware tools (Pathracer against the reconstructed
  GFA) recover genes fragmentation would otherwise split — at a scale (tens of
  thousands of SRA libraries) that de-novo assembly cannot touch. This is Logan's
  killer app, not a consolation prize.

One line: **read Logan as contigs and it disappoints; read it as a graph, or mine
it for genes at SRA scale, and it is transformative.**

---

## Fragmented contigs, intact graph -> graph-aware recovery (Pathracer)

The binning result measures the cost of reading Logan as *isolated contigs*. But
Logan contigs are not isolated: **every header carries the de Bruijn graph**.

```
>SRR2565980_0 ka:f:17.725   L:-:179254:-  L:+:13037:+ L:+:13038:+
```

- `ka:f:17.725` — average k-mer abundance (coverage) of the contig.
- `L:-:179254:-` — a **link** (edge) to contig `179254` in the compacted de
  Bruijn graph, with the orientation of each end (`+`/`-`). GFA `L`-line
  semantics, embedded per record; the target is the bare numeric contig id.

Contigs are the graph **nodes**; the `L:` tags are the **edges**. In our example
sample, **42% of contigs carry >=1 link** (195,092 directed edges). The sequence
continuity a de-novo assembler would have concatenated into one long contig is
not lost in Logan — it is stored as *adjacency* instead of *concatenation*. So a
gene split across two or three contigs by a coverage dip, a repeat, or strain
variation ("fragmented but significant" hits) is still a single **path** through
the graph.

**Pathracer** (Shlemov & Korobeynikov, RECOMB-Seq 2019) exploits exactly this: it
aligns a profile HMM against the assembly graph rather than against contigs,
"racing" high-scoring paths and returning the best traversals that match the
pHMM. That recovers the full-length gene (e.g. HPV **L1** for typing) across the
breaks that defeat a contig-level blastn/nhmmer screen.

Reconstruct the graph from the Logan FASTA, then race an HMM against it:

```bash
# 1. Logan contigs.fa -> GFA. Logan sequences are single-line; nodes are named by
#    the numeric contig id so they match the bare-number targets in the L: tags.
#    Set K to Logan's contig k-mer size (overlap = K-1; confirm from Logan docs).
K=31; OVL="$((K-1))M"
awk -v ovl="$OVL" '
  /^>/ {
    name=substr($1,2); sub(/^.*_/,"",name)          # numeric id, matches L: targets
    cov="0"; delete lk; n=0
    for (i=2;i<=NF;i++){
      if ($i ~ /^ka:f:/){ split($i,a,":"); cov=a[3] }
      else if ($i ~ /^L:/){ lk[++n]=$i }
    }
    getline seq                                      # sequence is the next line
    printf "S\t%s\t%s\tLN:i:%d\tKC:f:%s\n", name, seq, length(seq), cov
    for (j=1;j<=n;j++){ split(lk[j],p,":")            # L : end1 : target : end2
      printf "L\t%s\t%s\t%s\t%s\t%s\n", name, p[2], p[3], p[4], ovl }
  }
' SRR2565980.contigs.fa > SRR2565980.gfa

# 2. Recover the fragmented gene by racing an HMM through the graph. NOTE: a
#    hand-rolled GFA does NOT feed generic PathRacer -- Logan is a compacted dBG
#    (overlap k-1, even) while PathRacer wants a SPAdes graph (overlap = k, odd).
#    The working recipe (Shen et al., see demo/) uses logan-extract to rebuild
#    the subgraph, then PathRacer from SPAdes v4.0.0:
logan-extract -d 1000 -t 8   acc.contigs.fa  hit_contig_list.txt  -o sub.gfa
pathracer --length 0.95 --top 5 --threads 16 --rescore   jrHMM.hmm  sub.gfa  --output pr_out
# keep paths with >95% HMM coverage = recovered full-length L1 -> type vs PaVE
```

A **runnable demo** of this (the L1-fragmentation problem live, plus the exact
recovery recipe) is in [`demo/`](demo/).

**Why this is Logan's killer app.** De-novo assembling tens of thousands of SRA
libraries just to mine one gene is prohibitive; pulling their Logan contigs
(~25 MB each), rebuilding the graph from the `L:` links, and racing an L1 pHMM is
not. Shen et al. ("Petabase-scale Papillomavirus Discovery", bioRxiv
2026.04.21.719858) did exactly this: from **18,642** Logan libraries with
fragmented-but-significant L1, graph-aware recovery (`logan-extract` + PathRacer)
retrieved **7,275 additional full-length L1s (51.7%)** with zero assembly of their
own. It also rebuts the "is the signal just a fragmentation artifact?" objection,
since the gene is recovered *through* the graph rather than hoping the assembler
stitched it. See [`demo/`](demo/) for the recipe.

**Caveats.** Building a valid GFA needs Logan's k (for the overlap) and correct
handling of the `L:±` orientations. Tangled regions (high-coverage repeats,
co-occurring strains) can blow up the path search or return many near-ties, and a
recovered path can chimerise co-occurring strains — validate each path against
its `ka:f:` coverage consistency and, ideally, the reads before trusting a type
call. Note this rescues *targeted gene / marker* recovery, **not binning** — a MAG
needs whole-genome linkage and a graph-aware binner, a much harder problem, so the
1-vs-21 MAG result stands.

---

## Reproduce

On an HPC cluster (SLURM + singularity). Submit from your working dir so the
relative `logs/` path in each `#SBATCH -o` resolves; edit the `$HOME/...` DB,
read, container, and `-B /path/to/hostfs` bind paths to match your site:

```bash
bash scripts/stage.sh              # download Logan + decompress our contigs; writes inputs.tsv
sbatch scripts/sanity.sbatch       # container/DB smoke test
sbatch scripts/genomad_checkv.sbatch   # array 1-6: viral recovery
BIN=$(sbatch --parsable scripts/binning.sbatch)        # array 1-6: MAGs
sbatch --dependency=aftercorr:$BIN scripts/gtdbtk.sbatch  # array 1-6: taxonomy
python3 scripts/checkm2_compare.py   # MAG-quality table
python3 scripts/virome_compare.py    # viral-recovery table
python3 scripts/gtdb_compare.py      # taxonomy concordance
```

`inputs.tsv` is 6 rows (`source<TAB>acc<TAB>sample`); each array task indexes
into it by `SLURM_ARRAY_TASK_ID`.

## Gotchas found while building this

- **geNomad 1.11.0 needs DB >= v1.9.** With DB v1.7 it crashes in
  `get_marker_annotation` (`ValueError: invalid literal for int(): '1398618at2'`)
  after only the `annotate` step, producing no summary.
- **CoverM 0.7.0's container** can't self-detect its bundled `minimap2`
  (`which --tty-only` hits busybox `which`); call `minimap2 -ax sr` directly
  from the same container instead.

## Layout

```
scripts/    stage + sbatch harnesses + comparison scripts
results/    comparison tables + per-sample summary TSVs
            (per_sample/<src>_<acc>/ : checkm2, gtdbtk, genomad, checkv summaries)
```

Large inputs are **not** committed: our MEGAHIT contigs (~150 MB each) are
regenerable, and Logan contigs are re-downloadable from S3 (command above).
