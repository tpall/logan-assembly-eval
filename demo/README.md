# Demo: graph-aware L1 recovery from Logan (Pathracer)

A concrete, runnable illustration of the capstone claim in the parent README:
**Logan's fragmentation is recoverable for targeted genes because the de Bruijn
graph is intact in the contig headers.** The worked example is papillomavirus
**L1** -- the gene used for PV typing -- which fragmentation routinely splits
across contigs.

This demo follows the method of **Shen et al., "Petabase-scale Papillomavirus
Discovery"** (bioRxiv 2026.04.21.719858; code:
[github.com/syueqiao/pvs_in_logan](https://github.com/syueqiao/pvs_in_logan)),
who recovered **7,275 additional full-length L1s (51.7%)** from **18,642** Logan
libraries whose L1 was fragmented but significant.

## The problem (this part runs live)

`run_l1_demo.sh` pulls a Logan assembly for an HPV-rich library (SRR2565980),
predicts genes (prodigal) and searches them with an L1 HMM (Pfam **PF00500**,
`Late_protein_L1`). A full L1 is ~500 aa, but on the isolated contigs it comes
back **shattered into partial hits on different contigs**:

```
L1 hit: SRR2565980_155501_1  aa=68    env=1-68
L1 hit: SRR2565980_74276_1   aa=155   env=1-154
L1 hit: SRR2565980_94628_1   aa=184   env=1-184
L1 hit: SRR2565980_12123_1   aa=220   env=1-219
L1 hit: SRR2565980_22571_1   aa=239   env=6-239
L1 hit: SRR2565980_34731_1   aa=239   env=6-239
L1 hit: SRR2565980_9691_2    aa=476   env=4-476
L1 hit: SRR2565980_1558_2    aa=516   env=4-515
```

No single contig-level search reconstructs the full gene from the short
fragments; but they are adjacent **paths** in the assembly graph.

## The graph is in the headers

The demo also rebuilds a GFA from the `L:` header links (182,027 nodes /
195,092 edges for this library) to show the graph substrate is present. **This
illustrative GFA does not itself feed PathRacer** -- see the next section.

## Why generic PathRacer does not eat a hand-rolled Logan GFA

I tried the obvious `awk` FASTA->GFA then bioconda `pathracer`. It does not work,
and the reason is a real convention clash, not a formatting slip:

- PathRacer wants a **SPAdes-style** graph: overlap = k, and **k must be odd**,
  with segments that truly overlap by k.
- Logan is a **compacted unitig dBG** (Minia3): overlap = **k-1 = 30** (even).

A compacted dBG's overlap is always k-1 -- the opposite parity to the odd k
PathRacer requires -- so it is structurally a different graph type. Concretely,
bioconda PathRacer 3.16 gave: `30M` -> `k-mer length must be odd`; `31M` -> loads
the graph and **seeds 11 L1-matching edges** but then aborts extending paths
(`position()+1 == len+k`); `*` -> links dropped, "can't determine k".

## The recovery that works (Shen et al.) -- verified end-to-end here

The paper bridges the gap with a purpose-built tool, **`logan-extract`**, from the
SPAdes `graph-extract` branch, plus PathRacer from SPAdes v4.0.0. Both were built
and run for this demo (2026-07-27). Two build/usage gotchas cost real time, so they
are recorded here:

```bash
# --- build (from the graph-extract branch) ---
git clone -b graph-extract https://github.com/ablab/spades && cd spades
PREFIX="$PWD/inst" ./spades_compile.sh -j 8               # logan-extract (needs zstd)
# PathRacer is OFF by default; enable it explicitly at a FRESH configure
# (SPADES_ENABLE_PROJECTS is empty by default -> pathracer is skipped):
mkdir bp && cd bp
cmake -DSPADES_ENABLE_PROJECTS="pathracer" -DCMAKE_INSTALL_PREFIX="$PWD/../inst" ../src
make -j 8 pathracer

# --- run (per library) ---
# seed list = HMM hit names; logan-extract keys edges off split("_")[1]
logan-extract  acc.contigs.fa  hit_list.txt  sub.gfa  -d 1000 -t 8 --gfa   # NB: no -k !
pathracer  L1.hmm  sub.gfa  --hmm --length 0.95 --top 5 --rescore -t 8 --output pr_out
```

**Gotcha 1 -- no `-k`.** `logan-extract` defaults to k=29, which rebuilds a *linked*
subgraph. Passing `-k 31` yields a segments-only GFA (0 links) and PathRacer dies
with "Failed to determine k-mer length".
**Gotcha 2 -- enable pathracer.** `make pathracer` fails with "No rule to make
target" unless `SPADES_ENABLE_PROJECTS="pathracer"` is set at a clean configure.

**Result on SRR2565980** (`demo/example_output/`): from a 12,532-segment / 19,812-link
subgraph, PathRacer recovered a **508 aa L1 covering positions 1-492 of the 497-aa
model (~99%), score 813**, along a path `Edges=9391_7597_18139'` -- i.e. **stitched
across three graph edges**, which no single contig held. The protein opens
`WLPASGKVYLPP...`, the canonical papillomavirus L1 N-terminus. The paper uses a
bespoke **jrHMM** (jellyroll HMM over the six L1 sheets B/CD/E/F/GH/I); PF00500 is the
public stand-in here. `run_l1_demo.sh` runs this step automatically when
`logan-extract` and `pathracer` are on PATH (or via `LOGAN_EXTRACT=`/`PATHRACER=`),
and prints the build recipe otherwise.

## Run it

```bash
mamba create -n l1demo -c bioconda -c conda-forge hmmer prodigal seqkit zstd
conda activate l1demo
bash run_l1_demo.sh SRR2565980      # steps 1-3 run live; step 4 runs if logan-extract is built
```

## Files

- `run_l1_demo.sh` -- the demo driver (download + GFA + baseline + recovery).
- `PF00500.hmm.gz` -- Pfam papillomavirus L1 HMM (`Late_protein_L1`, CC0), bundled.
- `example_output/` -- the real recovered L1 from SRR2565980: `recovered_L1.faa`
  (protein, top 5 paths), `recovered_L1.fna` (nucleotide), `recovered_L1.domtblout`
  (HMM coverage of the 497-aa model).
