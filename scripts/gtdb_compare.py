import glob, os
BASE=os.path.expanduser("~/Projects/logan-test")

def rank(cl, pfx):
    for t in cl.split(";"):
        t=t.strip()
        if t.startswith(pfx) and len(t)>len(pfx): return t[len(pfx):]
    return None

def load_quality(src,acc):
    q=f"{BASE}/binning/{src}_{acc}/checkm2/quality_report.tsv"
    d={}
    if not os.path.exists(q): return d
    with open(q) as f:
        h=f.readline().rstrip("\n").split("\t"); ci=h.index("Completeness"); ki=h.index("Contamination")
        for line in f:
            p=line.rstrip("\n").split("\t")
            d[p[0]]=(float(p[ci]),float(p[ki]))
    return d

def load_gtdb(src,acc):
    d={}
    for s in glob.glob(f"{BASE}/gtdbtk/{src}_{acc}/**/gtdbtk.*.summary.tsv", recursive=True)+\
             glob.glob(f"{BASE}/gtdbtk/{src}_{acc}/gtdbtk.*.summary.tsv"):
        with open(s) as f:
            h=f.readline().rstrip("\n").split("\t"); gi=h.index("classification")
            for line in f:
                p=line.rstrip("\n").split("\t"); d[p[0]]=p[gi]
    return d

accs=["ERR2726406","ERR2726443","SRR2912786"]
print("%-20s%-7s%8s%9s%10s%11s"%("sample","src","MQ+bins","classfd","genera","species"))
tot={"logan":[set(),set()],"ours":[set(),set()]}
per={}
for a in accs:
    for src in ("logan","ours"):
        qual=load_quality(src,a); gt=load_gtdb(src,a)
        # medium+ bins only
        mq=[b for b,(c,k) in qual.items() if c>=50 and k<10]
        gens=set(); sps=set()
        for b in mq:
            cl=gt.get(b) or gt.get(b.replace(".fa",""))
            if not cl: continue
            g=rank(cl,"g__"); s=rank(cl,"s__")
            if g: gens.add(g)
            if s: sps.add(s)
        per[(a,src)]=(gens,sps)
        classified=sum(1 for b in mq if (gt.get(b) or gt.get(b.replace(".fa",""))))
        print("%-20s%-7s%8d%9d%10d%11d"%(a,src,len(mq),classified,len(gens),len(sps)))
        tot[src][0]|={a+"|"+x for x in gens}; tot[src][1]|={a+"|"+x for x in sps}
print("-"*70)
# per-sample genus overlap
print("\nGenus recovery among MQ+ MAGs (Logan vs ours), per sample:")
for a in accs:
    gl,_=per[(a,"logan")]; go,_=per[(a,"ours")]
    shared=gl&go
    print("  %-14s shared=%2d  logan_only=%2d  ours_only=%2d  (logan %d / ours %d genera)"
          %(a,len(shared),len(gl-go),len(go-gl),len(gl),len(go)))
    if go-gl: print("     genera ONLY in ours (Logan missed): "+", ".join(sorted(go-gl)))
    if gl-go: print("     genera ONLY in Logan: "+", ".join(sorted(gl-go)))
# aggregate species
sl=tot["logan"][1]; so=tot["ours"][1]
print("\nSpecies-level (union over 3 samples): logan=%d ours=%d shared=%d ours_only=%d logan_only=%d"
      %(len(sl),len(so),len(sl&so),len(so-sl),len(sl-so)))
