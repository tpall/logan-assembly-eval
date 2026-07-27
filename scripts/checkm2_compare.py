import glob, os, statistics as st
base=os.path.expanduser("~/Projects/logan-test/binning")
rows={}
for q in sorted(glob.glob(base+"/*/checkm2/quality_report.tsv")):
    d=os.path.basename(os.path.dirname(os.path.dirname(q)))
    src,acc=d.split("_",1)
    comp=[]; cont=[]
    with open(q) as f:
        h=f.readline().rstrip("\n").split("\t")
        ci=h.index("Completeness"); ki=h.index("Contamination")
        for line in f:
            p=line.rstrip("\n").split("\t")
            comp.append(float(p[ci])); cont.append(float(p[ki]))
    hq=sum(1 for c,k in zip(comp,cont) if c>=90 and k<5)
    mq=sum(1 for c,k in zip(comp,cont) if c>=50 and k<10)   # MIMAG medium+ (incl HQ)
    rows[(acc,src)]=(len(comp),hq,mq,st.median(comp),st.median(cont))
print("%-20s%-7s%5s%13s%15s%10s%10s"%("sample","src","bins","HQ(>=90/<5)","MQ+(>=50/<10)","med_comp","med_cont"))
accs=sorted({a for a,_ in rows})
agg={"logan":[0,0,0],"ours":[0,0,0]}
for a in accs:
    for src in ("logan","ours"):
        n,hq,mq,mc,mk=rows[(a,src)]
        print("%-20s%-7s%5d%13d%15d%10.1f%10.1f"%(a,src,n,hq,mq,mc,mk))
        agg[src][0]+=n; agg[src][1]+=hq; agg[src][2]+=mq
print("-"*80)
for src in ("logan","ours"):
    n,hq,mq=agg[src]
    print("%-20s%-7s%5d%13d%15d"%("TOTAL (3 samples)",src,n,hq,mq))
