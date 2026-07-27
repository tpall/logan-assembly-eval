import glob, os
BASE=os.path.expanduser("~/Projects/logan-test/virome")
accs=["ERR2726406","ERR2726443","SRR2912786"]

def genomad(src,acc):
    hits=glob.glob(f"{BASE}/{src}_{acc}/genomad/**/*_virus_summary.tsv", recursive=True)
    if not hits: return None
    lens=[]
    with open(hits[0]) as f:
        h=f.readline().rstrip("\n").split("\t"); li=h.index("length")
        for line in f:
            p=line.rstrip("\n").split("\t")
            if len(p)>li: lens.append(int(p[li]))
    return lens

def checkv(src,acc):
    q=f"{BASE}/{src}_{acc}/checkv/quality_summary.tsv"
    tiers={}
    if not os.path.exists(q): return None
    with open(q) as f:
        h=f.readline().rstrip("\n").split("\t"); qi=h.index("checkv_quality")
        for line in f:
            p=line.rstrip("\n").split("\t")
            if len(p)>qi: tiers[p[qi]]=tiers.get(p[qi],0)+1
    return tiers

order=["Complete","High-quality","Medium-quality","Low-quality","Not-determined"]
print("%-20s%-7s%8s%8s%9s%11s | %s"%("sample","src","nVir","n>=5kb","n>=10kb","viral_Mb","CheckV Comp/High/Med/Low/NA"))
agg={"logan":[0,0,0,0.0,[0,0,0,0,0]],"ours":[0,0,0,0.0,[0,0,0,0,0]]}
for a in accs:
    for src in ("logan","ours"):
        lens=genomad(src,a); tiers=checkv(src,a)
        if lens is None:
            print("%-20s%-7s  (geNomad not finished)"%(a,src)); continue
        n=len(lens); n5=sum(1 for L in lens if L>=5000); n10=sum(1 for L in lens if L>=10000); mb=sum(lens)/1e6
        tv=[ (tiers or {}).get(t,0) for t in order]
        print("%-20s%-7s%8d%8d%9d%11.2f | %d/%d/%d/%d/%d"%(a,src,n,n5,n10,mb,*tv))
        agg[src][0]+=n; agg[src][1]+=n5; agg[src][2]+=n10; agg[src][3]+=mb
        for i in range(5): agg[src][4][i]+=tv[i]
print("-"*95)
for src in ("logan","ours"):
    n,n5,n10,mb,tv=agg[src]
    print("%-20s%-7s%8d%8d%9d%11.2f | %d/%d/%d/%d/%d"%("TOTAL (3 samples)",src,n,n5,n10,mb,*tv))
