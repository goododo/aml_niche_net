# input : gene column of ws_s5_de_T.csv (names only), ambient table, HALLMARK v2024.1, repo MSigDB v2026.1 (Gavish sets),
#         Savage delay set, project IEG list (config_ws.R), 10x 2020-A GTF (gene_type)
# output: printed membership counts (outcome-independent: no logFC/P/FDR/hit/direction column is read)
# what it does: universe sizes, set coverage and handling fractions for the G07 draft
import csv, gzip, re
R='/FAST/gr10634/gaozy/aml_niche_net/results/tables/ws/'
with open(R+'ws_s5_de_T.csv') as f:
    rd=csv.reader(f); h=next(rd); gi=h.index('gene'); U=[r[gi] for r in rd]
U=set(U); print('T universe', len(U))
BL={}
with open(R+'decontx/GSE185381/ws_s1_ambient_blacklist.csv') as f:
    rd=csv.DictReader(f)
    for r in rd: BL[r['gene']]=int(float(r['n_units']))
ext={g for g,n in BL.items() if n>=2}; UA=U-ext
print('extended ambient', len(ext), 'in T universe', len(U&ext), 'ambient-free universe', len(UA))
def gmt(p, pat=None):
    S={}
    for line in open(p):
        x=line.rstrip('\n').split('\t')
        if pat is None or re.match(pat,x[0]): S[x[0]]=set(x[2:])
    return S
HM=gmt('/LARGE1/gr10634/gaozy/reference/h.all.v2024.1.Hs.symbols.gmt')
GV=gmt('/FAST/gr10634/gaozy/aml_niche_net/LCC_proj/scripts/panels/msigdb.v2026.1.Hs.symbols.gmt', r'GAVISH_3CA_METAPROGRAM_CD[48]_T_CELLS_')
SV=set()
with open('/LARGE1/gr10634/gaozy/reference/savage2021_processing_delay/savage2021_T_delay18h_consistent.csv') as f:
    rd=csv.DictReader(f); SV={r['gene'] for r in rd}
IEG={"FOS","FOSB","JUN","JUNB","JUND","EGR1","EGR2","EGR3","IER2","IER3","IER5","ATF3","DUSP1","ZFP36","NR4A1","KLF2","KLF4","KLF6","HSPA1A","HSPA1B","HSPA8","HSPH1","DNAJB1","PPP1R15A","SOCS3","CEBPD"}
H=(SV|IEG)&UA; print('Savage', len(SV), 'IEG', len(IEG), 'H in ambient-free universe', len(H), 'base rate %.3f'%(len(H)/len(UA)))
for name,S in (('HALLMARK',HM),('GAVISH_T',GV)):
    rows=[]
    for k,s in S.items():
        c=s&UA; ch=s&U
        rows.append((k, len(s), len(ch), len(c), len(c&H)/len(c) if c else float('nan')))
    n15=sum(1 for r in rows if r[3]>=15)
    print(f'{name}: {len(rows)} sets; >=15 genes in the ambient-free universe: {n15}')
    for r in sorted(rows, key=lambda r:-r[4]):
        if name=='GAVISH_T' or r[4]>=0.25 or r[3]<15:
            print('  %-60s size %3d inT %3d inUA %3d fracH %.2f%s'%(r[0],r[1],r[2],r[3],r[4],'  <15' if r[3]<15 else ''))
CLOCK=["ARNTL","CLOCK","NPAS2","PER1","PER2","PER3","CRY1","CRY2","NR1D1","NR1D2","DBP","TEF","HLF","CIART","BHLHE40","BHLHE41"]
GC=["FKBP5","TSC22D3","ZBTB16","KLF9","DDIT4","PER1","SGK1","DUSP1","TXNIP","CXCR4"]
print('clock in UA', [g for g in CLOCK if g in UA]); print('GC in UA', [g for g in GC if g in UA])
gt={}
with open('/LARGE1/gr10634/gaozy/reference/refdata-gex-GRCh38-2020-A/genes/genes.gtf') as f:
    for line in f:
        if line.startswith('#'): continue
        x=line.split('\t')
        if x[2]!='gene': continue
        m=re.search(r'gene_name "([^"]+)"',x[8]); t=re.search(r'gene_type "([^"]+)"',x[8])
        gt.setdefault(m.group(1), set()).add(t.group(1))
from collections import Counter
cU=Counter(); cA=Counter(); miss=0
for g in U:
    t=gt.get(g)
    if t is None: miss+=1; continue
    cU[ 'lncRNA' if t=={'lncRNA'} else ('protein_coding' if 'protein_coding' in t else 'other')]+=1
for g in UA:
    t=gt.get(g)
    if t is None: continue
    cA[ 'lncRNA' if t=={'lncRNA'} else ('protein_coding' if 'protein_coding' in t else 'other')]+=1
print('gene_type in T universe', dict(cU), 'unmatched names', miss); print('gene_type in ambient-free universe', dict(cA))
