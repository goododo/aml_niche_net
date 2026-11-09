# input : gene names only (ws_s5_de_T.csv gene column, ambient table, HALLMARK, Gavish, Savage with direction, IEG list)
# output: genes left per set after direction-matched removal (H_up = IEG + Savage-up; H_down = IEG + Savage-down) and after
#         full removal; outcome-independent (no logFC/P/FDR/hit/direction column of any S5 table is read)
import csv, re
R='/FAST/gr10634/gaozy/aml_niche_net/results/tables/ws/'
U={r['gene'] for r in csv.DictReader(open(R+'ws_s5_de_T.csv'))}
ext={r['gene'] for r in csv.DictReader(open(R+'decontx/GSE185381/ws_s1_ambient_blacklist.csv')) if float(r['n_units'])>=2}
UA=U-ext
def gmt(p,pat=None):
    S={}
    for l in open(p):
        x=l.rstrip('\n').split('\t')
        if pat is None or re.match(pat,x[0]): S[x[0]]=set(x[2:])
    return S
HM=gmt('/LARGE1/gr10634/gaozy/reference/h.all.v2024.1.Hs.symbols.gmt')
GV=gmt('/FAST/gr10634/gaozy/aml_niche_net/LCC_proj/scripts/panels/msigdb.v2026.1.Hs.symbols.gmt', r'GAVISH_3CA_METAPROGRAM_CD[48]_T_CELLS_')
SV={r['gene']:r['direction'] for r in csv.DictReader(open('/LARGE1/gr10634/gaozy/reference/savage2021_processing_delay/savage2021_T_delay18h_consistent.csv'))}
IEG={"FOS","FOSB","JUN","JUNB","JUND","EGR1","EGR2","EGR3","IER2","IER3","IER5","ATF3","DUSP1","ZFP36","NR4A1","KLF2","KLF4","KLF6","HSPA1A","HSPA1B","HSPA8","HSPH1","DNAJB1","PPP1R15A","SOCS3","CEBPD"}
Hup={g for g,d in SV.items() if d=='up'}|IEG; Hdn={g for g,d in SV.items() if d=='down'}|IEG; H=set(SV)|IEG
print('U_A',len(UA),'U_A-Hup',len(UA-Hup),'U_A-Hdown',len(UA-Hdn),'U_A-H',len(UA-H))
for name,S in (('HALLMARK',HM),('GAVISH',GV)):
    rows=[(k,len(s&UA),len(s&UA-Hup),len(s&UA-Hdn),len(s&UA-H)) for k,s in S.items() if len(s&UA)>=15]
    print(f'{name}: {len(rows)} sets with >=15 in U_A; min left: up-removal {min(r[2] for r in rows)}, down-removal {min(r[3] for r in rows)}, full {min(r[4] for r in rows)}')
    print('  sets <10 after up-removal:', [ (r[0],r[2]) for r in rows if r[2]<10]); print('  sets <10 after down-removal:', [(r[0],r[3]) for r in rows if r[3]<10])
    print('  sets <15 after full removal:', [(r[0].replace('GAVISH_3CA_METAPROGRAM_',''),r[4]) for r in rows if r[4]<15])
OPP=['ANXA1','APOBEC3G','CCL5','CD27','CKLF','GZMA','GZMH','GZMK','HOPX','KLRD1','KLRG1','NKG7','TIGIT','TRGC2']
print('opposite-sign class in U_A:', sum(g in UA for g in OPP), 'Savage-down:', sum(SV.get(g)=='down' for g in OPP))
