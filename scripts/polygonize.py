import numpy as np, pickle, collections
from rasterio import features
from shapely.geometry import shape
from shapely.ops import unary_union
X0,Y0=2300,2800
lab=np.load('/home/claude/work/labels_hi.npy'); area=np.load('/home/claude/work/area_hi.npy')
recs=[]
for i,l in enumerate(open('/home/claude/work/labels.tsv',encoding='utf-8')):
    l=l.strip()
    if l:
        c,t,x,y=l.split('|'); recs.append((i+1,c,t,float(x),float(y)))
H,W=lab.shape
def region_at(x,y):
    hx,hy=int(round(x*2-X0)),int(round(y*2-Y0))
    if not(0<=hx<W and 0<=hy<H): return 0
    v=lab[hy,hx]
    if v: return int(v)
    r=14; sub=lab[max(hy-r,0):hy+r+1,max(hx-r,0):hx+r+1]
    vals=sub[sub>0]
    return int(np.bincount(vals).argmax()) if len(vals) else 0
own=collections.defaultdict(list); hyd=collections.defaultdict(list)
for i,c,t,x,y in recs:
    if c=='OWN': own[region_at(x,y)].append(i)
    if c=='HYD' and t.startswith('Laguna'): hyd[region_at(x,y)].append(i)
own.pop(0,None)
keep=np.zeros(len(area),bool)
for r in range(1,len(area)):
    a=area[r]
    if r in own: keep[r]=a<=1_200_000 and a>=150
    else: keep[r]=400<=a<=150_000
print('owner regions',len(own),'kept',keep.sum(),'owners w/ region',sum(len(v) for k,v in own.items() if keep[k]))
lab2=np.where(keep[lab],lab,0).astype(np.int32); del lab
polys=collections.defaultdict(list)
for geom,val in features.shapes(lab2,mask=lab2>0,connectivity=4):
    polys[int(val)].append(shape(geom))
out=[]
for r,ps in polys.items():
    g=max(ps,key=lambda p:p.area) if len(ps)>1 else ps[0]
    g=g.simplify(1.6,preserve_topology=True)
    if g.is_empty or g.area<150: continue
    # to AGN px
    from shapely.affinity import affine_transform
    g=affine_transform(g,[0.5,0,0,0.5,X0/2,Y0/2])
    out.append(dict(region=r,geom=g,owners=own.get(r,[]),area_px=float(area[r])))
pickle.dump(out,open('/home/claude/work/parcels.pkl','wb'))
print('polygons',len(out),'with owner',sum(1 for o in out if o['owners']),'multi-owner',sum(1 for o in out if len(o['owners'])>1))
