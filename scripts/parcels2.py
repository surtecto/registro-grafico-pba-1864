import sys, numpy as np, cv2, gc, collections
sys.path.insert(0,'/home/claude/scripts')
from extend import extend
X0,Y0=2300,2800
B=np.load('/home/claude/work/barrier_hi.npy'); H,W=B.shape
for it in range(2):
    E=np.zeros_like(B); S=3000; O=260
    for y in range(0,H,S):
        for x in range(0,W,S):
            ya,xa=max(y-O,0),max(x-O,0); yb,xb=min(y+S+O,H),min(x+S+O,W)
            e=extend(np.ascontiguousarray(B[ya:yb,xa:xb]))
            E[y:min(y+S,H),x:min(x+S,W)]|=e[y-ya:y-ya+min(S,H-y),x-xa:x-xa+min(S,W-x)]
    print('pass',it,'added px',int(E.sum()),flush=True); B|=E; del E
np.save('/home/claude/work/barrier2_hi.npy',B)
def label(seed, maxgrow, restrict=None):
    n,lab=cv2.connectedComponents(seed,connectivity=4,ltype=cv2.CV_32S)
    grown=np.zeros_like(lab); S2=5000; O2=maxgrow+8
    for y in range(0,H,S2):
        for x in range(0,W,S2):
            ya,xa=max(y-O2,0),max(x-O2,0); yb,xb=min(y+S2+O2,H),min(x+S2+O2,W)
            sl=lab[ya:yb,xa:xb]; sd=(sl>0).astype(np.uint8)
            if not sd.any(): continue
            if sd.all(): g2=sl.copy()
            else:
                d2,idx=cv2.distanceTransformWithLabels(1-sd,cv2.DIST_L2,3,labelType=cv2.DIST_LABEL_PIXEL)
                lut=np.concatenate([[0],sl[sl>0]]).astype(np.int32); g2=lut[idx]; g2[d2>maxgrow]=0
            grown[y:min(y+S2,H),x:min(x+S2,W)]=g2[y-ya:y-ya+min(S2,H-y),x-xa:x-xa+min(S2,W-x)]
    return n,grown
dist=cv2.distanceTransform(1-B,cv2.DIST_L2,3); del B; gc.collect()
n1,lab=label((dist>3.5).astype(np.uint8),9)
area=np.bincount(lab.ravel(),minlength=n1)
own=collections.Counter()
for l in open('/home/claude/work/labels.tsv',encoding='utf-8'):
    l=l.strip()
    if not l: continue
    c,t,x,y=l.split('|')
    if c!='OWN': continue
    hx,hy=int(round(float(x)*2-X0)),int(round(float(y)*2-Y0))
    if 0<=hx<W and 0<=hy<H and lab[hy,hx]: own[int(lab[hy,hx])]+=1
prob=np.zeros(n1,bool); prob[area>300_000]=True
for r,c in own.items():
    if c>=2: prob[r]=True
prob[0]=False
print('regions',n1,'problem',int(prob.sum()),'problem area frac',float(area[prob].sum()/area.sum()),flush=True)
pm=prob[lab]
seed2=((dist>14)&pm).astype(np.uint8); del dist; gc.collect()
n2,lab2=label(seed2,32); del seed2; gc.collect()
lab2[~pm]=0
lab[pm]=0; m=lab2>0; lab[m]=lab2[m]+n1; del lab2,m,pm; gc.collect()
area=np.bincount(lab.ravel())
np.save('/home/claude/work/labels_hi.npy',lab); np.save('/home/claude/work/area_hi.npy',area)
print('final regions',len(area),'>=400:',int((area>=400).sum()),'>300k:',int((area[1:]>300_000).sum()))
