import sys, json, numpy as np, cv2, gc, pickle
sys.path.insert(0,'/home/claude/scripts')
from seg4 import barrier
from PIL import Image; Image.MAX_IMAGE_PIXELS=None
X0,Y0,X1,Y1=2300,2800,16800,21300   # hi-frame land bbox
g=np.array(Image.open('/home/claude/work/hi.tif').crop((X0,Y0,X1,Y1)).convert('L'))
H,W=g.shape; B=np.zeros((H,W),np.uint8); S=4000; O=120
for y in range(0,H,S):
    for x in range(0,W,S):
        ya,xa=max(y-O,0),max(x-O,0); yb,xb=min(y+S+O,H),min(x+S+O,W)
        b=barrier(g[ya:yb,xa:xb],6.0)
        B[y:min(y+S,H),x:min(x+S,W)]=b[y-ya:y-ya+min(S,H-y),x-xa:x-xa+min(S,W-x)]
    print('row',y,flush=True)
del g; gc.collect()
# mask non-map boxes (AGN px -> hi)
for (ax0,ay0,ax1,ay1) in [(6000,3850,8350,6350),(250,9050,2750,10650)]:
    B[max(ay0*2-Y0,0):max(ay1*2-Y0,0),max(ax0*2-X0,0):max(ax1*2-X0,0)]=1
free=(1-B).astype(np.uint8)
dist=cv2.distanceTransform(free,cv2.DIST_L2,3)
seed=(dist>3.5).astype(np.uint8); del dist; gc.collect()
n,lab=cv2.connectedComponents(seed,connectivity=4,ltype=cv2.CV_32S)
print('seeds',n,flush=True)
np.save('/home/claude/work/barrier_hi.npy',B); del B,free; gc.collect()
grown=np.zeros_like(lab); S2=5000; O2=24
for y in range(0,H,S2):
    for x in range(0,W,S2):
        ya,xa=max(y-O2,0),max(x-O2,0); yb,xb=min(y+S2+O2,H),min(x+S2+O2,W)
        sl=lab[ya:yb,xa:xb]; sd=(sl>0).astype(np.uint8)
        if sd.all(): g2=sl.copy()
        elif not sd.any(): continue
        else:
            d2,idx=cv2.distanceTransformWithLabels(1-sd,cv2.DIST_L2,3,labelType=cv2.DIST_LABEL_PIXEL)
            lut=np.concatenate([[0],sl[sl>0]]).astype(np.int32); g2=lut[idx]; g2[d2>9]=0
        grown[y:min(y+S2,H),x:min(x+S2,W)]=g2[y-ya:y-ya+min(S2,H-y),x-xa:x-xa+min(S2,W-x)]
del seed; gc.collect()
np.save('/home/claude/work/labels_hi.npy',grown)
area=np.bincount(grown.ravel()); print('regions',len(area),'>=400:',(area>=400).sum())
np.save('/home/claude/work/area_hi.npy',area)
