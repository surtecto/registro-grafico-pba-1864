import cv2, numpy as np, json, os, sys
from PIL import Image
Image.MAX_IMAGE_PIXELS=None
agn=cv2.imread('/home/claude/work/agn.tif',cv2.IMREAD_GRAYSCALE)
H,W=agn.shape
sift=cv2.SIFT_create(nfeatures=0,contrastThreshold=0.03)
# AGN features in blocks
kps=[];des=[]
B=2200
for y in range(0,H,B-200):
    for x in range(0,W,B-200):
        k,d=sift.detectAndCompute(agn[y:y+B,x:x+B],None)
        if d is None: continue
        kps+= [(p.pt[0]+x,p.pt[1]+y) for p in k]; des.append(d)
kps=np.array(kps,np.float32); des=np.vstack(des)
print('agn feats',len(kps),flush=True)
flann=cv2.FlannBasedMatcher(dict(algorithm=1,trees=6),dict(checks=64))
flann.add([des]); flann.train()
out={}
D='/mnt/user-data/uploads/pciabuenosaires1864/print/'
for r in 'ABCDEFGH':
    for c in range(1,9):
        f=f'{D}{r}{c}.jpg'
        if not os.path.exists(f): continue
        t=cv2.imread(f,cv2.IMREAD_GRAYSCALE); s=0.55
        ts=cv2.resize(t,None,fx=s,fy=s,interpolation=cv2.INTER_AREA)
        k,d=sift.detectAndCompute(ts,None)
        if d is None or len(k)<50: print(r,c,'few feats'); continue
        m=flann.knnMatch(d,k=2)
        good=[a for a,b in m if a.distance<0.72*b.distance]
        if len(good)<12: print(r,c,'few matches',len(good),flush=True); continue
        src=np.float32([k[g.queryIdx].pt for g in good])/s; dst=kps[[g.trainIdx for g in good]]
        Hm,mask=cv2.findHomography(src,dst,cv2.RANSAC,6.0)
        inl=int(mask.sum()) if mask is not None else 0
        if Hm is None: print(r,c,'no H'); continue
        cen=cv2.perspectiveTransform(np.float32([[[1200,t.shape[0]/2]]]),Hm)[0,0]
        sc=np.sqrt(abs(np.linalg.det(Hm[:2,:2])))
        print(r,c,'matches',len(good),'inl',inl,'center',cen.round(),'scale',round(float(sc),3),flush=True)
        out[f'{r}{c}']=dict(H=Hm.tolist(),inl=inl,h=t.shape[0])
json.dump(out,open('/home/claude/work/reg/H.json','w'))
