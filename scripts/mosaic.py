import cv2, numpy as np, json
S=2.0
agn=cv2.imread('/home/claude/work/agn.tif')
H0,W0=agn.shape[:2]; W,H=int(W0*S),int(H0*S)
hi=cv2.resize(agn,(W,H),interpolation=cv2.INTER_CUBIC); del agn
best=np.zeros((H,W),np.uint16)
src=np.zeros((H,W),np.uint8)
reg=json.load(open('/home/claude/work/reg/H.json'))
D='/mnt/user-data/uploads/pciabuenosaires1864/print/'
Sm=np.diag([S,S,1.0]); used=[]
for i,(name,r) in enumerate(sorted(reg.items())):
    Hm=np.array(r['H']); sc=np.sqrt(abs(np.linalg.det(Hm[:2,:2])))
    if r['inl']<60 or not(0.5<=sc<=0.585) or name in('E6','E7'): continue
    t=cv2.imread(f'{D}{name}.jpg'); th,tw=t.shape[:2]
    M=Sm@Hm
    cor=cv2.perspectiveTransform(np.float32([[[0,0],[tw,0],[tw,th],[0,th]]]),M)[0]
    x0,y0=np.floor(cor.min(0)).astype(int); x1,y1=np.ceil(cor.max(0)).astype(int)
    x0,y0=max(x0,0),max(y0,0); x1,y1=min(x1,W),min(y1,H)
    if x1<=x0 or y1<=y0: continue
    T=np.array([[1,0,-x0],[0,1,-y0],[0,0,1.0]])@M
    w=cv2.warpPerspective(t,T,(x1-x0,y1-y0),flags=cv2.INTER_LINEAR)
    yy,xx=np.mgrid[0:th,0:tw]
    d=np.minimum(np.minimum(xx,tw-1-xx),np.minimum(yy,th-1-yy)).astype(np.float32)
    dw=cv2.warpPerspective(d,T,(x1-x0,y1-y0),flags=cv2.INTER_LINEAR).astype(np.uint16)
    b=best[y0:y1,x0:x1]; m=dw>b
    hi[y0:y1,x0:x1][m]=w[m]; b[m]=dw[m]; src[y0:y1,x0:x1][m]=i+1
    used.append(name)
print('used',len(used),used)
cv2.imwrite('/home/claude/work/hi.tif',hi,[cv2.IMWRITE_TIFF_COMPRESSION,5])
cv2.imwrite('/home/claude/work/hi_src.png',src)
cv2.imwrite('/home/claude/work/hi_small.jpg',cv2.resize(hi,(W//16,H//16),interpolation=cv2.INTER_AREA))
