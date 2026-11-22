import numpy as np, cv2, sys
from PIL import Image
Image.MAX_IMAGE_PIXELS=None
def kernels(L=45,off=5,n=24):
    ks=[]
    for i in range(n):
        th=np.pi*i/n; c,s=np.cos(th),np.sin(th); R=L//2+off+2
        k=np.zeros((2*R+1,2*R+1),np.float32)
        def ln(o,val):
            ox,oy=-s*o,c*o
            p1=(int(round(R+ox-c*L/2)),int(round(R+oy-s*L/2))); p2=(int(round(R+ox+c*L/2)),int(round(R+oy+s*L/2)))
            m=np.zeros_like(k); cv2.line(m,p1,p2,1.0,1); return m/m.sum()*val
        k+=ln(0,1.0)+ln(off,-0.5)+ln(-off,-0.5); ks.append(k)
    return ks
KS=kernels()
def barrier(g,thr=7.0,minlen=60):
    g=g.astype(np.float32)
    bg=cv2.GaussianBlur(g,(0,0),12); dark=np.clip(bg-cv2.GaussianBlur(g,(0,0),0.8),0,60)
    R=np.zeros_like(dark)
    for k in KS: R=np.maximum(R,cv2.filter2D(dark,-1,k))
    B=((R>thr)&(dark>thr*0.6)).astype(np.uint8)
    n,lab,st,_=cv2.connectedComponentsWithStats(B,8)
    keep=np.maximum(st[:,2],st[:,3])>=minlen; keep[0]=False
    B=keep[lab].astype(np.uint8)
    return cv2.dilate(B,cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(5,5)))
if __name__=='__main__':
    x0,y0,x1,y1=[int(v)*2 for v in sys.argv[1:5]]
    im=Image.open('/home/claude/work/hi.tif').crop((x0,y0,x1,y1)); g=np.array(im.convert('L'))
    Ld=barrier(g,*[float(v) for v in sys.argv[5:6]])
    n2,lab2,st2,_=cv2.connectedComponentsWithStats(1-Ld,4)
    area=st2[:,4]; print('regions>1500:',(area>1500).sum(),'largest frac',area[1:].max()/area.sum())
    rng=np.random.default_rng(1); col=rng.integers(90,255,(n2,3)).astype(np.uint8); col[area<1500]=(255,255,255)
    out=col[lab2]; out[Ld>0]=0
    a=np.array(im); out=(0.6*out+0.4*a).astype(np.uint8)
    Image.fromarray(out).resize((1300,int(1300*(y1-y0)/(x1-x0)))).save('/home/claude/work/seg.jpg',quality=85)
