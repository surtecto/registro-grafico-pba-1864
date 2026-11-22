import numpy as np, cv2
def extend(B, G=120, skip=80, T=None):
    """Close gaps: prolong straight dangling line ends until they meet another line (max G px)."""
    H,W=B.shape; T=skip+G+2
    segs=cv2.HoughLinesP(B*255,1,np.pi/360,60,minLineLength=70,maxLineGap=6)
    if segs is None: return np.zeros_like(B)
    s=segs.reshape(-1,4).astype(np.float32)
    P=np.r_[s[:,2:4],s[:,0:2]]; Q=np.r_[s[:,0:2],s[:,2:4]]      # end, other end
    d=P-Q; L=np.hypot(d[:,0],d[:,1]); u=d/L[:,None]; n=np.c_[-u[:,1],u[:,0]]
    t=np.arange(0,T,dtype=np.float32)
    def samp(off_t,off_n):
        x=np.rint(P[:,0,None]+u[:,0,None]*off_t+n[:,0,None]*off_n).astype(np.int32)
        y=np.rint(P[:,1,None]+u[:,1,None]*off_t+n[:,1,None]*off_n).astype(np.int32)
        ok=(x>=0)&(x<W)&(y>=0)&(y<H)
        v=np.ones(x.shape,np.uint8); v[ok]=B[y[ok],x[ok]]   # outside = barrier
        return v
    c=samp(t[None,:],0)
    free=c==0
    has=free.any(1); i0=np.where(has,free.argmax(1),T)
    ok=has&(i0<=skip)&(i0>=1)
    # after i0: first barrier
    idx=np.arange(T)[None,:]
    after=(c==1)&(idx>i0[:,None])
    hit=after.any(1); j=np.where(hit,after.argmax(1),T)
    ok&=hit&((j-i0)<=G)&((j-i0)>=2)
    # dangling test: just before the end the line must be thin (no barrier at lateral ±6)
    for back in (2,4):
        tb=np.maximum(i0-back,0).astype(np.float32)[:,None]
        for lat in (-6,6):
            ok&=samp(tb,lat)[:,0]==0
    # forward corridor must be clear laterally for the first pixels (not running along another line)
    E=np.zeros_like(B)
    for k in np.nonzero(ok)[0]:
        a=P[k]+u[k]*(i0[k]-2); b=P[k]+u[k]*(j[k]+1)
        cv2.line(E,(int(round(a[0])),int(round(a[1]))),(int(round(b[0])),int(round(b[1]))),1,3)
    return E
