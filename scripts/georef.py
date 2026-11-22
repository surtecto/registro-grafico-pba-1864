import csv, numpy as np
from pyproj import Transformer
from scipy.interpolate import RBFInterpolator
T=Transformer.from_crs(4326,3857,always_xy=True); Ti=Transformer.from_crs(3857,4326,always_xy=True)
rows=list(csv.DictReader(open('/home/claude/work/gcps.csv')))
PX=np.array([[float(r['px']),float(r['py'])] for r in rows])
LL=np.array([[float(r['lon']),float(r['lat'])] for r in rows])
MX=np.c_[T.transform(LL[:,0],LL[:,1])]
def dp(p):
    x,y=(p[:,0]-4250)/4000,(p[:,1]-5400)/5000
    return np.c_[np.ones_like(x),x,y,x*x,x*y,y*y]
def dm(m):
    x,y=(m[:,0]+6.6e6)/4e5,(m[:,1]+4.3e6)/5e5
    return np.c_[np.ones_like(x),x,y,x*x,x*y,y*y]
class Model:
    def __init__(s,src,dst,D,smooth):
        s.D=D; s.c,*_=np.linalg.lstsq(D(src),dst,rcond=None)
        res=dst-D(src)@s.c
        s.r=RBFInterpolator(src,res,kernel='thin_plate_spline',smoothing=smooth) if smooth is not None else None
    def __call__(s,p):
        p=np.atleast_2d(p).astype(float); out=s.D(p)@s.c
        if s.r is not None: out=out+s.r(p)
        return out
def build(smooth=None, idx=None):
    i=np.arange(len(PX)) if idx is None else idx
    return Model(PX[i],MX[i],dp,smooth), Model(MX[i],PX[i],dm,smooth)
def px2ll(p,fwd):
    m=fwd(p); lon,lat=Ti.transform(m[:,0],m[:,1]); return np.c_[lon,lat]
if __name__=='__main__':
    for sm in (None,1e3,1e2,10,1,0):
        loo=[]
        for k in range(len(PX)):
            idx=np.array([j for j in range(len(PX)) if j!=k]); f,_=build(sm,idx)
            r=f(PX[k:k+1])[0]-MX[k]; loo.append(np.hypot(*r)*0.824)
        loo=np.array(loo)
        print(sm, 'LOO RMS km',round(np.sqrt((loo**2).mean())/1000,2),'median',round(np.median(loo)/1000,2))
