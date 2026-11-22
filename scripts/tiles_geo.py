import sys, os, json, math, numpy as np, cv2
sys.path.insert(0,'/home/claude/scripts')
from georef import build
from PIL import Image
Image.MAX_IMAGE_PIXELS=None
fwd,inv=build(0)
src=cv2.cvtColor(cv2.imread('/home/claude/work/hi.tif'),cv2.COLOR_BGR2BGRA)
K=2
H,W=src.shape[0]//K,src.shape[1]//K
mask=np.zeros(src.shape[:2],np.uint8); mask[90*K:10740*K,130*K:8370*K]=255
for poly in json.load(open('/home/claude/work/mask.json'))['mask_out']:
    cv2.fillPoly(mask,[np.array(poly,np.int32)*K],0)
src[:,:,3]=mask
OUT='/home/claude/site/tiles/geo'; R=6378137*math.pi
corners=fwd(np.array([[0,0],[W,0],[0,H],[W,H],[W/2,0],[W/2,H],[0,H/2],[W,H/2]]))
mx0,my0=corners.min(0); mx1,my1=corners.max(0)
ZMAX=int(sys.argv[1]) if len(sys.argv)>1 else 11
def tile_bounds(z,x,y):
    s=2*R/2**z; return -R+x*s, R-(y+1)*s, -R+(x+1)*s, R-y*s
z=ZMAX; s=2*R/2**z
tx0=int((mx0+R)//s); tx1=int((mx1+R)//s); ty0=int((R-my1)//s); ty1=int((R-my0)//s)
n=0; N=17
for tx in range(tx0,tx1+1):
    for ty in range(ty0,ty1+1):
        a,b,c,d=tile_bounds(z,tx,ty)
        gx,gy=np.meshgrid(np.linspace(a,c,N),np.linspace(d,b,N))
        p=inv(np.c_[gx.ravel(),gy.ravel()]).reshape(N,N,2).astype(np.float32)
        if p[...,0].max()<0 or p[...,0].min()>W or p[...,1].max()<0 or p[...,1].min()>H: continue
        # bilinear upsample of coarse map (pixel centers)
        mxp=cv2.resize(p[...,0],(256,256),interpolation=cv2.INTER_LINEAR)*K; myp=cv2.resize(p[...,1],(256,256),interpolation=cv2.INTER_LINEAR)*K
        t=cv2.remap(src,mxp,myp,cv2.INTER_LINEAR,borderMode=cv2.BORDER_CONSTANT,borderValue=(0,0,0,0))
        if t[:,:,3].max()==0: continue
        os.makedirs(f'{OUT}/{z}/{tx}',exist_ok=True)
        Image.fromarray(cv2.cvtColor(t,cv2.COLOR_BGRA2RGBA)).save(f'{OUT}/{z}/{tx}/{ty}.webp',quality=62,method=4)
        n+=1
print('z',z,'tiles',n)
for z in range(ZMAX-1,4,-1):
    cnt=0
    for tx in os.listdir(f'{OUT}/{z+1}'):
        pass
    parents=set()
    for tx in os.listdir(f'{OUT}/{z+1}'):
        for f in os.listdir(f'{OUT}/{z+1}/{tx}'):
            parents.add((int(tx)//2,int(f.split('.')[0])//2))
    for px,py in parents:
        canvas=Image.new('RGBA',(512,512),(0,0,0,0))
        for dx in (0,1):
            for dy in (0,1):
                f=f'{OUT}/{z+1}/{2*px+dx}/{2*py+dy}.webp'
                if os.path.exists(f): canvas.paste(Image.open(f),(dx*256,dy*256))
        os.makedirs(f'{OUT}/{z}/{px}',exist_ok=True)
        canvas.resize((256,256),Image.LANCZOS).save(f'{OUT}/{z}/{px}/{py}.webp',quality=62,method=4); cnt+=1
    print('z',z,cnt)
from georef import Ti
lon0,lat0=Ti.transform(mx0,my0); lon1,lat1=Ti.transform(mx1,my1)
json.dump(dict(bounds=[[lat0,lon0],[lat1,lon1]],zmin=5,zmax=ZMAX),open('/home/claude/site/data/tiles_geo.json','w'))
