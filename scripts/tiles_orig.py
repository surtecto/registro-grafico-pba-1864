import os, math
from PIL import Image
Image.MAX_IMAGE_PIXELS=None
im=Image.open('/home/claude/work/agn.tif'); W,H=im.size; OUT='/home/claude/site/tiles/orig'; ZM=6
cur=im
for z in range(ZM,-1,-1):
    w,h=cur.size; n=0
    for tx in range(math.ceil(w/256)):
        os.makedirs(f'{OUT}/{z}/{tx}',exist_ok=True)
        for ty in range(math.ceil(h/256)):
            t=cur.crop((tx*256,ty*256,min(tx*256+256,w),min(ty*256+256,h)))
            c=Image.new('RGB',(256,256),(20,20,20)); c.paste(t,(0,0)); c.save(f'{OUT}/{z}/{tx}/{ty}.jpg',quality=70); n+=1
    print(z,cur.size,n)
    cur=cur.resize((max(1,w//2),max(1,h//2)),Image.LANCZOS)
print(W,H)
