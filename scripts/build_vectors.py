import sys, json, re, collections, numpy as np
sys.path.insert(0,'/home/claude/scripts')
from georef import build, px2ll, PX, rows as GROWS, LL
from shapely.geometry import Point, MultiPoint, Polygon, mapping, box
from shapely.ops import unary_union, voronoi_diagram, transform as stransform
from pyproj import Transformer, Geod
fwd,inv=build(0)
geod=Geod(ellps='WGS84')
PARTIDOS={1:'General Belgrano',2:'San Isidro',3:'San Fernando',4:'General San Martín',5:'San José de Flores',6:'Barracas al Sur',7:'Quilmes',8:'Ensenada',9:'Magdalena',10:'Rivadavia',11:'Cañuelas',12:'Ranchos',13:'Chascomús',14:'Viedma',15:'Monte',16:'Lobos',17:'Navarro',18:'Las Flores',19:'Pila',20:'Castelli',21:'Tordillo',22:'Ajó',23:'Dolores',24:'Vecino',25:'Arenales',26:'Rauch',27:'Monsalvo',28:'Ayacucho',29:'Tandil',30:'Tuyú',31:'Mar Chiquita',32:'Patagones',36:'Las Conchas',37:'Morón',38:'Matanza',39:'Lomas de Zamora',40:'San Vicente',41:'Merlo',42:'General Las Heras',43:'Villa de Luján',44:'Moreno',45:'Pilar',46:'Exaltación de la Cruz',47:'Zárate',48:'San Antonio de Areco',49:'San Andrés de Giles',50:'Mercedes',51:'Suipacha',52:'Chivilcoy',53:'Chacabuco',54:'Carmen de Areco',55:'Baradero',56:'San Pedro',57:'Arrecifes',58:'Salto',59:'Junín',60:'Rojas',61:'Pergamino',62:'San Nicolás de los Arroyos',63:'Ramallo',64:'Bragado',65:'9 de Julio',66:'25 de Mayo',67:'Saladillo',68:'Tapalquén',69:'Azul',70:'Balcarce (Dr. Antonio)',71:'Lobería',72:'Necochea',73:'Tres Arroyos',74:'Bahía Blanca',75:'Lincoln'}
recs=[]
for i,l in enumerate(open('/home/claude/work/labels.tsv',encoding='utf-8')):
    l=l.strip()
    if not l: continue
    c,t,x,y=l.split('|'); recs.append(dict(id=i+1,cat=c,text=t,px=float(x),py=float(y)))
P=np.array([[r['px'],r['py']] for r in recs]); G=px2ll(P,fwd)
for r,g in zip(recs,G): r['lon'],r['lat']=round(float(g[0]),6),round(float(g[1]),6)
# ---------- envelope (pixel space) ----------
mask_out=[box(6000,3850,8350,6350), box(250,9050,2750,10650)]
land_pts=MultiPoint([(r['px'],r['py']) for r in recs if r['cat'] in('OWN','TOWN','PART','SEC','PNAME')])
from shapely import concave_hull
env=concave_hull(land_pts,ratio=0.12).buffer(110).buffer(-40)
env=env.difference(unary_union(mask_out))
# ---------- partidos: voronoi of seeds ----------
seeds=[(r['px'],r['py'],int(r['text'])) for r in recs if r['cat']=='PART']
vor=voronoi_diagram(MultiPoint([(x,y) for x,y,_ in seeds]),envelope=box(-2000,-2000,12000,14000))
cells=collections.defaultdict(list)
for poly in vor.geoms:
    for x,y,n in seeds:
        if poly.contains(Point(x,y)): cells[n].append(poly); break
partpoly={n:unary_union(ps).intersection(env) for n,ps in cells.items()}
def to_ll(geom):
    def f(x,y,z=None):
        g=px2ll(np.c_[np.atleast_1d(x),np.atleast_1d(y)],fwd); return g[:,0],g[:,1]
    return stransform(f,geom)
def area_ha(g):
    a,_=geod.geometry_area_perimeter(g); return abs(a)/1e4
def which_partido(x,y):
    p=Point(x,y)
    for n,g in partpoly.items():
        if g.contains(p): return n
    # nearest seed fallback
    return min(seeds,key=lambda s:(s[0]-x)**2+(s[1]-y)**2)[2]
for r in recs:
    n=which_partido(r['px'],r['py']); r['partido_n']=n; r['partido']=PARTIDOS.get(n,'')
# ---------- owner classification ----------
COMP=['Sáenz Valiente','Ramos Mejía','Ortiz Basualdo','Pérez Millán','Díaz Vélez','Real de Azúa','Montes de Oca','Peralta Ramos','Vélez Sarsfield','Plaza Montero','Díaz de Bedoya','Ponce de León','Martínez de Hoz','Blanco Casariego','Rodríguez Gaete','Pereyra Lucena','Ramírez de Castro','Fernández de Quiroga']
def classify(t):
    tl=t.lower()
    if tl.startswith('el estado') or tl.startswith('pedido por') or tl.startswith('solicitado'): return 'Estado / tierra pública'
    if 'testamentaria' in tl or 'herederos' in tl or 'sucesores' in tl: return 'Sucesión / herederos'
    if any(k in tl for k in ['sociedad','cía','compañía','hermanos','e hijos','e hijo','y cía','familia','varios']): return 'Sociedad / colectivo'
    return 'Particular'
def surname(t):
    s=re.sub(r'^(El Estado,? (conocido por de|reclamado por( los herederos de Dn\.)?|pedido por)\s*)','',t)
    s=re.sub(r'^(Testamentaria de (Dn\. )?|Herederos de |Solicitado por |Pedido por )','',s)
    s=re.split(r',| hoy | y | e hij',s)[0].strip()
    for c in COMP:
        if s.endswith(c): return c
    w=s.split()
    if not w: return ''
    sn=w[-1]
    if sn.lower() in ('varios','hijos','hijo','cía.','compañía','estado','argentina'): return ''
    return sn
own_cells=[]
owners=[r for r in recs if r['cat']=='OWN']
for r in owners:
    r['tipo']=classify(r['text']); r['apellido']=surname(r['text']) if r['tipo']!='Estado / tierra pública' or 'conocido' in r['text'] else ''
# ---------- approximate parcel zones: voronoi of owner points within partido ----------
ovor=voronoi_diagram(MultiPoint([(r['px'],r['py']) for r in owners]),envelope=box(-2000,-2000,12000,14000))
cellmap={}
for poly in ovor.geoms:
    c=poly.representative_point()
    for r in owners:
        if poly.contains(Point(r['px'],r['py'])): cellmap[r['id']]=poly; break
parcel_feats=[]
for r in owners:
    g=cellmap.get(r['id'])
    if g is None: continue
    g=g.intersection(partpoly.get(r['partido_n'],env)).intersection(Point(r['px'],r['py']).buffer(260))
    if g.is_empty: continue
    gl=to_ll(g); ha=area_ha(gl); r['zona_ha']=round(ha)
    parcel_feats.append(dict(type='Feature',geometry=mapping(gl),properties=dict(id=r['id'],propietario=r['text'],tipo=r['tipo'],apellido=r['apellido'],partido=r['partido'],partido_n=r['partido_n'],zona_aprox_ha=round(ha))))
def fc(feats): return dict(type='FeatureCollection',features=feats)
def ptfeat(r,extra=()):
    pr=dict(id=r['id'],nombre=r['text'],categoria=r['cat'],partido=r['partido'],partido_n=r['partido_n'],px=r['px'],py=r['py'])
    for k in extra: pr[k]=r.get(k)
    return dict(type='Feature',geometry=dict(type='Point',coordinates=[r['lon'],r['lat']]),properties=pr)
def hydtype(t):
    for k in ['Laguna','Arroyo','Río','Cañada','Ensenada','Brazo','Riacho','La Salada','Mar Chiquita','Bahía']:
        if t.startswith(k) or k in t: return {'La Salada':'Laguna','Mar Chiquita':'Laguna','Riacho':'Río','Brazo':'Río','Bahía':'Costa'}.get(k,k)
    return 'Otro'
def geotype(t):
    tl=t.lower()
    for k,v in [('fortín','Fortín'),('fuerte','Fortín'),('estancia','Estancia'),('sierra','Sierra / cerro'),('cerr','Sierra / cerro'),('médano','Médanos'),('punta','Costa'),('cabo','Costa'),('puerto','Puerto'),('isla','Isla'),('banco','Banco fluvial'),('ensenada','Costa'),('chacras','Chacras / ejido'),('paso','Paso'),('vestigios','Ruinas'),('montes','Monte'),('rincón','Paraje')]:
        if k in tl: return v
    return 'Paraje'
for r in recs:
    if r['cat']=='HYD': r['tipo']=hydtype(r['text'])
    if r['cat']=='GEO': r['tipo']=geotype(r['text'])
D='/home/claude/site/data/'
json.dump(fc([ptfeat(r,('tipo','apellido','zona_ha')) for r in owners]),open(D+'propietarios.geojson','w'),ensure_ascii=False)
json.dump(fc(parcel_feats),open(D+'parcelas_aprox.geojson','w'),ensure_ascii=False)
json.dump(fc([ptfeat(r) for r in recs if r['cat']=='TOWN']),open(D+'poblados.geojson','w'),ensure_ascii=False)
json.dump(fc([ptfeat(r,('tipo',)) for r in recs if r['cat']=='HYD']),open(D+'hidrografia.geojson','w'),ensure_ascii=False)
json.dump(fc([ptfeat(r,('tipo',)) for r in recs if r['cat'] in('GEO','OTHER')]),open(D+'geografia.geojson','w'),ensure_ascii=False)
json.dump(fc([ptfeat(r) for r in recs if r['cat'] in('SEC','PNAME','PART')]),open(D+'rotulos.geojson','w'),ensure_ascii=False)
pf=[]
stats_part={}
for n,g in sorted(partpoly.items()):
    if g.is_empty: continue
    gl=to_ll(g); ha=area_ha(gl)
    os_=[r for r in owners if r['partido_n']==n]
    stats_part[n]=dict(n=n,nombre=PARTIDOS.get(n,''),area_aprox_km2=round(ha/100),propietarios=len(os_),
        estado=sum(r['tipo'].startswith('Estado') for r in os_),poblados=[r['text'] for r in recs if r['cat']=='TOWN' and r['partido_n']==n],
        hidro=sum(1 for r in recs if r['cat']=='HYD' and r['partido_n']==n))
    pf.append(dict(type='Feature',geometry=mapping(gl),properties=dict(n=n,nombre=PARTIDOS.get(n,''),area_aprox_km2=round(ha/100),propietarios=len(os_),limite='aproximado (Voronoi de rótulos de partido)')))
json.dump(fc(pf),open(D+'partidos_aprox.geojson','w'),ensure_ascii=False)
# GCPs
fx=fwd(PX); from georef import Ti, MX
gc=[]
for r,p,m,mm in zip(GROWS,PX,fx,MX):
    lon,lat=Ti.transform(m[0],m[1])
    gc.append(dict(type='Feature',geometry=dict(type='Point',coordinates=[float(r['lon']),float(r['lat'])]),properties=dict(nombre=r['name'],px=float(p[0]),py=float(p[1]))))
json.dump(fc(gc),open(D+'puntos_control.geojson','w'),ensure_ascii=False)
# ---------- stats ----------
C=collections.Counter
st=dict(
 totales=dict(propietarios=len(owners),propietarios_distintos=len(set(r['text'] for r in owners)),poblados=sum(r['cat']=='TOWN' for r in recs),
   hidrografia=sum(r['cat']=='HYD' for r in recs),geografia=sum(r['cat'] in('GEO','OTHER') for r in recs),partidos_localizados=len(partpoly),partidos_lista=72,puntos_control=len(GROWS)),
 tipos_propietario=C(r['tipo'] for r in owners).most_common(),
 apellidos_top=[(a,c) for a,c in C(r['apellido'] for r in owners if r['apellido']).most_common(40)],
 propietarios_repetidos=[(a,c) for a,c in C(r['text'] for r in owners if not r['text'].startswith('El Estado')).most_common(25) if c>1],
 hidro_tipos=C(r['tipo'] for r in recs if r['cat']=='HYD').most_common(),
 geo_tipos=C(r['tipo'] for r in recs if r['cat']=='GEO').most_common(),
 partidos=list(stats_part.values()),
 lista_partidos=[dict(n=k,nombre=v,localizado=k in partpoly) for k,v in PARTIDOS.items()],
 apellidos_zona=[(a,round(sum(r.get('zona_ha',0) for r in owners if r['apellido']==a)/100)) for a,_ in C(r['apellido'] for r in owners if r['apellido']).most_common(200)],
)
st['apellidos_zona']=sorted(st['apellidos_zona'],key=lambda t:-t[1])[:25]
json.dump(st,open(D+'estadisticas.json','w'),ensure_ascii=False,indent=1)
print(st['totales']); print(st['tipos_propietario']); print(st['apellidos_top'][:15])
# env for raster mask (pixel space) 
json.dump(dict(mask_out=[list(m.exterior.coords) for m in mask_out]),open('/home/claude/work/mask.json','w'))
