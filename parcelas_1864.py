#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
parcelas_1864.py — Capa de parcelas del Registro Gráfico de la Provincia de Buenos Aires (1864)

Arq. Juan Patricio Addesso · Surtectura · Doctorado en Geografía (USAL)

Idea del método
---------------
No "inunda" desde un clic (lo que hacía que la varita se escapara por las líneas
débiles y dejara muescas en los textos). Trabaja sobre las LÍNEAS del plano:

  1. Realce de trazos finos oscuros (black-hat) sobre la lámina.
  2. Filtro orientado largo (35 px, 90 direcciones): las líneas de deslinde son
     rectas y largas; el texto responde parejo en todas las direcciones y se descarta.
  3. Orientaciones dominantes locales: un trazo débil se acepta sólo si es paralelo
     a los deslindes fuertes de su entorno (así aparecen las líneas tenues entre
     suertes y desaparecen las letras).
  4. Ríos y lagunas (rayado denso) se toman por su contorno.
  5. Cada línea que termina en el aire se prolonga en recta (hasta 50 px) para
     cerrar los cortes. Lo que sobra queda como "colgante" y no forma parcela.
  6. Poligonizado topológico: las parcelas vecinas comparten EXACTAMENTE el mismo
     borde (sin rendijas ni solapes).
  7. Limpieza de astillas y caras de texto; tus parcelas manuales quedan fijas y
     sus bordes se incorporan a la red.
  8. Propietario por los rótulos transcriptos que caen adentro; tus nombres mandan.

Flujo sugerido en QGIS
----------------------
  a) python parcelas_1864.py lamina   --teselas tiles/orig/7 --salida SIG
  b) python parcelas_1864.py lineas   --teselas tiles/orig/7 --salida SIG [--manual respaldo.json]
  c) Abrir en QGIS SIG/lamina_1864.jpg + SIG/lineas_deslinde.geojson y corregir LÍNEAS
     (borrar las falsas, dibujar las que faltan). Es mucho más rápido que corregir polígonos.
  d) python parcelas_1864.py parcelas --salida SIG [--manual respaldo.json]
     → SIG/parcelas.geojson (coordenadas de lámina) y SIG/parcelas_wgs84.geojson (georreferenciada),
       SIG/partidos_wgs84.geojson (disolución por partido asignado).
  Repetir c) y d) las veces que haga falta.

Coordenadas "de lámina": x = columna, y = −fila del mosaico de 17008 × 21670 px
(el signo negativo hace que la lámina y las capas coincidan en QGIS con el archivo .jgw).

Requisitos: Python 3.9+, numpy, opencv-python-headless, shapely>=2, scikit-image.
  pip install numpy opencv-python-headless shapely scikit-image
"""
import argparse, json, math, os, sys, time, csv
import numpy as np
import cv2
from skimage.morphology import skeletonize
from shapely.geometry import LineString, Polygon, MultiPolygon, Point, box, mapping, shape
from shapely.ops import unary_union, polygonize, linemerge
from shapely import STRtree

W_LAM, H_LAM = 17008, 21670          # mosaico de la lámina (px)
AGN_SCALE = 2.0                       # rótulos / puntos de control en px de la copia AGN (8504 × 10835)

PARAM = dict(
    k=7,            # black-hat: ancho máximo de trazo (px)
    L=35,           # largo del filtro orientado (px)
    nang=90,        # direcciones (2°)
    win=151,        # ventana para orientaciones dominantes (px)
    frac=0.35,      # una orientación es "dominante" si tiene ≥35 % del pico local
    t=6,            # umbral para trazos débiles paralelos a los dominantes
    t_hi=22,        # umbral para trazos fuertes (cualquier orientación)
    text=30,        # se descartan manchas de menos de 30 px de extensión
    spur=20,        # colgantes más cortos que esto se eliminan
    eps=1.0,        # simplificación de líneas (px)
    gap=50,         # prolongación máxima de una línea cortada (px)
    gapk=0.5,       # ... y como mucho la mitad de su propio largo
    extmin=25,      # sólo se prolongan líneas de al menos 25 px
    amin=150,       # caras menores (px²) se funden con la vecina
    wmin=3.0,       # caras de menos de 3 px de ancho medio se funden
    wf=0.76,
    norm=64,        # percentil 98 del realce al que se lleva cada ventana        # agua: oscuridad tras desenfoque fuerte (fracción de la mediana)
)

# ───────────────────────── lectura de la lámina ─────────────────────────
class Lamina:
    """Acceso por ventanas a la lámina: un archivo único o la carpeta de teselas z7 del repo."""
    def __init__(self, imagen=None, teselas=None):
        self.img = None; self.dir = teselas
        if imagen:
            cv2.setNumThreads(4)
            os.environ.setdefault('OPENCV_IO_MAX_IMAGE_PIXELS', str(2 ** 40))
            self.img = cv2.imread(imagen, cv2.IMREAD_COLOR)
            if self.img is None: sys.exit(f'No pude leer {imagen}')
            self.H, self.W = self.img.shape[:2]
        elif teselas:
            self.ts = 256; self.W, self.H = W_LAM, H_LAM
        else:
            sys.exit('Indicá --lamina o --teselas')
        self._cache = {}

    def ventana(self, x0, y0, x1, y1):
        x0, y0 = max(0, x0), max(0, y0); x1, y1 = min(self.W, x1), min(self.H, y1)
        if self.img is not None:
            return self.img[y0:y1, x0:x1].copy()
        out = np.full((y1 - y0, x1 - x0, 3), 255, np.uint8); T = self.ts
        for tx in range(x0 // T, (x1 - 1) // T + 1):
            for ty in range(y0 // T, (y1 - 1) // T + 1):
                p = os.path.join(self.dir, str(tx), f'{ty}.jpg')
                if not os.path.exists(p): continue
                t = cv2.imread(p)
                X, Y = tx * T, ty * T
                ax, ay = max(x0, X), max(y0, Y); bx, by = min(x1, X + t.shape[1]), min(y1, Y + t.shape[0])
                if bx > ax and by > ay:
                    out[ay - y0:by - y0, ax - x0:bx - x0] = t[ay - Y:by - Y, ax - X:bx - X]
        return out

# ───────────────────────── detección de líneas ─────────────────────────
NB = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]

def remove_small(mask, min_extent):
    n, lab, st, _ = cv2.connectedComponentsWithStats(mask, connectivity=8)
    keep = np.maximum(st[:, 2], st[:, 3]) >= min_extent; keep[0] = False
    return keep[lab].astype(np.uint8)

def _kern(L, th):
    K = np.zeros((L, L), np.float32); c = L // 2; dx, dy = math.cos(th), math.sin(th)
    cv2.line(K, (int(round(c - dx * c)), int(round(c - dy * c))), (int(round(c + dx * c)), int(round(c + dy * c))), 1, 1)
    return K / K.sum()

def line_response(bh, L, nang):
    bh = bh.astype(np.float32)
    R = np.full(bh.shape, -1, np.float32); A = np.zeros(bh.shape, np.int16); S = np.zeros_like(bh)
    for k in range(nang):
        r = cv2.filter2D(bh, -1, _kern(L, math.pi * k / nang))
        S += r; m = r > R; R[m] = r[m]; A[m] = k
    return R - S / nang, A

def nms(D, A, nang):
    th = np.pi * A.astype(np.float32) / nang
    nx, ny = (-np.sin(th)).astype(np.float32), np.cos(th).astype(np.float32)
    H, W = D.shape; yy, xx = np.mgrid[0:H, 0:W].astype(np.float32)
    a = cv2.remap(D, xx + nx, yy + ny, cv2.INTER_LINEAR)
    b = cv2.remap(D, xx - nx, yy - ny, cv2.INTER_LINEAR)
    return (D >= a) & (D >= b)

def dominant_ok(D, A, nang, t_hi, win, frac, f=4):
    """¿La orientación de cada píxel es una de las dominantes en su entorno? (histograma a 1/4 de resolución)"""
    strong = np.where(D > t_hi, D, 0).astype(np.float32)
    H, W = D.shape; h, w = (H + f - 1) // f, (W + f - 1) // f
    Hs = np.empty((nang, h, w), np.float32); kw = max(3, win // f)
    for k in range(nang):
        s = cv2.resize(strong * (A == k), (w, h), interpolation=cv2.INTER_AREA)
        Hs[k] = cv2.blur(s, (kw, kw))
    Hs = Hs + np.roll(Hs, 1, 0) + np.roll(Hs, -1, 0)
    mx = Hs.max(0) + 1e-9
    yy, xx = np.mgrid[0:H, 0:W]; yy //= f; xx //= f
    val = Hs[A.astype(np.int64), yy, xx]
    return val >= frac * mx[yy, xx]

def water_mask(bgr, f, med):
    g = cv2.cvtColor(bgr, cv2.COLOR_BGR2GRAY)
    b = cv2.GaussianBlur(g, (0, 0), 5)
    m = (b < med * f).astype(np.uint8)
    m = cv2.morphologyEx(m, cv2.MORPH_CLOSE, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (7, 7)))
    return remove_small(m, 200)

def line_map(bgr, P):
    g = cv2.cvtColor(bgr, cv2.COLOR_BGR2GRAY)
    bh = cv2.morphologyEx(g, cv2.MORPH_BLACKHAT, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (P['k'], P['k'])))
    if P.get('norm'):   # contraste parejo: sectores grises y desvaídos rinden como los nítidos
        ref = np.percentile(bh[g < 250], 98) if (g < 250).any() else 64
        bh = np.clip(bh.astype(np.float32) * float(np.clip(P['norm'] / max(ref, 1), 0.7, 2.5)), 0, 255)
    D, A = line_response(bh, P['L'], P['nang'])
    ok = dominant_ok(D, A, P['nang'], P['t_hi'], P['win'], P['frac'])
    m = nms(D, A, P['nang']) & ((D > P['t_hi']) | ((D > P['t']) & ok))
    m = remove_small(m.astype(np.uint8), P['text'])
    return skeletonize(m > 0).astype(np.uint8)

def skel_paths(sk):
    """Recorre el esqueleto y lo convierte en polilíneas entre extremos y nudos."""
    H, W = sk.shape
    P2 = np.pad(sk, 1)
    deg = sum(np.roll(np.roll(P2, dy, 0), dx, 1) for dy, dx in NB)[1:-1, 1:-1] * sk
    node = (deg != 2) & (sk > 0)
    S = sk.astype(bool); seen = np.zeros_like(S); used = set(); paths = []
    def nbrs(y, x):
        for dy, dx in NB:
            yy, xx = y + dy, x + dx
            if 0 <= yy < H and 0 <= xx < W and S[yy, xx]: yield yy, xx
    for y, x in zip(*np.nonzero(node)):
        for ny, nx in nbrs(y, x):
            key = ((y, x), (ny, nx)) if (y, x) < (ny, nx) else ((ny, nx), (y, x))
            if key in used: continue
            used.add(key); path = [(x, y)]; py, px, cy, cx = y, x, ny, nx
            while True:
                path.append((cx, cy)); seen[cy, cx] = True
                if node[cy, cx]: break
                nxt = None
                for t in nbrs(cy, cx):
                    if t != (py, px): nxt = t; break
                if nxt is None: break
                k2 = ((cy, cx), nxt) if (cy, cx) < nxt else (nxt, (cy, cx))
                if k2 in used: break
                used.add(k2); py, px, cy, cx = cy, cx, nxt[0], nxt[1]
            paths.append(path)
    rest = S & ~seen & ~node
    if rest.any():
        n, lab = cv2.connectedComponents(rest.astype(np.uint8), connectivity=8)
        for k in range(1, n):
            c, _ = cv2.findContours((lab == k).astype(np.uint8), cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
            if c and len(c[0]) > 8:
                pts = [tuple(p[0]) for p in c[0]]; paths.append(pts + [pts[0]])
    return paths, deg

def extend_end(ls, G, at_start):
    L = ls.length; t = min(10.0, L)
    if at_start: p0 = np.array(ls.coords[0]); p1 = np.array(ls.interpolate(t).coords[0])
    else: p0 = np.array(ls.coords[-1]); p1 = np.array(ls.interpolate(L - t).coords[0])
    d = p0 - p1; n = np.linalg.norm(d)
    if n == 0: return ls
    e = tuple(p0 + d / n * G); c = list(ls.coords)
    return LineString([e] + c) if at_start else LineString(c + [e])

def lines_from(bgr, P, med):
    """Líneas de deslinde de una ventana, en coordenadas de la ventana."""
    sk = line_map(bgr, P)
    paths, deg = skel_paths(sk)
    out = []
    wm = water_mask(bgr, P['wf'], med)
    cs, _ = cv2.findContours(wm, cv2.RETR_LIST, cv2.CHAIN_APPROX_NONE)
    for c in cs:
        if len(c) < 8: continue
        q = [tuple(map(float, v[0])) for v in c]
        out.append(('agua', LineString(q + [q[0]]).simplify(P['eps'])))
    for p in paths:
        if len(p) < 2: continue
        ls = LineString(p); a, b = p[0], p[-1]
        fa = deg[a[1], a[0]] <= 1; fb = deg[b[1], b[0]] <= 1
        if (fa or fb) and ls.length < P['spur']: continue
        ls = ls.simplify(P['eps'])
        if ls.length >= P['extmin']:
            G = min(P['gap'], P['gapk'] * ls.length)
            if fa: ls = extend_end(ls, G, True)
            if fb: ls = extend_end(ls, G, False)
        out.append(('auto', ls))
    return out, wm

# ───────────────────────── poligonizado y limpieza ─────────────────────────
def merge_small(polys, amin, wmin, keep=()):
    polys = [p.buffer(0) for p in polys]; keep = set(keep)
    changed = True; passes = 0
    while changed and passes < 4:
        changed = False; passes += 1; tree = STRtree(polys); dead = set()
        for i in sorted(range(len(polys)), key=lambda i: polys[i].area):
            if i in dead or i in keep: continue
            p = polys[i]
            if p.area >= amin and 2 * p.area / max(p.length, 1e-9) >= wmin: continue
            best, bl = None, 0
            for j in tree.query(p):
                if j == i or j in dead or j in keep: continue
                l = p.boundary.intersection(polys[j].boundary).length
                if l > bl: bl, best = l, j
            if best is None: continue
            polys[best] = unary_union([polys[best], p]).buffer(0); dead.add(i); changed = True
        idx = [k for k in range(len(polys)) if k not in dead]
        keep = {idx.index(k) for k in keep if k in idx}
        polys = [polys[k] for k in idx]
    return polys, keep

# ───────────────────────── georreferenciado ─────────────────────────
class Georef:
    """Polinomio de 2.º grado + thin-plate spline de residuos, desde los puntos de control."""
    def __init__(self, csv_path):
        rows = list(csv.DictReader(open(csv_path, encoding='utf-8')))
        px = np.array([[float(r['px']) * AGN_SCALE, float(r['py']) * AGN_SCALE] for r in rows])
        ll = np.array([[float(r['lon']), float(r['lat'])] for r in rows])
        self.src = px; self.c = self._mean_std(px)
        D = self._D(px); self.coef, *_ = np.linalg.lstsq(D, ll, rcond=None)
        res = ll - D @ self.coef
        n = len(px); K = self._U(px, px) + np.eye(n) * 1e-9
        Pm = np.c_[np.ones(n), self._n(px)]
        A = np.zeros((n + 3, n + 3)); A[:n, :n] = K; A[:n, n:] = Pm; A[n:, :n] = Pm.T
        b = np.zeros((n + 3, 2)); b[:n] = res
        self.w = np.linalg.solve(A, b)
    def _mean_std(self, p): return p.mean(0), p.std(0)
    def _n(self, p): return (p - self.c[0]) / self.c[1]
    def _D(self, p):
        x, y = self._n(p).T; return np.c_[np.ones_like(x), x, y, x * x, x * y, y * y]
    def _U(self, a, b):
        r = np.linalg.norm(self._n(a)[:, None, :] - self._n(b)[None, :, :], axis=2)
        with np.errstate(divide='ignore', invalid='ignore'): u = r * r * np.log(r)
        return np.nan_to_num(u)
    def __call__(self, pts):
        pts = np.atleast_2d(np.asarray(pts, float)); n = len(self.src)
        out = self._D(pts) @ self.coef
        out += self._U(pts, self.src) @ self.w[:n] + np.c_[np.ones(len(pts)), self._n(pts)] @ self.w[n:]
        return out

def geom_to_wgs(g, gr):
    def ring(r):
        return [tuple(v) for v in gr(np.array(r.coords))]
    if g.geom_type == 'Polygon':
        return Polygon(ring(g.exterior), [ring(i) for i in g.interiors])
    return MultiPolygon([geom_to_wgs(p, gr) for p in g.geoms])

def area_ha(gw):
    lat0 = gw.centroid.y; kx = 111320 * math.cos(math.radians(lat0)); ky = 110540
    from shapely.affinity import scale
    return scale(gw, kx, ky, origin=(0, 0)).area / 1e4

# ───────────────────────── E/S ─────────────────────────
def to_plano(g):            # (x, y_fila) → (x, −y)  para QGIS con el .jgw
    from shapely.affinity import scale
    return scale(g, 1, -1, origin=(0, 0))

def write_fc(path, feats):
    with open(path, 'w', encoding='utf-8') as f:
        json.dump({'type': 'FeatureCollection', 'features': feats}, f, ensure_ascii=False)

def read_manual(path):
    """Respaldo del editor («Copiar respaldo»): {id: {g:[x,y,...], o:propietario, p:partido}}."""
    if not path: return []
    d = json.load(open(path, encoding='utf-8'))
    if 'e' in d and isinstance(d['e'], dict): d = d['e']
    out = []
    for k, v in d.items():
        g = np.array(v['g'], float).reshape(-1, 2)
        if len(g) >= 3:
            out.append(dict(id=f'm{k}', geom=Polygon(g).buffer(0), o=v.get('o', ''), p=int(v.get('p') or 0)))
    return out

def read_rotulos(path):
    if not path or not os.path.exists(path): return []
    if path.endswith('.json'):
        return [dict(c=r['c'], t=r['t'], x=float(r['x']), y=float(r['y'])) for r in json.load(open(path, encoding='utf-8'))]
    out = []   # transcripcion_rotulos.tsv (coordenadas AGN)
    for line in open(path, encoding='utf-8'):
        a = line.rstrip('\n').split('|')
        if len(a) == 4: out.append(dict(c=a[0], t=a[1], x=float(a[2]) * AGN_SCALE, y=float(a[3]) * AGN_SCALE))
    return out

PARTIDOS = {1:'General Belgrano',2:'San Isidro',3:'San Fernando',4:'General San Martín',5:'San José de Flores',6:'Barracas al Sur',7:'Quilmes',8:'Ensenada',9:'Magdalena',10:'Rivadavia',11:'Cañuelas',12:'Ranchos',13:'Chascomús',14:'Viedma',15:'Monte',16:'Lobos',17:'Navarro',18:'Las Flores',19:'Pila',20:'Castelli',21:'Tordillo',22:'Ajó',23:'Dolores',24:'Vecino',25:'Arenales',26:'Rauch',27:'Monsalvo',28:'Ayacucho',29:'Tandil',30:'Tuyú',31:'Mar Chiquita',32:'Patagones',36:'Las Conchas',37:'Morón',38:'Matanza',39:'Lomas de Zamora',40:'San Vicente',41:'Merlo',42:'General Las Heras',43:'Villa de Luján',44:'Moreno',45:'Pilar',46:'Exaltación de la Cruz',47:'Zárate',48:'San Antonio de Areco',49:'San Andrés de Giles',50:'Mercedes',51:'Suipacha',52:'Chivilcoy',53:'Chacabuco',54:'Carmen de Areco',55:'Baradero',56:'San Pedro',57:'Arrecifes',58:'Salto',59:'Junín',60:'Rojas',61:'Pergamino',62:'San Nicolás de los Arroyos',63:'Ramallo',64:'Bragado',65:'9 de Julio',66:'25 de Mayo',67:'Saladillo',68:'Tapalquén',69:'Azul',70:'Balcarce',71:'Lobería',72:'Necochea',73:'Tres Arroyos',74:'Bahía Blanca',75:'Lincoln'}

# ───────────────────────── comandos ─────────────────────────
def _tile_job(args):
    lam_args, x0, y0, x1, y1, core, P, med = args
    lam = Lamina(**lam_args)
    img = lam.ventana(x0, y0, x1, y1)
    g = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    if (g < 250).mean() < 0.05: return [], None
    lines, wm = lines_from(img, P, med)
    cb = box(core[0] - x0, core[1] - y0, core[2] - x0, core[3] - y0)
    out = []
    for src, ls in lines:
        c = ls.intersection(cb)
        if c.is_empty: continue
        for part in getattr(c, 'geoms', [c]):
            if part.geom_type == 'LineString' and part.length > 0.5:
                a = np.array(part.coords) + [x0, y0]; out.append((src, a.round(2).tolist()))
    sub = wm[core[1] - y0:core[3] - y0, core[0] - x0:core[2] - x0]
    small = cv2.resize(sub, (max(1, sub.shape[1] // 4), max(1, sub.shape[0] // 4)), interpolation=cv2.INTER_AREA)
    return out, (core, small)

def cmd_lineas(a):
    P = dict(PARAM)
    cal = os.path.join(a.salida, 'calibracion.json')
    if os.path.exists(cal): P.update(json.load(open(cal))); print('Usando', cal, json.load(open(cal)))
    P.update(json.loads(a.param) if a.param else {})
    lam_args = dict(imagen=a.lamina, teselas=a.teselas); lam = Lamina(**lam_args)
    zx0, zy0, zx1, zy1 = (map(int, a.zona.split(','))) if a.zona else (0, 0, lam.W, lam.H)
    # mediana global del papel (muestra)
    smp = [lam.ventana(x, y, x + 256, y + 256) for x in range(zx0, zx1, max(256, (zx1 - zx0) // 8)) for y in range(zy0, zy1, max(256, (zy1 - zy0) // 8))]
    gs = np.concatenate([cv2.cvtColor(s, cv2.COLOR_BGR2GRAY).ravel() for s in smp]); gs = gs[gs < 250]
    med = float(np.median(gs)) if len(gs) else 135.0
    T, O = a.tesela, 96
    jobs = []
    for y in range(zy0, zy1, T):
        for x in range(zx0, zx1, T):
            core = (x, y, min(x + T, zx1), min(y + T, zy1))
            jobs.append((lam_args, x - O, y - O, core[2] + O, core[3] + O, core, P, med))
    print(f'{len(jobs)} ventanas de {T}px · mediana del papel {med:.0f}', flush=True)
    t0 = time.time(); res = []
    if a.procesos > 1:
        from multiprocessing import Pool
        with Pool(a.procesos, maxtasksperchild=4) as pool:
            for i, r in enumerate(pool.imap_unordered(_tile_job, jobs)):
                res.append(r); print(f'\r  {i + 1}/{len(jobs)}  {time.time() - t0:.0f}s', end='', flush=True)
    else:
        for i, j in enumerate(jobs):
            res.append(_tile_job(j)); print(f'\r  {i + 1}/{len(jobs)}  {time.time() - t0:.0f}s', end='', flush=True)
    print()
    os.makedirs(a.salida, exist_ok=True)
    feats = []
    for out, _ in res:
        for src, c in out:
            feats.append({'type': 'Feature', 'properties': {'fuente': src},
                          'geometry': {'type': 'LineString', 'coordinates': [[x, -y] for x, y in c]}})
    for m in read_manual(a.manual):
        c = [[x, -y] for x, y in m['geom'].exterior.coords]
        feats.append({'type': 'Feature', 'properties': {'fuente': 'manual'}, 'geometry': {'type': 'LineString', 'coordinates': c}})
    write_fc(os.path.join(a.salida, 'lineas_deslinde.geojson'), feats)
    # máscara de agua a 1/4 de resolución
    WM = np.zeros((lam.H // 4 + 1, lam.W // 4 + 1), np.uint8)
    for _, w in res:
        if w is None: continue
        core, small = w; X, Y = core[0] // 4, core[1] // 4
        h, wd = small.shape; WM[Y:Y + h, X:X + wd] = np.maximum(WM[Y:Y + h, X:X + wd], small * 255)
    cv2.imwrite(os.path.join(a.salida, 'agua_x4.png'), WM)
    json.dump(dict(P, mediana_papel=med, zona=[zx0, zy0, zx1, zy1]), open(os.path.join(a.salida, 'parametros.json'), 'w'), indent=1)
    print(f'{len(feats)} líneas → {a.salida}/lineas_deslinde.geojson')

def cmd_parcelas(a):
    P = dict(PARAM)
    pj = os.path.join(a.salida, 'parametros.json')
    if os.path.exists(pj): P.update({k: v for k, v in json.load(open(pj)).items() if k in PARAM})
    lp = a.lineas or os.path.join(a.salida, 'lineas_deslinde.geojson')
    L = json.load(open(lp, encoding='utf-8'))
    lines = []
    for f in L['features']:
        g = f.get('geometry')
        if not g: continue
        if (f.get('properties') or {}).get('fuente') == 'manual': continue   # se vuelven a tomar del respaldo
        parts = [g['coordinates']] if g['type'] == 'LineString' else g['coordinates'] if g['type'] == 'MultiLineString' else []
        for c in parts:
            if len(c) >= 2: lines.append(LineString([(x, -y) for x, y in c]))
    manual = read_manual(a.manual)
    for m in manual: lines.append(m['geom'].exterior)
    zona = json.load(open(pj)).get('zona') if os.path.exists(pj) else [0, 0, W_LAM, H_LAM]
    frame = box(*zona)
    print(f'{len(lines)} líneas · nodado y poligonizado…', flush=True); t0 = time.time()
    U = unary_union(lines + [frame.exterior])
    polys = [p for p in polygonize(U) if p.area > 0]
    print(f'  {len(polys)} caras ({time.time() - t0:.0f}s)', flush=True)
    # tus parcelas manuales mandan: se quitan las caras que caen adentro
    if manual:
        MU = unary_union([m['geom'] for m in manual]); keep = []
        for p in polys:
            ov = p.intersection(MU).area
            if ov > 0.5 * p.area: continue
            if ov > 0: p = p.difference(MU).buffer(0)
            if not p.is_empty: keep.extend(getattr(p, 'geoms', [p]))
        polys = [p for p in keep if p.geom_type == 'Polygon']
    # caras enormes = recintos que no cerraron (o el exterior del plano): se descartan antes de limpiar
    big = a.max_area
    nb = sum(1 for p in polys if p.area > big)
    polys = [p for p in polys if p.area <= big]
    print(f'  {nb} recintos abiertos descartados (> {big:.0f} px²)', flush=True)
    polys, _ = merge_small(polys, P['amin'], P['wmin'])
    print(f'  {len(polys)} caras tras la limpieza ({time.time() - t0:.0f}s)', flush=True)
    # agua
    WMp = os.path.join(a.salida, 'agua_x4.png'); WM = cv2.imread(WMp, 0) if os.path.exists(WMp) else None
    def agua_frac(p):
        if WM is None: return 0.0
        x0, y0, x1, y1 = [int(v // 4) for v in p.bounds]
        sub = WM[max(0, y0):y1 + 1, max(0, x0):x1 + 1]
        if sub.size == 0: return 0.0
        msk = np.zeros_like(sub)
        cv2.fillPoly(msk, [((np.array(p.exterior.coords) / 4) - [max(0, x0), max(0, y0)]).astype(np.int32)], 1)
        n = msk.sum(); return float((sub[msk > 0] > 127).sum() / n) if n else 0.0
    # rótulos
    rot = read_rotulos(a.rotulos)
    own = [r for r in rot if r['c'] == 'OWN']; part = [r for r in rot if r['c'] == 'PART']
    tree = STRtree(polys); hits = {i: [] for i in range(len(polys))}
    for r in own:
        for i in tree.query(Point(r['x'], r['y'])):
            if polys[i].contains(Point(r['x'], r['y'])): hits[i].append(r['t']); break
    gr = Georef(a.gcp) if a.gcp and os.path.exists(a.gcp) else None
    feats, featw = [], []
    def add(pid, g, origen, prop, nrot, todos, pn, calidad):
        props = dict(id=pid, propietario=prop, origen=origen, n_rotulos=nrot, rotulos=' / '.join(todos),
                     partido_n=pn or None, partido=PARTIDOS.get(pn, '') if pn else '', calidad=calidad, area_px=round(g.area))
        feats.append({'type': 'Feature', 'properties': props, 'geometry': mapping(to_plano(g))})
        if gr is not None:
            gw = geom_to_wgs(g, gr); props2 = dict(props, sup_ha=round(area_ha(gw)))
            featw.append({'type': 'Feature', 'properties': props2, 'geometry': mapping(gw)})
    # lo que cargaste en QGIS (propietario, partido) en la corrida anterior se conserva
    prev = []; pp = os.path.join(a.salida, 'parcelas.geojson')
    if os.path.exists(pp):
        for f in json.load(open(pp, encoding='utf-8'))['features']:
            pr = f['properties']
            if pr.get('origen') == 'manual': continue
            if pr.get('propietario') or pr.get('partido_n'):
                prev.append((to_plano(shape(f['geometry'])).buffer(0), pr.get('propietario') or '', int(pr.get('partido_n') or 0)))
        os.replace(pp, os.path.join(a.salida, 'parcelas_anterior.geojson'))
    ptree = STRtree([g for g, _, _ in prev]) if prev else None
    def heredado(p):
        if ptree is None: return '', 0
        q = p.representative_point()
        for j in ptree.query(q):
            if prev[j][0].contains(q): return prev[j][1], prev[j][2]
        return '', 0
    n = 0
    for m in manual:
        add(m['id'], m['geom'], 'manual', m['o'], 0, [], m['p'], 'manual'); n += 1
    for i, p in enumerate(polys):
        t = hits[i]; af = agua_frac(p)
        if af > 0.5: cal = 'agua'
        elif len(t) == 0: cal = 'sin_rotulo'
        elif len(t) == 1: cal = 'ok'
        else: cal = 'varios_rotulos'      # probablemente parcelas unidas: falta una línea
        ho, hp = heredado(p)
        prop = ho or (t[0] if len(t) == 1 else '')
        if ho and cal != 'agua': cal = 'ok'
        add(f'a{i + 1}', p, 'auto', prop, len(t), t, hp, cal); n += 1
    write_fc(os.path.join(a.salida, 'parcelas.geojson'), feats)
    if featw:
        write_fc(os.path.join(a.salida, 'parcelas_wgs84.geojson'), featw)
        # partidos = disolución de las parcelas con partido asignado
        byp = {}
        for f in featw:
            pn = f['properties']['partido_n']
            if pn: byp.setdefault(pn, []).append(shape(f['geometry']))
        write_fc(os.path.join(a.salida, 'partidos_wgs84.geojson'), [
            {'type': 'Feature', 'properties': {'partido_n': k, 'partido': PARTIDOS.get(k, ''), 'parcelas': len(v)},
             'geometry': mapping(unary_union([g.buffer(1e-6) for g in v]).buffer(-1e-6))} for k, v in sorted(byp.items())])
    from collections import Counter
    c = Counter(f['properties']['calidad'] for f in feats)
    print(f'{n} parcelas → {a.salida}/parcelas.geojson' + (' y parcelas_wgs84.geojson' if featw else ''))
    print('  ' + ' · '.join(f'{k}: {v}' for k, v in c.most_common()))

def cmd_lamina(a):
    lam = Lamina(imagen=a.lamina, teselas=a.teselas)
    os.makedirs(a.salida, exist_ok=True)
    out = os.path.join(a.salida, 'lamina_1864.jpg')
    big = np.full((lam.H, lam.W, 3), 255, np.uint8)
    for y in range(0, lam.H, 2048):
        big[y:y + 2048] = lam.ventana(0, y, lam.W, min(lam.H, y + 2048)); print(f'\r  fila {y}', end='', flush=True)
    cv2.imwrite(out, big, [cv2.IMWRITE_JPEG_QUALITY, 90])
    open(out[:-4] + '.jgw', 'w').write('1\n0\n0\n-1\n0.5\n-0.5\n')
    print(f'\n{out} (+ .jgw: 1 px = 1 unidad, y hacia abajo negativa)')

def cmd_calibrar(a):
    """Prueba combinaciones de parámetros en la zona de tus parcelas manuales y elige la que mejor las reproduce."""
    manual = read_manual(a.manual)
    if not manual: sys.exit('Hace falta --manual con al menos unas parcelas')
    lam = Lamina(imagen=a.lamina, teselas=a.teselas)
    MU = unary_union([m['geom'] for m in manual]); x0, y0, x1, y1 = [int(v) for v in MU.buffer(250).bounds]
    img = lam.ventana(x0, y0, x1, y1); g = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY); med = float(np.median(g[g < 250]))
    from shapely.affinity import translate
    M = [translate(m['geom'], -x0, -y0) for m in manual]
    grid = [dict(t=t, t_hi=th, gap=gp) for t in (5, 6, 8) for th in (18, 22, 26) for gp in (35, 50)]
    best = None
    for q in grid:
        P = dict(PARAM, **q)
        lines, _ = lines_from(img, P, med)
        polys = [p for p in polygonize(unary_union([l for _, l in lines] + [box(0, 0, x1 - x0, y1 - y0).exterior])) if p.area > 0]
        polys, _ = merge_small(polys, P['amin'], P['wmin'])
        tree = STRtree(polys); E = unary_union([p.boundary for p in polys]); d = []; io = []
        for m in M:
            io.append(max([p.intersection(m).area / p.union(m).area for p in (polys[i] for i in tree.query(m))], default=0))
            b = m.boundary; d += [b.interpolate(t).distance(E) for t in np.linspace(0, b.length, 40, endpoint=False)]
        s = float(np.mean(io)) + float((np.array(d) < 3).mean())
        print(f'  {q}  IoU medio {np.mean(io):.2f}  borde<3px {np.mean(np.array(d) < 3):.2f}', flush=True)
        if best is None or s > best[0]: best = (s, q)
    os.makedirs(a.salida, exist_ok=True)
    json.dump(best[1], open(os.path.join(a.salida, 'calibracion.json'), 'w'))
    print('Mejor:', json.dumps(best[1]), f'→ {a.salida}/calibracion.json (lo usa «lineas»)')

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    def common(p, lam=True):
        if lam:
            p.add_argument('--lamina', help='mosaico de la lámina (17008×21670)')
            p.add_argument('--teselas', help='carpeta tiles/orig/7 del repositorio')
        p.add_argument('--salida', default='SIG')
        p.add_argument('--manual', help='respaldo del editor (JSON de «Copiar respaldo»)')
    p = sp.add_parser('lineas'); common(p)
    p.add_argument('--zona', help='x0,y0,x1,y1 en px de la lámina (por defecto, toda)')
    p.add_argument('--tesela', type=int, default=1536); p.add_argument('--procesos', type=int, default=max(1, (os.cpu_count() or 2) - 1))
    p.add_argument('--param', help='JSON con parámetros a cambiar, p. ej. \'{"t":5,"gap":60}\'')
    p = sp.add_parser('parcelas'); common(p, False)
    p.add_argument('--lineas', help='líneas (editadas) en GeoJSON; por defecto SALIDA/lineas_deslinde.geojson')
    p.add_argument('--rotulos', default='data/transcripcion_rotulos.tsv')
    p.add_argument('--gcp', default='data/puntos_control.csv')
    p.add_argument('--max_area', type=float, default=400000, help='px²; recintos mayores se consideran abiertos (≈ 1300 km²)')
    p = sp.add_parser('lamina'); common(p)
    p = sp.add_parser('calibrar'); common(p)
    a = ap.parse_args()
    dict(lineas=cmd_lineas, parcelas=cmd_parcelas, lamina=cmd_lamina, calibrar=cmd_calibrar)[a.cmd](a)

if __name__ == '__main__':
    main()
