# Parcelas del Registro Gráfico 1864 · script para QGIS

`parcelas_1864.py` arma la capa de parcelas a partir de las **líneas** del plano, no inundando desde un clic. Por eso no se escapa a la parcela vecina, no deja muescas en los nombres y las parcelas vecinas comparten exactamente el mismo borde (sin rendijas ni solapes).

## Instalación (una vez)

En la **OSGeo4W Shell** de QGIS (o cualquier Python 3.9+):

```
pip install numpy opencv-python-headless shapely scikit-image
```

Bajá el repositorio `surtecto/registro-grafico-pba-1864` (Code → Download ZIP) y descomprimilo. El script lee la lámina desde `tiles/orig/7` y los puntos de control desde `data/`. Poné estos archivos en la carpeta del repositorio.

## Flujo

1. **Lámina para QGIS** (una vez):
   `python parcelas_1864.py lamina --teselas tiles/orig/7 --salida SIG`
   Genera `SIG/lamina_1864.jpg` y su `.jgw`. En QGIS: clic derecho → *Crear pirámides* para que ande fluido.

2. **(Opcional) Calibrar con tus parcelas manuales:**
   `python parcelas_1864.py calibrar --teselas tiles/orig/7 --manual respaldo_editor.json --salida SIG`
   Prueba combinaciones de parámetros en la zona de tus parcelas y guarda la mejor en `SIG/calibracion.json`.

3. **Líneas de deslinde:**
   `python parcelas_1864.py lineas --teselas tiles/orig/7 --manual respaldo_editor.json --salida SIG`
   Tarda unos 10 minutos para toda la lámina con 2 núcleos.

4. **Corregir LÍNEAS en QGIS** (no polígonos): abrí `lamina_1864.jpg` y `lineas_deslinde.geojson`, activá la edición y el autoensamblado (snapping) a vértices y segmentos.
   - Borrá las líneas falsas: rayas sobre textos y caminos que no son deslindes.
   - Dibujá las que faltan. Alcanza con que la línea cruce un poco la línea vecina; el sobrante se descarta solo.

5. **Parcelas:**
   `python parcelas_1864.py parcelas --manual respaldo_editor.json --rotulos rotulos_1864.json --salida SIG`
   Genera:
   - `parcelas.geojson`: en coordenadas de lámina, para seguir trabajando sobre la imagen.
   - `parcelas_wgs84.geojson`: georreferenciada, con `sup_ha`.
   - `partidos_wgs84.geojson`: disolución de las parcelas que tienen partido asignado.

Repetí los pasos 4 y 5 las veces que haga falta. Tus parcelas del editor (`--manual`) siempre mandan: quedan fijas y sus bordes entran a la red de líneas.

## Campo `calidad` (para pintar en QGIS)

| valor | qué significa | qué hacer |
|---|---|---|
| `manual` | la dibujaste vos | nada |
| `ok` | un solo rótulo transcripto adentro | revisar el nombre |
| `varios_rotulos` | dos o más nombres adentro: casi seguro faltan una o más líneas | dibujar el deslinde que falta |
| `sin_rotulo` | parcela cerrada sin nombre transcripto | cargar el propietario |
| `agua` | río, laguna o bañado | ignorar o borrar |

Las zonas grandes sin parcelas son recintos que no cerraron y se descartaron para no inventar polígonos. Se arreglan dibujando la línea que falta.

## Coordenadas

- Las capas de lámina usan x = columna e y = −fila del mosaico de 17008 × 21670 px. Asignales a todas el mismo SRC, cualquiera métrico; por ejemplo EPSG:3857, solo como "plano".
- El georreferenciado usa los 44 puntos de control: polinomio de 2.º grado y *thin-plate spline*. Su error es de unos 4 km de mediana, propio del plano de 1864.

## Partidos

La línea de rayas y puntos va pegada a los deslindes y no se puede seguir sola de forma confiable. El partido de cada parcela se asigna así:
- en el editor (las manuales ya lo traen);
- en QGIS: seleccionar las parcelas → calculadora de campos → `partido_n`.

`partidos_wgs84.geojson` es la disolución de esas parcelas.
