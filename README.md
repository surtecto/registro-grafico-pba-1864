# Registro Gráfico 1864 · Provincia de Buenos Aires · SIG

SIG histórico del **Registro Gráfico de las propiedades rurales de la Provincia de Buenos Aires** (Departamento Topográfico, 1864), georreferenciado sobre imagen satelital, con visor web para consulta y estadística.

**Visor:** abrir `index.html` (o la página de GitHub Pages del repositorio).

## Contenido

| Ruta | Qué es |
|---|---|
| `index.html` | Visor web (Leaflet): mapa georreferenciado sobre Esri / Google / OSM, lámina original navegable, búsqueda, filtros por partido y tipo de titular, estadística y método |
| `data/propietarios.geojson` · `.csv` | 1080 rótulos de propietario (978 titulares distintos), con tipo de titular, apellido y partido |
| `data/parcelas.geojson` | 7013 parcelas trazadas sobre los deslindes del plano; 458 con un propietario asignado, 119 con varios nombres (deslinde a revisar), el resto sin rótulo transcripto. Con superficie en ha |
| `data/partidos.geojson` | 70 partidos, cada uno como unión de sus parcelas |
| `data/poblados.geojson` | 57 pueblos y ciudades |
| `data/hidrografia.geojson` | 119 lagunas, arroyos, ríos y cañadas |
| `data/geografia.geojson` | Fortines, sierras, estancias, médanos, puertos, parajes |
| `data/rotulos.geojson` | Números y rótulos de partidos y secciones |
| `data/puntos_control.geojson` · `.csv` | 44 puntos de control del georreferenciado |
| `data/transcripcion_rotulos.tsv` | Transcripción cruda: `categoría|texto|x|y` en píxeles de la lámina AGN (8504 × 10835) |
| `data/estadisticas.json` | Estadística agregada |
| `tiles/geo/` | Teselas XYZ (WebP, z5–12) de la lámina georreferenciada, EPSG:3857 |
| `tiles/orig/` | Pirámide del mosaico de la lámina sin deformar (17008 × 21670 px) |
| `scripts/` | Georreferenciado, armado de capas y teselado (Python) |

Todas las capas vectoriales están en EPSG:4326. Cada registro conserva `px`/`py`, su posición en la lámina original.

## Método y precisión

- La lámina está dibujada con el sur hacia arriba. Se usaron 44 puntos de control (pueblos cabecera, Buenos Aires, Colonia, Cabo Corrientes, Bahía Blanca, etc.) con una transformación polinómica de 2.º grado más corrección *thin plate spline* de residuos.
- Error por validación cruzada: unos 7 km RMS, con mediana de 4 km. Es mayor hacia Bahía Blanca y en los bordes. Es la precisión esperable de un plano de 1864 compilado a partir de mensuras heterogéneas.
- **Parcelas:** recintos cerrados por las líneas de deslinde, detectadas sobre un mosaico de las 64 hojas de la digitalización original (Archivo Histórico Municipal de Olavarría, 2008) registradas contra la copia del AGN. El campo `estado` indica la confiabilidad de cada polígono. Los pliegues del papel generan algunas divisiones falsas.
- **Partidos:** unión de las parcelas de cada partido. La pertenencia de las parcelas de borde se asignó por cercanía al número del partido y debe revisarse contra la línea de puntos y rayas del plano.
- La transcripción no es exhaustiva: en los partidos del norte, con suertes muy angostas, quedan rótulos menores sin transcribir. Falta localizar el partido 25 (Arenales). El 32 (Patagones) está en la lámina inserta, fuera del georreferenciado.

## Fuente

Registro Gráfico, 1864 (Salas, Moreno, Kuhr, Benoit, Casagemas, Malaver). Copias digitales: Archivo General de la Nación y CD del Archivo Histórico Municipal de Olavarría (2008).

Arq. Juan Patricio Addesso · Surtectura · Doctorado en Geografía (USAL)
