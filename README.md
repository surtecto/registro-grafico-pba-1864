# Registro Gráfico 1864 · Provincia de Buenos Aires · SIG

SIG histórico del **Registro Gráfico de las propiedades rurales de la Provincia de Buenos Aires** (Departamento Topográfico, 1864), georreferenciado sobre imagen satelital, con visor web para consulta y estadística.

**Visor:** abrir `index.html` (o la página de GitHub Pages del repositorio).

## Contenido

| Ruta | Qué es |
|---|---|
| `index.html` | Visor web (Leaflet): mapa georreferenciado sobre Esri / Google / OSM, lámina original navegable, búsqueda, filtros por partido y tipo de titular, estadística y método |
| `data/propietarios.geojson` · `.csv` | 1080 rótulos de propietario (978 titulares distintos), con tipo de titular, apellido y partido |
| `data/parcelas_aprox.geojson` | Zonas de propiedad aproximadas (Voronoi por rótulo, recortado por partido) |
| `data/partidos_aprox.geojson` | 70 partidos localizados, con límite aproximado (Voronoi de los números de partido) |
| `data/poblados.geojson` | 57 pueblos y ciudades |
| `data/hidrografia.geojson` | 119 lagunas, arroyos, ríos y cañadas |
| `data/geografia.geojson` | Fortines, sierras, estancias, médanos, puertos, parajes |
| `data/rotulos.geojson` | Números y rótulos de partidos y secciones |
| `data/puntos_control.geojson` · `.csv` | 44 puntos de control del georreferenciado |
| `data/transcripcion_rotulos.tsv` | Transcripción cruda: `categoría|texto|x|y` en píxeles de la lámina AGN (8504 × 10835) |
| `data/estadisticas.json` | Estadística agregada |
| `tiles/geo/` | Teselas XYZ (WebP, z5–11) de la lámina georreferenciada, EPSG:3857 |
| `tiles/orig/` | Pirámide de la lámina original sin deformar |
| `scripts/` | Georreferenciado, armado de capas y teselado (Python) |

Todas las capas vectoriales están en EPSG:4326. Cada registro conserva `px`/`py`, su posición en la lámina original.

## Método y precisión

- La lámina está dibujada con el sur hacia arriba. Se usaron 44 puntos de control (pueblos cabecera, Buenos Aires, Colonia, Cabo Corrientes, Bahía Blanca, etc.) con una transformación polinómica de 2.º grado más corrección *thin plate spline* de residuos.
- Error por validación cruzada: unos 7 km RMS, con mediana de 4 km. Es mayor hacia Bahía Blanca y en los bordes. Es la precisión esperable de un plano de 1864 compilado a partir de mensuras heterogéneas.
- **Partidos y zonas de propiedad son aproximaciones**, no un calcado de los deslindes. Sirven para asignar cada rótulo a su partido y para análisis de densidad.
- La transcripción no es exhaustiva: en los partidos del norte, con suertes muy angostas, quedan rótulos menores sin transcribir. Falta localizar el partido 25 (Arenales). El 32 (Patagones) está en la lámina inserta, fuera del georreferenciado.

## Fuente

Reproducción digital del Archivo General de la Nación (Registro Gráfico, 1864; Salas, Moreno, Kuhr, Benoit, Casagemas, Malaver).

Arq. Juan Patricio Addesso · Surtectura · Doctorado en Geografía (USAL)
