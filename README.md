# Disección aórtica (HP:0002647) desde la biología de sistemas

Proyecto de la asignatura Biología de Sistemas (Universidad de Málaga, curso 2026/27).

Se parte del término de la Human Phenotype Ontology **HP:0002647, *Aortic dissection***, se obtienen los genes asociados a él, se construye su red de interacciones de proteínas con STRING y se analiza esa red: propiedades topológicas, detección de módulos y enriquecimiento funcional.

Todo el análisis es reproducible con dos órdenes (`setup.sh` y `launch.sh`), o con Docker (`docker build` y `docker run`). Los resultados no se guardan en el repositorio: se generan cada vez en `results/`.

## Flujo de trabajo

```
API de HPO  ->  genes del término  ->  STRING 12.5 (STRINGdb)  ->  red en igraph
            ->  topología  ->  clustering (Louvain)  ->  enriquecimiento funcional
```

El trabajo tiene dos partes, separadas en carpetas distintas:

- **Análisis mínimo** (`code/minimo/`, `results/minimo/`): el flujo obligatorio de la asignatura.
- **Ampliaciones** (`code/ampliaciones/`, `results/ampliaciones/`): análisis adicionales sobre el mínimo. Cada una tiene su script `A<n>_...` y su propia subcarpeta de resultados.

## Estructura del repositorio

```
project_template/
  README.md
  Dockerfile              receta de la imagen de Docker (alternativa a setup.sh + launch.sh)
  .dockerignore           lo que no entra en la imagen (resultados, software instalado, .git, informe)
  code/
    setup.sh              instala los paquetes de R y de Python dentro de software/
    launch.sh             ejecuta todo el análisis y escribe en results/
    minimo/               scripts del análisis mínimo (01_ a 05_)
    ampliaciones/         scripts de las ampliaciones (A1_, A2_, ...)
  report/
    report.tex            informe (plantilla BMC)
    report.pdf            informe compilado
    tex_files/            un archivo .tex por apartado
    bibliography.bib      bibliografía
  results/                lo crea launch.sh (no se sube)
  software/               lo crea setup.sh (no se sube)
```

## Requisitos

- Linux (probado en Ubuntu sobre WSL en Windows 11).
- R 4.6 (probado con R 4.6.2).
- Python 3 y `curl`, que ya vienen con Ubuntu (probado con el Python 3.12.3 de Ubuntu 24.04). Solo los usa la ampliación A8. No hace falta el paquete `python3-venv`: `setup.sh` crea el entorno de Python sin él.
- Conexión a internet: el análisis consulta la API de HPO y la de ChEMBL (A9) y descarga los datos de STRING, las anotaciones gen-fenotipo de HPO y la ontología de HPO (`hp.obo`, A7 y A10).
- Cinco librerías del sistema que necesitan algunos paquetes de R para compilar (`libharfbuzz-dev` y `libfribidi-dev`, para el paquete `textshaping`, del que dependen de forma indirecta `enrichplot` y `clusterProfiler`). Es lo único que requiere permisos de administrador, y solo hay que hacerlo una vez:

```bash
sudo apt install -y libcurl4-openssl-dev libuv1-dev libcairo2-dev libharfbuzz-dev libfribidi-dev
```

Los paquetes de R los instala `setup.sh` (no hace falta instalarlos a mano): `igraph`, `optparse`, `httr`, `jsonlite` y `BiocManager` (CRAN); `STRINGdb`, `clusterProfiler`, `org.Hs.eg.db` y `enrichplot` (Bioconductor). Los de Python también, en un entorno virtual dentro del repositorio (`software/venv_py`), con las versiones probadas: `numpy`, `scipy`, `scikit-learn`, `igraph` y `matplotlib`.

## Cómo ejecutarlo

Antes, una sola vez, hay que instalar las librerías del sistema del apartado anterior. Después, desde una terminal, en una carpeta cualquiera:

```bash
git clone https://github.com/antoniodoa9/project_template.git
cd project_template
./code/setup.sh
./code/launch.sh
```

- `setup.sh` instala los paquetes de R en `software/deps_R` y los de Python en `software/venv_py`, dentro del repositorio, sin `sudo` y sin escribir nada fuera de él. La primera vez tarda unos 7 minutos; si se vuelve a ejecutar, solo instala lo que falte. Termina con el mensaje `setup.sh terminado`.
- `launch.sh` ejecuta los scripts en orden y escribe todos los resultados en `results/`. En la terminal muestra un título antes de cada parte (`===== ANÁLISIS MÍNIMO =====`, `== 02_red.R: red de STRING ==`, `===== AMPLIACIONES =====`, `== A1. Barrido de umbrales de STRING ==`...) seguido del resumen que imprime cada script, y termina con `launch.sh terminado`. La primera vez descarga los datos de STRING (unos 25 MB) en `results/string_data/`.

Los dos scripts funcionan desde cualquier carpeta, porque cada uno se sitúa en `code/` al empezar.

### Con Docker (alternativa)

Si se prefiere no instalar nada en el sistema, el `Dockerfile` de la raíz construye una imagen con todo lo necesario: parte de R 4.6.2 sobre Ubuntu (imagen `rocker/r-ver:4.6.2`), instala las librerías del sistema, `python3` y `curl`, y ejecuta `setup.sh` dentro de la imagen. Al arrancar, el contenedor ejecuta `launch.sh`. Solo hace falta tener Docker instalado y poder usarlo con el usuario normal (en Windows, Docker Desktop abierto y con la integración con WSL activada). Probado en Ubuntu sobre WSL con Docker Desktop; no se ha probado en Mac con procesador ARM.

Desde una terminal, en una carpeta cualquiera:

```bash
git clone https://github.com/antoniodoa9/project_template.git
cd project_template
docker build -t proyecto_hp0002647 .
mkdir -p results
docker run --rm --user "$(id -u):$(id -g)" -v "$(pwd)/results:/proyecto/results" proyecto_hp0002647
```

- `docker build` construye la imagen (la primera vez, unos 11 minutos; ocupa unos 3,4 GB). Si después solo cambia algún script de `code/`, se reutiliza la instalación de los paquetes y tarda unos segundos.
- `docker run` ejecuta `launch.sh` dentro del contenedor (unos 9 minutos). La imagen contiene el código y los programas, pero no los resultados: `-v` conecta la carpeta `results/` del ordenador con la del contenedor, así que los resultados aparecen en `results/`, igual que sin Docker. Sin `-v`, se perderían al terminar.
- `--user "$(id -u):$(id -g)"` hace que los archivos de `results/` sean del usuario y no de `root`. `--rm` borra el contenedor al terminar (la imagen se conserva).

## Parámetros del análisis mínimo

| Parámetro | Valor | Motivo |
|---|---|---|
| Término HPO | HP:0002647 (*Aortic dissection*) | Término elegido de la lista de la asignatura |
| Versión de STRING | 12.5, *Homo sapiens* (9606) | Versión fija, para que el resultado no cambie |
| Umbral de puntuación de STRING | 700 (confianza alta) | Es el corte de confianza alta de STRING. Con 800 se quedan sin ninguna interacción 13 de los 41 genes; con 700, 10 (ver la ampliación A1) |
| Red | No dirigida, sin interacciones repetidas | Una interacción entre dos proteínas no tiene dirección, y STRINGdb puede devolver la misma pareja más de una vez |
| Semilla | `set.seed(123)` | Louvain y la disposición de las figuras tienen una parte aleatoria |

## Resultados que se generan

Todas las rutas son relativas a `results/`.

Además, `session_info.txt` guarda la versión de R y de cada paquete con la que se ha hecho el análisis (salida de `sessionInfo()`).

### Análisis mínimo

| Archivo | Contenido | Script |
|---|---|---|
| `minimo/tablas/genes_hpo.tsv` | Genes asociados al término HPO (`id`, `name`) | `01_genes.R` |
| `minimo/tablas/red.tsv` | Interacciones entre esos genes con puntuación de 700 o más (`from`, `to`) | `02_red.R` |
| `minimo/tablas/genes_string.tsv` | Genes encontrados en STRING y su identificador (`name`, `STRING_id`) | `02_red.R` |
| `string_data/` | Archivos de STRING que descarga STRINGdb (se reutilizan en las siguientes ejecuciones) | `02_red.R` |
| `minimo/tablas/topologia_resumen.tsv` | Propiedades globales: nodos, nodos aislados, aristas, grado medio, densidad y componentes | `03_topologia.R` |
| `minimo/tablas/topologia_genes.tsv` | Por gen: grado, cercanía dentro de la componente mayor y centralidad armónica | `03_topologia.R` |
| `minimo/figuras/red_fr.pdf` | Red con disposición de Fruchterman-Reingold (tamaño del nodo según su grado) | `03_topologia.R` |
| `minimo/figuras/red_circular.pdf` | Red con disposición circular | `03_topologia.R` |
| `minimo/tablas/clusters.tsv` | Cluster de cada gen según Louvain (`gen`, `cluster`) | `04_clustering.R` |
| `minimo/figuras/clusters_louvain.pdf` | Red coloreada por cluster | `04_clustering.R` |
| `minimo/tablas/GO_lista_completa.csv` | Enriquecimiento GO (procesos biológicos) de todos los genes del término | `05_enriquecimiento.R` |
| `minimo/tablas/GO_cluster_<n>.csv` | Enriquecimiento GO de cada cluster de 5 genes o más | `05_enriquecimiento.R` |
| `minimo/tablas/GO_clusters_resumen.tsv` | Por cluster enriquecido: genes, número de términos y término más significativo | `05_enriquecimiento.R` |
| `minimo/figuras/GO_lista_completa_dotplot.pdf`, `GO_lista_completa_barplot.pdf` | Los 10 términos más significativos de la lista completa | `05_enriquecimiento.R` |
| `minimo/figuras/GO_cluster_<n>_dotplot.pdf` | Los 10 términos más significativos de cada cluster con algún término | `05_enriquecimiento.R` |

Los genes que están en STRING pero no tienen ninguna interacción por encima del umbral se mantienen en la topología como nodos aislados, y `02_red.R` los muestra por pantalla. El clustering se hace solo con los genes conectados, porque cada gen aislado formaría un cluster de un solo gen.

En el enriquecimiento se usa como fondo, para la lista completa, el conjunto de genes con alguna anotación en HPO (`genes_to_phenotype.txt`), y para cada cluster, los genes de la red. Los p-valores se corrigen por Benjamini-Hochberg y solo se enriquecen los clusters de 5 genes o más.

### Ampliaciones

**A1. Barrido de umbrales de STRING** (`code/ampliaciones/A1_barrido_umbrales.R`). Responde a cómo cambia la red según la confianza que se exige a las interacciones. `launch.sh` construye la red con seis umbrales (400, 500, 600, 700, 800 y 900) usando el mismo `02_red.R`, y el script compara las seis redes con las mismas medidas que el análisis mínimo.

| Archivo | Contenido | Script |
|---|---|---|
| `ampliaciones/A1_umbrales/redes/umbral_<u>/` | La red (`red.tsv`) y los genes encontrados en STRING (`genes_string.tsv`) con cada umbral | `02_red.R`, llamado desde `launch.sh` |
| `ampliaciones/A1_umbrales/barrido_umbrales.tsv` | Por umbral: genes con y sin interacción, aristas, densidad, nodos de la componente mayor, clusters de Louvain, clusters de 5 genes o más y modularidad | `A1_barrido_umbrales.R` |
| `ampliaciones/A1_umbrales/genes_por_umbral.tsv` | Por gen: si tiene alguna interacción con cada umbral y el umbral más alto al que la tiene | `A1_barrido_umbrales.R` |
| `ampliaciones/A1_umbrales/barrido_umbrales.pdf` | Genes conectados, interacciones, clusters y modularidad frente al umbral | `A1_barrido_umbrales.R` |

**A2. Red real frente a redes aleatorias** (`code/ampliaciones/A2_redes_aleatorias.R`). Responde a si la organización en módulos de la red del análisis mínimo es mayor de lo que daría el azar. Genera 1000 redes aleatorias con los mismos genes y el mismo número de interacciones de dos maneras: colocando las interacciones al azar (Erdos-Renyi, `sample_gnm`) y barajando las parejas de cada gen sin cambiar su número de interacciones (grado conservado, `rewire`). En cada red calcula la modularidad de Louvain, la transitividad, la distancia media y la asortatividad por grado, y las compara con las de la red real con un valor z y un p empírico.

| Archivo | Contenido | Script |
|---|---|---|
| `ampliaciones/A2_redes_aleatorias/A2_nulos.tsv` | Las cuatro propiedades de cada una de las 1000 redes aleatorias de cada modelo | `A2_redes_aleatorias.R` |
| `ampliaciones/A2_redes_aleatorias/A2_resumen.tsv` | Por propiedad y modelo: valor real, media y desviación de las redes aleatorias, valor z y p empírico | `A2_redes_aleatorias.R` |
| `ampliaciones/A2_redes_aleatorias/A2_redes_aleatorias.pdf` | Un histograma por propiedad y modelo con el valor de la red real marcado en rojo | `A2_redes_aleatorias.R` |

**A3. Ataque a los hubs frente a fallo al azar** (`code/ampliaciones/A3_ataque_hubs.R`). Responde a qué genes sostienen la conexión de la red y a si la red se rompe antes cuando se quitan esos genes que cuando se quitan genes al azar. Calcula el grado y la intermediación (`betweenness`) de cada gen de la red del análisis mínimo y llama hubs a los del 20 % de mayor grado y cuellos de botella a los del 20 % de mayor intermediación. Después quita los genes de uno en uno en tres órdenes (de más a menos grado, de más a menos intermediación y 1000 órdenes al azar) y mide el tamaño de la componente mayor tras cada paso. Compara los ataques con el azar al quitar el 20 % de los genes y con el número de genes que hay que quitar para que la componente mayor baje a la mitad, con un p empírico.

| Archivo | Contenido | Script |
|---|---|---|
| `ampliaciones/A3_ataque_hubs/A3_hubs_cuellos.tsv` | Por gen: grado, intermediación (también normalizada), módulo de Louvain, número de módulos distintos entre sus vecinos y tipo (hub, cuello de botella, los dos o ninguno) | `A3_ataque_hubs.R` |
| `ampliaciones/A3_ataque_hubs/A3_curvas.tsv` | Tras quitar 0, 1, 2... genes: componente mayor con cada ataque, media e intervalo del 95 % de los órdenes al azar, y número de componentes | `A3_ataque_hubs.R` |
| `ampliaciones/A3_ataque_hubs/A3_resumen.tsv` | Por estrategia: genes quitados al llegar al 20 %, componente mayor, componentes, p empírico, genes necesarios para partir la componente mayor por la mitad y su p empírico | `A3_ataque_hubs.R` |
| `ampliaciones/A3_ataque_hubs/A3_ataque_hubs.pdf` | Curva de la componente mayor frente a los genes quitados con los dos ataques y el azar, y grado frente a intermediación con los hubs y los cuellos de botella | `A3_ataque_hubs.R` |

**A4. Otros algoritmos de clustering** (`code/ampliaciones/A4_clusterings.R`). Responde a si los módulos dependen del algoritmo elegido. Aplica a la red del análisis mínimo siete algoritmos (Louvain, Leiden, fast greedy, walktrap, infomap, label propagation y edge betweenness) y compara sus particiones con la información mutua normalizada (NMI: 1 si dos algoritmos agrupan igual). También prueba Louvain con tres valores del parámetro `resolution` (0,5, 1 y 1,5) y aplica los siete algoritmos a las redes de cada umbral de la A1, para ver si los módulos dependen del umbral además del algoritmo.

| Archivo | Contenido | Script |
|---|---|---|
| `ampliaciones/A4_clusterings/A4_resumen.tsv` | Por algoritmo: número de clusters, modularidad, NMI frente a Louvain y tamaños | `A4_clusterings.R` |
| `ampliaciones/A4_clusterings/A4_nmi.tsv` | NMI entre cada pareja de algoritmos | `A4_clusterings.R` |
| `ampliaciones/A4_clusterings/A4_particiones.tsv` | Cluster de cada gen con cada algoritmo | `A4_clusterings.R` |
| `ampliaciones/A4_clusterings/A4_resolution.tsv` | Louvain con cada valor de `resolution`: clusters, modularidad, NMI frente a Louvain y tamaños | `A4_clusterings.R` |
| `ampliaciones/A4_clusterings/A4_umbral_algoritmo.tsv` | Por umbral de la A1 y algoritmo: nodos, número de clusters y modularidad | `A4_clusterings.R` |
| `ampliaciones/A4_clusterings/A4_clusterings.pdf` | Una página por algoritmo, con la misma disposición de los nodos | `A4_clusterings.R` |

**A5. Enriquecimiento ampliado** (`code/ampliaciones/A5_enriquecimiento.R`). Responde a qué funciones moleculares y componentes celulares, además de procesos biológicos, comparten los genes del término. Repite el enriquecimiento de la lista completa con las tres ontologías de GO (procesos biológicos, BP; funciones moleculares, MF; componentes celulares, CC), con el mismo fondo de HPO que `05_enriquecimiento.R`, y quita con `simplify()` los términos casi repetidos por la jerarquía de GO.

| Archivo | Contenido | Script |
|---|---|---|
| `ampliaciones/A5_enriquecimiento/GO_<ont>.csv` | Enriquecimiento de cada ontología (`BP`, `MF`, `CC`) | `A5_enriquecimiento.R` |
| `ampliaciones/A5_enriquecimiento/GO_<ont>_simplificado.csv` | Lo mismo tras quitar los términos redundantes con `simplify()` | `A5_enriquecimiento.R` |
| `ampliaciones/A5_enriquecimiento/A5_resumen.tsv` | Por ontología: términos antes y después de `simplify()` y término más significativo | `A5_enriquecimiento.R` |
| `ampliaciones/A5_enriquecimiento/A5_enriquecimiento.pdf` | Un dotplot por ontología, la red genes-términos (`cnetplot`) y la red de términos que comparten genes (`emapplot`) de BP | `A5_enriquecimiento.R` |

**A6. Genes candidatos con random walk con reinicio** (`code/ampliaciones/A6_rwr_roc.R`). Responde a qué otros genes, fuera de los 41 del término, están más cerca de ellos en la red de interacciones y podrían estar relacionados con la disección aórtica. Parte de los genes de HPO (las semillas) y hace un random walk con reinicio (probabilidad de volver a las semillas en cada paso: 0,7) sobre una subred de STRING formada por las semillas y sus vecinos directos (a un paso), con la puntuación de cada interacción como peso. Se usa esta subred y no la red humana completa porque con la red completa el random walk tarda mucho y se reparte por zonas lejanas a la enfermedad. Lo repite con tres umbrales (500, 700 y 900) para ver qué candidatos no dependen del umbral. Para comprobar si el método funciona, esconde genes de HPO (uno a uno, leave-one-out, y en 5 bloques, 5-fold), repite el random walk con los demás y mira en qué puesto quedan los escondidos frente al resto de genes de la subred: de ahí salen el AUC y la curva ROC, calculados a mano. Con el umbral 700 repite la validación 100 veces con genes al azar de grado parecido al de las semillas, para ver si el AUC se debe solo a que las semillas tienen muchas interacciones. Además guarda la subred de 700, que usa la A8.

| Archivo | Contenido | Script |
|---|---|---|
| `ampliaciones/A6_rwr_roc/A6_candidatos.tsv` | Por umbral (500, 700 y 900): los 20 genes candidatos con más puntuación del random walk con reinicio, con su grado, con cuántos genes de HPO interaccionan directamente, el gen de HPO vecino que más puntuación le pasa en un paso y el grado de ese gen | `A6_rwr_roc.R` |
| `ampliaciones/A6_rwr_roc/A6_resumen_umbral.tsv` | Por umbral: genes e interacciones de la subred (semillas y vecinos directos), y genes de HPO que están y que no están en ella | `A6_rwr_roc.R` |
| `ampliaciones/A6_rwr_roc/A6_validacion.tsv` | Por umbral: el AUC de cada gen de HPO al esconderlo, con leave-one-out y con 5-fold | `A6_rwr_roc.R` |
| `ampliaciones/A6_rwr_roc/A6_auc.tsv` | Por umbral: el AUC medio con leave-one-out y con 5-fold | `A6_rwr_roc.R` |
| `ampliaciones/A6_rwr_roc/A6_roc.pdf` | Una página por umbral: la curva ROC con leave-one-out y con 5-fold, y la diagonal de un clasificador al azar (AUC 0,5); y una última página con el histograma de los AUC con genes al azar y el AUC real en rojo | `A6_rwr_roc.R` |
| `ampliaciones/A6_rwr_roc/A6_azar.tsv` | Con el umbral 700: el AUC medio (leave-one-out) de cada una de las 100 repeticiones con genes al azar de grado parecido al de las semillas | `A6_rwr_roc.R` |
| `ampliaciones/A6_rwr_roc/A6_subred_700.tsv` | La subred de 700, una fila por interacción: los dos genes (`from`, `to`) y la puntuación de STRING (`combined_score`). La usa la A8 | `A6_rwr_roc.R` |

**A7. Red bipartita gen-fenotipo** (`code/ampliaciones/A7_red_bipartita.R`). Responde a si los genes que interaccionan en STRING, o que están en el mismo módulo, también se parecen en los fenotipos que causan (relación entre genotipo y fenotipo). Construye una red bipartita con dos tipos de nodos, los 41 genes y sus fenotipos de HPO (`genes_to_phenotype.txt`), con una arista entre cada gen y cada fenotipo que tiene anotado. Se quitan el propio HP:0002647, que tienen todos, y los modos de herencia (la rama HP:0000005 de `hp.obo`, por ejemplo *Autosomal dominant inheritance*), que no son fenotipos. La red se proyecta sobre los genes (`bipartite_projection`) y el parecido entre dos genes se mide con el índice de Jaccard: fenotipos compartidos entre fenotipos de los dos juntos. Compara, frente a 1000 permutaciones de los nombres de los genes sobre la misma red y los mismos módulos (modelo al azar basado en los nodos), el Jaccard medio de las interacciones de STRING y el de los pares del mismo módulo de Louvain. Además agrupa la red de fenotipos con Louvain (con el Jaccard como peso) y compara esos módulos con los de STRING mediante la NMI, también frente a 1000 permutaciones.

| Archivo | Contenido | Script |
|---|---|---|
| `ampliaciones/A7_red_bipartita/A7_resumen.tsv` | Jaccard medio de las interacciones de STRING y del mismo módulo, y NMI entre los módulos de fenotipos y los de STRING: valor real, media y máximo al azar, p empírico, genes y pares | `A7_red_bipartita.R` |
| `ampliaciones/A7_red_bipartita/A7_bipartita.tsv` | Tamaño de la red bipartita (genes, fenotipos, aristas, anotaciones de herencia quitadas, pares de genes que comparten algún fenotipo) y versión de HPO | `A7_red_bipartita.R` |
| `ampliaciones/A7_red_bipartita/A7_genes.tsv` | Por gen: número de fenotipos, módulo de STRING y módulo de la red de fenotipos | `A7_red_bipartita.R` |
| `ampliaciones/A7_red_bipartita/A7_fenotipos_compartidos.tsv` | Los 30 fenotipos que comparten más genes | `A7_red_bipartita.R` |
| `ampliaciones/A7_red_bipartita/A7_pares_jaccard.tsv` | Cada par de genes que comparte algún fenotipo: fenotipos compartidos, Jaccard y si interaccionan en STRING | `A7_red_bipartita.R` |
| `ampliaciones/A7_red_bipartita/A7_azar.tsv` | Los tres valores de cada una de las 1000 permutaciones | `A7_red_bipartita.R` |
| `ampliaciones/A7_red_bipartita/A7_red_bipartita.pdf` | Los tres histogramas al azar con el valor real en rojo, y la matriz de Jaccard entre los genes de la red ordenados por módulo | `A7_red_bipartita.R` |

**A8. Embedding de la subred** (`code/ampliaciones/A8_embedding.py`, en Python). Responde a si una representación de cada gen como vector (embedding), aprendida de los paseos aleatorios por la red como en DeepWalk y node2vec, sirve para predecir interacciones y genes de la enfermedad. Parte de la subred de 700 de la A6. Hace 10 paseos aleatorios de 80 pasos desde cada gen, cuenta qué genes aparecen juntos (ventana de 10) y obtiene un vector de 128 números por gen con la PMI positiva y una SVD, que es lo que calcula word2vec skip-gram. Usa los vectores de tres formas, cada una comparada con una referencia: esconde el 10 % de las interacciones, rehace el embedding sin ellas y mira si esas parejas quedan más cerca que parejas al azar (frente al índice de Jaccard); esconde las semillas una a una (leave-one-out, como la A6) y las puntúa por su parecido con las demás (frente al random walk de la A6 y frente al grado); y agrupa los vectores con k-means y compara los grupos con Louvain (NMI, frente a 100 permutaciones). Los vectores se dibujan con PCA solo para verlos: las distancias se miden siempre con los vectores completos.

| Archivo | Contenido | Script |
|---|---|---|
| `ampliaciones/A8_embedding/A8_resumen.tsv` | Tamaño de la subred, AUC de las interacciones escondidas (embedding y Jaccard), AUC de las semillas escondidas (embedding con la media de las semillas o con la semilla más cercana, random walk de la A6 y grado), NMI entre k-means y Louvain y su media al azar, y candidatos comunes con la A6 | `A8_embedding.py` |
| `ampliaciones/A8_embedding/A8_validacion.tsv` | El AUC de cada gen de HPO al esconderlo, con las dos formas de puntuar | `A8_embedding.py` |
| `ampliaciones/A8_embedding/A8_candidatos.tsv` | Los 20 genes que no son de HPO más parecidos (coseno) a la media de las semillas, con su grado y si también son candidatos de la A6 con 700 | `A8_embedding.py` |
| `ampliaciones/A8_embedding/A8_pca.pdf` | Los vectores de los genes dibujados con PCA, con las semillas en rojo | `A8_embedding.py` |

**A9. Fármacos de ChEMBL** (`code/ampliaciones/A9_chembl.R`). Responde a si los genes del término y los candidatos de la A6 y de la A8 son diana de fármacos aprobados o en ensayos clínicos más a menudo de lo que daría el azar. Descarga de la API de ChEMBL todos los mecanismos de acción con fase clínica 1 o más (4 es aprobado) y las dianas humanas que son una sola proteína, con su símbolo de gen, y cuenta qué genes de cada conjunto (semillas, candidatos de la A6 con 700 y candidatos de la A8) tienen algún fármaco. Para saber si son más de lo esperado, repite el recuento 1000 veces con genes de la subred de 700 que no son de HPO y tienen un grado parecido (mismo decil), porque los fármacos suelen ir a genes muy estudiados, que tienen muchas interacciones (p empírico, como en la A6). Guarda la versión de ChEMBL y la fecha de descarga.

| Archivo | Contenido | Script |
|---|---|---|
| `ampliaciones/A9_chembl/A9_farmacos.tsv` | Una fila por conjunto, gen y fármaco: nombre e identificador de ChEMBL del fármaco, fase máxima, tipo de acción y mecanismo | `A9_chembl.R` |
| `ampliaciones/A9_chembl/A9_resumen.tsv` | Por conjunto: genes, genes con fármaco, media y máximo al azar, p empírico, versión de ChEMBL y fecha de descarga | `A9_chembl.R` |
| `ampliaciones/A9_chembl/A9_azar.tsv` | Genes con fármaco en cada una de las 1000 repeticiones al azar de cada conjunto | `A9_chembl.R` |
| `ampliaciones/A9_chembl/A9_azar.pdf` | Un histograma por conjunto con los recuentos al azar y el real en rojo | `A9_chembl.R` |

**A10. Evidencia de los genes y red previa** (`code/ampliaciones/A10_evidencia_genes.R`). Responde a si todos los genes del término están igual de relacionados con él. Un gen entra en HP:0002647 porque alguna enfermedad causada por ese gen tiene la disección aórtica entre sus fenotipos, y HPO indica con qué frecuencia la presentan los pacientes (`genes_to_phenotype.txt`: *Very frequent* 80-99 %, *Frequent* 30-79 %, *Occasional* 5-29 %, *Very rare* 1-4 %, o una fracción de pacientes n/m). El script toma, para cada gen, la frecuencia más alta en el término y sus descendientes (de `hp.obo`) y lo clasifica en evidencia frecuente (30 % o más), ocasional (5-29 %) o muy rara (menos del 5 %). Después construye una red previa solo con los genes de evidencia frecuente (sus interacciones de STRING con 700), repite Louvain y la compara con la red completa (modularidad y NMI de los módulos en los genes comunes). Por último cruza la evidencia con el AUC *leave-one-out* de cada semilla en la A6 (umbral 700), para ver si la red recupera peor los genes de evidencia débil.

| Archivo | Contenido | Script |
|---|---|---|
| `ampliaciones/A10_evidencia_genes/A10_evidencia_genes.tsv` | Por gen: frecuencia máxima, grupo de evidencia, enfermedades (OMIM u Orphanet) por las que está en el término, módulo y grado en la red del mínimo y AUC de la A6 con 700 | `A10_evidencia_genes.R` |
| `ampliaciones/A10_evidencia_genes/A10_resumen.tsv` | Versión de HPO, genes de cada grupo, y genes, interacciones, módulos y modularidad de la red completa y de la red previa, con la NMI entre sus módulos | `A10_evidencia_genes.R` |
| `ampliaciones/A10_evidencia_genes/A10_modulos_red_fuerte.tsv` | Módulo de cada gen en la red previa y en la red completa | `A10_evidencia_genes.R` |
| `ampliaciones/A10_evidencia_genes/A10_evidencia.pdf` | La red del mínimo con cada gen coloreado por su evidencia, y el AUC de cada semilla en la A6 según su evidencia (en rojo, las que no tienen interacciones en la red del mínimo) | `A10_evidencia_genes.R` |

## Ejecutar un script por separado

Cada script se puede ejecutar también por separado desde `code/`. Con `--help` se ven sus opciones, por ejemplo:

```bash
cd code
Rscript minimo/02_red.R --help
../software/venv_py/bin/python ampliaciones/A8_embedding.py --help
```

## Solución de problemas

| Problema | Causa | Solución |
|---|---|---|
| `setup.sh` falla al compilar un paquete y el error menciona `curl`, `uv`, `cairo`, `harfbuzz`, `fribidi` o `textshaping`, o termina con `Faltan paquetes: clusterProfiler, enrichplot` | Faltan librerías del sistema | `sudo apt install -y libcurl4-openssl-dev libuv1-dev libcairo2-dev libharfbuzz-dev libfribidi-dev` y volver a ejecutar `./code/setup.sh` |
| `setup.sh` se para con `Faltan paquetes: ...` | Se cortó una descarga (por ejemplo `org.Hs.eg.db`, que es grande) | Volver a ejecutar `./code/setup.sh`: solo instala lo que falte |
| `setup.sh` falla al crear el entorno de Python o al descargar `get-pip.py` | No hay conexión con `bootstrap.pypa.io` o el sistema no tiene `python3` o `curl` | Comprobar la conexión y que `python3 --version` y `curl --version` funcionan, y volver a ejecutar `./code/setup.sh` |
| `docker: Cannot connect to the Docker daemon` o `docker: command not found` | Docker no está en marcha (en Windows, Docker Desktop cerrado o sin la integración con WSL) | Abrir Docker Desktop y, en Settings, Resources, WSL integration, activar la distribución de Ubuntu; o en Linux, `sudo systemctl start docker` |
| `permission denied while trying to connect to the Docker daemon socket` | El usuario no está en el grupo `docker` | Pedir al administrador que lo añada (`sudo usermod -aG docker $USER`) y volver a entrar en la sesión |
| `Permission denied` al ejecutar `./code/setup.sh` o `./code/launch.sh` | El archivo perdió el permiso de ejecución (por ejemplo, al copiarlo desde Windows) | `chmod +x code/setup.sh code/launch.sh` |
| `A9_chembl.R` termina con un error HTTP (por ejemplo 500 o 503) y `launch.sh` muestra `AVISO: la A9 no ha terminado` | La API de ChEMBL del EBI no responde en ese momento (el script ya repite cada descarga hasta 5 veces). `launch.sh` sigue con el resto del análisis | Esperar unos minutos y repetir solo la A9 desde `code/`: `export R_LIBS_USER=$(pwd)/../software/deps_R` y después la línea de la A9 de `launch.sh` (`Rscript ampliaciones/A9_chembl.R ...`) |
| `01_genes.R` termina con un error HTTP | La API de HPO no responde en ese momento | Comprobar la conexión y volver a ejecutar `./code/launch.sh` |
| `05_enriquecimiento.R` o `A5_enriquecimiento.R` se paran con `Timeout of 60 seconds was reached` | La descarga de `genes_to_phenotype.txt` de HPO tarda más de lo que R espera por defecto (60 segundos) | `launch.sh` ya da hasta 10 minutos a cada descarga (`R_DEFAULT_INTERNET_TIMEOUT=600`). Si se ejecuta un script suelto, antes hacer `export R_DEFAULT_INTERNET_TIMEOUT=600`. Si aun así falla, la conexión es muy lenta: volver a ejecutar `./code/launch.sh` |
| Durante `05_enriquecimiento.R` y `A5_enriquecimiento.R` salen mensajes `'select()' returned 1:1 mapping...` y los avisos `0.06% of input gene IDs are fail to map` y `semData is not provided` | Son mensajes normales de `clusterProfiler`: al traducir los nombres de los genes a identificadores Entrez, unos pocos genes del fondo no tienen identificador; y `simplify()` calcula por su cuenta el parecido entre términos de GO | Nada: el análisis es correcto |
| `there is no package called ...` al ejecutar un script suelto | Se ha ejecutado sin la carpeta de paquetes del proyecto | Ejecutar con `./code/launch.sh`, o antes hacer `export R_LIBS_USER=$(pwd)/../software/deps_R` desde `code/` |

## Autores

| Persona | GitHub | Parte del trabajo |
|---|---|---|
| Antonio Sánchez Doña (coordinador) | [@antoniodoa9](https://github.com/antoniodoa9) | Código: instalación y ejecución (también con Docker), genes, red y topología, y sus ampliaciones |
| Nicolás Sánchez Alonso | [@NicoSA03](https://github.com/NicoSA03) | Código: clustering y enriquecimiento, y sus ampliaciones. Informe: Material y métodos y Resultados de A6, A8 y A9 |
| Eva Lucía López Porras | [@evalu5](https://github.com/evalu5) | Informe: Introducción y Discusión |
| Aitana Nieto Reus | [@aitananietoreus](https://github.com/aitananietoreus) | Informe: Material y métodos y Resultados (análisis mínimo y resto de ampliaciones) y Conclusiones |

Repositorio creado a partir de la plantilla [seoanezonjic/project_template](https://github.com/seoanezonjic/project_template).
