# Disección aórtica (HP:0002647) desde la biología de sistemas

Proyecto de la asignatura Biología de Sistemas (Universidad de Málaga, curso 2026/27).

Se parte del término de la Human Phenotype Ontology **HP:0002647, *Aortic dissection***, se obtienen los genes asociados a él, se construye su red de interacciones de proteínas con STRING y se analiza esa red: propiedades topológicas, detección de módulos y enriquecimiento funcional.

Todo el análisis es reproducible con dos órdenes (`setup.sh` y `launch.sh`). Los resultados no se guardan en el repositorio: se generan cada vez en `results/`.

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
  code/
    setup.sh              instala los paquetes de R dentro de software/
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
- Conexión a internet: el análisis consulta la API de HPO y descarga los datos de STRING y las anotaciones gen-fenotipo de HPO.
- Tres librerías del sistema que necesitan algunos paquetes de Bioconductor para compilar. Es lo único que requiere permisos de administrador, y solo hay que hacerlo una vez:

```bash
sudo apt install -y libcurl4-openssl-dev libuv1-dev libcairo2-dev
```

Los paquetes de R los instala `setup.sh` (no hace falta instalarlos a mano): `igraph`, `optparse`, `httr`, `jsonlite` y `BiocManager` (CRAN); `STRINGdb`, `clusterProfiler`, `org.Hs.eg.db` y `enrichplot` (Bioconductor).

## Cómo ejecutarlo

Antes, una sola vez, hay que instalar las librerías del sistema del apartado anterior. Después, desde una terminal, en una carpeta cualquiera:

```bash
git clone https://github.com/antoniodoa9/project_template.git
cd project_template
./code/setup.sh
./code/launch.sh
```

- `setup.sh` instala los paquetes de R en `software/deps_R`, dentro del repositorio, sin `sudo` y sin escribir nada fuera de él. La primera vez tarda unos 7 minutos; si se vuelve a ejecutar, solo instala lo que falte. Termina con el mensaje `setup.sh terminado`.
- `launch.sh` ejecuta los scripts en orden y escribe todos los resultados en `results/`. En la terminal muestra un título antes de cada parte (`===== ANÁLISIS MÍNIMO =====`, `== 02_red.R: red de STRING ==`, `===== AMPLIACIONES =====`, `== A1. Barrido de umbrales de STRING ==`...) seguido del resumen que imprime cada script, y termina con `launch.sh terminado`. La primera vez descarga los datos de STRING (unos 25 MB) en `results/string_data/`.

Los dos scripts funcionan desde cualquier carpeta, porque cada uno se sitúa en `code/` al empezar.

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

**A6. Genes candidatos con random walk con reinicio** (`code/ampliaciones/A6_rwr_roc.R`). Responde a qué otros genes, fuera de los 41 del término, están más cerca de ellos en la red de interacciones y podrían estar relacionados con la disección aórtica. Parte de los genes de HPO (las semillas) y hace un random walk con reinicio (probabilidad de volver a las semillas en cada paso: 0,7) sobre la red humana completa de STRING, con la puntuación de cada interacción como peso. Lo repite con tres umbrales (500, 700 y 900) para ver qué candidatos no dependen del umbral. Para comprobar si el método funciona, esconde genes de HPO (uno a uno, leave-one-out, y en 5 bloques, 5-fold), repite el random walk con los demás y mira en qué puesto quedan los escondidos frente al resto de genes de la red: de ahí salen el AUC y la curva ROC, calculados a mano.

| Archivo | Contenido | Script |
|---|---|---|
| `ampliaciones/A6_rwr_roc/A6_candidatos.tsv` | Por umbral (500, 700 y 900): los 20 genes candidatos con más puntuación del random walk con reinicio, con su grado, con cuántos genes de HPO interaccionan directamente, el gen de HPO vecino que más puntuación le pasa en un paso y el grado de ese gen | `A6_rwr_roc.R` |
| `ampliaciones/A6_rwr_roc/A6_resumen_umbral.tsv` | Por umbral: genes e interacciones de la red completa de STRING, y genes de HPO que están y que no están en ella | `A6_rwr_roc.R` |
| `ampliaciones/A6_rwr_roc/A6_validacion.tsv` | Por umbral: el AUC de cada gen de HPO al esconderlo, con leave-one-out y con 5-fold | `A6_rwr_roc.R` |
| `ampliaciones/A6_rwr_roc/A6_auc.tsv` | Por umbral: el AUC medio con leave-one-out y con 5-fold | `A6_rwr_roc.R` |
| `ampliaciones/A6_rwr_roc/A6_roc.pdf` | Una página por umbral: la curva ROC con leave-one-out y con 5-fold, y la diagonal de un clasificador al azar (AUC 0,5) | `A6_rwr_roc.R` |

## Ejecutar un script por separado

Cada script se puede ejecutar también por separado desde `code/`. Con `--help` se ven sus opciones, por ejemplo:

```bash
cd code
Rscript minimo/02_red.R --help
```

## Solución de problemas

| Problema | Causa | Solución |
|---|---|---|
| `setup.sh` falla al compilar un paquete y el error menciona `curl`, `uv` o `cairo` | Faltan librerías del sistema | `sudo apt install -y libcurl4-openssl-dev libuv1-dev libcairo2-dev` y volver a ejecutar `./code/setup.sh` |
| `setup.sh` se para con `Faltan paquetes: ...` | Se cortó una descarga (por ejemplo `org.Hs.eg.db`, que es grande) | Volver a ejecutar `./code/setup.sh`: solo instala lo que falte |
| `Permission denied` al ejecutar `./code/setup.sh` o `./code/launch.sh` | El archivo perdió el permiso de ejecución (por ejemplo, al copiarlo desde Windows) | `chmod +x code/setup.sh code/launch.sh` |
| `01_genes.R` termina con un error HTTP | La API de HPO no responde en ese momento | Comprobar la conexión y volver a ejecutar `./code/launch.sh` |
| `05_enriquecimiento.R` o `A5_enriquecimiento.R` se paran con `Timeout of 60 seconds was reached` | La descarga de `genes_to_phenotype.txt` de HPO tarda más de lo que R espera por defecto (60 segundos) | `launch.sh` ya da hasta 10 minutos a cada descarga (`R_DEFAULT_INTERNET_TIMEOUT=600`). Si se ejecuta un script suelto, antes hacer `export R_DEFAULT_INTERNET_TIMEOUT=600`. Si aun así falla, la conexión es muy lenta: volver a ejecutar `./code/launch.sh` |
| Durante `05_enriquecimiento.R` y `A5_enriquecimiento.R` salen mensajes `'select()' returned 1:1 mapping...` y los avisos `0.06% of input gene IDs are fail to map` y `semData is not provided` | Son mensajes normales de `clusterProfiler`: al traducir los nombres de los genes a identificadores Entrez, unos pocos genes del fondo no tienen identificador; y `simplify()` calcula por su cuenta el parecido entre términos de GO | Nada: el análisis es correcto |
| `there is no package called ...` al ejecutar un script suelto | Se ha ejecutado sin la carpeta de paquetes del proyecto | Ejecutar con `./code/launch.sh`, o antes hacer `export R_LIBS_USER=$(pwd)/../software/deps_R` desde `code/` |

## Autores

| Persona | GitHub | Parte del trabajo |
|---|---|---|
| Antonio Sánchez Doña (coordinador) | [@antoniodoa9](https://github.com/antoniodoa9) | Código: instalación y ejecución, genes, red y topología, y sus ampliaciones |
| Nicolás Sánchez Alonso | [@NicoSA03](https://github.com/NicoSA03) | Código: clustering y enriquecimiento, y sus ampliaciones |
| Eva Lucía López Porras | [@evalu5](https://github.com/evalu5) | Informe: Introducción y Discusión |
| Aitana Nieto Reus | [@aitananietoreus](https://github.com/aitananietoreus) | Informe: Material y métodos, Resultados y Conclusiones |

Repositorio creado a partir de la plantilla [seoanezonjic/project_template](https://github.com/seoanezonjic/project_template).
