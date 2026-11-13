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
- Conexión a internet: el análisis consulta la API de HPO y descarga los datos de STRING.
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
- `launch.sh` ejecuta los scripts en orden y escribe todos los resultados en `results/`. La primera vez descarga los datos de STRING (unos 25 MB) en `results/string_data/`.

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

Los genes que están en STRING pero no tienen ninguna interacción por encima del umbral se mantienen en la topología como nodos aislados, y `02_red.R` los muestra por pantalla. El clustering se hace solo con los genes conectados, porque cada gen aislado formaría un cluster de un solo gen.

### Ampliaciones

**A1. Barrido de umbrales de STRING** (`code/ampliaciones/A1_barrido_umbrales.R`). Responde a cómo cambia la red según la confianza que se exige a las interacciones. `launch.sh` construye la red con seis umbrales (400, 500, 600, 700, 800 y 900) usando el mismo `02_red.R`, y el script compara las seis redes con las mismas medidas que el análisis mínimo.

| Archivo | Contenido | Script |
|---|---|---|
| `ampliaciones/A1_umbrales/redes/umbral_<u>/` | La red (`red.tsv`) y los genes encontrados en STRING (`genes_string.tsv`) con cada umbral | `02_red.R`, llamado desde `launch.sh` |
| `ampliaciones/A1_umbrales/barrido_umbrales.tsv` | Por umbral: genes con y sin interacción, aristas, densidad, nodos de la componente mayor, clusters de Louvain, clusters de 5 genes o más y modularidad | `A1_barrido_umbrales.R` |
| `ampliaciones/A1_umbrales/genes_por_umbral.tsv` | Por gen: si tiene alguna interacción con cada umbral y el umbral más alto al que la tiene | `A1_barrido_umbrales.R` |
| `ampliaciones/A1_umbrales/barrido_umbrales.pdf` | Genes conectados, interacciones, clusters y modularidad frente al umbral | `A1_barrido_umbrales.R` |

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
| `there is no package called ...` al ejecutar un script suelto | Se ha ejecutado sin la carpeta de paquetes del proyecto | Ejecutar con `./code/launch.sh`, o antes hacer `export R_LIBS_USER=$(pwd)/../software/deps_R` desde `code/` |

## Autores

| Persona | GitHub | Parte del trabajo |
|---|---|---|
| Antonio Sánchez Doña (coordinador) | [@antoniodoa9](https://github.com/antoniodoa9) | Código: instalación y ejecución, genes, red y topología, y sus ampliaciones |
| Nicolás Sánchez Alonso | [@NicoSA03](https://github.com/NicoSA03) | Código: clustering y enriquecimiento, y sus ampliaciones |
| Eva Lucía López Porras | [@evalu5](https://github.com/evalu5) | Informe: Introducción y Discusión |
| Aitana Nieto Reus | [@aitananietoreus](https://github.com/aitananietoreus) | Informe: Material y métodos, Resultados y Conclusiones |

Repositorio creado a partir de la plantilla [seoanezonjic/project_template](https://github.com/seoanezonjic/project_template).
