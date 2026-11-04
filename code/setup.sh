#!/usr/bin/env bash
# Instala todos los paquetes de R dentro del repositorio (software/), sin sudo
set -e
cd "$(dirname "$0")"
SOFT="$(pwd)/../software"
mkdir -p "$SOFT/deps_R" "$SOFT/tmp"
export TMPDIR="$SOFT/tmp"
export R_LIBS_USER="$SOFT/deps_R"
Rscript -e "
options(timeout = 1200)
lib  <- Sys.getenv('R_LIBS_USER')
cran <- c('BiocManager', 'igraph', 'optparse', 'httr', 'jsonlite')
bioc <- c('STRINGdb', 'clusterProfiler', 'org.Hs.eg.db', 'enrichplot')
falta <- function(p) p[!p %in% rownames(installed.packages(lib.loc = lib))]
if (length(falta(cran)) > 0) install.packages(falta(cran), lib = lib, repos = 'https://cloud.r-project.org', Ncpus = 4)
if (length(falta(bioc)) > 0) BiocManager::install(falta(bioc), lib = lib, ask = FALSE, update = FALSE, Ncpus = 4)
if (length(falta(c(cran, bioc))) > 0) stop('Faltan paquetes: ', paste(falta(c(cran, bioc)), collapse = ', '))
"
echo "setup.sh terminado"
