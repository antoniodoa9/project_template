#!/usr/bin/env bash
# Instala todos los paquetes de R y de Python dentro del repositorio (software/), sin sudo
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
# Paquetes de Python de la A8 en un entorno virtual dentro del repositorio (software/venv_py), sin sudo.
# El entorno se crea sin pip, porque en Ubuntu crearlo con pip exige el paquete python3-venv (sudo apt), y pip se instala con get-pip.py, el instalador oficial de pip
if [ ! -x "$SOFT/venv_py/bin/pip" ]; then
  python3 -m venv --without-pip "$SOFT/venv_py"
  curl -fsSL https://bootstrap.pypa.io/get-pip.py -o "$SOFT/tmp/get-pip.py"
  "$SOFT/venv_py/bin/python" "$SOFT/tmp/get-pip.py" --no-cache-dir --disable-pip-version-check
fi
# Versiones fijas, las probadas; --no-cache-dir para no escribir nada fuera del repositorio. Si ya están, pip no hace nada
"$SOFT/venv_py/bin/pip" install --no-cache-dir --disable-pip-version-check numpy==2.5.5 scipy==1.18.1 scikit-learn==1.9.1 igraph==1.0.0 matplotlib==3.11.2
echo "setup.sh terminado"
