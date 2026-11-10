#!/usr/bin/env bash
# Ejecuta todo el análisis. Todo lo que genera va a results/
set -e
cd "$(dirname "$0")"
export R_LIBS_USER="$(pwd)/../software/deps_R"
export TMPDIR="$(pwd)/../software/tmp"
mkdir -p ../results/minimo/tablas ../results/minimo/figuras ../results/ampliaciones

# ===== ANÁLISIS MÍNIMO (flujo del ejemplo del profesor) =====
Rscript minimo/01_genes.R -p HP:0002647 -o ../results/minimo/tablas/genes_hpo.tsv
Rscript minimo/02_red.R -i ../results/minimo/tablas/genes_hpo.tsv -s 700 -d ../results/string_data -o ../results/minimo/tablas
Rscript minimo/03_topologia.R -i ../results/minimo/tablas/red.tsv -g ../results/minimo/tablas/genes_string.tsv -o ../results/minimo/tablas -f ../results/minimo/figuras
Rscript minimo/04_clustering.R -i ../results/minimo/tablas/red.tsv -o ../results/minimo/tablas/clusters.tsv -f ../results/minimo/figuras/clusters_louvain.pdf

# ===== AMPLIACIONES =====
# A1. Barrido de umbrales de STRING: una red por umbral con 02_red.R y su comparación
for umbral in 400 500 600 700 800 900; do
  Rscript minimo/02_red.R -i ../results/minimo/tablas/genes_hpo.tsv -s $umbral -d ../results/string_data -o ../results/ampliaciones/A1_umbrales/redes/umbral_$umbral
done
Rscript ampliaciones/A1_barrido_umbrales.R -r ../results/ampliaciones/A1_umbrales/redes -o ../results/ampliaciones/A1_umbrales
