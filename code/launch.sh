#!/usr/bin/env bash
# Ejecuta todo el análisis. Todo lo que genera va a results/
set -e
cd "$(dirname "$0")"
export R_LIBS_USER="$(pwd)/../software/deps_R"
export TMPDIR="$(pwd)/../software/tmp"
mkdir -p ../results/minimo ../results/ampliaciones

# ===== ANÁLISIS MÍNIMO (flujo del ejemplo del profesor) =====
Rscript 01_genes.R -p HP:0002647 -o ../results/minimo/genes_hpo.tsv

# ===== AMPLIACIONES =====
# (se irán añadiendo aquí)
