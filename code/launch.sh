#!/usr/bin/env bash
# Ejecuta todo el análisis. Todo lo que genera va a results/
set -e
cd "$(dirname "$0")"
export R_LIBS_USER="$(pwd)/../software/deps_R"
export TMPDIR="$(pwd)/../software/tmp"
# Hasta 10 minutos por descarga (R espera 60 segundos por defecto y no basta con una conexión lenta)
export R_DEFAULT_INTERNET_TIMEOUT=600
mkdir -p ../results/minimo/tablas ../results/minimo/figuras ../results/ampliaciones

# ===== ANÁLISIS MÍNIMO (flujo del ejemplo del profesor) =====
printf "\n===== ANÁLISIS MÍNIMO =====\n"
printf "\n== 01_genes.R: genes del término HPO ==\n"
Rscript minimo/01_genes.R -p HP:0002647 -o ../results/minimo/tablas/genes_hpo.tsv
printf "\n== 02_red.R: red de STRING ==\n"
Rscript minimo/02_red.R -i ../results/minimo/tablas/genes_hpo.tsv -s 700 -d ../results/string_data -o ../results/minimo/tablas
printf "\n== 03_topologia.R: propiedades de la red ==\n"
Rscript minimo/03_topologia.R -i ../results/minimo/tablas/red.tsv -g ../results/minimo/tablas/genes_string.tsv -o ../results/minimo/tablas -f ../results/minimo/figuras
printf "\n== 04_clustering.R: módulos con Louvain ==\n"
Rscript minimo/04_clustering.R -i ../results/minimo/tablas/red.tsv -o ../results/minimo/tablas/clusters.tsv -f ../results/minimo/figuras/clusters_louvain.pdf
printf "\n== 05_enriquecimiento.R: enriquecimiento funcional GO ==\n"
Rscript minimo/05_enriquecimiento.R -g ../results/minimo/tablas/genes_hpo.tsv -c ../results/minimo/tablas/clusters.tsv -t ../results/minimo/tablas -f ../results/minimo/figuras

# ===== AMPLIACIONES =====
printf "\n===== AMPLIACIONES =====\n"
# A1. Barrido de umbrales de STRING: una red por umbral con 02_red.R y su comparación
printf "\n== A1. Barrido de umbrales de STRING ==\n"
for umbral in 400 500 600 700 800 900; do
  Rscript minimo/02_red.R -i ../results/minimo/tablas/genes_hpo.tsv -s $umbral -d ../results/string_data -o ../results/ampliaciones/A1_umbrales/redes/umbral_$umbral
done
Rscript ampliaciones/A1_barrido_umbrales.R -r ../results/ampliaciones/A1_umbrales/redes -o ../results/ampliaciones/A1_umbrales

# A2. Modularidad y otras propiedades de la red frente a 1000 redes aleatorias de dos modelos (Erdos-Renyi y grado conservado)
printf "\n== A2. Red real frente a redes aleatorias ==\n"
Rscript ampliaciones/A2_redes_aleatorias.R -i ../results/minimo/tablas/red.tsv -n 1000 -o ../results/ampliaciones/A2_redes_aleatorias

# A3. Robustez: ataque dirigido a los hubs y a los cuellos de botella frente a fallo al azar
printf "\n== A3. Ataque a los hubs frente a fallo al azar ==\n"
Rscript ampliaciones/A3_ataque_hubs.R -i ../results/minimo/tablas/red.tsv -c ../results/minimo/tablas/clusters.tsv -n 1000 -o ../results/ampliaciones/A3_ataque_hubs

# A4. Otros algoritmos de clustering comparados con Louvain (modularidad y NMI)
printf "\n== A4. Otros algoritmos de clustering ==\n"
Rscript ampliaciones/A4_clusterings.R -i ../results/minimo/tablas/red.tsv -o ../results/ampliaciones/A4_clusterings

# A5. Enriquecimiento ampliado: GO BP, MF y CC, simplify() y figuras de enrichplot
printf "\n== A5. Enriquecimiento ampliado ==\n"
Rscript ampliaciones/A5_enriquecimiento.R -g ../results/minimo/tablas/genes_hpo.tsv -o ../results/ampliaciones/A5_enriquecimiento

# A6. Genes candidatos: random walk con reinicio desde los genes de HPO sobre la red ponderada de STRING con tres umbrales
printf "\n== A6. Genes candidatos (random walk con reinicio) ==\n"
Rscript ampliaciones/A6_rwr_roc.R -g ../results/minimo/tablas/genes_hpo.tsv -u 500,700,900 -d ../results/string_data -o ../results/ampliaciones/A6_rwr_roc

# A8. Embedding de la subred de 700 de la A6 (paseos aleatorios, PMI y SVD), usado para predecir interacciones y semillas escondidas
printf "\n== A8. Embedding de la subred ==\n"
../software/venv_py/bin/python ampliaciones/A8_embedding.py -i ../results/ampliaciones/A6_rwr_roc/A6_subred_700.tsv -g ../results/minimo/tablas/genes_hpo.tsv -a ../results/ampliaciones/A6_rwr_roc -o ../results/ampliaciones/A8_embedding

# A9. Fármacos de ChEMBL para las semillas y los candidatos de la A6 y de la A8, comparados con conjuntos al azar de la subred de 700
printf "\n== A9. Fármacos de ChEMBL ==\n"
Rscript ampliaciones/A9_chembl.R -i ../results/ampliaciones/A6_rwr_roc/A6_subred_700.tsv -g ../results/minimo/tablas/genes_hpo.tsv --a6 ../results/ampliaciones/A6_rwr_roc/A6_candidatos.tsv --a8 ../results/ampliaciones/A8_embedding/A8_candidatos.tsv -o ../results/ampliaciones/A9_chembl

# Versiones de R, de Python y de los paquetes usados, para poder citarlas en Métodos
printf "\n== Versiones de R, de Python y de los paquetes ==\n"
Rscript -e "suppressPackageStartupMessages({library(igraph); library(optparse); library(httr); library(jsonlite); library(STRINGdb); library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot); library(Matrix)}); sessionInfo()" > ../results/session_info.txt
# Y las de Python y sus paquetes (A8)
{ printf "\nPython de la A8: "; ../software/venv_py/bin/python --version; ../software/venv_py/bin/pip freeze; } >> ../results/session_info.txt
echo "Guardadas en results/session_info.txt"
printf "\nlaunch.sh terminado\n"
