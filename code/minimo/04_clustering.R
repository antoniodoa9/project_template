#!/usr/bin/env Rscript
# 04_clustering.R - Análisis mínimo, parte 5
# Agrupa los genes de la red en clusters con Louvain y guarda la tabla gen-cluster
# y la figura de la red coloreada por cluster.

suppressPackageStartupMessages({
  library(optparse)
  library(igraph)
})

opciones <- list(
  make_option(c("-i", "--red"), type = "character",
              default = "../results/minimo/tablas/red.tsv",
              help = "TSV de la red de 02_red.R [por defecto: %default]"),
  make_option(c("-o", "--salida"), type = "character",
              default = "../results/minimo/tablas/clusters.tsv",
              help = "TSV de salida con el cluster de cada gen [por defecto: %default]"),
  make_option(c("-f", "--figura"), type = "character",
              default = "../results/minimo/figuras/clusters_louvain.pdf",
              help = "PDF de salida con la red coloreada por cluster [por defecto: %default]")
)
opt <- parse_args(OptionParser(option_list = opciones))

dir.create(dirname(opt$salida), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$figura), recursive = TRUE, showWarnings = FALSE)

# Red no dirigida y sin genes aislados, como en el ejemplo del profesor
red <- read.delim(opt$red)
g <- graph_from_data_frame(red, directed = FALSE)

# Louvain, como en el ejemplo. El resultado puede cambiar entre ejecuciones,
# así que se fija la semilla
set.seed(123)
communities <- cluster_louvain(g)

# Tabla gen -> cluster (05_enriquecimiento.R la usa para el enriquecimiento de cada cluster)
clusters <- data.frame(gen = V(g)$name, cluster = as.vector(membership(communities)))
clusters <- clusters[order(clusters$cluster, clusters$gen), ]
write.table(clusters, opt$salida, sep = "\t", quote = FALSE, row.names = FALSE)

# Figura de la red coloreada por cluster, como en el ejemplo
# (la disposición de los nodos también tiene azar, por eso otra vez la semilla)
pdf(opt$figura)
set.seed(123)
plot(communities, g, vertex.size = 5, vertex.label.cex = 0.6,
     main = "Clustering con Louvain")
invisible(dev.off())

# Resumen para la salida de launch.sh
cat("Nodos de la red:", gorder(g), "\n")
cat("Número de clusters:", length(communities), "\n")
cat("Tamaños de los clusters:", sizes(communities), "\n")
cat("Modularidad:", round(modularity(communities), 3), "\n")
