#!/usr/bin/env Rscript
# A4_clusterings.R - Ampliación A4
# Compara Louvain (el del análisis mínimo) con otros algoritmos de clustering vistos en clase: modularidad, número de clusters y parecido entre particiones (NMI).
# También aplica los algoritmos a las redes de cada umbral de la A1 (tabla umbral × algoritmo).

suppressPackageStartupMessages({
  library(optparse)
  library(igraph)
})

opciones <- list(
  make_option(c("-i", "--red"), type = "character",
              default = "../results/minimo/tablas/red.tsv",
              help = "TSV de la red de 02_red.R [por defecto: %default]"),
  make_option(c("-o", "--salida"), type = "character",
              default = "../results/ampliaciones/A4_clusterings",
              help = "Carpeta de salida de las tablas y la figura [por defecto: %default]"),
  make_option(c("-r", "--redes"), type = "character",
              default = "../results/ampliaciones/A1_umbrales/redes",
              help = "Carpeta de la A1 con una subcarpeta umbral_<u> por umbral [por defecto: %default]"),
  make_option(c("-u", "--umbrales"), type = "character", default = "400,500,600,700,800,900",
              help = "Umbrales de STRING separados por comas, los mismos que en la A1 [por defecto: %default]")
)
opt <- parse_args(OptionParser(option_list = opciones))

dir.create(opt$salida, recursive = TRUE, showWarnings = FALSE)

# La misma red que en 04_clustering.R: no dirigida y sin genes aislados
red <- read.delim(opt$red)
g <- graph_from_data_frame(red, directed = FALSE)

# Los 7 algoritmos de compara_metodos() de Clustering (I) de Perkins. Louvain va el primero y justo después de set.seed(123), como en 04_clustering.R, para que salga
# la misma partición. Leiden, con la modularidad como objetivo
metodos <- list(
  "Louvain"           = cluster_louvain,
  "Leiden"            = function(g) cluster_leiden(g, objective_function = "modularity"),
  "Fast greedy"       = cluster_fast_greedy,
  "Walktrap"          = cluster_walktrap,
  "Infomap"           = cluster_infomap,
  "Label propagation" = cluster_label_prop,
  "Edge betweenness"  = cluster_edge_betweenness
)
set.seed(123)
resultados <- lapply(metodos, function(metodo) metodo(g))

# Tabla resumen, como compara_metodos() de Clustering (I). No hay comunidades reales, así que la NMI se calcula frente a Louvain (1 = misma partición)
resumen <- data.frame(
  metodo      = names(resultados),
  clusters    = sapply(resultados, length),
  modularidad = round(sapply(resultados, function(r) modularity(g, membership(r))), 3),
  NMI_louvain = round(sapply(resultados, compare, resultados[["Louvain"]], method = "nmi"), 2),
  tamanos     = sapply(resultados, function(r) paste(sort(sizes(r), decreasing = TRUE), collapse = ","))
)
write.table(resumen, file.path(opt$salida, "A4_resumen.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# NMI de cada pareja de algoritmos
nmi <- sapply(resultados, function(a) sapply(resultados, compare, a, method = "nmi"))
nmi <- round(nmi, 2)
write.table(data.frame(metodo = rownames(nmi), nmi, check.names = FALSE),
            file.path(opt$salida, "A4_nmi.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Cluster de cada gen con cada algoritmo, para ver qué genes cambian de módulo
particiones <- data.frame(gen = V(g)$name, sapply(resultados, membership), check.names = FALSE)
particiones <- particiones[order(particiones$Louvain, particiones$gen), ]
write.table(particiones, file.path(opt$salida, "A4_particiones.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Límite de resolución (Clustering I): la modularidad tiende a fusionar módulos
# pequeños. Louvain con varios valores de resolution (mayor = clusters más pequeños;
# 1 es el valor por defecto, el del análisis mínimo)
resoluciones <- c(0.5, 1, 1.5)
louvain_res <- lapply(resoluciones, function(r) {
  set.seed(123)
  cluster_louvain(g, resolution = r)
})
tabla_res <- data.frame(
  resolution  = resoluciones,
  clusters    = sapply(louvain_res, length),
  modularidad = round(sapply(louvain_res, function(r) modularity(g, membership(r))), 3),
  NMI_louvain = round(sapply(louvain_res, compare, resultados[["Louvain"]], method = "nmi"), 2),
  tamanos     = sapply(louvain_res, function(r) paste(sort(sizes(r), decreasing = TRUE), collapse = ","))
)
write.table(tabla_res, file.path(opt$salida, "A4_resolution.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Tabla umbral × algoritmo: los mismos 7 algoritmos en la red de cada umbral de
# la A1 (como en la A1 y en 04_clustering.R, sin los genes aislados). Sirve para
# ver si los módulos dependen del umbral elegido o del algoritmo
umbrales <- as.integer(strsplit(opt$umbrales, ",")[[1]])
tabla_umbral <- data.frame()
for (u in umbrales) {
  red_u <- read.delim(file.path(opt$redes, paste0("umbral_", u), "red.tsv"))
  g_u <- graph_from_data_frame(red_u, directed = FALSE)
  set.seed(123)
  res_u <- lapply(metodos, function(metodo) metodo(g_u))
  tabla_umbral <- rbind(tabla_umbral, data.frame(
    umbral      = u,
    metodo      = names(res_u),
    nodos       = gorder(g_u),
    clusters    = sapply(res_u, length),
    modularidad = round(sapply(res_u, function(r) modularity(g_u, membership(r))), 3)
  ))
}
write.table(tabla_umbral, file.path(opt$salida, "A4_umbral_algoritmo.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Una página por algoritmo, con la misma disposición de nodos para poder compararlas que luego no me entero
set.seed(123)
disposicion <- layout_with_fr(g)
pdf(file.path(opt$salida, "A4_clusterings.pdf"))
for (m in names(resultados)) {
  plot(resultados[[m]], g, layout = disposicion, vertex.size = 5, vertex.label.cex = 0.6,
       main = paste0(m, ": ", length(resultados[[m]]), " clusters, Q = ",
                     round(modularity(g, membership(resultados[[m]])), 3)))
}
invisible(dev.off())

# Resumen para la salida de launch.sh
cat("Nodos de la red:", gorder(g), "\n")
print(resumen, row.names = FALSE)
cat("\nNMI entre algoritmos:\n")
print(nmi)
cat("\nLouvain con distintos valores de resolution:\n")
print(tabla_res, row.names = FALSE)
cat("\nModularidad de cada algoritmo en cada umbral:\n")
print(xtabs(modularidad ~ umbral + metodo, data = tabla_umbral))
cat("\nNúmero de clusters de cada algoritmo en cada umbral:\n")
print(xtabs(clusters ~ umbral + metodo, data = tabla_umbral))
