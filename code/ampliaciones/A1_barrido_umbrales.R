#!/usr/bin/env Rscript
# A1_barrido_umbrales.R - Ampliación 1
# Compara las redes de STRING obtenidas con varios umbrales de puntuación:
# tamaño de la red, genes sin interacción, clusters de Louvain y a qué umbral
# deja de tener interacciones cada gen.

suppressPackageStartupMessages({
  library(optparse)
  library(igraph)
})

opciones <- list(
  make_option(c("-u", "--umbrales"), type = "character", default = "400,500,600,700,800,900",
              help = "Umbrales de STRING separados por comas [por defecto: %default]"),
  make_option(c("-r", "--redes"), type = "character",
              default = "../results/ampliaciones/A1_umbrales/redes",
              help = "Carpeta con una subcarpeta umbral_<u> por umbral, hecha con 02_red.R [por defecto: %default]"),
  make_option(c("-o", "--salida"), type = "character",
              default = "../results/ampliaciones/A1_umbrales",
              help = "Carpeta de salida [por defecto: %default]")
)
opt <- parse_args(OptionParser(option_list = opciones))

umbrales <- as.integer(strsplit(opt$umbrales, ",")[[1]])
dir.create(opt$salida, recursive = TRUE, showWarnings = FALSE)

barrido <- data.frame()   # una fila por umbral
conectados <- list()      # genes con alguna interacción en cada umbral

for (u in umbrales) {
  carpeta <- file.path(opt$redes, paste0("umbral_", u))
  red <- read.delim(file.path(carpeta, "red.tsv"))
  genes <- read.delim(file.path(carpeta, "genes_string.tsv"))

  # Igual que en 03_topologia.R: todos los genes de STRING, con los aislados
  g <- graph_from_data_frame(red, directed = FALSE, vertices = genes)
  grado <- degree(g)

  # Igual que en 04_clustering.R: Louvain sin los aislados y con la misma semilla
  g_conectada <- graph_from_data_frame(red, directed = FALSE)
  set.seed(123)
  louvain <- cluster_louvain(g_conectada)

  barrido <- rbind(barrido, data.frame(
    umbral = u,
    genes_conectados = sum(grado > 0),
    genes_sin_interaccion = sum(grado == 0),
    aristas = gsize(g),
    densidad = round(edge_density(g), 4),
    nodos_componente_mayor = max(components(g)$csize),
    clusters = length(louvain),
    clusters_5_o_mas = sum(sizes(louvain) >= 5),
    modularidad = round(modularity(louvain), 3)))

  conectados[[as.character(u)]] <- V(g)$name[grado > 0]
}

write.table(barrido, file.path(opt$salida, "barrido_umbrales.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Tabla por gen: si tiene alguna interacción en cada umbral ("si" o "no") y el
# umbral más alto al que todavía la tiene (NA si no tiene ninguna en ninguno)
todos <- sort(genes$name)
genes_umbral <- data.frame(gen = todos)
for (u in umbrales) {
  genes_umbral[[paste0("umbral_", u)]] <- ifelse(todos %in% conectados[[as.character(u)]], "si", "no")
}
genes_umbral$umbral_maximo <- sapply(todos, function(gen) {
  u_con <- umbrales[sapply(umbrales, function(u) gen %in% conectados[[as.character(u)]])]
  if (length(u_con) == 0) NA else max(u_con)
})
# Primero los genes que antes se quedan sin interacciones
genes_umbral <- genes_umbral[order(genes_umbral$umbral_maximo, genes_umbral$gen, na.last = FALSE), ]
write.table(genes_umbral, file.path(opt$salida, "genes_por_umbral.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Figura: cuatro paneles frente al umbral; la línea discontinua marca el 700
# que usa el análisis mínimo
pdf(file.path(opt$salida, "barrido_umbrales.pdf"), width = 8, height = 7)
par(mfrow = c(2, 2), mar = c(5, 4.5, 3, 1))
plot(barrido$umbral, barrido$genes_conectados, type = "b", pch = 19,
     xlab = "Umbral de STRING", ylab = "Genes con interacción", main = "Genes conectados")
abline(v = 700, lty = 2, col = "red")
plot(barrido$umbral, barrido$aristas, type = "b", pch = 19,
     xlab = "Umbral de STRING", ylab = "Interacciones", main = "Interacciones")
abline(v = 700, lty = 2, col = "red")
plot(barrido$umbral, barrido$clusters, type = "b", pch = 19, ylim = c(0, max(barrido$clusters)),
     xlab = "Umbral de STRING", ylab = "Número de clusters", main = "Clusters de Louvain")
lines(barrido$umbral, barrido$clusters_5_o_mas, type = "b", pch = 17, col = "blue")
legend("bottomleft", legend = c("Todos", "De 5 genes o más"),
       pch = c(19, 17), col = c("black", "blue"), bty = "n")
abline(v = 700, lty = 2, col = "red")
plot(barrido$umbral, barrido$modularidad, type = "b", pch = 19,
     xlab = "Umbral de STRING", ylab = "Modularidad (Q)", main = "Modularidad de Louvain")
abline(v = 700, lty = 2, col = "red")
invisible(dev.off())

# Resumen para la salida de launch.sh
print(barrido, row.names = FALSE)
