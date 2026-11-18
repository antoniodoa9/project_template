#!/usr/bin/env Rscript
# A2_redes_aleatorias.R - Ampliación 2
# Compara la red del mínimo con redes aleatorias del mismo tamaño para saber si
# su modularidad y sus otras propiedades son mayores de lo que daría el azar
# (Perkins, Topología II: "Un valor de Q por sí solo no dice mucho: hay que
# compararlo con el que se obtiene en redes aleatorias del mismo tamaño y densidad").

suppressPackageStartupMessages({
  library(optparse)
  library(igraph)
})

opciones <- list(
  make_option(c("-i", "--red"), type = "character",
              default = "../results/minimo/tablas/red.tsv",
              help = "TSV de la red de 02_red.R [por defecto: %default]"),
  make_option(c("-n", "--repeticiones"), type = "integer", default = 1000,
              help = "Redes aleatorias de cada modelo [por defecto: %default]"),
  make_option(c("-o", "--salida"), type = "character",
              default = "../results/ampliaciones/A2_redes_aleatorias",
              help = "Carpeta de salida [por defecto: %default]")
)
opt <- parse_args(OptionParser(option_list = opciones))
dir.create(opt$salida, recursive = TRUE, showWarnings = FALSE)

# La misma red que usa 04_clustering.R: no dirigida y sin los genes aislados
red <- read.delim(opt$red)
g <- graph_from_data_frame(red, directed = FALSE)

# Propiedades que se comparan (todas de Topología I y II de Perkins):
# modularidad de Louvain, transitividad global (cuántos triángulos se cierran),
# distancia media entre genes y asortatividad por grado (si los genes muy
# conectados se unen entre sí)
propiedades <- function(red) {
  c(modularidad   = modularity(cluster_louvain(red)),
    transitividad = transitivity(red, type = "global"),
    distancia     = mean_distance(red),
    asortatividad = assortativity_degree(red, directed = FALSE))
}

# Valores reales, con la misma semilla que 04_clustering.R para que la
# modularidad sea la misma (0.432)
set.seed(123)
real <- propiedades(g)

# Dos modelos de red aleatoria, los dos de clase:
# - Erdos-Renyi (sample_gnm): mismos genes y mismas interacciones, puestas al
#   azar. Es el que usa Perkins en Topología II
# - Rewire conservando el grado: cada gen mantiene su número de interacciones
#   y solo cambia con quién. Responde a si la estructura se debe solo a que
#   hay unos pocos genes muy conectados
n <- gorder(g)
m <- gsize(g)
set.seed(123)
nulos <- data.frame()
for (i in seq_len(opt$repeticiones)) {
  er <- sample_gnm(n, m)
  gr <- rewire(g, with = keeping_degseq(niter = 10 * m))
  nulos <- rbind(nulos,
                 data.frame(modelo = "Erdos-Renyi", repeticion = i, t(propiedades(er))),
                 data.frame(modelo = "Grado conservado", repeticion = i, t(propiedades(gr))))
}
write.table(nulos, file.path(opt$salida, "A2_nulos.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Para cada propiedad y modelo: media y desviación de las redes aleatorias,
# z = (real - media) / desviación, y p empírico. El p empírico cuenta cuántas
# redes aleatorias llegan al valor real, en la dirección en que se aleja:
# (1 + número de redes con un valor igual o más extremo) / (1 + repeticiones)
resumen <- data.frame()
for (modelo in unique(nulos$modelo)) {
  for (p in names(real)) {
    x <- nulos[nulos$modelo == modelo, p]
    x <- x[is.finite(x)]
    extremos <- if (real[p] >= mean(x)) sum(x >= real[p]) else sum(x <= real[p])
    resumen <- rbind(resumen, data.frame(
      propiedad = p, modelo = modelo,
      real = round(real[p], 3),
      media_azar = round(mean(x), 3),
      desviacion_azar = round(sd(x), 3),
      z = round((real[p] - mean(x)) / sd(x), 2),
      p_empirico = signif((1 + extremos) / (1 + length(x)), 2),
      redes_validas = length(x)))
  }
}
write.table(resumen, file.path(opt$salida, "A2_resumen.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Figura: un histograma por propiedad y modelo, con el valor real en rojo,
# como el de la modularidad en Topología II
pdf(file.path(opt$salida, "A2_redes_aleatorias.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
for (p in names(real)) {
  for (modelo in unique(nulos$modelo)) {
    x <- nulos[nulos$modelo == modelo, p]
    x <- x[is.finite(x)]
    hist(x, breaks = 30, col = "grey", border = "white",
         xlim = range(c(x, real[p])),
         main = paste0(p, " (", modelo, ")"),
         xlab = paste(p, "en", length(x), "redes aleatorias"))
    abline(v = real[p], col = "red", lwd = 2)
  }
}
invisible(dev.off())

# Resumen para la salida de launch.sh
cat("Red:", n, "genes y", m, "interacciones;", opt$repeticiones, "redes aleatorias de cada modelo\n")
print(resumen, row.names = FALSE)
