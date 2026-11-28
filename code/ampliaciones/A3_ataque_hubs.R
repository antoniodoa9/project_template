#!/usr/bin/env Rscript
# A3_ataque_hubs.R - Ampliación 3
# Robustez de la red: quita genes de la red del mínimo de dos formas
# (Perkins, Topología I, ejercicio "¿Qué pasa si atacamos los hubs?"):
# - Ataque dirigido: primero los genes con más interacciones (grado) o los
#   que más caminos cortos atraviesan (intermediación)
# - Fallo al azar: en un orden aleatorio, repetido muchas veces
# y mide cuánto se rompe la red (tamaño de la componente mayor).
# Además clasifica cada gen como hub (mucho grado) o cuello de botella
# (mucha intermediación), con su módulo de Louvain.

suppressPackageStartupMessages({
  library(optparse)
  library(igraph)
})

opciones <- list(
  make_option(c("-i", "--red"), type = "character",
              default = "../results/minimo/tablas/red.tsv",
              help = "TSV de la red de 02_red.R [por defecto: %default]"),
  make_option(c("-c", "--clusters"), type = "character",
              default = "../results/minimo/tablas/clusters.tsv",
              help = "TSV de los módulos de 04_clustering.R [por defecto: %default]"),
  make_option(c("-n", "--repeticiones"), type = "integer", default = 1000,
              help = "Órdenes al azar para el fallo aleatorio [por defecto: %default]"),
  make_option(c("-p", "--porcentaje"), type = "double", default = 20,
              help = "Porcentaje de genes que se quitan en el ejercicio de clase [por defecto: %default]"),
  make_option(c("-o", "--salida"), type = "character",
              default = "../results/ampliaciones/A3_ataque_hubs",
              help = "Carpeta de salida [por defecto: %default]")
)
opt <- parse_args(OptionParser(option_list = opciones))
dir.create(opt$salida, recursive = TRUE, showWarnings = FALSE)

# La misma red que usan 04_clustering.R, A2 y A4: no dirigida y sin los genes
# aislados (un gen sin interacciones no sostiene ninguna conexión, así que
# quitarlo no cambia nada)
red <- read.delim(opt$red)
g <- graph_from_data_frame(red, directed = FALSE)
clusters <- read.delim(opt$clusters)
genes <- V(g)$name
n <- gorder(g)

# ===== 1. Hubs y cuellos de botella =====
# Grado: número de interacciones del gen (degree, Topología I)
# Intermediación: número de caminos más cortos entre otros dos genes que
# pasan por él (betweenness, Topología I); normalizada, entre 0 y 1
grado <- degree(g)
intermediacion <- betweenness(g, directed = FALSE)
intermediacion_norm <- betweenness(g, directed = FALSE, normalized = TRUE)

# Hub = grado en el 20 % más alto; cuello de botella = intermediación en el
# 20 % más alto (el mismo corte para las dos, como Yu et al. 2007)
corte <- 1 - opt$porcentaje / 100
es_hub <- grado >= quantile(grado, corte)
es_cuello <- intermediacion >= quantile(intermediacion, corte)
tipo <- ifelse(es_hub & es_cuello, "hub y cuello de botella",
        ifelse(es_hub, "solo hub",
        ifelse(es_cuello, "solo cuello de botella", "ninguno")))

# Módulo de Louvain de cada gen (de 04_clustering.R) y cuántos módulos
# distintos hay entre sus vecinos: un gen que une módulos tiene más de 1
modulo <- clusters$cluster[match(genes, clusters$gen)]
names(modulo) <- genes
modulos_vecinos <- sapply(genes, function(x)
  length(unique(modulo[names(neighbors(g, x))])))

genes_tabla <- data.frame(gen = genes,
                          grado = grado,
                          intermediacion = round(intermediacion, 1),
                          intermediacion_normalizada = round(intermediacion_norm, 3),
                          modulo = modulo,
                          modulos_vecinos = modulos_vecinos,
                          tipo = tipo)
genes_tabla <- genes_tabla[order(-genes_tabla$intermediacion, -genes_tabla$grado,
                                 genes_tabla$gen), ]
write.table(genes_tabla, file.path(opt$salida, "A3_hubs_cuellos.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ===== 2. Ataque dirigido frente a fallo al azar =====
# Tamaño de la componente mayor y número de componentes tras quitar unos genes
tras_quitar <- function(quitar) {
  h <- delete_vertices(g, quitar)
  if (gorder(h) == 0) return(c(mayor = 0, componentes = 0))
  comp <- components(h)
  c(mayor = max(comp$csize), componentes = comp$no)
}
# Quita los genes de uno en uno en el orden dado, de 0 a todos
curva <- function(orden) {
  t(sapply(0:n, function(k) tras_quitar(orden[seq_len(k)])))
}

# Órdenes de ataque: de más a menos grado o intermediación (los empates, por
# nombre del gen, para que el resultado no cambie entre ejecuciones)
orden_grado <- genes[order(-grado, genes)]
orden_intermediacion <- genes[order(-intermediacion, genes)]
ataque_grado <- curva(orden_grado)
ataque_intermediacion <- curva(orden_intermediacion)

# Fallo al azar: el mismo proceso con órdenes aleatorios
set.seed(123)
azar <- sapply(seq_len(opt$repeticiones), function(r) curva(sample(genes))[, "mayor"])

mayor_inicial <- ataque_grado[1, "mayor"]
curvas <- data.frame(
  genes_quitados = 0:n,
  fraccion_quitada = round((0:n) / n, 3),
  mayor_ataque_grado = ataque_grado[, "mayor"],
  mayor_ataque_intermediacion = ataque_intermediacion[, "mayor"],
  mayor_azar_media = round(rowMeans(azar), 2),
  mayor_azar_p2.5 = apply(azar, 1, quantile, 0.025),
  mayor_azar_p97.5 = apply(azar, 1, quantile, 0.975),
  componentes_ataque_grado = ataque_grado[, "componentes"],
  componentes_ataque_intermediacion = ataque_intermediacion[, "componentes"])
write.table(curvas, file.path(opt$salida, "A3_curvas.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# El ejercicio de clase: quitar el 20 % de los genes (6 de 31) y comparar.
# p empírico, como en la A2: (1 + órdenes al azar que dejan una componente
# mayor igual o más pequeña que el ataque) / (1 + repeticiones)
k <- round(n * opt$porcentaje / 100)
# Genes que hay que quitar para que la componente mayor baje a la mitad o menos
# (resume la curva entera, no solo el punto del 20 %)
hasta_la_mitad <- function(mayor) min(which(mayor <= mayor_inicial / 2)) - 1
azar_k <- azar[k + 1, ]
azar_mitad <- apply(azar, 2, hasta_la_mitad)
mitad_grado <- hasta_la_mitad(ataque_grado[, "mayor"])
mitad_intermediacion <- hasta_la_mitad(ataque_intermediacion[, "mayor"])
resumen <- data.frame(
  estrategia = c("ataque por grado", "ataque por intermediacion", "fallo al azar"),
  genes_quitados = k,
  quitados = c(paste(orden_grado[1:k], collapse = ","),
               paste(orden_intermediacion[1:k], collapse = ","),
               paste0("al azar (", opt$repeticiones, " repeticiones)")),
  componente_mayor = c(ataque_grado[k + 1, "mayor"],
                       ataque_intermediacion[k + 1, "mayor"],
                       round(mean(azar_k), 2)),
  componentes = c(ataque_grado[k + 1, "componentes"],
                  ataque_intermediacion[k + 1, "componentes"], NA),
  p_empirico = c(signif((1 + sum(azar_k <= ataque_grado[k + 1, "mayor"])) /
                          (1 + opt$repeticiones), 2),
                 signif((1 + sum(azar_k <= ataque_intermediacion[k + 1, "mayor"])) /
                          (1 + opt$repeticiones), 2), NA),
  genes_para_la_mitad = c(mitad_grado, mitad_intermediacion,
                          round(mean(azar_mitad), 2)),
  # p empírico de la mitad: órdenes al azar que la parten con tan pocos genes
  # como el ataque o menos
  p_empirico_mitad = c(signif((1 + sum(azar_mitad <= mitad_grado)) /
                                (1 + opt$repeticiones), 2),
                       signif((1 + sum(azar_mitad <= mitad_intermediacion)) /
                                (1 + opt$repeticiones), 2), NA))
write.table(resumen, file.path(opt$salida, "A3_resumen.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Conectividad de la componente mayor (Topología I): cuántos genes
# (vertex_connectivity) o interacciones (edge_connectivity) como mínimo hay
# que quitar para partirla
mayor_comp <- largest_component(g)
conectividad_genes <- vertex_connectivity(mayor_comp)
conectividad_aristas <- edge_connectivity(mayor_comp)

# ===== 3. Figuras =====
pdf(file.path(opt$salida, "A3_ataque_hubs.pdf"), width = 7, height = 6)
# Curva de robustez: componente mayor (relativa a la inicial) según la
# fracción de genes quitados; banda gris = 95 % de los órdenes al azar
x <- curvas$fraccion_quitada
plot(x, curvas$mayor_ataque_grado / mayor_inicial, type = "n", ylim = c(0, 1),
     xlab = "Fracción de genes quitados",
     ylab = "Componente mayor (fracción de la inicial)",
     main = "Ataque dirigido frente a fallo al azar")
polygon(c(x, rev(x)),
        c(curvas$mayor_azar_p2.5, rev(curvas$mayor_azar_p97.5)) / mayor_inicial,
        col = "grey85", border = NA)
lines(x, curvas$mayor_azar_media / mayor_inicial, col = "black", lwd = 2)
lines(x, curvas$mayor_ataque_grado / mayor_inicial, col = "red", lwd = 2, type = "o", pch = 16, cex = 0.6)
lines(x, curvas$mayor_ataque_intermediacion / mayor_inicial, col = "blue", lwd = 2, type = "o", pch = 17, cex = 0.6)
abline(v = k / n, lty = 2)
legend("topright", bty = "n",
       legend = c("Ataque por grado", "Ataque por intermediación",
                  paste0("Fallo al azar (media de ", opt$repeticiones, ")"),
                  "95 % de los órdenes al azar",
                  paste0(opt$porcentaje, " % de los genes")),
       col = c("red", "blue", "black", "grey85", "black"),
       lwd = c(2, 2, 2, 10, 1), lty = c(1, 1, 1, 1, 2))

# Hubs frente a cuellos de botella: grado frente a intermediación, con los
# dos cortes del 20 %
colores <- c("hub y cuello de botella" = "purple", "solo hub" = "red",
             "solo cuello de botella" = "blue", "ninguno" = "grey50")
plot(genes_tabla$grado, genes_tabla$intermediacion,
     col = colores[genes_tabla$tipo], pch = 16,
     xlab = "Grado", ylab = "Intermediación",
     xlim = c(min(genes_tabla$grado), max(genes_tabla$grado) + 1),   # margen para que no se corten los nombres
     ylim = c(0, 1.1 * max(genes_tabla$intermediacion)),
     main = "Hubs y cuellos de botella")
# Etiquetas de los genes destacados: encima del punto, y debajo si otro gen
# destacado con el mismo grado está casi en el mismo sitio (para que no se pisen)
destacados <- genes_tabla[genes_tabla$tipo != "ninguno", ]
juntos <- duplicated(paste(destacados$grado,
                           round(destacados$intermediacion / 5)))
text(destacados$grado, destacados$intermediacion, destacados$gen,
     pos = ifelse(juntos, 1, 3), cex = 0.9)
abline(v = quantile(grado, corte), h = quantile(intermediacion, corte), lty = 2)
legend("topleft", bty = "n", legend = names(colores), col = colores, pch = 16)
invisible(dev.off())

# Resumen para la salida de launch.sh
cat("Red:", n, "genes y", gsize(g), "interacciones; componente mayor:",
    mayor_inicial, "genes\n")
cat("Para partir la componente mayor hay que quitar como mínimo",
    conectividad_genes, "gen(es) o", conectividad_aristas, "interacción(es)\n")
cat("Hubs y cuellos de botella (", opt$porcentaje, "% más alto):\n")
print(genes_tabla[genes_tabla$tipo != "ninguno", ], row.names = FALSE)
cat("\nQuitando", k, "de", n, "genes:\n")
print(resumen, row.names = FALSE)
