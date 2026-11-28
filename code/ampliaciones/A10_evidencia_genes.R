#!/usr/bin/env Rscript
# A10_evidencia_genes.R - Ampliación 10
# No todos los genes del término están igual de relacionados con él. Un gen entra en el término porque alguna
# enfermedad causada por ese gen tiene la disección aórtica entre sus fenotipos, y HPO indica con qué frecuencia
# la presentan los pacientes de esa enfermedad. Este script mide esa evidencia para cada gen, construye una red
# previa solo con los genes de evidencia fuerte y la compara con la red completa (módulos y modularidad), y cruza
# la evidencia con el AUC de cada semilla en la A6 (¿la red recupera peor los genes de evidencia débil?).

suppressPackageStartupMessages({
  library(optparse)
  library(igraph)
})

# genes_to_phenotype.txt y hp.obo son grandes: hasta 10 minutos por descarga
options(timeout = 600)

opciones <- list(
  make_option(c("-g", "--genes"), type = "character",
              default = "../results/minimo/tablas/genes_hpo.tsv",
              help = "TSV de genes de 01_genes.R [por defecto: %default]"),
  make_option(c("-i", "--red"), type = "character",
              default = "../results/minimo/tablas/red.tsv",
              help = "TSV de la red de 02_red.R [por defecto: %default]"),
  make_option(c("-c", "--clusters"), type = "character",
              default = "../results/minimo/tablas/clusters.tsv",
              help = "TSV de módulos de 04_clustering.R [por defecto: %default]"),
  make_option(c("-p", "--hpo"), type = "character", default = "HP:0002647",
              help = "Término HPO de partida [por defecto: %default]"),
  make_option(c("-a", "--auc"), type = "character",
              default = "../results/ampliaciones/A6_rwr_roc/A6_validacion.tsv",
              help = "TSV con el AUC de cada semilla de la A6 (si no existe, se omite el cruce) [por defecto: %default]"),
  make_option(c("-f", "--fuerte"), type = "double", default = 0.30,
              help = "Frecuencia mínima para considerar fuerte la evidencia de un gen [por defecto: %default]"),
  make_option(c("-o", "--salida"), type = "character",
              default = "../results/ampliaciones/A10_evidencia_genes",
              help = "Carpeta de salida [por defecto: %default]")
)
opt <- parse_args(OptionParser(option_list = opciones))
dir.create(opt$salida, recursive = TRUE, showWarnings = FALSE)

# 1. Anotaciones gen-fenotipo de HPO y jerarquía de la ontología (los mismos archivos que la A7)
url_hpo <- "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download"
g2p <- read.delim(file.path(url_hpo, "genes_to_phenotype.txt"))
obo <- readLines(file.path(url_hpo, "hp.obo"))
version_hpo <- sub("^data-version: ", "", obo[grepl("^data-version:", obo)][1])
genes <- sort(unique(read.delim(opt$genes)$name))

# El término y todos sus descendientes (por ejemplo "Ascending aortic dissection"), bajando por los is_a de hp.obo:
# la API de HPO devuelve también los genes anotados a los descendientes
id_actual <- cumsum(grepl("^id: HP:", obo))
ids <- sub("^id: ", "", obo[grepl("^id: HP:", obo)])
lineas_is_a <- grepl("^is_a: HP:", obo) & id_actual > 0
is_a <- data.frame(hijo = ids[id_actual[lineas_is_a]],
                   padre = sub("^is_a: (HP:[0-9]+).*", "\\1", obo[lineas_is_a]))
terminos <- opt$hpo
repeat {
  nuevos <- setdiff(is_a$hijo[is_a$padre %in% terminos], terminos)
  if (length(nuevos) == 0) break
  terminos <- c(terminos, nuevos)
}

# 2. Frecuencia de la disección en las enfermedades de cada gen. HPO la da como un término de frecuencia
# (Very frequent 80-99 %, Frequent 30-79 %, Occasional 5-29 %, Very rare 1-4 %) o como fracción de pacientes (n/m);
# los términos se convierten al centro de su intervalo y de cada gen se toma la frecuencia más alta
ev <- g2p[g2p$gene_symbol %in% genes & g2p$hpo_id %in% terminos, c("gene_symbol", "hpo_id", "frequency", "disease_id")]
centro <- c("HP:0040280" = 1, "HP:0040281" = 0.895, "HP:0040282" = 0.545, "HP:0040283" = 0.17,
            "HP:0040284" = 0.025, "HP:0040285" = 0)
a_numero <- function(f) {
  if (f %in% names(centro)) return(centro[[f]])
  if (grepl("^[0-9]+/[0-9]+$", f)) { x <- as.numeric(strsplit(f, "/")[[1]]); return(if (x[2] > 0) x[1] / x[2] else NA) }
  if (grepl("^[0-9.]+%$", f)) return(as.numeric(sub("%$", "", f)) / 100)
  NA
}
ev$frecuencia <- vapply(ev$frequency, a_numero, numeric(1))
maxima <- tapply(ev$frecuencia, ev$gene_symbol, function(x) if (all(is.na(x))) NA else max(x, na.rm = TRUE))

evidencia <- data.frame(gen = genes, frecuencia_maxima = round(as.vector(maxima[genes]), 3))
evidencia$grupo <- as.character(cut(evidencia$frecuencia_maxima, c(-Inf, 0.05, opt$fuerte, Inf), right = FALSE,
                                    labels = c("muy_rara", "ocasional", "frecuente")))
evidencia$enfermedades <- vapply(genes, function(x) paste(sort(unique(ev$disease_id[ev$gene_symbol == x])), collapse = ","), "")
red_completa <- simplify(graph_from_data_frame(read.delim(opt$red), directed = FALSE))
clusters <- read.delim(opt$clusters)
evidencia$modulo_string <- clusters$cluster[match(genes, clusters$gen)]
evidencia$grado_string <- ifelse(genes %in% V(red_completa)$name, degree(red_completa)[match(genes, V(red_completa)$name)], 0)
if (file.exists(opt$auc)) {
  auc <- read.delim(opt$auc)
  auc <- auc[auc$umbral == 700, ]
  evidencia$auc_loo_A6_700 <- auc$auc_loo[match(genes, auc$gen)]
}
evidencia <- evidencia[order(-evidencia$frecuencia_maxima, evidencia$gen), ]

# 3. Red previa solo con los genes de evidencia fuerte: las interacciones de STRING entre ellos, sin aislados,
# Louvain otra vez y comparación con los módulos de la red completa (NMI en los genes comunes)
fuertes <- evidencia$gen[evidencia$grupo %in% "frecuente"]
g_fuerte <- induced_subgraph(red_completa, intersect(V(red_completa)$name, fuertes))
g_fuerte <- delete_vertices(g_fuerte, degree(g_fuerte) == 0)
set.seed(123)
louvain_fuerte <- cluster_louvain(g_fuerte)
mod_fuerte <- setNames(as.vector(membership(louvain_fuerte)), V(g_fuerte)$name)
memb_completa <- clusters$cluster[match(V(red_completa)$name, clusters$gen)]
nmi_fuerte <- compare(mod_fuerte, clusters$cluster[match(names(mod_fuerte), clusters$gen)], method = "nmi")

resumen <- data.frame(
  medida = c("version_hpo", "genes", "genes_frecuente", "genes_ocasional", "genes_muy_rara", "genes_sin_frecuencia",
             "red_completa_genes", "red_completa_interacciones", "red_completa_modulos", "red_completa_modularidad",
             "red_fuerte_genes", "red_fuerte_interacciones", "red_fuerte_modulos", "red_fuerte_modularidad",
             "nmi_modulos_red_fuerte_vs_completa"),
  valor = c(version_hpo, length(genes), sum(evidencia$grupo %in% "frecuente"), sum(evidencia$grupo %in% "ocasional"),
            sum(evidencia$grupo %in% "muy_rara"), sum(is.na(evidencia$grupo)),
            vcount(red_completa), ecount(red_completa), length(unique(memb_completa)), round(modularity(red_completa, memb_completa), 3),
            vcount(g_fuerte), ecount(g_fuerte), length(louvain_fuerte), round(modularity(louvain_fuerte), 3), round(nmi_fuerte, 3))
)
modulos <- data.frame(gen = names(mod_fuerte), modulo_red_fuerte = mod_fuerte,
                      modulo_red_completa = clusters$cluster[match(names(mod_fuerte), clusters$gen)])
modulos <- modulos[order(modulos$modulo_red_fuerte, modulos$gen), ]
write.table(evidencia, file.path(opt$salida, "A10_evidencia_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(resumen, file.path(opt$salida, "A10_resumen.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(modulos, file.path(opt$salida, "A10_modulos_red_fuerte.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# 4. Figura: la red completa con cada gen coloreado por su evidencia (misma semilla que la figura de 04_clustering.R)
# y, si está la A6, el AUC de cada semilla según su evidencia
colores <- c(frecuente = "#7fbf7b", ocasional = "#fdb863", muy_rara = "#f4a582")
pdf(file.path(opt$salida, "A10_evidencia.pdf"), width = 9, height = 9)
grupo_nodo <- evidencia$grupo[match(V(red_completa)$name, evidencia$gen)]
set.seed(123)
plot(red_completa, vertex.color = colores[grupo_nodo], vertex.size = 9, vertex.label.cex = 0.9,
     vertex.label.color = "black", vertex.label.family = "sans", edge.color = "grey60",
     main = "Evidencia de cada gen: frecuencia de la disección en sus enfermedades")
legend("topright", bty = "n", fill = colores,
       legend = c("Frecuente (30 % o más)", "Ocasional (5-29 %)", "Muy rara (menos del 5 %)"))
if (!is.null(evidencia$auc_loo_A6_700)) {
  con_auc <- evidencia[!is.na(evidencia$auc_loo_A6_700) & !is.na(evidencia$grupo), ]
  con_auc$grupo <- factor(con_auc$grupo, levels = c("frecuente", "ocasional", "muy_rara"))
  # cada semilla es un punto, desplazado un poco al azar en horizontal para que no se tapen (solo afecta al dibujo)
  set.seed(123)
  x <- as.numeric(con_auc$grupo) + runif(nrow(con_auc), -0.15, 0.15)
  plot(x, con_auc$auc_loo_A6_700, xlim = c(0.5, 3.5), ylim = c(0, 1), xaxt = "n", pch = 19,
       col = ifelse(con_auc$grado_string > 0, "black", "red"),
       xlab = "Evidencia del gen", ylab = "AUC leave-one-out de la semilla (A6, umbral 700)",
       main = "Recuperación de cada semilla por el RWR según su evidencia")
  axis(1, at = 1:3, labels = c("Frecuente", "Ocasional", "Muy rara"))
  abline(h = 0.5, lty = 2, col = "grey50")
  legend("right", bty = "n", pch = 19, col = c("black", "red"),
         legend = c("Con interacciones en la red del mínimo", "Sin interacciones en la red del mínimo"))
}
invisible(dev.off())

# Resumen para la salida de launch.sh
cat("HPO", version_hpo, ":", length(genes), "genes; evidencia (frecuencia máxima de la disección en sus enfermedades):",
    paste(names(table(evidencia$grupo)), table(evidencia$grupo), collapse = ", "), "\n")
cat("Genes de cada grupo por módulo de la red completa:\n")
print(table(modulo = ifelse(is.na(evidencia$modulo_string), "sin_interacciones", evidencia$modulo_string), grupo = evidencia$grupo))
cat("Red completa:", vcount(red_completa), "genes,", ecount(red_completa), "interacciones,", length(unique(memb_completa)),
    "módulos, Q", round(modularity(red_completa, memb_completa), 3), "\n")
cat("Red previa (solo evidencia fuerte):", vcount(g_fuerte), "genes,", ecount(g_fuerte), "interacciones,",
    length(louvain_fuerte), "módulos, Q", round(modularity(louvain_fuerte), 3), "; NMI frente a la red completa", round(nmi_fuerte, 3), "\n")
if (!is.null(evidencia$auc_loo_A6_700)) {
  medias <- tapply(evidencia$auc_loo_A6_700, evidencia$grupo, mean, na.rm = TRUE)
  cat("AUC leave-one-out medio de la A6 (700) por grupo:", paste(names(medias), round(medias, 3), collapse = ", "),
      "; con interacciones en la red del mínimo:", round(mean(evidencia$auc_loo_A6_700[evidencia$grado_string > 0], na.rm = TRUE), 3),
      "; sin interacciones:", round(mean(evidencia$auc_loo_A6_700[evidencia$grado_string == 0], na.rm = TRUE), 3), "\n")
}
