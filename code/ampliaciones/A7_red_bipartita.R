#!/usr/bin/env Rscript
# A7_red_bipartita.R - Ampliación 7
# Red bipartita gen-fenotipo (Perkins, Tema 11, redes multipartitas): une cada gen del término con los fenotipos
# de HPO que tiene anotados, la proyecta sobre los genes con el índice de Jaccard (cuántos fenotipos comparten dos genes)
# y comprueba, frente al azar, si los genes que interaccionan en STRING o que están en el mismo módulo de Louvain
# se parecen más en sus fenotipos de lo esperado (cierre genotipo-fenotipo, Tema 12).

suppressPackageStartupMessages({
  library(optparse)
  library(igraph)
})

# genes_to_phenotype.txt y hp.obo son grandes: hasta 10 minutos por descarga si el script se ejecuta suelto
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
              help = "Término HPO de partida, que se quita de los fenotipos porque lo tienen todos los genes [por defecto: %default]"),
  make_option(c("-n", "--repeticiones"), type = "integer", default = 1000,
              help = "Permutaciones al azar [por defecto: %default]"),
  make_option(c("-o", "--salida"), type = "character",
              default = "../results/ampliaciones/A7_red_bipartita",
              help = "Carpeta de salida [por defecto: %default]")
)
opt <- parse_args(OptionParser(option_list = opciones))
dir.create(opt$salida, recursive = TRUE, showWarnings = FALSE)

# 1. Fenotipos de cada gen: las anotaciones gen-fenotipo de HPO (el mismo archivo que usa 05_enriquecimiento.R)
url_hpo <- "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download"
g2p <- read.delim(file.path(url_hpo, "genes_to_phenotype.txt"))
genes <- sort(unique(read.delim(opt$genes)$name))

# Solo cuentan las anomalías fenotípicas: se quitan los modos de herencia (la rama HP:0000005 de la ontología,
# por ejemplo "Autosomal dominant inheritance"), que no son fenotipos, y el propio término de partida.
# La rama se saca de hp.obo: cada término dice de quién es hijo (is_a), y se baja desde HP:0000005 por los hijos
obo <- readLines(file.path(url_hpo, "hp.obo"))
id_actual <- cumsum(grepl("^id: HP:", obo))
ids <- sub("^id: ", "", obo[grepl("^id: HP:", obo)])
lineas_is_a <- grepl("^is_a: HP:", obo) & id_actual > 0
is_a <- data.frame(hijo = ids[id_actual[lineas_is_a]],
                   padre = sub("^is_a: (HP:[0-9]+).*", "\\1", obo[lineas_is_a]))
herencia <- "HP:0000005"
repeat {
  nuevos <- setdiff(is_a$hijo[is_a$padre %in% herencia], herencia)
  if (length(nuevos) == 0) break
  herencia <- c(herencia, nuevos)
}
version_hpo <- sub("^data-version: ", "", obo[grepl("^data-version:", obo)][1])

anot <- unique(g2p[g2p$gene_symbol %in% genes, c("gene_symbol", "hpo_id", "hpo_name")])
quitadas_herencia <- sum(anot$hpo_id %in% herencia)
anot <- anot[!anot$hpo_id %in% c(herencia, opt$hpo), ]

# 2. Red bipartita: dos tipos de nodos (genes y fenotipos) y una arista entre un gen y cada fenotipo que tiene
b <- graph_from_data_frame(unique(anot[, c("gene_symbol", "hpo_id")]), directed = FALSE)
V(b)$type <- V(b)$name %in% anot$hpo_id   # FALSE = gen, TRUE = fenotipo
fenotipos_gen <- degree(b)[!V(b)$type]

# Proyección sobre los genes: dos genes se unen si comparten algún fenotipo, con el número de fenotipos compartidos como peso
proy <- bipartite_projection(b, which = "false", multiplicity = TRUE)
ar <- as_data_frame(proy, what = "edges")
ar$compartidos <- ar$weight
# Índice de Jaccard: fenotipos compartidos / fenotipos de los dos genes juntos (Tema 11, índices de asociación)
ar$jaccard <- ar$compartidos / (fenotipos_gen[ar$from] + fenotipos_gen[ar$to] - ar$compartidos)

# Matriz de Jaccard entre todos los genes con fenotipos (0 si no comparten ninguno)
gf <- sort(names(fenotipos_gen))
J <- matrix(0, length(gf), length(gf), dimnames = list(gf, gf))
J[cbind(ar$from, ar$to)] <- ar$jaccard
J[cbind(ar$to, ar$from)] <- ar$jaccard

# 3. ¿Se parecen más en sus fenotipos los genes que interaccionan en STRING?
# Red del mínimo (no dirigida, sin aislados) y módulos de Louvain; solo los genes que tienen fenotipos
red <- read.delim(opt$red)
g <- simplify(graph_from_data_frame(red, directed = FALSE))
g <- induced_subgraph(g, intersect(V(g)$name, gf))
clusters <- read.delim(opt$clusters)
modulo <- setNames(clusters$cluster, clusters$gen)[V(g)$name]
pares <- as_edgelist(g)

jaccard_aristas <- function(nombres) {
  # nombres: qué gen ocupa cada nodo de la red; con los reales, el Jaccard medio de las interacciones de STRING
  mean(J[cbind(nombres[match(pares[, 1], V(g)$name)], nombres[match(pares[, 2], V(g)$name)])])
}
jaccard_modulo <- function(mod) {
  # Jaccard medio entre los pares de genes del mismo módulo
  dentro <- outer(mod, mod, "==") & upper.tri(diag(length(mod)))
  mean(J[V(g)$name, V(g)$name][dentro])
}

# Modelo al azar basado en los nodos (Tema 11): se barajan los nombres de los genes sobre la misma red y los mismos módulos,
# así se conserva la forma de la red y el tamaño de los módulos y solo cambia qué gen está en cada sitio
set.seed(123)
real_aristas <- jaccard_aristas(V(g)$name)
real_modulo <- jaccard_modulo(modulo)
todos_pares <- J[V(g)$name, V(g)$name][upper.tri(diag(vcount(g)))]
azar_aristas <- replicate(opt$repeticiones, jaccard_aristas(sample(V(g)$name)))
azar_modulo <- replicate(opt$repeticiones, jaccard_modulo(sample(modulo)))

# 4. ¿Salen los mismos módulos con los fenotipos que con STRING?
# Louvain sobre la proyección, con el Jaccard como peso, y NMI frente a los módulos de STRING en los genes comunes
fen <- graph_from_data_frame(ar[, c("from", "to")], directed = FALSE, vertices = data.frame(name = gf))
E(fen)$weight <- ar$jaccard
louvain_fen <- cluster_louvain(fen, weights = E(fen)$weight)
mod_fen <- setNames(as.vector(membership(louvain_fen)), V(fen)$name)
comunes <- intersect(names(mod_fen), V(g)$name)
real_nmi <- compare(mod_fen[comunes], modulo[comunes], method = "nmi")
azar_nmi <- replicate(opt$repeticiones, compare(sample(mod_fen[comunes]), modulo[comunes], method = "nmi"))

# p empírico como en la A2: (1 + nº de repeticiones al azar con un valor igual o mayor que el real) / (1 + repeticiones)
p_emp <- function(real, azar) (1 + sum(azar >= real)) / (1 + length(azar))

resumen <- data.frame(
  medida = c("jaccard_interacciones_string", "jaccard_mismo_modulo", "nmi_modulos_fenotipo_vs_string"),
  real = round(c(real_aristas, real_modulo, real_nmi), 3),
  media_azar = round(c(mean(azar_aristas), mean(azar_modulo), mean(azar_nmi)), 3),
  maximo_azar = round(c(max(azar_aristas), max(azar_modulo), max(azar_nmi)), 3),
  p_empirico = signif(c(p_emp(real_aristas, azar_aristas), p_emp(real_modulo, azar_modulo), p_emp(real_nmi, azar_nmi)), 3),
  genes = c(vcount(g), vcount(g), length(comunes)),
  pares = c(nrow(pares), sum(outer(modulo, modulo, "==")[upper.tri(diag(vcount(g)))]), NA)
)

# Tablas de salida
bipartita <- data.frame(
  medida = c("genes_con_fenotipos", "fenotipos", "aristas_gen_fenotipo", "anotaciones_de_herencia_quitadas",
             "pares_de_genes_que_comparten_algun_fenotipo", "pares_de_genes_posibles", "version_hpo"),
  valor = c(length(gf), sum(V(b)$type), ecount(b), quitadas_herencia, nrow(ar), choose(length(gf), 2), version_hpo)
)
genes_fen <- data.frame(gen = gf, fenotipos = as.vector(fenotipos_gen[gf]),
                        modulo_string = as.vector(modulo[gf]), modulo_fenotipos = as.vector(mod_fen[gf]))
genes_fen <- genes_fen[order(-genes_fen$fenotipos, genes_fen$gen), ]
frec <- table(anot$hpo_id)
nombre_fen <- setNames(anot$hpo_name, anot$hpo_id)
fen_compartidos <- data.frame(hpo_id = names(frec), fenotipo = nombre_fen[names(frec)], genes = as.vector(frec))
fen_compartidos <- fen_compartidos[order(-fen_compartidos$genes, fen_compartidos$hpo_id), ]
pares_tabla <- ar[, c("from", "to", "compartidos", "jaccard")]
names(pares_tabla)[1:2] <- c("gen1", "gen2")
pares_tabla$interaccion_string <- paste(pares_tabla$gen1, pares_tabla$gen2) %in% c(paste(pares[, 1], pares[, 2]), paste(pares[, 2], pares[, 1]))
pares_tabla$jaccard <- round(pares_tabla$jaccard, 3)
pares_tabla <- pares_tabla[order(-pares_tabla$jaccard), ]

write.table(resumen, file.path(opt$salida, "A7_resumen.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(bipartita, file.path(opt$salida, "A7_bipartita.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(genes_fen, file.path(opt$salida, "A7_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(head(fen_compartidos, 30), file.path(opt$salida, "A7_fenotipos_compartidos.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pares_tabla, file.path(opt$salida, "A7_pares_jaccard.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(repeticion = seq_len(opt$repeticiones), jaccard_interacciones = azar_aristas,
                       jaccard_mismo_modulo = azar_modulo, nmi = azar_nmi),
            file.path(opt$salida, "A7_azar.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Figura: una página con los tres histogramas al azar y el valor real en rojo, y otra con la matriz de Jaccard
pdf(file.path(opt$salida, "A7_red_bipartita.pdf"), width = 11, height = 5)
par(mfrow = c(1, 3))
dibujar <- function(azar, real, titulo, eje) {
  hist(azar, breaks = 30, col = "grey80", xlim = range(c(azar, real)), main = titulo, xlab = eje, ylab = "Repeticiones al azar")
  abline(v = real, col = "red", lwd = 2)
}
dibujar(azar_aristas, real_aristas, paste0("Interacciones de STRING\np empírico ", signif(p_emp(real_aristas, azar_aristas), 2)), "Jaccard medio")
dibujar(azar_modulo, real_modulo, paste0("Mismo módulo de Louvain\np empírico ", signif(p_emp(real_modulo, azar_modulo), 2)), "Jaccard medio")
dibujar(azar_nmi, real_nmi, paste0("Módulos de fenotipos frente a STRING\np empírico ", signif(p_emp(real_nmi, azar_nmi), 2)), "NMI")
par(mfrow = c(1, 1))
orden <- V(g)$name[order(modulo, V(g)$name)]
image(J[orden, orden][, rev(seq_along(orden))], axes = FALSE, col = hcl.colors(30, "YlOrRd", rev = TRUE),
      main = "Jaccard entre los genes de la red, ordenados por módulo de Louvain")
axis(1, at = seq(0, 1, length.out = length(orden)), labels = orden, las = 2, cex.axis = 0.5)
axis(2, at = seq(0, 1, length.out = length(orden)), labels = rev(orden), las = 2, cex.axis = 0.5)
invisible(dev.off())

cat("HPO", version_hpo, ":", length(gf), "genes,", sum(V(b)$type), "fenotipos,", ecount(b), "aristas gen-fenotipo (quitadas",
    quitadas_herencia, "anotaciones de herencia y el término", opt$hpo, ")\n")
cat("Pares de genes que comparten algún fenotipo:", nrow(ar), "de", choose(length(gf), 2), "\n")
cat("Jaccard medio de las", nrow(pares), "interacciones de STRING:", round(real_aristas, 3), "; de todos los pares de la red:",
    round(mean(todos_pares), 3), "; al azar", round(mean(azar_aristas), 3), "(máx", round(max(azar_aristas), 3), "), p",
    signif(p_emp(real_aristas, azar_aristas), 3), "\n")
cat("Jaccard medio dentro de los módulos:", round(real_modulo, 3), "; al azar", round(mean(azar_modulo), 3), "(máx",
    round(max(azar_modulo), 3), "), p", signif(p_emp(real_modulo, azar_modulo), 3), "\n")
cat("Módulos de la red de fenotipos:", max(mod_fen), "; NMI frente a STRING en", length(comunes), "genes:", round(real_nmi, 3),
    "; al azar", round(mean(azar_nmi), 3), "(máx", round(max(azar_nmi), 3), "), p", signif(p_emp(real_nmi, azar_nmi), 3), "\n")
cat("Fenotipos compartidos por más genes:", paste0(head(fen_compartidos$fenotipo, 5), " (", head(fen_compartidos$genes, 5), ")", collapse = ", "), "\n")
