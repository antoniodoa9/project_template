#!/usr/bin/env Rscript
# A6_rwr_roc.R - Ampliación A6
# Prioriza genes candidatos con un random walk con reinicio (RWR) desde los genes de HPO (semillas) sobre la red humana completa de STRING con el mismo umbral.

suppressPackageStartupMessages({
  library(optparse)
  library(igraph)
  library(Matrix)
  library(STRINGdb)
})

opciones <- list(
  make_option(c("-g", "--genes"),
              type = "character",
              default = "../results/minimo/tablas/genes_hpo.tsv",
              help = "TSV de genes de 01_genes.R: las semillas [por defecto: %default]"),
  
  make_option(c("-s", "--umbral"),
              type = "integer", 
              default = 700,
              help = "Puntuación mínima de STRING, la misma que en el mínimo [por defecto: %default]"),
  
  make_option(c("-d", "--datos"), 
              type = "character",
              default = "../results/string_data",
              help = "Carpeta de descargas de STRINGdb, la misma que usa 02_red.R [por defecto: %default]"),
  
  make_option(c("-r", "--reinicio"),
              type = "double", 
              default = 0.7,
              help = "Probabilidad de volver a las semillas en cada paso [por defecto: %default]"),
  
  make_option(c("-n", "--candidatos"), 
              type = "integer", 
              default = 20,
              help = "Genes candidatos que se guardan [por defecto: %default]"),
  
  make_option(c("-o", "--salida"), 
              type = "character",
              default = "../results/ampliaciones/A6_rwr_roc",
              help = "Carpeta de salida [por defecto: %default]")
)

opt <- parse_args(OptionParser(option_list = opciones))

dir.create(opt$salida, recursive = TRUE, showWarnings = FALSE)
genes <- read.delim(opt$genes)

# La misma versión de STRING y el mismo umbral que en 02_red.R, pero ahora con la red humana entera: así los candidatos pueden ser genes que no están en la lista
string_db <- STRINGdb$new(version = "12.5", 
                          species = 9606,
                          score_threshold = opt$umbral, input_directory = opt$datos)

semillas <- string_db$map(genes, "name", removeUnmappedRows = TRUE)

# Red humana completa con interacciones de puntuación igual o mayor que el umbral, sin interacciones repetidas y conservando la puntuación de STRING de cada una
g <- igraph::simplify(string_db$get_graph(), edge.attr.comb = list(combined_score = "max"))

# Nombre de gen de cada proteína de STRING, para la tabla de candidatos
proteinas <- string_db$get_proteins()
V(g)$gen <- proteinas$preferred_name[match(V(g)$name, proteinas$protein_external_id)]

# Semillas que están en la red (las que no tienen ninguna interacción con el umbral no están en el grafo)
en_red <- semillas$STRING_id %in% V(g)$name
sem_ids <- semillas$STRING_id[en_red]
es_semilla <- V(g)$name %in% sem_ids

# Matriz de transición W: la matriz de adyacencia con la puntuación de STRING de cada interacción como peso y cada columna dividida por su suma, así cada columna suma 1 (está en Temas 7 y 8 de Ranea)
A <- as_adjacency_matrix(g, attr = "combined_score", sparse = TRUE)
W <- A %*% Diagonal(x = 1 / colSums(A))   # %*% multiplicaión de matrices

# Random walk con reinicio, a mano: p <- (1 - r) W p + r p0, hasta que p no cambie. p0 Reparte la probabilidad por igual entre las semillas
rwr <- function(W, p0, r, tol = 1e-10) {
  p <- p0
  cambio <- Inf
  
 while (cambio >= tol){
   p_nuevo <- (1 - r) * (W %*% p) + r * p0
   cambio <- max(abs(p_nuevo - p))
   p <- p_nuevo
 }
  as.numeric(p)
}
p0 <- as.numeric(es_semilla) / sum(es_semilla)
puntuacion <- rwr(W, p0, opt$reinicio)

# Candidatos: los genes que no son semilla con más puntuación. Para cada uno, su
# grado y cuántas semillas toca directamente (Tema 12: contar vecinos no basta)
orden <- order(puntuacion, decreasing = TRUE)
orden <- orden[!es_semilla[orden]][1:opt$candidatos]

candidatos <- data.frame(
  puesto           = seq_along(orden),
  gen              = V(g)$gen[orden],
  puntuacion       = signif(puntuacion[orden], 3),
  grado            = degree(g)[orden],
  semillas_vecinas = rowSums(A[orden, es_semilla, drop = FALSE] > 0)
)
write.table(candidatos, 
            file.path(opt$salida, "A6_candidatos.tsv"),
            sep = "\t", quote = FALSE, 
            row.names = FALSE)

# Resumen para la salida de launch.sh
cat("Red de STRING completa (umbral", opt$umbral, "):", gorder(g), "genes y", gsize(g), "interacciones\n")
cat("Semillas:", nrow(semillas), "genes de HPO en STRING;", sum(en_red), "en la red\n")

if (any(!en_red)) cat("Semillas sin ninguna interacción:", semillas$name[!en_red], "\n")

cat("\nLos", opt$candidatos, "genes candidatos con más puntuación (reinicio", opt$reinicio, "):\n")
print(candidatos, row.names = FALSE)
