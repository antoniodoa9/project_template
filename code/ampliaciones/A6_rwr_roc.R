#!/usr/bin/env Rscript
# A6_rwr_roc.R - Ampliación A6
# Prioriza genes candidatos con un random walk con reinicio (RWR) desde los genes de HPO (semillas) sobre la red humana completa de STRING, con varios umbrales, basada en Köhler.

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
  
  make_option(c("-u", "--umbrales"), 
              type = "character", 
              default = "500,700,900",
              help = "Puntuaciones mínimas de STRING separadas por comas; 700 es la del mínimo [por defecto: %default]"),
  
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
              help = "Genes candidatos que se guardan por umbral [por defecto: %default]"),
  
  make_option(c("-o", "--salida"), 
              type = "character",
              default = "../results/ampliaciones/A6_rwr_roc",
              help = "Carpeta de salida [por defecto: %default]")
)
opt <- parse_args(OptionParser(option_list = opciones))

dir.create(opt$salida, recursive = TRUE, showWarnings = FALSE)
genes <- read.delim(opt$genes)
umbrales <- as.integer(strsplit(opt$umbrales, ",")[[1]])

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

tabla <- data.frame()
resumen <- data.frame()

for (u in umbrales) {
  # La misma versión de STRING que en 02_red.R, pero con la red humana entera: así los candidatos pueden ser genes que no están en la lista
  string_db <- STRINGdb$new(version = "12.5", 
                            species = 9606,
                            score_threshold = u, input_directory = opt$datos)
  
  semillas <- string_db$map(genes, "name", removeUnmappedRows = TRUE)

  # Red humana completa con interacciones de puntuación igual o mayor que el umbral, sin interacciones repetidas y conservando la puntuación de STRING de cada una
  g <- igraph::simplify(string_db$get_graph(), edge.attr.comb = list(combined_score = "max"))

  # Nombre de gen de cada proteína de STRING, para la tabla de candidatos
  proteinas <- string_db$get_proteins()
  V(g)$gen <- proteinas$preferred_name[match(V(g)$name, proteinas$protein_external_id)]

  # Semillas que están en la red (las que no tienen ninguna interacción con el umbral no están en el grafo)
  en_red <- semillas$STRING_id %in% V(g)$name
  es_semilla <- V(g)$name %in% semillas$STRING_id[en_red]

  # Matriz de transición W: la matriz de adyacencia con la puntuación de STRING de cada interacción como peso y cada columna dividida por su suma, así cada columna suma 1 (está en Temas 7 y 8 de Ranea)
  A <- as_adjacency_matrix(g, attr = "combined_score", sparse = TRUE)
  W <- A %*% Diagonal(x = 1 / colSums(A))   # %*% multiplicaión de matrices
  p0 <- as.numeric(es_semilla) / sum(es_semilla)
  puntuacion <- rwr(W, p0, opt$reinicio)

  # Candidatos: los genes que no son semilla con más puntuación. Para cada uno, su grado y cuántas semillas toca directamente (Tema 12: contar vecinos no basta)
  orden <- order(puntuacion, decreasing = TRUE)
  orden <- orden[!es_semilla[orden]][1:opt$candidatos]
  
  tabla <- rbind(tabla, data.frame(
    umbral           = u,
    puesto           = seq_along(orden),
    gen              = V(g)$gen[orden],
    puntuacion       = signif(puntuacion[orden], 3),
    grado            = degree(g)[orden],
    semillas_vecinas = rowSums(A[orden, es_semilla, drop = FALSE] > 0)
  ))
  resumen <- rbind(resumen, data.frame(
    umbral = u, genes = gorder(g), interacciones = gsize(g),
    semillas = nrow(semillas), semillas_en_red = sum(en_red),
    semillas_fuera = paste(semillas$name[!en_red], collapse = ",")
  ))
}
write.table(tabla, 
            file.path(opt$salida, "A6_candidatos.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(resumen, 
            file.path(opt$salida, "A6_resumen_umbral.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Candidatos que salen con todos los umbrales: los que no dependen del umbral
comunes <- names(which(table(tabla$gen) == length(umbrales)))

# Resumen para la salida de launch.sh
cat("Redes de STRING completas por umbral:\n")
print(resumen, row.names = FALSE)

# Puesto de cada candidato con cada umbral; NA es que con ese umbral no está entre los primeros
puestos <- reshape(tabla[, c("gen", "umbral", "puesto")], 
                   idvar = "gen",
                   timevar = "umbral", 
                   direction = "wide")

cat("\nPuesto de cada candidato con cada umbral (reinicio", opt$reinicio, "):\n")
print(puestos, row.names = FALSE)
cat("\nCandidatos con todos los umbrales (", length(comunes), "):", comunes, "\n")
