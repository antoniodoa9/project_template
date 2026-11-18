#!/usr/bin/env Rscript
# A6_rwr_roc.R - Ampliación A6
# Prioriza genes candidatos con un random walk con reinicio (RWR) desde los genes de HPO (semillas) sobre la red humana completa de STRING, con varios umbrales, basada en Köhler.
# Valida el RWR escondiendo semillas, con leave-one-out y k-fold, y calcula a mano la curva ROC y el AUC.

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
  
  make_option(c("-k", "--bloques"), 
              type = "integer", 
              default = 5,
              help = "Bloques de la validación k-fold [por defecto: %default]"),
  
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

# Validación: se esconde cada bloque de semillas, se ejecuta el RWR con las demás y se mira dónde quedan las escondidas frente a los genes que no son de HPO
# El AUC de cada semilla escondida es la fracción de esos genes que quedan por debajo de ella: 1 si sale la primera, 0.5 como el azar (Tema 8 de Ranea)
validar <- function(W, es_semilla, bloque, r) {
  ids <- which(es_semilla)
  negativos <- which(!es_semilla)
  bloques <- sort(unique(bloque))

  # Un p0 por bloque, sin las semillas de ese bloque. Cada p0 es una columna de P0, así el RWR hace todos los bloques a la vez
  P0 <- sapply(bloques, function(b) {
    p0 <- as.numeric(es_semilla)
    p0[ids[bloque == b]] <- 0
    p0 / sum(p0)
  })
  P <- matrix(rwr(W, P0, r), ncol = length(bloques))

  sapply(seq_along(ids), function(i) {
    p <- P[, match(bloque[i], bloques)]
    mean(p[negativos] < p[ids[i]]) + 0.5 * mean(p[negativos] == p[ids[i]])
  })
}

tabla <- data.frame()
resumen <- data.frame()
validacion <- data.frame()

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
  # Ojo con los candidatos que tocan una sola semilla (SLC2A11-2, grado 5, y RIPK1, grado 102, con 700; GGCX con 900): salen por delante de TGFB1, que toca muchas más. No depende del grado del candidato sino del de la semilla que toca: una semilla con pocos vecinos reparte su probabilidad entre pocos genes y cada uno recibe mucho
  orden <- order(puntuacion, decreasing = TRUE)
  orden <- orden[!es_semilla[orden]][1:opt$candidatos]

  # Semilla vecina que más probabilidad pasa a cada candidato en un paso (columna de W por la puntuación de la semilla) y cuántos vecinos tiene esa semilla
  aporte <- as.matrix(W[orden, es_semilla, drop = FALSE] %*% Diagonal(x = puntuacion[es_semilla]))
  principal <- which(es_semilla)[max.col(aporte, ties.method = "first")]
  principal[rowSums(aporte) == 0] <- NA
  
  tabla <- rbind(tabla, data.frame(
    umbral = u,
    puesto = seq_along(orden),
    gen = V(g)$gen[orden],
    puntuacion = signif(puntuacion[orden], 3),
    grado = degree(g)[orden],
    semillas_vecinas = rowSums(A[orden, es_semilla, drop = FALSE] > 0),
    semilla_principal = V(g)$gen[principal],
    grado_semilla = degree(g)[principal]
  ))
  resumen <- rbind(resumen, data.frame(
    umbral = u, genes = gorder(g), interacciones = gsize(g),
    semillas = nrow(semillas), semillas_en_red = sum(en_red),
    semillas_fuera = if (all(en_red)) "(ninguna)" else paste(semillas$name[!en_red], collapse = ",")
  ))

  # Leave-one-out: cada semilla es su propio bloque. k-fold: las semillas repartidas al azar en k bloques
  n_semillas <- sum(es_semilla)
  auc_loo <- validar(W, es_semilla, seq_len(n_semillas), opt$reinicio)
  set.seed(123)
  auc_kfold <- validar(W, es_semilla, sample(rep(seq_len(opt$bloques), length.out = n_semillas)), opt$reinicio)
  validacion <- rbind(validacion, data.frame(
    umbral = u,
    gen = V(g)$gen[es_semilla],
    auc_loo = round(auc_loo, 4),
    auc_kfold = round(auc_kfold, 4)
  ))
}
write.table(tabla, 
            file.path(opt$salida, "A6_candidatos.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(resumen, 
            file.path(opt$salida, "A6_resumen_umbral.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(validacion, 
            file.path(opt$salida, "A6_validacion.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# AUC medio de cada umbral con los dos esquemas
auc <- aggregate(cbind(auc_loo, auc_kfold) ~ umbral, data = validacion, FUN = mean)
auc[, -1] <- round(auc[, -1], 3)
write.table(auc, 
            file.path(opt$salida, "A6_auc.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Curva ROC de cada umbral: para cada fracción de genes no HPO que se acepta (1 - especificidad), la fracción de semillas escondidas que se recupera (sensibilidad)
# Es la distribución acumulada de 1 - AUC de cada semilla escondida, y el área bajo ella es el AUC medio
pdf(file.path(opt$salida, "A6_roc.pdf"))
for (u in umbrales) {
  v <- validacion[validacion$umbral == u, ]
  plot(ecdf(1 - v$auc_loo), do.points = FALSE, verticals = TRUE, xlim = c(0, 1), ylim = c(0, 1),
    main = paste("Curva ROC del RWR, umbral", u), xlab = "1 - especificidad", ylab = "Sensibilidad")
  lines(ecdf(1 - v$auc_kfold), do.points = FALSE, verticals = TRUE, col = "red")
  abline(0, 1, lty = 2, col = "grey")
  legend("bottomright", lty = c(1, 1, 2), col = c("black", "red", "grey"),
    legend = c(paste("Leave-one-out, AUC", round(mean(v$auc_loo), 3)),
      paste(opt$bloques, "bloques, AUC", round(mean(v$auc_kfold), 3)),
      "Azar, AUC 0.5"))
}
invisible(dev.off())

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

cat("\nAUC medio escondiendo semillas (leave-one-out y", opt$bloques, "bloques):\n")
print(auc, row.names = FALSE)
