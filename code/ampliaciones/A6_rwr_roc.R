#!/usr/bin/env Rscript
# A6_rwr_roc.R - Ampliación A6
# Prioriza genes candidatos con un random walk con reinicio (RWR) desde los genes de HPO (semillas) sobre la subred de STRING formada por las semillas y sus vecinos directos, con varios umbrales, basada en Köhler.
# Valida el RWR escondiendo semillas, con leave-one-out y k-fold, y calcula a mano la curva ROC y el AUC.
# Con el umbral 700 compara el AUC con el de semillas al azar de grado parecido (p empírico) y guarda la subred para la A8.

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
  
  make_option(c("-a", "--azar"), 
              type = "integer", 
              default = 100,
              help = "Repeticiones del modelo con semillas al azar, con el umbral 700 [por defecto: %default]"),
  
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
auc_azar <- NULL

for (u in umbrales) {
  # La misma versión de STRING que en 02_red.R, pero se parte de la red humana entera para sacar de ella los vecinos de las semillas: así los candidatos pueden ser genes que no están en la lista
  string_db <- STRINGdb$new(version = "12.5", 
                            species = 9606,
                            score_threshold = u, 
                            input_directory = opt$datos)
  
  semillas <- string_db$map(genes, "name", removeUnmappedRows = TRUE)

  # Red humana completa con interacciones de puntuación igual o mayor que el umbral, sin interacciones repetidas y conservando la puntuación de STRING de cada una
  g <- igraph::simplify(string_db$get_graph(), edge.attr.comb = list(combined_score = "max"))

  # Nombre de gen de cada proteína de STRING, para la tabla de candidatos
  proteinas <- string_db$get_proteins()
  V(g)$gen <- proteinas$preferred_name[match(V(g)$name, proteinas$protein_external_id)]

  # Subred: las semillas y sus vecinos directos (1 paso), con todas las interacciones que hay entre ellos. La red humana entera tarda mucho y el RWR se dispersa por zonas poco fiables
  # ego() da los genes a 1 paso de cada semilla (la semilla incluida) e induced_subgraph() se queda con esos genes y las interacciones entre ellos, con su puntuación
  g <- induced_subgraph(g, unique(unlist(ego(g, order = 1, nodes = which(V(g)$name %in% semillas$STRING_id)))))

  # Semillas que están en la subred (las que no tienen ninguna interacción con el umbral no están en el grafo)
  en_red <- semillas$STRING_id %in% V(g)$name
  es_semilla <- V(g)$name %in% semillas$STRING_id[en_red]

  # Matriz de transición W: la matriz de adyacencia con la puntuación de STRING de cada interacción como peso y cada columna dividida por su suma, así cada columna suma 1 (está en Temas 7 y 8 de Ranea)
  A <- as_adjacency_matrix(g, attr = "combined_score", sparse = TRUE)
  W <- A %*% Diagonal(x = 1 / colSums(A))   # %*% multiplicaión de matrices
  p0 <- as.numeric(es_semilla) / sum(es_semilla)
  puntuacion <- rwr(W, p0, opt$reinicio)

  # Candidatos: los genes que no son semilla con más puntuación. Para cada uno, su grado y cuántas semillas toca directamente (Tema 12: contar vecinos no basta)
  # Ojo con los candidatos que tocan una sola semilla: pueden salir por delante de otros que tocan muchas más. No depende del grado del candidato sino del de la semilla que toca: una semilla con pocos vecinos reparte su probabilidad entre pocos genes y cada uno recibe mucho
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
    umbral = u, genes_subred = gorder(g), interacciones_subred = gsize(g),
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

  # Modelo al azar, solo con 700: ¿el AUC es alto solo porque las semillas tienen muchas interacciones? Se repite el leave-one-out con genes de la subred que no son de HPO y tienen grado parecido al de las semillas
  if (u == 700) {
    auc_real <- mean(auc_loo)
    grado <- degree(g)
    # Decil de grado de cada gen de la subred; unique() porque con muchos genes del mismo grado dos cortes pueden coincidir
    decil <- cut(grado, unique(quantile(grado, 0:10 / 10)), include.lowest = TRUE, labels = FALSE)
    set.seed(123)
    
    auc_azar <- sapply(seq_len(opt$azar), function(i) {
      # Por cada semilla, un gen que no es semilla del mismo decil de grado, sin repetir genes
      falsas <- unlist(lapply(unique(decil[es_semilla]), function(d) {
        sample(which(!es_semilla & decil == d), sum(es_semilla & decil == d))
      }))
      
      mean(validar(W, seq_len(gorder(g)) %in% falsas, seq_len(n_semillas), opt$reinicio))
    })
    # p empírico: fracción de repeticiones al azar con un AUC igual o mayor que el real, contando la real (por eso el 1 + ...)
    p_empirico <- (1 + sum(auc_azar >= auc_real)) / (1 + opt$azar)

    # La subred de 700 se guarda para la A8 (embedding): una fila por interacción, con el gen de cada extremo y la puntuación de STRING
    extremos <- ends(g, E(g), names = FALSE)
    write.table(data.frame(from = V(g)$gen[extremos[, 1]], to = V(g)$gen[extremos[, 2]], combined_score = E(g)$combined_score),
      file.path(opt$salida, "A6_subred_700.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}
write.table(tabla, 
            file.path(opt$salida, "A6_candidatos.tsv"),
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)

write.table(resumen, 
            file.path(opt$salida, "A6_resumen_umbral.tsv"),
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)

write.table(validacion, 
            file.path(opt$salida, "A6_validacion.tsv"),
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)

# AUC medio de cada umbral con los dos esquemas
auc <- aggregate(cbind(auc_loo, auc_kfold) ~ umbral, data = validacion, FUN = mean)
auc[, -1] <- round(auc[, -1], 3)
write.table(auc, 
            file.path(opt$salida, "A6_auc.tsv"),
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)

if (700 %in% umbrales) {
  write.table(data.frame(repeticion = seq_along(auc_azar), auc_medio = round(auc_azar, 4)),
              file.path(opt$salida, "A6_azar.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# Curva ROC de cada umbral: para cada fracción de genes no HPO que se acepta (1 - especificidad), la fracción de semillas escondidas que se recupera (sensibilidad)
# Es la distribución acumulada de 1 - AUC de cada semilla escondida, y el área bajo ella es el AUC medio
pdf(file.path(opt$salida, "A6_roc.pdf"))
for (u in umbrales) {
  v <- validacion[validacion$umbral == u, ]
  
  plot(ecdf(1 - v$auc_loo), 
       do.points = FALSE, 
       verticals = TRUE, 
       xlim = c(0, 1), 
       ylim = c(0, 1),
    main = paste("Curva ROC del RWR, umbral", u), xlab = "1 - especificidad", ylab = "Sensibilidad")
  
  lines(ecdf(1 - v$auc_kfold), do.points = FALSE, verticals = TRUE, col = "red")
  abline(0, 1, lty = 2, col = "grey")
  legend("bottomright", lty = c(1, 1, 2), col = c("black", "red", "grey"),
    legend = c(paste("Leave-one-out, AUC", round(mean(v$auc_loo), 3)),
      paste(opt$bloques, "bloques, AUC", round(mean(v$auc_kfold), 3)),
      "Azar, AUC 0.5"))
}

# Histograma de los AUC con semillas al azar y el AUC real en rojo, como la modularidad frente a redes aleatorias en Topología II (Perkins)
if (700 %in% umbrales) {
  hist(auc_azar, xlim = c(min(auc_azar, 0.5), 1), col = "grey80",
    main = paste0("AUC con semillas al azar de grado parecido, umbral 700\nEn rojo, los genes de HPO: AUC ", round(auc_real, 3), ", p empírico ", signif(p_empirico, 2)),
    xlab = "AUC medio (leave-one-out)", ylab = "Repeticiones")
  
  abline(v = auc_real, col = "red", lwd = 3)
}
invisible(dev.off())

# Candidatos que salen con todos los umbrales: los que no dependen del umbral
comunes <- names(which(table(tabla$gen) == length(umbrales)))

# Resumen para la salida de launch.sh
cat("Subredes de STRING por umbral (semillas y sus vecinos directos):\n")
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

if (700 %in% umbrales) {
  cat("\nModelo al azar con 700 (", opt$azar, "repeticiones): AUC real", round(auc_real, 3),
    "; AUC al azar medio", round(mean(auc_azar), 3), "y máximo", round(max(auc_azar), 3),
    "; p empírico", signif(p_empirico, 2), "\n")
}
