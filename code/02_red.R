#!/usr/bin/env Rscript
# 02_red.R - Análisis mínimo, parte 2
# Busca los genes en STRING (versión 12.5), saca las interacciones entre ellos
# con el umbral indicado y cuenta los genes que se pierden.

suppressPackageStartupMessages({
  library(optparse)
  library(STRINGdb)
})

opciones <- list(
  make_option(c("-i", "--genes"), type = "character",
              default = "../results/minimo/genes_hpo.tsv",
              help = "TSV de genes de 01_genes.R [por defecto: %default]"),
  make_option(c("-s", "--umbral"), type = "integer", default = 700,
              help = "Puntuación mínima de STRING, de 0 a 1000 [por defecto: %default]"),
  make_option(c("-d", "--datos"), type = "character",
              default = "../results/string_data",
              help = "Carpeta donde STRINGdb guarda sus descargas [por defecto: %default]"),
  make_option(c("-o", "--salida"), type = "character",
              default = "../results/minimo",
              help = "Carpeta de salida [por defecto: %default]")
)
opt <- parse_args(OptionParser(option_list = opciones))

genes <- read.delim(opt$genes)
dir.create(opt$datos, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$salida, recursive = TRUE, showWarnings = FALSE)

# Versión fija de STRING, como en el ejemplo del profesor
string_db <- STRINGdb$new(version = "12.5", species = 9606,
                          score_threshold = opt$umbral, input_directory = opt$datos)

# Cada gen, a su identificador de STRING
string_ids <- string_db$map(genes, "name", removeUnmappedRows = TRUE)

# Interacciones entre nuestros genes, sin duplicados
red <- string_db$get_interactions(string_ids$STRING_id)
red <- unique(red[, 1:2])

# Identificadores de STRING -> nombres de gen
red$from <- string_ids$name[match(red$from, string_ids$STRING_id)]
red$to   <- string_ids$name[match(red$to, string_ids$STRING_id)]

# Genes perdidos: no están en STRING, o están pero no tienen ninguna interacción
no_mapeados <- setdiff(genes$name, string_ids$name)
sin_interaccion <- setdiff(string_ids$name, c(red$from, red$to))

cat("Umbral de STRING:", opt$umbral, "\n")
cat("Genes de HPO:", nrow(genes), "\n")
cat("Genes encontrados en STRING:", nrow(string_ids), "\n")
cat("Genes no encontrados en STRING:", length(no_mapeados), ":", no_mapeados, "\n")
cat("Genes sin ninguna interacción:", length(sin_interaccion), ":", sin_interaccion, "\n")
cat("Interacciones:", nrow(red), "\n")

# Se guardan la red y todos los genes encontrados en STRING
# (03_topologia.R los usa para incluir los genes sin interacción como nodos aislados)
write.table(red, file.path(opt$salida, "red.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(string_ids[, c("name", "STRING_id")], file.path(opt$salida, "genes_string.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
