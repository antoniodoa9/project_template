#!/usr/bin/env Rscript
# 01_genes.R - Análisis mínimo, parte 1
# Descarga de la API de HPO los genes asociados a un término y los guarda en un TSV.

suppressPackageStartupMessages({
  library(optparse)
  library(httr)
  library(jsonlite)
})

opciones <- list(
  make_option(c("-p", "--hpo"), type = "character", default = "HP:0002647",
              help = "Término HPO [por defecto: %default]"),
  make_option(c("-o", "--salida"), type = "character",
              default = "../results/minimo/genes_hpo.tsv",
              help = "Archivo TSV de salida [por defecto: %default]")
)
opt <- parse_args(OptionParser(option_list = opciones))

# Petición a la API de HPO (igual que en el ejemplo mínimo del profesor)
url <- paste0("https://ontology.jax.org/api/network/annotation/", opt$hpo)
respuesta <- GET(url)
stop_for_status(respuesta)  # si la API falla, el script se para con error

# httr::content y no content: Biobase (de Bioconductor) tapa la función content
datos <- httr::content(respuesta, as = "text", encoding = "UTF-8")
genes <- fromJSON(datos)$genes  # columnas: id (NCBIGene:...) y name (símbolo del gen)

dir.create(dirname(opt$salida), recursive = TRUE, showWarnings = FALSE)
write.table(genes, opt$salida, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Término", opt$hpo, ":", nrow(genes), "genes guardados en", opt$salida, "\n")
