#!/usr/bin/env Rscript
# A5_enriquecimiento.R - Ampliación A5
# Amplía el enriquecimiento de 05_enriquecimiento.R: GO de las tres ontologías (BP, MF y CC), quita los términos repetidos con simplify() y dibuja las redes de genes y términos (cnetplot y emapplot).

suppressPackageStartupMessages({
  library(optparse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
})

opciones <- list(
  make_option(c("-g", "--genes"), type = "character",
              default = "../results/minimo/tablas/genes_hpo.tsv",
              help = "TSV de genes de 01_genes.R [por defecto: %default]"),
  make_option(c("-o", "--salida"), type = "character",
              default = "../results/ampliaciones/A5_enriquecimiento",
              help = "Carpeta de salida de las tablas y la figura [por defecto: %default]")
)
opt <- parse_args(OptionParser(option_list = opciones))

dir.create(opt$salida, recursive = TRUE, showWarnings = FALSE)

genes <- read.delim(opt$genes)

# El mismo fondo que en 05_enriquecimiento.R: los genes con anotación en HPO
g2p <- read.delim("https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/genes_to_phenotype.txt")
fondo_hpo <- unique(g2p$gene_symbol)

# Las tres ontologías de GO: procesos biológicos (BP), funciones moleculares (MF)
# y componentes celulares (CC). El resto de argumentos, como en 05
ontologias <- c("BP", "MF", "CC")
enrich_go <- list()
enrich_simple <- list()
for (ont in ontologias) {
  enrich_go[[ont]] <- enrichGO(
    gene          = genes$name,
    universe      = fondo_hpo,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
  # simplify() (Enriquecimiento funcional de Perkins): por la jerarquía de GO
  # salen muchos términos casi iguales, y se queda con el más significativo.
  # Si una ontología no da ningún término, no hay nada que simplificar
  if (nrow(as.data.frame(enrich_go[[ont]])) > 0) {
    enrich_simple[[ont]] <- simplify(enrich_go[[ont]])
  } else {
    enrich_simple[[ont]] <- enrich_go[[ont]]
  }
  write.csv2(as.data.frame(enrich_go[[ont]]),
             file.path(opt$salida, paste0("GO_", ont, ".csv")), row.names = FALSE)
  write.csv2(as.data.frame(enrich_simple[[ont]]),
             file.path(opt$salida, paste0("GO_", ont, "_simplificado.csv")), row.names = FALSE)
}

# Cuántos términos quedan antes y después de simplify() en cada ontología
resumen <- data.frame(
  ontologia     = ontologias,
  terminos      = sapply(enrich_go, function(x) nrow(as.data.frame(x))),
  simplificados = sapply(enrich_simple, function(x) nrow(as.data.frame(x))),
  top_term      = sapply(enrich_simple, function(x) {
    res <- as.data.frame(x)
    if (nrow(res) == 0) "(ninguno)" else res$Description[1]
  })
)
write.table(resumen, file.path(opt$salida, "A5_resumen.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Figuras de enrichplot. Con enrichplot hace falta print()
pdf(file.path(opt$salida, "A5_enriquecimiento.pdf"), width = 10, height = 8)
# Un dotplot por ontología, con los términos ya simplificados
for (ont in ontologias) {
  if (nrow(as.data.frame(enrich_simple[[ont]])) > 0) {
    print(dotplot(enrich_simple[[ont]], showCategory = 10,
                  title = paste0("Enriquecimiento GO (", ont, ") simplificado, lista completa")))
  }
}
# cnetplot: qué genes tiene cada uno de los 5 primeros términos BP
# emapplot: red de términos, unidos si comparten genes (necesita pairwise_termsim)
bp <- enrich_simple[["BP"]]
if (nrow(as.data.frame(bp)) > 0) {
  print(cnetplot(bp, showCategory = 5))
}
if (nrow(as.data.frame(bp)) > 1) {
  print(emapplot(pairwise_termsim(bp), showCategory = 20))
}
invisible(dev.off())

# Resumen para la salida de launch.sh
cat("Genes de la lista completa:", nrow(genes), "; fondo de HPO:", length(fondo_hpo), "genes\n")
cat("Términos GO antes y después de simplify():\n")
print(resumen, row.names = FALSE)
