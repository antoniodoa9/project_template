#!/usr/bin/env Rscript
# 05_enriquecimiento.R - Análisis mínimo, parte 6
# Enriquecimiento GO (procesos biológicos) con clusterProfiler de la lista
# completa de genes de HPO y de cada cluster de 5 genes o más.

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
  make_option(c("-c", "--clusters"), type = "character",
              default = "../results/minimo/tablas/clusters.tsv",
              help = "TSV gen-cluster de 04_clustering.R [por defecto: %default]"),
  make_option(c("-m", "--minimo"), type = "integer", default = 5,
              help = "Genes mínimos de un cluster para enriquecerlo [por defecto: %default]"),
  make_option(c("-t", "--tablas"), type = "character",
              default = "../results/minimo/tablas",
              help = "Carpeta de salida de las tablas [por defecto: %default]"),
  make_option(c("-f", "--figuras"), type = "character",
              default = "../results/minimo/figuras",
              help = "Carpeta de salida de las figuras [por defecto: %default]")
)
opt <- parse_args(OptionParser(option_list = opciones))

dir.create(opt$tablas, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$figuras, recursive = TRUE, showWarnings = FALSE)

genes <- read.delim(opt$genes)
clusters <- read.delim(opt$clusters)

# ===== Lista completa =====
# Fondo: todos los genes con anotación en HPO, como propone el ejemplo. Nuestra
# lista sale de HPO, así que solo podría contener genes de HPO
g2p <- read.delim("https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/genes_to_phenotype.txt")
fondo_hpo <- unique(g2p$gene_symbol)

# ORA con términos GO de procesos biológicos, como en el ejemplo
enrich_results <- enrichGO(
  gene          = genes$name,   # los 41 genes de HPO
  universe      = fondo_hpo,    # fondo: genes con anotación en HPO
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2
)
tabla_lista <- as.data.frame(enrich_results)
write.csv2(tabla_lista, file.path(opt$tablas, "GO_lista_completa.csv"), row.names = FALSE)

# Figuras del ejemplo (dotplot y barplot). Con enrichplot hace falta print()
if (nrow(tabla_lista) > 0) {
  pdf(file.path(opt$figuras, "GO_lista_completa_dotplot.pdf"), width = 9)
  print(dotplot(enrich_results, showCategory = 10, title = "Enriquecimiento GO (BP), lista completa"))
  invisible(dev.off())
  pdf(file.path(opt$figuras, "GO_lista_completa_barplot.pdf"), width = 9)
  print(barplot(enrich_results, showCategory = 10, title = "Enriquecimiento GO (BP), lista completa"))
  invisible(dev.off())
}

# ===== Cada cluster =====
# Solo los clusters con el mínimo de genes: con 2 o 3 genes salen términos
# "significativos" que no dicen nada (lo explica el ejemplo)
genes_clusters <- split(clusters$gen, clusters$cluster)
genes_clusters <- genes_clusters[lengths(genes_clusters) >= opt$minimo]

# Fondo: los genes de la red, porque un cluster solo puede tener genes de la red
# (como en el ejemplo)
fondo_red <- clusters$gen

enrichment_results_clusters <- list()
for (i in names(genes_clusters)) {
  enrichment_results_clusters[[i]] <- enrichGO(
    gene          = genes_clusters[[i]],
    universe      = fondo_red,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
  tabla_cluster <- as.data.frame(enrichment_results_clusters[[i]])
  write.csv2(tabla_cluster, file.path(opt$tablas, paste0("GO_cluster_", i, ".csv")), row.names = FALSE)
  if (nrow(tabla_cluster) > 0) {
    pdf(file.path(opt$figuras, paste0("GO_cluster_", i, "_dotplot.pdf")), width = 9)
    print(dotplot(enrichment_results_clusters[[i]], showCategory = 10,
                  title = paste("Enriquecimiento GO (BP), cluster", i)))
    invisible(dev.off())
  }
}

# Término más significativo de cada cluster, como en el ejemplo
cluster_summary <- data.frame(
  cluster  = names(genes_clusters),
  genes    = lengths(genes_clusters),
  terminos = sapply(enrichment_results_clusters, function(x) nrow(as.data.frame(x))),
  top_term = sapply(enrichment_results_clusters, function(x) {
    res <- as.data.frame(x)
    if (nrow(res) == 0) "(ninguno)" else res$Description[1]
  })
)
write.table(cluster_summary, file.path(opt$tablas, "GO_clusters_resumen.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Resumen para la salida de launch.sh
cat("Genes de la lista completa:", nrow(genes), "; fondo de HPO:", length(fondo_hpo), "genes\n")
cat("Términos GO (BP) de la lista completa:", nrow(tabla_lista), "\n")
if (nrow(tabla_lista) > 0) cat("Primeros:", paste(head(tabla_lista$Description, 3), collapse = "; "), "\n")
cat("Clusters con", opt$minimo, "genes o más:", length(genes_clusters), "\n")
print(cluster_summary, row.names = FALSE)
