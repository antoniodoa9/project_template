#!/usr/bin/env Rscript
# 03_topologia.R - Análisis mínimo, partes 3 y 4
# Pasa la red a igraph con los genes sin interacción como nodos aislados,
# calcula sus propiedades y la dibuja.

suppressPackageStartupMessages({
  library(optparse)
  library(igraph)
})

opciones <- list(
  make_option(c("-i", "--red"), type = "character",
              default = "../results/minimo/tablas/red.tsv",
              help = "TSV de la red de 02_red.R [por defecto: %default]"),
  make_option(c("-g", "--genes"), type = "character",
              default = "../results/minimo/tablas/genes_string.tsv",
              help = "TSV de los genes encontrados en STRING, de 02_red.R [por defecto: %default]"),
  make_option(c("-o", "--salida"), type = "character",
              default = "../results/minimo/tablas",
              help = "Carpeta de salida de las tablas [por defecto: %default]"),
  make_option(c("-f", "--figuras"), type = "character",
              default = "../results/minimo/figuras",
              help = "Carpeta de salida de las figuras [por defecto: %default]")
)
opt <- parse_args(OptionParser(option_list = opciones))

red <- read.delim(opt$red)
genes <- read.delim(opt$genes)
dir.create(opt$salida, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$figuras, recursive = TRUE, showWarnings = FALSE)

# Red no dirigida, como en el ejemplo del profesor. Con vertices= entran
# también los genes sin ninguna interacción, como nodos aislados
# (el ejemplo lo propone), para que cuenten en los nodos y en la densidad
g <- graph_from_data_frame(red, directed = FALSE, vertices = genes)

# Propiedades de la red, como en el ejemplo
grado <- degree(g)
componentes <- components(g)

# La red no es conexa, y la cercanía falla en componentes pequeñas (el ejemplo
# lo avisa). Se dan las dos soluciones que propone: la cercanía dentro de la
# componente más grande y la centralidad armónica, que vale para toda la red
cercania <- closeness(largest_component(g))
armonica <- harmonic_centrality(g)

# Tabla por gen, de más a menos conectado
# (los genes fuera de la componente más grande quedan sin cercanía: NA)
metricas <- data.frame(gen = V(g)$name,
                       grado = grado,
                       cercania_componente_mayor = round(cercania[V(g)$name], 4),
                       centralidad_armonica = round(armonica, 4))
metricas <- metricas[order(-metricas$grado, metricas$gen), ]
write.table(metricas, file.path(opt$salida, "topologia_genes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Tabla resumen de la red
resumen <- data.frame(
  propiedad = c("nodos", "nodos_con_interaccion", "nodos_aislados", "aristas",
                "grado_medio", "densidad", "componentes", "nodos_componente_mayor"),
  valor = c(gorder(g), sum(grado > 0), sum(grado == 0), gsize(g),
            round(mean(grado), 3), round(edge_density(g), 4),
            componentes$no, max(componentes$csize)))
write.table(resumen, file.path(opt$salida, "topologia_resumen.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Figuras, como en el ejemplo: tamaño del nodo según el grado (+3 para que
# los aislados, con grado 0, también se vean). Con 41 genes las etiquetas
# se leen, así que se dejan
pdf(file.path(opt$figuras, "red_fr.pdf"))
set.seed(123)
plot(g, layout = layout_with_fr(g), vertex.color = "skyblue",
     vertex.size = grado + 3, vertex.label.cex = 0.6, vertex.label.color = "black",
     edge.width = 0.5, main = "Red (Fruchterman-Reingold)")
invisible(dev.off())

pdf(file.path(opt$figuras, "red_circular.pdf"))
plot(g, layout = layout_in_circle(g), vertex.color = "lightgreen",
     vertex.size = grado + 3, vertex.label.cex = 0.6, vertex.label.color = "black",
     edge.width = 0.5, main = "Red (disposición circular)")
invisible(dev.off())

# Resumen para la salida de launch.sh
cat("Nodos:", gorder(g), "(", sum(grado == 0), "aislados )\n")
cat("Aristas:", gsize(g), "\n")
cat("Grado medio:", round(mean(grado), 3), "\n")
cat("Densidad:", round(edge_density(g), 4), "\n")
cat("Componentes:", componentes$no, "; la mayor tiene", max(componentes$csize), "nodos\n")
cat("Genes con más interacciones:", head(metricas$gen, 5), "\n")
