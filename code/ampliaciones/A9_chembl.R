#!/usr/bin/env Rscript
# A9_chembl.R - Ampliación A9
# Busca en ChEMBL qué genes de HPO (semillas) y qué candidatos de la A6 y de la A8 son diana de un fármaco aprobado o en ensayos clínicos.
# Usa los mecanismos de acción de ChEMBL, que están curados a mano para fármacos aprobados y en ensayos, y la diana humana de cada uno.
# Compara cuántos genes de cada conjunto tienen fármaco con conjuntos al azar de la subred de 700 con grado parecido (p empírico), como en la A6.

suppressPackageStartupMessages({
  library(optparse)
  library(httr)
  library(jsonlite)
})

opciones <- list(
  make_option(c("-i", "--subred"),
              type = "character",
              default = "../results/ampliaciones/A6_rwr_roc/A6_subred_700.tsv",
              help = "Subred de 700 que guarda la A6 [por defecto: %default]"),

  make_option(c("-g", "--genes"),
              type = "character",
              default = "../results/minimo/tablas/genes_hpo.tsv",
              help = "TSV de genes de 01_genes.R: las semillas [por defecto: %default]"),

  make_option(c("--a6"),
              type = "character",
              default = "../results/ampliaciones/A6_rwr_roc/A6_candidatos.tsv",
              help = "Candidatos de la A6; se usan los del umbral 700 [por defecto: %default]"),

  make_option(c("--a8"),
              type = "character",
              default = "../results/ampliaciones/A8_embedding/A8_candidatos.tsv",
              help = "Candidatos de la A8 [por defecto: %default]"),

  make_option(c("-f", "--fase"),
              type = "double",
              default = 1,
              help = "Fase clínica mínima del fármaco: 1, 2 y 3 son ensayos y 4 es aprobado [por defecto: %default]"),

  make_option(c("-a", "--azar"),
              type = "integer",
              default = 1000,
              help = "Conjuntos de genes al azar por cada conjunto real [por defecto: %default]"),

  make_option(c("--api"),
              type = "character",
              default = "https://www.ebi.ac.uk/chembl/api/data",
              help = "Dirección de la API de ChEMBL [por defecto: %default]"),

  make_option(c("-o", "--salida"),
              type = "character",
              default = "../results/ampliaciones/A9_chembl",
              help = "Carpeta de salida [por defecto: %default]")
)
opt <- parse_args(OptionParser(option_list = opciones))

dir.create(opt$salida, recursive = TRUE, showWarnings = FALSE)

# Pide a la API una consulta entera. ChEMBL devuelve como mucho 1000 filas por página y en page_meta$next dice dónde está la siguiente RETRY() repite la llamada si falla (hasta 5 veces), porque la API del EBI a veces tarda o corta
servidor <- sub("^(https?://[^/]+).*", "\\1", opt$api)
descargar <- function(recurso, consulta, clave) {
  filas <- list()
  siguiente <- paste0(opt$api, "/", recurso, ".json?limit=1000", consulta)

  while (!is.null(siguiente)) {
    respuesta <- RETRY("GET", siguiente, timeout(120), times = 5, quiet = TRUE)
    stop_for_status(respuesta)
    datos <- fromJSON(content(respuesta, as = "text", encoding = "UTF-8"), simplifyVector = FALSE)
    filas <- c(filas, datos[[clave]])
    siguiente <- datos$page_meta[["next"]]
    if (!is.null(siguiente)) siguiente <- paste0(servidor, siguiente)
  }
  filas
}

# Saca un campo de cada fila; si falta o viene vacío (null en el JSON) pone NA
campo <- function(filas, nombre) {
  sapply(filas, function(x) if (is.null(x[[nombre]])) NA else as.character(x[[nombre]]))
}

# Versión de ChEMBL, para el informe y para que se sepa de cuándo son los datos
# Si el servidor de estado de ChEMBL falla (a veces devuelve un error 500 con una página web en vez de JSON), el análisis sigue igual y la versión queda como no disponible
estado <- RETRY("GET", paste0(opt$api, "/status.json"), timeout(60), times = 3, quiet = TRUE)
version <- if (status_code(estado) == 200 && grepl("json", headers(estado)[["content-type"]])) {
  fromJSON(content(estado, as = "text", encoding = "UTF-8"))$chembl_db_version
} else {
  "no disponible"
}
# La fecha de la descarga sirve para saber de cuándo son los datos aunque no se sepa la versión
fecha <- format(Sys.Date())

# 1. Todos los mecanismos de acción. max_phase es la fase máxima del fármaco para cualquier enfermedad, no solo para la disección aórtica
# Se usa la molécula madre (parent): las sales de un mismo fármaco (por ejemplo, losartán y losartán potásico) cuentan como uno
mecanismos <- descargar("mechanism", "", "mechanisms")
mecanismos <- data.frame(
  molecula = campo(mecanismos, "parent_molecule_chembl_id"),
  diana = campo(mecanismos, "target_chembl_id"),
  accion = campo(mecanismos, "action_type"),
  mecanismo = campo(mecanismos, "mechanism_of_action"),
  fase = as.numeric(campo(mecanismos, "max_phase"))
)
mecanismos <- mecanismos[!is.na(mecanismos$fase) & mecanismos$fase >= opt$fase & !is.na(mecanismos$diana), ]

# 2. Dianas humanas que son una sola proteína, con su símbolo de gen. Los complejos y las familias se dejan fuera
# porque repartirían el fármaco entre todas sus subunidades aunque solo actúe sobre una
dianas <- descargar("target", "&organism=Homo%20sapiens&target_type=SINGLE%20PROTEIN", "targets")
dianas <- do.call(rbind, lapply(dianas, function(d) {
  sinonimos <- unlist(lapply(d$target_components, function(componente) {
    sapply(componente$target_component_synonyms, function(s) if (identical(s$syn_type, "GENE_SYMBOL")) s$component_synonym else NA)
  }))

  sinonimos <- unique(sinonimos[!is.na(sinonimos)])
  if (length(sinonimos) == 0) return(NULL)
  data.frame(diana = d$target_chembl_id, gen = sinonimos)
}))

# Un fármaco por gen y mecanismo: unique() quita las filas repetidas que salen cuando un fármaco tiene varios registros del mismo mecanismo
farmacos <- unique(merge(mecanismos, dianas, by = "diana"))

# 3. Subred de 700 y grado de cada gen, contando las interacciones de la tabla de la A6
subred <- read.delim(opt$subred)
grado <- table(c(subred$from, subred$to))
genes_red <- names(grado)
grado <- as.numeric(grado)

semillas <- intersect(read.delim(opt$genes)$name, genes_red)
a6 <- read.delim(opt$a6)
conjuntos <- list(
  semillas = semillas,
  A6 = a6$gen[a6$umbral == 700],
  A8 = read.delim(opt$a8)$gen
)

con_farmaco <- genes_red %in% farmacos$gen[farmacos$fase >= opt$fase]

# 4. Nombre de los fármacos que tocan algún gen de los conjuntos. Se piden de 50 en 50 para que la dirección no sea demasiado larga
farmacos <- farmacos[farmacos$gen %in% unlist(conjuntos), ]
ids <- unique(farmacos$molecula)

moleculas <- unlist(lapply(split(ids, ceiling(seq_along(ids) / 50)), function(trozo) {
  descargar("molecule", paste0("&molecule_chembl_id__in=", paste(trozo, collapse = ",")), "molecules")
}), recursive = FALSE)

nombres <- setNames(campo(moleculas, "pref_name"), campo(moleculas, "molecule_chembl_id"))
farmacos$farmaco <- nombres[farmacos$molecula]

# Si un fármaco no tiene nombre en ChEMBL se deja su identificador
farmacos$farmaco[is.na(farmacos$farmaco)] <- farmacos$molecula[is.na(farmacos$farmaco)]

# Tabla: una fila por conjunto, gen y fármaco. Un gen puede salir en dos conjuntos (por ejemplo, un candidato de la A6 y de la A8)
tabla <- do.call(rbind, lapply(names(conjuntos), function(nombre) {
  f <- farmacos[farmacos$gen %in% conjuntos[[nombre]], ]

  if (nrow(f) == 0) return(NULL)
  data.frame(conjunto = nombre, gen = f$gen, farmaco = f$farmaco, molecula = f$molecula,
    fase = f$fase, accion = f$accion, mecanismo = f$mecanismo)
}))
if (!is.null(tabla)) tabla <- tabla[order(tabla$conjunto, tabla$gen, -tabla$fase, tabla$farmaco), ]

# 5. Modelo al azar: ¿hay más genes con fármaco en cada conjunto que en genes de la subred de grado parecido?
# Los fármacos suelen ir a genes muy estudiados, que tienen muchas interacciones, así que comparar con genes al azar sin mirar el grado daría un p demasiado bueno
# Por cada gen del conjunto se coge un gen que no es de HPO del mismo decil de grado, sin repetir, como en el modelo al azar de la A6
decil <- cut(grado, unique(quantile(grado, 0:10 / 10)), include.lowest = TRUE, labels = FALSE)
es_semilla <- genes_red %in% semillas
set.seed(123)

azar <- sapply(conjuntos, function(conjunto) {
  en_conjunto <- genes_red %in% conjunto
  replicate(opt$azar, {
    falsos <- unlist(lapply(unique(decil[en_conjunto]), function(d) {
      posibles <- which(!es_semilla & decil == d)
      # sample.int y no sample(posibles, ...): si en un decil solo quedara un gen, sample(5, 1) cogería un número entre 1 y 5
      posibles[sample.int(length(posibles), sum(en_conjunto & decil == d))]
    }))
    sum(con_farmaco[falsos])
  })
})

real <- sapply(conjuntos, function(conjunto) sum(con_farmaco[genes_red %in% conjunto]))
# p empírico: fracción de conjuntos al azar con tantos genes con fármaco como el real o más, contando el real (por eso el 1 + ...)
p_empirico <- sapply(names(conjuntos), function(nombre) (1 + sum(azar[, nombre] >= real[nombre])) / (1 + opt$azar))

resumen <- data.frame(
  conjunto = names(conjuntos),
  genes = sapply(conjuntos, function(conjunto) sum(genes_red %in% conjunto)),
  con_farmaco = real,
  azar_medio = round(colMeans(azar), 2),
  azar_maximo = apply(azar, 2, max),
  p_empirico = signif(p_empirico, 3),
  chembl = version,
  fecha_descarga = fecha,
  row.names = names(conjuntos)
)

write.table(tabla,
            file.path(opt$salida, "A9_farmacos.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(resumen,
            file.path(opt$salida, "A9_resumen.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(data.frame(repeticion = seq_len(opt$azar), azar),
            file.path(opt$salida, "A9_azar.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Figura: para cada conjunto, el histograma de los conjuntos al azar y en rojo el real, como en la A6
pdf(file.path(opt$salida, "A9_azar.pdf"), width = 11, height = 4)
par(mfrow = c(1, 3))
for (nombre in names(conjuntos)) {
  hist(azar[, nombre], breaks = seq(-0.5, max(azar[, nombre], real[nombre]) + 0.5, 1), col = "grey80",
    xlim = c(0, max(azar[, nombre], real[nombre]) + 1),
    main = paste0(nombre, ": ", real[nombre], " de ", resumen[nombre, "genes"], " con fármaco\np empírico ", signif(p_empirico[nombre], 2)),
    xlab = "Genes con fármaco en el conjunto al azar", ylab = "Conjuntos al azar")
  abline(v = real[nombre], col = "red", lwd = 2)
}
invisible(dev.off())

cat("ChEMBL", version, "(descargado el", fecha, "):", nrow(mecanismos), "mecanismos con fase", opt$fase, "o más y", length(unique(dianas$diana)), "dianas humanas\n")
cat("Subred:", length(genes_red), "genes;", sum(con_farmaco), "con fármaco (", round(100 * mean(con_farmaco), 1), "% )\n")
for (nombre in names(conjuntos)) {
  cat(nombre, ":", real[nombre], "de", resumen[nombre, "genes"], "con fármaco; al azar", resumen[nombre, "azar_medio"],
    "de media y", resumen[nombre, "azar_maximo"], "como máximo; p empírico", signif(p_empirico[nombre], 3), "\n")
  cat("  ", paste(sort(unique(tabla$gen[tabla$conjunto == nombre])), collapse = ", "), "\n")
}
