#!/usr/bin/env python3
# A8_embedding.py - Ampliación A8
# Embedding de la subred de 700 de la A6 (semillas y sus vecinos directos) con paseos aleatorios, como DeepWalk y node2vec (Tema 13 de Perkins).
# Los vectores salen de la matriz de PMI de los genes que aparecen juntos en los paseos, reducida con SVD: es lo que calcula word2vec skip-gram (Levy y Goldberg 2014; Qiu et al. 2018).
# El embedding se usa para predecir: se esconden interacciones y semillas y se mira si los vectores las recuperan, comparado con el azar y con el random walk de la A6.
# Además agrupa los vectores con k-means y lo compara con Louvain (NMI), y dibuja los vectores con PCA.

import argparse
import csv
import os
import random

import igraph
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import normalized_mutual_info_score
from sklearn.utils.extmath import randomized_svd

parser = argparse.ArgumentParser(description="Ampliación A8: embedding de la subred de la A6")
parser.add_argument("-i", "--subred", default="../results/ampliaciones/A6_rwr_roc/A6_subred_700.tsv",
                    help="TSV de la subred de 700 de la A6 (from, to, combined_score) [por defecto: %(default)s]")

parser.add_argument("-g", "--genes", default="../results/minimo/tablas/genes_hpo.tsv",
                    help="TSV de genes de 01_genes.R: las semillas [por defecto: %(default)s]")

parser.add_argument("-a", "--a6", default="../results/ampliaciones/A6_rwr_roc",
                    help="Carpeta de resultados de la A6, para comparar [por defecto: %(default)s]")

parser.add_argument("-d", "--dimensiones", type=int, default=128,
                    help="Números de cada vector [por defecto: %(default)s]")

parser.add_argument("-n", "--candidatos", type=int, default=20,
                    help="Genes candidatos que se guardan [por defecto: %(default)s]")

parser.add_argument("-o", "--salida", default="../results/ampliaciones/A8_embedding",
                    help="Carpeta de salida [por defecto: %(default)s]")

opt = parser.parse_args()
os.makedirs(opt.salida, exist_ok=True)

# Todo lo que tiene azar sale de la semilla 123: random para igraph y numpy para elegir interacciones y para la SVD
random.seed(123)
rng = np.random.default_rng(123)

# Subred no dirigida, con la puntuación de STRING como peso de cada interacción
with open(opt.subred) as f:
    aristas = [(fila["from"], fila["to"], float(fila["combined_score"])) for fila in csv.DictReader(f, delimiter="\t")]
g = igraph.Graph.TupleList(aristas, directed=False, weights=True)
genes = g.vs["name"]
n = g.vcount()

with open(opt.genes) as f:
    nombres_hpo = {fila["name"] for fila in csv.DictReader(f, delimiter="\t")}
es_semilla = np.array([gen in nombres_hpo for gen in genes])


# Embedding como node2vec con p = q = 1 (DeepWalk): paseos aleatorios por la red, que hacen de "frases" con los genes como "palabras"
# 10 paseos de 80 pasos por gen, ventana de 10 y 128 dimensiones, los valores de node2vec y del ejemplo del Tema 13
def embedding(red):
    n_red = red.vcount()
    paseos = []

    for vuelta in range(10):
        orden = list(range(n_red))
        random.shuffle(orden)

        for v in orden:
            paseos.append(np.array(red.random_walk(v, 80, stuck="return", weights="weight")))

    # Cuántas veces aparecen juntos dos genes a 10 pasos o menos dentro de un mismo paseo (la "ventana" de word2vec)
    filas, columnas = [], []
    for camino in paseos:
        for d in range(1, 11):
            filas += [camino[:-d], camino[d:]]
            columnas += [camino[d:], camino[:-d]]

    filas, columnas = np.concatenate(filas), np.concatenate(columnas)
    C = sparse.coo_matrix((np.ones(len(filas)), (filas, columnas)), shape=(n_red, n_red)).tocsr()
    # PMI positiva: log(C_ij * total / (veces de i * veces de j)), con los valores negativos a 0. Mide si dos genes salen juntos más de lo esperado por azar
    total = C.sum()
    fila = np.asarray(C.sum(axis=1)).ravel()
    C = C.tocoo()
    pmi = np.log(C.data * total / (fila[C.row] * fila[C.col]))
    M = sparse.coo_matrix((np.maximum(pmi, 0), (C.row, C.col)), shape=(n_red, n_red)).tocsr()
    # SVD: se queda con las 128 direcciones que mejor resumen la matriz. Vector de cada gen = U * raíz(S), como Levy y Goldberg
    U, S, Vt = randomized_svd(M, opt.dimensiones, random_state=123)
    vectores = U * np.sqrt(S)
    # Cada vector dividido por su longitud: así el producto de dos vectores es su coseno (un gen sin interacciones se queda con el vector 0)
    largo = np.linalg.norm(vectores, axis=1, keepdims=True)

    return vectores / np.where(largo == 0, 1, largo)


# AUC a mano, como en la A6: fracción de negativos que quedan por debajo de cada positivo, con los empates a 0.5
def auc(positivos, negativos):
    positivos = np.asarray(positivos)[:, None]
    negativos = np.asarray(negativos)[None, :]
    return float(np.mean((negativos < positivos) + 0.5 * (negativos == positivos)))


X = embedding(g)
# =============================================================================
# ========================= 1. Esconder interacciones =========================
# =============================================================================
# Se esconde el 10 % de las interacciones, se rehace el embedding sin ellas y se mira si esas parejas quedan más cerca (coseno) que parejas al azar sin interacción
escondidas = rng.choice(g.ecount(), size=g.ecount() // 10, replace=False)
pares_pos = [g.es[e].tuple for e in escondidas]
g_resto = g.copy()
g_resto.delete_edges(escondidas)
X_resto = embedding(g_resto)

# Parejas sin interacción, tantas como escondidas
pares_neg = []
while len(pares_neg) < len(pares_pos):
    a, b = rng.choice(n, size=2, replace=False)
    if not g.are_adjacent(a, b):
        pares_neg.append((a, b))


def coseno(V, pares):
    return [float(V[a] @ V[b]) for a, b in pares]


auc_aristas = auc(coseno(X_resto, pares_pos), coseno(X_resto, pares_neg))
# Referencia sencilla: el índice de Jaccard (vecinos comunes / vecinos totales) en la red sin las escondidas
auc_jaccard = auc(g_resto.similarity_jaccard(pairs=pares_pos), g_resto.similarity_jaccard(pairs=pares_neg))

# ====================================================================================================
# ========================= 2. Esconder semillas (leave-one-out, como la A6) =========================
# ====================================================================================================
# Puntuación de cada gen: su coseno con la media de los vectores de las demás semillas. La semilla escondida se compara con los genes que no son de HPO
ids = np.where(es_semilla)[0]
negativos = np.where(~es_semilla)[0]
auc_semillas = []

for s in ids:
    centro = X[ids[ids != s]].mean(axis=0)
    puntuacion = X @ centro
    auc_semillas.append(auc([puntuacion[s]], puntuacion[negativos]))
auc_loo = float(np.mean(auc_semillas))
# Otra forma de puntuar: el coseno con la semilla más parecida en vez de con la media. Las semillas están en módulos distintos y la media de sus vectores los mezcla
cos_semillas = X @ X[ids].T
auc_semillas_cercana = []
for j, s in enumerate(ids):
    # Se quita la columna de la semilla escondida y cada gen se queda con su coseno más alto con las demás
    puntuacion = np.delete(cos_semillas, j, axis=1).max(axis=1)
    auc_semillas_cercana.append(auc([puntuacion[s]], puntuacion[negativos]))
auc_cercana = float(np.mean(auc_semillas_cercana))
# Referencia: ordenar los genes solo por su grado. En la subred de 1 paso las semillas son los centros, así que el grado ya las distingue del resto
grado = np.array(g.degree())
auc_grado = auc(grado[ids], grado[negativos])

# AUC leave-one-out del random walk de la A6 con 700, para compararlo
auc_rwr = None
ruta_a6 = os.path.join(opt.a6, "A6_auc.tsv")

if os.path.exists(ruta_a6):
    with open(ruta_a6) as f:
        for fila in csv.DictReader(f, delimiter="\t"):
            if fila["umbral"] == "700":
                auc_rwr = float(fila["auc_loo"])

# ===============================================================================
# ========================= 3. k-means frente a Louvain =========================
# ===============================================================================
# Louvain sobre la subred con los pesos (igraph, como en el mínimo) y k-means con tantos grupos como módulos de Louvain
louvain = g.community_multilevel(weights="weight")
k = len(louvain)
grupos = KMeans(n_clusters=k, n_init=10, random_state=123).fit_predict(X)
nmi = normalized_mutual_info_score(louvain.membership, grupos)
# Azar: la NMI con los módulos de Louvain barajados (100 veces)
nmi_azar = [normalized_mutual_info_score(rng.permutation(louvain.membership), grupos) for i in range(100)]

# =================================================================
# ========================= 4. Candidatos =========================
# =================================================================
# Genes que no son de HPO más cercanos (coseno) a la media de todas las semillas
puntuacion = X @ X[ids].mean(axis=0)
orden = [i for i in np.argsort(-puntuacion) if not es_semilla[i]][:opt.candidatos]
candidatos_a6 = set()
ruta_cand = os.path.join(opt.a6, "A6_candidatos.tsv")

if os.path.exists(ruta_cand):
    with open(ruta_cand) as f:
        candidatos_a6 = {fila["gen"] for fila in csv.DictReader(f, delimiter="\t") if fila["umbral"] == "700"}

# =============================================================
# ========================= 5. Tablas =========================
# =============================================================
with open(os.path.join(opt.salida, "A8_candidatos.tsv"), "w") as f:
    f.write("puesto\tgen\tcoseno\tgrado\ten_A6\n")
    for puesto, i in enumerate(orden, start=1):
        f.write(f"{puesto}\t{genes[i]}\t{puntuacion[i]:.4f}\t{g.degree(i)}\t{'si' if genes[i] in candidatos_a6 else 'no'}\n")

with open(os.path.join(opt.salida, "A8_validacion.tsv"), "w") as f:
    f.write("gen\tauc_loo\tauc_loo_cercana\n")
    for s, valor, cercana in zip(ids, auc_semillas, auc_semillas_cercana):
        f.write(f"{genes[s]}\t{valor:.4f}\t{cercana:.4f}\n")

resumen = [
    ("genes", n), ("interacciones", g.ecount()), ("semillas", int(es_semilla.sum())),
    ("interacciones_escondidas", len(pares_pos)), ("auc_aristas_embedding", round(auc_aristas, 3)),
    ("auc_aristas_jaccard", round(auc_jaccard, 3)), ("auc_semillas_embedding", round(auc_loo, 3)),
    ("auc_semillas_cercana", round(auc_cercana, 3)),
    ("auc_semillas_rwr_A6", "NA" if auc_rwr is None else auc_rwr), ("auc_semillas_grado", round(auc_grado, 3)),
    ("modulos_louvain", k), ("nmi_kmeans_louvain", round(nmi, 3)), ("nmi_azar_media", round(float(np.mean(nmi_azar)), 3)),
    ("candidatos_comunes_A6", sum(genes[i] in candidatos_a6 for i in orden)),
]
with open(os.path.join(opt.salida, "A8_resumen.tsv"), "w") as f:
    f.write("medida\tvalor\n")
    for medida, valor in resumen:
        f.write(f"{medida}\t{valor}\n")

# =============================================================
# ========================= 6. Figura =========================
# =============================================================
# PCA de los vectores a 2 dimensiones, solo para dibujar: las distancias se miden siempre con los vectores completos, no en el dibujo
xy = PCA(n_components=2, random_state=123).fit_transform(X)
plt.figure(figsize=(7, 6))
plt.scatter(xy[~es_semilla, 0], xy[~es_semilla, 1], s=8, c="grey", alpha=0.5, label="Vecinos (no HPO)")
plt.scatter(xy[es_semilla, 0], xy[es_semilla, 1], s=30, c="red", label="Semillas (HPO)")
for i in ids:
    plt.annotate(genes[i], xy[i], fontsize=6)
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("Embedding de la subred de 700 (PCA de los vectores)")
plt.legend()
# Sin fecha de creación en el PDF, para que dos ejecuciones den el mismo archivo
plt.savefig(os.path.join(opt.salida, "A8_pca.pdf"), metadata={"CreationDate": None})
plt.close()

# Resumen para la salida de launch.sh
print(f"Subred: {n} genes, {g.ecount()} interacciones, {int(es_semilla.sum())} semillas")
print(f"Interacciones escondidas ({len(pares_pos)}): AUC con el embedding {auc_aristas:.3f}; con Jaccard {auc_jaccard:.3f}")
print(f"Semillas escondidas (leave-one-out): AUC con el embedding {auc_loo:.3f} (media de las semillas) y {auc_cercana:.3f} (semilla más cercana); con el random walk de la A6 {auc_rwr}; solo con el grado {auc_grado:.3f}")
print(f"k-means frente a Louvain ({k} módulos): NMI {nmi:.3f}; al azar {np.mean(nmi_azar):.3f} de media")
print("Candidatos:", ", ".join(genes[i] for i in orden))
