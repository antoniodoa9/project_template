# Imagen para reproducir todo el análisis con Docker (alternativa a ejecutar setup.sh y launch.sh directamente).
# Construir (desde la raíz del repositorio):  docker build -t proyecto_hp0002647 .
# Ejecutar (los resultados quedan en results/ del ordenador):
#   docker run --rm --user "$(id -u):$(id -g)" -v "$(pwd)/results:/proyecto/results" proyecto_hp0002647

# R 4.6.2 sobre Ubuntu, imagen oficial del proyecto Rocker (versión de R fija)
FROM rocker/r-ver:4.6.2

# Librerías del sistema que necesitan los paquetes de R para compilarse (las mismas que el apt install del ejemplo
# del profesor, más las de igraph, XML y las fuentes de las figuras), python3 para la A8 y curl para get-pip.py
RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev libssl-dev libuv1-dev libcairo2-dev libxml2-dev \
        libglpk-dev libgmp-dev libfontconfig1-dev libfreetype-dev libharfbuzz-dev libfribidi-dev \
        libpng-dev libjpeg-dev libtiff-dev libwebp-dev zlib1g-dev libbz2-dev liblzma-dev \
        python3 curl ca-certificates \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /proyecto

# Primero solo setup.sh, para que Docker guarde en caché los paquetes instalados:
# si después cambia un script de análisis, no hay que volver a instalar nada
COPY code/setup.sh code/setup.sh
RUN bash code/setup.sh \
    && chmod -R a+rwX /proyecto/software/tmp

# Después el resto del código
COPY code/ code/

# Con --user el contenedor no tiene carpeta personal: se le da /tmp, en la que puede escribir cualquier usuario,
# para R y para la configuración de matplotlib (A8)
ENV HOME=/tmp \
    MPLCONFIGDIR=/tmp/matplotlib

# Al arrancar el contenedor se ejecuta el análisis completo; results/ se conecta con el ordenador con -v
CMD ["bash", "code/launch.sh"]
