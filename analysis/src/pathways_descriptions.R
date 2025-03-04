# Instalar y cargar el paquete dplyr si no está instalado
if (!requireNamespace('dplyr', quietly = TRUE)) install.packages('dplyr', dependencies=TRUE)
library(dplyr)

# Definir los argumentos de línea de comandos
args <- commandArgs(trailingOnly = TRUE)

# Verificar que se proporcionen los nombres de archivo de entrada y salida
if (length(args) != 4) {
  stop("Uso: script.R -i archivo_base -o nombre_archivo_output")
}

# Obtener los nombres de archivo de entrada y salida desde los argumentos
archivo_base <- args[2]
archivo_output <- args[4]

# URL del archivo que ingresas manualmente (comprimido)
#url_manual <- "https://github.com/picrust/picrust2/raw/master/picrust2/default_files/description_mapfiles/metacyc_pathways_info.txt.gz"
url_manual <- "https://github.com/picrust/picrust2/raw/master/picrust2/default_files/description_mapfiles/ko_info.tsv.gz"
# Descargar el archivo .gz localmente
temp_file <- tempfile(fileext = ".gz")
download.file(url_manual, temp_file, mode = "wb")

# Leer el archivo comprimido con read.table
datos_manual <- read.table(gzfile(temp_file), header=FALSE, sep="\t", quote="", comment.char="")
colnames(datos_manual) <- c("#OTU ID", "description")

# Cargar datos desde el archivo base
datos_base <- read.delim(archivo_base, check.names = F, sep="\t", skip = 1)

# Realizar el join con dplyr
resultado <- dplyr::inner_join(datos_manual, datos_base)

# Imprimir el resultado
print(resultado)

# Guardar el resultado en un nuevo archivo
write.table(resultado, archivo_output, sep='\t', quote=FALSE, row.names=FALSE)

# Limpiar el archivo temporal
unlink(temp_file)
