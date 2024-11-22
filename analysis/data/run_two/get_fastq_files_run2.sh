#!/bin/bash
#obtener solo fastq de sedimentos
#con ayuda de chatgpt
dir="/axolote/mvazquez/DATOS_MVMM_MirnaVRL_2/Data/Intensities/BaseCalls"

ids_file="sediment_samples.txt"

while IFS= read -r ids
do
  for file in "$dir"/$ids*; do
    # Verificar si el archivo existe para evitar intentar enlazar archivos no existentes.
    if [ -e "$file" ]; then
      # Extraer parte del nombre del archivo para determinar si es R1 o R2.
      if [[ "$file" =~ _R1_ ]]; then
        type="R1"
      elif [[ "$file" =~ _R2_ ]]; then
        type="R2"
      fi

      # Crear liga simbólica con el nuevo nombre.
      #ls -lh  "$file"
      ln -s "$file" "./${ids}${type}.fastq.gz"
    fi
  done
done < "$ids_file"
