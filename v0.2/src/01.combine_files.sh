#!/bin/bash
#Concat fastq files
#DianaOaxaca

#define dirs
dir1="data/run_one"
dir2="data/run_two"
output_dir="data"

# Concat
for file in "$dir1"/*.gz; do
    filename=$(basename "$file")
    if [ -e "$dir2/$filename" ]; then
        #ls "$dir1/$filename" "$dir2/$filename"
	cat "$dir1/$filename" "$dir2/$filename" > "$output_dir/$filename"
        echo "Combine: $output_dir/$filename"
    else
        echo "There isn't $filename in $dir2, skip."
    fi
done

echo "Finish concat!"
