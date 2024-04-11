#DianaOaxaca
#Remove primers

mkdir -p results/02.cutadapt
out="results/02.cutadapt"
FASTQ=$(ls data/*.gz | sed 's/_.*//' | sed 's/data\///' | sort -u)

date

for FILE in ${FASTQ[@]}; do
	echo -e "Run cutadapt to $FILE sample"
        CUTADAPT='cutadapt -m 200 --pair-filter any --no-indels -g CCTACGGGNGGCWGCAG -G GACTACHVGGGTATCTAATCC -Z -j 80 -o '$out'/'$FILE'_1.fastq.gz -p '$out'/'$FILE'_2.fastq.gz data/'$FILE'_R1.fastq.gz  data/'$FILE'_R2.fastq.gz'
        echo -e $CUTADAPT "\n"
        $CUTADAPT
done

date
