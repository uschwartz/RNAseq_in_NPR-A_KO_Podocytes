echo 'starting alignment for project'  $1

ID=$1
Anno=$2

cd $ID

mkdir -p ../alignment

outDir='../alignment'

### set limit of threads
ulimit -n 10000

for entry in *".fastq.gz"; do

echo 'mapping for sample' $entry


patientID=$(echo $entry| cut -d'_' -f 2-5)

echo $patientID

mkdir $outDir/$patientID/

STAR --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --alignIntronMin 10 --alignIntronMax 1000000 --outFilterMismatchNoverReadLmax 0.04 --runThreadN 12 --outSAMtype BAM SortedByCoordinate --outWigType wiggle --outSAMmultNmax 1 --outMultimapperOrder Random  --outFileNamePrefix $outDir/$patientID/ --genomeDir $Anno --readFilesCommand gunzip -c --readFilesIn $entry

samtools index $outDir/$patientID/*.bam

for wig in $outDir/$patientID/*.wig; do

echo ${wig::${#wig}-3}bw
wigToBigWig  $wig $Anno'/chrNameLength.txt' ${wig::${#wig}-7}$patientID.bw
rm $wig

done


done



