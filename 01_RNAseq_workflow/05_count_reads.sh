echo 'starting counting for project'  $1

ID=$1
GTF=$2

cd $ID

mkdir -p counts
mkdir -p bams
mkdir -p counts/bams

cd $ID

for entry in alignment/*_*; do

echo $entry
mv $entry/dup_report/*bam  bams

done

cd counts

featureCounts -T 10 -s 2 -a $GTF -o count_table.txt ../bams/*.bam &>count_info.txt

mkir multiQC
multiQC -o multiQC .



