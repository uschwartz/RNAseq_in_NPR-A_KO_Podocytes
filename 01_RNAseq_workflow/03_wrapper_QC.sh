echo 'starting QC for project'  $1

ID=$1
GTF=$2
R_script=$3

cd $ID


for entry in *_*; do

echo 'QC for sample' $entry

### set limit of threads
ulimit -n 10000

mkdir $entry/quali-report
qualimap rnaseq --java-mem-size=8G -a uniquely-mapped-reads -bam $entry/*.bam -gtf $GTF -outdir $entry/quali-report

mkdir $entry/dup_report

java -Xmx20g -jar ~/Tools/picard.jar MarkDuplicates M=$entry/dup_report/dupstats.txt REMOVE_DUPLICATES=FALSE I=$entry"/Aligned.sortedByCoord.out.bam" O=$entry/dup_report/$entry.Aligned.out.dupmark.bam

Rscript $R_script $entry/dup_report $GTF

samtools idxstats $entry/Aligned.sortedByCoord.out.bam >$entry/chr_stats.txt

done



