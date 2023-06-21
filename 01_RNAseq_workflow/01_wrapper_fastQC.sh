
echo 'starting FastQC for project'  $1

cd $1

mkdir -p FastQC

for entry in *.gz; do

echo 'fastqc for sample' $entry

fastqc -o FastQC  --extract $entry

done

multiqc  -o FastQC FastQC/*.zip 
