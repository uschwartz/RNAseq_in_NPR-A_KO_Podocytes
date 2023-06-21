echo 'multiqc in path'  $1

ID=$1


cd $ID
mkdir multiQC

multiqc -o multiQC -d .
