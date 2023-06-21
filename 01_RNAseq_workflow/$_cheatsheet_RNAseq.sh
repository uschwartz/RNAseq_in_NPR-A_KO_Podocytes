AnalysisPath="/Users/admin/Analysis/17_20191202_Mm_RNAseq"
AnnoPath="/Users/admin/Annotation/MusMusculus"

############### Genome Index ####################
cd $AnnoPath

STAR --runThreadN 12 --runMode genomeGenerate --genomeDir STARidx --genomeFastaFiles Mus_musculus.GRCm38.dna.primary_assembly.fa --sjdbGTFfile Mus_musculus.GRCm38.98.gtf


################################## FastQC ###############################
$AnalysisPath/script/01_wrapper_fastQC.sh $AnalysisPath/data


################################## alignment ###############################
$AnalysisPath/script/02_wrapper_STAR_align.sh $AnalysisPath/data $AnnoPath"/STARidx"




################################## QC ###############################
$AnalysisPath/script/03_wrapper_QC.sh  $AnalysisPath/alignment $AnnoPath/protein_coding.gtf $AnalysisPath/script/support_scripts/dupRadar_script.R



$AnalysisPath/script/04_multiqc.sh  $AnalysisPath/alignment


################################## counting  ###############################
$AnalysisPath/script/05_count_reads.sh  $AnalysisPath $AnnoPath"/protein_coding_and_lincRNA.gtf"


################################## statistics ###############################
Rscript $AnalysisPath/script/06_count_statistics.R $AnalysisPath/counts





