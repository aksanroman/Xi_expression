firstPassDir=/YOUR_PATH
genomeDir=YOUR_PATH/STAR/hg38_STAR_PARYMASK

samplename=${1}
read1_file=${2}
read2_file=${3}
runDir=$firstPassDir/$samplename
echo $runDir
echo ${read1_file} 
echo ${read2_file}
mkdir $runDir
cd $runDir

#readFilesIn                 Read1 Read2
#for paired end: 
STAR --genomeDir $genomeDir --readFilesIn ${read1_file} ${read2_file} --runThreadN 8 -twopassMode Basic --readFilesCommand zcat --outSAMtype BAM Unsorted
#small subset of the files are not gzipped 
#STAR --genomeDir $genomeDir --readFilesIn ${read1_file} ${read2_file} --runThreadN 8 -twopassMode Basic --outSAMtype BAM Unsorted
#|samtools view -bo $runDir/aligned.bam

#STAR --genomeDir $genomeDir --readFilesIn $file --runThreadN 8 -twopassMode Basic --readFilesCommand zcat
#--readFilesCommand zcat so that i dont have to gunzip each file.... 
#command rm Aligned.out.sam
