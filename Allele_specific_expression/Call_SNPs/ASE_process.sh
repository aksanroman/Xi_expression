# process alele specific expression at all embryo-level het snps for each cell 

samplefiles="YOUR_PATH/ASE_filemap.txt"


#for individual file in `cat $samplefiles`;do
#for line in `cat $samplefiles`;do
while read -r line;do
#	echo $line
	individual=$(echo $line | cut -f1 -d" ")
	file=$(echo $line | cut -f2 -d" ")
	bsub -n 1 -e ASE.err -o ASE.out ./bsub_ASE.sh $file $individual
done < $samplefiles


