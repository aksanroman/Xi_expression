# submit GATK processing and variant calling file

samplefiles=YOUR_PATH/samplenames.txt

for file in `cat $samplefiles`;do
	bsub -n 1 -e processSingleFilespicardGATK.err -o processSingleFilespicard.out ./bsub_GATK_variantcalls_singleFiles.sh $file
	sleep 1
done
