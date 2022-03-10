# submit GATK processing and variant calling file

samplefiles= YOUR_PATH/star_merged_bams/star_merged_bams.txt


for file in `cat $samplefiles`;do
	bsub -n 1 -e starmergedGATK.err -o starmergedGATK.out ./bsub_GATK_variantcalls_starmergedFiles.sh $file
done
