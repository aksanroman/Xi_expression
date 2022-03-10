#make individual genome directory and then align to this. delete individual genome! 

#TEXT FILE WITH PATHS TO FASTQ FILES
samplefiles= #Add path to your sample file 


#IFS=$'\n'
#for line in `cat $samplefiles`;do
#	echo $line
#  	samplename=$(echo ${line} | cut -f1)
#	file=$(echo ${line} | cut -f2)
#	echo $samplename
#	echo $file
#	bsub -n 8 -e star_twopass.err -o star_twopass.out ./bsub_star_twopassmode.sh $samplename $file
#	./bsub_star_twopassmode.sh $line
#	sleep 60
#done


while IFS= read -r line; do
	name=$(echo ${line} | cut -d " " -f1) >> /dev/null
	R1=$(echo ${line} | cut -d " " -f2) >> /dev/null
	R2=$(echo ${line} | cut -d " " -f3) >> /dev/null

	#call processingScript
	#./processingScript.sh ${name} ${location}
	echo $name
	echo $R1
	echo $R2

	#if alignment file doesnt already exist, submit alignment job  
	#	if [ ! -f #YOUR_PATH/$name/Aligned.out.bam# ] && [ ! -f #YOUR_PATH/$name/recal_reads.bam ];then
	bsub -n 8 -e star_twopass.err -o star_twopass.out ./bsub_star_twopassmode.sh $name $R1 $R2

#		sleep 15
#	fi
done < ${samplefiles}
