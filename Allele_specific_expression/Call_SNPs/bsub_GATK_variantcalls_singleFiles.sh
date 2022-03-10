#CHECK THAT DEDUPPED FILE EXISTS! IF NOT, EXIT! and THROW ERROR!!
#otherwise continue, and delete aligned_out.sam first (also, genome and anything else if its left) 

file=$1
genomeFasta=YOUR_PATH/GRCh38_chrYPAR_masked.fa
orig_location=YOUR_PATH/alignments_PARY_mask
scratch_space=YOUR_PATH/single_file_alignment_processing/
dbsnp=YOUR_PATH/dbSNP/common_all_20180418_chr.vcf.gz #DOWNLOAD HERE: ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/common_all_20180418.vcf.gz and update chromosome names for GATK compatibility

#if [ -f $scratch_space/$file ];then
#java -jar /usr/local/share/picard-tools/picard.jar AddOrReplaceReadGroups I=$orig_location/$file/Aligned.out.bam \
#O=$scratch_space/$file.rg_added_sorted.bam SO=coordinate RGID=1 RGLB=$file RGPL=illumina RGPU=some_machine RGSM=$file
#fi

##only dedup if readgroup file exists
#if [ -f $scratch_space/$file.rg_added_sorted.bam ];then
#	java -jar /usr/local/share/picard-tools/picard.jar MarkDuplicates I=$scratch_space/$file.rg_added_sorted.bam O=$scratch_space/$file.dedupped.bam \
#CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics
#fi

##SplitNtrim
#if [ -f $scratch_space/$file.dedupped.bam ] && [ ! -f $scratch_space/$file.split.bam ];then
	YOUR_PATH/gatk-4.1.2.0/gatk SplitNCigarReads -R $genomeFasta \
-I $scratch_space/$file.dedupped.bam -O $scratch_space/$file.split.bam #mapping quality remap is now default?
#fi


###########################

#base recalibation
#https://software.broadinstitute.org/gatk/documentation/article?id=2801
#if [ -f $scratch_space/$file.split.bam ] && [ ! -f $scratch_space/$file.recal_data.table ];then
	YOUR_PATH/gatk-4.1.2.0/gatk BaseRecalibrator -R $genomeFasta -I $scratch_space/$file.split.bam --known-sites $dbsnp -O $scratch_space/$file.recal_data.table
#fi
#their docs are THE WORST
#https://gatkforums.broadinstitute.org/gatk/discussion/23296/second-pass-base-recalibration-gatk4

#apply recalibration
#if [ -f $scratch_space/$file.recal_data.table ] && [ ! -f $scratch_space/$file.recal_reads.bam ];then
        YOUR_PATH/gatk-4.1.2.0/gatk ApplyBQSR -R $genomeFasta -I $scratch_space/$file.split.bam -bqsr $scratch_space/$file.recal_data.table -O $scratch_space/$file.recal_reads.bam
#fi

samtools index $scratch_space/$file.recal_reads.bam
