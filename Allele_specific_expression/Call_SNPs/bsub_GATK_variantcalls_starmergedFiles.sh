#CHECK THAT DEDUPPED FILE EXISTS! IF NOT, EXIT! and THROW ERROR!!
#otherwise continue, and delete aligned_out.sam first (also, genome and anything else if its left) 

file=$1
genomeFasta=YOUR_PATH/GRCh38_chrYPAR_masked.fa #FASTA FILE WITH Y PAR MASKED
scratch_space=YOUR_PATH/alignments_PARY_mask
vcfParent=YOUR_PATH/starmerged_VCFs
dbsnp=YOUR_PATH/dbSNP/common_all_20180418_chr.vcf.gz #DOWNLOAD HERE: ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/common_all_20180418.vcf.gz and update chromosome names for GATK compatibility
vcfdir=$vcfParent/$file
mkdir $vcfdir
echo $vcfdir
#mkdir $scratch_space/$file

if [ -f $scratch_space/$file ];then
       java -jar /usr/local/share/picard-tools/picard.jar AddOrReplaceReadGroups I=$scratch_space/$file \
O=$scratch_space/$file.rg_added_sorted.bam SO=coordinate RGID=1 RGLB=$file RGPL=illumina RGPU=some_machine RGSM=$file
fi

##only dedup if readgroup file exists
if [ -f $scratch_space/$file.rg_added_sorted.bam ];then
	java -jar /usr/local/share/picard-tools/picard.jar MarkDuplicates I=$scratch_space/$file.rg_added_sorted.bam O=$scratch_space/$file.dedupped.bam \
CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics
fi

##SplitNtrim
if [ -f $scratch_space/$file.dedupped.bam ] && [ ! -f $scratch_space/$file.split.bam ];then
        YOUR_PATH/gatk-4.1.2.0/gatk SplitNCigarReads -R $genomeFasta \
-I $scratch_space/$file.dedupped.bam -O $scratch_space/$file.split.bam #mapping quality remap is now default?
fi


###########################

#base recalibation
#https://software.broadinstitute.org/gatk/documentation/article?id=2801
if [ -f $scratch_space/$file.split.bam ] && [ ! -f $scratch_space/$file.recal_data.table ];then
	YOUR_PATH/gatk-4.1.2.0/gatk BaseRecalibrator -R $genomeFasta -I $scratch_space/$file.split.bam --known-sites $dbsnp -O $scratch_space/$file.recal_data.table
fi

#https://gatkforums.broadinstitute.org/gatk/discussion/23296/second-pass-base-recalibration-gatk4

#apply recalibration
if [ -f $scratch_space/$file.recal_data.table ] && [ ! -f $scratch_space/$file.recal_reads.bam ];then
	YOUR_PATH/gatk-4.1.2.0/gatk ApplyBQSR -R $genomeFasta -I $scratch_space/$file.split.bam -bqsr $scratch_space/$file.recal_data.table -O $scratch_space/$file.recal_reads.bam
fi

samtools index $scratch_space/$file.recal_reads.bam

##call variants
YOUR_PATH/gatk-4.1.2.0/gatk HaplotypeCaller -R $genomeFasta -I $scratch_space/$file.recal_reads.bam --dont-use-soft-clipped-bases -stand-call-conf 20.0 -O $vcfdir/output.g.vcf.gz -ERC GVCF


##filter
YOUR_PATH/gatk-4.1.2.0/gatk VariantFiltration -R $genomeFasta -V $vcfdir/output.g.vcf.gz -window 35 -cluster 3 -filter-name "FSgreaterthan30" -filter "FS > 30.0" \
-filter-name "QDlessthan2" -filter "QD < 2.0" -O $vcfdir/filtered_genotyped_output.vcf
vcftools --vcf $vcfdir/filtered_genotyped_output.vcf --remove-filtered-all --recode --recode-INFO-all --out $vcfdir/filtered_output

bgzip -f $vcfdir/filtered_output.recode.vcf
tabix -f -p vcf $vcfdir/filtered_output.recode.vcf.gz
vcffilter -g "GT = 0/1" -r "chrX" $vcfdir/filtered_output.recode.vcf.gz | vcffixup - | vcffilter -f "AC > 0 & DP > 4" > $vcfdir/chr_X_hets.vcf
YOUR_PATH/gatk-4.1.2.0/gatk VariantAnnotator -V $vcfdir/chr_X_hets.vcf -O $vcfdir/annotated_chrX_hets.vcf --dbsnp $dbsnp
YOUR_PATH/gatk-4.1.2.0/gatk SelectVariants -select-type SNP -V $vcfdir/annotated_chrX_hets.vcf -O $vcfdir/annotated_chrX_SNPs_filtered.vcf
bgzip -f $vcfdir/annotated_chrX_SNPs_filtered.vcf
tabix -f -p vcf $vcfdir/annotated_chrX_SNPs_filtered.vcf.gz
