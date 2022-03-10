#CHECK THAT DEDUPPED FILE EXISTS! IF NOT, EXIT! and THROW ERROR!!
#otherwise continue, and delete aligned_out.sam first (also, genome and anything else if its left) 

file=$1
individual=$2
genomeFasta=YOUR_PATH/GRCh38_chrYPAR_masked.fa

scratch_space=YOUR_PATH/single_file_alignment_processing
bam=$scratch_space/$file.recal_reads.bam

vcfParent=YOUR_PATH/starmerged_VCFs
vcfdir=$vcfParent/$individual.bam

outdir=YOUR_PATH/star_ASE/$file
mkdir $outdir
echo $outdir

#annotated version
#if [ ! -f $outdir/chrX.annotated.allsnps.ASE.output.table ];then
#	YOUR_PATH/gatk-4.1.2.0/gatk ASEReadCounter -R $genomeFasta -I $bam -V $vcfdir/annotated_chrX_SNPs_filtered.vcf.gz -O $outdir/chrX.annotated.allsnps.ASE.output.table
#fi

YOUR_PATH/gatk-4.1.2.0/gatk ASEReadCounter -R $genomeFasta -I $bam -V $vcfdir/annotated_chr8_SNPs_filtered.vcf.gz -O $outdir/chr8.annotated.allsnps.ASE.output.table
