#FOR EACH SAMPLE, MERGE BAMS. EXAMPLE BELOW

#160A
bsub -J 160A -e merge160A_star.err -o merge160A_star.out samtools merge YOUR_PATH/star_merged_bams/160A.bam YOUR_PATH/alignments_PARY_mask/160A_XXXX/Aligned.out.bam
