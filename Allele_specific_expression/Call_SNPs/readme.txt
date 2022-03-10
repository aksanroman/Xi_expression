How to run these scripts: 

Step A: create VCFs for individuals: 
1. Align your raw data (./parallel_alignSTAR_twopassmode.sh)

2. Merge all bams from each individual: use makeMergeScript_starbams.py to generate the mergebams script: mergebams_star_Script.sh (example output provided)

3. Run GATK pipe on merged bams to call variants from RNAseq and filter using dbsnp for expressed het snps (./parallel_starmergedGATK.sh). 

Step B: process individual sample bams and run ASE
4. Run GATK pipe on individual bams to get sample-level processed bams for calling alelle specific expression (./parallel_GATKvariants.sh)

5. Using filtered VCF and individual processed bams, run ASE (./ASE_process.sh)

