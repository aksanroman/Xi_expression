
#from RNASeq.misc import pp

# Setup
project_root = #DIRECTORY FOR ALIGNED FILES NAME
outroot= #OUTPUT DIRECTORY


########################
# Alignment
########################

#Fastq info
fastq_file_handle = open(project_root+"scripts/"+"Aligned_FILES.txt",'r') #GET ALIGNED FILES
headerVals = fastq_file_handle.readline().rstrip().split("\t")
#print(headerVals)
samples = {}
jobinfo= {}
#print(headerVals)

#Samples
for line in fastq_file_handle:
	line = dict(zip(headerVals,line.rstrip().split("\t")))
	#remove .fastq.gz from "sample_id"
	new_sampleID=line['sample_id'][:4]
	print(line['sample_id']+"\t"+new_sampleID)
	output_file=outroot+new_sampleID
	line['sample_id']=new_sampleID
	#update each sample to contain full path to that bam 
	file=line['filename'];
	#full_file_path=outroot+cell_id+"/Aligned.out.bam";
	jobinfo[line['sample_id']]=output_file
	samples.setdefault(line['sample_id'],[]).append(file);

#Write alignment script
outHandle = open(project_root+"scripts/mergebams_star_Script.sh",'w')
sampleNames = [x for x in samples.keys()]
sampleNames.sort()
for sample in sampleNames:
        print("#"+sample, file=outHandle)
        print("bsub -J "+sample+" -e merge"+sample+"_star.err -o merge"+sample+"_star.out samtools merge "+jobinfo[sample]+".bam "+" ".join(samples[sample]),file=outHandle)
