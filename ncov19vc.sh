#!/usr/bin/env bash
#set -e

#directory hierrachy: pwd:adapter; exe; fastqc; output; primer; ref; sra
#scripts are executed from pwd

printf "loading modules...\n\ttrimmomatic/0.38\n\tsamtools/1.11\n\tbwa/0.7.17\n\tpython/3.7\n\tbedtools/2.29.0 \
\n\tivar/1.3\n\tsratoolkit/2.10.0\n"

module load trimmomatic/0.38 samtools/1.11 bwa/0.7.17 python/3.7 bedtools/2.29.0 ivar/1.3 \
sratoolkit/2.10.0 

echo "####################################################################################################"

###################################################################################################
#retrieving slurm output file name
slurm_out=$(ls -t | head -n1)

#initialise sample count
count=0

###################################################################################################
echo "indexing refcov genome..."
#### DB Indexing
#bwa: index genome
bwa index -a bwtsw -p ref/refcov ref/refcov.fa

echo "finished genome indexing"
echo "####################################################################################################"

###################################################################################################
echo "trimming, aligning & sorting files..."
#### Trimming(trm), Aligning(aln) & Sort(srt)
for f in sra/*_R1_001.fastq.gz
do
	count=$(( count + 1 ))

	base=$(basename ${f} _R1_001.fastq.gz)
	#trimmomatic: trim seq adapters
    trimmomatic PE ${f} sra/${base}_R2_001.fastq.gz \
    ${base}_1.trm1.fastq ${base}_1un.trm1.fastq \
    ${base}_2.trm1.fastq ${base}_2un.trm1.fastq \
    #ILLUMINACLIP:adapter/adapters.fa:2:40:15:2:True \
    TRAILING:19 \
    HEADCROP:8

	#bwmem: align reads to indexed genome (output .sam)
	#sview: convert aln-sam to aln-bam (output .bam)
	#ssort: srt aln-bam (output .bam)
	bwa mem -t 6 ref/refcov ${base}_1.fastq \
	${base}_2.fastq | samtools view -hu -F 4 -F 2048\
	| samtools sort -o ${base}_trm1.aln.srt.bam

	#itrim: (input aln-srt.bam)trim pcr primers (requres aln-srt & primer bed)
	#(output is .bam)
	ivar trim -m 20 -i ${base}_trm1.aln.srt.bam -b primer/nCoV-2019.bed \
	-p ${base}_trm2
	
	#ssort: (input is aln.bam)sort trm2 bam (output is .bam)
	samtools sort -o ${base}_trm2.srt.bam ${base}_trm2.bam

	echo "finished adapter and primer trimming!"
	echo "####################################################################################################"

###################################################################################################
	#create output directories for every sample
	dir=$(mkdir -p output/${base}_outputs)
###################################################################################################
	echo "extracting primer read coverage..."
	#cat outputs the content of current slurm output
	#grep finds the primer coverage data (^n[C] pattern), plus its header line (-B)
	#sed replaces spaces with underscore in the header line
	cat $slurm_out | grep -B 1 '^n[C]' | sed '1 s/ /_/g' > ${base}_primers.dph

###################################################################################################
	echo "computing read depths before & after primer trim..."
	#### Compute the read depths at each position for comparison
	#sdepth: (input is aln-srt.bam)compare depths after primer trim (output is .tsv)
	samtools depth -aa ${base}_trm1.aln.srt.bam > ${base}_trm1.dph

	#sdepth: (input is aln-srt.bam)compare depths after primer trim (output is .tsv)
	samtools depth -aa ${base}_trm2.srt.bam > ${base}.trm2.dph

	echo "finished computing depths!"
	echo "####################################################################################################"

###################################################################################################
	echo "identifying mismatched primers..."
	#### ID mismatched primer to consensus
	#smpileup: (input is .bam, no ref here)pile mrg reps (output is .bam)
	#iconsensus: (input is .bam, piped) create consensus from reps mrg
	#(output is .fasta cons seq & .txt ave base Q)
	#disabled BAQ (-B) because there is no ref
	samtools mpileup -A -d 0 -Q 20 ${base}_trm2.srt.bam \
	| ivar consensus -p ${base}_trm2.cnscov

	#bindex: (input .fa)index consensus (output mtlp)
	bwa index -a bwtsw -p ${base}_trm2.cnscov ${base}_trm2.cnscov.fa
	
	#bmem: (input .fa)aln primers to consensus ref (output .sam, piped)
	#sview|ssort: convert aln-sam to aln-bam | sort aln-bam
	bwa mem -k 5 -T 16 ${base}_trm2.cnscov primer/nCoV-2019.fa \
	| samtools view -hu -F 4 | samtools sort -o \
	${base}_nCoV-2019.cns.aln.srt.bam
	
	#smpileup: pileup primers against cns ref
	#ivariants: call primer variants from consensus (output .tsv)
	samtools mpileup -AB -d 0 -Q 20 --reference ${base}_trm2.cnscov.fa \
	${base}_nCoV-2019.cns.aln.srt.bam | ivar variants \
	-p ${base}_nCoV-2019.cns.var -t 0.25
	
	#get indices for mismatched primers
	bedtools bamtobed -i ${base}_nCoV-2019.cns.aln.srt.bam \
	> ${base}_nCoV-2019.cns.aln.srt.bed
	
	#get primers with mismatches to cns ref (output .txt)
	ivar getmasked -i ${base}_nCoV-2019.cns.var.tsv -b \
	${base}_nCoV-2019.cns.aln.srt.bed -f \
	primer/nCoV-2019.tsv -p ${base}_nCoV-2019.msk

	echo "finished identifying mismatched primer reads!"
	echo "####################################################################################################"

###################################################################################################
	echo "trimming reads with mismatched primers..."
	#### Trim reads from mismatched primers

	#iremovereads: (input aln-srt-trm2(ivar trm) .bam)remove reads 
	#associated with mismatched primer indeces (why not on the merged?) (output .bam)*_rep${c}.trm2.srt.bam
	ivar removereads -i ${base}_trm2.srt.bam -t ${base}_nCoV-2019.msk.txt \
	-b ${base}_nCoV-2019.cns.aln.srt.bed -p ${base}_msk.trm3
	
	#sort trm3 bam
	samtools sort -o ${base}_msk.trm3.srt.bam ${base}_msk.trm3.bam
	
	echo "computing read depth after trm3..."
	#sdepth: (input is aln-srt.bam)compare depths after primer trim (output is .tsv)
	samtools depth -aa ${base}_msk.trm3.srt.bam > ${base}_trm3.dph

	echo "finished removing reads with mismatched primers!"
	echo "####################################################################################################"

####################################################################################################
	echo "calling variants..."
	#### Variant Calling

	#smpileup: pile trim3 bam
	samtools mpileup -AB -d 0 -Q 20 --reference ref/refcov.fa \
	${base}_msk.trm3.srt.bam | ivar variants -m 10 -p ${base}_variants -t 0.25 \
	-r ref/refcov.fa -g ref/refcov.gff

	#ifiltervariants: (input .tsv)filter variants (output .tsv)
#	ivar filtervariants -p ${base}_flt.variants ${base}_variants.tsv

	echo "finished variant calling!"
	echo "####################################################################################################"

####################################################################################################
	echo "converting iVar .tsv to .vcf (with spike protein subsetting)..."
    ### Generating .vcf

    #ivar_variants_to_vcf.py: converting variants.tsv to variants.vcf
    python ivar_variants_s_to_vcf.py  -po ${base}_variants.tsv ${base}_variants.vcf

####################################################################################################
	echo "####################################################################################################"
	echo "annotating variants using snpEff..."
	### Annotate variants
	
	#snpeff: (input .vcf) (output ann.vcf)
	snpeff="java -Xmx1g -Xmx8g -jar exe/snpEff/snpEff.jar -c exe/snpEff/snpEff.config"
	${snpeff} covHu.1 ${base}_variants.vcf > ${base}_variants.ann.vcf
	
	echo "finished variant annotation!"
	echo "####################################################################################################"

####################################################################################################
	echo "logging sample analysis report (slurm.out)..."
	
	#filter out iVar-specific warnings and store the slurm.out report for the sample
	cat $slurm_out | grep -vw '^WARNING\|^iVar\|^It' | sed -r '/^\s*$/d' > ${base}_slurm.out

####################################################################################################
	echo "archiving sample analysis output..."
	
	cp *variants.tsv *.vcf /home/douso/variant_files

	cp *.dph /home/douso/depth_files

	#move all output associated with every sample to a respective directory
	mv ${base}* snpEff_genes.txt snpEff_summary.html output/${base}_outputs
	
	echo "finished processing sample ${count} (${base} read pair)"
####################################################################################################
	#clean the slurm.out for the next sample
	cat $slurm_out >> all_slurm.out

	> $slurm_out
	
done
echo "####################################################################################################"
echo "finished processing ${count} samples!"

####################################################################################################

#ghapcall: (input is .bam) use gatk haplotypecaller to call variants
#	gatk HaplotypeCaller -ploidy 1 -I ${base}_msk.trm3.bam \
#	-R refcov.fa -O ${base}_g.vcf.gz  

#../bbmap/bbduk.sh -Xmx8g in1=COVC1849_S1_L001_R1_001.fastq.gz in2=COVC1849_S1_L001_R2_001.fastq.gz \
#out=COVC1849_S1_matched.fq outm=COVC1849_S1_umatched.fq ref=hum.fa.gz k=25 hdist=0 stats=stats.txt