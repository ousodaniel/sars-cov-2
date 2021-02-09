set -e
echo "loading modules..."
module load trimmomatic/0.38 samtools/1.11 bwa/0.7.17 python/3.7 bedtools/2.29.0 ivar/1.3 \
sratoolkit/2.10.0 #bcftools/1.8

###################################################################################################
echo "indexing refcov genome..."
#### DB Indexing
#bwa: index genome
bwa index -a bwtsw -p refcov refcov.fa

echo "finished genome indexing"
echo "########"

###################################################################################################
echo "trimming, aligning & sorting files..."
#### Trimming(trm), Aligning(aln) & Sort(srt)
for f in *_1.fastq
do
	base=$(basename ${f} _1.fastq)
	#trimmomatic: trim seq adapters
	#trimmomatic PE ${f} ${base}_2.fastq \
	#${base}_1.trm1.fastq ${base}_1un.trm1.fastq \
	#${base}_2.trm1.fastq ${base}_2un.trm1.fastq \
	#ILLUMINACLIP:NexteraPE-PE.fa:2:40:15

	#bwmem: align reads to indexed genome (output .sam)
	#sview: convert aln-sam to aln-bam (output .bam)
	#ssort: srt aln-bam (output .bam)
	bwa mem -tM 4 refcov ${base}_1.fastq \
	${base}_2.fastq | samtools view -hu \
	| samtools sort -o ${base}_trm1.aln.srt.bam

	#itrim: (input aln-srt.bam)trim pcr primers (requres aln-srt & primer bed)
	#(output is .bam)
	ivar trim -m 20 -i ${base}_trm1.aln.srt.bam -b data/primers/nCoV-2019.bed \
	-p ${base}_trm2
	
	#ssort: (input is aln.bam)sort trm2 bam (output is .bam)
	samtools sort -o ${base}_trm2.srt.bam ${base}_trm2.bam

	echo "finished adapter and primer trimming!"
	echo "########"

###################################################################################################
	echo "computing read depths before & after primer trim..."
	#### Compute the read depths at each position for comparison
	#sdepth: (input is aln-srt.bam)compare depths after primer trim (output is .tsv)
	samtools depth -aa ${base}_trm1.aln.srt.bam > ${base}_trm1.dph

	#sdepth: (input is aln-srt.bam)compare depths after primer trim (output is .tsv)
	samtools depth -aa ${base}_trm2.srt.bam > ${base}.trm2.dph

	echo "finished computing depth!"
	echo "########"

###################################################################################################
	echo "identifying mismatched primers..."
	#### ID mismatched primer to consensus
	#smpileup: (input is .bam, no ref here)pile mrg reps (output is .bam)
	#iconsensus: (input is .bam, piped) create consensus from reps mrg
	#(output is .fasta cons seq & .txt ave base Q)
	#disabled BAQ (-B) because there is no ref
	samtools mpileup -AB -d 0 -Q 20 ${base}_trm2.srt.bam \
	| ivar consensus -p ${base}_trm2.cnscov

	#bindex: (input .fa)index consensus (output mtlp)
	bwa index -a bwtsw -p ${base}_trm2.cnscov ${base}_trm2.cnscov.fa
	
	#bmem: (input .fa)aln primers to consensus ref (output .sam, piped)
	#sview|ssort: convert aln-sam to aln-bam | sort aln-bam
	bwa mem -k 5 -T 16 ${base}_trm2.cnscov data/primers/nCoV-2019.fa \
	| samtools view -hb -F 4 | samtools sort -o \
	data/primers/${base}_nCoV-2019.cns.aln.srt.bam
	
	#smpileup: pileup primers against cns ref
	#ivariants: call primer variants from consensus (output .tsv)
	samtools mpileup -A -d 0 -Q 20 --reference ${base}_trm2.cnscov.fa \
	data/primers/${base}_nCoV-2019.cns.aln.srt.bam | ivar variants \
	-p data/primers/${base}_nCoV-2019.cns.var -t 0.25
	
	#get indices for mismatched primers
	bedtools bamtobed -i data/primers/${base}_nCoV-2019.cns.aln.srt.bam \
	> data/primers/${base}_nCoV-2019.cns.aln.srt.bed
	
	#get primers with mismatches to cns ref (output .txt)
	ivar getmasked -i data/primers/${base}_nCoV-2019.cns.var.tsv -b \
	data/primers/${base}_nCoV-2019.cns.aln.srt.bed -f \
	data/primers/nCoV-2019.tsv -p data/primers/${base}_nCoV-2019.msk

	echo "finished step two!"
	echo "########"

###################################################################################################
	echo "trimming reads with mismatched primers..."
	#### Trim reads from mismatched primers

	#iremovereads: (input aln-srt-trm2(ivar trm) .bam)remove reads 
	#associated with mismatched primer indeces (why not on the merged?) (output .bam)*_rep${c}.trm2.srt.bam
	ivar removereads -i ${base}_trm2.srt.bam -t data/primers/${base}_nCoV-2019.msk.txt \
	-b data/primers/${base}_nCoV-2019.cns.aln.srt.bed \
	-p ${base}_msk.trm3
	
	#sort trm3 bam
	samtools sort -o ${base}_msk.trm3.srt.bam ${base}_msk.trm3.bam
	
	#sdepth: (input is aln-srt.bam)compare depths after primer trim (output is .tsv)
	samtools depth -aa ${base}_msk.trm3.srt.bam > ${base}_trm3.dph

	echo "finished step three!"
	echo "########"

####################################################################################################
	echo "calling variants..."
	#### Variant Calling

	#smpileup: pile trim3 bam
	samtools mpileup -A -d 0 -Q 20 --reference refcov.fa \
	${base}_msk.trm3.srt.bam | ivar variants -m 10 -p ${base}_variants -t 0.25 \
	-r refcov.fa -g refcov.gff

	#ifiltervariants: (input .tsv)filter variants (output .tsv)
	ivar filtervariants -p ${base}_flt.variants ${base}_variants.tsv
done
echo "########"
echo "finished the pipeline!"

###################################################################################################

#ghapcall: (input is .bam) use gatk haplotypecaller to call variants
	gatk HaplotypeCaller -ploidy 1 -I ${base}_msk.trm3.bam \
	-R refcov.fa -O ${base}_g.vcf.gz  