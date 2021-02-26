# sars-cov-2
This folder is associated with Sars-CoV-2 variant calling and annotation. It can be downloaded as a semi-standalone to perform analyses on a distributed cluster, based on `module` package software module manager, installed with the required modules. It was tested on the [ILRI cluster](https://hpc.ilri.cgiar.org/using-the-cluster) and its functioning cannot be garanteed elsewhere, especially, without adjusting or setting variables accordingly. It is still under active development.

## Pipeline
The pipeline script was created based on [gkarthik's](https://github.com/andersen-lab/paper_2018_primalseq-ivar/blob/master/cookbook/CookBook.ipynb) iVar works. It entails variant calling all the way to annotation. It is still in active development.

## Directory Organisation
### slurm.sbatch
`SLURM` scheduler script: for running the variant calling script (`sbatch -w <computeNode> slurm.sbatch`) from `pwd`. You will have to modify the email field approriately.

### adapter
 - Trimming adapter `.fa` files.

### exe
 - [**ivar_variants_to_vcf.py**](https://raw.githubusercontent.com/nf-core/viralrecon/master/bin/ivar_variants_to_vcf.py): two copies, one (`*_s_*`) of which is modified to filter out only the Spike gene variants.
 - **ncov19vc.sh**: variant calling script based on `slurm job scheduler`. Be sure to modify read file **suffixes** and **input**/**output** directories accordingly. Outputs of the script are dumped in the `pwd` while script is running, but moved to respective output directories at the end of every sample processing cycle. Some outputs of interest are copied to separate directories for consolidation, especially for downstream analyses and reporting on local client.
 - **snpEff dir**: `snpEff` _installation_ configured with sars-cov-2 references 

### output
 - Will contain samples analyses output organised in directories per sample, including assssociated `slurm.out` files

### primer
 - Artic `V3` primer `.bed`, `.fa` and `.tsv` files

### ref
 - Sars-CoV-2 reference genome [**NC_045512.2**](https://www.ncbi.nlm.nih.gov/genome/?term=NC_045512.2) `.fa` sequence and `.gff` feature files

### sra
 - Contains read file pairs (`.fatq`/`.fastq.gz`)
