# ind_comparisons_trop_cal_mello_2024

I am going to cal FST, Admixture and PI individual wise - removing 1-20mb in chr7 for trop individuals using bam files

# FST

working in beluga

first I copied all the bamfiles and bai files into the directory

Then created a file with all bam file names

```bash
ls ../*.bam>all_samples.txt
```

Then I did spli in into different files including only 1 sample in each file and naming then pop1 , pop2 etc using following command

```bash
split --numeric-suffixes=1 -l 1 all_samples.txt "pop"
```
This adds a 0 in front of numbers 1-9. This complicates using arrays. Therefore i am removing 0 in file names from 1 to 9

```bash
for i in {1..9} ; do mv pop0${i} pop${i}; done
```

Now I want to remove Chr7 1-20 mb region. I am creating a file with information on what to include in the analysis

```txt
Chr1:                   \\all of chr1
Chr2:
Chr3:
Chr4:
Chr5:
Chr6:
Chr7:20000000-          \\chr7 but exclude the first 20000000 bases
Chr8:
Chr9:
Chr10:
```

# Admixture

Made a sub folder inside previous folder with bam files and started Admixture there

first created the text file bam file list

bam.filelist

```txt
../F_Ghana_WZ_BJE4687_combined__sorted.bam_rg_rh.bam
../F_IvoryCoast_xen228_combined__sorted.bam_rg_rh.bam
../F_Nigeria_EUA0331_combined__sorted.bam_rg_rh.bam
../F_Nigeria_EUA0333_combined__sorted.bam_rg_rh.bam
../F_SierraLeone_AMNH17272_combined__sorted.bam_rg_rh.bam
../F_SierraLeone_AMNH17274_combined__sorted.bam_rg_rh.bam
../all_ROM19161_sorted.bam
../XT10_WZ_no_adapt._sorted.bam_rg_rh.bam
../XT11_WW_trim_no_adapt_scafconcat_sorted.bam_rg_rh.bam
../M_Ghana_WY_BJE4362_combined__sorted.bam_rg_rh.bam
../M_Ghana_ZY_BJE4360_combined__sorted.bam_rg_rh.bam
../M_Nigeria_EUA0334_combined__sorted.bam_rg_rh.bam
../M_Nigeria_EUA0335_combined__sorted.bam_rg_rh.bam
../M_SierraLeone_AMNH17271_combined__sorted.bam_rg_rh.bam
../M_SierraLeone_AMNH17273_combined__sorted.bam_rg_rh.bam
../XT1_ZY_no_adapt._sorted.bam_rg_rh.bam
../XT7_WY_no_adapt__sorted.bam_rg_rh.bam
```
then prepared the files for Admixture analysis with

prep_admix.sh

```bash
#!/bin/sh
#SBATCH --job-name=fst
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=512gb
#SBATCH --output=abba.%J.out
#SBATCH --error=abba.%J.err
#SBATCH --account=def-ben

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

# Load modules
module load nixpkgs/16.09  gcc/7.3.0
module load angsd
module load gsl/2.5
module load htslib

angsd -GL 1 -out genolike -nThreads 10 -doGlf 2 -doMajorMinor 1 -doMaf 0 -doCounts 1 -setMinDepthInd 5 -setMaxDepth 100 -minMapQ 20 -bam bam.filelist
```
cal Admix

```bash
#!/bin/sh
#SBATCH --job-name=bwa_505
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=168:00:00
#SBATCH --mem=256gb
#SBATCH --output=bwa505.%J.out
#SBATCH --error=bwa505.%J.err
#SBATCH --account=def-ben
#SBATCH --array=2-5

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

# Load modules
module load nixpkgs/16.09  gcc/7.3.0
module load angsd
module load gsl/2.5
module load htslib

NGSadmix -likes genolike.beagle.gz -K ${SLURM_ARRAY_TASK_ID} -P 4 -o myoutfiles${SLURM_ARRAY_TASK_ID} -minMaf 0.05
```
# PI/Nucleotide diversity for individuals

```bash
#!/bin/sh
#SBATCH --job-name=fst
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=64gb
#SBATCH --output=abba.%J.out
#SBATCH --error=abba.%J.err
#SBATCH --account=def-ben
#SBATCH --array=1-20

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

# Load modules
module load nixpkgs/16.09  gcc/7.3.0
module load angsd
module load gsl/2.5
module load htslib

angsd -bam ind${SLURM_ARRAY_TASK_ID} -doMaf 0 -doCounts 1 -setMinDepthInd 5 -setMaxDepth 100 -minMapQ 20 -doSaf 1 -anc ../../reference_genome/XENTR_10.0_genome_scafconcat_goodnamez.fasta -GL 1 -P 24 -out out_${SLURM_ARRAY_TASK_ID}

realSFS out_${SLURM_ARRAY_TASK_ID}.saf.idx -P 24 > out_${SLURM_ARRAY_TASK_ID}.sfs
realSFS saf2theta out_${SLURM_ARRAY_TASK_ID}.saf.idx -sfs out_${SLURM_ARRAY_TASK_ID}.sfs -outname out_${SLURM_ARRAY_TASK_ID}
thetaStat do_stat out_${SLURM_ARRAY_TASK_ID}.thetas.idx -win 5000000 -step 5000000  -outnames theta.thetasWindow_${SLURM_ARRAY_TASK_ID}.gz
```
