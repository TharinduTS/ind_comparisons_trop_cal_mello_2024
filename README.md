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
After downloading the files,

First I renamed files with the format

Filename_sex_popularion

eg: allROM19161_F_Liberia

and then I did plot this data with following R script

```R
library("rstudioapi") 
setwd(dirname(getActiveDocumentContext()$path))

library(ggplot2)
#install.packages("Hmisc")
library(Hmisc)

# Delete previous plot before starting if it is in the same folder*****************


#remove scientific notation
options(scipen=999)

pop_name<-""

file_list<-grep(list.files(path="./individual_pi",full.names = TRUE), pattern='\\PI$', value=TRUE)

print(file_list)

# combine files adding pop and sex from file name
df <- do.call(rbind, lapply(file_list, function(x) cbind(read.table(x), file_name=tail(strsplit(x, "/")[[1]],n=1))))

# give new column names
col_names<-c('info','chr','loc','tW','tP','tF','tH','tL','Tajima','fuf','fud','fayh','zeng','nSites','file_name')
colnames(df) <- col_names

# Split name column into firstname and last name
require(stringr)
df[c( 'ind','sex','pop','info')] <- str_split_fixed(df$file_name, '_', 4)

#remove first 20mb of chr7
twenty_mb_removed<-df[df$chr!="Chr7" | df$loc>20000000,]

#remove JBL sample
twenty_mb_removed<-twenty_mb_removed[twenty_mb_removed$ind!="JBL052",]

# renaming all samples
twenty_mb_removed[twenty_mb_removed == "AMNH17272"] <- "SL_F1"
twenty_mb_removed[twenty_mb_removed == "AMNH17274"] <- "SL_F2"
twenty_mb_removed[twenty_mb_removed == "AMNH17271"] <- "SL_M1"
twenty_mb_removed[twenty_mb_removed == "AMNH17273"] <- "SL_M2"
twenty_mb_removed[twenty_mb_removed == "allROM19161"] <- "LB_F1"
twenty_mb_removed[twenty_mb_removed == "xen228"] <- "IC_F1"
twenty_mb_removed[twenty_mb_removed == "BJE4687"] <- "GH_F1"
twenty_mb_removed[twenty_mb_removed == "BJE4362"] <- "GH_M1"
twenty_mb_removed[twenty_mb_removed == "BJE4360"] <- "GH_M2"
twenty_mb_removed[twenty_mb_removed == "XT10"] <- "LT_F1"
twenty_mb_removed[twenty_mb_removed == "XT11"] <- "LT_F2"
twenty_mb_removed[twenty_mb_removed == "XT1"] <- "LT_M1"
twenty_mb_removed[twenty_mb_removed == "XT7"] <- "LT_M2"
twenty_mb_removed[twenty_mb_removed == "EUA0331"] <- "NG_F1"
twenty_mb_removed[twenty_mb_removed == "EUA0333"] <- "NG_F2"
twenty_mb_removed[twenty_mb_removed == "EUA0334"] <- "NG_M1"
twenty_mb_removed[twenty_mb_removed == "EUA0335"] <- "NG_M2"

twenty_mb_removed[twenty_mb_removed == "cal"] <- "Xcal"
twenty_mb_removed[twenty_mb_removed == "mello"] <- "Xmel"

#change pop names

twenty_mb_removed[twenty_mb_removed == "Cal"] <- "Cameroon"
twenty_mb_removed[twenty_mb_removed == "Mello"] <- "Gabon"


# use this section to set sample order

sample_list<-c('SL_F1',
               'SL_F2',
               'SL_M1',
               'SL_M2',
               'LB_F1',
               'IC_F1',
               'GH_F1',
               'GH_M1',
               'GH_M2',
               'LT_F1',
               'LT_F2',
               'LT_M1',
               'LT_M2',
               'NG_F1',
               'NG_F2',
               'NG_M1',
               'NG_M2',
               'Xcal',
               'Xmel')

#use sample list order as levels

twenty_mb_removed$ind=factor(twenty_mb_removed$ind,levels = sample_list)

ND_plot<-ggplot(twenty_mb_removed,aes(x=ind,y=tP/nSites,color=pop))+
  geom_violin()+
  theme_classic()+
  ylab(expression(pi))+
  xlab("Sample")+ 
  guides(color = guide_legend(title = "Origin"))+
  theme(axis.text=element_text(size=10), axis.title=element_text(size=18))

ND_plot<-ND_plot + stat_summary(fun.data="mean_sdl", mult=1, 
                 geom="pointrange", width=0.2 )

# Use custom color palettes
ND_plot+scale_color_manual(breaks = c("Sierra", "Liberia", "Ivory","Ghana","Tad","Nigeria","Cameroon","Gabon"),values=c("lightblue", "orange", "purple","red","forestgreen",'green','grey','black'))

ggsave("individual_PI_plot.pdf",width = 12, height = 3)
```
