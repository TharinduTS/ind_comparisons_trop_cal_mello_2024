# ind_comparisons_trop_cal_mello_2024

I am going to cal FST, Admixture and PI individual wise - removing 1-20mb in chr7 for trop individuals using bam files

# ******************
PI and Admix was calculated inside the old directory new_project_Apr_2023
and FST was in ind_comparisons_trop_cal_mello_2024 
# *****************
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
Then I calculated admixture with following script

prep_and_cal_admix.sh

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
#SBATCH --array=2-5

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

# Load modules
module load nixpkgs/16.09  gcc/7.3.0
module load StdEnv/2023
module load angsd
module load gsl/2.5
module load htslib

angsd -GL 1 -out genolike -nThreads 10 -doGlf 2 -doMajorMinor 1 -doMaf 0 -doCounts 1 -setMinDepthInd 5 -setMaxDepth 100 -minMapQ 20 -bam bam.filelist -rf ../ref_file_to_remove_7_20.txt
NGSadmix -likes genolike.beagle.gz -K ${SLURM_ARRAY_TASK_ID} -P 4 -o myoutfiles${SLURM_ARRAY_TASK_ID} -minMaf 0.05
```
Then I did the same in a seperate directory for first 20mb of chr 7

created a new directory in the directory with bamfiles

```bash
mkdir admix_for_first_20mb
```
and copied bam.filelist there
and ran admix changing the region

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
#SBATCH --array=2-5

#SBATCH --mail-user=premacht@mcmaster.ca
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=REQUEUE
#SBATCH --mail-type=ALL

# Load modules
module load nixpkgs/16.09  gcc/7.3.0
module load StdEnv/2023
module load angsd
module load gsl/2.5
module load htslib

angsd -GL 1 -out genolike -nThreads 10 -doGlf 2 -doMajorMinor 1 -doMaf 0 -doCounts 1 -setMinDepthInd 5 -setMaxDepth 100 -minMapQ 20 -bam bam.filelist -r Chr7:1-20000000
NGSadmix -likes genolike.beagle.gz -K ${SLURM_ARRAY_TASK_ID} -P 4 -o myoutfiles${SLURM_ARRAY_TASK_ID} -minMaf 0.05
```
Then I
1)Downloaded all the data
2)put them in a directory named outfiles inside the directory with my R script to plot
3)Created pop_info.txt file with the same order as bam.filelist

pop_info.txt
```txt
F_Ghana_WZ_BJE4687_combined__sorted.bam	Ghana
F_IvoryCoast_xen228_combined__sorted.bam	Ivory_coast
F_Nigeria_EUA0331_combined__sorted.bam	Nigeria
F_Nigeria_EUA0333_combined__sorted.bam	Nigeria
F_SierraLeone_AMNH17272_combined__sorted.bam	Sierra_Leone
F_SierraLeone_AMNH17274_combined__sorted.bam	Sierra_Leone
all_ROM19161_sorted.bam	Liberia
XT10_WZ_no_adapt._sorted.bam	Lab_tads
XT11_WW_trim_no_adapt_scafconcat_sorted.bam	Lab_tads
M_Ghana_WY_BJE4362_combined__sorted.bam	Ghana
M_Ghana_ZY_BJE4360_combined__sorted.bam	Ghana
M_Nigeria_EUA0334_combined__sorted.bam	Nigeria
M_Nigeria_EUA0335_combined__sorted.bam	Nigeria
M_SierraLeone_AMNH17271_combined__sorted.bam	Sierra_Leone
M_SierraLeone_AMNH17273_combined__sorted.bam	Sierra_Leone
XT1_ZY_no_adapt._sorted.bam	Lab_tads
XT7_WY_no_adapt__sorted.bam	Lab_tads
```
4) and copied this R script in to the same directory and ran it
```R
#!/usr/bin/Rscript

# Usage: plotADMIXTURE.r -p <prefix> -i <info file, 2-column file with ind name and population/species name> 
#                        -k <max K value> -l <comma-separated list of populations/species in the order to be plotted>
# This R script makes barplots for K=2 and all other K values until max K (specified with -k). It labels the individuals 
# and splits them into populations or species according to the individual and population/species names in the 2-column file specified with -i.
# The order of populations/species follows the list of populations/species given with -l.
# Usage example: plotADMIXTURE.r -p fileXY -i file.ind.pop.txt -k 4 -pop pop1,pop2,pop3
# In this example, the script would use the files fileXY.2.Q, fileXY.3.Q, fileXY.4.Q to make barplots for the three populations.
# file.ind.pop.txt should contain one line for each individual in the same order as in the admixture files e.g.
# ind1 pop1
# ind2 pop1
# ind3 pop2
# ind4 pop3

# Author: Joana Meier, September 2019

# I am customizing this script. Making it easy to use with R studio-TS
# ***keeping this script in current directory for clarity. But setting the working directory to the directory with all output files***


#set working directory to all_outputs inside the current path
setwd(paste(dirname(rstudioapi::getSourceEditorContext()$path),"/outfiles",sep=""))

#set default values here so it can run on R studio without parsing options
#These will come into action if you do not define these options like you do in bash

#change prefix here
p_input<-"myoutfiles"
#change sample list file here
# ******************* nKEEP THE SAME ORDE AS IN bam.filelist ************
i_input<-paste(dirname(rstudioapi::getSourceEditorContext()$path),"/pop_info.txt",sep="")
# change maximum k value to plot here
k_input<-5
# change minimum k value to plot here
m_input<-2

#************** Change your pop order by changing the list below*******************#
#add a list of populations seperated by commas here. This should be exactly similiar to the populations in your sample list file. plots will be created according to this population order
# you will have to change this every time you edit sample list
l_input<-"Sierra_Leone,Liberia,Ivory_coast,Ghana,Lab_tads,Nigeria"


#add the location and file name for the plots here. I am setting this to the directory I am creating in the next line
o_input<-paste(dirname(rstudioapi::getSourceEditorContext()$path),"/plot_outs/plot",sep="")

#create a directory for plot output if it doesn't already exist in the directory with the script
dir.create(file.path(dirname(rstudioapi::getSourceEditorContext()$path), "plot_outs"))

# Read in the arguments
library("optparse")
option_list = list(
  make_option(c("-p", "--prefix"), type="character", default=p_input, 
              help="prefix name (with path if not in the current directory)", metavar="character"),
  make_option(c("-i", "--infofile"), type="character", default=i_input, 
              help="info text file containing for each individual the population/species information", metavar="character"),
  make_option(c("-k", "--maxK"), type="integer", default=k_input, 
              help="maximum K value", metavar="integer"),
  make_option(c("-m", "--minK"), type="integer", default=m_input, 
              help="minimum K value", metavar="integer"),
  make_option(c("-l", "--populations"), type="character", default=l_input, 
              help="comma-separated list of populations/species in the order to be plotted", metavar="character"),
  make_option(c("-o", "--outPrefix"), type="character", default=o_input, 
              help="output prefix (default: name provided with prefix)", metavar="character")
) 
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# Check that all required arguments are provided
if (is.null(opt$prefix)){
  print_help(opt_parser)
  stop("Please provide the prefix", call.=FALSE)
}else if (is.null(opt$infofile)){
  print_help(opt_parser)
  stop("Please provide the info file", call.=FALSE)
}else if (is.null(opt$maxK)){
  print_help(opt_parser)
  stop("Please provide the maximum K value to plot", call.=FALSE)
}else if (is.null(opt$populations)){
  print_help(opt_parser)
  stop("Please provide a comma-separated list of populations/species", call.=FALSE)
}

# If no output prefix is given, use the input prefix
if(opt$outPrefix=="default") opt$outPrefix=opt$prefix

# Assign the first argument to prefix
prefix=opt$prefix

# Get individual names in the correct order
labels<-read.table(opt$infofile)

# Name the columns
names(labels)<-c("ind","pop")

# Add a column with population indices to order the barplots
# Use the order of populations provided as the fourth argument (list separated by commas)
labels$n<-factor(labels$pop,levels=unlist(strsplit(opt$populations,",")))
levels(labels$n)<-c(1:length(levels(labels$n)))
labels$n<-as.integer(as.character(labels$n))

# read in the different admixture output files
minK=opt$minK
maxK=opt$maxK
tbl<-lapply(minK:maxK, function(x) read.table(paste0(prefix,x,".qopt")))

# Prepare spaces to separate the populations/species
rep<-as.vector(table(labels$n))
spaces<-0
for(i in 1:length(rep)){spaces=c(spaces,rep(0,rep[i]-1),0.15)}
spaces<-spaces[-length(spaces)]

#change the space between individuals
spaces<-replace(spaces, spaces==0, 0.02)

#change the space between populations- only change pop_space

pop_space<-0.10
spaces<-replace(spaces, spaces==0.15, pop_space)


# Plot the cluster assignments as a single bar for each individual for each K as a separate row
tiff(file=paste0(opt$outPrefix,".tiff"),width = 800, height = 500,res=200)

#*************** Change font size by changing cex.axis here****************#
# change 3rd value in "oma" if you wanna add labels to gain enough space. Use mai to adjust space between different k valued plots
par(mfrow=c(maxK-1,1),mar=c(0,1,0.1,1),oma=c(2,4,4,1),mgp=c(0,0.2,0),mai=c(0.05,0,0,0),xaxs="i",cex.lab=0.5,cex.axis=0.55)

#Create a list for plots
plot_list<-list()
#assign colors to populations at once

col1<-"red"
col2<-"green"
col3<-"orange"
col4<-"skyblue"
col5<-"darkblue"
col6<-"yellow"
col7<-"pink"
col8<-"purple"
col9<-"grey"
col10<-"black"
col11<-"forestgreen"
col12<-"brown"

#create color palettes for each k value
col_palette_k2<-c(col1,col2)
col_palette_k3<-c(col3,col2,col1)
col_palette_k4<-c(col4,col1,col2,col3)
col_palette_k5<-c(col2,col1,col4,col5,col3)


col_palette_k6<-c(col1,col2,col6,col5,col4,col3)
col_palette_k7<-c(col1,col7,col5,col6,col4,col3,col2)
col_palette_k8<-c(col5,col8,col1,col7,col6,col3,col4,col2)
col_palette_k9<-c(col8,col6,col7,col2,col5,col4,col3,col1,col9)
col_palette_k10<-c(col8,col7,col4,col10,col2,col9,col6,col1,col3,col5)
col_palette_k11<-c(col6,col10,col2,col4,col1,col7,col9,col8,col3,col11,col5)
col_palette_k12<-c(col1,col2,col3,col4,col5,col6,col7,col8,col9,col10,col11,col12)

# make a variable with full sample names 
full_sample_names<-labels$ind[order(labels$n)]

# get simpler names  for those complex sample names
library(plyr)
#first renaming the sample names which does not follow the specific format manually**DO NOT USE"_" HERE AS SPLIT WILL USE THIS IN NEXT LINE
simple_sample_list_changing<-mapvalues(full_sample_names, from = c("all_ROM19161_sorted.bam","XT10_WZ_no_adapt._sorted.bam","XT11_WW_trim_no_adapt_scafconcat_sorted.bam","XT1_ZY_no_adapt._sorted.bam","XT7_WY_no_adapt__sorted.bam"), 
                                       to = c("o_Liberia1","F_Lab1","F_Lab2","M_Lab3","M_Lab4"))

#convert facter list into chrs to rename
s_list_chr<-as.character(simple_sample_list_changing)


#then remove the parts after the first"_" from other samples
shortened_sample_list<-sapply(strsplit(s_list_chr,split = "_"),`[`, 2)

#I am changing samples names here to keep it consistant with other analysis
# renaming all samples
#****** THIS CANCELS OUT PREVIOUS STEPS**********
# renaming all samples
full_sample_names[full_sample_names == "F_SierraLeone_AMNH17272_combined__sorted.bam"] <- "SL_F1"
full_sample_names[full_sample_names == "F_SierraLeone_AMNH17274_combined__sorted.bam"] <- "SL_F2"
full_sample_names[full_sample_names == "M_SierraLeone_AMNH17271_combined__sorted.bam"] <- "SL_M1"
full_sample_names[full_sample_names == "M_SierraLeone_AMNH17273_combined__sorted.bam"] <- "SL_M2"
full_sample_names[full_sample_names == "all_ROM19161_sorted.bam"] <- "LB_F1"
full_sample_names[full_sample_names == "F_IvoryCoast_xen228_combined__sorted.bam"] <- "IC_F1"
full_sample_names[full_sample_names == "F_Ghana_WZ_BJE4687_combined__sorted.bam"] <- "GH_F1"
full_sample_names[full_sample_names == "M_Ghana_WY_BJE4362_combined__sorted.bam"] <- "GH_M1"
full_sample_names[full_sample_names == "M_Ghana_ZY_BJE4360_combined__sorted.bam"] <- "GH_M2"
full_sample_names[full_sample_names == "XT10_WZ_no_adapt._sorted.bam"] <- "LT_1"
full_sample_names[full_sample_names == "XT11_WW_trim_no_adapt_scafconcat_sorted.bam"] <- "LT_2"
full_sample_names[full_sample_names == "XT1_ZY_no_adapt._sorted.bam"] <- "LT_3"
full_sample_names[full_sample_names == "XT7_WY_no_adapt__sorted.bam"] <- "LT_4"
full_sample_names[full_sample_names == "F_Nigeria_EUA0331_combined__sorted.bam"] <- "NG_F1"
full_sample_names[full_sample_names == "F_Nigeria_EUA0333_combined__sorted.bam"] <- "NG_F2"
full_sample_names[full_sample_names == "M_Nigeria_EUA0334_combined__sorted.bam"] <- "NG_M1"
full_sample_names[full_sample_names == "M_Nigeria_EUA0335_combined__sorted.bam"] <- "NG_M2"

shortened_sample_list<-full_sample_names
#**********************************************

#re oreder samples

#shortened_sample_list<-c("Liberia1","SierraLeone","SierraLeone","SierraLeone","SierraLeone","IvoryCoast","Ghana","Ghana","Ghana","Lab1","Lab2","Lab3","Lab4","Nigeria","Nigeria","Nigeria","Nigeria" )


# following plots are written in a way you can just copy paste only changing k_val for the different number of 'k's
#paste whats inside * marks for different k values and then change k_val
#**********
# Plot k=2
# change only here
k_val<-2

points<-c(0.5,1.5,2.5,3.5,4.65,5.75,6.85,7.85,8.85,9.95,10.95,11.95,12.95,14.25,15.25,16.25,17.25)

#uncomment and change values of this and line after bp to manually change label placing
label_points<-c(2,4.3,5.95,8.00,11.80,15.72)

#this changes K value label position
p2<-c(0.60)

bp<-barplot(t(as.matrix(tbl[[k_val-1]][order(labels$n),])), col=get(paste("col_palette_k",k_val,sep="")),xaxt="n", border=NA,ylab=paste0("K=",k_val,sep=""),yaxt="n",space=spaces)
#using following line to use manual measures
axis(3,at=label_points,las=1,cex=6,tick = F,
     labels=unlist(strsplit(opt$populations,",")))
axis(2, at =p2,labels=c("K=2"),las=2,tick=F,cex=2)
  




#***********

#commwnt out this section to automate coloring. Keep commented to manually assign colours for each pop with the next section( you can do it once by selecting lines and ctrl+shift+c)
# Plot higher K values
#  if(maxK>minK)lapply(2:(maxK-1), function(x) barplot(t(as.matrix(tbl[[x]][order(labels$n),])), col=rainbow(n=x+1),xaxt="n", border=NA,ylab=paste0("K=",x+1),yaxt="n",space=spaces))
#  axis(1,at=c(which(spaces==0.15),bp[length(bp)])-diff(c(1,which(spaces==pop_space),bp[length(bp)]))/2,
#      labels=unlist(strsplit(opt$populations,",")))
# dev.off()

#**********
# Plot k
# change only here
k_val<-3


#this changes K value label position
p3<-c(0.5)

bp<-barplot(t(as.matrix(tbl[[k_val-1]][order(labels$n),])), col=get(paste("col_palette_k",k_val,sep="")),xaxt="n", border=NA,ylab=paste0("K=",k_val,sep=""),yaxt="n",space=spaces)
axis(2, at =p3,labels=c("K=3"),las=2,tick=F,cex=2)




#***********

#**********
# Plot k
# change only here
k_val<-4


#this changes K value label position
p4<-c(0.50)

bp<-barplot(t(as.matrix(tbl[[k_val-1]][order(labels$n),])), col=get(paste("col_palette_k",k_val,sep="")),xaxt="n", border=NA,ylab=paste0("K=",k_val,sep=""),yaxt="n",space=spaces)
axis(2, at =p4,labels=c("K=4"),las=2,tick=F,cex=2)



#***********

#**********
# Plot k
# change only here
k_val<-5

#this changes K value label position
p4<-c(0.45)

bp<-barplot(t(as.matrix(tbl[[k_val-1]][order(labels$n),])), col=get(paste("col_palette_k",k_val,sep="")),xaxt="n", border=NA,ylab=paste0("K=",k_val,sep=""),yaxt="n",space=spaces)
#axis(1,at=c(which(spaces==pop_space),bp[length(bp)])-diff(c(1,which(spaces==pop_space),bp[length(bp)]))/2,
axis(1, at =points,labels=shortened_sample_list,las=2,tick=F,cex=2)
axis(2, at =p4,labels=c("K=5"),las=2,tick=F,cex=2)


# #***********
# #**********
# # Plot k
# # change only here
# k_val<-12
# 
# 
# 
# 
# bp<-barplot(t(as.matrix(tbl[[k_val-1]][order(labels$n),])), col=get(paste("col_palette_k",k_val,sep="")),xaxt="n", border=NA,ylab=paste0("K=",k_val,sep=""),yaxt="n",space=spaces)
# axis(1,at=c(which(spaces==0.15),bp[length(bp)])-diff(c(1,which(spaces==0.15),bp[length(bp)]))/2,
#      labels=unlist(strsplit(opt$populations,",")))

dev.off()
```

Then I did plot first 20mb data with the same way

I did repeat this two more times with different seed values and outfile names to check whether results are reproducable

# PI/Nucleotide diversity for individuals

#** Here I removed Chr7 first 20 mb in the R script********

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
twenty_mb_removed[twenty_mb_removed == "XT10"] <- "LT_1"
twenty_mb_removed[twenty_mb_removed == "XT11"] <- "LT_2"
twenty_mb_removed[twenty_mb_removed == "XT1"] <- "LT_3"
twenty_mb_removed[twenty_mb_removed == "XT7"] <- "LT_4"
twenty_mb_removed[twenty_mb_removed == "EUA0331"] <- "NG_F1"
twenty_mb_removed[twenty_mb_removed == "EUA0333"] <- "NG_F2"
twenty_mb_removed[twenty_mb_removed == "EUA0334"] <- "NG_M1"
twenty_mb_removed[twenty_mb_removed == "EUA0335"] <- "NG_M2"

twenty_mb_removed[twenty_mb_removed == "cal"] <- "Xcal"
twenty_mb_removed[twenty_mb_removed == "mello"] <- "Xmel"

#change pop names

twenty_mb_removed[twenty_mb_removed == "Cal"] <- "Cameroon"
twenty_mb_removed[twenty_mb_removed == "Mello"] <- "Gabon"

#Get rid of label tad

twenty_mb_removed[twenty_mb_removed == "Tad"] <- "Ghana"

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
               'LT_1',
               'LT_2',
               'LT_3',
               'LT_4',
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
  theme(axis.text=element_text(size=10), axis.title=element_text(size=16))

ND_plot<-ND_plot + stat_summary(fun.data="mean_sdl", mult=1, 
                 geom="pointrange", width=0.2 )

# Use custom color palettes
ND_plot+scale_color_manual(breaks = c("Sierra", "Liberia", "Ivory","Ghana","Tad","Nigeria","Cameroon","Gabon"),values=c("lightblue", "orange", "purple","red","red",'green','grey','black'))

ggsave("individual_PI_plot.pdf",width = 12, height = 3)
  
```
