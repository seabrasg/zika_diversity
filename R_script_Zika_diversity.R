#####################################################################################
#  Authors: Kristof Theys and Sofia Seabra (sgseabra@ihmt.unl.pt)

# R script for the analyses of the manuscript:

# Genome-wide diversity of Zika virus: exploring spatio-temporal dynamics to guide a new nomenclature proposal
#
# (Seabra SG, Theys K, Libin PJK, Zhukova A, Potter BI, Nebenzahl-Guimaraes H, Gorbalenya AE, Sidorov IA, 
# Pimentel V, Pingarilho M, de Vasconcelos, ATR, Dellicour S, Khouri R, Gascuel O, Vandamme AM, Baele G, 
# Cuypers L, Abecasis AB)

# Input files:
# AA_mutation_table_from_VIRULIGN.csv
# NT_mutation_table_from_VIRULIGN.csv
# Metadata_505_with_results.tsv 
# Zika_Aligned_505.fasta
# Zika_Aligned_437.fasta

# datamonkey-table_505.csv
# datamonkey-table_437.csv
# datamonkey-table_83.csv

# INDEX:
# 1. FILTER MUTATION FILES
# 2. DIVERSITY ANALYSES
# 3. GENETIC DISTANCES
# 4. AMOVA - analysis of molecular variance
# 5. RhierBAPS - bayesian clustering
# 6. PLOTS RESULTS OF FUBAR method for detection of selection

# Source file with functions used in this script:
source("./R_functions_Zika.R")

# Packages needed
library(ggplot2)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(dplyr)
library(tidyr)
library(plyr)
library(stringr)
library(phytools)
library(phangorn)
library(data.table)
library(ape)
library(vegan)
library(apex)
library(adegenet)
library(pegas)
library(mmod)
library(poppr)
library(ggtree)
library(rhierbaps)

setwd("path_to_working_directory")

proteinorder<-c("C","PR","M","E","NS1","NS2A","NS2B","NS3","NS4A","pep2K","NS4B","NS5")

######################################################################################
######################################################################################
# 1. FILTER MUTATION FILES
######################################################################################
######################################################################################

## AA MUTATION file - AMINOACIDS
inputaa <-read.csv("./AA_mutation_table_from_VIRULIGN.csv",sep=',',na.strings=c("-","","*","ACDEFGHIKLMNPQRSTVWY*",'X','NA'),
                   header=TRUE,  colClasses = c("character"))
# get index of column for the position in genome (first column is seqid, so the index is +1 in relation to original sequence)
begin_genome_aa<-match('C_1',names(inputaa))
end_genome_aa<-match('NS5_903',names(inputaa))

# AA check and filtering - many with multiple values (e.g. ABCD...) due to ambiguities and Ns
######################################
sort(table(unlist(as.data.frame(inputaa[,begin_genome_aa:end_genome_aa]))))
barplot(tail(sort(prop.table(table(unlist(as.data.frame(inputaa[,begin_genome_aa:end_genome_aa]))))),21))

# check multiple values 
testlevels<-apply(inputaa[,begin_genome_aa:end_genome_aa],2,function(x) levels(factor(x)))  # each positions, levels 
multiplevalues<-NULL;i<-1;for (k in 1:length(testlevels)){if(length(testlevels[names(testlevels)[k]][[1]])>1){multiplevalues<-c(multiplevalues,(testlevels[k][[1]]));i<-i+1}}
length(unique(multiplevalues))
multiplevalues[nchar(multiplevalues)!=1] # check single letters
sort(unique(multiplevalues[nchar(multiplevalues)==1]))

# Convert all aminoacids with "multiple values" (more than 1 aminoacid per position per sample) in NA (missing values)
for(z in 1:nrow(inputaa)) {
  inputaa[z,which(nchar(inputaa[z,])>1)[names(which(nchar(inputaa[z,])>1)) %in% names(inputaa)[begin_genome_aa:end_genome_aa]]] <-NA
  inputaa[z,grep('\\*',inputaa[z,])] <-NA
  #w_qc[z,which(nchar(w_qc[z,])>4)]<-'bad_codon'
} 

# Check if there are no multiple values after the filtering
sort(table(unlist(as.data.frame(inputaa[,begin_genome_aa:end_genome_aa])),useNA="always"))
barplot(tail(sort(prop.table(table(unlist(as.data.frame(inputaa[,begin_genome_aa:end_genome_aa])),useNA="always"))),21))
barplot(prop.table(table(unlist(as.data.frame(inputaa[,begin_genome_aa:end_genome_aa])),useNA="always")))

# Write filtered inputaa
write.table(inputaa,"./inputaa_total.csv",sep=",", quote = FALSE, row.names = FALSE)


######################################################################################
######################################################################################
## NT MUTATION file - NUCLEOTIDES
inputnt <-read.csv("./NT_mutation_table_from_VIRULIGN.csv",sep=',',na.strings=c("-","","N","K","Y","W","R","M","S","B","NA"),
                   header=TRUE,  colClasses = c("character"))
dim(inputnt)

# get index of column for the position in genome (first column is seqid, so the index is +1 in relation to original sequence)
begin_genome_nt<-match('C_1_1',names(inputnt))
end_genome_nt<-match('NS5_903_3',names(inputnt))

# NT check and filtering
######################################3
unique(unlist(as.data.frame(inputnt[,begin_genome_nt:end_genome_nt])))
sort(table(unlist(as.data.frame(inputnt[,begin_genome_nt:end_genome_nt])),useNA="always"))
barplot(sort(prop.table(table(unlist(as.data.frame(inputnt[,begin_genome_nt:end_genome_nt])),useNA="always"))))

# Write filtered inputnt
write.table(inputnt,"./inputnt_total.csv",sep=",", quote = FALSE, row.names = FALSE)

######################################################################################
######################################################################################

#################################
# CHECK INPUT FILES 
############################
inputaa <-read.csv("./inputaa_total.csv",sep=',', header=TRUE,  colClasses = c("character"))
dim(inputaa)
begin_genome_aa<-match('C_1',names(inputaa))
end_genome_aa<-match('NS5_903',names(inputaa))
table(unlist(as.data.frame(inputaa[,begin_genome_aa:end_genome_aa])),useNA = "always")


inputnt <-read.csv("./inputnt_total.csv",sep=',', header=TRUE,  colClasses = c("character"))
dim(inputnt)
begin_genome_nt<-match('C_1_1',names(inputnt))
end_genome_nt<-match('NS5_903_3',names(inputnt))
table(unlist(as.data.frame(inputnt[,begin_genome_nt:end_genome_nt])),useNA = "always")

######################################################################################
# FIND INDEX OF START AND END for each protein in the input file (start at 2 in the file)
# column index 1 is name; column index 2 is the first position in genome

# inputaa 
############################################################################
protein_length_aa<-table(gsub("_.*","",names(inputaa[begin_genome_aa:end_genome_aa])))
protein_index_aa<- data.frame(protein=character(), index_start=integer(), index_end=integer(), stringsAsFactors=FALSE)
for (prot in 1:(length(names(protein_length_aa))))  {
  protein_index_aa[prot,1] <-  names(protein_length_aa)[prot]
  protein_index_aa[prot,2] <- min(grep(paste(names(protein_length_aa)[prot],'_',sep=''),names(inputaa))) 
  protein_index_aa[prot,3] <- max(grep(paste(names(protein_length_aa)[prot],'_',sep=''),names(inputaa)))
}
protein_index_aa<-protein_index_aa[order(protein_index_aa[,2]),]  # order by index
protein_index_aa

# Start and end postions of each protein in genome
protein_position_aa<-protein_index_aa
protein_position_aa$start_pos<-protein_index_aa$index_start-1
protein_position_aa$end_pos<-protein_index_aa$index_end-1
protein_position_aa$length<-protein_position_aa$end_pos-protein_position_aa$start_pos+1
sum(protein_position_aa$length)
protein_position_aa
write.table(protein_position_aa,"./Protein_start_end_positions_lenght_aa.tsv",sep="\t", quote = FALSE, row.names = FALSE)


# inputnt
############################################################################
protein_length_nt<-table(gsub("_.*","",names(inputnt[begin_genome_nt:end_genome_nt])))
protein_index_nt<- data.frame(protein=character(), index_start=integer(), index_end=integer(), stringsAsFactors=FALSE)
for (prot in 1:(length(names(protein_length_nt)))) {
  protein_index_nt[prot,1] <-  names(protein_length_nt)[prot]
  protein_index_nt[prot,2] <- min(grep(paste(names(protein_length_nt)[prot],'_',sep=''),names(inputnt))) #-1 since C should start at pos 1, not pos 2 since pos is in front
  protein_index_nt[prot,3] <- max(grep(paste(names(protein_length_nt)[prot],'_',sep=''),names(inputnt)))
}
protein_index_nt<-protein_index_nt[order(protein_index_nt[,2]),]  # order 

# Start and end postions of each protein in genome
protein_position_nt<-protein_index_nt
protein_position_nt$start_pos<-protein_position_nt$index_start-1
protein_position_nt$end_pos<-protein_position_nt$index_end-1
protein_position_nt$length<-protein_position_nt$end_pos-protein_position_nt$start_pos+1
sum(protein_position_nt$length)
protein_position_nt
write.table(protein_position_nt,"./Protein_start_end_positions_lenght_nt.tsv",sep="\t", quote = FALSE, row.names = FALSE)


############################################################################
# MAKE A FULL GENOME MAPPING - start at 1 instead of 2 (which was due to index 2 in file)
## aa
xx<-NULL;
for (x in 1:nrow(protein_index_aa)) { 
  xx<-c(xx,(protein_index_aa$index_start[x]:protein_index_aa$index_end[x])-protein_index_aa$index_start[x]+1)
}

length_region<-length(min(begin_genome_aa:end_genome_aa):max(begin_genome_aa:end_genome_aa))
genome_protein_positions_aa<-cbind.data.frame(gsub("_.*","",names(inputaa)[begin_genome_aa:end_genome_aa]), 1:length_region, xx)
colnames(genome_protein_positions_aa)<-c('Protein','Pos_genome_aa','Pos_protein_aa')
write.table(genome_protein_positions_aa,"./Mapping_genome_protein_positions_aa.tsv",sep="\t", quote = FALSE, row.names = FALSE)

#same for nt 
xx_nt<-NULL;
for (x in 1:nrow(protein_index_nt)) {
  xx_nt<-c(xx_nt,(protein_index_nt$index_start[x]:protein_index_nt$index_end[x])-protein_index_nt$index_start[x]+1)
}

length_region_nt<-length(min(begin_genome_nt:end_genome_nt):max(begin_genome_nt:end_genome_nt))
genome_protein_positions_nt<-cbind.data.frame(gsub("_.*","",names(inputnt)[begin_genome_nt:end_genome_nt]), 1:length_region_nt, xx_nt)
colnames(genome_protein_positions_nt)<-c('Protein','Pos_genome_nt','Pos_protein_nt')
write.table(genome_protein_positions_nt,"./Mapping_genome_protein_positions_nt.tsv",sep="\t", quote = FALSE, row.names = FALSE)

#######################
# REFERENCE SEQUENCES
inputref_aa <-inputaa[inputaa$seqid=='KJ776791',]
dim(inputref_aa)
write.table(inputref_aa,"./inputref_aa.csv",sep=",", quote = FALSE, row.names = FALSE)

inputNATALref_aa <-inputaa[inputaa$seqid=='KU527068',]
dim(inputNATALref_aa)
write.table(inputNATALref_aa,"./inputNATALref_aa.csv",sep=",", quote = FALSE, row.names = FALSE)


inputref_nt <-inputnt[inputnt$seqid=='KJ776791',]
dim(inputref_nt)
write.table(inputref_nt,"./inputref_nt.csv",sep=",", quote = FALSE, row.names = FALSE)

inputNATALref_nt <-inputnt[inputnt$seqid=='KU527068',]
dim(inputNATALref_nt)
write.table(inputNATALref_nt,"./inputNATALref_nt.csv",sep=",", quote = FALSE, row.names = FALSE)



######################################################################################
######################################################################################
## 2. DIVERSITY ANALYSES
######################################################################################

# INPUT FILES (matrix samples vs positions)
###############################
inputaa <-read.csv("./inputaa_total.csv",sep=',', header=TRUE,  colClasses = c("character"))
dim(inputaa)
# get index of column for the position in genome (first column is seqid, so the index is +1 in relation to original sequence)
begin_genome_aa<-match('C_1',names(inputaa))
end_genome_aa<-match('NS5_903',names(inputaa))

table(unlist(as.data.frame(inputaa[,begin_genome_aa:end_genome_aa])),useNA = "always")

##############################
inputnt <-read.csv("./inputnt_total.csv",sep=',', header=TRUE,  colClasses = c("character"))
dim(inputnt)
# get index of column for the position in genome (first column is seqid, so the index is +1 in relation to original sequence)
begin_genome_nt<-match('C_1_1',names(inputnt))
end_genome_nt<-match('NS5_903_3',names(inputnt))

table(unlist(as.data.frame(inputnt[,begin_genome_nt:end_genome_nt])),useNA = "always")

# POSITIONS FILES
#############################
genome_protein_positions_aa<-read.csv("./Mapping_genome_protein_positions_aa.tsv",sep="\t",header=TRUE,  colClasses = c("character"))

protein_position_aa<-read.csv("./Protein_start_end_positions_lenght_aa.tsv",sep='\t', header=TRUE,  colClasses = c("character"))

#############################
genome_protein_positions_nt<-read.csv("./Mapping_genome_protein_positions_nt.tsv",sep="\t",header=TRUE,  colClasses = c("character"))

protein_position_nt<-read.csv("./Protein_start_end_positions_lenght_nt.tsv",sep='\t', header=TRUE,  colClasses = c("character"))


# REFERENCES
#############################
inputref_aa <-read.csv("./inputref_aa.csv",sep=',', header=TRUE,  colClasses = c("character"))
dim(inputref_aa)

inputNATALref_aa <-read.csv("./inputNATALref_aa.csv",sep=',', header=TRUE,  colClasses = c("character"))
dim(inputNATALref_aa)

#############################
inputref_nt <-read.csv("./inputref_nt.csv",sep=',', header=TRUE,  colClasses = c("character"))
dim(inputref_nt)

inputNATALref_nt <-read.csv("./inputNATALref_nt.csv",sep=',', header=TRUE,  colClasses = c("character"))
dim(inputNATALref_aa)

##########################################################################
# Get frequency matrix of AA (number of samples harboring each variant in each position -
##########################################################################
source("./R_functions_Zika.R")
freq_matrixaa<-makeMatrixAA(inputaa[,begin_genome_aa:end_genome_aa])
dim(freq_matrixaa) #rows: position (3423); columns: aminoacids (23 - 21+ 2 columns))
colnames(freq_matrixaa)

# check particular position
freq_matrixaa[which(freq_matrixaa[,1]=="E_66"),]

write.table(freq_matrixaa, "./freq_matrixaa.tsv", sep="\t",quote=F, col.names=TRUE, row.names=F)

# Open file
freq_matrixaa<-read.csv("./freq_matrixaa.tsv",sep='\t', header=TRUE,  colClasses = c("character"))
dim(freq_matrixaa)

##########################################################################
# Consensus from the AA distribution frequency file 
##########################################################################
# list of genome position (name of protein _ position in protein)
WTpos_aa<-freq_matrixaa[,1]
# frequency values of the aminoacid with highest frequency for each position
WTvalue_aa<-as.numeric(as.character(apply(freq_matrixaa[,-c(1,2)],1,function(x) x[which(x==max(as.numeric(x)))])))
# aminoacid with highest frequency for each position
WTaa<-apply(freq_matrixaa[,-c(1,2)],1,function(x) names(x)[which(x==max(as.numeric(x)))])

print(length(WTpos_aa));print(length(WTvalue_aa));print(length(WTaa));
# create dataframe
consensus_aa = as.data.frame(cbind(WTpos_aa,WTvalue_aa,WTaa))

# rownames to colum (it is the position in genome)
consensus_aa<-tibble::rownames_to_column(consensus_aa,"Pos")
# convert tibble to dataframe
consensus_aa<-as.data.frame(consensus_aa)
dim(consensus_aa)

write.table(consensus_aa, "./consensus_aa.tsv", sep="\t",quote=F, col.names=TRUE, row.names=F)

consensus_aa<-read.csv("./consensus_aa.tsv",sep='\t', header=TRUE,  colClasses = c("character"))
dim(consensus_aa)

plot(as.character(consensus_aa$WTvalue_aa))
max(as.character(consensus_aa$WTvalue_aa))

####################################################################3
# Plot names of aminoacids with colours
####################################################################3

color1<-ifelse(as.numeric(as.character(consensus_aa$WTvalue_aa)) == 100 ,"lightgrey",
               ifelse(as.numeric(as.character(consensus_aa$WTvalue_aa)) > 99, "black",
                      ifelse(as.numeric(as.character(consensus_aa$WTvalue_aa)) > 95,"blue",
                             ifelse(as.numeric(as.character(consensus_aa$WTvalue_aa)) > 90,"green",
                                    ifelse(as.numeric(as.character(consensus_aa$WTvalue_aa)) > 75,"orange","red")))))

jpeg(file = "./consensusAA_color_mutations.jpg", 
     bg="white", antialias = "default", width = 20, height =30, units = "in", res = 300)
par(xpd=T, mar=par()$mar+c(0,0,0,0))
length<-25
ypos <- 1*length*floor(nrow(consensus_aa)/length)
plot(as.numeric(consensus_aa$WTvalue_aa[1:length]),type='n',ylim=c(0,ypos),xaxt='n',yaxt='n', 
     main="",ylab="",xlab="",bty='n',axes=FALSE)
round<-0
for (i in 1:nrow(consensus_aa) ){#print(i)
  print(i)
  text(i-length*round,ypos-length*round,paste(i,consensus_aa$WTaa[i],sep = ""),col=color1[i],cex=1.2)
  if (i == length+length*round)
  {round<-round+1}
}

legend(length/2-3,-20, inset=.05, c("0","]0-1]","]1-5]","]5-10]","]10-25]","]25-100]"), fill=c("lightgrey","black","blue","green","orange","red"),title=("Frequency of non-consensus amino acids (%)"),cex=0.8,horiz=TRUE)
dev.off()

#############################################
# Compare consensus with REFERENCE (get position, frequency and nucleotides that differ between consensus and ref)
#############################################

source("./R_functions_Zika.R")
consensus_vs_reference_aa<-function_consensus_vs_ref(consensus_aa, inputref_aa[,begin_genome_aa:end_genome_aa], genome_protein_positions_aa)

consensus_vs_referenceNATAL_aa<-function_consensus_vs_ref(consensus_aa, inputNATALref_aa[,begin_genome_aa:end_genome_aa], genome_protein_positions_aa)

# View positions in inputaa that have different aminoacids in consensus and reference - Check the number of sequences for each aminoacid at those positions
apply(consensus_vs_reference_aa,1,function(x)  table(inputaa[,(x[2])]))

for (d in 1:nrow(consensus_vs_reference_aa)) {
  print(as.character(consensus_vs_reference_aa$WTpos_aa)[d])
  print(table(inputaa[,as.character(consensus_vs_reference_aa$WTpos_aa)[d]]))  # verschil met refseq
}

# NATAL ref
apply(consensus_vs_referenceNATAL_aa,1,function(x)  table(inputaa[,(x[2])]))

for (d in 1:nrow(consensus_vs_referenceNATAL_aa)) {
  print(as.character(consensus_vs_referenceNATAL_aa$WTpos_aa)[d])
  print(table(inputaa[,as.character(consensus_vs_referenceNATAL_aa$WTpos_aa)[d]]))  # verschil met refseq
}

######################################################################3
#Get positions with consensus frequency less than 90%
positions_under_90 <- consensus_aa[as.numeric(as.character(consensus_aa[,3])) < 90,]

# View positions in inputaa 
apply(positions_under_90,1,function(x)  table(inputaa[,(x[2])]))

for (d in 1:nrow(positions_under_90)) {
  print(as.character(positions_under_90$WTpos_aa)[d])
  print(table(inputaa[,as.character(positions_under_90$WTpos_aa)[d]]))  # # variable positions  # check hoeveel aa hier
}

################################################################
# Proportion of positions with major aminoacid frequencies in the categories '<50','50<75','75<95','95<99','99<100'
################################################################
binaa=cut(as.numeric(as.character(consensus_aa[,3])), breaks=c(0,50,75,95,99,100),  labels=c('<50','50<75','75<95','95<99','99<100'))
table(binaa);round(prop.table(table(binaa)),3)*100

# Per protein
aa_cuts<-NULL
for(z in 1:12) {
  print(unique(genome_protein_positions_aa[,3])[z])
  temp_x<-as.numeric(as.character(consensus_aa[gsub('_.*','',consensus_aa[,2])==unique(genome_protein_positions_aa[,1])[z],3]))
  temp_cut<-as.numeric(round(prop.table(table(cut(temp_x,breaks=c(0,50,75,95,99,100),  labels=c('<50','50<75','75<95','95<99','99<100')))),3)*100)
  aa_cuts<-rbind(aa_cuts,c(as.character(unique(genome_protein_positions_aa[,1])[z]),temp_cut))
}
colnames(aa_cuts)<-c('Protein','<50','50<75','75<95','95<99','99<100')
write.table(aa_cuts, "./aa_cuts_freq_per_protein.tsv", sep="\t",quote=F, col.names=TRUE, row.names=F)

aa_cuts<-read.csv("./aa_cuts_freq_per_protein_total.tsv",sep='\t', header=TRUE,  colClasses = c("character"))
colnames(aa_cuts)<-c('Protein','<50','50<75','75<95','95<99','99<100')

barplot(t(aa_cuts[,2:6]),names.arg=c(aa_cuts[,1]),col=rainbow(5),legend=T,las=2)

jpeg(file = "./proportion_positions_freq_aa.jpg", 
     bg="white", antialias = "default", width = 6, height = 3, units = "in", res = 100)
barplot(t(aa_cuts[,2:6]),names.arg=c(aa_cuts[,1]),col=rainbow(5),legend=T,las=2)
dev.off()

#list of positions with major variant frequency <95
c95<-consensus_aa[as.numeric(as.character(consensus_aa[,3])) <95,]
c95
#___write.table(c95, paste0("./2_Mutations/output/variable_positions_aa_freq_lower95",dataset,".tsv"), sep="\t",quote=F, col.names=TRUE, row.names=F)

#ordered according to freq
consensus_aa[as.numeric(as.character(consensus_aa[,3])) <95,] [order(consensus_aa[as.numeric(as.character(consensus_aa[,3])) <95,][,2]),]


## STATISTICAL TEST OF  DIFFERENCE BETWEEN FREQUENCY OF CONSENSUS in each protein
# Kruskal Wallis test 
kruskal.test(as.numeric(consensus_aa[,3])~gsub('_.*','',consensus_aa[,2]))

# ANOVA
summary(aov(as.numeric(consensus_aa[,3])~gsub('_.*','',consensus_aa[,2])))

pairwise.wilcox.test(as.numeric(consensus_aa[,3]),gsub('_.*','',consensus_aa[,2]),
                     p.adjust.method = "BH")
boxplot(as.numeric(consensus_aa[,3])~gsub('_.*','',consensus_aa[,2]))

# ANOVA
summary(aov(as.numeric(consensus_aa[,3])~gsub('_.*','',consensus_aa[,2])))

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
##ENTROPY** AA level
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# Calculate entropy (R package entropy) using the function in source file /R_functions_Zika.R 

library(entropy)
proteinorder<-c("C","PR","M","E","NS1","NS2A","NS2B","NS3","NS4A","pep2K","NS4B","NS5")
source("./R_functions_Zika.R")
entropy_aa<-f_entropy(inputaa[begin_genome_aa:end_genome_aa])
colnames(entropy_aa)

sum(as.numeric(entropy_aa$Entropy))/length(entropy_aa$Entropy) #average protease entropy

write.table(entropy_aa, "./entropy_aa.tsv", sep="\t",quote=F, col.names=TRUE, row.names=F)

entropy_aa<-read.csv("./entropy_aa.tsv",sep='\t', header=TRUE,  colClasses = c("character"))
colnames(entropy_aa)


# Kruskal Wallis test - difference in entropy between proteins 
kruskal.test(entropy_aa$Entropy~entropy_aa$Protein)


pairwise.wilcox.test(entropy_aa$Entropy,entropy_aa$Protein,
                     p.adjust.method = "BH")
# ANOVA
anova_test<-aov(entropy_aa$Entropy~entropy_aa$Protein)
summary(anova_test)
# Check homogeneity of variance
LeveneTest(anova_test) 

# Check normality of residuals
anova_residuals<-residuals(anova_test)
shapiro.test(anova_residuals)

hist(entropy_aa$Entropy)
qqnorm(entropy_aa$Entropy)

# Check independence
DurbinWatsonTest(anova_test)

# Statistics for each protein
library(dplyr)
library(gmodels)
output_ci<-ci(entropy_aa$Entropy, na.rm = TRUE) # from package gmodels
output_ci

output_ci<-group_by(entropy_aa,Protein) %>%
  dplyr::summarise(
    mean_ci = ci(Entropy, na.rm = TRUE)[1],
    CI95_low = ci(Entropy, na.rm = TRUE)[2],
    CI95_high = ci(Entropy, na.rm = TRUE)[3],
    SE = ci(Entropy, na.rm = TRUE)[4])
output_ci


############################################################################################
##########################################################################
# NUCLEOTIDES
##########################################################################
# Get frequency matrix of NT
##########################################################################
source("./R_functions_Zika.R")
freq_matrixnt<-makeMatrixNT(inputnt[,begin_genome_nt:end_genome_nt])
dim(freq_matrixnt) #rows: position (3423); columns: nucleotides (23))
colnames(freq_matrixnt)
freq_matrixnt[1,2]

# check if sum of frequencies is 100%
sum_df<-NULL
for(i in 1:nrow(freq_matrixnt)){
  sum_df<-rbind(sum_df,sum(as.numeric(freq_matrixnt[i,3]),as.numeric(freq_matrixnt[i,4]),
                           as.numeric(freq_matrixnt[i,5]),as.numeric(freq_matrixnt[i,6])))
}
sum_df<-cbind(freq_matrixnt[,1],sum_df)
sum_df[which(as.numeric(sum_df[,2])<100),]

# check particular position
freq_matrixnt[which(freq_matrixnt[,1]=="PR_18_3"),]
freq_matrixnt[which(freq_matrixnt[,1]=="NS1_324_2"),]
freq_matrixnt[which(freq_matrixnt[,1]=="NS1_324_3"),]

write.table(freq_matrixnt, "./freq_matrixnt.tsv", sep="\t",quote=F, col.names=TRUE, row.names=F)

# Open file
freq_matrixnt<-read.csv("./freq_matrixnt.tsv",sep='\t', header=TRUE,  colClasses = c("character"))
dim(freq_matrixnt)


##########################################################################
# Consensus from the NT distribution frequency file 
##########################################################################
# list of genome position (name of gene _ position in gene)
WTpos_nt<-freq_matrixnt[,1]
# frequency values of the nucleotides with highest frequency for each position
WTvalue_nt<-as.numeric(as.character(apply(freq_matrixnt[,-c(1,2)],1,function(x) x[which(x==max(as.numeric(x)))])))
# nucleotide with highest frequency for each position
WTnt<-apply(freq_matrixnt[,-c(1,2)],1,function(x) names(x)[which(x==max(as.numeric(x)))])
print(length(WTpos_nt));print(length(WTvalue_nt));print(length(WTnt));
consensus_nt = as.data.frame(cbind(WTpos_nt,WTvalue_nt,WTnt))
# rownames to colum (it is the position in genome)
consensus_nt<-tibble::rownames_to_column(consensus_nt,"Pos")

write.table(consensus_nt, "./consensus_nt.tsv", sep="\t",quote=F, col.names=TRUE, row.names=F)

consensus_nt<-read.csv("./consensus_nt.tsv",sep='\t', header=TRUE,  colClasses = c("character"))
dim(consensus_nt)

plot(as.character(consensus_nt$WTvalue_nt))
max(as.character(consensus_nt$WTvalue_nt))

hist(sort((round(100-as.numeric(as.character(consensus_nt$WTvalue_nt)),4))),breaks=1000)
plot(sort((round(100-as.numeric(as.character(consensus_nt$WTvalue_nt)),4))),type='l')


########################3
# REFERENCE
# Compare consensus with reference (get position, frequency and nucleotides that differ between consensus and ref)

source("./R_functions_Zika.R")
consensus_vs_reference_nt<-function_consensus_vs_ref(consensus_nt, inputref_nt[,begin_genome_nt:end_genome_nt], genome_protein_positions_nt)

consensus_vs_referenceNATAL_nt<-function_consensus_vs_ref(consensus_nt, inputNATALref_nt[,begin_genome_nt:end_genome_nt], genome_protein_positions_nt)


# check particular position
inputref_nt[,which(names(inputref_nt)=="C_107_3")]
freq_matrixnt[which(freq_matrixnt[,1]=="C_107_3"),]
genome_protein_positions_nt[which(genome_protein_positions_nt[,1]=="PR" & genome_protein_positions_nt[,3]=="18"),]

# View positions in inputnt that have different nucleotides in consensus and reference - Check the number of sequences for each aminoacid at those positions
apply(consensus_vs_reference_nt,1,function(x)  table(inputnt[,(x[2])]))

for (d in 1:nrow(consensus_vs_reference_nt)) {
  print(as.character(consensus_vs_reference_nt$WTpos_nt)[d])
  print(table(inputnt[,as.character(consensus_vs_reference_nt$WTpos_nt)[d]]))  # verschil met refseq
}

######################################################################3
#Get positions with consensus frequency less than 90%
positions_under_90 <- consensus_nt[as.numeric(as.character(consensus_nt[,3])) < 90,]

apply(positions_under_90,1,function(x)  table(inputnt[,(x[2])]))

for (d in 1:nrow(positions_under_90)) {
  print(as.character(positions_under_90$WTpos_nt)[d])
  print(table(inputnt[,as.character(positions_under_90$WTpos_nt)[d]]))  # # variable positions  # check hoeveel aa hier
}

################################################################
# Proportion of positions with major nucleotide frequencies in the categories '<50','50<75','75<95','95<99','99<100'
################################################################
binnt<-cut(as.numeric(as.character(consensus_nt[,3])), breaks=c(0,50,75,95,99,100), labels=c('<50','50<75','75<95','95<99','99<100'))
table(binnt);round(prop.table(table(binnt)),2)*100

# 1 protein, verschillende thresholds 
nt_cuts<-NULL
for(z in 1:12){
  print(unique(genome_protein_positions_nt[,3])[z])
  temp_x<-as.numeric(as.character(consensus_nt[gsub('_.*','',consensus_nt[,2])==unique(genome_protein_positions_nt[,1])[z],3]))
  temp_cut<-as.numeric(round(prop.table(table(cut(temp_x,breaks=c(0,50,75,95,99,100),  labels=c('<50','50<75','75<95','95<99','99<100')))),3)*100)
  nt_cuts<-rbind(nt_cuts,c(as.character(unique(genome_protein_positions_nt[,1])[z]),temp_cut))
}
colnames(nt_cuts)<-c('Prot','<50','50<75','75<95','95<99','99<100')
write.table(nt_cuts, "./nt_cuts_freq_per_protein.tsv", sep="\t",quote=F, col.names=TRUE, row.names=F)

barplot(t(nt_cuts[,2:6]),names.arg=c(proteinorder),col=rainbow(5),legend=T,las=2,ylab='% of positions')

jpeg(file = "./proportion_positions_freq_nt.jpg", 
     bg="white", antialias = "default", width = 6, height = 3, units = "in", res = 100)
barplot(t(nt_cuts[,2:6]),names.arg=c(proteinorder),col=rainbow(5),legend=T,las=2,ylab='% of positions')
dev.off()


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## NT frequency of consensus nucleotide
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
consensus_threshold<-95 # 25 sequences in 505
# List which positions have frequency lower than threshold
consensus_nt_variable95<-consensus_nt[as.numeric(as.character(consensus_nt[,3])) < consensus_threshold,]
consensus_nt_variable95
nrow(consensus_nt_variable95)
# Add column with position of aminoacids
consensus_nt_variable95$AA_pos<-ceiling(as.numeric(as.character(consensus_nt_variable95$Pos))/3)

write.table(consensus_nt_variable95, "./variable_positions_nt_freq_lower95.tsv", sep="\t",quote=F, col.names=TRUE, row.names=F)


#####################################################################################3
# Get list of nucleotide positions at codons that have a change of aminoacid
# From  the list of variable nucleotide positions (consensus_nt_maf) 
df_nt_nonsyn<-data.frame(stringsAsFactors = FALSE)
for(i in consensus_nt_variable95$WTpos_nt){
  aa_pos<-gsub("_[^_]*$","",i) # remove last _number
  freq_aa<-as.numeric(as.character(consensus_aa[which(consensus_aa[,2]== aa_pos),"WTvalue_aa"]))
  if(freq_aa<100){
    df_nt_nonsyn<-rbind(df_nt_nonsyn,freq_matrixnt[which(freq_matrixnt[,1]==paste0(aa_pos,"_1")),],stringsAsFactors = FALSE)
    df_nt_nonsyn<-rbind(df_nt_nonsyn,freq_matrixnt[which(freq_matrixnt[,1]==paste0(aa_pos,"_2")),],stringsAsFactors = FALSE)
    df_nt_nonsyn<-rbind(df_nt_nonsyn,freq_matrixnt[which(freq_matrixnt[,1]==paste0(aa_pos,"_3")),],stringsAsFactors = FALSE)
      }
  
}
colnames(df_nt_nonsyn)<-colnames(freq_matrixnt)

df_nt_nonsyn

write.table(df_nt_nonsyn, "./codons_nonsynonymous_positions.tsv",  sep="\t",quote=F, col.names=TRUE, row.names=F)

# STATISTICAL TEST OF  DIFFERENCE BETWEEN FREQUENCY OF CONSENSUS IN EACH protein
# Kruskal Wallis test 
kruskal.test(as.numeric(consensus_nt[,3])~gsub('_.*','',consensus_nt[,2]))
pairwise.wilcox.test(as.numeric(consensus_nt[,3]),gsub('_.*','',consensus_nt[,2]),
                     p.adjust.method = "BH")
boxplot(as.numeric(consensus_nt[,3])~gsub('_.*','',consensus_nt[,2]))

# ANOVA
summary(aov(as.numeric(consensus_nt[,3])~gsub('_.*','',consensus_nt[,2])))

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
##ENTROPY** NT level
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
# Calculate entropy (R package entropy) using the function in source file /R_functions_Zika.R 
library(entropy)
proteinorder<-c("C","PR","M","E","NS1","NS2A","NS2B","NS3","NS4A","pep2K","NS4B","NS5")

source("./R_functions_Zika.R")
entropy_nt<-f_entropy(inputnt[begin_genome_nt:end_genome_nt])
colnames(entropy_nt)

sum(as.numeric(as.character(entropy_nt$Entropy)))/length(entropy_nt$Entropy) #average protease entropy

write.table(entropy_nt, "./entropy_nt.tsv", sep="\t",quote=F, col.names=TRUE, row.names=F)

# Statistics for each protein
library(dplyr)
library(gmodels)
output_ci<-ci(entropy_nt$Entropy, na.rm = TRUE) # from package gmodels
output_ci

output_ci<-group_by(entropy_nt,Protein) %>%
  dplyr::summarise(
    mean_ci = ci(Entropy, na.rm = TRUE)[1],
    CI95_low = ci(Entropy, na.rm = TRUE)[2],
    CI95_high = ci(Entropy, na.rm = TRUE)[3],
    SE = ci(Entropy, na.rm = TRUE)[4])
output_ci

# Kruskal Wallis test - difference in entropy between proteins 
kruskal.test(entropy_nt$Entropy~entropy_nt$Protein)
pairwise.wilcox.test(entropy_nt$Entropy,entropy_nt$Protein,
                     p.adjust.method = "BH")


######################################################################################################################################################
######################################################################################################################################################
# CHECK PARTICULAR MUTATIONS
#####################################################################################################################################################################
######################################################################################################################################################

# Open files
freq_matrixnt<-read.csv("./freq_matrixnt.tsv",sep='\t', header=TRUE,  colClasses = c("character"))
dim(freq_matrixnt)

freq_matrixaa<-read.csv("./freq_matrixaa.tsv",sep='\t', header=TRUE,  colClasses = c("character"))
dim(freq_matrixaa)

consensus_aa<-read.csv("./consensus_aa.tsv",sep='\t', header=TRUE,  colClasses = c("character"))
dim(consensus_aa)

consensus_nt<-read.csv("./consensus_nt.tsv",sep='\t', header=TRUE,  colClasses = c("character"))
dim(consensus_nt)

inputnt <-read.csv("./inputnt_total.csv",sep=',', header=TRUE,  colClasses = c("character"))
dim(inputnt)

metadata <-read.csv('./Metadata_505_with_results.tsv',sep='\t')
colnames(metadata)

# Get tables for particular mutations
source("./R_functions_Zika.R")
for(i in c("C_107",
           "PR_17",
           "NS1_100",
           "NS1_194",
           "NS1_349",
           "NS3_572",
           "NS3_584",
           "NS5_114",
           "NS5_267",
           "NS5_322",
           "NS5_525",
           "NS5_525",
           "NS5_642",
           "NS5_808",
           "NS5_878")){ 
f_tables_mutation("NAMING_PROPOSAL",i)

}

source("./R_functions_Zika.R")
f_tables_mutation("NAMING_PROPOSAL","C_110")

f_tables_mutation("VirusType","M_47")

#######################################################################################
#####################################################################################
# get aminoacid changes in specific proteins
##########################################################################################

# Protein PR
df_protein<-consensus_aa[which((gsub('_.*','',consensus_aa[,2])=="PR") & as.numeric(as.character(consensus_aa[,3])) < 100) ,]
table(df_protein[,3])

nrow(df_protein)

ids_aa<-inputaa[,1]
gr_aa<-metadata[match(ids_aa,as.character(metadata[,"accession"])),"NAMING_PROPOSAL"] # Get list of group names (from metadata file), in the order found in input
d_aa<-cbind(inputaa,gr_aa)
tail(colnames(d_aa))
d_aa$gr_aa<-as.factor(d_aa$gr_aa)

sink(file="./proteinPR_AA_mutations.txt",type="output")
for(i in df_protein[,2]){
  print(freq_matrixaa[which(freq_matrixaa[,1]==i),])
  print(as.matrix(table(d_aa$gr_aa,d_aa[,i])))
}
sink()



#########################################################################
#########################################################################
# Check mutations differentiating two groups - exclusively present in one group and not the other, and vice-versa
#########################################################################
colnames(metadata)   
f_aa_mutation_different_group("VirusType")
f_nt_mutation_different_group("VirusType")

# Host Mosquito vs Human
table(metadata$Host)
metadata$Host_human_mosquito<-as.character(metadata$Host)
metadata$Host_human_mosquito[metadata$Host_human_mosquito == "Monkey"]<-NA
table(metadata$Host_human_mosquito)

source("./R_functions_Zika.R")
f_aa_mutation_different_group("Host_human_mosquito")
f_nt_mutation_different_group("Host_human_mosquito")

# PreAm vs Am
table(metadata$NAMING_PROPOSAL)
metadata$Am_PreAm<-as.character(metadata$NAMING_PROPOSAL)
metadata$Am_PreAm[(metadata$Am_PreAm != "ZB.2.0") & (metadata$Am_PreAm != "ZB.2.1") & 
                    (metadata$Am_PreAm != "ZB.2.2")]<-"PreAm"
metadata$Am_PreAm[(metadata$Am_PreAm == "ZB.2.0") | (metadata$Am_PreAm == "ZB.2.1") | 
                    (metadata$Am_PreAm == "ZB.2.2")]<-"Am"
table(metadata$Am_PreAm)

source("./R_functions_Zika.R")
f_aa_mutation_different_group("Am_PreAm")
f_nt_mutation_different_group("Am_PreAm")

# Continent
metadata$American_continent<-as.character(metadata$Continent)
metadata$American_continent[(metadata$American_continent != "NAM") & (metadata$American_continent != "SA")]<-"Others"
metadata$American_continent[(metadata$American_continent == "NAM") | (metadata$American_continent == "SA")]<-"America"

table(metadata$American_continent)
source("./R_functions_Zika.R")
f_aa_mutation_different_group("American_continent")
f_nt_mutation_different_group("American_continent")



####################################################################################
####################################################################################
# PLOT MUTATION FREQUENCIES
####################################################################################
########################################################################
consensus_aa<-read.csv("./consensus_aa.tsv",sep='\t', header=TRUE,  colClasses = c("character"))
dim(consensus_aa)
colnames(consensus_aa)

positions_protein_aa<-data.frame(xstart = c(1,123,216,291,795,1147,1373,1503,2120,2247,2270,2521), 
                                 xend = c(122,215,290,794,1146,1372,1502,2119,2246,2269,2520,3423), col = rep(c("odd","even"),6))

consensus_nt<-read.csv("./consensus_nt.tsv",sep='\t', header=TRUE,  colClasses = c("character"))
dim(consensus_nt)
colnames(consensus_nt)

positions_protein_nt<-data.frame(xstart = c(1,367,646,871,2383,3439,4117,4507,6358,6739,6808,7561), 
                                 xend = c(366,645,870,2382,3438,4116,4506,6357,6738,6807,7560,10269), col = rep(c("odd","even"),6))

# Freq matrix
freq_matrixaa<-read.csv("./freq_matrixaa.tsv",sep='\t', header=TRUE,  colClasses = c("character"))
dim(freq_matrixaa)


# Metadata
metadata <-read.csv('./Metadata_505_with_results.csv',sep='\t')
colnames(metadata)


# AA

# plot frequency aa
p_aa<-ggplot() +
  geom_rect(data = positions_protein_aa, aes(xmin = xstart, xmax = xend, ymin = 0, ymax = Inf,fill = col),  color="grey", alpha = 0.4) +
  geom_tile() + # remove the margins
  scale_x_continuous(expand=c(0,0)) + # remove the margins
  scale_y_continuous(expand=c(0,0)) + # remove the margins
  scale_fill_manual(values=c("grey","white"),guide="none") +
  geom_bar(data = consensus_aa, aes(x=as.numeric(as.character(Pos)),
                                    y=100-as.numeric(as.character(WTvalue_aa))), 
           color="grey50",stat="identity") +
  ylim(0,45)+
  geom_text(data=consensus_aa %>% filter(as.numeric(as.character(WTvalue_aa))<95), # Choose max frequency of major variant to plot label of aminoacid + position
            aes(label=paste0(Pos,""),x=as.numeric(as.character(Pos)),
                y=100-as.numeric(as.character(WTvalue_aa))+1),angle=0,hjust=1) +
  labs(x="",y="Amino acid\nfrequency (%)",tag="C")  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position="none")

# map protein 
p_aa_map<-ggplot() +
  geom_rect(data = positions_protein_aa, aes(xmin = xstart, xmax = xend, ymin = -3, ymax = -0, fill = col), color="black", alpha = 0.4) +
  geom_tile() + # remove the margins
  scale_x_continuous(expand=c(0,0)) + # remove the margins
  scale_y_continuous(expand=c(0,0)) + # remove the margins
  annotate("text",label=c("C","Pr","M","E","NS1","NS2A","NS2B","NS3","NS4A","","NS4B","NS5"),
           x=positions_protein_aa$xstart+(positions_protein_aa$xend-positions_protein_aa$xstart)/2,
           y=rep(-1.5,12)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.position="none") +
  scale_fill_manual(values=c("grey","white"),guide="none")

# NT

# plot frequency nt 
p_nt<-ggplot() +
  geom_rect(data = positions_protein_nt, aes(xmin = xstart, xmax = xend, ymin = 0, ymax = Inf,fill = col),  color="grey", alpha = 0.4) +
  geom_tile() + # remove the margins
  scale_x_continuous(expand=c(0,0)) + # remove the margins
  scale_y_continuous(expand=c(0,0)) + # remove the margins
  scale_fill_manual(values=c("grey","white"),guide="none") +
  geom_bar(data = consensus_nt, aes(x=as.numeric(as.character(Pos)),
                                    y=100-as.numeric(as.character(WTvalue_nt))), 
           color="grey50",stat="identity") +
  labs(x="",y="Nucleotide\nfrequency (%)",tag="B")  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position="none")

p_nt_map<-ggplot() +
  geom_rect(data = positions_protein_nt, aes(xmin = xstart, xmax = xend, ymin = -3, ymax = -0, fill = col), color="black", alpha = 0.4) +
  geom_tile() + # remove the margins
  scale_x_continuous(expand=c(0,0)) + # remove the margins
  scale_y_continuous(expand=c(0,0)) + # remove the margins
  annotate("text",label=c("C","Pr","M","E","NS1","NS2A","NS2B","NS3","NS4A","","NS4B","NS5"),
           x=positions_protein_nt$xstart+(positions_protein_nt$xend-positions_protein_nt$xstart)/2,
           y=rep(-1.5,12)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.position="none") +
  scale_fill_manual(values=c("grey","white"),guide="none")


##############################################
# COVERAGE per position NT (proportion of sequences with nucleotide information - not missing)
length_region<-length(min(begin_genome_nt:end_genome_nt):max(begin_genome_nt:end_genome_nt))
coverage<-apply(inputnt[,begin_genome_nt:end_genome_nt],2,function(x) prop.table(table(!is.na(x)))[names(prop.table(table(!is.na(x))))==TRUE])
boe<-data.frame(cbind(1:length_region,unlist(coverage)))
colnames(boe)<-c('pos','value')

min(boe$value)
max(boe$value)
1-min(boe$value)

p_coverage <-ggplot() +
  geom_rect(data = positions_protein_nt, aes(xmin = xstart, xmax = xend, ymin = 0, ymax = Inf,fill = col),  color="grey", alpha = 0.4) +
  geom_tile() + # remove the margins
  scale_x_continuous(expand=c(0,0)) + # remove the margins
  scale_y_continuous(expand=c(0,0)) + # remove the margins
  scale_fill_manual(values=c("grey","white"),guide="none") +
  geom_line(data = boe, aes(x=pos, y=value*100), 
            color="grey50") +
  labs(x="",y="Coverage (%)",tag="A")  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position="none") 


jpeg(file = "./MutationPrevalence_GenomeWide_AA_NT.jpg", 
     bg="white", antialias = "default", width = 14, height =7, units = "in", res = 100)

plot <- plot_grid(p_aa_map,p_coverage,p_nt,p_aa, align = "v", ncol = 1, axis = "tb", rel_heights = c(0.5,3,3,3))
plot

dev.off()

####################################################################################
####################################################################################
# PLOT ENTROPY
####################################################################################
####################################################################################

entropy_aa<-read.csv("./entropy_aa.tsv",sep='\t', header=TRUE,  colClasses = c("character"))
colnames(entropy_aa)

positions_protein_aa<-data.frame(xstart = c(1,123,216,291,795,1147,1373,1503,2120,2247,2270,2521), 
                                 xend = c(122,215,290,794,1146,1372,1502,2119,2246,2269,2520,3423), col = rep(c("odd","even"),6))

entropy_nt<-read.csv("./entropy_nt.tsv",sep='\t', header=TRUE,  colClasses = c("character"))
colnames(entropy_nt)

positions_protein_nt<-data.frame(xstart = c(1,367,646,871,2383,3439,4117,4507,6358,6739,6808,7561), 
                                 xend = c(366,645,870,2382,3438,4116,4506,6357,6738,6807,7560,10269), col = rep(c("odd","even"),6))


# plot entropy aa
p_aa<-ggplot() +
  geom_rect(data = positions_protein_aa, aes(xmin = xstart, xmax = xend, ymin = 0, ymax = Inf,fill = col),  color="grey", alpha = 0.4) +
  geom_tile() + # remove the margins
  scale_x_continuous(expand=c(0,0)) + # remove the margins
  scale_y_continuous(expand=c(0,0)) + # remove the margins
  scale_fill_manual(values=c("grey","white"),guide="none") +
  geom_bar(data = entropy_aa, aes(x=as.numeric(as.character(Pos)),
                                  y=as.numeric(as.character(Entropy))), 
           color="grey50",stat="identity") +
  ylim(0,0.75)+
  geom_text(data=entropy_aa %>% filter(as.numeric(as.character(Entropy))>0.2), # Choose max frequency of major variant to plot label of aminoacid + position
            aes(label=paste0(Pos,""),x=as.numeric(as.character(Pos)),
                y=as.numeric(as.character(Entropy))+0.01),angle=0,hjust=1) +
  labs(x="",y="Entropy",tag="B")  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position="none")

# map protein 
p_aa_map<-ggplot() +
  geom_rect(data = positions_protein_aa, aes(xmin = xstart, xmax = xend, ymin = -3, ymax = -0, fill = col), color="black", alpha = 0.4) +
  geom_tile() + # remove the margins
  scale_x_continuous(expand=c(0,0)) + # remove the margins
  scale_y_continuous(expand=c(0,0)) + # remove the margins
  annotate("text",label=c("C","Pr","M","E","NS1","NS2A","NS2B","NS3","NS4A","","NS4B","NS5"),
           x=positions_protein_aa$xstart+(positions_protein_aa$xend-positions_protein_aa$xstart)/2,
           y=rep(-1.5,12)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.position="none") +
  scale_fill_manual(values=c("grey","white"),guide="none")

# plot entropy nt
p_nt<-ggplot() +
  geom_rect(data = positions_protein_nt, aes(xmin = xstart, xmax = xend, ymin = 0, ymax = Inf,fill = col),  color="grey", alpha = 0.4) +
  geom_tile() + # remove the margins
  scale_x_continuous(expand=c(0,0)) + # remove the margins
  scale_y_continuous(expand=c(0,0)) + # remove the margins
  scale_fill_manual(values=c("grey","white"),guide="none") +
  geom_bar(data = entropy_nt, aes(x=as.numeric(as.character(Pos)),
                                  y=as.numeric(as.character(Entropy))), 
           color="grey50",stat="identity") +
  labs(x="",y="Entropy",tag="A")  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position="none")

# map protein 
p_nt_map<-ggplot() +
  geom_rect(data = positions_protein_nt, aes(xmin = xstart, xmax = xend, ymin = -3, ymax = -0, fill = col), color="black", alpha = 0.4) +
  geom_tile() + # remove the margins
  scale_x_continuous(expand=c(0,0)) + # remove the margins
  scale_y_continuous(expand=c(0,0)) + # remove the margins
  annotate("text",label=c("C","Pr","M","E","NS1","NS2A","NS2B","NS3","NS4A","","NS4B","NS5"),
           x=positions_protein_nt$xstart+(positions_protein_nt$xend-positions_protein_nt$xstart)/2,
           y=rep(-1.5,12)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.position="none") +
  scale_fill_manual(values=c("grey","white"),guide="none")


jpeg(file = "./Entropy_AA_NT.jpg", 
     bg="white", antialias = "default", width = 14, height = 5, units = "in", res = 100)

plot <- plot_grid(p_aa_map,p_nt,p_aa, align = "v", ncol = 1, axis = "tb", rel_heights = c(0.5, 3, 3))
plot

dev.off()



############################################################################
############################################################################
# 3. GENETIC DISTANCES *******************************************************
############################################################################
############################################################################

## Obtain matrix of pairwise genetic distances between groups
fasta_file_path="./Zika_Aligned_505.fasta"
metadata_file_path="./Metadata_505_with_results.tsv"
acession_name="accession"
group_name="NAMING_PROPOSAL"
output_file="./mean_distance_matrix_505_"

# check files
fasta_505<-read.FASTA(fasta_file_path,type="DNA")
fasta_505
metadata_505<-read.csv(metadata_file_path,sep='\t')
colnames(metadata_505)
prop.table(table(metadata_505$VirusType))

source("./R_functions_Zika.R") 
f_genetic_distances_groups(fasta_file_path,metadata_file_path,acession_name,group_name,output_file)


# Distribution of Pairwise Genetic distances for each group

f_plot_hist_distances<-function(group){
  names_gr<-metadata_505[which(metadata_505[,"NAMING_PROPOSAL"]==group),"accession"]
  table(metadata_505$NAMING_PROPOSAL)
  new_fasta<-fasta_505[names(fasta_505) %in% names_gr]
  # Get distance matrix (lower matrix)
  matrix_dist<-dist.dna(new_fasta, model="K80",pairwise.deletion = TRUE)
  df_distances<-na.omit(reshape2::melt(as.matrix(matrix_dist)))
  # Plot histogram of distances
  ggplot(data= df_distances, aes(x=value)) + 
    geom_histogram(bins=40) +
    xlab("Genetic distance") + ylab("Frequency") +
    xlim(-0.001,0.010) +
    labs(title=group)
} 

p1<-f_plot_hist_distances("ZA.1.1")
p2<-f_plot_hist_distances("ZA.1.2")
p3<-f_plot_hist_distances("ZB.1.0")
p4<-f_plot_hist_distances("ZB.1.1")
p5<-f_plot_hist_distances("ZB.1.2")
p6<-f_plot_hist_distances("ZB.2.0")
p7<-f_plot_hist_distances("ZB.2.1")
p8<-f_plot_hist_distances("ZB.2.2")

jpeg(file = "./hist_genetic_distances_groups_American.jpeg", bg="white", 
     antialias = "default", width = 6, height = 6, units = "in", res = 300)
grid.newpage()
grid.draw(rbind(ggplotGrob(p6),
                ggplotGrob(p7), ggplotGrob(p8),
                size = "last"))

dev.off() 

###############################################################
# 4. AMOVA - analysis of molecular variance
###############################################################
# based on https://popgen.nescent.org/PopDiffSequenceData.html

# Test genetic differentiation between groups
# 505 dataset without sample LC219720 (ZB1.1/ZB2.0)
################
fasta_file_path="./Zika_Aligned_505_without_LC219720.fasta"
metadata_file_path="./Metadata_505_with_results_without_LC219720.csv"
acession_name="accession"
group_name="NAMING_PROPOSAL"

source("./R_functions_Zika.R") 
f_amova(fasta_file_path,metadata_file_path,acession_name,group_name)

# Test genetic differentiation between groups of ZB clade only
# ZB dataset without sample LC219720 (ZB1.1/ZB2.0)
################
fasta_file_path="./Zika_Aligned_437_withoutLC219720.fasta"
metadata_file_path="./Metadata_NAMING_ZB.csv"
acession_name="fasta_id"
group_name="NAMING_PROPOSAL"

source("./R_functions_Zika.R") 
f_amova(fasta_file_path,metadata_file_path,acession_name,group_name)


###############################################################
# 5. RhierBAPS - bayesian clustering
###############################################################

snp.matrix <- load_fasta("./Zika_Aligned_505.fasta")
dimnames(snp.matrix)

######################################
# RUN ANALYSIS
# 10 independent runs
#############################################
set.seed(NULL)
snp.matrix <- load_fasta("./Zika_Aligned_505.fasta")
df<-data.frame(matrix(nrow=505,ncol=0))
for(i in seq(1,10)){
  hb.results <- hierBAPS(snp.matrix, max.depth = 2, n.pops = 20, quiet = TRUE)
  df<-cbind(df,hb.results$partition.df)
}
write.table(df,"./hb_results_505_2levels_runs.tsv", sep="\t", quote=F, row.names = F)

set.seed(NULL)
snp.matrix_437 <- load_fasta("./Zika_Aligned_437.fasta")
df<-data.frame(matrix(nrow=437,ncol=0))
for(i in seq(1,10)){
  hb.results <- hierBAPS(snp.matrix_437, max.depth = 2, n.pops = 20, quiet = TRUE)
  df<-cbind(df,hb.results$partition.df)
}
write.table(df,"./hb_results_437_2levels_runs.tsv", sep="\t", quote=F, row.names = F)


##########################################################
##########################################################
# 6. PLOTS RESULTS OF FUBAR method for detection of selection
##########################################################
# FUBAR  - Hyphy - Datamonkey   
# alpha = dS # synonymous substitution rate
# beta = dN # non-synonymous substitution rate
# beta.alpha = beta - alpha
# Prob.alpha_higher_beta = Prob[alpha>beta] # NEGATIVE SELECTION
# Prob.alpha_lower_beta = Prob[alpha<beta] # POSITIVE DIVERSIFYING SELECTION
# BayesFactor.alpha_lower_beta = BayesFactor[alpha<beta] # POSITIVE

fubar_505<-read.csv('./datamonkey-table_505.csv')
colnames(fubar_505)<-c("Site","Partition","alpha","beta","beta.alpha",
                       "Prob.alpha_higher_beta","Prob.alpha_lower_beta","BayesFactor.alpha_lower_beta")
fubar_437<-read.csv('./datamonkey-table_437.csv')
colnames(fubar_437)<-c("Site","Partition","alpha","beta","beta.alpha",
                       "Prob.alpha_higher_beta","Prob.alpha_lower_beta","BayesFactor.alpha_lower_beta")
fubar_83<-read.csv('./datamonkey-table_83.csv')
colnames(fubar_83)<-c("Site","Partition","alpha","beta","beta.alpha",
                      "Prob.alpha_higher_beta","Prob.alpha_lower_beta","BayesFactor.alpha_lower_beta")

col<-ifelse(fubar_505$Prob.alpha_lower_beta > 0.9,'#D55E00',
            ifelse(fubar_505$Prob.alpha_higher_beta> 0.9,"#56B4E9",'black'))
table(col)
p_505<-ggplot() +
  geom_rect(data = positions_protein_aa, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf,fill = col),  color="grey", alpha = 0.4) +
  geom_tile() + # remove the margins
  scale_x_continuous(expand=c(0,0)) + # remove the margins
  scale_y_continuous(expand=c(0,0)) + # remove the margins
  scale_fill_manual(values=c("grey","white"),guide="none") +
  geom_bar(data = fubar_505, 
           aes(x=as.numeric(as.character(Site)),y=as.numeric(as.character(beta.alpha))), 
           stat="identity",fill=col, width=3) +
  geom_text(data=fubar_505 %>% filter(as.numeric(as.character(Prob.alpha_lower_beta))>0.9), # Choose max frequency of major variant to plot label of aminoacid + position
            aes(label=paste0(Site,""),x=as.numeric(as.character(Site)),
                y=as.numeric(as.character(beta.alpha))-1),angle=0,hjust=1, color="#D55E00",cex=3) +
  labs(x="",y="Mean posterior\nbeta - alpha",tag="A")  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position="none")

col<-ifelse(fubar_437$Prob.alpha_lower_beta > 0.9,'#D55E00',
            ifelse(fubar_437$Prob.alpha_higher_beta> 0.9,"#56B4E9",'black'))
p_437<-ggplot() +
  geom_rect(data = positions_protein_aa, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf,fill = col),  color="grey", alpha = 0.4) +
  geom_tile() + # remove the margins
  scale_x_continuous(expand=c(0,0)) + # remove the margins
  scale_y_continuous(expand=c(0,0)) + # remove the margins
  scale_fill_manual(values=c("grey","white"),guide="none") +
  geom_bar(data = fubar_437, 
           aes(x=as.numeric(as.character(Site)),y=as.numeric(as.character(beta.alpha))), 
           stat="identity",fill=col, width=3) +
  geom_text(data=fubar_437 %>% filter(as.numeric(as.character(Prob.alpha_lower_beta))>0.9), # Choose max frequency of major variant to plot label of aminoacid + position
            aes(label=paste0(Site,""),x=as.numeric(as.character(Site)),
                y=as.numeric(as.character(beta.alpha))-0.5),angle=0,hjust=1, color="#D55E00",cex=3) +
  labs(x="",y="Mean posterior\nbeta - alpha",tag="B")  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position="none")

col<-ifelse(fubar_83$Prob.alpha_lower_beta > 0.9,'#D55E00',
            ifelse(fubar_83$Prob.alpha_higher_beta> 0.9,"#56B4E9",'black'))
p_83<-ggplot() +
  geom_rect(data = positions_protein_aa, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf,fill = col),  color="grey", alpha = 0.4) +
  geom_tile() + # remove the margins
  scale_x_continuous(expand=c(0,0)) + # remove the margins
  scale_y_continuous(expand=c(0,0)) + # remove the margins
  scale_fill_manual(values=c("grey","white"),guide="none") +
  geom_bar(data = fubar_83, 
           aes(x=as.numeric(as.character(Site)),y=as.numeric(as.character(beta.alpha))), 
           stat="identity",fill=col, width=3) +
  geom_text(data=fubar_83 %>% filter(as.numeric(as.character(Prob.alpha_lower_beta))>0.9), # Choose max frequency of major variant to plot label of aminoacid + position
            aes(label=paste0(Site,""),x=as.numeric(as.character(Site)),
                y=as.numeric(as.character(beta.alpha))-2),angle=0,hjust=1, color="#D55E00",cex=3) +
  labs(x="",y="Mean posterior\nbeta - alpha",tag="C")  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position="none")

# map protein 
p_aa_map<-ggplot() +
  geom_rect(data = positions_protein_aa, aes(xmin = xstart, xmax = xend, ymin = -3, ymax = -0, fill = col), color="black", alpha = 0.4) +
  geom_tile() + # remove the margins
  scale_x_continuous(expand=c(0,0)) + # remove the margins
  scale_y_continuous(expand=c(0,0)) + # remove the margins
  annotate("text",label=c("C","Pr","M","E","NS1","NS2A","NS2B","NS3","NS4A","","NS4B","NS5"),
           x=positions_protein_aa$xstart+(positions_protein_aa$xend-positions_protein_aa$xstart)/2,
           y=rep(-1.5,12)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.position="none") +
  scale_fill_manual(values=c("grey","white"),guide="none")

jpeg(file = "./plot_FUBAR_datasets.jpg", 
     bg="white", antialias = "default", width = 12, height =6, units = "in", res = 300)
plot <- plot_grid(p_aa_map,p_505,p_437,p_83, align = "v", ncol = 1, axis = "tb", rel_heights = c(0.5, 3, 3,3))
plot

dev.off()

################################################
################################################
# PLOT Frequency of amino acids per group - positive selection
################################################
################################################
ids_aa<-inputaa[,1]
gr<-metadata[match(ids_aa,as.character(metadata[,"accession"])),"NAMING_PROPOSAL"] # Get list of group names (from metadata file), in the order found in input
d_aa<-cbind(inputaa,gr)
dim(d_aa)

d_aa$gr<-as.character(d_aa$gr)
d_aa<-d_aa[-which(d_aa[,"gr"]=="ZB1.1 / ZB2.0"),] # remove Vietnam  sample ZB1.1 / ZB2.0
d_aa$gr<-as.factor(d_aa$gr)

plot_function<-function(d_aa,pos){
  m<-as.data.frame(table(d_aa$gr,d_aa[,pos]))
  pos_aa<-which(freq_matrixaa[,"position"]==pos)
  ggplot(data= m, aes(y=Freq, x=Var1)) + 
    geom_bar(aes(fill=Var2), colour="black",position="fill", stat="identity") +
    xlab("") + ylab("") + ggtitle(paste0(pos," (position ",pos_aa,")"))+
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(values=c("#a6611a",
                               "#dfc27d","#80cdc1","#018571"),name="Amino acid") 
}
list_positions<-c("NS1_100","NS4B_24","NS4B_186","NS5_74","NS5_287",
                  "NS5_642","NS5_703","NS5_808","NS5_878")

allplots<-lapply(list_positions,plot_function,d_aa=d_aa)

jpeg(file = "./Frequency_AA_positive_selection.jpg", 
     bg="white", antialias = "default", width = 6, height =12, units = "in", res = 100)
grid.arrange(grobs=allplots, nrow=length(list_positions),ncol=1)
dev.off()



