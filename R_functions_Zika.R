######################################################################################################
# Source file with functions used in R_script_Zika_diversity.R
######################################################################################################

####################################################
# Function to make frequency matrix of AA - use dataset inputaa[,begin_genome_aa:end_genome_aa]
####################################################
makeMatrixAA<-function(dataset) {
  names<-c("position", "support", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "\\*")
  #serotype_numb<-matrix(ncol=length(names),nrow=(length(names(dataset))-1),dimnames =list(1:(ncol(dataset)-1),names),data=0)
  serotype_perc<-matrix(ncol=length(names),nrow=(length(names(dataset))),dimnames =list(1:(ncol(dataset)),names),data=0)
  
  for (k in 1:ncol(dataset)){ # per column # start at 2 to ignore "seqid" name
    serotype_perc[k,1]<-names(dataset)[k]
    serotype_perc[k,2]<-length(na.omit(dataset[,k]))   
    tab<-round(prop.table(table(dataset[,k]))*100,2)   # table 
    namestab<-names(tab)                # table headers 
    #print(k)
    #print(namestab)
    for(l in 1:dim(tab))
    {
      sep_aa<-unlist(strsplit(namestab[l],""))
      #per single 
      for (su_aa in sep_aa)
      {
        serotype_perc[k,match(su_aa,colnames(serotype_perc))] <-as.numeric(serotype_perc[k,match(su_aa,colnames(serotype_perc))]) + (tab[l]/length(sep_aa))[[1]]}
    }
    # as.numeric(lapply(tabl[l]/length(sep_aa)*100,formatC,digits=5,format="f")[[1]])
  }
  return(serotype_perc)
}


####################################################
# Function to make frequency matrix of NT - use dataset inputnt[,begin_genome_nt:end_genome_nt]
####################################################
makeMatrixNT<-function(dataset) {
  k <- 1 
  namesnt<-c("position", "support", "A","C","G","T")
  serotype_percnt<-matrix(ncol=length(namesnt),nrow=(length(names(dataset))),dimnames =list(1:(ncol(dataset)),namesnt),data=0)
  
  for (we in 1:ncol(dataset)) { # per column # take car seqid
    
    if(we==500*k){print(500*k)}
    #print(we)
    serotype_percnt[we,1]<-names(dataset)[we]   # first column name 
    serotype_percnt[we,2]<-length(na.omit(dataset[,we]))    # second column support 
    # following columns 
    tab<-prop.table(table(dataset[,we]))*100   # table 
    namestab<-names(tab)                # table headers 
    for(l in 1:dim(tab)) {
      serotype_percnt[we,match(namestab[l],colnames(serotype_percnt))]<-tab[l][[1]]
    }
    k<-k+1
  }
  return(serotype_percnt)
}


####################################################
# Function to Compare the consensus with reference
####################################################
function_consensus_vs_ref <- function(consensus, inputref, genome_protein_positions) { 
  nrow(consensus)
  ncol(inputref)
  rr<-cbind(consensus,as.character(inputref))
  names(rr)[5]<-'WTref'
  rrmismatch<-rr[which(rr[,4]!=rr[,5]),]
 
  print(rrmismatch)
  
}


####################################################
# Function to calculate entropy per position, format and add protein id
####################################################

f_entropy<-function(input){
    counter=0
    entropmix =NULL;info=NULL # nieuw mixtures
    for (c in 1:(ncol(input))) {
        countsmix<-table(na.omit(input[,c]))
        entropmix = c(entropmix, entropy.ChaoShen(countsmix))
        counter=counter+1
        info<-rbind(info,c(counter,(colnames(input))[c],entropy.ChaoShen(countsmix),
                                 length(input[,c][!is.na(input[,c])])))
    }
      colnames(info)<-c("Pos","Protein_pos","Entropy","N")
      # Prepare the dataframe to analyse by protein
      df<-as.data.frame(info)
      df$Protein_pos<-as.character(df$Protein_pos)
      df$Entropy<-as.numeric(as.character(df$Entropy))
      df$Protein<-NULL
      for(i in 1:nrow(df)){
        df[i,"Protein"]<-unlist(strsplit(df[i,"Protein_pos"],"_"))[1]
      } 
      df$Protein<-factor(df$Protein, levels=proteinorder)
      
  print(df)

}


###############################################################3
# FUNCTION to Calculate Genetic distances matrix
###############################################################3
f_genetic_distances_groups<-function(fasta_file_path,metadata_file_path,acession_name,group_name,output_file){
  
  # Import FASTA file
  fasta  <-read.FASTA(fasta_file_path,type="DNA")
  
  ## Import METADATA file
  metadata <-read.csv(metadata_file_path,sep='\t')
  
  # Get distance matrix
  matrix_dist_seqs<-dist.dna(fasta, model="K80",pairwise.deletion = TRUE)
  
  # Get list of group names (from metadata file), in the order found in the distance matrix
  ids<-labels(matrix_dist_seqs)
  gr<-metadata[match(ids,as.character(metadata[,acession_name])),group_name]
  
  # Get mean distance matrix within and between groups
  matrix_dist_groups<-meandist(matrix_dist_seqs,gr)
  
  print(matrix_dist_groups)
  
  # write distance matrix to file
  write.table(matrix_dist_groups, paste0(output_file,noquote(group_name),".tsv"),
              sep="\t",quote=F, col.names=TRUE, row.names=T)  
}


########################################3
# FUNCTION to perform AMOVA analysis
########################################3

f_amova<-function(fasta_file_path,metadata_file_path,acession_name,group_name){
  
  # Import fasta
  fasta  <-read.FASTA(fasta_file_path,type="DNA")
  f_dist<-dist.dna(fasta, model="K80",pairwise.deletion = TRUE) # K80 - Kimura 2-parameter model
  #f<-read.multiFASTA(fasta_file_path)
  #f_dist <- dist.multidna(f, pool = TRUE)
  
  # Create genind object
  f<-read.multiFASTA(fasta_file_path)
  g<-multidna2genind(f, mlst = FALSE)
  
  # Import metadata
  metadata <-read.csv(metadata_file_path,sep='\t')
  
  # Set the population (groupings) file
  ids<-getSequenceNames(f)
  
  groups<-metadata[match(ids[[1]],as.character(metadata[,acession_name])),group_name]
  gr<-data.frame(groups)
  
  strata(g)<-data.frame(gr)
  
  # AMOVA - test hypothesis about differentiation between subgroups
  setPop(g)<-~groups # name of column with population group in the data.frame imported to @strata
  g
  print(amova(f_dist ~ groups, data = strata(g), nperm = 1000))
  
}


########################################
# FUNCTION to Get tables for particular mutations
########################################
f_tables_mutation<-function(group_name,pos,acession_name){
  ids_aa<-inputaa[,1]
  gr_aa<-metadata[match(ids_aa,as.character(metadata[,acession_name])),group_name] # Get list of group names (from metadata file), in the order found in input
  d_aa<-cbind(inputaa,gr_aa)
  tail(colnames(d_aa))
  d_aa$gr_aa<-as.factor(d_aa$gr_aa)
  
  ids_nt<-inputnt[,1]
  gr_nt<-metadata[match(ids_nt,as.character(metadata[,acession_name])),group_name] # Get list of group names (from metadata file), in the order found in inputnt
  d_nt<-cbind(inputnt,gr_nt)
  tail(colnames(d_nt))
  d_nt$gr_nt<-as.factor(d_nt$gr_nt)
  
  # print(consensus_aa[which(consensus_aa$WTpos_aa==pos),])
  print(freq_matrixaa[which(freq_matrixaa[,1]==pos),])
  print(as.matrix(table(d_aa$gr_aa,d_aa[,pos])))
  #print(consensus_nt[which(consensus_nt$WTpos_nt==paste0(pos,"_1")),])
  #print(consensus_nt[which(consensus_nt$WTpos_nt==paste0(pos,"_2")),])
  #print(consensus_nt[which(consensus_nt$WTpos_nt==paste0(pos,"_3")),])
  
  print(freq_matrixnt[which(freq_matrixnt[,1]==paste0(pos,"_1")),])
  print(freq_matrixnt[which(freq_matrixnt[,1]==paste0(pos,"_2")),])
  print(freq_matrixnt[which(freq_matrixnt[,1]==paste0(pos,"_3")),])
  
  print(table(d_nt$gr,d_nt[,paste0(pos,"_1")]))
  print(table(d_nt$gr,d_nt[,paste0(pos,"_2")]))
  print(table(d_nt$gr,d_nt[,paste0(pos,"_3")]))
}

########################################
# FUNCTION to find Aminoacid mutation differentiating groups
########################################
f_aa_mutation_different_group<-function(group_name,acession_name){
  ids_aa<-inputaa[,1]
  gr_aa<-metadata[match(ids_aa,as.character(metadata[,acession_name])),group_name] # Get list of group names (from metadata file), in the order found in input
  d_aa<-cbind(inputaa,gr_aa)
  tail(colnames(d_aa))
  d_aa$gr_aa<-as.factor(d_aa$gr_aa)
  
  ids_nt<-inputnt[,1]
  gr_nt<-metadata[match(ids_nt,as.character(metadata[,acession_name])),group_name] # Get list of group names (from metadata file), in the order found in inputnt
  d_nt<-cbind(inputnt,gr_nt)
  tail(colnames(d_nt))
  d_nt$gr_nt<-as.factor(d_nt$gr_nt)
  
  counter<-0
  for(i in colnames(d_aa)[-length(colnames(d_aa))][-1]){ # exclude seqid and gr name column
    m<-as.matrix(table(d_aa$gr_aa,d_aa[,i]))
    if(ncol(m)==2){
      if((m[1,1]==0 & m[1,2]>0 & m[2,1]>0 & m[2,2]==0) | 
         (m[1,1]>0 & m[1,2]==0 & m[2,1]==0 & m[2,2]>0)){
        print(i)
        print(m)
        
        print(freq_matrixaa[which(freq_matrixaa[,1]==i),])
        
        print(freq_matrixnt[which(freq_matrixnt[,1]==paste0(i,"_1")),])
        print(freq_matrixnt[which(freq_matrixnt[,1]==paste0(i,"_2")),])
        print(freq_matrixnt[which(freq_matrixnt[,1]==paste0(i,"_3")),])
        
        print(table(d_nt$gr_nt,d_nt[,paste0(i,"_1")]))
        print(table(d_nt$gr_nt,d_nt[,paste0(i,"_2")]))
        print(table(d_nt$gr_nt,d_nt[,paste0(i,"_3")]))
        
        counter<-counter+1
      }
    }
    else if(ncol(m)==3){
      if((m[1,1]==0 & m[1,2]==0 & m[1,3]>0 & 
          m[2,1]>0 & m[2,2]>0 & m[2,3]==0) | 
         
         (m[1,1]==0 & m[1,2]>0 & m[1,3]==0 &
          m[2,1]>0 & m[2,2]==0 & m[2,3]>0) |
         
         (m[1,1]>0 & m[1,2]==0 & m[1,3]==0 & 
          m[2,1]==0 & m[2,2]>0 & m[2,3]>0) | 
         
         (m[1,1]==0 & m[1,2]>0 & m[1,3]>0 & 
          m[2,1]>0 & m[2,2]==0 & m[2,3]==0) | 
         
         (m[1,1]>0 & m[1,2]==0 & m[1,3]>0 & 
          m[2,1]==0 & m[2,2]>0 & m[2,3]==0) | 
         
         (m[1,1]>0 & m[1,2]>0 & m[1,3]==0 & 
          m[2,1]==0 & m[2,2]==0 & m[2,3]>0) ){
        
        print(i)
        print(m)
        
        print(freq_matrixaa[which(freq_matrixaa[,1]==i),])
        
        print(freq_matrixnt[which(freq_matrixnt[,1]==paste0(i,"_1")),])
        print(freq_matrixnt[which(freq_matrixnt[,1]==paste0(i,"_2")),])
        print(freq_matrixnt[which(freq_matrixnt[,1]==paste0(i,"_3")),])
        
        print(table(d_nt$gr_n,d_nt[,paste0(i,"_1")]))
        print(table(d_nt$gr_n,d_nt[,paste0(i,"_2")]))
        print(table(d_nt$gr_n,d_nt[,paste0(i,"_3")]))
        
        counter<-counter+1
      }   
    } 
  }
  print(counter)
}

########################################
# FUNCTION to find  Nucleotide mutation differentiating groups
########################################
f_nt_mutation_different_group<-function(group_name,acession_name){
  ids_aa<-inputaa[,1]
  gr_aa<-metadata[match(ids_aa,as.character(metadata[,acession_name])),group_name]
  d_aa<-cbind(inputaa,gr_aa)
  tail(colnames(d_aa))
  d_aa$gr_aa<-as.factor(d_aa$gr_aa)
  
  ids_nt<-inputnt[,1]
  gr_nt<-metadata[match(ids_nt,as.character(metadata[,acession_name])),group_name] # Get list of group names (from metadata file), in the order found in inputnt
  d_nt<-cbind(inputnt,gr_nt)
  tail(colnames(d_nt))
  d_nt$gr_nt<-as.factor(d_nt$gr_nt)
  
  counter<-0
  for(i in colnames(d_nt)[-length(colnames(d_nt))][-1]){ # exclude seqid and gr name column
    m<-as.matrix(table(d_nt$gr,d_nt[,i]))
    if(ncol(m)==2){
      if((m[1,1]==0 & m[1,2]>0 & m[2,1]>0 & m[2,2]==0) | 
         (m[1,1]>0 & m[1,2]==0 & m[2,1]==0 & m[2,2]>0)){
        print(i)
        print(m)
        
        print(freq_matrixnt[which(freq_matrixnt[,1]==paste0(i)),])
        
        counter<-counter+1
      }
    } 
    else if(ncol(m)==3){
      if((m[1,1]==0 & m[1,2]==0 & m[1,3]>0 & 
          m[2,1]>0 & m[2,2]>0 & m[2,3]==0) | 
         
         (m[1,1]==0 & m[1,2]>0 & m[1,3]==0 &
          m[2,1]>0 & m[2,2]==0 & m[2,3]>0) |
         
         (m[1,1]>0 & m[1,2]==0 & m[1,3]==0 & 
          m[2,1]==0 & m[2,2]>0 & m[2,3]>0) | 
         
         (m[1,1]==0 & m[1,2]>0 & m[1,3]>0 & 
          m[2,1]>0 & m[2,2]==0 & m[2,3]==0) | 
         
         (m[1,1]>0 & m[1,2]==0 & m[1,3]>0 & 
          m[2,1]==0 & m[2,2]>0 & m[2,3]==0) | 
         
         (m[1,1]>0 & m[1,2]>0 & m[1,3]==0 & 
          m[2,1]==0 & m[2,2]==0 & m[2,3]>0) ){
        
        print(i)
        print(m)
        
        print(freq_matrixnt[which(freq_matrixnt[,1]==paste0(i)),])
        
        counter<-counter+1
      }   
    } 
  }
  print(counter)
}

