library(dplyr)
library(tm)
library(ggplot2)
library(forcats)
library(plotrix)
library(matrixStats)
library(remotes)
library("multcompView")
#remotes::install_github("snandi/RFunctionsSN")
library(RFunctionsSN)
options(dplyr.summarise.inform = FALSE)


path <- "~/storage/qPCR Analysis/Phenanthrene qPCR"
output_dir = paste0(path,"/output/")
unlink(output_dir, recursive = TRUE)
reps <- list.dirs(path,recursive = FALSE)
rep <-0
CTdf <- data.frame()


for (reppath in reps){ 
  rep <- rep+1
  plates <- list.dirs(reppath,recursive = FALSE) 
  rep_result_df <- data_frame("chemical" = character(0), "gene" = character(0), "RQ"  = numeric(0))
  
  #for each of the plates available 
  for (platepath in plates) {
    files <- list.files(platepath, pattern=NULL, all.files=FALSE,full.names=FALSE)
    platesoutput <- data.frame(value = character(),well = character())
    
    #match value in file type to well to be able to link together later
    for (file in files) { 
      dftemp = data_frame()
      table = read.csv(paste0(platepath,"/",file),sep=",",header=FALSE,fill = TRUE)
      colnames(table) <- table[1,]
      for(row in 2:nrow(table)){ 
        for (col in 2:ncol(table)) { 
          dftemp = rbind(dftemp, c(table[row,col], paste0(table[row,1],table[1,col])))
        }
      }
      
      #generate plate data given file type inputs 
      if (grepl("-CTmap", file, ignore.case = TRUE)){ 
        colnames(dftemp) <- c("CTvalue", "well")
      } else if (grepl("-wellmap", file, ignore.case = TRUE)) { 
        colnames(dftemp) <- c("chemical", "well")
      } else if (grepl("-genemap", file, ignore.case = TRUE)) { 
        colnames(dftemp) <- c("gene", "well")
      }
      if (nrow(platesoutput) ==0){ 
        platesoutput <- dftemp
      } else { 
        platesoutput <- left_join(platesoutput,dftemp, by="well")
      }
      
    }
    #find technical replicates average for plate 
    platesoutput <- platesoutput[complete.cases(platesoutput), ]
    platesoutput$CTvalue<- as.numeric(platesoutput$CTvalue)
    average_df <- platesoutput %>%
      group_by(gene, chemical) %>%
      summarize(average_CT = mean(CTvalue, na.rm = TRUE))
    
    #find control gene (usually b-actin) and calculate DeltaCT
    bactin <- filter(average_df, gene == 'b-actin')
    gene_of_interest <- filter(average_df, gene != 'b-actin')
    joined <- left_join(bactin, gene_of_interest, by ='chemical')
    joined$DeltaCT <- joined$average_CT.y-joined$average_CT.x
    
    #find control treatment samples and calculate Delta Delta CT and RQ
    control_to_check <- "MIR + veh con "
    joined <- joined %>%
      mutate(control = ifelse(chemical ==control_to_check, "C", ""))
    joined <- joined %>% 
      group_by(gene.y) %>%
      mutate(DDCT = DeltaCT-DeltaCT[control=="C"]) 
    joined$RQ <- round(2^(-joined$DDCT), digits=2) 
    names(joined)[names(joined) == "gene.y"] <- "gene"
    rep_result_df <- rbind(rep_result_df, joined[,c("chemical","gene","RQ")])
  }
   
   rep_result_df <- rep_result_df %>% 
     group_by(chemical,gene) %>%
     summarize(RQ = mean(RQ))
  if(nrow(CTdf)==0){ 
    CTdf <- bind_rows(CTdf,rep_result_df)
  }
  else{ 
    #checks to see if there are any genes tested for in the rep that do not have a previous rep 
    #if it is new genes and not a rep then it adds to the table 
    RQ_rep<- paste0("RQ",as.character(rep))
    names(rep_result_df)[names(rep_result_df) == "RQ"] <- RQ_rep
    nonmatching <-anti_join(rep_result_df, CTdf, by= c("chemical", "gene"))
    if(nrow(nonmatching)!=0) {
      CTdf <- bind_rows(CTdf,nonmatching)
      rep_result_df <- rep_result_df[!rep_result_df$chemical %in% nonmatching$chemical, ]
    }
    #if this is a new rep it adds the column 
    if(all(rep_result_df$gene %in% CTdf$gene) && all(rep_result_df$chemical %in% CTdf$chemical)){
      CTdf<- left_join(CTdf, rep_result_df, by= c("chemical", "gene"))
    }
    columns_to_merge <- grep(paste0("^",RQ_rep), names(CTdf), value = TRUE)
    # Merge selected columns into a new column
     
    # Remove the original columns if needed
    if(length(columns_to_merge)>1){ 
      RQ_toadd<- coalesce(!!!c(CTdf[columns_to_merge]))
      RQ_tobind <- data.frame(value = RQ_toadd)
      colnames(RQ_tobind) <- RQ_rep
      CTdf <- CTdf[, !names(CTdf) %in% columns_to_merge]
      CTdf<-bind_cols(CTdf,RQ_tobind) 
    }
    
    
  }
   
}

#group by chemical and dose and calculate avg between biological replicates and standard deviation
CTdf<- CTdf[!CTdf$chemical=="MIR + veh con ",]

CTdf$chemicalgroup <- gsub(control_to_check, "", removeNumbers(CTdf$chemical))
CTdf$chemicalgroup <- gsub("\\.", "", CTdf$chemicalgroup)
CTdf$chemicalgroup <- trimws(CTdf$chemicalgroup)
dir.create(output_dir)

CTdf$average_RQ <- round(rowMeans(CTdf[, grepl("^RQ", names(CTdf))], na.rm = TRUE), digits=2) 
CTdf$SEM <- apply(CTdf[, grepl("^RQ", names(CTdf))], 1, function(row) {
  std.error(na.omit(row))
})
unique_chemicalgroups <- unique(CTdf$chemicalgroup)
unique_genes <- unique(CTdf$gene)


write.table(CTdf,paste0(output_dir,"resall.txt"),sep="\t", row.names=FALSE)


