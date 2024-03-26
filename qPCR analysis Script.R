library(dplyr)
library(tm)
library(ggplot2)
library(forcats)
library(plotrix)
library(matrixStats)
library(remotes)
library("multcompView")
library("plotly")
library("crosstalk")
library("scatterD3")
library(tidyr)
#remotes::install_github("snandi/RFunctionsSN")
library(RFunctionsSN)
#devtools::install_github("jcheng5/d3scatter")
library(tibble)
options(dplyr.summarise.inform = FALSE)

library(dplyr)
library("DescTools")
library("CATT")
library("coin")
library(plotrix)
library(multcomp)
library("aod")
library(outliers)

path <- "~/storage/qPCR Analysis/Phenanthrene qPCR2"
output_dir = paste0(path,"/output/")
unlink(output_dir, recursive = TRUE)
reps <- list.dirs(path,recursive = FALSE)
rep <-0
plate <-0
CTdf <- data.frame()

chemicalsymbols = c("phe","9p","bap") #if using full names already leave empty 
chemicalnames= c("Phenanthrene", "9-Chloro-Phenanthrene","BaP")
unit = "uM"

control_to_check <- "mir+vehcon"
analysistype = "qPCR"

for (reppath in reps){ 
  rep <- rep+1
  plates <- list.dirs(reppath,recursive = FALSE) 
  rep_result_df <- data_frame("chemical" = character(0), "gene" = character(0), "RQ"  = numeric(0))
  
  #for each of the plates available 
  for (platepath in plates) {
    plate<- plate+1
    #print(rep)
    #print(plate)
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
      } else if (grepl("-genemap", file, ignore.case = TRUE) || grepl("-proteinmap", file, ignore.case = TRUE)) { 
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
    platesoutput$chemical <- tolower(gsub(" ", "", platesoutput$chemical))
    platesoutput$CTvalue<- as.numeric(platesoutput$CTvalue)
    average_df <- platesoutput %>%
      group_by(gene, chemical) %>%
      summarize(average_CT = mean(CTvalue, na.rm = TRUE))
    
    #find control gene (usually b-actin) and calculate DeltaCT
    bactin <- filter(average_df, gene == 'b-actin')
    gene_of_interest <- filter(average_df, gene != 'b-actin')
    #print(platesoutput)
    joined <- left_join(bactin, gene_of_interest, by ='chemical')
    
    
    #find control treatment samples and calculate Delta Delta CT and RQ
    joined <- joined %>%
      mutate(control = ifelse(chemical ==control_to_check, "C", ""))
    
    #western processing mostly uses qPCR terminology during process just so I dont have to redo it all 
    #equaton is updated tho 
    if(analysistype=="western"){
      joined$DeltaCT <- joined$average_CT.y/joined$average_CT.x
      joined <- joined %>% 
        group_by(gene.y) %>%
        mutate(RQ = round(DeltaCT / DeltaCT[control == "C"], 2))
    }
    else{ 
      joined$DeltaCT <- joined$average_CT.y-joined$average_CT.x
      joined <- joined %>% 
        group_by(gene.y) %>%
        mutate(DDCT = DeltaCT-DeltaCT[control=="C"])
      joined$RQ <- round(2^(-joined$DDCT), digits=2) 
    }
    names(joined)[names(joined) == "gene.y"] <- "gene"
    #print(joined)
    rep_result_df <- rbind(rep_result_df, joined[,c("chemical","gene","DeltaCT","RQ")])
  }
  
  #calculate average RQ from technical replicates
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
      rep_result_df <- anti_join(rep_result_df, nonmatching, by = c("chemical","gene"))
    }
    
    #if this is a new rep it adds the column 
    CTdf<- left_join(CTdf, rep_result_df, by= c("chemical", "gene"))
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

#right now controls are removed after calculate
CTdf<- CTdf[!CTdf$chemical==control_to_check,]
dir.create(output_dir)

#outlier test
for (i in 1:nrow(CTdf)) {
  rq_columns <- names(CTdf)[startsWith(names(CTdf), "RQ")]
  row <- CTdf[i,rq_columns]
  if(rowSums(!is.na(row))>2){ 
    values_to_check_outliers = unlist(row)
    grubbs_result <- grubbs.test(values_to_check_outliers)
    matches <- gregexpr("\\b\\d+\\.\\d+\\b", grubbs_result$alternative, perl = TRUE)
    numbers <- regmatches(grubbs_result$alternative, matches)[[1]]
    outlier_to_remove = as.numeric(numbers)
    if(grubbs_result$p.value<0.05){ 
      CTdf[i,which(CTdf[i,] == outlier_to_remove)] <- NA
    }
  }
  
}

#round SEM 
CTdf$average_RQ <- round(rowMeans(CTdf[, grepl("^RQ", names(CTdf))], na.rm = TRUE), digits=2) 
CTdf$SEM <- apply(CTdf[, grepl("^RQ", names(CTdf))], 1, function(row) {
  std.error(na.omit(row))
})

#concentration
extract_concentration <- function(string) {
  result <- regmatches(string, regexpr("^\\d*\\.?\\d+", string))
  if (length(result) == 0) {
    return(NA)
  } else {
    return(result)
  }
}
CTdf$concentration <- unlist(lapply(CTdf$chemical,extract_concentration))

# Function to replace chemical symbols with names
replace_symbols <- function(string) {
  index <- match(chemicalsymbols, chemicalsymbols)
  i =0 
  for(chemicalsym in chemicalsymbols){ 
    i<-i+1
    if(grepl(chemicalsym,string)){ 
      return(as.character(chemicalnames[i]))
    }
  }
  return(string)
}



CTdf$chemicalgroup <- unlist(lapply(CTdf$chemical,replace_symbols))
CTdf$concentration <- as.double(CTdf$concentration)

if(analysistype=="western"){
  gene_or_protein_label = "protein"
  names(CTdf)[names(CTdf) == "gene"] <- gene_or_protein_label
  CTdf[[gene_or_protein_label]] <- tolower(CTdf[[gene_or_protein_label]])
  CTdf$graphing_sort <-paste(CTdf$chemicalgroup, CTdf$protein, sep = "_")
}else{ 
  gene_or_protein_label = "gene"
  CTdf[[gene_or_protein_label]] <- tolower(CTdf[[gene_or_protein_label]])
  CTdf$graphing_sort <-paste(CTdf$chemicalgroup, CTdf$gene, sep = "_")
}



CTdf$chemical <- paste0(format(CTdf$concentration, nsmall = 1),unit," ",as.character(CTdf$chemicalgroup)," ",CTdf[[gene_or_protein_label]])
CTdf$concentration <- as.character(CTdf$concentration)

unique_chemicalgroups <- unique(CTdf$chemicalgroup)
unique_genes <- unique(CTdf$gene)

#move concentration and chemical columsn to the front of dataframe 
columns_to_move <- c("chemical", "concentration","chemicalgroup")
column_indices <- which(names(CTdf) %in% columns_to_move)
CTdf <- CTdf[, c(column_indices, setdiff(seq_along(CTdf), column_indices))]



write.table(CTdf,paste0(output_dir,"resall.txt"),sep="\t", row.names=FALSE)





for(chemicals in unique (CTdf$chemicalgroup)) {
  for(genes in unique(CTdf[[gene_or_protein_label]])) {
    df2 <- CTdf[CTdf$chemicalgroup==chemicals & CTdf[[gene_or_protein_label]]==genes, ]
    if(length(unique(df2$concentration))>1){ 
      #print(chemicals)
      #print(genes)
      #df2 <- rbind(df2, CTdf[CTdf$chemicalgroup=="mir+vehcon" & CTdf$protein==genes, ])
      df2[is.na(df2$concentration), "concentration"] <- 0
      
      #only include concentrations with more than 2 reps for statistical test
      rq_columns <- names(df2)[startsWith(names(df2), "RQ")]
      filled_values <- rowSums(!is.na(df2[, rq_columns]))
      more_than_2_filled <- filled_values > 1
      df2 <- df2[more_than_2_filled, ]
      
      if(nrow(df2)>0){ 
        table_for_anova<- df2 %>% 
          pivot_longer(cols = starts_with("RQ"), names_to = "RQ_number", values_to = "RQs")
        
        
        #add control row to compare, because it is fold changes it will always be 1, probably need to properly add SEM
        new_row <- table_for_anova[1,]
        new_row$concentration <- 0 
        new_row$SEM <- 0 
        new_row$RQs <-1
        new_row$chemical <- "control"
        table_for_anova <- rbind(table_for_anova, new_row)
        table_for_anova <- rbind(table_for_anova, new_row)
        table_for_anova <- rbind(table_for_anova, new_row)
        table_for_anova <- table_for_anova[complete.cases(table_for_anova[, "RQs"]), ]
        
        
        anova <- aov(RQs~factor(concentration), data=table_for_anova)
        summary <- summary(glht(anova, linfct = mcp('factor(concentration)' = "Dunnett")))
        p_values <- summary$test$pvalues
        tukey<-TukeyHSD(anova)
    
        sigLetters <- multcompLetters4(anova,tukey)
        sigLetters<-as.data.frame.list(sigLetters$`factor(concentration)`)
        sigLetters$chemicalgroup <- chemicals
        sigLetters[[gene_or_protein_label]] <- genes 
        sigLetters <- rownames_to_column(sigLetters, var ="concentration")
        sigLetters <- sigLetters[, !(names(sigLetters) == "a")]
        sigLetters <- sigLetters[, !grepl("LetterMatrix", names(sigLetters))]
        sigLetters <- sigLetters[sigLetters$concentration != 0, ]
        if (nrow(sigLetters)== length(p_values)){ 
          sigLetters$pvalue <- p_values
        }
        
        CTdf <- left_join(CTdf, sigLetters, by = c("chemicalgroup", gene_or_protein_label, "concentration"))
        columns_to_coalesce <- grep("Letters", names(CTdf), value = TRUE)
        
        # Coalesce selected columns
        CTdf <- CTdf %>% 
          mutate(L = coalesce(!!!syms(columns_to_coalesce)))
        
        # Remove the original columns if needed
        CTdf <- CTdf[, !names(CTdf) %in% columns_to_coalesce]
        names(CTdf)[names(CTdf) == "L"] <- "Letters"
        #print(CTdf)
        
        columns_to_coalesce <- grep("pvalue", names(CTdf), value = TRUE)
        
        # Coalesce selected columns
        CTdf <- CTdf %>% 
          mutate(p = coalesce(!!!syms(columns_to_coalesce)))
        
        # Remove the original columns if needed
        CTdf <- CTdf[, !names(CTdf) %in% columns_to_coalesce]
        names(CTdf)[names(CTdf) == "p"] <- "pvalue"
        }
    }
      }
    
} 

CTdf$concentration <- as.double(CTdf$concentration)
CTdf$concentration[is.na(CTdf$concentration)] <- 0

CTdf$chemical <- paste0(format(CTdf$concentration, nsmall = 1),unit," ",as.character(CTdf$chemicalgroup)," ",CTdf[[gene_or_protein_label]])
CTdf$concentration <- as.character(CTdf$concentration)

CTdf <- CTdf[order(CTdf$chemical), ]
CTdf$chemical <- factor(CTdf$chemical, levels = unique(CTdf$chemical))

rownames(CTdf) <- NULL
