library(dplyr)
library(tm)
library(ggplot2)
library(forcats)
library(plotrix)
path <- "~/storage/qPCR Analysis/Phenanthrene qPCR"
output_dir = paste0(path,"/output/")
unlink(output_dir, recursive = TRUE)
reps <- list.dirs(path,recursive = FALSE)
rep<-0
CTdf <- data.frame()

for (reppath in reps){ 
  print(reppath)
  rep <- rep+1
  plates <- list.dirs(reppath,recursive = FALSE) 
  plate <-0
  for (platepath in plates) { 
    files <- list.files(platepath, pattern=NULL, all.files=FALSE, 
                        full.names=FALSE)
    plate<-plate+1
    df <- data.frame(value = character(),well = character())
    colnames(df) <- c("value", "well")
    for (file in files) { 
      dftemp = data_frame()
      table = read.csv(paste0(platepath,"/",file),sep=",",header=FALSE,fill = TRUE)
      colnames(table) <- table[1,]
      for(row in 2:nrow(table)){ 
        for (col in 2:ncol(table)) { 
          dftemp = rbind(dftemp, c(table[row,col], paste0(table[row,1],table[1,col])))
        }
      }
      if (grepl("-CTmap", file, ignore.case = TRUE)){ 
        colnames(dftemp) <- c("CT value", "well")
      } else if (grepl("-wellmap", file, ignore.case = TRUE)) { 
        colnames(dftemp) <- c("chemical", "well")
      } else if (grepl("-genemap", file, ignore.case = TRUE)) { 
        colnames(dftemp) <- c("gene", "well")
      }
      
      if (nrow(df) ==0){ 
        df = rbind(df,dftemp)
        colnames(df) <- colnames(dftemp)
      } else { 
        df <- left_join(df,dftemp, by="well")
      }
    }
    df <- df[complete.cases(df), ]
    df$`CT value`<- as.numeric(df$`CT value`)
    average_df <- df %>%
      group_by(gene, chemical) %>%
      summarize(average_CT = mean(`CT value`, na.rm = TRUE))
    bactin <- filter(average_df, gene == 'b-actin')
    gene_of_interest <- filter(average_df, gene != 'b-actin')
    joined <- inner_join(bactin, gene_of_interest, by ='chemical')
    joined
    joined$DeltaCT <- joined$average_CT.y-joined$average_CT.x
    control_to_check <- "MIR + veh con "
    joined$chemicalgroup <- gsub(control_to_check, "", removeNumbers(joined$chemical))
    joined$chemicalgroup <- gsub("\\.", "", joined$chemicalgroup)
    joined$chemicalgroup <- trimws(joined$chemicalgroup)
    joined <- joined %>%
      mutate(control = ifelse(chemical ==control_to_check, "C", ""))
    joined <- joined %>% 
      group_by(gene.y) %>%
      mutate(DDCT = DeltaCT-DeltaCT[control=="C"]) 
    joined$RQ <- round(2^(-joined$DDCT), digits=2) 
    joined$rep <- rep
    joined$plate <- plate
    CTdf <- rbind(CTdf,joined %>% select("gene" = gene.y,RQ,chemical, chemicalgroup,rep,plate))
  }
  
}
dir.create(output_dir)
rep_avg <- CTdf %>%
  group_by(chemical, gene, chemicalgroup) %>%
  summarize(
    average_RQ = round(mean(RQ, na.rm = TRUE), digits = 2),
    SEM = std.error(RQ, na.rm = TRUE),
  )
rep_avg <- rep_avg[is.na(rep_avg$SEM)==FALSE,]
unique_chemicalgroups <- unique(rep_avg$chemicalgroup)
unique_genes <- unique(rep_avg$gene)

for (chemgroup in unique_chemicalgroups){ 
  for (g in unique_genes) { 
    subset_table <- rep_avg %>%
      filter(chemicalgroup == chemgroup, gene == g)
    subset_table$concentration <- as.numeric(gsub("[^0-9]", "", subset_table$chemical))
    subset_table <- subset_table[order(subset_table$concentration), ]
    print(subset_table)
      #write.table(subset_table,paste0(output_dir,chemgroup,"-",g,".txt"),sep="\t", row.names=FALSE)
      plot <- ggplot(subset_table, aes(x = fct_inorder(chemical), y = average_RQ)) +
        geom_bar(stat = "identity", fill = "grey", color = "black", width = 0.7) +
        geom_errorbar(aes(ymin = average_RQ - SEM, ymax = average_RQ + SEM), width = 0.25, position = position_dodge(0.7)) +
        labs(title =g, x = "Chemical Concentration", y = "Fold Change") + 
        theme_bw()
      ggsave(paste0(output_dir,chemgroup,"-",g,".png"), plot = plot, width = 6, height = 4, units = "in", dpi = 300)
    }
  }


write.table(CTdf,paste0(output_dir,"resall.txt"),sep="\t", row.names=FALSE)
write.table(rep_avg,paste0(output_dir,"ressummary.txt"),sep="\t", row.names=FALSE)


