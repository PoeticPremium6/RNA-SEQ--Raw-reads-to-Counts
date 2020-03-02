Tuxedo_Jonathan <- read.csv(file.path("/home/jonathan/Documents/Tuxedo_Counts_(Jon).csv"), header=TRUE)


head(Tuxedo_Jonathan)

#Extract RowNames
rownames(Tuxedo_Jonathan) <- Tuxedo_Jonathan$gene.AF8.id

#Create a function for quantile normalization
quantile_normalization <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)

  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}
    
Normalized_Tuxedo_Jonathan <- quantile_normalization(Tuxedo_Jonathan[,2:70])
library(xlsx)

#Export data
write.table(Normalized_Tuxedo_Jonathan, "/home/jonathan/Documents/Comp_Micro/iMAT/Normalized_Tuxedo_Jonathan.txt", sep="\t")
 
