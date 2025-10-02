library(readxl)



setwd("D:/Projects/3-My_Projects/3-Students/mohsen bagheri/")
# Reading methilation and DEG excel files

Methilation <- read_excel("data_processed/Diff_(un)annotated_methlome_hypo_fold_processed.xlsx",sheet = 'Sheet2')
class(Methilation)
Methilation <- as.data.frame(Methilation)
head(Methilation)


DEGs_1.5 <- as.data.frame(read_excel("data_processed/Final_DEGsSC_Processed.xlsx" , sheet = c(2)))
for (i in c(3,4,5)) {
  DEGs_1.5 <- rbind(DEGs_1.5 , as.data.frame(read_excel("data_processed/Final_DEGsSC_Processed.xlsx" , sheet = c(i))))
}

head(DEGs_1.5)
dim(DEGs_1.5)



colnames(Methilation) -> temp
temp[82] <- "Gene Symbol"
temp[82]
colnames(Methilation) <- temp

temp_merged_1.5DEGs_Methy <-merge(Methilation , DEGs_1.5 , by = "Gene Symbol")
temp_merged_1.5DEGs_Methy$`Gene Symbol`









DEGs_1 <- as.data.frame(read_excel("data_processed/Final_DEGsSC_Processed.xlsx" , sheet = c(1)))
DEGs_1 <- DEGs_1[abs(DEGs_1$logFC) > 1,]
DEGs_1 <- DEGs_1[DEGs_1$P.Value < 0.05,]
dim(DEGs_1)
head(DEGs_1)
colnames(DEGs_1) -> temp
temp[1] <- "Gene Symbol"
colnames(DEGs_1) <- temp


temp_merged_1DEGs_Methy <-merge(Methilation , DEGs_1 , by = "Gene Symbol")
temp_merged_1DEGs_Methy$`Gene Symbol`
unique((temp_merged_1DEGs_Methy$`Gene Symbol`))
dim(DEGs_1)






