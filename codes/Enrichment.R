library(enrichR)
library(readxl)
library(writexl)
library(GOplot)
library(stringr)  
# For error EnrichR is not responding 
#install.packages("devtools") 



setwd("D:/Projects/3-My_Projects/3-Students/mohsen bagheri/")
C <- as.data.frame(read_excel("Networks/1.5_DEGs_all_DEGs_Subnetwork1/Net_info.xlsx",sheet=1))

list_of_genes <- as.data.frame(C$`Differentialy Expressed Genes`[!is.na(C$`Differentialy Expressed Genes`)])
colnames(list_of_genes) <- 'DEGs'


listEnrichrSites()
setEnrichrSite("Enrichr") #

websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

dbs <- c("KEGG_2021_Human","GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015")



DEG_Enr <- enrichr(list_of_genes$DEGs, dbs)

plotEnrich(DEG_Enr[[1]] , showTerms = 20, y="Caount")
write_xlsx(DEG_Enr , "Final_Results/C1_Enrichment.xlsx")




for (Cluster in seq(1,10)) {
  print(Cluster)
  C <- as.data.frame(read_excel("Networks/1.5_DEGs_all_DEGs_Subnetwork1/Net_info.xlsx",sheet=Cluster))
  
  list_of_genes <- as.data.frame(C$`Differentialy Expressed Genes`[!is.na(C$`Differentialy Expressed Genes`)])
  colnames(list_of_genes) <- 'DEGs'
  
  DEG_Enr <- enrichr(list_of_genes$DEGs, dbs)
  
  pdf(file = file.path('Final_Results' , paste('enrichment_', Cluster, '.pdf', sep = '')))
  plotEnrich(DEG_Enr[[1]] , showTerms = 20, y="Caount")
  dev.off()
  write_xlsx(DEG_Enr ,  file.path('Final_Results' , paste('enrichment_', Cluster, '.xlsx',sep = '')))
  
  
  
  
}









#### This section is not needed.... 
#_______________________________________________________________________________
# NS
NS <- read_excel("N_S/Gene_list/Node_loist_clustered.xlsx",sheet=1)

if (websiteLive) {
  NS_Enr_C1 <- enrichr(NS[NS$Clustering ==1 , ]$ORF, dbs)
  NS_Enr_C2 <- enrichr(NS[NS$Clustering ==2 , ]$ORF, dbs)
  NS_Enr_C3 <- enrichr(NS[NS$Clustering ==3 , ]$ORF, dbs)

}
plotEnrich(NS_Enr_C1[[1]] , showTerms = 10 , y="Caount")
plotEnrich(NS_Enr_C2[[1]] , showTerms = 10 , y="Caount")
plotEnrich(NS_Enr_C3[[1]] , showTerms = 10 , y="Caount")
write_xlsx(NS_Enr_C1 , "N_S/Enrich_Results_C1.xlsx")
write_xlsx(NS_Enr_C2 , "N_S/Enrich_Results_C2.xlsx")
write_xlsx(NS_Enr_C3 , "N_S/Enrich_Results_C3.xlsx")

#_______________________________________________________________________________
#NM
NM <- as.data.frame(read_excel("N_M/Gene_list/Node_list_clustered.xlsx",sheet=1))
L <- list()
n <- 1
Clusters <- c(2,4,5,9,11,12,13)

NM_Enr_C2 <- enrichr(NM[NM$Cluster ==2 , ]$ORF, dbs)
NM_Enr_C4 <- enrichr(NM[NM$Cluster ==4 , ]$ORF, dbs)
NM_Enr_C5 <- enrichr(NM[NM$Cluster ==5 , ]$ORF, dbs)
NM_Enr_C9 <- enrichr(NM[NM$Cluster ==9 , ]$ORF, dbs)
NM_Enr_C11 <- enrichr(NM[NM$Cluster ==11 , ]$ORF, dbs)
NM_Enr_C12 <- enrichr(NM[NM$Cluster ==12 , ]$ORF, dbs)
NM_Enr_C13 <- enrichr(NM[NM$Cluster ==13 , ]$ORF, dbs)
write_xlsx(NM_Enr_C2 , "N_M/Enrich_Results_C2.xlsx")
write_xlsx(NM_Enr_C4 , "N_M/Enrich_Results_C4.xlsx")
write_xlsx(NM_Enr_C5 , "N_M/Enrich_Results_C5.xlsx")
write_xlsx(NM_Enr_C9 , "N_M/Enrich_Results_C9.xlsx")
write_xlsx(NM_Enr_C11 , "N_M/Enrich_Results_C11.xlsx")
write_xlsx(NM_Enr_C12 , "N_M/Enrich_Results_C12.xlsx")
write_xlsx(NM_Enr_C13 , "N_M/Enrich_Results_C13.xlsx")

#_______________________________________________________________________________
# MS
MS <- as.data.frame(read_excel("M_S//Gene_list/Node_list_clustered.xlsx",sheet=1))
MS_Enr_C1 <- enrichr(MS$ORF, dbs)
write_xlsx(MS_Enr_C1 , "M_S/Enrich_Results_C1.xlsx")





if (websiteLive) plotEnrich(enriched[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")

mod = mod4
colnames(mod) = c("logFC","ID","degree","sum","sum2")
mod_enriched = mod4_enriched

genes= data.frame(ID=mod[,2],
                  logFC=mod[,1],
                  adj.P.Val = rep(0.002,times=39),
                  sum=mod[,4],
                  sum2=mod[,5],
                  degree=mod[,3])
modcirclist = list(
  kegg= data.frame(Category= rep("a",times=102),
                   Term=mod_enriched[["KEGG_2021_Human"]][["Term"]],
                   adj_pval= mod_enriched[["KEGG_2021_Human"]][["Adjusted.P.value"]],
                   Genes= str_replace_all( mod_enriched[["KEGG_2021_Human"]][["Genes"]],";",", ")),
  genes= genes)


modcirc <- circle_dat(modcirclist[["kegg"]], modcirclist[["genes"]])

genes = genes[genes$sum>7,]
genes = genes[abs(genes$logFC)>0,]
genes = genes[order(abs(genes$logFC), decreasing = TRUE),]
                          
chord <- chord_dat(data = modcirc, genes = genes[,1:2],process =colnames(n))
chord <- chord_dat(data = modcirc, genes = genes[,1:2])
pdf("m4_chord.pdf",         # File name
    width = 22, height = 20) 
GOChord(chord,space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5,limit = c(2, 2))
dev.off()

write_xlsx(mod0_enriched,
           "mod0EnrichR.xlsx")

write_xlsx(mod1_enriched,
           "mod1EnrichR.xlsx")
write_xlsx(mod2_enriched,
           "mod2EnrichR.xlsx")
write_xlsx(mod3_enriched,
           "mod3EnrichR.xlsx")
write_xlsx(mod4_enriched,
           "mod4EnrichR.xlsx")
write_xlsx(mod5_enriched,
           "mod5EnrichR.xlsx")
