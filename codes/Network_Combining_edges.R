library(readxl)


setwd("D:/Projects/3-My_Projects/3-Students/mohsen bagheri/")

# This section is for constructing the first network. 3 gene with 3 methylation

nodes <- read.table("Networks/1.5_DEG/DEG_To_Methabolits/MethaboAnalysis/Download (1)/orig_node_list.csv" , sep = ',' , header = T)
edges <- read.table("Networks/1.5_DEG/DEG_To_Methabolits/MethaboAnalysis/Download (1)/orig_edge_list.csv" , sep = ',' , header = T)



nodes[nodes$Id==edges[,'Source'],]
nodes[nodes$Id==edges[,'Target'],]


temp1 <- edges
temp2 <- nodes

head(temp1)

for (i in seq(length(temp1$Target))) {
  for (j in seq(dim(nodes)[1])) {
    for (j in seq(dim(nodes)[1])) {
      print(i)
      print(nodes$Label[j])
      temp1[i,2] = nodes$Label[j]
    }
  }
}
      
tail(temp1)  

t <- temp[,c("Gene Symbol" , "annot.id")] # Merged DEGs with methylation
colnames(t) <- colnames(temp1)

tt <- rbind(temp1 , t)
tail(tt)
write.table(tt , "Networks/1.5_DEG/DEG_To_Methabolits/MethaboAnalysis/Download (1)/final_edge_list.txt" ,quote = F,row.names = F,)





# this section is for constructing the second network, subnetwork 1 with methylations and methabol
Subnet_Node <- as.data.frame(read_excel("Networks/1.5_DEGs_all_DEGs_Subnetwork1/Subnet1_Nodes.xlsx"))
Subnet_edges <- as.data.frame(read_excel("Networks/1.5_DEGs_all_DEGs_Subnetwork1/Subnet1_Edges.xlsx"))

head(Subnet_edges)
head(Subnet_Node)

temp_nodes <- Subnet_Node[,c("name" , "display name")]
colnames(temp_nodes) <- c("Source" , "display name")
temp_out <-merge(temp_nodes,Subnet_edges , by = "Source")
dim(temp_out) 
dim(Subnet_edges)
temp_out
colnames(temp_out) <- c("Source_id","Source","EdgeBetweenness","Target")
Subnet_edges <- temp_out


colnames(temp_nodes) <- c("Target" , "display name")
temp_out <-merge(temp_nodes,Subnet_edges , by = "Target")
dim(temp_out) 
dim(Subnet_edges)
temp_out
colnames(temp_out) <- c("Target_id","Target","Source_id","Source" , "EdgeBetweenness")
Subnet_edges <- temp_out
Subnet_edges




metha_edges <- as.data.frame(read_excel("Networks/1.5_DEGs_all_DEGs_Subnetwork1/1.5_DEGs,AllEnriched_Subnetwork1_withmethabolits/Edges_File1.xlsx" ,))


t <- temp_merged_1.5DEGs_Methy[,c("Gene Symbol" , "annot.id")]
t
colnames(t) <- c("Source","Target")
final_edge <- rbind(Subnet_edges[,c("Source","Target")],t, metha_edges)
dim(final_edge)


write.table(final_edge, "Networks/1.5_DEGs_all_DEGs_Subnetwork1/final_edge_list.txt" ,quote = F,row.names = F)









temp2 <- unique(final_edge$Source)
s <- ''
for (i in temp2) {
    s <- paste(s , i , sep = ' OR ')
}
s


t

temp2 <- unique(t$Target)
s <- ''
for (i in temp2) {
  s <- paste(s , i , sep = ' OR ')
}
s
