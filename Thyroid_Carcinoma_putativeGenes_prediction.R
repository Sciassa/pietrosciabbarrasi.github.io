if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks, quietly = T)
library(SummarizedExperiment,quietly = T)
library(dplyr, quietly = T)
library(ggplot2,quietly = T)
BiocManager::install("DESeq2")
library(DESeq2, quietly = T)

install.packages("ggrepel")
library(ggrepel)
# Define project ID
proj <- "TCGA-THCA"
dir.create(file.path(proj), showWarnings = FALSE)

load("my_data.RData")

###########################################################

# Query for Tumor Samples
rna.query.Tumor <- GDCquery(
  project = proj,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor"
)
# Download Tumor Samples
GDCdownload(query = rna.query.Tumor, directory = "GDCdata", method = "api")
rna.data.Tumor <- GDCprepare(rna.query.Tumor, directory = "GDCdata")

# Query for Normal Samples
rna.query.Normal <- GDCquery(
  project = proj,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Solid Tissue Normal"
)
# Download Normal Samples
GDCdownload(query = rna.query.Normal, directory = "GDCdata", method = "api")
rna.data.Normal <- GDCprepare(rna.query.Normal, directory = "GDCdata")
###########################################################

# represent the RNA-seq data for normal and tumor samples,
head(rna.data.Tumor)  
head(rna.data.Normal)


#  extracting the expression matrices 
# row: genes ; col: samples
rna.expr.data.Tumor <- assay(rna.data.Tumor)
rna.expr.data.Normal <- assay(rna.data.Normal)


# converting the genomic information from both the tumor and normal RNA-seq datasets into data frames.
genes.info <- BiocGenerics::as.data.frame(rowRanges(rna.data.Tumor))
genes.info2 <- BiocGenerics::as.data.frame(rowRanges(rna.data.Normal))
head(genes.info,1)
all(na.omit(genes.info2) == na.omit(genes.info))

# Data Cleaning
ncol(rna.expr.data.Normal) # 59
head(colnames(rna.expr.data.Normal))
head(substr(colnames(rna.expr.data.Normal), 1,12))

dim(rna.expr.data.Tumor) # 60660 x 505
length(unique(substr(colnames(rna.expr.data.Tumor), 1,12))) #no duplicates - > 505
dim(rna.expr.data.Tumor)

patients.Tumor <- substr(colnames(rna.expr.data.Tumor), 1,12) # Get first 12 characters that uniquely identify patient
## sort(table(patients.Tumor))
unique.patients.Tumor <- names(which(table(patients.Tumor) == 1)) # Get unique patients (we don't have duplicates)
# index of patients only keeping patients with one sample
idx.unique.pats <- match(unique.patients.Tumor, substr(colnames(rna.expr.data.Tumor), 1,12) )
expr.Tumor <- as.data.frame(rna.expr.data.Tumor[,idx.unique.pats])
expr.Normal <- as.data.frame(rna.expr.data.Normal)

#let's rename patients in a shorter way
colnames(expr.Tumor) <- substr(colnames(expr.Tumor), 1,12)
unique(colnames(expr.Tumor))

colnames(expr.Normal) <- substr(colnames(expr.Normal), 1,12)
unique(colnames(expr.Normal))
intersect(colnames(expr.Normal), colnames(expr.Tumor)) # still 59
setdiff(colnames(expr.Normal), colnames(expr.Tumor))
match(setdiff(colnames(expr.Normal), colnames(expr.Tumor)), colnames(expr.Normal)) #idx to remove (we don't need to remove any)
length(intersect(colnames(expr.Normal), colnames(expr.Tumor))) # 59

#let's check the actual counts
typeof(expr.Tumor[1,1]) # integer -> ok
any(is.na(expr.Tumor)) # False -> ok
any(is.nan(as.matrix(expr.Tumor))) # False -> ok

typeof(expr.Normal[1,1]) # integer -> ok
any(is.na(expr.Normal)) # False -> ok
any(is.nan(as.matrix(expr.Normal))) # False -> ok

#let's consider only patients for which we have both normal and cancer samples
expr.Tumor <- expr.Tumor[, colnames(expr.Normal)]



#### Normalizing With DESeq2
# Concatenate normal and tumor into a matrix for normalization
all(rownames(expr.Tumor) == rownames(expr.Normal)) # TRUE -> ok
full.data <- cbind(expr.Normal, expr.Tumor)

dim(full.data) # 60660 x 118
full.data <- data.frame(full.data)

metad <- rep("cancer", 118)
metad[1:59] <- "normal"
metad
metad <- data.frame(metad)
rownames(metad) <- colnames(full.data)
colnames(metad)[1] <- "condition"
metad[,1] <- as.factor(metad[,1]) # Converts the column to a factor data type which is a categorical variable with a fixed number of unique levels.
full.data <- cbind(rownames(full.data), full.data)

# Normalization
dds <- DESeqDataSetFromMatrix(countData=full.data,colData=metad,design= ~condition,tidy=TRUE)

View(counts(dds))
dim(counts(dds)) # 60660 x 118

# filtering: at least ten counts on 90% of patients?
( 118*90 )/100
keep <- rowSums(counts(dds) >= 10) >= 106 #ensures that only genes with at least 10 counts in 90% of the samples (i.e., 106 out of 118 samples) are retained.

dds <- dds[keep,]
dim(counts(dds)) # 17001 x 118

dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
sum(rowSums(normalized_counts == 0) == 118) #no null rows

filtr.expr.n <- as.data.frame(normalized_counts[, 1:59])
filtr.expr.c <- as.data.frame(normalized_counts[, 60:118])
colnames(filtr.expr.c) <- substr(colnames(filtr.expr.c), 1,12)
colnames(filtr.expr.n) <- substr(colnames(filtr.expr.n), 1,12)

colnames(filtr.expr.c) <-gsub("\\.", "-", colnames(filtr.expr.c))
colnames(filtr.expr.n) <-gsub("\\.", "-", colnames(filtr.expr.n))

###############

# Calculate Fold Change
fc <- log2(rowMeans(filtr.expr.c) / rowMeans(filtr.expr.n))
names(fc) <- rownames(filtr.expr.c)
head(fc)
# ENSG00000000003.15 ENSG00000000419.13 ENSG00000000457.14 ENSG00000000460.17 ENSG00000000938.13 ENSG00000000971.16 
# 0.11627942        -0.05581529        -0.16456876         0.30007278         0.60682520         0.33657513 

# Perform Paired T-tests
pval.fc <- sapply(1:nrow(filtr.expr.c), function(i) (t.test(as.numeric(filtr.expr.c[i,]), as.numeric(filtr.expr.n[i,]), paired = T ))$p.value)

# Adjust P-values for Multiple Comparisons (FDR)
pval.fc.fdr <- p.adjust(pval.fc, method = "fdr")

# Create DEG Table
expr.table <- data.frame(fc = fc, pval.fc.fdr = pval.fc.fdr)
expr.table <- na.omit(expr.table)
deg.genes <- rownames(expr.table[abs(expr.table$fc) >= 2 & expr.table$pval.fc.fdr <= 0.01, ])
length(deg.genes) # 430

deg.table=expr.table[deg.genes,]
deg.table=deg.table[order(-abs(deg.table$fc)),] 




################ VOLCANO PLOT ! https://biostatsquid.com/volcano-plots-r-tutorial/


deg.table=expr.table[deg.genes,]
deg.table=deg.table[order(-abs(deg.table$fc)),]
top_degs=genes.info[rownames(deg.table),"gene_name"]
top3_degs=genes.info[c("ENSG00000145864.14","ENSG00000174460.4","ENSG00000137648.19"),"gene_name"]
bottom3_degs=genes.info[c("ENSG00000183775.11","ENSG00000160180.15","ENSG00000138722.10"),"gene_name"]

expr.table$gene_symbol=genes.info[rownames(expr.table),"gene_name"]
expr.table$delabel <- ifelse(expr.table$gene_symbol %in% top3_degs | expr.table$gene_symbol %in% bottom3_degs, expr.table$gene_symbol, NA)
expr.table$diffexpressed <- "NO"
expr.table$diffexpressed[expr.table$fc >= 2 & expr.table$pval.fc.fdr <= 0.01] <- "UP"
expr.table$diffexpressed[expr.table$fc <= -2 & expr.table$pval.fc.fdr <= 0.01] <- "DOWN"



# top3 up regulated e down 3 regulated fare dopo tabella
volcanoplot =
  ggplot(data = expr.table, aes(x = fc, y = -log10(pval.fc.fdr), col = diffexpressed, label = delabel)) +
  geom_hline(yintercept = -log10(0.01), col = "white", linetype = 'dashed') + 
  geom_vline(xintercept = c(-2, 2), col = "white", linetype = 'dashed') +
  geom_point(size = 0.8) + 
  scale_color_manual(values = c("green", "grey", "orange"), # to set the colours of our variable  
                     labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(-1, 23), xlim = c(-5, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis
  geom_label_repel(max.overlaps = Inf, key_glyph = draw_key_point) + # To show all labels 
  theme_dark() +
  theme(legend.position = "none") # Remove the legend

volcanoplot # Execute the plot


###########
install.packages("pacman")
pacman::p_load(psych, network, sna, NetworkToolbox, ggplot2, GGally) # efficient way to insall and load packages

# Spearman correlation matrix for cancer samples (only DEGs)
cor_mat_c = corr.test(t(log2(filtr.expr.c[deg.genes,] + 1)), use = "pairwise", method = "spearman", adjust = "fdr", ci = FALSE)
# Spearman correlation matrix for normal samples (only DEGs)
cor_mat_n = corr.test(t(log2(filtr.expr.n[deg.genes,] + 1) ), use = "pairwise", method = "spearman", adjust = "fdr", ci = FALSE)

rho_c=cor_mat_c$r # correlation matrix
rho_n=cor_mat_n$r

diag(rho_c) <- 0  # Set diagonal to 0 (no self-correlation)
diag(rho_n) <- 0

# Define threshold for adjacency matrix
threshold <- 0.6
adj_mat_c = (abs(rho_c)>=threshold)*1
adj_mat_n = (abs(rho_n)>=threshold)*1

net_c=network(adj_mat_c,matrix.type="adjacency",directed = F,ignore.eval = FALSE, names.eval = "weights")
round(network.density(net_c),2) # 0.11
round(network.size(net_c),2) # 430 
round(network.edgecount(net_c),2) # 9708
round(clustcoeff(adj_mat_c,weighted = F)$CC,2) # 0.51

degree_c=rowSums(adj_mat_c!=0)
degree_c=sort(degree_c,decreasing = T)

degree_c=degree_c[degree_c>0]
length(degree_c)
x_c=quantile(degree_c,0.95)

net_n=network(adj_mat_n,matrix.type="adjacency",directed = F,ignore.eval = FALSE, names.eval = "weights")
round(network.density(net_n),2) # 0.02
round(network.size(net_n),2) # 430
round(network.edgecount(net_n),2) # 2086
round(clustcoeff(adj_mat_n,weighted = F)$CC,2) # 0.28


degree_n=rowSums(adj_mat_n!=0)
degree_n=sort(degree_n,decreasing = T)

degree_n=degree_n[degree_n>0]
length(degree_n) #285
x_n=quantile(degree_n,0.95)

par(mfrow=c(1,2))
hist(degree_n,col=NA,border="blue",lwd=2,breaks=20,main="",xlab="Degree Normal Co-Expresson")
hist(degree_c,col=NA,border="red",lwd=2,breaks=20,main="",xlab="Degree Cancer Co-Expresson")

hubs_n <- degree_n[degree_n >= x_n]
genes.info[names(hubs_n),"gene_name"]
# TMC6 CYP2S1 RUNX1 TRIM46 DUSP4 MYO1G KCNN4 APOE\\ PLEKHN1 BID SPOCK2 ALDH3B1 ALOX5 TIMP1 GPRIN1 CEROX1
hubs_c
hubs_c <- degree_c[degree_c >= x_c]
genes.info[names(hubs_c),"gene_name"]
genes.info["ENSG00000134042.14","gene_name"] 
cbind(genes.info[names(hubs_c),"gene_name"], as.numeric(hubs_c))
# SLC4A4 MRO CYP2S1 TMPRSS4 KRT19 SERPINA1 MPPED2 NECTIN4 LAMB3 GRB7 ALDH3B1 FN1 KCNN4 PTPRE PDLIM4 STAC LGALS3 RUNX1 LIPH EPHB3

intersect(genes.info[names(hubs_n),"gene_name"],genes.info[names(hubs_c),"gene_name"]) # "CYP2S1"  "RUNX1"   "KCNN4"   "ALDH3B1"

########## Compute Centrality Index (CI) and Check Overlap #############
bet_c=sna::betweenness(net_c)
top_bet_c=bet_c[which(bet_c>=quantile(bet_c, 0.95))]
top_bet_c_names=genes.info[network.vertex.names(net_c)[which(bet_c >= quantile(bet_c, 0.95))],"gene_name"]
d_bet_c=as.data.frame(cbind(top_bet_c_names,top_bet_c))
intersect(top_bet_c_names,genes.info[names(hubs_c),"gene_name"]) # "ALDH3B1" "MPPED2"  "SLC4A4"  "MRO" 

bet_n=sna::betweenness(net_n)
top_bet_n=bet_n[which(bet_n>=quantile(bet_n, 0.95))]
top_bet_n_names=genes.info[network.vertex.names(net_n)[which(bet_n >= quantile(bet_n, 0.95))],"gene_name"]
d_bet_n=as.data.frame(cbind(top_bet_n_names,top_bet_n))
intersect(top_bet_n_names,genes.info[names(hubs_n),"gene_name"]) # "ALDH3B1" "BID"     "TIMP1"   "SPOCK2"  "DUSP4"   "TMC6"    "TRIM46"  "PLEKHN1" "CEROX1" 



# enrichment
study=c("FN1","TFF3","ALDH3B1","SLC4A4")

databases <- c("GO_Biological_Process_2023", "KEGG_2021_Human")

# GSEA for cancer hubs
enriched_s <- enrichr(study, databases)
print(enriched_s[["GO_Biological_Process_2023"]])

top_enriched_s <- enriched_s[["GO_Biological_Process_2023"]][1:10, ]
ggplot(top_enriched_s, aes(x = reorder(Term, -log10(P.value)), y = -log10(P.value))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Enriched GO Terms for Cancer Hubs",
       x = "GO Term", y = "-log10(P-value)") +
  theme_minimal()
cbind(top_enriched_s$Genes[c(1,2,4,10)],cbind(top_enriched_s$Term,top_enriched_s$Genes)[c(1,2,4,10)],round(-log10(top_enriched_s$P.value)[c(1,2,4,10)],3))

#######

# Fisher Z-transformation
z_transform <- function(rho) {
  return(0.5 * log((1 + rho) / (1 - rho)))
}

# Apply Z-transformation
z_cancer <- z_transform(rho_c)
z_normal <- z_transform(rho_n)
length(z_normal)


# Calculate Z-difference using Z-score formula
n_cancer <- ncol(filtr.expr.c) # Number of cancer samples
n_normal <- ncol(filtr.expr.n) # Number of normal samples


z_diff <- (z_cancer - z_normal) / sqrt(1 / (n_cancer - 3) + 1 / (n_normal - 3))

# Binary adjacency matrix |Z| < 3
adj_diff <- ifelse(abs(z_diff) >= 3, 1, 0)
diag(adj_diff) <- 0 # Ensure no self-loops

# Create the differential network
net_diff <- network(adj_diff,matrix.type="adjacency",ignore.eval = FALSE,names.eval = "weights", directed = F)
network.density(net_diff) # 0.1222529
network.size(net_diff) # 430 
network.edgecount(net_diff) #11276
round(clustcoeff(adj_diff,weighted = F)$CC,2) #0.26

degree_diff=rowSums(adj_diff!=0)
degree_diff=sort(degree_diff,decreasing = T)

degree_diff=degree_diff[degree_diff>0]
length(degree_diff) #429

par(mfrow=c(1,3))
hist(degree_n,col=NA,border="blue",lwd=2,breaks=20,main="",xlab="Degree Normal Co-Expresson Network")
hist(degree_c,col=NA,border="red",lwd=2,breaks=20,main="",xlab="Degree Cancer Co-Expresson Network")

hist(degree_diff,col=NA,border="darkviolet",lwd=2,breaks=30,main="",xlab="Degree Differential Co-Expresson Network")

# Identify hubs (top 5%)
threshold_hub <- quantile(degree_diff, 0.95)
hubs_diff <- names(degree_diff[degree_diff >= threshold_hub])
genes.info[hubs_diff,"gene_name"]
#MEX3A MET KCP RSPO4 FOXQ1 AHNAK2 MXRA8 EVA1A PDE5A FAXC GRB7 AC119427.1 TMEM184A GALNT7 PLPP2 TFF3 CD55 CAMK2N1 ITGA2 SDC4 DOCK9-DT IGF2BP2

# Compare with Task 3 Hubs
common_hubs_normal <- intersect(hubs_diff, names(hubs_n))
genes.info[common_hubs_normal,"gene_name"] # NA

common_hubs_cancer <- intersect(hubs_diff, names(hubs_c))
genes.info[common_hubs_cancer,"gene_name"] # GRB7


#######
bet_diff=sna::betweenness(net_diff)
top_bet_diff=bet_diff[which(bet_diff>=quantile(bet_diff, 0.95))]
top_bet_diff_names=genes.info[network.vertex.names(net_diff)[which(bet_diff >= quantile(bet_diff, 0.95))],"gene_name"]
d_bet_diff=as.data.frame(cbind(top_bet_diff_names,top_bet_diff))
d_bet_diff
intersect(intersect(top_bet_diff_names,genes.info[hubs_diff,"gene_name"]),head(top_degs,100))
intersect(genes.info[hubs_diff,"gene_name"],head(top_degs,100))

######

hubs_diff_ids <- vector("integer",length(hubs_diff))
for (i in 1:length(hubs_diff)){
  hubs_diff_ids[i] <- match(hubs_diff[i],rownames(adj_diff))
}

hubs_diff_neigh <- c()
for (f in hubs_diff_ids){
  hubs_diff_neigh <- append(hubs_diff_neigh, get.neighborhood(net_diff, f))
}

hubs_diff_neigh <- unique(hubs_diff_neigh)
hubs_diff_neigh_names <- rownames(adj_diff[hubs_diff_neigh,])
subnet <- unique(c(names(hubs_diff), hubs_diff_neigh_names))

hub_diff_adj <- adj_diff[subnet, subnet]
dim(hub_diff_adj) # 374 x 374

net_diff_hub <- network(hub_diff_adj, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights")

net_diff_hub %v% "type" = ifelse(network.vertex.names(net_diff_hub) %in% hubs_diff,"hub", "non-hub")
net_diff_hub %v% "color" <- ifelse(network.vertex.names(net_diff_hub) %in% hubs_diff,"purple","green") # hubs

node_types_hub_diff <- net_diff_hub %v% "type"
edge_colors_hub_diff <- apply(as.matrix(net_diff_hub, matrix.type = "edgelist"), 1, function(edge) {
  node1 <- as.numeric(edge[1])
  node2 <- as.numeric(edge[2])
  if (node_types_hub_diff[node1] == "hub" || node_types_hub_diff[node2] == "hub") {"violet" } else {"palegreen"} })

network::set.edge.attribute(net_diff_hub, "ecolor", edge_colors_hub_diff)

table(net_diff_hub %v% "color") # green 352 ; purpl  22 

base_plot2=ggnet2(net_diff_hub,color = "color",alpha = 0.9,size = 2,
                  edge.color = "ecolor",edge.alpha = 0.9,edge.size = 0.15,
                  node.label = NA,
) + guides(size = "none")#+ theme(panel.background = element_rect(fill="black"))

base_plot2

########
head(rownames(filtr.expr.c))
head(colnames(filtr.expr.c))
# Compute patient similarity matrix using Spearman correlation
similarity_matrix <- cor(log2(filtr.expr.c[deg.genes, ]+1), method = "spearman")


# Visualize the similarity matrix
diag(similarity_matrix)=0
head(similarity_matrix)

# Threshold the similarity matrix to create a binary adjacency matrix
threshold <- quantile(similarity_matrix, 0.70)
adj_matrix_psn <- ifelse(similarity_matrix >= threshold, 1, 0)
diag(adj_matrix_psn) <- 0  # Remove self-loops

# Create Patient Similarity Network
psn_network <- network(adj_matrix_psn, directed = FALSE)

# Calculate basic network properties
cat("Number of Nodes (Patients):", network.size(psn_network), "\n") # 59
cat("Number of Edges:", network.edgecount(psn_network), "\n") # 523
cat("Network Density:", network.density(psn_network), "\n") # 0.3056 sparse network

degree_psn=sort(rowSums(adj_matrix_psn!=0),decreasing = T)
degree_psn=degree_psn[degree_psn>0]

x_psn=quantile(degree_psn,0.95)
hubs_psn <- degree_psn[degree_psn >= x_psn]
hubs_psn_name <- names(degree_psn[degree_psn >= x_psn])


######
# B]

psn_comp=component.largest(psn_network,result="graph")
psn_comp <- adj_matrix_psn[rownames(psn_comp), rownames(psn_comp)]
head(psn_comp)

write.csv2(psn_comp, "input-matrix.csv")
View(as.data.frame(read.csv2("input-matrix.csv", row.names = 1) ) )

# let's open the terminal 
# pip install bctpy
# python3 btc-community.py input-matrix.csv

psn_comm <- read.csv2("output.txt", header = FALSE)
rownames(psn_comm) <- rownames(psn_comp)
length(table(psn_comm[,1]))
table(psn_comm[,1])

psn_final <- network(psn_comp, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights", directed = F)

all(psn_final %v% "vertex.names" == rownames(psn_comm)) 
psn_final  %v% "community" <-  as.character(psn_comm[,1])

# Map community colors to the network nodes
node_colors <- rainbow(length(unique(psn_comm[, 1])),start=0.1)
names(node_colors) <- unique(psn_comm[, 1])
psn_final %v% "color" <- node_colors[psn_final %v% "community"]

# Assign edge colors based on node attributes
network::set.edge.attribute(psn_final, "edgecolor",
                            ifelse(as.numeric(psn_final %e% "in.community") == 1,"color", "black")  # Edges within communities = "color", others = "grey50"
)


# Find centroids of each community
centroids_psn <- sapply(unique(psn_comm[, 1]), function(community) {
  community_nodes <- which(psn_final %v% "community" == community)
  degrees <- rowSums(as.matrix(psn_final, matrix.type = "adjacency")[community_nodes, ])
  community_nodes[which.max(degrees)]  # Node with the highest degree in the community
})

# Get the names of the centroids
centroid_labels_psn <- network.vertex.names(psn_final)[centroids_psn]

# Visualize the network
# https://briatte.github.io/ggnet/
ggnet2(psn_final,color = "color",size = 4,edge.color = c("color", "black"),  # Use community colors for edges or default grey
       edge.alpha = 0.75)#,label = centroid_labels_psn,  # Label only centroids
#  label.color = "black",label.size = 4) 



#### Survival analysis and Enrichment analysis (OPTIONAL) #################
## IT DIDN'T FIT IN THE REPORT
install.packages("rjson")
install.packages("WriteXLS")
install.packages("ggsurvfit")
install.packages("survminer")
install.packages("enrichR_3.2.tar.gz", repos = NULL, type = "source")
library(enrichR)
library(survival)
library(ggsurvfit)
library(dplyr)
library(survminer)

clinical_data <- read.delim("clinical.tsv", header = TRUE, sep = "\t")

head(clinical_data)
colnames(clinical_data)
table(clinical_data$vital_status) # Alive 982, Dead 32
clinical_data$status <- ifelse(clinical_data$vital_status == "Alive", 0, 1)

psn_nodes <- network.vertex.names(psn_final) # Nodes of psn which are patients IDs
sum(duplicated(clinical_data$case_submitter_id)) # 507 duplicates
clinical_data <- clinical_data[!duplicated(clinical_data$case_submitter_id), ]
patients.clinical_data <- clinical_data[clinical_data$case_submitter_id %in% psn_nodes, ]
patients.clinical_data$community <- psn_final %v% "community" # Add community
length(psn_nodes) == length(psn_final %v% "community") # TRUE-> ok

setdiff(psn_nodes, patients.clinical_data$case_submitter_id)  # IDs in PSN but not in clinical data
setdiff(patients.clinical_data$case_submitter_id, psn_nodes)  # IDs in clinical data but not in PSN
table(patients.clinical_data$community) # 59 patiens




# Overall survival:
# Checks whether the patient's vital_status is "Alive".
# If TRUE (patient is alive), the survival time is taken as the days_to_last_follow_up.
# If FALSE (patient is deceased), the survival time is taken as the days_to_death.
patients.clinical_data$os <- ifelse(patients.clinical_data$vital_status == "Alive", 
                                    patients.clinical_data$days_to_last_follow_up,
                                    patients.clinical_data$days_to_death)
patients.clinical_data$os <- as.numeric(patients.clinical_data$os)

fit <- survfit(Surv(os, status) ~ community, data = patients.clinical_data)
table(patients.clinical_data$community, patients.clinical_data$status)
# With only 2 and 2 deaths in the two communities, 
#the number of observed events is insufficient for reliable survival analysis.

# Plot survival curves for each community
ggsurvplot(
  fit,
  data = patients.clinical_data,
  conf.int = TRUE,             
  # risk.table = TRUE,           
  # surv.median.line = "hv",     
  palette = "Dark2",           
  title = "Kaplan-Meier Survival Curves by PSN Community",
  ggtheme = theme_bw()  
)


# Survival probability plot
survfit2(Surv(os, status) ~ 1, data = patients.clinical_data) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability")

#CI
survfit2(Surv(os, status) ~ 1, data = patients.clinical_data) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability")+ 
  add_confidence_interval()

#risk table
survfit2(Surv(os, status) ~ 1, data = patients.clinical_data) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability")+ 
  add_confidence_interval()+
  add_risktable()

# Log-rank test to compare survival between communities
logrank_test <- survdiff(Surv(os, status) ~ community, data = patients.clinical_data)
print(logrank_test)

# Since p=0.7, there is no evidence to reject the null hypothesis, 
# which states that the survival distributions of the two communities are the same.

summary(fit)

# Enrichment analysis on clinical data 
colnames(patients.clinical_data)
patients.clinical_data$ajcc_pathologic_stage # Tumor stage
patients.clinical_data$gender # gender
patients.clinical_data$status # Alive or dead
patients.clinical_data$age_at_diagnosis # Age at diagnosis

# Tumor stage
table_stage <- table(patients.clinical_data$community, patients.clinical_data$ajcc_pathologic_stage)
chisq_test <- chisq.test(table_stage)
fisher_test <- fisher.test(table_stage)
print(chisq_test)
print(fisher_test)

# gender
table_gender <- table(patients.clinical_data$community, patients.clinical_data$gender)
chisq_test_gender <- chisq.test(table_gender)
fisher_test <- fisher.test(table_gender)
print(chisq_test_gender)
print(fisher_test)

# Perform Kruskal-Wallis Test on age_at_diagnosis
kruskal_test_age <- kruskal.test(age_at_diagnosis ~ community, data = patients.clinical_data)
print(kruskal_test_age)

# Status
table_status <- table(patients.clinical_data$community, patients.clinical_data$status)
print(table_status)
chisq_test_status <- chisq.test(table_status)
print(chisq_test_status)


ggplot(patients.clinical_data, aes(x = community, fill = ajcc_pathologic_stage)) +
  geom_bar(position = "fill") +
  labs(title = "Tumor Stage Distribution Across Communities",
       x = "Community",
       y = "Proportion",
       fill = "Tumor Stage") +
  theme_minimal()
ggplot(patients.clinical_data, aes(x = community, y = age_at_diagnosis, fill = community)) +
  geom_boxplot() +
  labs(title = "Age Distribution by Community",
       x = "Community",
       y = "Age at Diagnosis") +
  theme_minimal()
######################


# C]
load("my_data2.RData") #mut.data
head(mut.data)
# create binary mutation matrix
BiocManager::install("maftools")
library(maftools,quietly = T)
install.packages("SNFtool")
library(SNFtool,quietly = T)
library(DT,quietly = T)
library(vegan,quietly = T)

###
proj <- "TCGA-THCA"  
quetab2query.mut <- GDCquery(
  project = proj, 
  data.category = "Simple Nucleotide Variation", 
  data.type = "Masked Somatic Mutation"
)
GDCdownload(quetab2query.mut)
mut.data <- GDCprepare(quetab2query.mut) 
# This data contains filtered mutation calls for somatic variations.
###

mut.data$Tumor_Sample_Barcode <- substr(mut.data$Tumor_Sample_Barcode, 1, 12)
colnames(filtr.expr.c) <- gsub("\\.", "-", colnames(filtr.expr.c))
mut.data <- mut.data[mut.data$Tumor_Sample_Barcode %in% colnames(filtr.expr.c),]

maf1 <- mut.data %>% read.maf


d=datatable(getSampleSummary(maf1), filter = 'top',
            options = list(scrollX = TRUE, keys = TRUE, pageLength = 10),  rownames = FALSE)
plotmafSummary(maf = maf1, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
oncoplot(maf = maf1, top = 10, removeNonMutated = TRUE) # displaying the top 10 most frequently mutated genes
lollipopPlot( maf = maf1,gene = 'BRAF',  showMutationRate = F)
lollipopPlot( maf = maf1,gene = 'NRAS',  showMutationRate = F)
somaticInteractions(maf = maf1, top = 20, pvalue = c(0.05, 0.01))
###########################################

# Create a binary mutation matrix (patients as columns, genes as rows)
mut_binary_matrix <- table( mut.data$Tumor_Sample_Barcode,mut.data$Hugo_Symbol) > 0 # each cell indicates whether a specific gene is mutated in a patient
mut_binary_matrix <- as.matrix(mut_binary_matrix)

similarity_matrix_mut <- 1 - as.matrix(vegdist(mut_binary_matrix, method = "jaccard", binary = TRUE))
diag(similarity_matrix_mut) <- 0  # Remove self-loops
head(similarity_matrix_mut)

# Normalize both similarity matrices
diag(similarity_matrix)=0
normalized_expr <- standardNormalization(similarity_matrix)
normalized_mut <- standardNormalization(similarity_matrix_mut)

common_ids <- intersect(rownames(normalized_expr), rownames(normalized_mut))
normalized_expr <- normalized_expr[common_ids, common_ids]
normalized_mut <- normalized_mut[common_ids, common_ids]

range(normalized_expr)
range(normalized_mut)

# Perform SNF
fused_network <- SNF(list(normalized_expr, normalized_mut), K = 30, t = 50)
diag(fused_network)=0
range(fused_network)


quantile(fused_network,0.8)
threshold <- 0.013  # Similarity threshold
fused_adjacency <- ifelse(fused_network >= threshold, 1, 0)

psn_2layer <- network(fused_adjacency, directed = FALSE)



psn_2l_comp=component.largest(psn_2layer,result="graph")
head(psn_2l_comp,1)

write.csv2(psn_2l_comp, "input-matrix.csv")
# python3 btc-community.py input-matrix.csv
psn_2l_comm <- read.csv2("output.txt", header = FALSE)
rownames(psn_2l_comm) <- rownames(psn_2l_comp)
length(table(psn_2l_comm[,1]))
table(psn_2l_comm[,1])

psn_2l_final <- network(psn_2l_comp, matrix.type="adjacency",ignore.eval = FALSE, names.eval = "weights", directed = F)

all(psn_2l_final %v% "vertex.names" == rownames(psn_2l_comm)) 
psn_2l_final  %v% "community" <-  as.character(psn_2l_comm[,1])

# Map community colors to the network nodes
node_colors_2l <- rainbow(length(unique(psn_2l_comm[, 1])),start=0.1)
names(node_colors_2l) <- unique(psn_2l_comm[, 1])
psn_2l_final %v% "color" <- node_colors_2l[psn_2l_final %v% "community"]

# Assign edge colors based on node attributes
network::set.edge.attribute(psn_2l_final, "edgecolor",
                            ifelse(as.numeric(psn_final %e% "in.community") == 1,"color", "white")  # Edges within communities = "color", others = "grey50"
)

# Find centroids of each community
centroids_psn_2l <- sapply(unique(psn_2l_comm[, 1]), function(community) {
  community_nodes <- which(psn_2l_final %v% "community" == community)
  degrees <- rowSums(as.matrix(psn_2l_final, matrix.type = "adjacency")[community_nodes, ])
  community_nodes[which.max(degrees)]  # Node with the highest degree in the community
})

# Get the names of the centroids
centroid_labels_psn_2l <- network.vertex.names(psn_2l_final)[centroids_psn_2l]

# Visualize the network
ggnet2(psn_2l_final,color = "color",size = 4,edge.color = c("color", "black"),  
       edge.alpha = 0.75)
#label = centroid_labels_psn_2l,label.color = "black",label.size = 4)






########################################

install.packages('intergraph')
install.packages('mclust')
library(mclust)
library(intergraph)
library(igraph)
# PSN 

if (all(rownames(psn_comm) == psn_final %v% "vertex.names")) {
  # Assign community information to psn_final
  psn_final %v% "community" <- as.character(psn_comm[, 1])
} else {
  stop("Vertex names in psn_final and psn_comm do not match.")
}
table(psn_final %v% "community")

igraph_network_1 <- asIgraph(psn_final)
modularity_value_1 <- modularity(igraph_network_1, as.numeric(psn_final %v% "community"))
cat("Modularity of the network:", modularity_value_1, "\n") # In this case, the modularity is 0.1544717 , suggesting a reasonably a weak community structure


# Calculate edge density for each community in both structures
density_3_communities_1 <- edge_density(induced_subgraph(igraph_network_1, V(igraph_network_1)$community == 1))  # Community 1 density 0.3181818
density_3_communities_2 <- edge_density(induced_subgraph(igraph_network_1, V(igraph_network_1)$community == 2))  # Community 2 density 0.9157895
density_3_communities_3 <- edge_density(induced_subgraph(igraph_network_1, V(igraph_network_1)$community == 3))  # Community 3 density 0.6666667


centroid_labels_psn # Centroid
centroid_expression <- filtr.expr.c[, centroid_labels_psn, drop = FALSE] # genes of centroid patients
top_genes <- apply(centroid_expression, 2, function(expr_values) { # top genes
  gene_index <- which.max(expr_values)  
  rownames(centroid_expression)[gene_index]  
})
names(top_genes) <- centroid_labels_psn 
print(top_genes)
# TCGA-BJ-A28R         TCGA-BJ-A2NA         TCGA-EL-A3ZS 
# "ENSG00000115414.21" "ENSG00000115414.21" "ENSG00000042832.12" 
#     (TFF3)               (TFF3)                 (IL6)



#########
# SNF
psn_2l_final
psn_2l_comm

if (all(rownames(psn_2l_comm) == psn_2l_final %v% "vertex.names")) {
  # Assign community information to psn_final
  psn_2l_final %v% "community" <- as.character(psn_2l_comm[, 1])
} else {
  stop("Vertex names in psn_final and psn_comm do not match.")
}
table(psn_2l_final %v% "community")

igraph_network_2 <- asIgraph(psn_2l_final)
modularity_value_2 <- modularity(igraph_network_2, as.numeric(psn_2l_final %v% "community"))
cat("Modularity of the network:", modularity_value_2, "\n") # In this case, the modularity is 0.4909765



density_4_communities <- sapply(1:4, function(x) {
  edge_density(induced_subgraph(igraph_network_2, V(igraph_network_2)$community == x))
})
cat("Density for 4-community structure:\n", density_4_communities, "\n") #  0.6623377 0.7628458 0.8 1 

centroid_labels_psn_2l # centroid
centroid_expression <- filtr.expr.c[, centroid_labels_psn_2l, drop = FALSE] # genes of centroid patients
top_genes <- apply(centroid_expression, 2, function(expr_values) { # top genes
  gene_index <- which.max(expr_values)  
  rownames(centroid_expression)[gene_index]  
})
names(top_genes) <- centroid_labels_psn_2l 
print(top_genes) 
#   TCGA-H2-A2K9         TCGA-EL-A3ZO         TCGA-H2-A3RI         TCGA-EL-A3MX 
# "ENSG00000042832.12" "ENSG00000115414.21" "ENSG00000042832.12"  "ENSG00000198804.2" 
#     (TG)               (FN1)                (TG)                   (MT-CO1)


# COMPARISON

# Gene Expression Only (PSN 1):
# Modularity = 0.1544717 Suggests a weak community structure
# 3 communities identified. The 3-community structure provides a lower-level partitioning with a very moderate dense communities
# Community 1 density 0.3181818
# Community 2 density 0.9157895
# Community 3 density 0.6666667

# SNF (Gene Expression + Mutation Profiles):
# Modularity = 0.491 Indicates a stronger but still weak community structure when integrating mutation profiles
# 4 communities identified. The 4-community structure reveals finer subgroups, with two highly cohesive groups (Communities 3 and 4) and two moderately cohesive ones (Communities 1 and 2).
# Community 1 0.6623377 
# Community 2 0.7628458 
# Community 3 0.8 
# Community 4 1 
# The SNF network creates a more selective and interpretable similarity structure, reducing the number of edges while retaining key connections.

# SNF Network has 52 nodes and 371 edges while PSN has 51 nodes and 523 edges
# The gene ENSG00000042832.12  is consistently the top-expressed gene for both networks in multiple centroid patients:
#PSN: It appears for TCGA-EL-A3ZS
# SNF-PSN: It remains a top gene for TCGA-H2-A2K9 and TCGA-H2-A3RI.
# Also the gene ENSG00000115414.21  is consistently top-expressed

# Let's investigate more about this gene
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Retrieve gene information
gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),
                   filters = "ensembl_gene_id",
                   values = "ENSG00000042832",
                   mart = ensembl)
print(gene_info)
# 1 ENSG00000042832 -> TG thyroglobulin [Source:HGNC Symbol;Acc:HGNC:11764]

# High levels of thyroglobulin are often observed in thyroid-related disorders, such as Thyroid cancers
# It is widely used as a tumor marker for monitoring thyroid cancer recurrence post-treatment.
# Reference https://academic.oup.com/ejendo/article/189/2/R11/7251404?


# Retrieve gene information
gene_info_2 <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),
                     filters = "ensembl_gene_id",
                     values = "ENSG00000115414",
                     mart = ensembl)
print(gene_info_2)
# 1 ENSG00000115414 -> FN1 fibronectin 1 [Source:HGNC Symbol;Acc:HGNC:3778]
# FN1 promotes cell movement and is critical in processes such as embryogenesis, angiogenesis, and metastasis of cancer cells.
# FN1 is often overexpressed in cancer tissues, promoting tumor progression, angiogenesis, and invasion.
# source: https://jeccr.biomedcentral.com/articles/10.1186/s13046-021-01908-8



