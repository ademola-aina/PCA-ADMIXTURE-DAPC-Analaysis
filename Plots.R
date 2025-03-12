#RVERSION 4.3.2 ################################## DEALING WITH ADMIXTURE DATA SETS #############################################################################
################################### Remember to always set your working directory appropriately ####################
##After you run admix and get the Q matrix results, attach the samples names to the first column and filter based on 70 or whatever % you prefer
##A new column Grp with samples assigned into unique grps based on your chosen % will be created
################EXAMPLE: Group Assignment Based on 0.7% Threshold from ADMIX Data #############################
# Assign groups based on the highest value > 0.7 among Q1:Q7
data <- read.csv("metadata@K7.csv", header = FALSE)
data$Grp <- apply(data[, 3:9], 1, function(row) {
  grp <- which.max(row)
  if (row[grp] > 0.7) {
    return(grp)
  } else {
    return(NA)
  }
})

# Print the data frame with the assigned groups
#print(data)
# Write out and save data to your local computer for future use
write.csv(data, "admix@K7data.csv", row.names = FALSE)

# To Create a data frame without NA (non-missing) in the Grp column
K7data_no_na <- data[!is.na(data$Grp), ]

# To Create a data frame with only rows that have NA (missing data) in the Grp column
K7data_na <- data[is.na(data$Grp), ]

# Write the data frames to CSV files saved to your working directory folder
write.csv(K7data_no_na, "K7data_no_na.csv", row.names = FALSE)
write.csv(K7data_na, "K7data_na.csv", row.names = FALSE)

# Check to see if truly those files were created and save properly 
cat("CSV files created:\n - K7data_no_na.csv (without NA)\n - K7data_na.csv (with only NA)\n")


##################################################################################################################################
#Quick CV plot from ADMIXTURE ANALYSIS
##################################################################################################################################
##
CV <- c(0.71025,0.69185,0.68317,0.67658,0.67369,0.67178,0.67238,0.67096,0.67307,0.67134,0.67354,0.67416)
plot(CV,type = "o", xlab = "K", ylab = "Cross-validation error")

#########################PCA plot###########################################################
#########################PCA FROM ADMIX RESULTS WITH ADMIXED INDVS REMOVED 
#Eigenvec and eigenval data from pca run from unix command line
#Install required libraries and run the scrip below
library(tidyverse)
library(ggplot2)
plinkPCA <- read_table("SNPs.QCnoadmix@k70.prune_pca.eigenvec", col_names = F)
plinkPCA <- plinkPCA[c(-1,-2)]
EigenValue <- scan("SNPs.QCnoadmix@k70.prune_pca.eigenval")
view(EigenValue)
##set column names
names(plinkPCA)[1:ncol(plinkPCA)] <- paste0("PC", 1:(ncol(plinkPCA)))
##Percentage variance explained
pve <- data.frame(PC = 1:20, pve = EigenValue/sum(EigenValue)*100)
##PCA plot
plinkPCA2 <- plinkPCA[,1:2]
mypop <- read.csv("Admixdata@K5_70Thresholdmini.csv", header = FALSE)
mypop$Grp3 <- as.factor(mypop$Grp3)
plinkPCA2$pop <- mypop$Grp3
Groups <- plinkPCA2$pop
####color by population use 5 colors 
rcolor <- c("gold", "darkmagenta", "red", "chocolate", "springgreen", "darkgreen", "blue")
####or color by states use 12 colors. remember to switch line 53 to column St for states
rcolor <- c("red", "chocolate", "cyan", "blue", "darkgreen", "darkmagenta","springgreen", "gold", "maroon2", "moccasin", "#FF7F00", "darkturquoise")
##PLOT PCA###
pca_plot <- ggplot(plinkPCA2, aes(x = PC1, y = PC2, color = Groups)) +
  geom_point(size = 2.5) +
  scale_color_manual(values = rcolor) + coord_equal()+
  theme_light()+ xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)"))+ ylab(paste0("PC2 (",signif(pve$pve[2], 3), "%)")) +
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(face = "bold"),
        axis.title = element_text(face = "bold"))
ggsave("pca_admix.pdf", plot = pca_plot, width = 7.5, height = 5.7, units = "in", dpi = 300)

######################################################################################################
##DAPC Analysis
######################################################################################################
library(adegenet)
library(vcfR)
library(ggplot2)
library(poppr)
vcf <- read.vcfR("SNPs.merged@K7Final.vcf", verbose = FALSE )
gid <- vcfR2genind(vcf, return.alleles = TRUE)
class(gid)
gid
grp <- find.clusters(gid, max.n.clust=40)
names(grp)
head(grp$Kstat, 8)
grp$stat
head(grp$grp, 30)
ind <- grp$grp
ind <- as.data.frame(ind)
write.csv(ind, "ind.csv", row.names = TRUE)
grp$size
#table(pop(gid), grp$grp)
#table.value(table(pop(gid), grp$grp), col.lab=paste("inf", 1:6),row.lab=paste("ori", 1:6))
dapc1 <- dapc(gid, grp$grp)
dapc1
scatter(dapc1)
scatter(dapc1, posi.da="bottomright", bg="white", pch=17:22)
myCol <- c("blue", "darkgreen", "red", "darkmagenta", "gold", "chocolate", "springgreen")
scatter(dapc1, posi.da="bottomright", bg="white",pch=17:22, cstar=0, col=myCol, scree.pca=FALSE,scree.da=FALSE,clabel = 0,
        posi.pca="bottomleft")

############################## MAKE A BETTER DAPC PLOT #############################################
pdf("DAPC_K7final")
scatter(dapc1, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=myCol, solid=.4,cex=1.7,clab=0, leg=TRUE, txt.leg=paste("Group",1:7))
dev.off()

####################################################################################################
##NJ TREE Analysis
#############################################################################################################
# Dist.matrix already generated in Linux PLINK 
#load distance matrix from PLINK
library(vcfR)
library(ape)
library(phangorn)
library(phyclust)
IBS<-read.table("Ferals.mdist") 
ID <-read.table("Ferals.mdist.id", row.names = 2)
# convert to Numeric
IBS.num <- sapply(IBS, as.numeric)
rownames(IBS.num) <- rownames(ID)
#save(IBS.num, file = "IBS.num.QC0.Rdata")
IBS.num.dist<-as.dist(IBS.num) # QC0
NJ.IBS<-nj(IBS.num.dist) # QC0
plot(NJ.IBS, cex = 0.1, type = "u", lab4ut = "axial")
class(NJ.IBS)
###Save to Newick file#####
nj_tree <- as.phylo(NJ.IBS)
# Read the group information from the CSV file
grp_info <- read.csv("Admixdata@K5_70Thresholdmini.csv", header = TRUE)

# Assign colors to groups (1 to 7)
#grp_colors <- c("gold", "darkmagenta", "red", "chocolate", "springgreen", "darkgreen", "blue")
grp_colors <- c("darkgreen", "gold", "darkmagenta", "red", "blue")
names(grp_colors) <- 1:5
# Extract tip labels from the tree
tip_labels <- nj_tree$tip.label
# Create a vector of colors for the tip labels based on the group information
tip_colors <- grp_colors[as.character(grp_info$Grp[match(tip_labels, grp_info$Acc)])]

# Plot the unrooted nj_tree with assigned colored tip labels
pdf("nj_tree_unrooted_colored.pdf")
plot(nj_tree, type = "unrooted", tip.color = tip_colors)
tiplabels(col = tip_colors, pch = 19, cex = 0.1)
dev.off()

write.tree(nj_tree, file = "nj_tree_unrooted_colored.nwk") ## Writes a Newick file

#Jump here from line 169 to color branches only
# Map colors to branches based on the tip colors
# Find the corresponding edges for each tip label
edge_colors <- rep("black", nrow(nj_tree$edge))
for (i in 1:length(tip_labels)) {
  tip_index <- which(nj_tree$tip.label == tip_labels[i])
  edge_indices <- which(nj_tree$edge[, 2] == tip_index)
  edge_colors[edge_indices] <- tip_colors[i]
}

# Plot the unrooted NJtree with colored branches and no tip labels
plot(nj_tree, type = "unrooted", edge.color = edge_colors, show.tip.label = FALSE)

# Optionally, you can also save the tree to a file
write.tree(nj_tree, file = "nj_tree_unrooted_colored.nwk")

############PUBLICATION READY ADMIXTURE PLOTS##################
library(tidyverse)
library(ggplot2)
tbl<- read.csv("K7data_no_na copy.csv", header = TRUE)
tbl2 <- tbl[,3:9]
plot_data <- tbl2 %>% 
  mutate(id=row_number()) %>% 
  gather("pop", "prob", Q1:Q7) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))
custom_colors <- c("Q1" = "gold", "Q2" = "darkmagenta", "Q3" = "red", "Q4" = "chocolate", "Q5" = "springgreen", "Q6" ="darkgreen", "Q7"="blue")
p <- ggplot(plot_data, aes(id, prob, fill = pop)) +
  geom_col() +
  scale_fill_manual(values = custom_colors) +
  theme_classic()  
  
P <- ggplot(plot_data, aes(id, prob, fill = pop)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = custom_colors) +
  facet_grid(~likely_assignment, scales = 'free', space = 'free') 

P + theme(strip.text.x.top = element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        rect = element_blank()) + 
  labs(x=NULL, y="Membership")
  

#####CHECKING ADMIXTURE DIRECTLY WITHOUT INITIAl GROUPING#######
library(tidyverse)
setwd("/Users/macbookair/Desktop/Rexercise")
tbl <- read.delim("Arl@K6.txt", header = FALSE)
plot_data <- tbl %>% 
  mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V4) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

####with GGplot#####################################################
ggplot(plot_data, aes(id, prob, fill = pop)) +
  geom_col() +
  theme_classic()

#####With Facets#####################################################
ggplot(plot_data, aes(id, prob, fill = pop)) +
  geom_col() +
  facet_grid(~likely_assignment, scales = 'free', space = 'free')

######################################################################
#######Extract samples names from Groups##############################
library(dplyr)
# Read the CSV file into a data frame
data <- read.csv("K7data_no_na.csv")
# Filter the data to keep only the required rows
grp_data <- data %>%
  filter(Grp == 7)
# Print the filtered data
print(grp_data)
# Save the filtered data to a new CSV file
write.csv(grp_data, "grp_7_samples.csv", row.names = FALSE)




  
