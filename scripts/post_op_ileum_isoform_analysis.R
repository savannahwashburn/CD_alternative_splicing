#Script for analysis post-op dataset 

#load in libraries 
library(lme4)
library(readxl)
library(dplyr)
library(glmmTMB)
library(ggvenn)
library(ggplot2)
library(ggrepel)
library(SummarizedExperiment)
library(pvca)
library(multcomp)
library(readr)


#calcualte the mean difference of IF (dIF)
calc_dIF <- function(df){
  mean_IF <- aggregate(df[,1],list(df$affected), FUN=mean)
  dIF_each <- mean_IF[mean_IF$Group.1=='R1','x'] - mean_IF[mean_IF$Group.1=='R0','x']
  dIF <- c(names(df)[1],dIF_each)
  
  
  for (m in 2:(ncol(df)-1)){
    mean_IF <- aggregate(df[,m],list(df$affected), FUN=mean)
    dIF_each <- mean_IF[mean_IF$Group.1=='R1','x'] - mean_IF[mean_IF$Group.1=='R0','x']
    dIF <- rbind(dIF,c(names(df)[m],dIF_each))
    
  }
  dIF <- as.data.frame(dIF)
  colnames(dIF) <- c('transcript_id','dIF')
  row.names(dIF) <- dIF$transcript_id
  dIF$dIF <- as.numeric(dIF$dIF)
  return(dIF)
}

#------read in phenotype data--------# 

pheno <- read_excel("/Users/swashburn30/Desktop/isoform/validate_kates_results/PhenoSummary.xlsx",sheet = 'Data')
pheno <- as.data.frame(pheno)
rownames(pheno) <- make.names(pheno$RNAid, unique = TRUE)
#row.names(pheno)<-pheno$RNAid
pheno$Batch <- as.factor(pheno$Batch)
pheno <- pheno[,c("affected","PhenoSex", "Batch",'Smoking','race',"AgeAtSample","consortium_id")]
pheno <-na.omit(pheno) #334 samples instead of 335


#----------read in isoform fraction-----------#


IF_df <- read.table('/Users/swashburn30/Desktop/isoform/validate_kates_results/washburn_isoform_fraction.tsv',sep='\t', header=TRUE)
row.names(IF_df) <- IF_df$transcript_id
gene_mapping <- IF_df %>% dplyr::select(transcript_id,gene_id)
IF_df <- IF_df%>% dplyr::select(row.names(pheno)) #select samples - now there are the same in pheno and iso table
IF_df <- na.omit(IF_df)  #remove NA (where the tpm for gene is 0)
IF_df <- IF_df[rowSums(IF_df != 0) > (0.05*ncol(IF_df)), ] #keep transcript only no more than 95% ids have non-zero IF
IF_df <- IF_df[rowSums(IF_df < 1) > (0.95*ncol(IF_df)), ] #keep transcript only 95% individuals have it IF < 1


IF_df <- IF_df[,row.names(pheno)] #reorder - pheno and if_df samples are in same order
IF_trans <- as.data.frame(t(IF_df)) #transpose


#combine pheno and IF_df tables - test for just one IF 
IF_design <- cbind(pheno[,c("affected","PhenoSex", "Batch",'Smoking','race',"AgeAtSample","consortium_id")],IF_trans[,1])
names(IF_design) <- c("affected","PhenoSex", "Batch",'Smoking','race',"AgeAtSample","consortium_id","IF")


IF_design <- na.omit(IF_design) #remove na values

############## linear mixed model, gaussian distribution#######

model.full = lmer(IF ~ affected + PhenoSex + Batch + Smoking + race + AgeAtSample + (1|consortium_id),  data=IF_design ,REML = FALSE) #full model

model.affected = lmer(IF ~ PhenoSex + Batch + Smoking + race + AgeAtSample + (1|consortium_id),  data=IF_design, REML = FALSE ) #model missing affected

model.Batch = lmer(IF ~ affected + PhenoSex + Smoking + race + AgeAtSample + (1|consortium_id),  data=IF_design, REML = FALSE ) #full model w/o batch

model.race = lmer(IF ~ affected + PhenoSex + Batch + Smoking + AgeAtSample + (1|consortium_id), data=IF_design,REML = FALSE ) 


#extract pvalue using maximum likelyhood
prob_affected <- anova(model.full,model.affected)[2,8]
prob_batch <- anova(model.full,model.Batch)[2,8]
prob_race <- anova(model.full,model.race)[2,8]
pvals <- matrix(c(prob_affected,prob_batch,prob_race),nrow=1,ncol=3)

#loop through the rest of isoforms
for (m in 2:ncol(IF_trans)){
  IF_design <- cbind(pheno[,c("affected","PhenoSex", "Batch",'Smoking','race',"AgeAtSample","consortium_id")],IF_trans[,m])
  names(IF_design) <- c("affected","PhenoSex", "Batch",'Smoking','race',"AgeAtSample","consortium_id","IF")
  
  IF_design <- na.omit(IF_design) #remove na values
  
  
  model.full = lmer(IF ~ affected + PhenoSex + Batch + Smoking + race + AgeAtSample + (1|consortium_id),  data=IF_design, REML = FALSE ) #full model
  
  model.affected = lmer(IF ~ PhenoSex + Batch + Smoking + race + AgeAtSample + (1|consortium_id),  data=IF_design, REML = FALSE ) #model missing affected
  
  model.Batch = lmer(IF ~ affected + PhenoSex + Smoking + race + AgeAtSample + (1|consortium_id), data=IF_design, REML = FALSE ) #full model w/o batch
  
  model.race = lmer(IF ~ affected + PhenoSex + Batch + Smoking + AgeAtSample + (1|consortium_id),  data=IF_design, REML = FALSE ) 
  
  
  #extract pvalue using maximum likelyhood
  prob_affected <- anova(model.full,model.affected)[2,8]
  prob_batch <- anova(model.full,model.Batch)[2,8]
  prob_race <- anova(model.full,model.race)[2,8]
  pvals <- rbind(pvals, c(prob_affected,prob_batch,prob_race))
}

#change the rownames 
row.names(pvals) <- row.names(IF_df)
colnames(pvals) <- c('affected.p','batch.p','race.p')

#adjust p value - calcualte 
pvals_gaussian <- as.data.frame(pvals)
pvals_gaussian$affected.adjp <- p.adjust(pvals_gaussian[,1], method='fdr')
pvals_gaussian$batch.adjp <- p.adjust(pvals_gaussian[,2], method='fdr')
pvals_gaussian$race.adjp <- p.adjust(pvals_gaussian[,3], method='fdr')
pvals_gaussian$transcript_id<- row.names(pvals_gaussian)
pvals_gaussian<- left_join(pvals_gaussian,gene_mapping,by='transcript_id')


###########calculate dIF: mean (IF(recurrence)) - mean (IF(non-recurrence))##########
IF_trans$affected <- pheno$affected
#mean of just the first isoform, for each group
mean_IF <- aggregate(IF_trans[,1],list(IF_trans$affected), FUN=mean)
dIF_each <- mean_IF[mean_IF$Group.1=='R1','x'] - mean_IF[mean_IF$Group.1=='R0','x']
dIF <- c(names(IF_trans)[1],dIF_each)


for (m in 2:(ncol(IF_trans)-1)){
  mean_IF <- aggregate(IF_trans[,m],list(IF_trans$affected), FUN=mean)
  dIF_each <- mean_IF[mean_IF$Group.1=='R1','x'] - mean_IF[mean_IF$Group.1=='R0','x']
  dIF <- rbind(dIF,c(names(IF_trans)[m],dIF_each))
  
}
dIF <- as.data.frame(dIF)
colnames(dIF) <- c('transcript_id','dIF')
row.names(dIF) <- dIF$transcript_id
dIF$dIF <- as.numeric(dIF$dIF)



###############find overlapping genes in DEG and DTU#########

#DTU cut-off: affected adjp by FDR <0.05 and abs(dIF) >0.1
pvals_gaussian <- left_join(pvals_gaussian,dIF,by='transcript_id')
DTU <- pvals_gaussian[pvals_gaussian$affected.adjp<0.05 & (pvals_gaussian$dIF>0.1| pvals_gaussian$dIF< (-0.1)),]
DTU_relaxed <- pvals_gaussian[pvals_gaussian$affected.adjp<0.05,]


#DEG cut-off: adjp by FDR <0.05 and log2FC > log2(1.5) or log2FC < log2(-1.5)
DEG <- read.table('/Users/swashburn30/Desktop/isoform/validate_kates_results/recur0vs1.csv',sep=',',header=TRUE)
names(DEG)[1] <- 'gene_id'
DEG <- DEG[(DEG$log2FoldChange<log2(1/1.5) | DEG$log2FoldChange>log2(1.5/1)) & DEG$padj<0.05,]

#draw venn diagram
df_venn <- data.frame('gene_id' = unique(c(DTU$gene_id, DEG$gene_id)))
names(df_venn)[1] <- 'gene_id'
df_venn$DEG <- df_venn$gene_id %in% DEG$gene_id
df_venn$DTU <- df_venn$gene_id %in% DTU$gene_id

venn_gaussian <- ggplot(df_venn, # Apply geom_venn function
                        aes( A = DEG , B = DTU)) +
  geom_venn() +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle('overlapping genes between differential expression gene and differential transcript usage gene')

#draw venn diagram for relaxed DTU criteria
df_venn_relaxed <- data.frame('gene_id' = unique(c(DTU_relaxed$gene_id, DEG$gene_id)))
names(df_venn_relaxed)[1] <- 'gene_id'
df_venn_relaxed$DEG <- df_venn_relaxed$gene_id %in% DEG$gene_id
df_venn_relaxed$DTU <- df_venn_relaxed$gene_id %in% DTU_relaxed$gene_id

venn_gaussian_relaxed <- ggplot(df_venn_relaxed, # Apply geom_venn function
                                aes( A = DEG , B = DTU)) +
  geom_venn() +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ggtitle('overlapping genes between differential expression gene and differential transcript usage gene')


#save table with pvals and dIF
write.table(pvals_gaussian,'/Users/swashburn30/Desktop/isoform/validate_kates_results/DTU_R1R0.tsv',sep='\t',quote=FALSE,dec='.',row.names = FALSE)


###########volcano plot############
df_volcano <- pvals_gaussian
df_volcano$group <- 'NS'
df_volcano$group[df_volcano$dIF< (-0.1) & df_volcano$affected.adjp<0.05] <-'DOWN'
df_volcano$group[df_volcano$dIF > (0.1) & df_volcano$affected.adjp<0.05] <- 'UP'
df_volcano$gene_transcript <- paste0(df_volcano$gene_id,': ',df_volcano$transcript_id)
df_volcano$label[df_volcano$group !='NS'] <- df_volcano$gene_transcript[df_volcano$group != 'NS'] #gene label

volcano_gaussian <- ggplot(data=df_volcano,mapping=aes(x=dIF,y=-log10(affected.adjp),color=group,label=label))+
  geom_point()+
  xlim(-0.5,0.5) +
  labs(title='volcano plot of significance against difference in isoform fraction between recurrence and non-recurrence') +
  geom_label_repel(na.rm = TRUE, show.legend = FALSE, box.padding = 0.5) +
  theme_minimal()+ 
  theme(legend.position="none")+
  geom_vline(xintercept=c(-0.1, 0.1), col="red", linetype='longdash') +
  geom_hline(yintercept=-log10(0.05), col="red", linetype='longdash') +
  geom_vline(xintercept = 0, col="grey", linetype='longdash') +
  scale_color_manual(values=c("DOWN"= "#56B4E9", "NS"="#999999", "UP"="#E69F00"))


############PCA score#########
library(ggpubr)
library(glue)
library(dplyr)
library(ggpubr)
#stringent DTU filter: adjp < 0.05 and dIF > ± 0.1
mat <- IF_trans%>% select(DTU$transcript_id)
mat <- as.matrix(mat)
res.pca.DTU_stringent <- prcomp(mat)
df_pca.DTU_stringent <- as.data.frame(res.pca.DTU_stringent$x)
df_pca.DTU_stringent <- df_pca.DTU_stringent %>% select(PC1)

unique(row.names(df_pca.DTU_stringent) == row.names(pheno)) #checking rna id matches


df_pca.DTU_stringent$affected <- pheno$affected
df_pca.DTU_stringent$PhenoSex <- pheno$PhenoSex

#24 DTUs
filt <- IF_df[rownames(IF_df) %in% c("NM_001244583.1", "XM_005245796.2", "XM_006717970.1", "NM_000084.4", "XM_005251975.1", "XM_005274390.2", "rna85636", "XM_006716789.1", "NM_019101.2", "NM_201446.2", "NM_001286688.1", "NM_015577.2", "NM_001127216.1", "NM_001190796.1", "NM_012252.3", "NM_173157.2", "NM_024756.2", "NM_001145522.1", "XM_005259384.2", "XM_005258201.1", "NM_170782.2", "NM_001127215.1", "NM_001135940.1", "NM_001006638.2"), ]
gene_filt <- gene_mapping[gene_mapping$transcript_id %in% rownames(IF_df), ]
IF_df$gene_id <- "NA"
IF_df$gene_id <- gene_filt$gene_id
IF_df$transcript_id <- gene_filt$transcript_id
transcript_filt <- gene_filt[gene_filt$transcript_id %in% rownames(filt), ]
filt$transcript <- transcript_filt$transcript_id
filt$gene <- transcript_filt$gene_id

filt <- filt[, -c(335,336)]
filt_t <- t(filt)
df_pca <- prcomp(filt_t)

all(row.names(df_pca) == row.names(pheno))
df_pca$affected <- pheno$affected

#variance explained by each PC
pc_eigenvalues <- df_pca$sdev^2

pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_eigenvalues

# The PC scores are stored in the "x" value of the prcomp object
pc_scores <- df_pca$x

pc_scores <- pc_scores %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")

pc_scores

#add metadata information 
#read in sample and recurrence vs non-recurrence data 
pheno <- read.csv("/Users/swashburn30/Desktop/isoform/kiera_analysis/phenodata.csv", header = F)

colnames(pheno) <- c("sra_id", "geo_id", "sex", "tissue", "sample")

pc_score_metadata<- full_join(pc_scores, pheno, by = ("sample"))
test_pc_score_metadata <- na.omit(pc_score_metadata)

pheno$sample <- rownames(pheno)
pc_score_metadata <- full_join(pc_scores, pheno, by = ("sample"))

pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(affected))) +
  geom_point() 
#geom_text(aes(label = sample))
pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(Batch))) +
  geom_point() 
pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(PhenoSex))) +
  geom_point() 
pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(race))) +
  geom_point() 



boxplot_stringent_affected <- ggplot(df_pca.DTU_stringent, aes(x=affected, y=PC1, fill= affected)) + 
  geom_boxplot()+
  labs(subtitle="DTU PC score in all samples (n=335)",x="Recurrence", y = "PC1")+
  theme_minimal() +
  scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  stat_compare_means(method = "wilcox.test", label.x = 2) +
  theme(legend.position="none",axis.title.x=element_blank(),text = element_text(size=12),axis.text = element_text(size=12))+
  scale_x_discrete(labels=c("R0" = glue("R0 (n={nrow(df_pca.DTU_stringent[df_pca.DTU_stringent$affected=='R0',])})"), "R1" = glue("R1 (n={nrow(df_pca.DTU_stringent[df_pca.DTU_stringent$affected=='R1',])})")))



boxplot_stringent_sex <- ggplot(df_pca.DTU_stringent, aes(x=PhenoSex, y=PC1, fill= PhenoSex)) + 
  geom_boxplot()+
  labs(title="PC1 score of differential transcript usage transcripts in male and female",x="Recurrence", y = "PC1")+
  theme_minimal() +
  scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  stat_compare_means(method = "wilcox.test", label.x = 2) +
  theme(legend.position="none",axis.title.x=element_blank(),text = element_text(size=12),axis.text = element_text(size=12))









#relaxed DTU filter: adjp < 0.05 
mat <- IF_trans%>% select(DTU_relaxed$transcript_id)
mat <- as.matrix(mat)
res.pca.DTU_relaxed <- prcomp(mat)
df_pca.DTU_relaxed <- as.data.frame(res.pca.DTU_relaxed$x)
df_pca.DTU_relaxed <- df_pca.DTU_relaxed %>% select(PC1)

unique(row.names(df_pca.DTU_relaxed) == row.names(pheno)) #checking rna id matches


#relaxed DTU filter: adjp < 0.05 
df_pca.DTU_relaxed$affected <- pheno$affected
df_pca.DTU_relaxed$PhenoSex <- pheno$PhenoSex


boxplot_relaxed_affected <- ggplot(df_pca.DTU_relaxed, aes(x=affected, y=PC1, fill= affected)) + 
  geom_boxplot()+
  labs(title="PC1 score of differential transcript usage in recurrence and non-recurrence group",x="Recurrence", y = "PC1")+
  theme_minimal() +
  scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  stat_compare_means(method = "wilcox.test", label.x = 2) +
  theme(legend.position="none",axis.title.x=element_blank(),text = element_text(size=12),axis.text = element_text(size=12))



boxplot_relaxed_sex <- ggplot(df_pca.DTU_relaxed, aes(x=PhenoSex, y=PC1, fill= PhenoSex)) + 
  geom_boxplot()+
  labs(title="PC1 score of differential transcript usage transcripts in male and female",x="Recurrence", y = "PC1")+
  theme_minimal() +
  scale_color_manual(values=c("#E69F00", "#56B4E9"))+
  stat_compare_means(method = "wilcox.test", label.x = 2) +
  theme(legend.position="none",axis.title.x=element_blank(),text = element_text(size=12),axis.text = element_text(size=12))


#--------------------rectum vs. ileum DTUs stratify r0 vs. r1?-------------------#

#------read in phenotype data--------# 
library(readxl)
pheno <- read_excel("/Users/swashburn30/Desktop/isoform/validate_kates_results/PhenoSummary.xlsx",sheet = 'Data')
pheno <- as.data.frame(pheno)
rownames(pheno) <- make.names(pheno$RNAid, unique = TRUE)
#row.names(pheno)<-pheno$RNAid
pheno$Batch <- as.factor(pheno$Batch)
pheno <- pheno[,c("affected","PhenoSex", "Batch",'Smoking','race',"AgeAtSample","consortium_id")]
pheno <-na.omit(pheno) #334 samples instead of 335


#----------read in isoform fraction-----------#


IF_df <- read.table('/Users/swashburn30/Desktop/isoform/validate_kates_results/washburn_isoform_fraction.tsv',sep='\t', header=TRUE)
row.names(IF_df) <- IF_df$transcript_id
gene_mapping <- IF_df %>% dplyr::select(transcript_id,gene_id)
IF_df <- IF_df%>% dplyr::select(row.names(pheno)) #select samples - now there are the same in pheno and iso table
IF_df <- na.omit(IF_df)  #remove NA (where the tpm for gene is 0)
#for this part - did not filter out the IFs (just filter based on transcript_id in tissue_DTU)

IF_df <- IF_df[,row.names(pheno)] #reorder - pheno and if_df samples are in same order


#load in ileum vs. rectum DTUs 
tissue_DTU <- read_csv(file = "/Users/swashburn30/Desktop/isoform/kiera_analysis/sig_DTUs_2_9_24.csv")

#subset IF_df/IF_trans based on tissue DTUs 
IF_df <- IF_df[rownames(IF_df) %in% tissue_DTU$transcript_id, ] #only 201 DTU from tissue instead of 203
IF_trans <- as.data.frame(t(IF_df)) #transpose

#perform PCA using the tissue DTUs 
df_pca <- prcomp(IF_trans)

all(row.names(df_pca) == row.names(pheno)) #check that samples match ##TRUE

#variance explained by each PC
pc_eigenvalues <- df_pca$sdev^2

pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_eigenvalues #PC1: 11.6% of variance; PC2: 9.06% of variance 

# The PC scores are stored in the "x" value of the prcomp object
pc_scores <- df_pca$x

pc_scores <- pc_scores %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")

pc_scores

#add sample name column 
pheno$sample <- rownames(pheno)
pc_score_metadata <- full_join(pc_scores, pheno, by = ("sample"))

#create PCA plot 
pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(affected))) +
  geom_point() 

#create PCA plot - PC1 is just batch 
pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(Batch))) +
  geom_point() 


#save the PC scores with metadata
write.csv(pc_score_metadata, file = "/Users/swashburn30/Desktop/isoform/validate_kates_results/tissue_DTU_pca_results_r0_vs_r1.csv")

#which PC is associated with recurrence status?
pc1 <- lm(pc_score_metadata$PC1 ~ affected, data = pc_score_metadata)
pc2 <- lm(pc_score_metadata$PC2 ~ affected, data = pc_score_metadata)
pc3 <- lm(pc_score_metadata$PC3 ~ affected, data = pc_score_metadata)
pc4 <- lm(pc_score_metadata$PC4 ~ affected, data = pc_score_metadata)
pc5 <- lm(pc_score_metadata$PC5 ~ affected, data = pc_score_metadata)
pc6 <- lm(pc_score_metadata$PC6 ~ affected, data = pc_score_metadata)
pc7 <- lm(pc_score_metadata$PC7 ~ affected, data = pc_score_metadata)
pc8 <- lm(pc_score_metadata$PC8 ~ affected, data = pc_score_metadata)
pc9 <- lm(pc_score_metadata$PC9 ~ affected, data = pc_score_metadata)
pc10 <- lm(pc_score_metadata$PC10 ~ affected, data = pc_score_metadata)

#PC2 and PC5 are associated with recurrence status

#Create Boxplot of PC2 score 
p1 <- pc_score_metadata %>%
  ggplot( aes(x=affected, y=PC2, fill=affected)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=20),
    text = element_text(size = 20)
  ) +
  ggtitle("DTU PC1 score in all samples (n=334)") +
  stat_compare_means(method = "wilcox.test", label.x = 2) +
  xlab("Recurrence Status") +
  ylab("PC2") +
  scale_x_discrete(labels=c("R0" = "R0 (n=238)", "R1" = "R1 (n=96)"))

p1 + stat_compare_means(method = "wilcox.test", label.x = 2)


#extract top 100 genes from PC2 

pca_gene <- df_pca$rotation
pca_gene <- as.data.frame(pca_gene)
pca_gene$gene <- rownames(pca_gene)

#select top genes 
top_gene_PC2 <- pca_gene %>% 
  select(gene, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  #group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>%
  slice(1:100) %>% 
  # pull the gene column as a vector
  #pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

gene_pc2 <- top_gene_PC2$gene

#match isoform to gene
gene_pc2 <- data.frame(gene_pc2)
colnames(gene_pc2) <- "isoform"
rownames(gene_pc2) <- gene_pc2$isoform

iso_gene <- iso[, 1, drop = FALSE]
iso_gene$isoform <- rownames(iso_gene)

iso_gene_sub <- subset(iso_gene, isoform %in% gene_pc2$isoform)

# Get the row names of df1
row_names_pc2 <- rownames(gene_pc2)

# Use the row names of df1 to reorder df2
iso_gene_sub <- iso_gene_sub[row_names_pc2, ]

unique_pc2 <- unique(iso_gene_sub$gene_id)

#look at top 200 isoforms
#select top genes 
top_gene_PC2_200 <- pca_gene %>% 
  select(gene, PC2) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  #group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>%
  slice(1:200) %>% 
  # pull the gene column as a vector
  #pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

gene_pc2_200 <- top_gene_PC2_200$gene

#match isoform to gene
gene_pc2_200 <- data.frame(gene_pc2_200)
colnames(gene_pc2_200) <- "isoform"
rownames(gene_pc2_200) <- gene_pc2_200$isoform

iso_gene <- iso[, 1, drop = FALSE]
iso_gene$isoform <- rownames(iso_gene)

iso_gene_sub_200 <- subset(iso_gene, isoform %in% gene_pc2_200$isoform)

# Get the row names of df1
row_names_pc2_200 <- rownames(gene_pc2_200)

# Use the row names of df1 to reorder df2
iso_gene_sub_200 <- iso_gene_sub_200[row_names_pc2_200, ]

unique_pc2_200 <- unique(iso_gene_sub_200$gene_id)


#save the pc2 genes 
write.csv(unique_pc2, file = "/Users/swashburn30/Desktop/isoform/pc2_genes.csv")

#save top 200 pc2 genes
write.csv(unique_pc2_200, file = "/Users/swashburn30/Desktop/isoform/pc2_genes_top_200.csv")


#save top 

#identify pathways in PC2 
library(DOSE)
library(pathview)
library(clusterProfiler)

#boferroni - more stringent threshold 
ego <- enrichGO(gene = unique_pc2, 
                #universe = allOE_genes,
                keyType = "SYMBOL",
                OrgDb = "org.Hs.eg.db", 
                ont = "BP", 
                pAdjustMethod = "bonferroni", 
                qvalueCutoff = 0.05) 

#save as data table 
go_summary_pc2 <- data.frame(ego)

#visualize results
dotplot(ego)
barplot(ego)

#pathways from ToppFunn - PC2
pathways_pc2 <- read.table(file = "/Users/swashburn30/Desktop/isoform/pc2_pathways_4_18_24.txt", sep = '\t', header = TRUE)

#test making ggplot barplot 
pathways_pc2$q.value.Bonferroni <- as.numeric(pathways_pc2$q.value.Bonferroni)
pathways_pc2$log10 <- -log10(pathways_pc2$q.value.Bonferroni)

ggplot(pathways_pc2, aes(x = reorder(Name, +Hit.Count.in.Query.List), y = Hit.Count.in.Query.List, fill = q.value.Bonferroni)) +
  geom_bar(stat = 'identity') +
  xlab("Pathways") +
  ylab("Count") +
  scale_fill_gradient(low = "blue", high = "red") +
  theme(
    text = element_text(size = 15)
  ) +
  coord_flip()
test <- scale_fill_brewer(palette="Dark2") 

#top 200 pc2 pathways 
pathways_pc2_200 <- read.table(file = "/Users/swashburn30/Desktop/isoform/pc2_top200_pathways.txt", sep = '\t', header = TRUE)
#test making ggplot barplot 
pathways_pc2_200$q.value.Bonferroni <- as.numeric(pathways_pc2_200$q.value.Bonferroni)

ggplot(pathways_pc2_200, aes(x = reorder(Name, +Hit.Count.in.Query.List), y = Hit.Count.in.Query.List, fill = q.value.Bonferroni)) +
  geom_bar(stat = 'identity') +
  xlab("Pathways") +
  ylab("Count") +
  scale_fill_gradient(low = "blue", high = "red") +
  theme(
    text = element_text(size = 15)
  ) +
  coord_flip()


#---------analysis at gene level------------#

#read in pheno data
pheno <- read_excel("/Users/swashburn30/Desktop/isoform/validate_kates_results/PhenoSummary.xlsx",sheet = 'Data')
pheno <- as.data.frame(pheno)
row.names(pheno)<-pheno$RNAid
pheno$Batch <- as.factor(pheno$Batch)
pheno <- pheno[,c("affected","PhenoSex", "Batch",'Smoking','race',"AgeAtSample","consortium_id")]
pheno <-na.omit(pheno)

#read in gene expression data 
gene_df <- read.table('/Users/swashburn30/Desktop/isoform/validate_kates_results/washburn_rsem.merged.gene_tpm.tsv', sep = '\t', header = TRUE)
row.names(gene_df) <- gene_df$gene_id
gene_df <- gene_df[, -1] #remove the rownames (it was just numbers)
gene_mapping <- gene_df %>% dplyr::select(transcript_id.s.,gene_id) 
gene_df <- gene_df%>% dplyr::select(row.names(pheno)) #select samples - 334 samples
gene_df <- na.omit(gene_df)  #remove NA (where the tpm for gene is 0) #29972 genes
gene_df <- gene_df[ rowMeans(gene_df) > 5, ] #mean coverage > 5 reads - 5505 genes

gene_df <- gene_df[,row.names(pheno)] #reorder
gene_trans <- as.data.frame(t(gene_df)) #transpose

#scale the log2 normalize the gene expression level 
gene_trans <- log2(gene_trans + 1)

#create the design matrix for linear model 
gene_design <- cbind(pheno[,c("affected","PhenoSex", "Batch",'Smoking','race',"AgeAtSample","consortium_id")],gene_trans[,1])
names(gene_design) <- c("affected","PhenoSex", "Batch",'Smoking','race',"AgeAtSample","consortium_id","gene")


gene_design <- na.omit(gene_design) #remove na values

#try removing row 3956 and see if this helps
gene_trans <- gene_trans[ ,-3956] #RNU5A-1
gene_df <- gene_df[-3956, ]

############## linear mixed model, gaussian distribution#######

model.full = lmer(gene ~ affected + PhenoSex + Batch + Smoking + race + AgeAtSample + (1|consortium_id),  data=gene_design ,REML = FALSE) #full model

model.affected = lmer(gene ~ PhenoSex + Batch + Smoking + race + AgeAtSample + (1|consortium_id),  data=gene_design, REML = FALSE ) #model missing affected

model.Batch = lmer(gene ~ affected + PhenoSex + Smoking + race + AgeAtSample + (1|consortium_id),  data=gene_design, REML = FALSE ) #full model w/o batch

model.race = lmer(gene ~ affected + PhenoSex + Batch + Smoking + AgeAtSample + (1|consortium_id), data=gene_design,REML = FALSE ) 


#extract pvalue using maximum likelyhood
prob_affected <- anova(model.full,model.affected)[2,8]
prob_batch <- anova(model.full,model.Batch)[2,8]
prob_race <- anova(model.full,model.race)[2,8]
pvals <- matrix(c(prob_affected,prob_batch,prob_race),nrow=1,ncol=3)



#loop through the rest of isoforms
for (m in 2:ncol(gene_trans)){
  gene_design <- cbind(pheno[,c("affected","PhenoSex", "Batch",'Smoking','race',"AgeAtSample","consortium_id")],gene_trans[,m])
  names(gene_design) <- c("affected","PhenoSex", "Batch",'Smoking','race',"AgeAtSample","consortium_id","gene")
  
  gene_design <- na.omit(gene_design) #remove na values
  
  
  model.full = lmer(gene ~ affected + PhenoSex + Batch + Smoking + race + AgeAtSample + (1|consortium_id),  data=gene_design, REML = FALSE ) #full model
  
  model.affected = lmer(gene ~ PhenoSex + Batch + Smoking + race + AgeAtSample + (1|consortium_id),  data=gene_design, REML = FALSE ) #model missing affected
  
  model.Batch = lmer(gene ~ affected + PhenoSex + Smoking + race + AgeAtSample + (1|consortium_id), data=gene_design, REML = FALSE ) #full model w/o batch
  
  model.race = lmer(gene ~ affected + PhenoSex + Batch + Smoking + AgeAtSample + (1|consortium_id),  data=gene_design, REML = FALSE ) 
  
  
  #extract pvalue using maximum likelyhood
  prob_affected <- anova(model.full,model.affected)[2,8]
  prob_batch <- anova(model.full,model.Batch)[2,8]
  prob_race <- anova(model.full,model.race)[2,8]
  pvals <- rbind(pvals, c(prob_affected,prob_batch,prob_race))
}

row.names(pvals) <- row.names(gene_df)
colnames(pvals) <- c('affected.p','batch.p','race.p')




#adjust p value
pvals_gaussian <- as.data.frame(pvals)
pvals_gaussian$affected.adjp <- p.adjust(pvals_gaussian[,1], method='fdr')
pvals_gaussian$batch.adjp <- p.adjust(pvals_gaussian[,2], method='fdr')
pvals_gaussian$race.adjp <- p.adjust(pvals_gaussian[,3], method='fdr')
pvals_gaussian$gene_id<- row.names(pvals_gaussian)
pvals_gaussian<- left_join(pvals_gaussian,gene_df,by=rownames(gene_df))

#-------calculate log 2 FC------------#

# Subset the dataframe based on disease status
#combine pheno data with gene expression data
meta_gene <- cbind(gene_trans,pheno)

r1_samples <- meta_gene[meta_gene$affected == "R1", ]
r0_samples <- meta_gene[meta_gene$affected == "R0", ]

# Calculate the mean expression for each gene in each condition
mean_expression_r1 <- as.data.frame(colMeans(r1_samples[, -c(5504, 5505,5506,5507,5508,5509,5510,5511)], na.rm = TRUE))
colnames(mean_expression_r1) <- "R1"
mean_expression_r0 <- as.data.frame(colMeans(r0_samples[, -c(5504, 5505,5506,5507,5508,5509,5510,5511)], na.rm = TRUE))
colnames(mean_expression_r0) <- "R0"

# Calculate the log2-fold change - already log2 normalized 

test <- (mean_expression_r1 - mean_expression_r0)
test$gene_id <- rownames(test)
colnames(test) <- c("log2FC", "gene_id")

###########volcano plot############

pvals_gaussian <- left_join(pvals_gaussian,test,by='gene_id')
df_volcano <- pvals_gaussian
df_volcano$group <- 'NS'
df_volcano$group[df_volcano$log2FC< (-0.25) & df_volcano$affected.adjp<0.05] <-'DOWN'
df_volcano$group[df_volcano$log2FC > (0.25) & df_volcano$affected.adjp<0.05] <- 'UP'
#df_volcano$gene_transcript <- paste0(df_volcano$gene_id,': ',df_volcano$transcript_id)
df_volcano$label[df_volcano$group !='NS'] <- df_volcano$gene_id[df_volcano$group != 'NS'] #gene label

volcano_gaussian <- ggplot(data=df_volcano,mapping=aes(x=log2FC,y=-log10(affected.adjp),color=group,label=label))+
  geom_point()+
  #xlim(-0.228,0.228) +
  labs(title='volcano plot of significance against DEG recurrence and non-recurrence') +
  geom_label_repel(na.rm = TRUE, show.legend = T, box.padding = 0.5) +
  theme_minimal()+ 
  theme(legend.position="none")+
  geom_vline(xintercept=c(-1, 1), col="red", linetype='longdash') +
  geom_hline(yintercept=-log10(0.05), col="red", linetype='longdash') +
  geom_vline(xintercept = 0, col="grey", linetype='longdash') +
  labs(color = df_volcano$group) + 
  scale_color_manual(values=c("DOWN"= "#56B4E9", "NS"="#999999", "UP"="#E69F00"))

volcano_gaussian

#save the df_volcano - has the log2FC 
write.csv(df_volcano, file = "/Users/swashburn30/Desktop/isoform/validate_kates_results/prelim_DEG_results_2_15_2024.csv")

############PCA score#########
library(ggpubr)
library(glue)
library(dplyr)
library(ggpubr)
#stringent DEG filter: adjp < 0.05 and Log2FC > ± 0.25 - 762 DEGs
DEG <- pvals_gaussian[pvals_gaussian$affected.adjp<0.05 & (pvals_gaussian$log2FC>0.25| pvals_gaussian$log2FC< (-0.25)),]

#now test a more stringent filter of log2FC of > 1 - use the similar workflow below (just replaced DEG with DEG_filt in r console)
DEG_filt <- DEG[DEG$affected.adjp<0.05 & (DEG$log2FC>1| DEG$log2FC< (-1)),] #23 DEGs 


#select the sig up/down genes 
mat <- gene_trans%>% select(DEG$gene_id)
mat <- as.matrix(mat)
gene_pca <- prcomp(mat)
gene_pca_df <- as.data.frame(gene_pca$x)
#gene_pca <- gene_pca %>% select(PC1)

unique(row.names(gene_pca) == row.names(pheno)) #checking rna id matches

gene_pca_df$sample <- rownames(gene_pca_df)
pheno$sample <- rownames(pheno)
pc_score_metadata <- full_join(gene_pca_df, pheno, by = ("sample"))

#create PCA plot 
pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(affected))) +
  geom_point() 
#
pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(Batch))) +
  geom_point() 

pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(PhenoSex))) +
  geom_point() 

pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(race))) +
  geom_point() 



#variance explained by each PC
pc_eigenvalues <- gene_pca$sdev^2

pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))
pc_eigenvalues

#--------Boxplots - PC1---------#
p <- pc_score_metadata %>%
  ggplot( aes(x=affected, y=PC1, fill=affected)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("DEG PC1 score in all samples (n=334)") +
  xlab("Recurrence Status") +
  ylab("PC1") +
  scale_x_discrete(labels=c("R0" = "R0 (n=238)", "R1" = "R1 (n=96)"))

p + stat_compare_means(method = "wilcox.test", label.x = 2)

#wilcoxon test
new_df <- select(test_pc_score_metadata, "sample", "PC1", "tissue")

#calculate summary statistics
group_by(new_df, tissue) %>%
  summarise(
    count = n(),
    median = median(PC1, na.rm = TRUE),
    IQR = IQR(PC1, na.rm = TRUE)
  )


r0 <- subset(pc_score_metadata, affected == 'R0')
r1 <- subset(pc_score_metadata, affected == 'R1')

res <- wilcox.test(r0$PC1, r1$PC1)
res
#p-value = 8.148e-11

mean(r0$PC1) #-3.561835
mean(r1$PC1) #8.830382

#DEG PC1 score female 

#extract females from DF
female_df <- select(pc_score_metadata, "sample", "PC1", "PhenoSex", "affected", )
female_df <- female_df[female_df$PhenoSex == 'Female', ]

fp <- female_df %>%
  ggplot( aes(x=affected, y=PC1, fill=affected)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("DEG PC1 score in all female (n=164)") +
  xlab("Recurrence Status") +
  ylab("PC1") +
  scale_x_discrete(labels=c("R0" = "R0 (n=123)", "rectum" = "rectum (n=41)"))

fp + stat_compare_means(method = "wilcox.test", label.x = 2)


f_r0 <- subset(female_df, affected == 'R0')
f_r1 <- subset(female_df, affected == 'R1')

fres <- wilcox.test(f_r0$PC1, f_r1$PC1)
fres
#p-value = 0.000104

mean(f_r0$PC1) #-3.565999
mean(f_r1$PC1) #8.464961

#DEG PC1 score male 

#extract male from DF
male_df <- select(pc_score_metadata, "sample", "PC1", "PhenoSex", "affected", )
male_df <- male_df[male_df$PhenoSex == 'Male', ]

mp <- male_df %>%
  ggplot( aes(x=affected, y=PC1, fill=affected)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("DEG PC1 score in all male (n=170)") +
  xlab("Recurrence Status") +
  ylab("PC1") +
  scale_x_discrete(labels=c("R0" = "R0 (n=115)", "R1" = "R1 (n=55)"))

mp + stat_compare_means(method = "wilcox.test", label.x = 2)

m_r0 <- subset(male_df, affected == 'R0')
m_r1 <- subset(male_df, affected == 'R1')

mres <- wilcox.test(m_r0$PC1, m_r1$PC1)
mres
#p-value = 1.73e-07

mean(m_r0$PC1) #-3.557381
mean(m_r1$PC1) # 9.102787


#---------use stringent DEG threshold from rectum vs. ileum (validation) to compare r1 vs. r0-----------#

#load in pheno data 
pheno <- read_excel("/Users/swashburn30/Desktop/isoform/validate_kates_results/PhenoSummary.xlsx",sheet = 'Data')
pheno <- as.data.frame(pheno)
row.names(pheno)<-pheno$RNAid
pheno$Batch <- as.factor(pheno$Batch)
pheno <- pheno[,c("affected","PhenoSex", "Batch",'Smoking','race',"AgeAtSample","consortium_id")]
pheno <-na.omit(pheno)

#read in gene expression data 
gene_df <- read.table('/Users/swashburn30/Desktop/isoform/validate_kates_results/washburn_rsem.merged.gene_tpm.tsv', sep = '\t', header = TRUE)
row.names(gene_df) <- gene_df$gene_id
gene_df <- gene_df[, -1] #remove the rownames (it was just numbers)
gene_mapping <- gene_df %>% dplyr::select(transcript_id.s.,gene_id) 
gene_df <- gene_df%>% dplyr::select(row.names(pheno)) #select samples - 334 samples
gene_df <- na.omit(gene_df)  #remove NA (where the tpm for gene is 0) #29972 genes
#gene_df <- gene_df[ rowMeans(gene_df) > 5, ] #mean coverage > 5 reads - 5505 genes

gene_df <- gene_df[,row.names(pheno)] #reorder
#filter gene_df to contain only the top DEGs - did not filter out the mean expression (most likely won't be the top DEGs anyway)
gene_df <- gene_df[rownames(gene_df) %in% DEG_filt$gene_id, ]

gene_trans <- as.data.frame(t(gene_df)) #transpose

#scale the log2 normalize the gene expression level 
gene_trans <- log2(gene_trans + 1)


#------PCA-------#
mat <- gene_trans
mat <- as.matrix(mat)
gene_pca <- prcomp(mat)
gene_pca_df <- as.data.frame(gene_pca$x)
gene_pca_df$sample <- rownames(gene_pca_df)
#gene_pca <- gene_pca %>% select(PC1)

unique(row.names(gene_pca) == row.names(pheno)) #checking rna id matches

gene_pca_df$sample <- rownames(gene_pca_df)
pheno$sample <- rownames(pheno)
pc_score_metadata <- full_join(gene_pca_df, pheno, by = ("sample"))

#create PCA plot 
pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(affected))) +
  geom_point() 
#

pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(PhenoSex))) +
  geom_point() 

pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(race))) +
  geom_point() 

pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(Batch))) +
  geom_point() 

#pc plot of PC1 vs. PC4 
#create PCA plot 
pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC4, colour = factor(affected))) +
  xlab("PC1 (41.8% variance)") +
  ylab("PC4 (4.43% variance)") +
  geom_point() +
  guides(color = guide_legend(title = "Recurrence Status")) +
  theme(
    text = element_text(size = 20),
    
    
  )
#theme(text = element_text(size=20) 



#variance explained by each PC
pc_eigenvalues <- gene_pca$sdev^2

pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))
pc_eigenvalues

#which PCs are associated with recurrence status?

#pc1
pc1 <- lm(pc_score_metadata$PC1 ~ affected, data = pc_score_metadata)

#pc2
pc2 <- lm(pc_score_metadata$PC2 ~ affected, data = pc_score_metadata)

#pc3 
pc3 <- lm(pc_score_metadata$PC3 ~ affected, data = pc_score_metadata)

#pc4
pc4 <- lm(pc_score_metadata$PC4 ~ affected, data = pc_score_metadata)

#pc5
pc5 <- lm(pc_score_metadata$PC5 ~ affected, data = pc_score_metadata)

#pc6
pc6 <- lm(pc_score_metadata$PC6 ~ affected, data = pc_score_metadata)

#pc7
pc7 <- lm(pc_score_metadata$PC7 ~ affected, data = pc_score_metadata)

#pc8 
pc8 <- lm(pc_score_metadata$PC8 ~ affected, data = pc_score_metadata)

#pc9
pc9 <- lm(pc_score_metadata$PC9 ~ affected, data = pc_score_metadata)

#pc10
pc10 <- lm(pc_score_metadata$PC10 ~ affected, data = pc_score_metadata)



#--------Boxplots - PC1---------#
p1 <- pc_score_metadata %>%
  ggplot( aes(x=affected, y=PC1, fill=affected)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=20),
    text = element_text(size = 20)
  ) +
  ggtitle("DEG PC1 score in all samples (n=334)") +
  stat_compare_means(method = "wilcox.test", label.x = 2) +
  xlab("Recurrence Status") +
  ylab("PC1") +
  scale_x_discrete(labels=c("R0" = "R0 (n=238)", "R1" = "R1 (n=96)"))

p1 + stat_compare_means(method = "wilcox.test", label.x = 2)


r0 <- subset(pc_score_metadata, affected == 'R0')
r1 <- subset(pc_score_metadata, affected == 'R1')

res <- wilcox.test(r0$PC1, r1$PC1)
res
#p-value = 1.873e-10

mean(r0$PC1) #-2.80646
mean(r1$PC1) #6.957681

#DEG PC1 score female 

#extract females from DF
female_df <- select(pc_score_metadata, "sample", "PC1", "PhenoSex", "affected", )
female_df <- female_df[female_df$PhenoSex == 'Female', ]

fp1 <- female_df %>%
  ggplot( aes(x=affected, y=PC1, fill=affected)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("DEG PC1 score in all female (n=164)") +
  xlab("Recurrence Status") +
  ylab("PC1") +
  scale_x_discrete(labels=c("R0" = "R0 (n=123)", "R1" = "R1 (n=41)"))

fp1 + stat_compare_means(method = "wilcox.test", label.x = 2)


f_r0 <- subset(female_df, affected == 'R0')
f_r1 <- subset(female_df, affected == 'R1')

fres <- wilcox.test(f_r0$PC1, f_r1$PC1)
fres
#p-value = 0.0001626

mean(f_r0$PC1) #-2.778652
mean(f_r1$PC1) #6.89145

#DEG PC1 score male 

#extract male from DF
male_df <- select(pc_score_metadata, "sample", "PC1", "PhenoSex", "affected", )
male_df <- male_df[male_df$PhenoSex == 'Male', ]

mp1 <- male_df %>%
  ggplot( aes(x=affected, y=PC1, fill=affected)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("DEG PC1 score in all male (n=170)") +
  xlab("Recurrence Status") +
  ylab("PC1") +
  scale_x_discrete(labels=c("R0" = "R0 (n=115)", "R1" = "R1 (n=55)"))

mp1 + stat_compare_means(method = "wilcox.test", label.x = 2)

m_r0 <- subset(male_df, affected == 'R0')
m_r1 <- subset(male_df, affected == 'R1')

mres <- wilcox.test(m_r0$PC1, m_r1$PC1)
mres
#p-value = 2.346e-07

mean(m_r0$PC1) #-2.836202
mean(m_r1$PC1) #7.007053

#-------PC4 associated with recurrence status?--------#

p4 <- pc_score_metadata %>%
  ggplot( aes(x=affected, y=PC4, fill=affected)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("DEG PC4 score in all samples (n=334)") +
  xlab("Recurrence Status") +
  ylab("PC4") +
  scale_x_discrete(labels=c("R0" = "R0 (n=238)", "R1" = "R1 (n=96)"))

p4 + stat_compare_means(method = "wilcox.test", label.x = 2)


res4 <- wilcox.test(r0$PC4, r1$PC4)
res4
#p-value = 4.41e-09

mean(r0$PC4) #0.901859
mean(r1$PC4) #-2.235859

#extract females from DF
female_df4 <- select(pc_score_metadata, "sample", "PC4", "PhenoSex", "affected", )
female_df4 <- female_df4[female_df4$PhenoSex == 'Female', ]


fp4 <- female_df4 %>%
  ggplot( aes(x=affected, y=PC4, fill=affected)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("DEG PC4 score in all female (n=164)") +
  xlab("Recurrence Status") +
  ylab("PC4") +
  scale_x_discrete(labels=c("R0" = "R0 (n=123)", "R1" = "R1 (n=41)"))

fp4 + stat_compare_means(method = "wilcox.test", label.x = 2)


f_r04 <- subset(female_df4, affected == 'R0')
f_r14 <- subset(female_df4, affected == 'R1')

fres4 <- wilcox.test(f_r04$PC4, f_r14$PC4)
fres4
#p-value = 1.88e-06

mean(f_r04$PC4) #1.4226
mean(f_r14$PC4) #-2.388859

#extract male from DF
male_df4 <- select(pc_score_metadata, "sample", "PC4", "PhenoSex", "affected", )
male_df4 <- male_df4[male_df4$PhenoSex == 'Male', ]

mp4 <- male_df4 %>%
  ggplot( aes(x=affected, y=PC4, fill=affected)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("DEG PC4 score in all male (n=170)") +
  xlab("Recurrence Status") +
  ylab("PC4") +
  scale_x_discrete(labels=c("R0" = "R0 (n=115)", "R1" = "R1 (n=55)"))

mp4 + stat_compare_means(method = "wilcox.test", label.x = 2)

m_r04 <- subset(male_df4, affected == 'R0')
m_r14 <- subset(male_df4, affected == 'R1')

mres4 <- wilcox.test(m_r04$PC4, m_r14$PC4)
mres4
#p-value = 0.0004524

mean(m_r04$PC4) #0.3448928
mean(m_r14$PC4) #-2.121804

#----------use top 100 genes from PC1 to perform pathway enrichment analysis-------------#

DEG <- read_csv(file = "/Users/swashburn30/Desktop/isoform/kiera_analysis/deseq_deg_tissue_4_18_24.csv")

gene_df <- gene_df[,row.names(pheno)] #reorder

test_merge$gene_id <- rownames(test_merge)
#filter gene_df to contain only the top DEGs - did not filter out the mean expression (most likely won't be the top DEGs anyway)
gene_df <- gene_df[rownames(gene_df) %in% test_merge$gene_id, ]

gene_trans <- as.data.frame(t(gene_df)) #transpose
gene_trans <- "numeric"

#scale the log2 normalize the gene expression level 
gene_trans <- log2(gene_trans + 1.0)

#

#------PCA-------#
mat <- gene_trans
mat <- as.matrix(mat)
gene_pca <- prcomp(mat) #tested scale
gene_pca_df <- as.data.frame(gene_pca$x)
gene_pca_df$sample <- rownames(gene_pca_df)
#gene_pca <- gene_pca %>% select(PC1)

unique(row.names(gene_pca) == row.names(pheno)) #checking rna id matches

pheno$sample <- rownames(pheno)
pc_score_metadata <- full_join(gene_pca_df, pheno, by = ("sample"))

#Variance explained by PC
pc_eigenvalues <- gene_pca$sdev^2

pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))
pc_eigenvalues

#pc1
pc1 <- lm(pc_score_metadata$PC1 ~ affected, data = pc_score_metadata)

#pc2
pc2 <- lm(pc_score_metadata$PC2 ~ affected, data = pc_score_metadata)

#pc3 
pc3 <- lm(pc_score_metadata$PC3 ~ affected, data = pc_score_metadata)

#pc4
pc4 <- lm(pc_score_metadata$PC4 ~ affected, data = pc_score_metadata)

#pc5
pc5 <- lm(pc_score_metadata$PC5 ~ affected, data = pc_score_metadata)

#pc6
pc6 <- lm(pc_score_metadata$PC6 ~ affected, data = pc_score_metadata)

#pc7
pc7 <- lm(pc_score_metadata$PC7 ~ affected, data = pc_score_metadata)

#pc8 
pc8 <- lm(pc_score_metadata$PC8 ~ affected, data = pc_score_metadata)

#pc9
pc9 <- lm(pc_score_metadata$PC9 ~ affected, data = pc_score_metadata)

#pc10
pc10 <- lm(pc_score_metadata$PC10 ~ affected, data = pc_score_metadata)


#create PCA plot 
pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(affected))) +
  geom_point() 
#

pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC4, colour = factor(affected))) +
  geom_point() 

pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC3, colour = factor(affected))) +
  geom_point() 

#box plot of PC1 

pc1_violin <- pc_score_metadata %>%
  ggplot( aes(x=affected, y=PC1, fill=affected)) +
  geom_violin(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=10),
    text = element_text(size = 10)
  ) +
  ggtitle("DEG PC1 score in all samples (n=334)") +
  stat_compare_means(method = "wilcox.test", label.x.npc = "middle", size = 4, label = paste0("P" = "p.format")) +
  xlab("Recurrence Status") +
  ylab("PC1") +
  scale_x_discrete(labels=c("R0" = "R0 (n=238)", "R1" = "R1 (n=96)"))
pc1_violin + geom_boxplot(width=0.1)

pc3_violin <- pc_score_metadata %>%
  ggplot( aes(x=affected, y=PC3, fill=affected)) +
  geom_violin(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=10),
    text = element_text(size = 10)
  ) +
  ggtitle("DEG PC3 score in all samples (n=334)") +
  stat_compare_means(method = "wilcox.test", label.x.npc = "middle", size = 4, label = paste0("P" = "p.format")) +
  xlab("Recurrence Status") +
  ylab("PC3") +
  scale_x_discrete(labels=c("R0" = "R0 (n=238)", "R1" = "R1 (n=96)"))
pc3_violin + geom_boxplot(width=0.1)

#save PCA results 
write.csv(pc_score_metadata, file = "/Users/swashburn30/Desktop/isoform/validate_kates_results/pc_results_postop_tissueDEG_4_18_24.csv")

#extract top 100 genes from PCA 
gene_pca <- gene_pca$rotation
gene_pca <- as.data.frame(gene_pca)
gene_pca$gene <- rownames(gene_pca)

sorted_genes <- order(gene_pca$PC1, decreasing = TRUE)
top_100_genes <- rownames(gene_pca)[sorted_genes[1:100]]

top_gene_pc1 <- gene_pca %>% 
  select(gene, PC1) %>%
  # convert to a "long" format
  pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  # for each PC
  #group_by(PC) %>% 
  # arrange by descending order of loading
  arrange(desc(abs(loading))) %>%
  slice(1:100) %>% 
  # pull the gene column as a vector
  #pull(gene) %>% 
  # ensure only unique genes are retained
  unique()

#boferroni - more stringent threshold 
ego <- enrichGO(gene = top_gene_pc1$gene, 
                #universe = allOE_genes,
                keyType = "SYMBOL",
                OrgDb = "org.Hs.eg.db", 
                ont = "BP", 
                pAdjustMethod = "bonferroni", 
                qvalueCutoff = 0.05) 

#save as data table 
go_summary_pc1 <- data.frame(ego)

#visualize results
dotplot(ego)
barplot(ego)


#Figures for Manuscript 

#--------------A PC plot of PC1 vs. PC2 isoform fraction--------# 
#create PCA plot 
pdf("/Users/swashburn30/Desktop/isoform/manuscript_figures/preliminary_figures/PCA_p1_vs_p2_IF.pdf", width = 6.5, height = 6)
plot(pc_score_metadata %>% 
       # create the plot
       ggplot(aes(x = PC1, y = PC2, colour = factor(affected))) +
       xlab("PC1 (11.6% variance)") +
       ylab("PC2 (9.06% variance)") +
       geom_point() +
       guides(color = guide_legend(title = "Recurrence Status")) +
       theme(
         text = element_text(size = 10),
         
       ))
#theme(text = element_text(size=20) 
dev.off()

#--------------B Boxplot of PC2 score for recurrence status--------------#
pdf("/Users/swashburn30/Desktop/isoform/manuscript_figures/preliminary_figures/violinplot_pc2_tissue_predict_r0_vs_r1_IF.pdf", width = 6, height = 6)
pc2_violin <- pc_score_metadata %>%
  ggplot( aes(x=affected, y=PC2, fill=affected)) +
  geom_violin(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=10),
    text = element_text(size = 10)
  ) +
  ggtitle("DTU PC2 score in all samples (n=334)") +
  stat_compare_means(method = "wilcox.test", label.x.npc = "middle", size = 4, label = paste0("P" = "p.format")) +
  xlab("Recurrence Status") +
  ylab("PC2") +
  scale_x_discrete(labels=c("R0" = "R0 (n=238)", "R1" = "R1 (n=96)"))

plot(pc2_violin + geom_boxplot(width=0.1))
dev.off()

#-----------C pc plot of PC1 vs. PC4-------------#

#created using the DEGs comparing rectum vs. ileum 
#shows that PC1 adn PC4 are orthoganal 

#scaled version 
pdf("/Users/swashburn30/Desktop/isoform/manuscript_figures/preliminary_figures/PCA_p1_vs_p4_scaled.pdf", width = 6.5, height = 6)
plot(pc_score_metadata %>% 
       # create the plot
       ggplot(aes(x = PC1, y = PC4, colour = factor(affected))) +
       xlab("PC1 (30.4% variance)") +
       ylab("PC4 (4.05% variance)") +
       geom_point() +
       guides(color = guide_legend(title = "Recurrence Status")) +
       theme(
         text = element_text(size = 10),
         
       ))
#theme(text = element_text(size=20) 
dev.off()


#--------D Violin Plot of recurrence status predicted by tissue---------#
#scaled version 
pdf("/Users/swashburn30/Desktop/isoform/manuscript_figures/preliminary_figures/violinplot_pc1_tissue_predict_r0_vs_r1_scale.pdf", width = 6, height = 6)
test_violin <- pc_score_metadata %>%
  ggplot( aes(x=affected, y=PC1, fill=affected)) +
  geom_violin(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=10),
    text = element_text(size = 10)
  ) +
  ggtitle("DEG PC1 score in all samples (n=334)") +
  stat_compare_means(method = "wilcox.test", label.x.npc = "middle", size = 4, label = paste0("P" = "p.format")) +
  xlab("Recurrence Status") +
  ylab("PC1") +
  scale_x_discrete(labels=c("R0" = "R0 (n=238)", "R1" = "R1 (n=96)"))

plot(test_violin + geom_boxplot(width=0.1))
dev.off()

#------------E pathway enrichment in top 100 PC1 genes---------------#
#load in PC1 genes 

pc1_gene <- read.csv(file = "/Users/swashburn30/Desktop/isoform/validate_kates_results/pc1_genes.csv")
pc1_gene <- pc1_gene[, -1]
pc1_gene <- data.frame(pc1_gene)

library(DOSE)
library(pathview)
library(clusterProfiler)

#boferroni - more stringent threshold 
ego <- enrichGO(gene = pc1_gene$pc1_gene, 
                #universe = allOE_genes,
                keyType = "SYMBOL",
                OrgDb = "org.Hs.eg.db", 
                ont = "BP", 
                pAdjustMethod = "bonferroni", 
                qvalueCutoff = 0.05) 

#save as data table 
go_summary_pc1 <- data.frame(ego)

#save file 
write.csv(go_summary_pc1, file = "/Users/swashburn30/Desktop/isoform/enrich_go_pc1_gene_pathway.csv")

#visualize results
dotplot(ego)
barplot(ego)


#--------------F pathway enrichment from top 100 PC4 genes-------------#

#load in pc4 pathway results from ToppFun
pathways_pc4 <- read.table(file = "/Users/swashburn30/Desktop/isoform/PC4_topp_gene_pathway.txt", sep = '\t', header = TRUE)

#test making ggplot barplot 
pathways_pc4$q.value.Bonferroni <- as.numeric(pathways_pc4$q.value.Bonferroni)
pathways_pc4["8", "Name"] <- "REACTOME ANTIMICROBIAL PEPTIDES"
pathways_pc4$Name <- str_wrap(pathways_pc4$Name, width = 10)


#save this barplot - will have to edit the scale to flip in adobe 
pdf("/Users/swashburn30/Desktop/isoform/manuscript_figures/preliminary_figures/pc4_pathway_barplot.pdf", width = 7, height = 5.5)
plot(ggplot(pathways_pc4, aes(x = reorder(Name, +Hit.Count.in.Query.List), y = Hit.Count.in.Query.List, fill = q.value.Bonferroni)) +
       geom_bar(stat = 'identity') +
       #guides(colour = guide_legend(reverse=T)) +
       ylab("Count") +
       scale_fill_gradient(low = "red", high = "blue", name = "p.adjust") +
       #guides(colour = guide_legend(reverse=T)) +
       theme_bw(base_size = 15) +
       theme(axis.title.y=element_blank()) + 
       coord_flip())
dev.off()

#-----------------create supplemental figures for paper------------# 

#CBX3
pdf("/Users/swashburn30/Desktop/isoform/manuscript_figures/preliminary_figures/violinPlot_cbx3.pdf", width = 6, height = 6)
cbx3 <- gene_comb %>%
  ggplot( aes(x=affected, y=CBX3, fill=affected)) +
  geom_violin(alpha=0.3) + 
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=10)
  ) +
  stat_compare_means(method = "wilcox.test", label.x.npc = "middle", size = 4, label = paste0("P" = "p.format")) +
  #ggtitle("expression of gene") +
  xlab("Recurrence Status") +
  ylab("CBX3") +
  scale_x_discrete(labels=c("R0" = "R0 (n= 241)", "R1" = "R1 (n=96)")
  )
plot(cbx3 + geom_boxplot(width=0.1))
dev.off()

#CBX5
pdf("/Users/swashburn30/Desktop/isoform/manuscript_figures/preliminary_figures/violinPlot_cbx5.pdf", width = 6, height = 6)
cbx5 <- gene_comb %>%
  ggplot( aes(x=affected, y=CBX5, fill=affected)) +
  geom_violin(alpha=0.3) + 
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=10)
  ) +
  stat_compare_means(method = "wilcox.test", label.x.npc = "middle", size = 4, label = paste0("P" = "p.format")) +
  #ggtitle("expression of gene") +
  xlab("Recurrence Status") +
  ylab("CBX5") +
  scale_x_discrete(labels=c("R0" = "R0 (n=238)", "R1" = "R1 (n=96)")
  )
plot(cbx5 + geom_boxplot(width=0.1))
dev.off()

#MKI67
pdf("/Users/swashburn30/Desktop/isoform/manuscript_figures/preliminary_figures/violinPlot_mki67.pdf", width = 6, height = 6)
mki67 <- gene_comb %>%
  ggplot( aes(x=affected, y=MKI67, fill=affected)) +
  geom_violin(alpha=0.3) + 
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=10)
  ) +
  stat_compare_means(method = "wilcox.test", label.x.npc = "middle", size = 4, label = paste0("P" = "p.format")) +
  #ggtitle("expression of gene") +
  xlab("Recurrence Status") +
  ylab("MKI67") +
  scale_x_discrete(labels=c("R0" = "R0 (n=238)", "R1" = "R1 (n=96)")
  )
plot(mki67 + geom_boxplot(width=0.1))
dev.off()

#GSEA 
#+ is up in r0 and - is up in r1
DEG <- read.csv(file = "/Users/swashburn30/Desktop/isoform/validate_kates_results/recur0vs1.csv")
DEG$padj <- as.numeric(DEG$padj)
DEG <- DEG %>%
  filter(padj <= 0.05)

DEG$group <- 'NS'
DEG$group[DEG$log2FoldChange< (-1) & DEG$padj<0.05] <-'DOWN'
DEG$group[DEG$log2FoldChange > (1) & DEG$padj<0.05] <- 'UP'
#df_volcano$gene_transcript <- paste0(df_volcano$gene_id,': ',df_volcano$transcript_id)
#df_volcano$label[df_volcano$group !='NS'] <- df_volcano$gene_id[df_volcano$group != 'NS'] #gene label
DEG$label <- ifelse(DEG$group == "NS", NA, DEG$group)
DEG$label[DEG$group !='NS'] <- DEG$X[DEG$group != 'NS']
rank <- DEG$log2FoldChange
names(rank) <- DEG$X
barplot(sort(rank, decreasing = T))
rank = sort(rank, decreasing = TRUE)


res2 <- DEG %>% 
  dplyr::select(X, log2FoldChange) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(X) %>% 
  summarize(stat=mean(log2FoldChange))
res2

#remove the weird gene numbers 
res2 <- res2[-c(1:7), ]
ranks <- deframe(res2)


#kegg pathways 
BiocManager::install("KEGGREST")
library(KEGGREST)
library(org.Hs.eg.db)
library(fgsea)
data(KEGGREST)



go <- gmtPathways("/Users/swashburn30/Downloads/c5.go.v2023.2.Hs.symbols.gmt")

go_res <- fgsea(pathways = go, 
                stats    = ranks,
                minSize  = 15,
                maxSize  = 500)

go_res <- go_res[go_res$padj < 0.05, ]

gsea_plot <- plotEnrichment(go[["GOBP_REGULATION_OF_CELL_CYCLE_PROCESS"]],
                            ranks) + labs(title="Regulation of Cell Cycle Process") #adj p 0.02


pdf("/Users/swashburn30/Desktop/isoform/manuscript_figures/preliminary_figures/gsea_cellcycle.pdf", width = 7.5, height = 6)
plot(gsea_plot)
dev.off()

#---------------p53 investigation------------#

library(readxl)
pheno_postop <- read_excel("/Users/swashburn30/Desktop/isoform/validate_kates_results/PhenoSummary.xlsx",sheet = 'Data')
pheno_postop <- as.data.frame(pheno_postop)
rownames(pheno_postop) <- make.names(pheno_postop$RNAid, unique = TRUE)
#row.names(pheno)<-pheno$RNAid
pheno_postop$Batch <- as.factor(pheno_postop$Batch)
pheno_postop <- pheno_postop[,c("affected","PhenoSex", "Batch",'Smoking','race',"AgeAtSample","consortium_id")]
pheno_postop <-na.omit(pheno_postop) #334 samples instead of 335


#----------read in isoform fraction-----------#


IF_df <- read.table('/Users/swashburn30/Desktop/isoform/validate_kates_results/washburn_isoform_fraction.tsv',sep='\t', header=TRUE)
row.names(IF_df) <- IF_df$transcript_id
#gene_mapping <- IF_df %>% dplyr::select(transcript_id,gene_id)
IF_df <- IF_df%>% dplyr::select(row.names(pheno_postop)) #select samples - now there are the same in pheno and iso table
IF_df <- na.omit(IF_df)  #remove NA (where the tpm for gene is 0)
#for this part - did not filter out the IFs (just filter based on transcript_id in tissue_DTU)

IF_df <- IF_df[,row.names(pheno_postop)] #reorder - pheno and if_df samples are in same order

gene_mapping_filt <- subset(gene_mapping, transcript_id %in% rownames(IF_df))

all(rownames(IF_df) == gene_mapping_filt$transcript_id)

IF_df$gene_ID <- gene_mapping_filt$gene_id

#subset IF_df/IF_trans based on p53 genes 
IF_df_hall <- IF_df[IF_df$gene_ID %in% p53_gene_hall[[1]], ] #30 transcripts

#genes in IF_df_hall
gene_hall_if_df <- IF_df_hall$gene_ID

gene_hall_if_df <- as.data.frame(gene_hall_if_df)

rownames(gene_hall_if_df) <- rownames(IF_df_hall)

#remove gene IDs from IF df hall 
IF_df_hall$gene_ID <- NULL

IF_trans_hall <- as.data.frame(t(IF_df_hall)) #transpose

#perform PCA using the tissue DTUs 
df_pca <- prcomp(IF_trans_hall)

all(row.names(df_pca) == row.names(pheno_postop)) #check that samples match ##TRUE

#variance explained by each PC
pc_eigenvalues <- df_pca$sdev^2

pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_eigenvalues #PC1: 55.4% of variance; PC2: 14.4% of variance 

# The PC scores are stored in the "x" value of the prcomp object
pc_scores <- df_pca$x

pc_scores <- pc_scores %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")

pc_scores

#add sample name column 
pheno_postop$sample <- rownames(pheno_postop)
pc_score_metadata <- full_join(pc_scores, pheno_postop, by = ("sample"))

#create PCA plot 
pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(affected))) +
  geom_point() 

#create PCA plot - PC1 is just batch 
pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(Batch))) +
  geom_point() 


#save the PC scores with metadata
write.csv(pc_score_metadata, file = "/Users/swashburn30/Desktop/isoform/validate_kates_results/P53_ileum_vs_rectum_enrichment_genes_r0_vs_r1_11_12_24.csv")

#which PC is associated with recurrence status?
pc1 <- lm(pc_score_metadata$PC1 ~ affected, data = pc_score_metadata)
pc2 <- lm(pc_score_metadata$PC2 ~ affected, data = pc_score_metadata)
pc3 <- lm(pc_score_metadata$PC3 ~ affected, data = pc_score_metadata)
pc4 <- lm(pc_score_metadata$PC4 ~ affected, data = pc_score_metadata)
pc5 <- lm(pc_score_metadata$PC5 ~ affected, data = pc_score_metadata)
pc6 <- lm(pc_score_metadata$PC6 ~ affected, data = pc_score_metadata)
pc7 <- lm(pc_score_metadata$PC7 ~ affected, data = pc_score_metadata)
pc8 <- lm(pc_score_metadata$PC8 ~ affected, data = pc_score_metadata)
pc9 <- lm(pc_score_metadata$PC9 ~ affected, data = pc_score_metadata)
pc10 <- lm(pc_score_metadata$PC10 ~ affected, data = pc_score_metadata)

#PC6 associated with recurrence status

#Create Boxplot of PC6 score 
p1 <- pc_score_metadata %>%
  ggplot( aes(x=affected, y=PC6, fill=affected)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=20),
    text = element_text(size = 20)
  ) +
  ggtitle("PC6 score (n=334)") +
  stat_compare_means(method = "wilcox.test", label.x = 2) +
  xlab("Recurrence Status") +
  ylab("PC6 score") +
  scale_x_discrete(labels=c("R0" = "R0 (n=238)", "R1" = "R1 (n=96)"))

p1 + stat_compare_means(method = "wilcox.test", label.x = 2)


#Create Boxplot of PC7 score 
p2 <- pc_score_metadata %>%
  ggplot( aes(x=affected, y=PC7, fill=affected)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=20),
    text = element_text(size = 20)
  ) +
  ggtitle("PC7 score (n=334)") +
  stat_compare_means(method = "wilcox.test", label.x = 2) +
  xlab("Recurrence Status") +
  ylab("PC7 score") +
  scale_x_discrete(labels=c("R0" = "R0 (n=238)", "R1" = "R1 (n=96)"))

p2 + stat_compare_means(method = "wilcox.test", label.x = 2)

#---------p53 oncogene pathway----------#

IF_df <- read.table('/Users/swashburn30/Desktop/isoform/validate_kates_results/washburn_isoform_fraction.tsv',sep='\t', header=TRUE)
row.names(IF_df) <- IF_df$transcript_id
#gene_mapping <- IF_df %>% dplyr::select(transcript_id,gene_id)
IF_df <- IF_df%>% dplyr::select(row.names(pheno_postop)) #select samples - now there are the same in pheno and iso table
IF_df <- na.omit(IF_df)  #remove NA (where the tpm for gene is 0)
#for this part - did not filter out the IFs (just filter based on transcript_id in tissue_DTU)

IF_df <- IF_df[,row.names(pheno_postop)] #reorder - pheno and if_df samples are in same order

gene_mapping_filt <- subset(gene_mapping, transcript_id %in% rownames(IF_df))

all(rownames(IF_df) == gene_mapping_filt$transcript_id)

IF_df$gene_ID <- gene_mapping_filt$gene_id

#subset IF_df/IF_trans based on p53 genes 
IF_df_onco <- IF_df[IF_df$gene_ID %in% p53_gene_onco[[1]], ] #69 transcripts

#genes in IF_df_hall
gene_onco_if_df <- IF_df_onco$gene_ID

gene_onco_if_df <- as.data.frame(gene_onco_if_df)

rownames(gene_onco_if_df) <- rownames(IF_df_onco)

#remove gene IDs from IF df hall 
IF_df_onco$gene_ID <- NULL

IF_trans_onco <- as.data.frame(t(IF_df_onco)) #transpose

#perform PCA using the tissue DTUs 
df_pca_onco <- prcomp(IF_trans_onco)

all(row.names(df_pca_onco) == row.names(pheno_postop)) #check that samples match ##TRUE

#variance explained by each PC
pc_eigenvalues_onco <- df_pca_onco$sdev^2

pc_eigenvalues_onco <- tibble(PC = factor(1:length(pc_eigenvalues_onco)), 
                              variance = pc_eigenvalues_onco) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_eigenvalues_onco #PC1: 24% of variance; PC2: 17.8% of variance 

# The PC scores are stored in the "x" value of the prcomp object
pc_scores_onco <- df_pca_onco$x

pc_scores_onco <- pc_scores_onco %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")

pc_scores_onco

#add sample name column 
pheno_postop$sample <- rownames(pheno_postop)
pc_score_metadata_onco <- full_join(pc_scores_onco, pheno_postop, by = ("sample"))

#create PCA plot 
pc_score_metadata_onco %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(affected))) +
  geom_point() 

#create PCA plot - PC1 is just batch 
pc_score_metadata_onco %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(Batch))) +
  geom_point() 


#save the PC scores with metadata
write.csv(pc_score_metadata_onco, file = "/Users/swashburn30/Desktop/isoform/validate_kates_results/P53_ileum_vs_rectum_enrichment_onco_genes_r0_vs_r1_11_12_24.csv")

#which PC is associated with recurrence status?
pc1 <- lm(pc_score_metadata_onco$PC1 ~ affected, data = pc_score_metadata_onco)
pc2 <- lm(pc_score_metadata_onco$PC2 ~ affected, data = pc_score_metadata_onco)
pc3 <- lm(pc_score_metadata_onco$PC3 ~ affected, data = pc_score_metadata_onco)
pc4 <- lm(pc_score_metadata_onco$PC4 ~ affected, data = pc_score_metadata_onco)
pc5 <- lm(pc_score_metadata_onco$PC5 ~ affected, data = pc_score_metadata_onco)
pc6 <- lm(pc_score_metadata_onco$PC6 ~ affected, data = pc_score_metadata_onco)
pc7 <- lm(pc_score_metadata_onco$PC7 ~ affected, data = pc_score_metadata_onco)
pc8 <- lm(pc_score_metadata_onco$PC8 ~ affected, data = pc_score_metadata_onco)
pc9 <- lm(pc_score_metadata_onco$PC9 ~ affected, data = pc_score_metadata_onco)
pc10 <- lm(pc_score_metadata_onco$PC10 ~ affected, data = pc_score_metadata_onco)

#PC4 PC7 and PC10 associated with recurrence status

#Create Boxplot of PC6 score 
p1 <- pc_score_metadata_onco %>%
  ggplot( aes(x=affected, y=PC4, fill=affected)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=20),
    text = element_text(size = 20)
  ) +
  ggtitle("PC4 score (n=334)") +
  stat_compare_means(method = "wilcox.test", label.x = 2) +
  xlab("Recurrence Status") +
  ylab("PC4 score") +
  scale_x_discrete(labels=c("R0" = "R0 (n=238)", "R1" = "R1 (n=96)"))

p1 + stat_compare_means(method = "wilcox.test", label.x = 2)


#Create Boxplot of PC7 score 
p2 <- pc_score_metadata_onco %>%
  ggplot( aes(x=affected, y=PC7, fill=affected)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=20),
    text = element_text(size = 20)
  ) +
  ggtitle("PC7 score (n=334)") +
  stat_compare_means(method = "wilcox.test", label.x = 2) +
  xlab("Recurrence Status") +
  ylab("PC7 score") +
  scale_x_discrete(labels=c("R0" = "R0 (n=238)", "R1" = "R1 (n=96)"))

p2 + stat_compare_means(method = "wilcox.test", label.x = 2)

p3 <- pc_score_metadata_onco %>%
  ggplot( aes(x=affected, y=PC10, fill=affected)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=20),
    text = element_text(size = 20)
  ) +
  ggtitle("PC10 score (n=334)") +
  stat_compare_means(method = "wilcox.test", label.x = 2) +
  xlab("Recurrence Status") +
  ylab("PC10 score") +
  scale_x_discrete(labels=c("R0" = "R0 (n=238)", "R1" = "R1 (n=96)"))

p3 + stat_compare_means(method = "wilcox.test", label.x = 2)

#----------test on gene level-----------#

gene_df_postop <- read.table('/Users/swashburn30/Desktop/isoform/validate_kates_results/washburn_rsem.merged.gene_tpm.tsv', sep = '\t', header = TRUE)
row.names(gene_df_postop) <- gene_df_postop$gene_id
gene_df_postop <- gene_df_postop[, -1] #remove the rownames (it was just numbers)
#gene_mapping <- gene_df %>% dplyr::select(transcript_id.s.,gene_id) 
gene_df_postop <- gene_df_postop%>% dplyr::select(row.names(pheno_postop)) #select samples - 334 samples
gene_df_postop <- na.omit(gene_df_postop)  #remove NA (where the tpm for gene is 0) #29972 genes
#gene_df <- gene_df[ rowMeans(gene_df) > 5, ] #mean coverage > 5 reads - 5505 genes

gene_df_postop <- gene_df_postop[,row.names(pheno_postop)] #reorder
#filter gene_df to contain only the top DEGs - did not filter out the mean expression (most likely won't be the top DEGs anyway)
gene_postop_hall <- gene_df_postop[rownames(gene_df_postop) %in% p53_gene_hall[[1]], ]

gene_postop_hall_trans <- as.data.frame(t(gene_postop_hall)) #transpose

#scale the log2 normalize the gene expression level 
gene_postop_hall_trans <- log2(gene_postop_hall_trans + 1)

#perform PCA using the tissue DTUs 
df_pca_gene_hall <- prcomp(gene_postop_hall_trans)

all(row.names(df_pca_gene_hall) == row.names(pheno_postop)) #check that samples match ##TRUE

#variance explained by each PC
pc_eigenvalues_gene_hall <- df_pca_gene_hall$sdev^2

pc_eigenvalues_gene_hall <- tibble(PC = factor(1:length(pc_eigenvalues_gene_hall)), 
                                   variance = pc_eigenvalues_gene_hall) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_eigenvalues_gene_hall #PC1: 40.1% of variance; PC2: 20.8% of variance 

# The PC scores are stored in the "x" value of the prcomp object
pc_scores_gene_hall <- df_pca_gene_hall$x

pc_scores_gene_hall <- pc_scores_gene_hall %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")

pc_scores_gene_hall

#add sample name column 
pheno_postop$sample <- rownames(pheno_postop)
pc_score_metadata_gene_hall <- full_join(pc_scores_gene_hall, pheno_postop, by = ("sample"))

#create PCA plot 
pc_score_metadata_gene_hall %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(affected))) +
  geom_point() 

#create PCA plot - PC1 is just batch 
pc_score_metadata_gene_hall %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(Batch))) +
  geom_point() 


#save the PC scores with metadata
write.csv(pc_score_metadata_gene_hall, file = "/Users/swashburn30/Desktop/isoform/validate_kates_results/P53_ileum_vs_rectum_enrichment_hallmark_genes_r0_vs_r1_gene_level_11_12_24.csv")

#which PC is associated with recurrence status?
pc1 <- lm(pc_score_metadata_gene_hall$PC1 ~ affected, data = pc_score_metadata_gene_hall)
pc2 <- lm(pc_score_metadata_gene_hall$PC2 ~ affected, data = pc_score_metadata_gene_hall)
pc3 <- lm(pc_score_metadata_gene_hall$PC3 ~ affected, data = pc_score_metadata_gene_hall)
pc4 <- lm(pc_score_metadata_gene_hall$PC4 ~ affected, data = pc_score_metadata_gene_hall)
pc5 <- lm(pc_score_metadata_gene_hall$PC5 ~ affected, data = pc_score_metadata_gene_hall)
pc6 <- lm(pc_score_metadata_gene_hall$PC6 ~ affected, data = pc_score_metadata_gene_hall)
pc7 <- lm(pc_score_metadata_gene_hall$PC7 ~ affected, data = pc_score_metadata_gene_hall)
pc8 <- lm(pc_score_metadata_gene_hall$PC8 ~ affected, data = pc_score_metadata_gene_hall)
pc9 <- lm(pc_score_metadata_gene_hall$PC9 ~ affected, data = pc_score_metadata_gene_hall)
pc10 <- lm(pc_score_metadata_gene_hall$PC10 ~ affected, data = pc_score_metadata_gene_hall)

#PC2, PC3, PC9, and PC10 associated with recurrence status

#Create Boxplot of PC2 score 
p1 <- pc_score_metadata_gene_hall %>%
  ggplot( aes(x=affected, y=PC2, fill=affected)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=20),
    text = element_text(size = 20)
  ) +
  ggtitle("PC2 score (n=334)") +
  stat_compare_means(method = "wilcox.test", label.x = 2) +
  xlab("Recurrence Status") +
  ylab("PC2 score") +
  scale_x_discrete(labels=c("R0" = "R0 (n=238)", "R1" = "R1 (n=96)"))

p1 + stat_compare_means(method = "wilcox.test", label.x = 2)


#Create Boxplot of PC3 score 
p2 <- pc_score_metadata_gene_hall %>%
  ggplot( aes(x=affected, y=PC3, fill=affected)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=20),
    text = element_text(size = 20)
  ) +
  ggtitle("PC3 score (n=334)") +
  stat_compare_means(method = "wilcox.test", label.x = 2) +
  xlab("Recurrence Status") +
  ylab("PC3 score") +
  scale_x_discrete(labels=c("R0" = "R0 (n=238)", "R1" = "R1 (n=96)"))

p2 + stat_compare_means(method = "wilcox.test", label.x = 2)

p3 <- pc_score_metadata_gene_hall %>%
  ggplot( aes(x=affected, y=PC9, fill=affected)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=20),
    text = element_text(size = 20)
  ) +
  ggtitle("PC9 score (n=334)") +
  stat_compare_means(method = "wilcox.test", label.x = 2) +
  xlab("Recurrence Status") +
  ylab("PC9 score") +
  scale_x_discrete(labels=c("R0" = "R0 (n=238)", "R1" = "R1 (n=96)"))

p3 + stat_compare_means(method = "wilcox.test", label.x = 2)

p4 <- pc_score_metadata_gene_hall %>%
  ggplot( aes(x=affected, y=PC10, fill=affected)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=20),
    text = element_text(size = 20)
  ) +
  ggtitle("PC10 score (n=334)") +
  stat_compare_means(method = "wilcox.test", label.x = 2) +
  xlab("Recurrence Status") +
  ylab("PC10 score") +
  scale_x_discrete(labels=c("R0" = "R0 (n=238)", "R1" = "R1 (n=96)"))

p4 + stat_compare_means(method = "wilcox.test", label.x = 2)

#--------oncogene pathway--------#
gene_postop_onco <- gene_df_postop[rownames(gene_df_postop) %in% p53_gene_onco[[1]], ]

gene_postop_onco_trans <- as.data.frame(t(gene_postop_onco)) #transpose

#scale the log2 normalize the gene expression level 
gene_postop_onco_trans <- log2(gene_postop_onco_trans + 1)

#perform PCA using the tissue DTUs 
df_pca_gene_onco <- prcomp(gene_postop_onco_trans)

all(row.names(df_pca_gene_onco) == row.names(pheno_postop)) #check that samples match ##TRUE

#variance explained by each PC
pc_eigenvalues_gene_onco <- df_pca_gene_onco$sdev^2

pc_eigenvalues_gene_onco <- tibble(PC = factor(1:length(pc_eigenvalues_gene_onco)), 
                                   variance = pc_eigenvalues_gene_onco) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_eigenvalues_gene_onco #PC1: 44.8% of variance; PC2: 14.5% of variance 

# The PC scores are stored in the "x" value of the prcomp object
pc_scores_gene_onco <- df_pca_gene_onco$x

pc_scores_gene_onco <- pc_scores_gene_onco %>% 
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample")

pc_scores_gene_onco

#add sample name column 
pheno_postop$sample <- rownames(pheno_postop)
pc_score_metadata_gene_onco <- full_join(pc_scores_gene_onco, pheno_postop, by = ("sample"))

#create PCA plot 
pc_score_metadata_gene_onco %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(affected))) +
  geom_point() 

#create PCA plot - PC1 is just batch 
pc_score_metadata_gene_onco %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(Batch))) +
  geom_point() 


#save the PC scores with metadata
write.csv(pc_score_metadata_gene_onco, file = "/Users/swashburn30/Desktop/isoform/validate_kates_results/P53_ileum_vs_rectum_enrichment_onco_genes_r0_vs_r1_gene_level_11_12_24.csv")

#which PC is associated with recurrence status?
pc1 <- lm(pc_score_metadata_gene_onco$PC1 ~ affected, data = pc_score_metadata_gene_onco)
pc2 <- lm(pc_score_metadata_gene_onco$PC2 ~ affected, data = pc_score_metadata_gene_onco)
pc3 <- lm(pc_score_metadata_gene_onco$PC3 ~ affected, data = pc_score_metadata_gene_onco)
pc4 <- lm(pc_score_metadata_gene_onco$PC4 ~ affected, data = pc_score_metadata_gene_onco)
pc5 <- lm(pc_score_metadata_gene_onco$PC5 ~ affected, data = pc_score_metadata_gene_onco)
pc6 <- lm(pc_score_metadata_gene_onco$PC6 ~ affected, data = pc_score_metadata_gene_onco)
pc7 <- lm(pc_score_metadata_gene_onco$PC7 ~ affected, data = pc_score_metadata_gene_onco)
pc8 <- lm(pc_score_metadata_gene_onco$PC8 ~ affected, data = pc_score_metadata_gene_onco)
pc9 <- lm(pc_score_metadata_gene_onco$PC9 ~ affected, data = pc_score_metadata_gene_onco)
pc10 <- lm(pc_score_metadata_gene_onco$PC10 ~ affected, data = pc_score_metadata_gene_onco)

#PC1, PC2, PC3, PC7, and PC10 associated with recurrence status

#Create Boxplot of PC2 score 
p1 <- pc_score_metadata_gene_onco %>%
  ggplot( aes(x=affected, y=PC1, fill=affected)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=20),
    text = element_text(size = 20)
  ) +
  ggtitle("PC1 score (n=334)") +
  stat_compare_means(method = "wilcox.test", label.x = 2) +
  xlab("Recurrence Status") +
  ylab("PC1 score") +
  scale_x_discrete(labels=c("R0" = "R0 (n=238)", "R1" = "R1 (n=96)"))

p1 + stat_compare_means(method = "wilcox.test", label.x = 2)


#Create Boxplot of PC3 score 
p2 <- pc_score_metadata_gene_onco %>%
  ggplot( aes(x=affected, y=PC2, fill=affected)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=20),
    text = element_text(size = 20)
  ) +
  ggtitle("PC2 score (n=334)") +
  stat_compare_means(method = "wilcox.test", label.x = 2) +
  xlab("Recurrence Status") +
  ylab("PC2 score") +
  scale_x_discrete(labels=c("R0" = "R0 (n=238)", "R1" = "R1 (n=96)"))

p2 + stat_compare_means(method = "wilcox.test", label.x = 2)

p3 <- pc_score_metadata_gene_onco %>%
  ggplot( aes(x=affected, y=PC3, fill=affected)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=20),
    text = element_text(size = 20)
  ) +
  ggtitle("PC3 score (n=334)") +
  stat_compare_means(method = "wilcox.test", label.x = 2) +
  xlab("Recurrence Status") +
  ylab("PC3 score") +
  scale_x_discrete(labels=c("R0" = "R0 (n=238)", "R1" = "R1 (n=96)"))

p3 + stat_compare_means(method = "wilcox.test", label.x = 2)

p4 <- pc_score_metadata_gene_onco %>%
  ggplot( aes(x=affected, y=PC7, fill=affected)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=20),
    text = element_text(size = 20)
  ) +
  ggtitle("PC7 score (n=334)") +
  stat_compare_means(method = "wilcox.test", label.x = 2) +
  xlab("Recurrence Status") +
  ylab("PC7 score") +
  scale_x_discrete(labels=c("R0" = "R0 (n=238)", "R1" = "R1 (n=96)"))

p4 + stat_compare_means(method = "wilcox.test", label.x = 2)

p5 <- pc_score_metadata_gene_onco %>%
  ggplot( aes(x=affected, y=PC10, fill=affected)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=20),
    text = element_text(size = 20)
  ) +
  ggtitle("PC10 score (n=334)") +
  stat_compare_means(method = "wilcox.test", label.x = 2) +
  xlab("Recurrence Status") +
  ylab("PC10 score") +
  scale_x_discrete(labels=c("R0" = "R0 (n=238)", "R1" = "R1 (n=96)"))

p5 + stat_compare_means(method = "wilcox.test", label.x = 2)
