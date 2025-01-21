#replication from GEO GSE158952

#load libraries 
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

#load in pheno 
pheno <- read_csv("/Users/swashburn30/Desktop/isoform/kiera_analysis/phenodata_2_9_24.csv", col_names = F)
pheno <- as.data.frame(pheno)
colnames(pheno) <- c("donor", 'geo_id', 'sex', 'tissue', 'sample', 'age', 'race')
row.names(pheno)<-pheno$sample
pheno <-na.omit(pheno)
pheno$donor_label <- c("donor1", "donor1", "donor2", "donor2", "donor3", "donor3", "donor4", "donor4", 'donor5', "donor5", "donor6", "donor6", "donor7", "donor7", "donor8", "donor8", "donor9", "donor9", "donor10", "donor10")

#load in IF  - calculated the same as post-op ileal dataset
#read in the isoform fraction tsv file 
IF_df <- read_tsv("/Users/swashburn30/Desktop/isoform/kiera_analysis/isoform_fraction_2_5_24.tsv")
IF_df <- as.data.frame(IF_df)

row.names(IF_df) <- IF_df$transcript_id
gene_mapping <- IF_df %>% dplyr::select(transcript_id,gene_id)
#IF_df <- IF_df%>% dplyr::select(row.names(pheno))#select samples - now there are the same in pheno and iso table
IF_df <- IF_df[, -c(1,2) , drop = FALSE]
IF_df <- na.omit(IF_df)  #remove NA (where the tpm for gene is 0)
IF_df <- IF_df[rowSums(IF_df != 0) > (0.05*ncol(IF_df)), ] #keep transcript only no more than 95% ids have non-zero IF
IF_df <- IF_df[rowSums(IF_df < 1) > (0.95*ncol(IF_df)), ] #keep transcript only 95% individuals have it IF < 1


IF_df <- IF_df[,row.names(pheno)] #reorder - pheno and if_df samples are in same order
IF_trans <- as.data.frame(t(IF_df)) #transpose


#combine pheno and IF_df tables - test for just one IF 
IF_design <- cbind(pheno[,c("donor","geo_id", "sex",'tissue','sample',"age","race", "donor_label")],IF_trans[,1])
names(IF_design) <- c("donor","geo_id", "sex",'tissue','sample',"age","race", "donor_label","IF")


IF_design <- na.omit(IF_design) #remove na values

############## linear mixed model, gaussian distribution#######
model.full = lmer(IF ~ tissue + sex + race + (1|donor_label),  data=IF_design ,REML = FALSE) #full model

model.tissue = lmer(IF ~ sex + race + (1|donor_label),  data=IF_design ,REML = FALSE)

#extract pvalue using maximum likelyhood
prob_tissue <- anova(model.full,model.tissue)[2,8]
pvals <- matrix(c(prob_tissue),nrow=1,ncol=1)

#loop through the rest of isoforms
for (m in 2:ncol(IF_trans)){
  IF_design <- cbind(pheno[,c("donor","geo_id", "sex",'tissue','sample',"age","race", "donor_label")],IF_trans[,m])
  names(IF_design) <- c("donor","geo_id", "sex",'tissue','sample',"age","race", "donor_label","IF")
  
  IF_design <- na.omit(IF_design) #remove na values
  
  model.full = lmer(IF ~ tissue + sex + race + (1|donor_label),  data=IF_design ,REML = FALSE) #full model
  
  model.tissue = lmer(IF ~ sex + race + (1|donor_label),  data=IF_design ,REML = FALSE)
  
  
  #extract pvalue using maximum likelyhood
  
  prob_tissue <- anova(model.full,model.tissue)[2,8]
  pvals <- rbind(pvals, c(prob_tissue))
  
}

#change the rownames 
row.names(pvals) <- row.names(IF_df)
colnames(pvals) <- c('tissue.p')

#adjust p value - calcualte 
pvals_gaussian <- as.data.frame(pvals)
pvals_gaussian$tissue.adjp <- p.adjust(pvals_gaussian[,1], method='fdr')
pvals_gaussian$transcript_id<- row.names(pvals_gaussian)
pvals_gaussian<- left_join(pvals_gaussian,gene_mapping,by='transcript_id')


###########calculate dIF: mean (IF(recurrence)) - mean (IF(non-recurrence))##########
IF_trans$tissue <- pheno$tissue
#mean of just the first isoform, for each group
mean_IF <- aggregate(IF_trans[,1],list(IF_trans$tissue), FUN=mean)
dIF_each <- mean_IF[mean_IF$Group.1=='rectum','x'] - mean_IF[mean_IF$Group.1=='ileum','x']
dIF <- c(names(IF_trans)[1],dIF_each)

#i dont think this worked
colnames(dIF) <- c('transcript_id','dIF')

for (m in 2:(ncol(IF_trans)-1)){
  mean_IF <- aggregate(IF_trans[,m],list(IF_trans$tissue), FUN=mean)
  dIF_each <- mean_IF[mean_IF$Group.1=='rectum','x'] - mean_IF[mean_IF$Group.1=='ileum','x']
  dIF <- rbind(dIF,c(names(IF_trans)[m],dIF_each))
  
}
dIF <- as.data.frame(dIF)
colnames(dIF) <- c('transcript_id','dIF')
row.names(dIF) <- dIF$transcript_id
dIF$dIF <- as.numeric(dIF$dIF)
pvals_gaussian <- left_join(pvals_gaussian,dIF,by='transcript_id')


###########volcano plot############
df_volcano <- pvals_gaussian
df_volcano$group <- 'NS'
df_volcano$group[df_volcano$dIF< (-0.1) & df_volcano$tissue.adjp<0.05] <-'DOWN'
df_volcano$group[df_volcano$dIF > (0.1) & df_volcano$tissue.adjp<0.05] <- 'UP'
df_volcano$gene_transcript <- paste0(df_volcano$gene_id,': ',df_volcano$transcript_id)
df_volcano$label[df_volcano$group !='NS'] <- df_volcano$gene_transcript[df_volcano$group != 'NS'] #gene label

volcano_gaussian <- ggplot(data=df_volcano,mapping=aes(x=dIF,y=-log10(tissue.adjp),color=group,label=label))+
  geom_point()+
  #xlim(-0.228,0.228) +
  labs(title='volcano plot of significance against difference in isoform fraction between recurrence and non-recurrence') +
  geom_label_repel(na.rm = TRUE, show.legend = FALSE, box.padding = 0.5) +
  theme_minimal()+ 
  theme(legend.position="none")+
  geom_vline(xintercept=c(-0.1, 0.1), col="red", linetype='longdash') +
  geom_hline(yintercept=-log10(0.05), col="red", linetype='longdash') +
  geom_vline(xintercept = 0, col="grey", linetype='longdash') +
  scale_color_manual(values=c("DOWN"= "#56B4E9", "NS"="#999999", "UP"="#E69F00"))

volcano_gaussian


#filter to only include the DEGs from r0 vs. r1 - this is from df_volcano computed earlier
DEG <- df_volcano[df_volcano$affected.adjp<0.05 & (df_volcano$log2FC>0.25| df_volcano$log2FC< (-0.25)),]

#read in pheno data
pheno <- read_csv("/Users/swashburn30/Desktop/isoform/kiera_analysis/phenodata_2_9_24.csv", col_names = F)
pheno <- as.data.frame(pheno)
colnames(pheno) <- c("donor", 'geo_id', 'sex', 'tissue', 'sample', 'age', 'race')
row.names(pheno)<-pheno$sample
pheno <-na.omit(pheno)

#read in gene expression data 
gene_df <- read.table('/Users/swashburn30/Desktop/isoform/kiera_analysis/rsem.merged.gene_tpm.tsv', sep = '\t', header = TRUE)
row.names(gene_df) <- gene_df$gene_id
gene_mapping <- gene_df %>% dplyr::select(transcript_id.s.,gene_id) 
gene_df <- gene_df%>% dplyr::select(row.names(pheno)) #select samples - 334 samples
gene_df <- na.omit(gene_df)  #remove NA (where the tpm for gene is 0) #29972 genes

gene_df <- gene_df[,row.names(pheno)] #reorder

#filter genes based on DEG from r0 vs r1
gene_df <- gene_df[rownames(gene_df) %in% DEG$gene_id, ]

#transpose 
gene_trans <- as.data.frame(t(gene_df)) #transpose

#scale the log2 normalize the gene expression level 
gene_trans <- log2(gene_trans + 1)

#perform PCA with these genes 

df_pca <- prcomp(gene_trans)

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


pc_score_metadata <- full_join(pc_scores, pheno, by = ("sample"))
test_pc_score_metadata <- na.omit(pc_score_metadata)

pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(tissue))) +
  geom_point() 
#geom_text(aes(label = sample))

p <- pc_score_metadata %>%
  ggplot( aes(x=tissue, y=PC1, fill=tissue)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("DEG PC1 score in all samples (n=20)") +
  xlab("tissue") +
  ylab("PC1") +
  scale_x_discrete(labels=c("ileum" = "ileum (n=10)", "rectum" = "rectum (n=10)"))

p + stat_compare_means(method = "wilcox.test", label.x = 2)

ileum <- subset(pc_score_metadata, tissue == 'ileum')
rectum <- subset(pc_score_metadata, tissue == 'rectum')

res <- wilcox.test(ileum$PC1, rectum$PC1)
res
#p-value = 0.01854

mean(ileum$PC1) #-20.75314
mean(rectum$PC1) #20.75314

#DTU PC1 score female 

#extract females from DF
female_df <- pc_score_metadata[pc_score_metadata$sex == 'female', ]

fp <- female_df %>%
  ggplot( aes(x=tissue, y=PC1, fill=tissue)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("DEG PC1 score in all female (n=10)") +
  xlab("tissue") +
  ylab("PC1") +
  scale_x_discrete(labels=c("ileum" = "ileum (n=5)", "rectum" = "rectum (n=5)"))

fp + stat_compare_means(method = "wilcox.test", label.x = 2)

f_ileum <- subset(female_df, tissue == 'ileum')
f_rectum <- subset(female_df, tissue == 'rectum')

fres <- wilcox.test(f_ileum$PC1, f_rectum$PC1)
fres
#p-value = 0.6905

mean(f_ileum$PC1) #-12.44206
mean(f_rectum$PC1) #16.79153

#DTU PC1 score male 

#extract male from DF
male_df <- pc_score_metadata[pc_score_metadata$sex == 'male', ]

mp <- male_df %>%
  ggplot( aes(x=tissue, y=PC1, fill=tissue)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("DEG PC1 score in all male (n=10)") +
  xlab("tissue") +
  ylab("PC1") +
  scale_x_discrete(labels=c("ileum" = "ileum (n=5)", "rectum" = "rectum (n=5)"))

mp + stat_compare_means(method = "wilcox.test", label.x = 2)


m_ileum <- subset(male_df, tissue == 'ileum')
m_rectum <- subset(male_df, tissue == 'rectum')

mres <- wilcox.test(m_ileum$PC1, m_rectum$PC1)
mres
#p-value = 0.007937

mean(m_ileum$PC1) #-29.06422
mean(m_rectum$PC1) #24.71476

#------------use recurring vs. non recurring DEGs in validation cohort-------------#

#filter to only include the DEGs from r0 vs. r1 - this is from df_volcano computed earlier - from post-op dataset 
DEG <- df_volcano[df_volcano$affected.adjp<0.05 & (df_volcano$log2FC>0.25| df_volcano$log2FC< (-0.25)),]

#read in pheno data
pheno <- read_csv("/Users/swashburn30/Desktop/isoform/kiera_analysis/phenodata_2_9_24.csv", col_names = F)
pheno <- as.data.frame(pheno)
colnames(pheno) <- c("donor", 'geo_id', 'sex', 'tissue', 'sample', 'age', 'race')
row.names(pheno)<-pheno$sample
pheno <-na.omit(pheno)

#read in gene expression data 
gene_df <- read.table('/Users/swashburn30/Desktop/isoform/kiera_analysis/rsem.merged.gene_tpm.tsv', sep = '\t', header = TRUE)
row.names(gene_df) <- gene_df$gene_id
gene_mapping <- gene_df %>% dplyr::select(transcript_id.s.,gene_id) 
gene_df <- gene_df%>% dplyr::select(row.names(pheno)) #select samples - 334 samples
gene_df <- na.omit(gene_df)  #remove NA (where the tpm for gene is 0) #29972 genes

gene_df <- gene_df[,row.names(pheno)] #reorder

#filter genes based on DEG from r0 vs r1
gene_df <- gene_df[rownames(gene_df) %in% DEG$gene_id, ]

#transpose 
gene_trans <- as.data.frame(t(gene_df)) #transpose

#scale the log2 normalize the gene expression level 
gene_trans <- log2(gene_trans + 1)

#perform PCA with these genes 

df_pca <- prcomp(gene_trans)

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


pc_score_metadata <- full_join(pc_scores, pheno, by = ("sample"))
test_pc_score_metadata <- na.omit(pc_score_metadata)

pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(tissue))) +
  geom_point() 
#geom_text(aes(label = sample))

p <- pc_score_metadata %>%
  ggplot( aes(x=tissue, y=PC1, fill=tissue)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("DEG PC1 score in all samples (n=20)") +
  xlab("tissue") +
  ylab("PC1") +
  scale_x_discrete(labels=c("ileum" = "ileum (n=10)", "rectum" = "rectum (n=10)"))

p + stat_compare_means(method = "wilcox.test", label.x = 2)

ileum <- subset(pc_score_metadata, tissue == 'ileum')
rectum <- subset(pc_score_metadata, tissue == 'rectum')

res <- wilcox.test(ileum$PC1, rectum$PC1)
res
#p-value = 0.01854

mean(ileum$PC1) #-20.75314
mean(rectum$PC1) #20.75314

#DTU PC1 score female 

#extract females from DF
female_df <- pc_score_metadata[pc_score_metadata$sex == 'female', ]

fp <- female_df %>%
  ggplot( aes(x=tissue, y=PC1, fill=tissue)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("DEG PC1 score in all female (n=10)") +
  xlab("tissue") +
  ylab("PC1") +
  scale_x_discrete(labels=c("ileum" = "ileum (n=5)", "rectum" = "rectum (n=5)"))

fp + stat_compare_means(method = "wilcox.test", label.x = 2)

f_ileum <- subset(female_df, tissue == 'ileum')
f_rectum <- subset(female_df, tissue == 'rectum')

fres <- wilcox.test(f_ileum$PC1, f_rectum$PC1)
fres
#p-value = 0.6905

mean(f_ileum$PC1) #-12.44206
mean(f_rectum$PC1) #16.79153

#DTU PC1 score male 

#extract male from DF
male_df <- pc_score_metadata[pc_score_metadata$sex == 'male', ]

mp <- male_df %>%
  ggplot( aes(x=tissue, y=PC1, fill=tissue)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("DEG PC1 score in all male (n=10)") +
  xlab("tissue") +
  ylab("PC1") +
  scale_x_discrete(labels=c("ileum" = "ileum (n=5)", "rectum" = "rectum (n=5)"))

mp + stat_compare_means(method = "wilcox.test", label.x = 2)


m_ileum <- subset(male_df, tissue == 'ileum')
m_rectum <- subset(male_df, tissue == 'rectum')

mres <- wilcox.test(m_ileum$PC1, m_rectum$PC1)
mres
#p-value = 0.007937

mean(m_ileum$PC1) #-29.06422
mean(m_rectum$PC1) #24.71476



#----------now compute DEG comparing rectum vs. ileum in validation dataset------------# 
#read in pheno data
pheno <- read_csv("/Users/swashburn30/Desktop/isoform/kiera_analysis/phenodata_2_9_24.csv", col_names = F)
pheno <- as.data.frame(pheno)
colnames(pheno) <- c("donor", 'geo_id', 'sex', 'tissue', 'sample', 'age', 'race')
row.names(pheno)<-pheno$sample
pheno <-na.omit(pheno)
pheno$donor_label <- c("donor1", "donor1", "donor2", "donor2", "donor3", "donor3", "donor4", "donor4", 'donor5', "donor5", "donor6", "donor6", "donor7", "donor7", "donor8", "donor8", "donor9", "donor9", "donor10", "donor10")

#read in gene expression data 
gene_df <- read.table('/Users/swashburn30/Desktop/isoform/kiera_analysis/rsem.merged.gene_tpm.tsv', sep = '\t', header = TRUE)
row.names(gene_df) <- gene_df$gene_id
gene_mapping <- gene_df %>% dplyr::select(transcript_id.s.,gene_id) 
gene_df <- gene_df%>% dplyr::select(row.names(pheno)) #select samples - 334 samples
gene_df <- na.omit(gene_df)  #remove NA (where the tpm for gene is 0) #29972 genes
gene_df <- gene_df[ rowMeans(gene_df) > 5, ] #mean coverage > 5 reads - 5505 genes

gene_df <- gene_df[,row.names(pheno)] #reorder
gene_trans <- as.data.frame(t(gene_df)) #transpose

#scale the log2 normalize the gene expression level 
gene_trans <- log2(gene_trans + 1)

#create the design matrix for linear model 
gene_design <- cbind(pheno[,c("tissue", "sex", "race", "donor_label")],gene_trans[,1])
names(gene_design) <- c("tissue", "sex", "race", "donor_label", "gene")


gene_design <- na.omit(gene_design) #remove na values

#try removing row 3956 and see if this helps
gene_trans <- gene_trans[ ,-3956] #RNU5A-1
gene_df <- gene_df[-3956, ]

############## linear mixed model, gaussian distribution#######

model.full = lmer(gene ~ tissue + sex + race + (1|donor_label),  data=gene_design ,REML = FALSE) #full model

model.tissue = lmer(gene ~ sex + race + (1|donor_label),  data=gene_design, REML = FALSE ) #model missing affected


#extract pvalue using maximum likelyhood
prob_tissue <- anova(model.full,model.tissue)[2,8]

pvals <- matrix(c(prob_tissue),nrow=1,ncol=1)



#loop through the rest of isoforms
for (m in 2:ncol(gene_trans)){
  gene_design <- cbind(pheno[,c("tissue", "sex", "race", "donor_label")],gene_trans[,m])
  names(gene_design) <- c("tissue", "sex", "race", "donor_label", "gene")
  
  gene_design <- na.omit(gene_design) #remove na values
  
  
  model.full = lmer(gene ~ tissue + sex + race + (1|donor_label),  data=gene_design, REML = FALSE ) #full model
  
  model.tissue = lmer(gene ~ sex + race + (1|donor_label),  data=gene_design, REML = FALSE ) #model missing affected
  
  
  #extract pvalue using maximum likelyhood
  prob_tissue <- anova(model.full,model.tissue)[2,8]
  pvals <- rbind(pvals, c(prob_tissue))
}

row.names(pvals) <- row.names(gene_df)
colnames(pvals) <- c('tissue.p')




#adjust p value - 11,137 genes
pvals_gaussian <- as.data.frame(pvals)
pvals_gaussian$tissue.adjp <- p.adjust(pvals_gaussian[,1], method='fdr')
pvals_gaussian$gene_id<- row.names(pvals_gaussian)
pvals_gaussian<- left_join(pvals_gaussian,gene_df,by=rownames(gene_df))

#calculate log 2 FC 
meta_gene <- cbind(gene_trans,pheno)

rectum_samples <- meta_gene[meta_gene$tissue == "rectum", ]
ileum_samples <- meta_gene[meta_gene$tissue == "ileum", ]

# Calculate the mean expression for each gene in each condition
mean_expression_rectum <- as.data.frame(colMeans(rectum_samples[, -c(11138, 11139, 11140, 11141, 11142, 11143 ,11144, 11145)], na.rm = TRUE))
colnames(mean_expression_rectum) <- "rectum"
mean_expression_ileum <- as.data.frame(colMeans(ileum_samples[, -c(11138, 11139, 11140, 11141, 11142, 11143 ,11144, 11145)], na.rm = TRUE))
colnames(mean_expression_ileum) <- 'ileum'

# Calculate the log2-fold change - already log2 normalized 
test <- (mean_expression_rectum - mean_expression_ileum)
test$gene_id <- rownames(test)
colnames(test) <- c("log2FC", "gene_id")


###########volcano plot############

pvals_gaussian <- left_join(pvals_gaussian,test,by='gene_id')
df_volcano <- pvals_gaussian
df_volcano$group <- 'NS'
df_volcano$group[df_volcano$log2FC< (-0.25) & df_volcano$tissue.adjp<0.05] <-'DOWN'
df_volcano$group[df_volcano$log2FC > (0.25) & df_volcano$tissue.adjp<0.05] <- 'UP'
#df_volcano$gene_transcript <- paste0(df_volcano$gene_id,': ',df_volcano$transcript_id)
df_volcano$label[df_volcano$group !='NS'] <- df_volcano$gene_id[df_volcano$group != 'NS'] #gene label

volcano_gaussian <- ggplot(data=df_volcano,mapping=aes(x=log2FC,y=-log10(tissue.adjp),color=group,label=label))+
  geom_point()+
  #xlim(-0.228,0.228) +
  labs(title='volcano plot of significance against DEG rectum and ileum') +
  geom_label_repel(na.rm = TRUE, show.legend = FALSE, box.padding = 0.5) +
  theme_minimal()+ 
  theme(legend.position="none")+
  geom_vline(xintercept=c(-0.25, 0.25), col="red", linetype='longdash') +
  geom_hline(yintercept=-log10(0.05), col="red", linetype='longdash') +
  geom_vline(xintercept = 0, col="grey", linetype='longdash') +
  scale_color_manual(values=c("DOWN"= "#785EF0", "NS"="#999999", "UP"="#DC267F"))

volcano_gaussian


png(file = "/Users/swashburn30/Desktop/isoform/validate_kates_results/DEG_volcano_plot_rectum_vs_ileum.png", width = 1750, height = 1000)
plot(volcano_gaussian)
dev.off()


#save the df_volcano - has the log2FC 
write.csv(pvals_gaussian, file = "/Users/swashburn30/Desktop/isoform/kiera_analysis/DEG_rectum_vs_ileum_results_2_19_2024.csv")

############PCA score#########
library(ggpubr)
library(glue)
library(dplyr)
library(ggpubr)
#stringent DEG filter: adjp < 0.05 and DEG > ± 0.25 - 762 DEGs
DEG <- pvals_gaussian[pvals_gaussian$tissue.adjp<0.05 & (pvals_gaussian$log2FC>0.25| pvals_gaussian$log2FC< (-0.25)),]

#stringent filter of +- 1
DEG_filt <- DEG[DEG$tissue.adjp<0.05 & (DEG$log2FC>1| DEG$log2FC< (-1)),]


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
  ggplot(aes(x = PC1, y = PC2, colour = factor(tissue))) +
  geom_point() 
#

pc_score_metadata %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(sex))) +
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
  ggplot( aes(x=tissue, y=PC1, fill=tissue)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("DEG PC1 score in all samples (n=20)") +
  xlab("Tissue") +
  ylab("PC1") +
  scale_x_discrete(labels=c("ileum" = "ileum (n=10)", "rectum" = "rectum (n=10)"))

p + stat_compare_means(method = "wilcox.test", label.x = 2)

ileum <- subset(pc_score_metadata, tissue == 'ileum')
rectum <- subset(pc_score_metadata, tissue == 'rectum')

res <- wilcox.test(ileum$PC1, rectum$PC1)
res
#p-value =1.083e-05

mean(ileum$PC1) #-36.02912
mean(rectum$PC1) #36.02912

#DEG PC1 score female 

#extract females from DF
female_df <- select(pc_score_metadata, "sample", "PC1", "sex", "tissue", )
female_df <- female_df[female_df$sex == 'female', ]

fp <- female_df %>%
  ggplot( aes(x=tissue, y=PC1, fill=tissue)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("DEG PC1 score in all female (n=10)") +
  xlab("Tissue") +
  ylab("PC1") +
  scale_x_discrete(labels=c("ileum" = "ileum (n=5)", "rectum" = "rectum (n=5)"))

fp + stat_compare_means(method = "wilcox.test", label.x = 1)


f_ileum <- subset(female_df, tissue == 'ileum')
f_rectum <- subset(female_df, tissue == 'rectum')

fres <- wilcox.test(f_ileum$PC1, f_rectum$PC1)
fres
#p-value = 0.007937

mean(f_ileum$PC1) #-32.50896
mean(f_rectum$PC1) #34.76156

#DEG PC1 score male 

#extract male from DF
male_df <- select(pc_score_metadata, "sample", "PC1", "sex", "tissue", )
male_df <- male_df[male_df$sex == 'male', ]

mp <- male_df %>%
  ggplot( aes(x=tissue, y=PC1, fill=tissue)) +
  geom_boxplot(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("DEG PC1 score in all male (n=10)") +
  xlab("Tissue") +
  ylab("PC1") +
  scale_x_discrete(labels=c("ileum" = "ileum (n=5)", "rectum" = "rectum (n=5)"))

mp + stat_compare_means(method = "wilcox.test", label.x = 1)

m_ileum <- subset(male_df, tissue == 'ileum')
m_rectum <- subset(male_df, tissue == 'rectum')

mres <- wilcox.test(m_ileum$PC1, m_rectum$PC1)
mres
#p-value = 0.007937

mean(m_ileum$PC1) #-39.54928
mean(m_rectum$PC1) #37.29668

#-----------p53 investigation-------------#

#read in pheno data
pheno <- read_csv("/Users/swashburn30/Desktop/isoform/kiera_analysis/phenodata_2_9_24.csv", col_names = F)
pheno <- as.data.frame(pheno)
colnames(pheno) <- c("donor", 'geo_id', 'sex', 'tissue', 'sample', 'age', 'race')
row.names(pheno)<-pheno$sample
pheno <-na.omit(pheno)
pheno$donor_label <- c("donor1", "donor1", "donor2", "donor2", "donor3", "donor3", "donor4", "donor4", 'donor5', "donor5", "donor6", "donor6", "donor7", "donor7", "donor8", "donor8", "donor9", "donor9", "donor10", "donor10")

#read in gene expression data 
gene_df <- read.table('/Users/swashburn30/Desktop/isoform/kiera_analysis/rsem.merged.gene_tpm.tsv', sep = '\t', header = TRUE)
row.names(gene_df) <- gene_df$gene_id
gene_mapping <- gene_df %>% dplyr::select(transcript_id.s.,gene_id) 
gene_df <- gene_df%>% dplyr::select(row.names(pheno)) #select samples - 334 samples
gene_df <- na.omit(gene_df)  #remove NA (where the tpm for gene is 0) #29972 genes
gene_df <- gene_df[ rowMeans(gene_df) > 5, ] #mean coverage > 5 reads - 5505 genes

gene_df <- gene_df[,row.names(pheno)] #reorder
gene_trans <- as.data.frame(t(gene_df)) #transpose

#scale the log2 normalize the gene expression level 
gene_trans <- log2(gene_trans + 1)

#load in DEGs 
tissue_DEG <- read_csv(file = "/Users/swashburn30/Desktop/isoform/kiera_analysis/DEG_filt_rectum_vs_ileum_results_2_19_2024.csv")


#merge log 2 normalized gene expression level with pheno data 
df <- cbind(gene_trans, pheno)
#Violin Plot of TP53INP2 expression in ileum vs. rectum 

TP53INP2_violin <- df %>%
  ggplot( aes(x=tissue, y=TP53INP2, fill=tissue)) +
  geom_violin(alpha=0.3) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  #theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=10),
    text = element_text(size = 10)
  ) +
  ggtitle("TP53INP2") +
  stat_compare_means(method = "wilcox.test", label.x.npc = "middle", size = 4, label = paste0("P" = "p.format")) +
  xlab("Recurrence Status") +
  ylab("Normalized gene expression") +
  scale_x_discrete(labels=c("ileum" = "ileum (n=10)", "rectum" = "rectum (n=10)"))

plot(TP53INP2_violin + geom_boxplot(width=0.1))


#------test GSEA for p53--------#
library(fgsea)

rank <- tissue_DEG$log2FC
names(rank) <- tissue_DEG$gene_id
barplot(sort(rank, decreasing = T))
rank = sort(rank, decreasing = TRUE)

#kegg pathways 
BiocManager::install("KEGGREST")
library(KEGGREST)
library(org.Hs.eg.db)
data(KEGGREST)

pathways <- gmtPathways("/Users/swashburn30/Downloads/h.all.v2023.2.Hs.symbols.gmt")



fgseaRes <- fgsea(pathways = pathways, 
                  stats    = rank,
                  minSize  = 15,
                  maxSize  = 500)

plotEnrichment(pathways[["HALLMARK_P53_PATHWAY"]],
               rank)

cell <- gmtPathways("/Users/swashburn30/Downloads/c8.all.v2023.2.Hs.symbols.gmt")

cell_res <- fgsea(pathways = cell, 
                  stats    = rank,
                  minSize  = 15,
                  maxSize  = 500)

#plotEnrichment(cell[["BUSSLINGER_DUODENAL_DIFFERENTIATING_STEM_CELLS"]],
#ranks)

go <- gmtPathways("/Users/swashburn30/Downloads/c5.go.v2023.2.Hs.symbols.gmt")

go_res <- fgsea(pathways = go, 
                stats    = rank,
                minSize  = 15,
                maxSize  = 500)

go_res <- go_res[go_res$padj < 0.05, ]

plotEnrichment(go[["GOBP_REGULATION_OF_CELL_CYCLE_PROCESS"]],
               ranks)

p53 <- gmtPathways("/Users/swashburn30/Downloads/P53_DN.V1_UP.v2024.1.Hs.gmt")

p53_res <- fgsea(pathways = p53, 
                 stats    = rank,
                 minSize  = 15,
                 maxSize  = 500)

plotEnrichment(p53[["P53_DN.V1_UP"]],
               rank)

#genes from hallmark p53 pathway

# Replace "Your Pathway Name" with the name of the pathway you are interested in
p53_gene_hall <- fgseaRes %>%
  filter(pathway == "HALLMARK_P53_PATHWAY") %>%         # Filter for the pathway of interest
  pull(leadingEdge)                                  # Extract the leadingEdge genes for that pathway


#genes from p53 oncogene pathway 
p53_gene_onco <- p53_res %>%
  filter(pathway == "P53_DN.V1_UP") %>%         # Filter for the pathway of interest
  pull(leadingEdge)                                  # Extract the leadingEdge genes for that pathway

