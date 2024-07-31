library(ggplot2)
library(dplyr)
library(png)
library(grid)
library(gridExtra)
library(ecodist)
library(vegan)
library(factoextra)
library(lemon)

AUCs <- read.csv("./csv_files/semi_leaky_AUCs.csv", colClasses = c("DEBIAS" = "factor"))
post_DEBIAS <- list.files("./cluster_runs/DEBIAS-M_semi_leaky/csv_files/DEBIAS-M_runs/semi_leaky/",
                          pattern = "debiased_lognorm", full.names = TRUE)
debias_weights <- list.files("./cluster_runs/DEBIAS-M_semi_leaky/csv_files/DEBIAS-M_runs/semi_leaky/",
                             pattern = "debias_weights", full.names = TRUE)
phenos <- read.csv("./csv_files/phenotypes.csv")
pval.df <- read.csv("./csv_files/AUC_pvals_semi_leaky.csv")
neg_pval <- read.csv("./csv_files/DEBIAS-M_AUC_pvals.csv")  %>% mutate(DEBIAS = if_else(DEBIAS == 1, 0.5, 0))
neg_auc <- read.csv("./cluster_runs/DEBIAS-M_100_perm_fixed/csv_files/AUCs.csv") %>% mutate(DEBIAS = 
                                                                                              if_else(DEBIAS == 1, 0.5, 0))
IDs <- distinct(AUCs, Study_ID)$Study_ID

a = list.files("./cluster_runs/DEBIAS-M_semi_leaky/csv_files/AUCs/semi_leaky/",
               pattern = "post_DEBIAS", full.names = TRUE)

AUCs = read.csv(a[1])
for(i in 2:length(a)){
  d = read.csv(a[i])
  AUCs = rbind(AUCs, d) 
}
write.csv(AUCs, "./csv_files/semi_leaky_AUCs.csv", row.names = FALSE)
##### histograms and p-val df ####
# can make pval.df or histograms

pval.df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(pval.df) = c("Study_ID", "pval", "DEBIAS")

for(i in 1:length(IDs)){
  
  phen <- phenos[phenos$ID == IDs[i],]$phenotype
  
  auc.df <- filter(AUCs, Study_ID == IDs[i])
  
  a <- auc.df[auc.df$Permutation == 0,]$AUC
  samp <- auc.df[auc.df$Permutation == 1,]$AUC
  z = (a-mean(samp))/(sd(samp)/sqrt(1))
  for.pval = pnorm(z, lower.tail = FALSE)
  
  pval.df[nrow(pval.df) + 1,] = c(IDs[i], for.pval, TRUE)
  
  # png(paste("./cluster_runs/DEBIAS-M_100_perm_fixed/output/pre_DEBIAS-M_RF_lognorm_histogram_fixed_", IDs[i], ".png", sep=""), width = 480, height = 480)
  # g <- ggplot() + geom_histogram(data = filter(auc.df, Permutation == TRUE), aes(x = AUC), bins = 40) +
  #   geom_vline(filter(auc.df, Permutation == FALSE), mapping = aes(xintercept=AUC), color = "cornflowerblue") +
  #   labs(title = paste("Training without: ", IDs[i], " (", phen, ")", sep=""), y = "count") +
  #   annotate("label", x=min(auc.df$AUC)+.01, y=5, size = 3, label = paste("p= ", signif(for.pval, digits=3), sep="")) +
  #   scale_y_continuous(expand = expansion(mult = c(0, .1)))
  # 
  # 
  # print(g)
  # dev.off()
}

pval.df <- merge(pval.df, phenos, by.x = "Study_ID", by.y = "ID", all.x = TRUE)
pval.df$phenotype[grep("cancer", pval.df$phenotype)] <- "cancer"
write.csv(pval.df, "./csv_files/AUC_pvals_semi_leaky.csv", row.names = FALSE)



#### save histograms as a pdf ####

all_images <- list.files("./cluster_runs/DEBIAS-M_semi_leaky/output/", full.names = TRUE)
ind <- grep("post_DEBIAS-M_RF_lognorm_hist", all_images)
all_images <- all_images[ind]

plots <- lapply(ll <- all_images, function(x){
  img <- as.raster(readPNG(x))
  rasterGrob(img, interpolate = FALSE)
})

ggsave("./output/semi_leaky/hist_100_perm_semi_leaky.pdf", marrangeGrob(grobs=plots, nrow=2, ncol=2), width = 4, height = 4, dpi = 300)

#### box plots ####

# auc box plot
a = AUCs
AUCs = a

# t-test between leaky and semi-leaky: t = 0.57981, df = 3027.8, p-value = 0.5621
neg_auc$DEBIAS = as.factor(neg_auc$DEBIAS)
AUCs <- rbind(AUCs, neg_auc)
t_auc = AUCs %>% select(DEBIAS, AUC) %>% filter(DEBIAS != 0)
with(t_auc, t.test(AUC[DEBIAS == 1], AUC[DEBIAS == 0.5]), alternative = "two.sided")

AUCs$DEBIAS <- as.numeric(levels(AUCs$DEBIAS))[(AUCs$DEBIAS)]
selected_AUCs <- AUCs %>% filter(Permutation == 0) %>% merge(y = phenos, by.x = "Study_ID", by.y = "ID", all.x = TRUE) %>% arrange(desc(DEBIAS))
selected_AUCs$phenotype[grep("cancer", selected_AUCs$phenotype)] <- "cancer"
auc_b <- ggplot(selected_AUCs, mapping = aes(x=as.factor(DEBIAS), y=AUC)) + 
  geom_boxplot() + 
  geom_point(aes(col = as.factor(phenotype)), alpha = 2) +
  scale_color_brewer(name = "pehnotype", palette = "Paired", guide = "none") +
  scale_x_discrete(labels = c("Before DEBIAS-M", "After DEBIAS-M (leaky)", "After DEBIAS-M (semi-leaky)")) +
  labs(x = "", title = "AUCs before and after DEBIAS-M") +
  ylim(0.4, 1.0) +
  theme_dark()

auc_b

# p-val box plot
pval.df = p
# p= pval.df

# t-test between leaky and semi-leaky: t = 2.1986, df = 21.152, p-value = 0.03916


neg_pval$DEBIAS = as.factor(neg_pval$DEBIAS)
pval.df$DEBIAS = as.factor(pval.df$DEBIAS)
pval.df <- rbind(pval.df, neg_pval)
t_pval = pval.df %>% select(DEBIAS, pval) %>% filter(DEBIAS != 0)
with(t_pval, t.test(pval[DEBIAS == 1], pval[DEBIAS == 0.5]), alternative = "two.sided")

pval.df$DEBIAS <- as.numeric(levels(pval.df$DEBIAS))[(pval.df$DEBIAS)]
pval.df <- pval.df %>% arrange(desc(DEBIAS))
pval_b <- ggplot(pval.df, mapping = aes(x = as.factor(DEBIAS), y = log(pval, base = 10))) + 
  geom_boxplot() +
  geom_point(aes(col = as.factor(phenotype))) +
  scale_color_brewer(name = "phenotype", palette = "Paired") +
  scale_x_discrete(labels = c("Before DEBIAS-M", "After DEBIAS-M (leaky)", "After DEBIAS-M (semi-leaky)")) +
  labs(x = "", y = "log10 p-value", title = "p-values before and after DEBIAS-M") +
  theme_dark()

pval_b

png("./output/semi_leaky/box_AUC_pval_semi_leaky.png", width = 1000, height = 480)
grid.arrange(auc_b, pval_b, nrow = 1)
dev.off()


#### pcoa ####

pcoa <- function(df, chosen_title, pc = NULL){
  bray <- vegdist(df, method = "bray")
  pcoa_val <- pco(bray, negvals = "zero", dround = 0)
  
  if(!is.null(pc)){
    
    cat("example to access elements in the list: \n
    l = pcoa(df)\n
    l[1] # returns bray-curtis values\n
    l[2] # returns the pcoa results\n
    ls[2][[1]] # returns vectors and values which can be accessed via $\n
    ex: ls[2][[1]]$vectors")
    
    return(list(bray, pcoa_val))
  }
  
  pco.labels = lapply(row.names(pcoa_val$vectors), function(x) unlist(strsplit(x, split = ".", fixed=TRUE))[1]) # remove .numbers from the end of the row names
  pcoa_val.df = data.frame(Study_ID = unlist(pco.labels),
                           PCoA1 = pcoa_val$vectors[,1], 
                           PCoA2 = pcoa_val$vectors[,2])
  
  pcoa_val.df = merge(pcoa_val.df, phenos, by.x = "Study_ID", by.y = "ID", all.x = TRUE)
  pcoa_val.df$phenotype[grep("cancer", pcoa_val.df$phenotype)] = "cancer"
  
  pco.plot <- ggplot(data = pcoa_val.df, mapping = aes(x = PCoA1, y = PCoA2)) + 
    geom_point(aes(col = as.factor(phenotype)), alpha = 0.7) +
    scale_color_brewer(name = "phenotype", palette = "Paired") +
    labs(title = chosen_title)
  return(pco.plot)
}


# after debias-m count pca

pd = post_DEBIAS

for(i in 1:length(pd)){
  post_DEBIAS <- read.csv(pd[i])
  rownames(post_DEBIAS) <- make.names(post_DEBIAS$Study_ID, unique = TRUE)
  rownames(post_DEBIAS) <- gsub("^X", "", rownames(post_DEBIAS))
  
  if((i+1)%%2 == 0 & (i+1)%%4 != 0){
    p1 = pcoa(post_DEBIAS[4:length(post_DEBIAS)], paste("PCoA on ", IDs[i], " count data after DEBIAS-M", sep = ""))
    print(paste(i, "p1"))}

  if(i%%2 == 0 & i%%4 != 0){
    p2 = pcoa(post_DEBIAS[4:length(post_DEBIAS)], paste("PCoA on ", IDs[i], " count data after DEBIAS-M", sep = ""))
    print(paste(i, "p2"))}

  if((i+1)%%4 == 0){
    p3 = pcoa(post_DEBIAS[4:length(post_DEBIAS)], paste("PCoA on ", IDs[i], " count data after DEBIAS-M", sep = ""))
    print(paste(i, "p3"))}

  if(i%%4 == 0){
    p4 = pcoa(post_DEBIAS[4:length(post_DEBIAS)], paste("PCoA on ", IDs[i], " count data after DEBIAS-M", sep = ""))
    print(paste(i, "p4"))
    
  
    png(paste("./output/semi_leaky/pcoa_post_DEBIAS",i-3,"-",i ,".png", sep=""))
    p1
    grid_arrange_shared_legend(p1, p2, p3, p4, ncol = 2, nrow = 2)
    print(p)
    dev.off()
    }
  print(paste(i, " of ", length(pd), " done ", sep = ""))
}

#### PCA weights ####
# debias-m weights
grid_arrange_shared_legend(p1, p2, p3, p4, ncol = 2, nrow = 2)

dw = debias_weights
for(i in 1:length(dw)){
  debias_weights = t(debias_weights)
  debias_weights = data.frame(debias_weights)
  # rownames(debias_weights) <- debias_weights$Study_ID
  # debias_weights = debias_weights[,2:length(debias_weights)]
  scaled <- as.data.frame(scale(debias_weights))
  # verify mean == 0 and sd == 1
  # sapply(scaled, mean)
  # sapply(scaled, sd)
  
  
  dw_pca <- prcomp(scaled)
  var = get_pca_var(dw_pca)
  dwpcacontr.df = as.data.frame(var$contrib)
  pca_contr = summary(dw_pca)$importance[2,]
  
  dwpca.df = as.data.frame(dw_pca$x)
  dwpca.df$Study_ID = row.names(dwpca.df)
  dwpca.df = merge(dwpca.df, phenos, by.x = "Study_ID", by.y = "ID", all.x = TRUE)
  dwpca.df$phenotype[grep("cancer", dwpca.df$phenotype)] = "cancer"
  
  dwpca.plot = ggplot(dwpca.df, mapping = aes(x = PC1, y = PC2)) +
    geom_point(aes(col = as.factor(phenotype))) +
    scale_color_brewer(name = "phenotype", palette = "Paired") +
    labs(title = "PCA of weights", x = paste("PC1 (", pca_contr[1], ")", sep = ""), y = paste("PC2 (", pca_contr[2],")", sep=""))
  
  png(paste("./output/semi_leaky/PCA_weights_PC1-2_leaky",IDs[i],".png"))
  print(dwpca.plot)
  dev.off()
  print(paste(i, " of ", length(dw), " done.", sep = ""))
  
}

# same plot but with study labels instead
dwpca_labels.plot = ggplot(dwpca.df, mapping = aes(x = PC1, y = PC2)) +
  geom_point(color = "white") +
  geom_text(mapping = aes(label = Study_ID)) +
  labs(title = "PCA of weights")
dwpca_labels.plot

png("./output/test_permutated/PCA_weights_PC1-2_leaky.png")
dwpca.plot
dev.off()
