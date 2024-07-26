library(ggplot2)
library(dplyr)
library(png)
library(grid)
library(gridExtra)
library(ecodist)

AUCs <- read.csv("./cluster_runs/pre_DEBIAS-M_100_perm/csv_files/AUCs.csv", colClasses = c("DEBIAS" = "factor"))
phenos <- read.csv("./csv_files/phenotypes.csv")
IDs <- distinct(AUCs, Study_ID)$Study_ID

dim(auc.df)

# pval.df <- data.frame(matrix(ncol = 3, nrow = 0))
# colnames(pval.df) = c("Study_ID", "pval", "DEBIAS")
selected_AUCs <- filter(AUCs, DEBIAS == 0)

for(i in 1:length(IDs)){
  
  phen <- phenos[phenos$ID == IDs[i],]$phenotype
  
  auc.df <- filter(selected_AUCs, Study_ID == IDs[i])
  
  a <- auc.df[auc.df$Permutation == 0,]$AUC
  samp <- auc.df[auc.df$Permutation == 1,]$AUC
  z = (a-mean(samp))/(sd(samp)/sqrt(1))
  for.pval = pnorm(z, lower.tail = FALSE)
  
  # pval.df[nrow(pval.df) + 1,] = c(IDs[i], for.pval, FALSE)
  
  png(paste("./cluster_runs/pre_DEBIAS-M_100_perm/output/by_hand/post_DEBIAS-M_RF_lognorm_histogram_", IDs[i], ".png", sep=""), width = 480, height = 480)
  g <- ggplot() + geom_histogram(data = filter(auc.df, Permutation == TRUE), aes(x = AUC), bins = 40) +
    geom_vline(filter(auc.df, Permutation == FALSE), mapping = aes(xintercept=AUC), color = "cornflowerblue") +
    labs(title = paste("Training without: ", IDs[i], " (", phen, ")", sep=""), y = "count") +
    annotate("label", x=min(auc.df$AUC)+.01, y=5, size = 3, label = paste("p= ", signif(for.pval, digits=3), sep="")) +
    scale_y_continuous(expand = expansion(mult = c(0, .1)))

  
  print(g)
  dev.off()
}



#### save histograms as a pdf ####

all_images <- list.files("./cluster_runs/pre_DEBIAS-M_100_perm/output/by_hand", full.names = TRUE)
ind <- grep("post_DEBIAS-M_RF_lognorm_histogram_", all_images)
all_images <- all_images[ind]

plots <- lapply(ll <- all_images, function(x){
  img <- as.raster(readPNG(x))
  rasterGrob(img, interpolate = FALSE)
})

ggsave("./output/test_permutated/hist_100_perm_post_by_hand.pdf", marrangeGrob(grobs=plots, nrow=2, ncol=2), width = 4, height = 4, dpi = 300)

#### box plots ####

selected_AUCs <- AUCs %>% filter(Permutation == 0) %>% merge(y = phenos, by.x = "Study_ID", by.y = "ID", all.x = TRUE) 
selected_AUCs$phenotype[grep("cancer", selected_AUCs$phenotype)] <- "cancer"

auc_b <- ggplot(selected_AUCs, mapping = aes(x=DEBIAS, y=AUC)) + 
  geom_boxplot() + 
  geom_point(aes(col = as.factor(phenotype)), alpha = 2) +
  scale_color_brewer(name = "pehnotype", palette = "Paired") +
  scale_x_discrete(labels = c("Before DEBIAS-M", "After DEBIAS-M")) +
  labs(x = "", title = "AUCs before and after DEBIAS-M") +
  theme_dark()

auc_b

pval.df <- merge(pval.df, phenos, by.x = "Study_ID", by.y = "ID", all.x = TRUE)
pval.df$phenotype[grep("cancer", pval.df$phenotype)] <- "cancer"
# write.csv(pval.df, "./csv_files/DEBIAS-M_AUC_pvals.csv")

pval_b <- ggplot(pval.df, mapping = aes(x = as.factor(DEBIAS), y = pval)) + 
  geom_boxplot() +
  geom_point(aes(col = as.factor(phenotype))) +
  scale_color_brewer(name = "phenotype", palette = "Paired") +
  scale_x_discrete(labels = c("Before DEBIAS-M", "After DEBIAS-M")) +
  labs(x = "", title = "pvalues before and after DEBIAS-M") +
  theme_dark()

pval_b

png("./output/test_permutated/box_AUC_pval.png", width = 1000, height = 480)
grid.arrange(auc_b, pval_b, nrow = 1)
dev.off()


#### pcoa ####

varespec.bray <- vegdist(varespec, method = "bray")




