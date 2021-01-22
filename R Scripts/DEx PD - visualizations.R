# some visualizations

#### Plotting results
# MA plot
par(mfrow=c(1,3), mar=c(4,4,2,1))
plotMA(res.DE.mRNA.new.ordered,main = "mRNA mean count")
plotMA(res.DE.mRNA.new_shrink_apeglm,main = "shrinked by apeglm")
plotMA(res.DE.mRNA.new_shrink_ashr,main = "shrinked by ashr")

# Volcano plot
# -> p-value<0.1
alpha <- 0.05 # Threshold on the adjusted p-value
volcano_dist_plot(res.DE.mRNA.new,alpha,"Volcano plot mRNA with p-value<0.1")

# Counts transformation into log2-like & PCA analysis for default
# -> VSD
vsd.mRNA <- varianceStabilizingTransformation(DE.mRNA.new,blind = FALSE)
normed_mRNA_VSD <- assay(vsd.mRNA)
meanSdPlot(normed_mRNA_VSD, xlab = "Mean in normalized mRNA VSD") 
dists_mRNA_VSD <- dist(t(normed_mRNA_VSD))
plot_dist_mRNA_VSD <- plot(hclust(dists_mRNA_VSD))
# -> PCA
pca_mRNA <- plotPCA(vsd.mRNA)
pca_data <- pca_mRNA$data
# example of how to label the outlier
idx_control_remove <- which(pca_data$PC1 == max(pca_data$PC1))
point <- pca_data[idx_control_remove,]
ggplot(pca_data, aes(x=PC1,y=PC2,color = condition)) + geom_point() + 
  geom_point(data = point, aes(x=PC1,y=PC2), colour = "red",size = 4) +
  annotate("text",label = point$name, x= point$PC1, y = point$PC2) +
  ggtitle("PCA mRNA by VSD")

# Comparison between GO analysis of frameworks A and B
enrichment_GO_common_UP <- intersect(enrichment_GO_B$`Functional Category`,
                                     enrichment_GO_A$`Functional Category`)
enrichment_GO_only_B <- enrichment_GO_B[!enrichment_GO_B$`Functional Category` %in% enrichment_GO_common_UP,]
enrichment_GO_only_A <- enrichment_GO_A[!enrichment_GO_A$`Functional Category` %in% enrichment_GO_common_UP,]

enrichment_GO_A_B <- merge(enrichment_GO_A,
                           enrichment_GO_B, 
                           by.x = "Functional Category",
                           by.y = "Functional Category", all = TRUE)
enrichment_GO_A_B$`Enrichment FDR.x.log10` <- -log10(enrichment_GO_A_B$`Enrichment FDR.x`)
enrichment_GO_A_B$`Enrichment FDR.y.log10` <- -log10(enrichment_GO_A_B$`Enrichment FDR.y`)
enrichment_GO_A_B <- enrichment_GO_A_B %>%
  select(`Functional Category`,`Enrichment FDR.x.log10`,`Enrichment FDR.y.log10`) %>%
  distinct()
enrichment_GO_A_B_new <- enrichment_GO_A_B
enrichment_GO_A_B_new$`Enrichment FDR.x.log10` <- ifelse(is.na(enrichment_GO_A_B_new$`Enrichment FDR.x.log10`),
                                                         0,
                                                         enrichment_GO_A_B_new$`Enrichment FDR.x.log10`)
enrichment_GO_A_B_new$`Enrichment FDR.y.log10` <- ifelse(is.na(enrichment_GO_A_B_new$`Enrichment FDR.y.log10`),
                                                         0,
                                                         enrichment_GO_A_B_new$`Enrichment FDR.y.log10`)
enrichment_GO_A_B_new$label <- ifelse(enrichment_GO_A_B_new$`Enrichment FDR.x.log10` > 8.7 &
                                        enrichment_GO_A_B_new$`Enrichment FDR.y.log10` > 3.3,
                                      "both",
                                      ifelse(enrichment_GO_A_B_new$`Enrichment FDR.y.log10` > 3.3,
                                             "B",
                                             ifelse(enrichment_GO_A_B_new$`Enrichment FDR.x.log10` > 8.7,
                                                    "A","no")))

ggplot(enrichment_GO_A_B_new,aes(x=`Enrichment FDR.x.log10`,
                                 y=`Enrichment FDR.y.log10`)) +
  geom_point(colour = "grey",size = 3,alpha = 0.75) +  
  geom_abline(colour = "grey",alpha = 0.8,type="l", lty=2) +
  geom_point(data=subset(enrichment_GO_A_B_new,label=="both"),
             aes(colour = "Both"),size = 3,alpha = 0.7) +
  geom_point(data=subset(enrichment_GO_A_B_new,label=="A"),
             aes(colour = "Framework A"),size = 3,alpha = 0.7) +
  geom_point(data=subset(enrichment_GO_A_B_new,label=="B"),
             aes(colour = "Framework B"),size = 3,alpha = 0.7) +
  geom_text_repel(data=subset(enrichment_GO_A_B_new,label!="no"),
                  aes(label=`Functional Category`),size = 3.2,
                  min.segment.length = unit(0.2, "lines")) +
  scale_color_manual(values = c("#3F5DF5","#EDB64E","#6AD082"),
                     breaks = c("Framework A","Framework B","Both")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),plot.background = element_blank()) + 
  labs(title = "GO Enrichment Framework A vs Framework B",
       x = "-log10(FDR) - Framework A",
       y = "-log10(FDR) - Framework B",
       color = "Top Enriched \nGO Categories") + 
  xlim(0,max(enrichment_GO_A_B_new$`Enrichment FDR.x.log10`)) + 
  ylim(0,max(enrichment_GO_A_B_new$`Enrichment FDR.y.log10`)) +
  theme(aspect.ratio = 1) +theme_classic()

# density plots to put on the side
density_plot_A <- enrichment_GO_A_B[!is.na(enrichment_GO_A_B$`Enrichment FDR.x.log10`),] %>%
  mutate(FDR = `Enrichment FDR.x.log10`,type = rep("A",length(which(!is.na(enrichment_GO_A_B$`Enrichment FDR.x.log10`))))) %>%
  select(FDR,type)

density_plot_B <- enrichment_GO_A_B[!is.na(enrichment_GO_A_B$`Enrichment FDR.y.log10`),] %>%
  mutate(FDR = `Enrichment FDR.y.log10`,type = rep("B",length(which(!is.na(enrichment_GO_A_B$`Enrichment FDR.y.log10`))))) %>%
  select(FDR,type)

density_plot_both <- rbind(density_plot_A,density_plot_B)
density_plot_both$type <- rep("both",nrow(density_plot_both))
density_plot_A_new <- rbind(density_plot_A,density_plot_both)
density_plot_B_new <- rbind(density_plot_B,density_plot_both)
density_plot_B_new <- density_plot_B_new[density_plot_B_new$FDR<max(enrichment_GO_A_B_new$`Enrichment FDR.y.log10`),]

ggplot(density_plot_A_new, aes(FDR, fill = type)) +
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = c("#3F5DF5","#6AD082"),breaks = c("A","both"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),plot.background = element_blank()) + 
  labs(title = "density framework A and both",
       x = "",
       y = "",
       color = "") +
  xlim(0,max(enrichment_GO_A_B_new$`Enrichment FDR.x.log10`)) + 
  ylim(0,1) +
  theme_classic()

ggplot(density_plot_B_new, aes(FDR, fill = type)) +
  geom_density(alpha = 0.6) + 
  scale_fill_manual(values = c("#EDB64E","#6AD082"),breaks = c("B","both"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),plot.background = element_blank()) + 
  labs(title = "density framework B and both",
       x = "",
       y = "",
       color = "") +
  xlim(0,5.2) + 
  ylim(0,1.3) +
  theme_classic()