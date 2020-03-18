library(rfPermute)
library(ggplot2)
library(vegan)

cf146<-read.csv("CF146 all data for stats.csv")
d<-as.matrix(cf146[,2:19])
rownames(d)<-cf146$Day


nmds1<-metaMDS(d, distance = "bray", autotransform =TRUE)

stressplot(nmds1)

plot(nmds1, type = "t")

data.scores<- as.data.frame(scores(nmds1))
data.scores$site <- row.names(data.scores)
data.scores$grp <- c("Colistin", "Colistin","Colistin", "No antibiotics", "No antibiotics", "No antibiotics", "No antibiotics", "Clindamycin", "Clindamycin")
head(data.scores)
species.scores <- as.data.frame(scores(nmds1, "species"))
species.scores$species <- rownames(species.scores)

plot1<-ggplot()+
  geom_text(data = species.scores, aes (x = NMDS1, y = NMDS2, label = species), alpha = 0.4, size = 1.25) +
  geom_point (data = data.scores, aes (x = NMDS1, y = NMDS2, shape = grp, colour = grp), size = 3) +
  geom_text(data = data.scores, aes(x = NMDS1, y = NMDS2, label = site), size = 1, vjust = 1)+
  coord_equal() +
  theme_bw() +
  theme(axis.line = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor=element_blank(), 
        panel.border = element_rect(colour="black"),
        axis.text=element_text(size=6),
        axis.title=element_text(size=6),
        #axis.title.y = element_text(margin = margin(t = 0, r = 2, b = 0, l = 0)),
        #axis.title.x = element_text("Days since hospital admission and IV antibiotic treatment"),
        legend.position=c(0.2,0.2),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.text=element_text(size=2),
        legend.key.size = unit(0.1, 'lines'),
        legend.spacing.x = unit(0.1, 'cm'))+
  ylim(c(-0.7, 0.6))+
  xlim(c(-1.1, 0.6))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggsave("NMDS_antibiotics.pdf",width = 8, height = 8, units = "cm", dpi = 300)
plot1


#Unsupervised random forest
library(randomForest)
set.seed(100)
rf1<-randomForest(d, proximity = TRUE, ntree=5000, mtry = 3)
rf1$votes
summary(rf1)
rf1

omics.pc <- prcomp(d[,1:18], center = FALSE, scale. = FALSE)
km.cluster <- kmeans(d[,1:18], centers = 3, iter.max = 20, nstart = 2)
omics.pc$kmeans.cluster <- km.cluster$cluster
hclust.rf <- hclust(as.dist(1-rf1$proximity), method = "ward.D2")
rf.cluster = cutree(hclust.rf, k=3)
omics.pc$rf.clusters <- rf.cluster
table(rf.cluster, cf146$Day)


#Random forest supervised by antibiotic treatment
grp<-c("No antibiotics", "Colistin","Colistin", "No antibiotics", "No antibiotics", "No antibiotics", "No antibiotics", "Clindamycin", "Clindamycin")
rf4<-rfPermute(factor(data.scores$grp)~., scale = TRUE, proximity = TRUE, data = d[,1:18], importance = TRUE, ntree = 5000, nrep = 10000)
rf4
plot(rp.importance(rf4, scale = TRUE))
proximityPlot(rf4)
confusionMatrix(rf4)
impHeatmap(rf4)
rf4$pval


######
#Plot of transcriptome data by gene normalized by transcripts per ml
genes<-read.csv("mt_genes_norm.csv")
library(reshape2)
mgenes<-melt(genes[,2:9], id.vars="Time", variable.name="Gene", value.name="Transc_per_ml")
plot6<-ggplot(mgenes, aes(x=Time, y = log10(Transc_per_ml), colour = Gene))+
  #geom_boxplot(lwd=0.3, outlier.shape = NA) +
  annotate("rect", xmin = 10, xmax = 27, ymin = -Inf, ymax = Inf,
           alpha = .2)+
  annotate("rect", xmin = 10, xmax = 13, ymin = -Inf, ymax = Inf,
           alpha = .3)+
  annotate("rect", xmin = 36, xmax = 40, ymin = -Inf, ymax = Inf,
           alpha = .2)+
  geom_point(size = 4)+
  geom_line(size = 1.5)+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor=element_blank(), 
        panel.border = element_blank(),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10),     
        axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position=c(0.2,0.7),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.text=element_text(size=8),
        legend.key.size = unit(0.3, 'lines'),
        legend.spacing.x = unit(0.2, 'cm'))+
  xlim(c(10,40))+
  ylab("Trancription (transcripts/ml)")+
  xlab("Days since hospital admission")+
  scale_colour_manual(limits = c("budA", "budB","budC"), values = c("black", "darkslategray", "slategray"))+
  ggsave("Trancrip_per_ml_genes_bud.eps",width = 8, height = 8, units = "cm", dpi = 300)+
  ggsave("Trancrip_per_ml_genes_bud.pdf",width = 8, height = 8, units = "cm", dpi = 300)

plot6

plot7<-ggplot(mgenes, aes(x=Time, y = log10(Transc_per_ml), colour = Gene))+
  annotate("rect", xmin = 10, xmax = 27, ymin = -Inf, ymax = Inf,
           alpha = .2)+
  annotate("rect", xmin = 10, xmax = 13, ymin = -Inf, ymax = Inf,
           alpha = .3)+
  annotate("rect", xmin = 36, xmax = 40, ymin = -Inf, ymax = Inf,
           alpha = .2)+
  geom_point(size = 4)+
  geom_line(size = 1.5)+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor=element_blank(), 
        panel.border = element_blank(),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        #axis.title.y = element_text(margin = margin(t = 0, r = 2, b = 0, l = 0)),
        #axis.title.x = element_text("Days since hospital admission and IV antibiotic treatment"),
        axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position=c(0.2,0.7),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.text=element_text(size=8),
        legend.key.size = unit(0.3, 'lines'),
        legend.spacing.x = unit(0.2, 'cm'))+
  xlim(c(10,40))+
  ylab("Trancription (transcripts/ml)")+
  xlab("Days since hospital admission")+
  scale_colour_manual(limits = c("acoA","acoB","Phenazine"), values = c("black", "darkslategray", "slategray" ))+
  ggsave("Trancrip_per_ml_genes_aco.eps",width = 8, height = 8, units = "cm", dpi = 300)+
  ggsave("Trancrip_per_ml_genes_aco.pdf",width = 8, height = 8, units = "cm", dpi = 300)

plot7


#####
#Plot clinical data

fev<-read.csv("CF146_fev.csv")

plot8<-ggplot(fev, aes(x=Day, y = FEV1p...))+
  annotate("rect", xmin = 0, xmax = 90, ymin = -Inf, ymax = Inf,
           alpha = .2)+
  geom_point(size = 2)+
  geom_vline(xintercept = 0, linetype = "dashed")+
  geom_vline(xintercept = 27, linetype = "dashed")+
  geom_vline(xintercept = 36, linetype = "dashed")+
  geom_vline(xintercept = 54, linetype = "dashed")+
  stat_smooth(method = "lm", formula = y ~ splines::bs(x,27))+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor=element_blank(), 
        panel.border = element_blank(),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10),
        axis.text.x=element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position=c(0.2,0.7),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.text=element_text(size=8),
        legend.key.size = unit(0.3, 'lines'),
        legend.spacing.x = unit(0.2, 'cm'))+
  ylab("FEV1p(%)")+
  xlab("Days since hospital admission ")+
  ggsave("FEVb.eps",width = 16, height = 8, units = "cm", dpi = 300)+
  ggsave("FEVb.pdf",width = 16, height = 8, units = "cm", dpi = 300)

plot8

library(fields)
tps1<-Tps(fev$FEV1p..., fev$Day)
summary(tps1)
plot(tps1)
surface(tps1)
tps1$eff.df
tps1$best.model

reg1<-qsreg(fev$Day, fev$FEV1p...)
summary(reg1)
plot(reg1)
predict(reg1)
reg1_25<-qsreg(fev$Day, fev$FEV1p..., alpha = 0.25)
summary(reg1_25)
reg1_75<-qsreg(fev$Day, fev$FEV1p..., alpha = 0.75)
summary(reg1_75)
plot(reg1_75)
