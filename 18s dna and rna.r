

setwd("C:/Users/XZM/Desktop/Coldseep 2nd study/DNA RNA togather analysis/R codes and files")


#  18S DNA-------------------------------------------------------

library(tidyverse)
library(vegan)

asv_18s = read.table("asv_18s.txt",header = T,check.names = F) %>% t()
tax_18s = read.delim("taxonomy-pr2.tsv")


## remove metazoan and land plants
tax_pro = tax_18s[!c(str_detect(tax_18s$Taxon, "Metazoa")),]
tax_pro = tax_pro[!c(str_detect(tax_pro$Taxon, "Streptophyta")),]
asv_18s = asv_18s[,match(tax_pro$Feature.ID, colnames(asv_18s))]
asv_18s = rrarefy(asv_18s,min(rowSums(asv_18s)))

asv18s_dna = asv_18s[asv_18s %>% rownames() %>% str_detect("dna"), ]
asv18s_dna %>% rownames()
asv18s_dna = asv18s_dna[,which(colSums(asv18s_dna)>0)]
asv18s_dna %>% dim
asv18s_dna = rrarefy(asv18s_dna, min(rowSums(asv18s_dna)))
asv18s_dna = asv18s_dna[,which(colSums(asv18s_dna)>0)]

grp1 <-
  asv18s_dna %>% rownames() %>% as_tibble() %>% 
  rename(sample.id = value) %>% 
  mutate(ROV = sample.id %>% 
           str_split("_") %>% sapply('[', 1) %>% str_sub(5,8)) %>% 
  mutate(Depth_cmbs = sample.id %>% 
           str_split("_") %>% sapply('[', 2)) %>% 
  mutate(Ribosome = sample.id %>% 
           str_split("_") %>% sapply('[', 3)) %>% 
  mutate(Nucleic = sample.id %>% 
           str_split("_") %>% sapply('[', 4))

## NMDS and 3D plot
d = vegdist(asv18s_dna %>% sqrt)
nmds = monoMDS(d,k=3)
NMDS = data.frame(MDS1 = nmds$points[,1], MDS2 = nmds$points[,2], MDS3=nmds$points[,3])
NMDS$ROV = as.factor(grp1$ROV)
#NMDS$ROV = as.factor(grp1)
library(rgl)
library(car)
scatter3d(x = NMDS$MDS1, y = NMDS$MDS2, z = NMDS$MDS3, groups = NMDS$ROV,
          surface=FALSE, ellipsoid = TRUE,grid = FALSE,
          xlab = "NMDS1", ylab = "NMDS2", zlab = "NMDS3",
          surface.col = c("red", "hotpink", "darkgoldenrod1", "cornflowerblue","darkcyan"))
rgl.snapshot("NMDS-18s DNA community.png", "png")


anosim(d, grp1$Depth_cmbs)
anosim(d, grp1$ROV)

# 18S RNA -----------------------------------------------------------------


asv18s_rna = asv_table_18s[asv_table_18s %>% rownames() %>% str_detect("rna"), ]
asv18s_rna = asv18s_rna[asv18s_rna %>% rownames() %>% str_detect("_15_", negate = T), ]# remove a few samples with depth = 15cm
asv18s_rna %>% rownames()
asv18s_rna = asv18s_rna[,which(colSums(asv18s_rna)>0)]
asv18s_rna %>% dim
asv18s_rna = rrarefy(asv18s_rna, min(rowSums(asv18s_rna))) 
asv18s_rna = asv18s_rna[,which(colSums(asv18s_rna)>0)]

asv18s_rna_alpha = alpha(asv18s_rna, tree_18s)

grp2 <-
  asv18s_rna %>% rownames() %>% as_tibble() %>% 
  rename(sample.id = value) %>% 
  mutate(ROV = sample.id %>% 
           str_split("_") %>% sapply('[', 1) %>% str_sub(5,8)) %>% 
  mutate(Depth_cmbs = sample.id %>% 
           str_split("_") %>% sapply('[', 2)) %>% 
  mutate(Ribosome = sample.id %>% 
           str_split("_") %>% sapply('[', 3)) %>% 
  mutate(Nucleic = sample.id %>% 
           str_split("_") %>% sapply('[', 4))


d2 = vegdist(asv18s_rna)
nmds = monoMDS(d2)
NMDS = data.frame(MDS1 = nmds$points[,1], MDS2 = nmds$points[,2])


NMDS %>% as_tibble() %>% 
  mutate(Depth_cmbs = factor(grp2$Depth_cmbs, levels = c("0", "5", "10")),
         ROV = grp2$ROV %>% as.factor) %>%
  ggplot(aes(x = MDS1,y = MDS2))+theme_bw()+
  geom_point(aes(shape = Depth_cmbs, color = ROV),size=5)+theme(panel.grid=element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=20),
        axis.text.y = element_text(color="black", size=20))+
  theme(legend.title=element_text(size=20,face = "bold"), legend.text=element_text(size=20,face = "bold"),
        axis.title=element_text(size=20))+
  scale_color_manual(values = c("firebrick1", "hotpink", "darkgoldenrod1", "cornflowerblue","darkcyan"))+
  ggtitle("18S rRNA community")

anosim(d2, grp2$Depth_cmbs)# ANOSIM statistic R: 0.07666; Significance: 0.006 
anosim(d2, grp2$ROV) #ANOSIM statistic R: 0.4114; Significance: 0.001 

d = vegdist(asv18s_rna %>% sqrt)
nmds = monoMDS(d,k=3)
NMDS = data.frame(MDS1 = nmds$points[,1], MDS2 = nmds$points[,2], MDS3=nmds$points[,3])
NMDS$ROV = as.factor(grp2$ROV)
#NMDS$ROV = as.factor(grp1)
library(rgl)
library(car)
scatter3d(x = NMDS1$MDS1, y = NMDS1$MDS2, z = NMDS1$MDS3, groups = NMDS1$ROV,
          surface=FALSE, ellipsoid = TRUE,grid = FALSE,
          xlab = "NMDS1", ylab = "NMDS2", zlab = "NMDS3",
          surface.col = c("red", "hotpink", "darkgoldenrod1", "cornflowerblue","darkcyan"))
rgl.snapshot("NMDS-18s RNA community.png", "png")



# 16S DNA---------------------------------------------------------------------
asv_16s = read.table("asv_table_16s.txt",header = T,check.names = F) %>% t()
asv_16s %>% dim

asv16s_dna = asv_16s[asv_16s %>% rownames() %>% str_detect("dna"), ]
asv16s_dna %>% rownames()
asv16s_dna = asv16s_dna[,which(colSums(asv16s_dna)>0)]
asv16s_dna %>% dim
asv16s_dna = rrarefy(asv16s_dna, min(rowSums(asv16s_dna)))
asv16s_dna = asv16s_dna[,which(colSums(asv16s_dna)>0)]

grp1 <-
  asv16s_dna %>% rownames() %>% as_tibble() %>% 
  rename(sample.id = value) %>% 
  mutate(ROV = sample.id %>% 
           str_split("_") %>% sapply('[', 1) %>% str_sub(5,8)) %>% 
  mutate(Depth_cmbs = sample.id %>% 
           str_split("_") %>% sapply('[', 2)) %>% 
  mutate(Ribosome = sample.id %>% 
           str_split("_") %>% sapply('[', 3)) %>% 
  mutate(Nucleic = sample.id %>% 
           str_split("_") %>% sapply('[', 4))

## NMDS and 3D plot
d1 = vegdist(asv16s_dna %>% sqrt)
nmds1 = monoMDS(d1,k=3)
NMDS1 = data.frame(MDS1 = nmds1$points[,1], MDS2 = nmds1$points[,2], MDS3=nmds1$points[,3])
NMDS1$ROV = as.factor(grp1$ROV)
#NMDS$ROV = as.factor(grp1)
library(rgl)
library(car)
scatter3d(x = NMDS1$MDS1, y = NMDS1$MDS2, z = NMDS1$MDS3, groups = NMDS1$ROV,
          surface=FALSE, ellipsoid = TRUE,grid = FALSE,
          xlab = "NMDS1", ylab = "NMDS2", zlab = "NMDS3",
          surface.col = c("red", "hotpink", "darkgoldenrod1", "cornflowerblue","darkcyan"))
rgl.snapshot("NMDS-16s DNA community.png", "png")


anosim(d1, grp1$Depth_cmbs)
anosim(d1, grp1$ROV)

# 16S RNA -----------------------------------------------------------------
asv16s_rna = asv_16s[asv_16s %>% rownames() %>% str_detect("rna"), ]
asv16s_rna %>% rownames()
asv16s_rna = asv16s_rna[,which(colSums(asv16s_rna)>0)]
asv16s_rna %>% dim
asv16s_rna = rrarefy(asv16s_rna, min(rowSums(asv16s_rna)))
asv16s_rna = asv16s_rna[,which(colSums(asv16s_rna)>0)]

grp1 <-
  asv16s_rna %>% rownames() %>% as_tibble() %>% 
  rename(sample.id = value) %>% 
  mutate(ROV = sample.id %>% 
           str_split("_") %>% sapply('[', 1) %>% str_sub(5,8)) %>% 
  mutate(Depth_cmbs = sample.id %>% 
           str_split("_") %>% sapply('[', 2)) %>% 
  mutate(Ribosome = sample.id %>% 
           str_split("_") %>% sapply('[', 3)) %>% 
  mutate(Nucleic = sample.id %>% 
           str_split("_") %>% sapply('[', 4))

## NMDS and 3D plot
d1 = vegdist(asv16s_rna %>% sqrt)
nmds1 = monoMDS(d1,k=3)
NMDS1 = data.frame(MDS1 = nmds1$points[,1], MDS2 = nmds1$points[,2], MDS3=nmds1$points[,3])
NMDS1$ROV = as.factor(grp1$ROV)
#NMDS$ROV = as.factor(grp1)
library(rgl)
library(car)
scatter3d(x = NMDS1$MDS1, y = NMDS1$MDS2, z = NMDS1$MDS3, groups = NMDS1$ROV,
          surface=FALSE, ellipsoid = TRUE,grid = FALSE,
          xlab = "NMDS1", ylab = "NMDS2", zlab = "NMDS3",
          surface.col = c("red", "hotpink", "darkgoldenrod1", "cornflowerblue","darkcyan"))
rgl.snapshot("NMDS-16s RNA community.png", "png")


anosim(d1, grp1$Depth_cmbs)
anosim(d1, grp1$ROV)

# d1 = vegdist(asv18s_dna)
# nmds = monoMDS(d1)
# NMDS = data.frame(MDS1 = nmds$points[,1], MDS2 = nmds$points[,2])
# NMDS %>% as_tibble() %>% 
#   mutate(Depth_cmbs = grp1$Depth_cmbs %>% as.factor,
#          ROV = grp1$ROV %>% as.factor) %>%
#   ggplot(aes(x = MDS1,y = MDS2))+theme_bw()+
#   geom_point(aes(color = ROV),size=5)+theme(panel.grid=element_blank())+
#   theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
#   theme(axis.text.x = element_text(color="black", size=15),
#         axis.text.y = element_text(color="black", size=15))+
#   theme(legend.title=element_text(size=15,face = "bold"), legend.text=element_text(size=12,face = "bold"),
#         axis.title=element_text(size=15))+
#   scale_color_manual(values = c("firebrick1", "hotpink", "darkgoldenrod1", "cornflowerblue","darkcyan"))


# alpha diversity -----------------------------------------------------
# define alpha() function
alpha <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  Richness <- est[1, ]
  Chao1 <- est[2, ]
  ACE <- est[4, ]
  Shannon <- diversity(x, index = 'shannon', base = base)
  Simpson <- diversity(x, index = 'simpson')    #Gini-Simpson
  Pielou <- Shannon / log(Richness, base)
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  result <- data.frame(Richness, Shannon, Simpson, Pielou, Chao1, ACE, goods_coverage)
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
  }
  result
}

library(picante)




# depth vs diversity -------------------------------------------------------

# alpha vs depth   #ALT + L
library(picante)
tree_18s = read.tree("rooted_tree_18s.nwk")
tree_16s = read.tree("tree-16s-rooted.nwk")


#18s-dna

asv18s_dna = asv_table_16s[rownames(asv_table_16s) %>% str_detect("dna"),]
asv18s_dna_alpha = alpha(asv18s_dna)



asv18s_dna_alpha %>% 
  rownames_to_column(., var = "sample_id") %>% 
  as_tibble() %>%  
  mutate(ROV = sample_id %>% str_sub(5,8)) %>%
  mutate(Depth_cmbs = sample_id %>%
           str_split("_") %>% sapply('[', 2) %>% 
            as.numeric()) %>% 
  # filter(Depth_cmbs < 30) %>% 
  ggplot(aes(x = Depth_cmbs, y = Richness)) +
  geom_point(aes(color = ROV), size = 3)+
 # geom_boxplot()+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12))+
  theme(legend.title=element_text(size=15), legend.text=element_text(size=12),
        axis.title=element_text(size=15))+
  geom_smooth(linewidth = 2, method = "lm")+
  scale_color_manual(values = c("firebrick1", "hotpink", "darkgoldenrod1", "cornflowerblue","darkcyan"))+
  ggtitle("18S DNA community") -> r1

asv18s_dna_alpha %>% 
  rownames_to_column(., var = "sample_id") %>% 
  as_tibble() %>%  
  mutate(ROV = sample_id %>% str_sub(5,8)) %>%
  mutate(Depth_cmbs = sample_id %>%
           str_split("_") %>% sapply('[', 2) %>% 
           as.numeric()) %>% 
  ggplot(aes(x = Depth_cmbs, y = Chao1)) +
  geom_point(aes(color = ROV), size = 3)+
  # geom_boxplot()+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12))+
  theme(legend.title=element_text(size=15), legend.text=element_text(size=12),
        axis.title=element_text(size=15))+
  geom_smooth(linewidth = 2)+
  scale_color_manual(values = c("firebrick1", "hotpink", "darkgoldenrod1", "cornflowerblue","darkcyan"))+
  ggtitle("18S DNA community") -> r2

asv18s_dna_alpha %>% 
  rownames_to_column(., var = "sample_id") %>% 
  as_tibble() %>%  
  mutate(ROV = sample_id %>% str_sub(5,8)) %>%
  mutate(Depth_cmbs = sample_id %>%
           str_split("_") %>% sapply('[', 2) %>% 
           as.numeric()) %>% 
  ggplot(aes(x = Depth_cmbs, y = PD_whole_tree)) +
  geom_point(aes(color = ROV), size = 3)+
  # geom_boxplot()+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12))+
  theme(legend.title=element_text(size=15), legend.text=element_text(size=12),
        axis.title=element_text(size=15))+
  geom_smooth(linewidth = 2)+
  scale_color_manual(values = c("firebrick1", "hotpink", "darkgoldenrod1", "cornflowerblue","darkcyan"))+
  ggtitle("18S DNA community") -> r3


#16s-dna
tree_16s = read.tree("tree-16s-rooted.nwk")
asv16s_dna_alpha = alpha(asv16s_dna, tree_16s)

asv16s_dna_alpha %>% 
  rownames_to_column(., var = "sample_id") %>% 
  as_tibble() %>%  
  mutate(ROV = sample_id %>% str_sub(5,8)) %>%
  mutate(Depth_cmbs = sample_id %>%
           str_split("_") %>% sapply('[', 2) %>% as.numeric()) %>% 
  ggplot(aes(x = Depth_cmbs, y = Richness)) +
  geom_point(aes(color = ROV), size = 3)+
  # geom_boxplot()+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12))+
  theme(legend.title=element_text(size=15), legend.text=element_text(size=12),
        axis.title=element_text(size=15))+
  geom_smooth(linewidth = 2)+
  scale_color_manual(values = c("firebrick1", "hotpink", "darkgoldenrod1", "cornflowerblue","darkcyan"))+
  ggtitle("16S DNA community") -> r4

asv16s_dna_alpha %>% 
  rownames_to_column(., var = "sample_id") %>% 
  as_tibble() %>%  
  mutate(ROV = sample_id %>% str_sub(5,8)) %>%
  mutate(Depth_cmbs = sample_id %>%
           str_split("_") %>% sapply('[', 2) %>% as.numeric()) %>% 
  ggplot(aes(x = Depth_cmbs, y = Chao1)) +
  geom_point(aes(color = ROV), size = 3)+
  # geom_boxplot()+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12))+
  theme(legend.title=element_text(size=15), legend.text=element_text(size=12),
        axis.title=element_text(size=15))+
  geom_smooth(linewidth = 2)+
  scale_color_manual(values = c("firebrick1", "hotpink", "darkgoldenrod1", "cornflowerblue","darkcyan"))+
  ggtitle("16S DNA community") ->r5


asv16s_dna_alpha %>% 
  rownames_to_column(., var = "sample_id") %>% 
  as_tibble() %>%  
  mutate(ROV = sample_id %>% str_sub(5,8)) %>%
  mutate(Depth_cmbs = sample_id %>%
           str_split("_") %>% sapply('[', 2) %>% as.numeric()) %>% 
  ggplot(aes(x = Depth_cmbs, y = PD_whole_tree)) +
  geom_point(aes(color = ROV), size = 3)+
  # geom_boxplot()+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12))+
  theme(legend.title=element_text(size=15), legend.text=element_text(size=12),
        axis.title=element_text(size=15))+
  geom_smooth(linewidth = 2)+
  scale_color_manual(values = c("firebrick1", "hotpink", "darkgoldenrod1", "cornflowerblue","darkcyan"))+
  ggtitle("16S DNA community") -> r6

library(ggpubr)
ggarrange(r1,r4,r2,r5,r3,r6, ncol = 2, nrow = 3)
##compare alpha diversity [0-10cm]

#18S DNA
asv_18s_apha = alpha(asv_18s, tree_18s)
asv_16s_apha = alpha(asv_16s, tree_16s)


r1<-
asv_18s_apha %>% 
  rownames_to_column(., var = "sample_id") %>% 
  as_tibble() %>%  
  mutate(ROV = sample_id %>% str_sub(5,8)) %>%
  mutate(Depth_cmbs = sample_id %>%
           str_split("_") %>% sapply('[', 2) %>% 
           str_replace("^5$","05")) %>%
  mutate(Ribosome = sample_id %>%
           str_split("_") %>% sapply('[', 3)) %>%
  mutate(Nucleic = sample_id %>%
           str_split("_") %>% sapply('[', 4)) %>% 
  filter(Nucleic == "dna" &
           Depth_cmbs < 15) %>% 
  ggplot(aes(x = ROV, y = Shannon)) +
  geom_boxplot(outlier.colour = "white", outlier.fill = "white", linewidth =2, color = "grey66")+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
             size = 4,aes(color= Depth_cmbs), alpha = 0.7,show.legend = T)+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=12,face = "bold"),
        axis.text.y = element_text(color="black", size=12,face = "bold"))+
  #annotate("text",label="ANOSIM_transplant: R = 0.19, p = 0.027",x=-0.72,y=0.58, size=5,fontface="bold")+
  #annotate("text",label="ANOSIM_tissue: R = 0.38, p = 0.001",x=-0.72,y=0.50, size=5,fontface="bold")+
  theme(legend.title=element_text(size=15,face = "bold"), legend.text=element_text(size=12,face = "bold"),
        axis.title=element_text(size=15,face = "bold"))+
  ggtitle("18S DNA")+
  xlab("")
r1

r2<-
asv_18s_apha %>% 
  rownames_to_column(., var = "sample_id") %>% 
  as_tibble() %>%  
  mutate(ROV = sample_id %>% str_sub(5,8)) %>%
  mutate(Depth_cmbs = sample_id %>%
           str_split("_") %>% sapply('[', 2) %>% 
           str_replace("^5$","05")) %>%
  mutate(Ribosome = sample_id %>%
           str_split("_") %>% sapply('[', 3)) %>%
  mutate(Nucleic = sample_id %>%
           str_split("_") %>% sapply('[', 4)) %>% 
  filter(Nucleic == "dna" &
           Depth_cmbs < 15) %>% 
  ggplot(aes(x = ROV, y  = Chao1)) +
  geom_boxplot(outlier.colour = "white", outlier.fill = "white", size =2, color = "grey66")+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
             size = 4,aes(color= Depth_cmbs), alpha = 0.7, show.legend = T)+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=12,face = "bold"),
        axis.text.y = element_text(color="black", size=12,face = "bold"))+
  #annotate("text",label="ANOSIM_transplant: R = 0.19, p = 0.027",x=-0.72,y=0.58, size=5,fontface="bold")+
  #annotate("text",label="ANOSIM_tissue: R = 0.38, p = 0.001",x=-0.72,y=0.50, size=5,fontface="bold")+
  theme(legend.title=element_text(size=15,face = "bold"), legend.text=element_text(size=12,face = "bold"),
        axis.title=element_text(size=15,face = "bold"))+
  ggtitle("18S DNA")

r3<-
asv_18s_apha %>% 
  rownames_to_column(., var = "sample_id") %>% 
  as_tibble() %>%  
  mutate(ROV = sample_id %>% str_sub(5,8)) %>%
  mutate(Depth_cmbs = sample_id %>%
           str_split("_") %>% sapply('[', 2) %>% 
           str_replace("^5$","05")) %>%
  mutate(Ribosome = sample_id %>%
           str_split("_") %>% sapply('[', 3)) %>%
  mutate(Nucleic = sample_id %>%
           str_split("_") %>% sapply('[', 4)) %>% 
  filter(Nucleic == "dna" &
           Depth_cmbs < 15) %>% 
  ggplot(aes(x = ROV, y = PD_whole_tree)) +
  geom_boxplot(outlier.colour = "white", outlier.fill = "white", size =2, color = "grey66")+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
             size = 4,aes(color= Depth_cmbs),  alpha = 0.7,show.legend = T)+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=12,face = "bold"),
        axis.text.y = element_text(color="black", size=12,face = "bold"))+
  #annotate("text",label="ANOSIM_transplant: R = 0.19, p = 0.027",x=-0.72,y=0.58, size=5,fontface="bold")+
  #annotate("text",label="ANOSIM_tissue: R = 0.38, p = 0.001",x=-0.72,y=0.50, size=5,fontface="bold")+
  theme(legend.title=element_text(size=15,face = "bold"), legend.text=element_text(size=12,face = "bold"),
        axis.title=element_text(size=15,face = "bold"))+
  ggtitle("18S DNA")

#18s RNA
r4<-
asv_18s_apha %>% 
  rownames_to_column(., var = "sample_id") %>% 
  as_tibble() %>%  
  mutate(ROV = sample_id %>% str_sub(5,8)) %>%
  mutate(Depth_cmbs = sample_id %>%
           str_split("_") %>% sapply('[', 2) %>% 
           str_replace("^5$","05")) %>%
  mutate(Ribosome = sample_id %>%
           str_split("_") %>% sapply('[', 3)) %>%
  mutate(Nucleic = sample_id %>%
           str_split("_") %>% sapply('[', 4)) %>% 
  filter(Nucleic == "rna" &
           Depth_cmbs < 15) %>% 
  ggplot(aes(x = ROV, y = Richness)) +
  geom_boxplot(outlier.colour = "white", outlier.fill = "white", size =2, color = "grey66")+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
             size = 4,aes(color= Depth_cmbs), alpha = 0.7, show.legend = T)+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=12,face = "bold"),
        axis.text.y = element_text(color="black", size=12,face = "bold"))+
  #annotate("text",label="ANOSIM_transplant: R = 0.19, p = 0.027",x=-0.72,y=0.58, size=5,fontface="bold")+
  #annotate("text",label="ANOSIM_tissue: R = 0.38, p = 0.001",x=-0.72,y=0.50, size=5,fontface="bold")+
  theme(legend.title=element_text(size=15,face = "bold"), legend.text=element_text(size=12,face = "bold"),
        axis.title=element_text(size=15,face = "bold"))+
  ggtitle("18S RNA")

r5<-
asv_18s_apha %>% 
  rownames_to_column(., var = "sample_id") %>% 
  as_tibble() %>%  
  mutate(ROV = sample_id %>% str_sub(5,8)) %>%
  mutate(Depth_cmbs = sample_id %>%
           str_split("_") %>% sapply('[', 2) %>% 
           str_replace("^5$","05")) %>%
  mutate(Ribosome = sample_id %>%
           str_split("_") %>% sapply('[', 3)) %>%
  mutate(Nucleic = sample_id %>%
           str_split("_") %>% sapply('[', 4)) %>% 
  filter(Nucleic == "rna" &
           Depth_cmbs < 15) %>% 
  ggplot(aes(x = ROV, y  = Chao1)) +
  geom_boxplot(outlier.colour = "white", outlier.fill = "white", size =2, color = "grey66")+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
             size =4,aes(color= Depth_cmbs), alpha = 0.7, show.legend = T)+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=12,face = "bold"),
        axis.text.y = element_text(color="black", size=12,face = "bold"))+
  #annotate("text",label="ANOSIM_transplant: R = 0.19, p = 0.027",x=-0.72,y=0.58, size=5,fontface="bold")+
  #annotate("text",label="ANOSIM_tissue: R = 0.38, p = 0.001",x=-0.72,y=0.50, size=5,fontface="bold")+
  theme(legend.title=element_text(size=15,face = "bold"), legend.text=element_text(size=12,face = "bold"),
        axis.title=element_text(size=15,face = "bold"))+
  ggtitle("18S RNA")

r6<-
asv_18s_apha %>% 
  rownames_to_column(., var = "sample_id") %>% 
  as_tibble() %>%  
  mutate(ROV = sample_id %>% str_sub(5,8)) %>%
  mutate(Depth_cmbs = sample_id %>%
           str_split("_") %>% sapply('[', 2) %>% 
           str_replace("^5$","05")) %>%
  mutate(Ribosome = sample_id %>%
           str_split("_") %>% sapply('[', 3)) %>%
  mutate(Nucleic = sample_id %>%
           str_split("_") %>% sapply('[', 4)) %>% 
  filter(Nucleic == "rna" &
           Depth_cmbs < 15) %>% 
  ggplot(aes(x = ROV, y = PD_whole_tree)) +
  geom_boxplot(outlier.colour = "white", outlier.fill = "white", size =2, color = "grey66")+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
             size =4,aes(color= Depth_cmbs),  alpha = 0.7,show.legend = T)+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=12,face = "bold"),
        axis.text.y = element_text(color="black", size=12,face = "bold"))+
  #annotate("text",label="ANOSIM_transplant: R = 0.19, p = 0.027",x=-0.72,y=0.58, size=5,fontface="bold")+
  #annotate("text",label="ANOSIM_tissue: R = 0.38, p = 0.001",x=-0.72,y=0.50, size=5,fontface="bold")+
  theme(legend.title=element_text(size=15,face = "bold"), legend.text=element_text(size=12,face = "bold"),
        axis.title=element_text(size=15,face = "bold"))+
  ggtitle("18S RNA")


ggarrange(r1,r4,r2,r5,r3,r6,ncol = 2, nrow = 3)

##16s DNA and RNA

asv_table_16s = read.table("asv_table_16s.txt",header = T,check.names = F) %>% t()
tree_16s = read.tree("tree-16s-rooted.nwk")
asv_16s_apha = asv_table_16s %>% alpha(tree_16s)

#16S DNA
r7<-
asv_16s_apha %>% 
  rownames_to_column(., var = "sample_id") %>% 
  as_tibble() %>%  
  mutate(ROV = sample_id %>% str_sub(5,8)) %>%
  mutate(Depth_cmbs = sample_id %>%
           str_split("_") %>% sapply('[', 2) %>% 
           str_replace("^5$","05")) %>%
  mutate(Ribosome = sample_id %>%
           str_split("_") %>% sapply('[', 3)) %>%
  mutate(Nucleic = sample_id %>%
           str_split("_") %>% sapply('[', 4)) %>% 
  filter(Nucleic == "dna" &
           Depth_cmbs < 15) %>% 
  ggplot(aes(x = ROV, y = Richness)) +
  geom_boxplot(outlier.colour = "white", outlier.fill = "white", size =2, color = "grey66")+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
             size = 4,aes(color= Depth_cmbs),alpha = 0.7, show.legend = T)+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=12,face = "bold"),
        axis.text.y = element_text(color="black", size=12,face = "bold"))+
  #annotate("text",label="ANOSIM_transplant: R = 0.19, p = 0.027",x=-0.72,y=0.58, size=5,fontface="bold")+
  #annotate("text",label="ANOSIM_tissue: R = 0.38, p = 0.001",x=-0.72,y=0.50, size=5,fontface="bold")+
  theme(legend.title=element_text(size=15,face = "bold"), legend.text=element_text(size=12,face = "bold"),
        axis.title=element_text(size=15,face = "bold"))+
  ggtitle("16S DNA")

r8<-
asv_16s_apha %>% 
  rownames_to_column(., var = "sample_id") %>% 
  as_tibble() %>%  
  mutate(ROV = sample_id %>% str_sub(5,8)) %>%
  mutate(Depth_cmbs = sample_id %>%
           str_split("_") %>% sapply('[', 2) %>% 
           str_replace("^5$","05")) %>%
  mutate(Ribosome = sample_id %>%
           str_split("_") %>% sapply('[', 3)) %>%
  mutate(Nucleic = sample_id %>%
           str_split("_") %>% sapply('[', 4)) %>% 
  filter(Nucleic == "dna" &
           Depth_cmbs < 15) %>% 
  ggplot(aes(x = ROV, y  = Chao1)) +
  geom_boxplot(outlier.colour = "white", outlier.fill = "white", size =2, color = "grey66")+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
             size = 4,aes(color= Depth_cmbs), alpha = 0.7,show.legend = T)+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=12,face = "bold"),
        axis.text.y = element_text(color="black", size=12,face = "bold"))+
  #annotate("text",label="ANOSIM_transplant: R = 0.19, p = 0.027",x=-0.72,y=0.58, size=5,fontface="bold")+
  #annotate("text",label="ANOSIM_tissue: R = 0.38, p = 0.001",x=-0.72,y=0.50, size=5,fontface="bold")+
  theme(legend.title=element_text(size=15,face = "bold"), legend.text=element_text(size=12,face = "bold"),
        axis.title=element_text(size=15,face = "bold"))+
  ggtitle("16S DNA")

r9<-
asv_16s_apha %>% 
  rownames_to_column(., var = "sample_id") %>% 
  as_tibble() %>%  
  mutate(ROV = sample_id %>% str_sub(5,8)) %>%
  mutate(Depth_cmbs = sample_id %>%
           str_split("_") %>% sapply('[', 2) %>% 
           str_replace("^5$","05")) %>%
  mutate(Ribosome = sample_id %>%
           str_split("_") %>% sapply('[', 3)) %>%
  mutate(Nucleic = sample_id %>%
           str_split("_") %>% sapply('[', 4)) %>% 
  filter(Nucleic == "dna" &
           Depth_cmbs < 15) %>% 
  ggplot(aes(x = ROV, y = PD_whole_tree)) +
  geom_boxplot(outlier.colour = "white", outlier.fill = "white", size =2, color = "grey66")+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
             size = 4,aes(color= Depth_cmbs),alpha = 0.7, show.legend = T)+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=12,face = "bold"),
        axis.text.y = element_text(color="black", size=12,face = "bold"))+
  #annotate("text",label="ANOSIM_transplant: R = 0.19, p = 0.027",x=-0.72,y=0.58, size=5,fontface="bold")+
  #annotate("text",label="ANOSIM_tissue: R = 0.38, p = 0.001",x=-0.72,y=0.50, size=5,fontface="bold")+
  theme(legend.title=element_text(size=15,face = "bold"), legend.text=element_text(size=12,face = "bold"),
        axis.title=element_text(size=15,face = "bold"))+
  ggtitle("16S DNA")

#16S RNA
r10<-
asv_16s_apha %>% 
  rownames_to_column(., var = "sample_id") %>% 
  as_tibble() %>%  
  mutate(ROV = sample_id %>% str_sub(5,8)) %>%
  mutate(Depth_cmbs = sample_id %>%
           str_split("_") %>% sapply('[', 2) %>% 
           str_replace("^5$","05")) %>%
  mutate(Ribosome = sample_id %>%
           str_split("_") %>% sapply('[', 3)) %>%
  mutate(Nucleic = sample_id %>%
           str_split("_") %>% sapply('[', 4)) %>% 
  filter(Nucleic == "rna" &
           Depth_cmbs < 15) %>% 
  ggplot(aes(x = ROV, y = Richness)) +
  geom_boxplot(outlier.colour = "white", outlier.fill = "white", size =2, color = "grey66")+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
             size = 4,aes(color= Depth_cmbs), alpha = 0.7,show.legend = T)+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=12,face = "bold"),
        axis.text.y = element_text(color="black", size=12,face = "bold"))+
  #annotate("text",label="ANOSIM_transplant: R = 0.19, p = 0.027",x=-0.72,y=0.58, size=5,fontface="bold")+
  #annotate("text",label="ANOSIM_tissue: R = 0.38, p = 0.001",x=-0.72,y=0.50, size=5,fontface="bold")+
  theme(legend.title=element_text(size=15,face = "bold"), legend.text=element_text(size=12,face = "bold"),
        axis.title=element_text(size=15,face = "bold"))+
  ggtitle("16S RNA")

r11<-
asv_16s_apha %>% 
  rownames_to_column(., var = "sample_id") %>% 
  as_tibble() %>%  
  mutate(ROV = sample_id %>% str_sub(5,8)) %>%
  mutate(Depth_cmbs = sample_id %>%
           str_split("_") %>% sapply('[', 2) %>% 
           str_replace("^5$","05")) %>%
  mutate(Ribosome = sample_id %>%
           str_split("_") %>% sapply('[', 3)) %>%
  mutate(Nucleic = sample_id %>%
           str_split("_") %>% sapply('[', 4)) %>% 
  filter(Nucleic == "rna" &
           Depth_cmbs < 15) %>% 
  ggplot(aes(x = ROV, y  = Chao1)) +
  geom_boxplot(outlier.colour = "white", outlier.fill = "white", size =2, color = "grey66")+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
             size = 4,aes(color= Depth_cmbs), alpha = 0.7,show.legend = T)+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=12,face = "bold"),
        axis.text.y = element_text(color="black", size=12,face = "bold"))+
  #annotate("text",label="ANOSIM_transplant: R = 0.19, p = 0.027",x=-0.72,y=0.58, size=5,fontface="bold")+
  #annotate("text",label="ANOSIM_tissue: R = 0.38, p = 0.001",x=-0.72,y=0.50, size=5,fontface="bold")+
  theme(legend.title=element_text(size=15,face = "bold"), legend.text=element_text(size=12,face = "bold"),
        axis.title=element_text(size=15,face = "bold"))+
  ggtitle("16S RNA")

r12<-
asv_16s_apha %>% 
  rownames_to_column(., var = "sample_id") %>% 
  as_tibble() %>%  
  mutate(ROV = sample_id %>% str_sub(5,8)) %>%
  mutate(Depth_cmbs = sample_id %>%
           str_split("_") %>% sapply('[', 2) %>% 
           str_replace("^5$","05")) %>%
  mutate(Ribosome = sample_id %>%
           str_split("_") %>% sapply('[', 3)) %>%
  mutate(Nucleic = sample_id %>%
           str_split("_") %>% sapply('[', 4)) %>% 
  filter(Nucleic == "rna" &
           Depth_cmbs < 15) %>% 
  ggplot(aes(x = ROV, y = PD_whole_tree)) +
  geom_boxplot(outlier.colour = "white", outlier.fill = "white", size =2, color = "grey66")+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
             size = 4,aes(color= Depth_cmbs), alpha = 0.7,show.legend = T)+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=12,face = "bold"),
        axis.text.y = element_text(color="black", size=12,face = "bold"))+
  #annotate("text",label="ANOSIM_transplant: R = 0.19, p = 0.027",x=-0.72,y=0.58, size=5,fontface="bold")+
  #annotate("text",label="ANOSIM_tissue: R = 0.38, p = 0.001",x=-0.72,y=0.50, size=5,fontface="bold")+
  theme(legend.title=element_text(size=15,face = "bold"), legend.text=element_text(size=12,face = "bold"),
        axis.title=element_text(size=15,face = "bold"))+
  ggtitle("16S RNA")


ggarrange(r7,r10,r8,r11,r9,r12, ncol = 2, nrow = 3)


## wilcox test
dat <-
asv_16s_apha %>% 
  rownames_to_column(., var = "sample_id") %>% 
  as_tibble() %>%  
  mutate(ROV = sample_id %>% str_sub(5,8)) %>%
  mutate(Depth_cmbs = sample_id %>%
           str_split("_") %>% sapply('[', 2) %>% 
           str_replace("^5$","05")) %>%
  mutate(Ribosome = sample_id %>%
           str_split("_") %>% sapply('[', 3)) %>%
  mutate(Nucleic = sample_id %>%
           str_split("_") %>% sapply('[', 4)) %>% 
  filter(Nucleic == "rna" &
           Depth_cmbs < 15) 

 
 wilcox.test(dat$Richness[str_detect(dat$ROV,"ROV1|ROV2|ROV3")],
             dat$Richness[str_detect(dat$ROV,"ROV4|ROV5")])
 
 
 wilcox.test(dat$Chao1[str_detect(dat$ROV,"ROV1|ROV2|ROV3")],
             dat$Chao1[str_detect(dat$ROV,"ROV4|ROV5")]) 
 
 wilcox.test(dat$PD_whole_tree[str_detect(dat$ROV,"ROV1|ROV2|ROV3")],
             dat$PD_whole_tree[str_detect(dat$ROV,"ROV4|ROV5")]) 
 


# sub-habitat and alpha diversity -----------------------------------------


# # sub-habitat
# met = read.table("sample-metadata-with habitat.txt",header = T)
# 
# 
# #18s dna
# rownames(met) <- paste0(rownames(met), "_18s_dna")
# 
# alpha_18s_dna = asv_18s_apha[str_detect(rownames(asv_18s_apha), "dna"),]
# 
# identical(alpha_18s_dna %>% rownames(), met %>% row.names())# true
# 
# alpha_18s_dna$sub_habitat = met$Habitat %>% as.factor()
# 
# 
# alpha_18s_dna %>% rownames_to_column(., var = "sample_id") %>%
#   as_tibble() %>%
#   mutate(ROV = sample_id %>% str_sub(5,8)) %>%
#   mutate(Depth_cmbs = sample_id %>%
#            str_split("_") %>% sapply('[', 2) %>%
#            str_replace("^5$","05") ) %>%
#   mutate(Ribosome = sample_id %>%
#            str_split("_") %>% sapply('[', 3)) %>%
#   ggplot(aes(x= ROV, y=Richness))+
#   geom_point(aes(color = sub_habitat), size =3)


# RNA_alpha diversity -----------------------------------------------------------

alpha_18s_rna = alpha(asv18s_rna)

#depth compare
alpha_18s_rna %>% rownames_to_column(., var = "sample_id") %>% 
  as_tibble() %>%  
  mutate(ROV = sample_id %>% str_sub(5,8)) %>%
  mutate(Depth_cmbs = sample_id %>%
           str_split("_") %>% sapply('[', 2) %>% 
           str_replace("^5$","05") ) %>%
  mutate(Ribosome = sample_id %>%
           str_split("_") %>% sapply('[', 3)) %>% 
  ggplot(aes(x= Depth_cmbs, y=Richness, group = Depth_cmbs))+
  geom_point(aes(color = ROV),size=2)+
  geom_boxplot()

#ROV compare

alpha_18s_rna %>% rownames_to_column(., var = "sample_id") %>% 
  as_tibble() %>%  
  mutate(ROV = sample_id %>% str_sub(5,8)) %>%
  mutate(Depth_cmbs = sample_id %>%
           str_split("_") %>% sapply('[', 2) %>% 
           str_replace("^5$","05") ) %>%
  mutate(Ribosome = sample_id %>%
           str_split("_") %>% sapply('[', 3)) %>% 
  ggplot(aes(x= ROV, y=Richness))+
  geom_point(aes(color = Depth_cmbs))+
  geom_boxplot()





# RNA:DNA  relative activity ---18s------------------------------------------

## Use Protists-only ASV; remove metazoan and land plants 

asv_table_18s = read.table("asv_18s.txt",header = T,check.names = F) %>% t()
tax_18s = read.delim("taxonomy-pr2.tsv")


## remove metazoan
tax_pro = tax_18s[!c(str_detect(tax_18s$Taxon, "Metazoa")),]
tax_pro = tax_pro[!c(str_detect(tax_pro$Taxon, "Streptophyta")),]
asv_table_18s = asv_table_18s[,match(tax_pro$Feature.ID, colnames(asv_table_18s))]
asv_table_18s = rrarefy(asv_table_18s,
                        min(rowSums(asv_table_18s)))

# rarefy
asv_table_18s = rrarefy(asv_table_18s,
                            min(rowSums(asv_table_18s)))

# relative abundance
asv_table_18s_ra = asv_table_18s / min(rowSums(asv_table_18s))

## divide to DNA and RNA samples

asv18s_ra_dna = 
  asv_table_18s_ra[asv_table_18s_ra %>% rownames() %>% str_detect("dna"), ]

asv18s_ra_rna = 
  asv_table_18s_ra[asv_table_18s_ra %>% rownames() %>% str_detect("rna"), ]

# Make DNA samples are equal to RNA samples, as there are much more DNA samples
asv18s_ra_rna %>% rownames()
m1 = asv18s_ra_rna %>% rownames() %>% str_remove("_rna")
asv18s_ra_dna %>% rownames()
m2 = asv18s_ra_dna %>% rownames() %>% str_remove("_dna")

asv18s_ra_dna = asv18s_ra_dna[match(m1,m2),]
identical(rownames(asv18s_ra_dna) %>% str_remove("_dna"),
          rownames(asv18s_ra_rna) %>% str_remove("_rna")) # True

identical(colnames(asv18s_ra_dna) %>% str_remove("_dna"),
          colnames(asv18s_ra_rna) %>% str_remove("_rna")) # True

# relative activity: %RNA : %DNA
rel_act_18s = t(asv18s_ra_rna) / t(asv18s_ra_dna)

sum_ra <- function(x)  x[!is.nan(x) & !is.infinite(x)] %>% mean


options(scipen = 200)## close Scientific notation


sp = c("Fungi","Labyrinthulomycetes",
       "Amoebozoa","Apicomplexa", "Breviatea", "Syndiniales",
       "Cercozoa", "Ciliophora","Radiolaria",
       "Dinophyceae", "Bacillariophyceae","Chlorophyta")


ra_18s = list()

for (i in seq_along(sp)){
  
  p <- 
  rel_act_18s %>% sqrt() %>%  data.frame %>% 
  rownames_to_column(., var = "ASV_ID") %>%
  as_tibble() %>% 
  mutate(Taxonomy = tax_18s$Taxon[match(colnames(asv18s_ra_dna),
                                        tax_18s$Feature.ID)]) %>% 
  filter(Taxonomy %>% str_detect(sp[i])) %>%  ######## variable here
  select(-Taxonomy) %>% 
  pivot_longer(!ASV_ID, names_to = "Sample", values_to = "RD_ratio") %>% 
  mutate(ROV = Sample %>% str_sub(5,8))  %>%
  mutate(Depth_cmbs = Sample %>% str_split("_") %>% sapply('[', 2)) %>% 
  filter(RD_ratio != "NaN") %>%
  filter(RD_ratio %>% is.finite()) %>% 
  ggplot(aes(x = ROV, y = RD_ratio, group = ROV))+
  # geom_boxplot(outlier.colour = "white", outlier.fill = "white", size =1, color = "grey66")+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
             size = 3,aes(color= ROV), show.legend = T)+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12))+
  theme(legend.title=element_text(size=12), legend.text=element_text(size=12),
        axis.title=element_text(size=12))+
  ylab("RNA:DNA ratio (Sqrt)")+xlab("Habitat")+
  #ylim(0,15)+
  scale_color_manual(values = c("firebrick1", "hotpink", "darkgoldenrod1", "cornflowerblue","darkcyan"))+
  ggtitle(sp[i])+                            ######## variable here
    theme(legend.position="none")+
    geom_boxplot(outlier.shape = NA)+
    geom_hline(yintercept=1, linetype="dashed", 
               color = "grey66", linewidth=1)+
    scale_x_discrete(limits = c("ROV1", "ROV2","ROV3","ROV4","ROV5"))
  
  ra_18s[[i]] <- p
}

ggarrange(ra_18s[[1]],ra_18s[[2]],ra_18s[[3]],ra_18s[[4]],ra_18s[[5]],ra_18s[[6]],ra_18s[[7]],
          ra_18s[[8]],ra_18s[[9]],ra_18s[[10]],ra_18s[[11]],ra_18s[[12]],
          nrow = 3,ncol = 4)

## calculate number of zero value
rel_act_18s %>% data.frame %>% rownames_to_column(., var = "ASV") %>% as_tibble() %>% 
  pivot_longer(!ASV, names_to = "Sample", values_to = "RA") %>% 
  filter(RA != "Inf") %>% 
  filter(RA != "NaN") %>%    ### remain 14518 ASVs
  filter(RA == 0) ### remain 13395 ASVs
### ONE WAY ANOVA TEST

df <-
ra_18s[[1]]$data %>% 
  mutate(Habitat = ROV %>% 
           str_replace_all(c("ROV1" = "Active Seep", 
                             "ROV2" = "Active Seep",
                             "ROV3" = "Seep",
                             "ROV4" = "Non-seep",
                             "ROV5" = "Non-seep")) %>% as.factor())
aov(RD_ratio ~ Habitat, df) %>% summary() 


### Wilcoxon test

ra_18s[[3]]$data %>% filter(ROV == "ROV1" | ROV == "ROV2" | ROV == "ROV3") %>% pull(RD_ratio) -> v1
ra_18s[[3]]$data %>% filter(ROV == "ROV4" | ROV == "ROV5") %>% pull(RD_ratio) -> v2
t.test(v1, v2)
wilcox.test(v1,v2)          

mean(v1)
mean(v2)
v2


# RNA:DNA  relative activity - 16s ----------------------------------------

asv_table_16s = read.table("asv_table_16s.txt",header = T,check.names = F) %>% t()
asv_table_16s %>% dim
tax_16s = read.delim("taxonomy-16s.tsv")


asv_table_16s = rrarefy(asv_table_16s,
                        min(rowSums(asv_table_16s))) #46886


# relative abundance
asv_table_16s_ra = asv_table_16s / min(rowSums(asv_table_16s))


## divide to DNA and RNA samples

asv16s_ra_dna = 
  asv_table_16s_ra[asv_table_16s_ra %>% rownames() %>% str_detect("dna"), ]

asv16s_ra_rna = 
  asv_table_16s_ra[asv_table_16s_ra %>% rownames() %>% str_detect("rna"), ]

# Make DNA samples are equal to RNA samples, as there are much more DNA samples
asv16s_ra_rna %>% rownames()
m1 = asv16s_ra_rna %>% rownames() %>% str_remove("_rna")
asv16s_ra_dna %>% rownames()
m2 = asv16s_ra_dna %>% rownames() %>% str_remove("_dna")

asv16s_ra_dna = asv16s_ra_dna[match(m1,m2),]
identical(rownames(asv16s_ra_dna) %>% str_remove("_dna"),
          rownames(asv16s_ra_rna) %>% str_remove("_rna")) # True

identical(colnames(asv16s_ra_dna) %>% str_remove("_dna"),
          colnames(asv16s_ra_rna) %>% str_remove("_rna")) # True

# relative activity: %RNA : %DNA
rel_act_16s = t(asv16s_ra_rna) / t(asv16s_ra_dna)

options(scipen = 200)## close Scientific notation


spp = c("Gammaproteobacteria", "ANME-1", "Desulfobacteria", "Campylobacteria", "JS1",
        "Methanosarcinia", "Desulfobulbia", "Bacilli", "Alphaproteobacteria")

ra_16s = list()

for (i in seq_along(spp)){
  
  p <- 
    rel_act_16s %>% sqrt  %>%  data.frame %>% 
    rownames_to_column(., var = "ASV_ID") %>%
    as_tibble() %>% 
    mutate(Taxonomy = tax_16s$Taxon[match(colnames(asv16s_ra_dna),
                                          tax_16s$Feature.ID)]) %>% 
    filter(Taxonomy %>% str_detect(spp[i])) %>%  ######## variable here
    select(-Taxonomy) %>% 
    pivot_longer(!ASV_ID, names_to = "Sample", values_to = "RD_ratio") %>% 
    mutate(ROV = Sample %>% str_sub(5,8))  %>%
    mutate(Depth_cmbs = Sample %>% str_split("_") %>% sapply('[', 2)) %>% 
    filter(RD_ratio != "NaN") %>%
    filter(RD_ratio %>% is.finite())%>%
    filter(Depth_cmbs %in% c("0","5","10")) %>% 
    ggplot(aes(x = ROV, y = RD_ratio, group = ROV))+
    # geom_boxplot(outlier.colour = "white", outlier.fill = "white", size =1, color = "grey66")+
    geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
               size = 3,aes(color= ROV), show.legend = T)+
    theme_classic()+
    theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
    theme(axis.text.x = element_text(color="black", size=15),
          axis.text.y = element_text(color="black", size=15))+
    theme(legend.title=element_text(size=15), legend.text=element_text(size=15),
          axis.title=element_text(size=15))+
    ylab("RNA:DNA ratio (Sqrt)")+xlab("Habitat")+
    scale_color_manual(values = c("firebrick1", "hotpink", "darkgoldenrod1", "cornflowerblue","darkcyan"))+
    theme(legend.position="none")+
    geom_boxplot(outlier.shape = NA)+
    geom_hline(yintercept=1, linetype="dashed", 
               color = "grey66", linewidth=1)+
    scale_x_discrete(limits = c("ROV1", "ROV2","ROV3","ROV4","ROV5"))+
    ggtitle(spp[i])
  
  
  ra_16s[[i]] <- p
}

ggarrange(ra_16s[[1]],ra_16s[[2]],ra_16s[[3]],ra_16s[[4]],ra_16s[[5]],ra_16s[[6]],ra_16s[[7]],
          ra_16s[[8]],ra_16s[[9]],
          nrow = 3,ncol = 3)

#  TEST
ra_16s[[6]]$data %>% filter(ROV == "ROV1" | ROV == "ROV2" | ROV == "ROV3") %>% pull(RD_ratio) -> v1
ra_16s[[6]]$data %>% filter(ROV == "ROV4" | ROV == "ROV5") %>% pull(RD_ratio) -> v2
t.test(v1, v2)
wilcox.test(v1,v2)          

mean(v1)
mean(v2)
v2


df <-
  ra_16s[[2]]$data %>% 
  mutate(Habitat = ROV %>% 
           str_replace_all(c("ROV1" = "Active Seep", 
                             "ROV2" = "Active Seep",
                             "ROV3" = "Seep",
                             "ROV4" = "Non-seep",
                             "ROV5" = "Non-seep")) %>% as.factor())
aov(RD_ratio ~ Habitat, df) %>% summary() 


## calculate number of zero value
rel_act_16s %>% data.frame %>% rownames_to_column(., var = "ASV") %>% as_tibble() %>% 
  pivot_longer(!ASV, names_to = "Sample", values_to = "RA") %>% 
  filter(RA != "Inf") %>% 
  filter(RA != "NaN") %>%    ### remain 14518 ASVs
  filter(RA == 0) ### remain 13395 ASVs

# Rna Dna ratio 1:1 PLOT -------------------------------------------------------------
#18s

asv18s_ra_rna %>% t() %>% data.frame() %>% 
  rownames_to_column(., var = "asv_id") %>% 
  as_tibble() %>% 
  pivot_longer(!asv_id, names_to = "rna_sample", values_to = "rna_ratio") %>% 
  mutate(sample_uni = rna_sample %>% str_remove("_rna")) -> df_rna


asv18s_ra_dna %>% t() %>% data.frame() %>% 
  rownames_to_column(., var = "asv_id") %>% 
  as_tibble() %>% 
  pivot_longer(!asv_id, names_to = "dna_sample", values_to = "dna_ratio") %>% 
  mutate(sample_uni = dna_sample %>% str_remove("_dna")) -> df_dna

identical(df_rna$asv_id, df_dna$asv_id)
identical(df_rna$sample_uni, df_dna$sample_uni)

df =
  cbind(df_rna[,c(1,4,3)], df_dna[,3]) %>% 
  mutate(ROV = sample_uni %>% str_sub(5,8))  %>%
  mutate(Depth_cmbs = sample_uni %>% str_split("_") %>% sapply('[', 2)) 

df$Taxonomy = tax_18s$Taxon[match(df$asv_id, tax_18s$Feature.ID)]

## for silva database
df %>% as_tibble %>% 
  filter(!(dna_ratio == "0" & dna_ratio == "0")) %>% 
  mutate(Protists = Taxonomy %>% 
           str_replace_all(c(".*Fungi.*" = "Fungi",
                             ".*Labyrinthulomycetes.*" = "Labyrinthulomycetes",
                             ".*Amoebozoa.*" = "Amoebozoa",
                             ".*Apicomplexa.*" = "Apicomplexa",
                             ".*Breviatea.*" = "Breviatea",
                             ".*Protalveolata.*" = "Protalveolata",
                             ".*Cercozoa.*" = "Cercozoa",
                             ".*Ciliophora.*" = "Ciliophora",
                             ".*Dinoflagellata.*" = "Dinoflagellata",
                             ".*Ochrophyta.*" = "Ochrophyta",
                             ".*Charophyta.*" = "Charophyta",
                             ".*Chlorophyta.*" = "Chlorophyta"))) %>% 
  mutate(Protists = Protists %>% str_replace(".*D_0__Eukaryota.*", "Others")) %>% 
  filter(!Protists == "Others") %>% 
  mutate(Habitat = ROV %>%  str_replace_all(c("ROV1" = "Seep",
                                              "ROV2" = "Seep",
                                              "ROV3" = "Seep",
                                              "ROV4" = "Non-seep",
                                              "ROV5" = "Non-seep"))) %>% 
  ggplot(aes(x = dna_ratio, y = rna_ratio))+
  geom_point(aes(color = Protists, shape = Habitat),size =3)+
  scale_y_log10(limits=c(0.0005,1))+
  scale_x_log10(limits=c(0.0005,1))+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0))+
  theme(axis.text.x = element_text(color="black", size=15),
        axis.text.y = element_text(color="black", size=15))+
  theme(legend.title=element_text(size=15), legend.text=element_text(size=15),
        axis.title=element_text(size=15))+
  ylab("%RNA (18S)")+xlab("%DNA (18S)")
  # scale_color_manual(values = c("firebrick1", "hotpink", "darkgoldenrod1", "cornflowerblue","darkcyan"))


# RA vs 13C ---------------------------------------------------------------

RA_18s = asv18s_rna / asv18s_dna

sum_ra <- function(x)  x[!is.nan(x) & !is.infinite(x)] %>% mean

data.frame(RA = apply(RA_18s, 1, sum_ra)) %>% 
  rownames_to_column(., var= "Sample") %>% as_tibble() %>% 
  mutate(ROV = Sample %>% str_sub(5,8)) %>% 
  ggplot(aes(x = ROV, y=RA, group = ROV))+
  geom_boxplot()


RA_16s = asv16s_rna / asv16s_dna

data.frame(RA = apply(RA_16s, 1, sum_ra)) %>% 
  rownames_to_column(., var= "Sample") %>% as_tibble() %>% 
  mutate(ROV = Sample %>% str_sub(5,8)) %>% 
  ggplot(aes(x = ROV, y=RA, group = ROV))+
  geom_boxplot()

### for PR2 database
df %>% as_tibble %>% 
  filter(!(dna_ratio == "0" & dna_ratio == "0")) %>% 
  mutate(Protists = Taxonomy %>% 
           str_replace_all(c(".*Fungi.*" = "Fungi",
                             ".*Labyrinthulomycetes.*" = "Labyrinthulomycetes",
                             ".*Amoebozoa.*" = "Amoebozoa",
                             ".*Apicomplexa.*" = "Apicomplexa",
                             ".*Breviatea.*" = "Breviatea",
                             ".*Syndiniales.*" = "Dinoflagellata (Syndiniales)",
                             ".*Cercozoa.*" = "Cercozoa",
                             ".*Ciliophora.*" = "Ciliophora",
                             ".*Radiolaria.*" = "Radiolaria",
                             ".*Chlorophyta.*" = "Chlorophyta",
                             ".*Haptophy.*" = "Haptophyta",
                             ".*Bacillariophyceae.*" = "Diatoms",
                             ".*Mediophyceae.*" = "Diatoms"
                         ))) %>% 
  filter(Protists %>% str_detect(
    
    "Fungi|Labyrinthulomycetes|Amoebozoa|Apicomplexa|Breviatea|Syndiniales|
    Cercozoa|Ciliophora|Radiolaria|Chlorophyta|Haptophyta|Diatoms")) %>% 
  
  # pull(Protists) %>% as.factor() %>% levels()
  
  mutate(Habitat = ROV %>%  str_replace_all(c("ROV1" = "Seep",
                                              "ROV2" = "Seep",
                                              "ROV3" = "Seep",
                                              "ROV4" = "Non-seep",
                                              "ROV5" = "Non-seep"))) %>% 
  ggplot(aes(x = dna_ratio, y = rna_ratio))+
  geom_point(aes(color = Protists, shape = Habitat),size = 3)+
  scale_y_log10(limits=c(0.0005,1))+
  scale_x_log10(limits=c(0.0005,1))+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=1))+
  theme(axis.text.x = element_text(color="black", size=15),
        axis.text.y = element_text(color="black", size=15))+
  theme(legend.title=element_text(size=15), legend.text=element_text(size=15),
        axis.title=element_text(size=15))+
  ylab("%RNA (18S)")+xlab("%DNA (18S)")
# scale_color_manual(values = c("firebrick1", "hotpink", "darkgoldenrod1", "cornflowerblue","darkcyan"))


#16s
asv16s_ra_rna %>% t() %>% data.frame() %>% 
  rownames_to_column(., var = "asv_id") %>% 
  as_tibble() %>% 
  pivot_longer(!asv_id, names_to = "rna_sample", values_to = "rna_ratio") %>% 
  mutate(sample_uni = rna_sample %>% str_remove("_rna")) -> df_rna


asv16s_ra_dna %>% t() %>% data.frame() %>% 
  rownames_to_column(., var = "asv_id") %>% 
  as_tibble() %>% 
  pivot_longer(!asv_id, names_to = "dna_sample", values_to = "dna_ratio") %>% 
  mutate(sample_uni = dna_sample %>% str_remove("_dna")) -> df_dna

identical(df_rna$asv_id, df_dna$asv_id)
identical(df_rna$sample_uni, df_dna$sample_uni)

df =
  cbind(df_rna[,c(1,4,3)], df_dna[,3]) %>% 
  mutate(ROV = sample_uni %>% str_sub(5,8))  %>%
  mutate(Depth_cmbs = sample_uni %>% str_split("_") %>% sapply('[', 2)) 

df$Taxonomy = tax_16s$Taxon[match(df$asv_id, tax_16s$Feature.ID)]


df %>% as_tibble %>% 
  filter(!(dna_ratio == "0" & dna_ratio == "0")) %>% 
  mutate(Tax = Taxonomy %>% str_replace_all(c(".*Gammaproteobacteria.*" = "Gammaproteobacteria",
                                              ".*ANME-1.*" = "ANME-1",
                                              ".*Desulfobacteria.*" = "Desulfobacteria",
                                              ".*Campylobacteria.*" = "Campylobacteria",
                                              ".*JS1.*" = "JS1",
                                              ".*Methanosarcinia.*" = "Methanosarcinia",
                                              ".*Desulfobulbia.*" = "Desulfobulbia",
                                              ".*Bacilli.*" = "Bacilli",
                                              ".*Alphaproteobacteria.*" = "Alphaproteobacteria",
                                              ".*Anaerolineae.*" = "Anaerolineae"
                                              ))) %>% 
  filter(!Tax %>% str_detect(".*d__Euk.*")) %>% 
  filter(!Tax %>% str_detect(".*Unassigned.*")) %>% 
  filter(!Tax %>% str_detect(".*d__Bacteria.*")) %>% 
  filter(!Tax %>% str_detect(".*d__Archaea.*")) %>% 
  mutate(Habitat = ROV %>%  str_replace_all(c("ROV1" = "Seep",
                                              "ROV2" = "Seep",
                                              "ROV3" = "Seep",
                                              "ROV4" = "Non-seep",
                                              "ROV5" = "Non-seep"))) %>% 
  ggplot(aes(x = dna_ratio, y = rna_ratio))+
  geom_point(aes(color = Tax, shape = Habitat),size =2)+
  scale_y_log10(limits=c(0.00001,1))+
  scale_x_log10(limits=c(0.00001,1))+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=15),
        axis.text.y = element_text(color="black", size=15))+
  theme(legend.title=element_text(size=15), legend.text=element_text(size=15),
        axis.title=element_text(size=15))+
  ylab("%RNA (16S)")+xlab("%DNA (16S)")




# Community composition ---------------------------------------------------
#18s
com = read.csv("level-4-18s.csv",check.names = F)
com = com[,1:194] %>% column_to_rownames(., var = "index")
com = rrarefy(com, min(rowSums(com)))
com = com/min(rowSums(com))
# extract the most abundant 15 taxa, with sum of remaining as "Others"
com1 = com %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:10] %>% data.frame() %>% rownames()
m = match(com1, colnames(com))
rownames(com) = gsub("_5_", "_05_",rownames(com))
# barplot of community composition-relative abundance
df = 
  cbind(com[,m], data.frame(apply(com[,-m],1,sum)))%>% 
  rename(Others = `apply.com....m...1..sum.`) %>% t() %>% 
  data.frame() %>% 
  rownames_to_column(., var = "Microeukaryotes") %>% as_tibble() %>% 
  pivot_longer(!Microeukaryotes, names_to = "Sample", values_to = "Seqs") %>% 
  mutate(ROV = Sample %>% str_sub(5,8)) %>% 
  mutate(Sample = Sample %>% str_replace("_5_", "_05_")) %>% 
  mutate(Nucleic = Sample %>% str_split("_") %>% sapply('[', 4)) %>% 
  mutate(Depth_cmbs = Sample %>% str_split("_") %>% sapply('[', 2)) %>% 
  filter(Nucleic == "rna")

df$Sample = factor(df$Sample, levels = rownames(com)[str_detect(rownames(com), "rna")])

g =  
  ggplot(df, aes(x = Sample, y = Seqs, fill = Microeukaryotes))+
  geom_bar(stat="identity", width = 1)+
  theme_classic()+
  theme(axis.text.y = element_text(color="black", size=15))+
  theme(legend.title=element_text(size=15), legend.text=element_text(size=5),
        axis.title=element_text(size=15))+
  ylab("Relative abundance")+
  ylim(0,1)+
  theme(legend.position="bottom")+
  ggtitle("18S rna community")+
  #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
  theme(axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

g  


## extract x-labels and make depth file
ggplot_build(g)$layout$panel_params[[1]]$x$breaks %>% 
  as_tibble() %>% mutate(Depth_cmbs = value %>% str_split("_") %>% 
                           sapply('[', 2) %>% as.numeric()) %>% 
  mutate(ROV = value %>% str_sub(5,8) %>% as.factor()) ->df1

df1$value = factor(df1$value, levels = rownames(com)[str_detect(rownames(com), "rna")])

ggplot(df1, aes(x = value, y = Depth_cmbs))+
  geom_bar(stat="identity", width = 1,aes(fill=ROV))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
  scale_fill_manual(values = c("blue", "dodgerblue", "cyan", "gray33","gray66"))


#16s
com = read.csv("level-3-16s.csv",check.names = F)
com = com[,1:261] %>% column_to_rownames(., var = "index")
com = rrarefy(com, min(rowSums(com)))
com = com/46886
# extract the most abundant 15 taxa, with sum of remaining as "Others"
com1 = com %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:10] %>% data.frame() %>% rownames()
m = match(com1, colnames(com))
rownames(com) = gsub("_5_", "_05_",rownames(com))
# barplot of community composition-relative abundance
df = 
cbind(com[,m], data.frame(apply(com[,-m],1,sum)))%>% 
  rename(Others = `apply.com....m...1..sum.`) %>% t() %>% 
  data.frame() %>% 
  rownames_to_column(., var = "Taxa_level_3") %>% as_tibble() %>% 
  pivot_longer(!Taxa_level_3, names_to = "Sample", values_to = "Seqs") %>% 
  mutate(ROV = Sample %>% str_sub(5,8)) %>% 
  mutate(Sample = Sample %>% str_replace("_5_", "_05_")) %>% 
  mutate(Nucleic = Sample %>% str_split("_") %>% sapply('[', 4)) %>% 
  mutate(Depth_cmbs = Sample %>% str_split("_") %>% sapply('[', 2)) %>% 
  filter(Nucleic == "rna")

df$Sample = factor(df$Sample, levels = rownames(com)[str_detect(rownames(com), "rna")])
 
g =  
  ggplot(df, aes(x = Sample, y = Seqs, fill = Taxa_level_3))+
  geom_bar(stat="identity", width = 1)+
  theme_classic()+
  theme(axis.text.y = element_text(color="black", size=15))+
  theme(legend.title=element_text(size=15), legend.text=element_text(size=5),
        axis.title=element_text(size=15))+
  ylab("Relative abundance")+
  ylim(0,1)+
  theme(legend.position="bottom")+
  ggtitle("16S rna community")+
 #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
  theme(axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
  
g  


## extract x-labels and make depth file
ggplot_build(g)$layout$panel_params[[1]]$x$breaks %>% 
  as_tibble() %>% mutate(Depth_cmbs = value %>% str_split("_") %>% 
                           sapply('[', 2) %>% as.numeric()) %>% 
  mutate(ROV = value %>% str_sub(5,8) %>% as.factor()) ->df1

df1$value = factor(df1$value, levels = rownames(com)[str_detect(rownames(com), "rna")])

  ggplot(df1, aes(x = value, y = Depth_cmbs))+
  geom_bar(stat="identity", width = 1,aes(fill=ROV))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+
    scale_fill_manual(values = c("blue", "dodgerblue", "cyan", "gray33","gray66"))
  
  
  

# alpha diversity and C N -------------------------------------------------
setwd("C:/Users/XZM/Desktop/Coldseep 2nd study/2021 C13 N15 measured by HKU")
df = read.csv("2021coldseepCN.csv",check.names = F)
df$Sample_name = df$Sample_name %>%  str_replace("-","_")

df %>% as_tibble() %>% 
  mutate(ROV = str_extract(Sample_name, "ROV\\d")) %>%
  mutate(Depth_m = str_split(Sample_name,"_") %>% sapply('[', 2)) %>%
  mutate(`C/N` = `%C` / `%N`) %>% 
  # filter(Depth_m %in% c(0,5,10)) %>%
  mutate(DNA_18s_richness = 
           asv18s_dna_alpha$Richness[match(df$Sample_name, 
                                           rownames(asv18s_dna_alpha) %>% 
                                             str_split("_18s_") %>% sapply('[', 1))]) %>% 
  ggplot(aes(x = `C/N`, y = d13C,color = ROV))+
  geom_point(size = 5)+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth = 2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=30),
        axis.text.y = element_text(color="black", size=30))+
  theme(legend.title=element_text(size=30), legend.text=element_text(size=30),
        axis.title=element_text(size=40,face = "bold"))+
  xlab(expression(paste(C/N)))+
  ylab(expression(paste(delta^13*C)))+
  ggtitle("18S-DNA")

df %>% as_tibble() %>% 
  mutate(ROV = str_extract(Sample_name, "ROV\\d")) %>%
  mutate(Depth_m = str_split(Sample_name,"_") %>% sapply('[', 2)) %>%
  mutate(`C/N` = `%C` / `%N`) %>% 
  # filter(Depth_m %in% c(0,5,10)) %>%
  mutate(RNA_18s_richness = 
           asv18s_rna_alpha$Richness[match(df$Sample_name, 
                                           rownames(asv18s_rna_alpha) %>% 
                                             str_split("_18s_") %>% sapply('[', 1))]) %>% 
  ggplot(aes(x = `C/N`, y = d13C,color = ROV))+
  geom_point(aes(size = RNA_18s_richness))+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth = 2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=15),
        axis.text.y = element_text(color="black", size=15))+
  theme(legend.title=element_text(size=20), legend.text=element_text(size=15),
        axis.title=element_text(size=20))+
  xlab(expression(paste(C/N)))+
  ylab(expression(paste(delta^13*C)))+
  ggtitle("18S-RNA")
 

df %>% as_tibble() %>% 
  mutate(ROV = str_extract(Sample_name, "ROV\\d")) %>%
  mutate(Depth_m = str_split(Sample_name,"_") %>% sapply('[', 2)) %>%
  mutate(`C/N` = `%C` / `%N`) %>% 
  # filter(Depth_m %in% c(0,5,10)) %>%
  mutate(DNA_16s_richness = 
           asv16s_dna_alpha$Richness[match(df$Sample_name, 
                                           rownames(asv16s_dna_alpha) %>% 
                                             str_split("_16s_") %>% sapply('[', 1))]) %>% 
  ggplot(aes(x = `C/N`, y = d13C,color = ROV))+
  geom_point(aes(size = DNA_16s_richness))+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth = 2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=15),
        axis.text.y = element_text(color="black", size=15))+
  theme(legend.title=element_text(size=20), legend.text=element_text(size=15),
        axis.title=element_text(size=20))+
  xlab(expression(paste(C/N)))+
  ylab(expression(paste(delta^13*C)))+
  ggtitle("16S-DNA")

# asv16s_rna_alpha= alpha(asv16s_rna, tree_16s)

df %>% as_tibble() %>% 
  mutate(ROV = str_extract(Sample_name, "ROV\\d")) %>%
  mutate(Depth_m = str_split(Sample_name,"_") %>% sapply('[', 2)) %>%
  mutate(`C/N` = `%C` / `%N`) %>% 
  # filter(Depth_m %in% c(0,5,10)) %>%
  mutate(DNA_16s_richness = 
           asv16s_rna_alpha$Richness[match(df$Sample_name, 
                                           rownames(asv16s_dna_alpha) %>% 
                                             str_split("_16s_") %>% sapply('[', 1))]) %>% 
  ggplot(aes(x = `C/N`, y = d13C,color = ROV))+
  geom_point(aes(size = DNA_16s_richness))+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth = 2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=15),
        axis.text.y = element_text(color="black", size=15))+
  theme(legend.title=element_text(size=20), legend.text=element_text(size=15),
        axis.title=element_text(size=20))+
  xlab(expression(paste(C/N)))+
  ylab(expression(paste(delta^13*C)))+
  ggtitle("16S-RNA")

# beta diversity ---------------------------------------------------------



# Ecological process ------------------------------------------------------
library(scales)

euk_nti = read.csv("euk-bnti.csv")
euk_nti %>% head()

pro_nti = read.csv("pro-bnti.csv")
pro_nti %>% head  

df1 = data.frame(bNTI = c(euk_nti$euk.dna.r1 %>% na.omit(),
                          euk_nti$euk.dna.r2 %>% na.omit(),
                          euk_nti$euk.dna.r3 %>% na.omit(),
                          euk_nti$euk.dna.r4 %>% na.omit(),
                          euk_nti$euk.dna.r5 %>% na.omit(),
                          euk_nti$euk.rna.r1 %>% na.omit(),
                          euk_nti$euk.rna.r2 %>% na.omit(),
                          euk_nti$euk.rna.r3 %>% na.omit(),
                          euk_nti$euk.rna.r4 %>% na.omit(),
                          euk_nti$euk.rna.r5 %>% na.omit()),
                 Name = c(rep("dna_ROV1", euk_nti$euk.dna.r1 %>% na.omit() %>% length()),
                          rep("dna_ROV2", euk_nti$euk.dna.r2 %>% na.omit() %>% length()),
                          rep("dna_ROV3", euk_nti$euk.dna.r3 %>% na.omit() %>% length()),
                          rep("dna_ROV4", euk_nti$euk.dna.r4 %>% na.omit() %>% length()),
                          rep("dna_ROV5", euk_nti$euk.dna.r5 %>% na.omit() %>% length()),
                          rep("rna_ROV1", euk_nti$euk.rna.r1 %>% na.omit() %>% length()),
                          rep("rna_ROV2", euk_nti$euk.rna.r2 %>% na.omit() %>% length()),
                          rep("rna_ROV3", euk_nti$euk.rna.r3 %>% na.omit() %>% length()),
                          rep("rna_ROV4", euk_nti$euk.rna.r4 %>% na.omit() %>% length()),
                          rep("rna_ROV5", euk_nti$euk.rna.r5 %>% na.omit() %>% length())))
p1<-
df1 %>% as_tibble() %>% mutate(ROV = Name %>% str_split("_") %>% sapply('[', 2)) %>% 
  ggplot(aes(x = Name, y = bNTI))+
  geom_boxplot(outlier.colour = "white", outlier.fill = "white", size =1.5, aes(color = ROV))+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
             size = 2,aes(color = ROV), alpha = 0.5,show.legend = T)+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=15),
        axis.text.y = element_text(color="black", size=15))+
  theme(legend.title=element_text(size=15), legend.text=element_text(size=12),
        axis.title=element_text(size=15))+
  geom_hline(aes(yintercept = 2), color = "grey66", linetype = "dashed", linewidth = 1)+
  geom_hline(aes(yintercept = -2), color = "grey66", linetype = "dashed", linewidth = 1)+
  scale_color_manual(values = c("firebrick1", "hotpink", "darkgoldenrod1", "cornflowerblue","darkcyan"))



pro_nti %>% head  

df2 = data.frame(bNTI = c(pro_nti$asv16s_dna_rov1_bnti %>% na.omit(),
                          pro_nti$asv16s_dna_rov2_bnti %>% na.omit(),
                          pro_nti$asv16s_dna_rov3_bnti %>% na.omit(),
                          pro_nti$asv16s_dna_rov4_bnti %>% na.omit(),
                          pro_nti$asv16s_dna_rov5_bnti %>% na.omit(),
                          pro_nti$asv16s_rna_rov1_bnti %>% na.omit(),
                          pro_nti$asv16s_rna_rov2_bnti %>% na.omit(),
                          pro_nti$asv16s_rna_rov3_bnti %>% na.omit(),
                          pro_nti$asv16s_rna_rov4_bnti %>% na.omit(),
                          pro_nti$asv16s_rna_rov5_bnti %>% na.omit()),
                 Name = c(rep("dna_ROV1", pro_nti$asv16s_dna_rov1_bnti %>% na.omit() %>% length()),
                          rep("dna_ROV2", pro_nti$asv16s_dna_rov2_bnti %>% na.omit() %>% length()),
                          rep("dna_ROV3", pro_nti$asv16s_dna_rov3_bnti %>% na.omit() %>% length()),
                          rep("dna_ROV4", pro_nti$asv16s_dna_rov4_bnti %>% na.omit() %>% length()),
                          rep("dna_ROV5", pro_nti$asv16s_dna_rov5_bnti %>% na.omit() %>% length()),
                          rep("rna_ROV1", pro_nti$asv16s_rna_rov1_bnti %>% na.omit() %>% length()),
                          rep("rna_ROV2", pro_nti$asv16s_rna_rov2_bnti %>% na.omit() %>% length()),
                          rep("rna_ROV3", pro_nti$asv16s_rna_rov3_bnti %>% na.omit() %>% length()),
                          rep("rna_ROV4", pro_nti$asv16s_rna_rov4_bnti %>% na.omit() %>% length()),
                          rep("rna_ROV5", pro_nti$asv16s_rna_rov5_bnti %>% na.omit() %>% length())))


p2 <-
  df2 %>% as_tibble() %>% mutate(ROV = Name %>% str_split("_") %>% sapply('[', 2)) %>% 
  ggplot(aes(x = Name, y = bNTI))+
  geom_boxplot(outlier.colour = "white", outlier.fill = "white", size =1.5, aes(color = ROV))+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
              size = 2,aes(color = ROV), alpha = 0.5,show.legend = T)+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=15),
        axis.text.y = element_text(color="black", size=15))+
  theme(legend.title=element_text(size=15), legend.text=element_text(size=12),
        axis.title=element_text(size=15))+
  geom_hline(aes(yintercept = 2), color = "grey66", linetype = "dashed", linewidth = 1)+
  geom_hline(aes(yintercept = -2), color = "grey66", linetype = "dashed", linewidth = 1)+
  ylim(-22,5)+
  scale_color_manual(values = c("firebrick1", "hotpink", "darkgoldenrod1", "cornflowerblue","darkcyan"))
    
ggarrange(p1,p2,nrow=2)    



#18s community
asv_table_18s = read.table("asv_18s.txt",header = T,check.names = F) %>% t()
tax_18s = read.delim("taxonomy-pr2.tsv")

## remove metazoan
tax_pro = tax_18s[!c(str_detect(tax_18s$Taxon, "Metazoa")),]
tax_pro = tax_pro[!c(str_detect(tax_pro$Taxon, "Streptophyta")),]
asv_table_18s = asv_table_18s[,match(tax_pro$Feature.ID, colnames(asv_table_18s))]
asv_table_18s = rrarefy(asv_table_18s,
                        min(rowSums(asv_table_18s)))

# rarefy
asv_table_18s = rrarefy(asv_table_18s,
                        min(rowSums(asv_table_18s)))


## divide to DNA and RNA samples

asv18s_dna = 
  asv_table_18s[asv_table_18s %>% rownames() %>% str_detect("dna"), ]

asv18s_rna = 
  asv_table_18s[asv_table_18s %>% rownames() %>% str_detect("rna"), ]

# Make DNA samples are equal to RNA samples, as there are much more DNA samples
asv18s_rna %>% rownames()
m1 = asv18s_rna %>% rownames() %>% str_remove("_rna")
asv18s_dna %>% rownames()
m2 = asv18s_dna %>% rownames() %>% str_remove("_dna")

asv18s_dna = asv18s_dna[match(m1,m2),]
identical(rownames(asv18s_dna) %>% str_remove("_dna"),
          rownames(asv18s_rna) %>% str_remove("_rna")) # True

identical(colnames(asv18s_dna) %>% str_remove("_dna"),
          colnames(asv18s_rna) %>% str_remove("_rna")) # True


#16s community
asv_table_16s = read.table("asv_table_16s.txt",header = T,check.names = F) %>% t()
asv_table_16s %>% dim
tax_16s = read.delim("taxonomy-16s.tsv") # 108166 taxon
tax_16s = tax_16s[!c(str_detect(tax_16s$Taxon, "Chloroplast")),] # 108061 taxon remained
tax_16s = tax_16s[!c(str_detect(tax_16s$Taxon, "Mitochondria")),] #107915 taxon remained

asv_table_16s = asv_table_16s[,match(tax_16s$Feature.ID, colnames(asv_table_16s))]
asv_table_16s = rrarefy(asv_table_16s,
                        min(rowSums(asv_table_16s)))


## divide to DNA and RNA samples
asv16s_dna = 
  asv_table_16s[asv_table_16s %>% rownames() %>% str_detect("dna"), ]

asv16s_rna = 
  asv_table_16s[asv_table_16s %>% rownames() %>% str_detect("rna"), ]

# Make DNA samples are equal to RNA samples, as there are much more DNA samples
asv16s_rna %>% rownames()
m1 = asv16s_rna %>% rownames() %>% str_remove("_rna")
asv16s_dna %>% rownames()
m2 = asv16s_dna %>% rownames() %>% str_remove("_dna")

asv16s_dna = asv16s_dna[match(m1,m2),]
identical(rownames(asv16s_dna) %>% str_remove("_dna"),
          rownames(asv16s_rna) %>% str_remove("_rna")) # True

identical(colnames(asv16s_dna) %>% str_remove("_dna"),
          colnames(asv16s_rna) %>% str_remove("_rna")) # True

#separate by ROVs
asv18s_dna_rov1 = asv18s_dna[asv18s_dna %>% rownames() %>% str_detect("ROV1"), ]
asv18s_dna_rov2 = asv18s_dna[asv18s_dna %>% rownames() %>% str_detect("ROV2"), ]
asv18s_dna_rov3 = asv18s_dna[asv18s_dna %>% rownames() %>% str_detect("ROV3"), ]
asv18s_dna_rov4 = asv18s_dna[asv18s_dna %>% rownames() %>% str_detect("ROV4"), ]
asv18s_dna_rov5 = asv18s_dna[asv18s_dna %>% rownames() %>% str_detect("ROV5"), ]


asv18s_rna_rov1 = asv18s_rna[asv18s_rna %>% rownames() %>% str_detect("ROV1"), ]
asv18s_rna_rov2 = asv18s_rna[asv18s_rna %>% rownames() %>% str_detect("ROV2"), ]
asv18s_rna_rov3 = asv18s_rna[asv18s_rna %>% rownames() %>% str_detect("ROV3"), ]
asv18s_rna_rov4 = asv18s_rna[asv18s_rna %>% rownames() %>% str_detect("ROV4"), ]
asv18s_rna_rov5 = asv18s_rna[asv18s_rna %>% rownames() %>% str_detect("ROV5"), ]

asv16s_dna_rov1 = asv16s_dna[asv16s_dna %>% rownames() %>% str_detect("ROV1"), ]
asv16s_dna_rov2 = asv16s_dna[asv16s_dna %>% rownames() %>% str_detect("ROV2"), ]
asv16s_dna_rov3 = asv16s_dna[asv16s_dna %>% rownames() %>% str_detect("ROV3"), ]
asv16s_dna_rov4 = asv16s_dna[asv16s_dna %>% rownames() %>% str_detect("ROV4"), ]
asv16s_dna_rov5 = asv16s_dna[asv16s_dna %>% rownames() %>% str_detect("ROV5"), ]


asv16s_rna_rov1 = asv16s_rna[asv16s_rna %>% rownames() %>% str_detect("ROV1"), ]
asv16s_rna_rov2 = asv16s_rna[asv16s_rna %>% rownames() %>% str_detect("ROV2"), ]
asv16s_rna_rov3 = asv16s_rna[asv16s_rna %>% rownames() %>% str_detect("ROV3"), ]
asv16s_rna_rov4 = asv16s_rna[asv16s_rna %>% rownames() %>% str_detect("ROV4"), ]
asv16s_rna_rov5 = asv16s_rna[asv16s_rna %>% rownames() %>% str_detect("ROV5"), ]


d = list(asv18s_dna_rov1,asv18s_dna_rov2,asv18s_dna_rov3,asv18s_dna_rov4,asv18s_dna_rov5,
          asv18s_rna_rov1,asv18s_rna_rov2,asv18s_rna_rov3,asv18s_rna_rov4,asv18s_rna_rov5,
          asv16s_dna_rov1,asv16s_dna_rov2,asv16s_dna_rov3,asv16s_dna_rov4,asv16s_dna_rov5,
          asv16s_rna_rov1,asv16s_rna_rov2,asv16s_rna_rov3,asv16s_rna_rov4,asv16s_rna_rov5)

a <- lapply(d, function(x){x = x[,which(colSums(x)>0)]}) %>% 
  lapply(function(x){x = rrarefy(x, min(rowSums(x)))})


asv18s_dna_rov1 = a[[1]]
asv18s_dna_rov2 = a[[2]]
asv18s_dna_rov3 = a[[3]]
asv18s_dna_rov4 = a[[4]]
asv18s_dna_rov5 = a[[5]]
asv18s_rna_rov1 = a[[6]]
asv18s_rna_rov2 = a[[7]]
asv18s_rna_rov3 = a[[8]]
asv18s_rna_rov4 = a[[9]]
asv18s_rna_rov5 = a[[10]]
asv16s_dna_rov1 = a[[11]]
asv16s_dna_rov2 = a[[12]]
asv16s_dna_rov3 = a[[13]]
asv16s_dna_rov4 = a[[14]]
asv16s_dna_rov5 = a[[15]]
asv16s_rna_rov1 = a[[16]]
asv16s_rna_rov2 = a[[17]]
asv16s_rna_rov3 = a[[18]]
asv16s_rna_rov4 = a[[19]]
asv16s_rna_rov5 = a[[20]]


df = list(asv18s_dna_rov1,asv18s_dna_rov2,asv18s_dna_rov3,asv18s_dna_rov4,asv18s_dna_rov5,
          asv18s_rna_rov1,asv18s_rna_rov2,asv18s_rna_rov3,asv18s_rna_rov4,asv18s_rna_rov5,
          asv16s_dna_rov1,asv16s_dna_rov2,asv16s_dna_rov3,asv16s_dna_rov4,asv16s_dna_rov5,
          asv16s_rna_rov1,asv16s_rna_rov2,asv16s_rna_rov3,asv16s_rna_rov4,asv16s_rna_rov5)

names(df)  = c( "asv18s_dna_rov1","asv18s_dna_rov2","asv18s_dna_rov3","asv18s_dna_rov4","asv18s_dna_rov5",
                "asv18s_rna_rov1","asv18s_rna_rov2","asv18s_rna_rov3","asv18s_rna_rov4","asv18s_rna_rov5",
                "asv16s_dna_rov1","asv16s_dna_rov2","asv16s_dna_rov3","asv16s_dna_rov4","asv16s_dna_rov5",
                "asv16s_rna_rov1","asv16s_rna_rov2","asv16s_rna_rov3","asv16s_rna_rov4","asv16s_rna_rov5")

tree_18s  =read.tree("rooted_tree_18s.nwk");
tree_16s = read.tree("tree-16s-rooted.nwk");
tree = list(tree_18s,tree_18s,tree_18s,tree_18s,tree_18s,
            tree_18s,tree_18s,tree_18s,tree_18s,tree_18s,
            tree_16s,tree_16s,tree_16s,tree_16s,tree_16s,
            tree_16s,tree_16s,tree_16s,tree_16s,tree_16s)

#Calculate bNTI using Stegen et al. 2013 methods;
df_bnti = list()

for (k in seq_along(df)) {
  
  match.phylo.otu = match.phylo.data(tree[[k]], t(df[[k]]));
  str(match.phylo.otu);
  
  ## calculate empirical betaMNTD
  beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T));
  dim(beta.mntd.weighted);
  beta.mntd.weighted[1:5,1:5];
  #write.csv(beta.mntd.weighted,'betaMNTD_weighted.csv',quote=F);
  
  identical(colnames(match.phylo.otu$data),colnames(beta.mntd.weighted)); # just a check, should be TRUE
  identical(colnames(match.phylo.otu$data),rownames(beta.mntd.weighted)); # just a check, should be TRUE
  
  # calculate randomized betaMNTD
  
  beta.reps = 999; # number of randomizations
  
  rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.otu$data),ncol(match.phylo.otu$data),beta.reps));
  dim(rand.weighted.bMNTD.comp);
  
  for (rep in 1:beta.reps) {
    
    rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),taxaShuffle(cophenetic(match.phylo.otu$phy)),abundance.weighted=T,exclude.conspecifics = F));
    
    print(c(date(),rep));
    
  }
  
  weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data));
  dim(weighted.bNTI);
  
  for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
    for (rows in (columns+1):ncol(match.phylo.otu$data)) {
      
      rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
      weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
      rm("rand.vals");
      
    };
  };
  
  rownames(weighted.bNTI) = colnames(match.phylo.otu$data);
  colnames(weighted.bNTI) = colnames(match.phylo.otu$data);
  bnti = as.matrix(weighted.bNTI)
  bnti = bnti[lower.tri(bnti)]
  df_bnti[[k]] = bnti
}

asv18s_dna_rov1_bnti = df_bnti[[1]]
asv18s_dna_rov2_bnti = df_bnti[[2]]
asv18s_dna_rov3_bnti = df_bnti[[3]]
asv18s_dna_rov4_bnti = df_bnti[[4]]
asv18s_dna_rov5_bnti = df_bnti[[5]]
asv18s_rna_rov1_bnti = df_bnti[[6]]
asv18s_rna_rov2_bnti = df_bnti[[7]]
asv18s_rna_rov3_bnti = df_bnti[[8]]
asv18s_rna_rov4_bnti = df_bnti[[9]]
asv18s_rna_rov5_bnti = df_bnti[[10]]
asv16s_dna_rov1_bnti = df_bnti[[11]]
asv16s_dna_rov2_bnti = df_bnti[[12]]
asv16s_dna_rov3_bnti = df_bnti[[13]]
asv16s_dna_rov4_bnti = df_bnti[[14]]
asv16s_dna_rov5_bnti = df_bnti[[15]]
asv16s_rna_rov1_bnti = df_bnti[[16]]
asv16s_rna_rov2_bnti = df_bnti[[17]]
asv16s_rna_rov3_bnti = df_bnti[[18]]
asv16s_rna_rov4_bnti = df_bnti[[19]]
asv16s_rna_rov5_bnti = df_bnti[[20]]

write.csv(asv18s_dna_rov1_bnti,"asv18s_dna_rov1_bnti.csv",quote = F)
write.csv(asv18s_dna_rov2_bnti,"asv18s_dna_rov2_bnti.csv",quote = F)
write.csv(asv18s_dna_rov3_bnti,"asv18s_dna_rov3_bnti.csv",quote = F)
write.csv(asv18s_dna_rov4_bnti,"asv18s_dna_rov4_bnti.csv",quote = F)
write.csv(asv18s_dna_rov5_bnti,"asv18s_dna_rov5_bnti.csv",quote = F)
write.csv(asv18s_rna_rov1_bnti,"asv18s_rna_rov1_bnti.csv",quote = F)
write.csv(asv18s_rna_rov2_bnti,"asv18s_rna_rov2_bnti.csv",quote = F)
write.csv(asv18s_rna_rov3_bnti,"asv18s_rna_rov3_bnti.csv",quote = F)
write.csv(asv18s_rna_rov4_bnti,"asv18s_rna_rov4_bnti.csv",quote = F)
write.csv(asv18s_rna_rov5_bnti,"asv18s_rna_rov5_bnti.csv",quote = F)


write.csv(asv16s_dna_rov1_bnti,"asv16s_dna_rov1_bnti.csv",quote = F)
write.csv(asv16s_dna_rov2_bnti,"asv16s_dna_rov2_bnti.csv",quote = F)
write.csv(asv16s_dna_rov3_bnti,"asv16s_dna_rov3_bnti.csv",quote = F)
write.csv(asv16s_dna_rov4_bnti,"asv16s_dna_rov4_bnti.csv",quote = F)
write.csv(asv16s_dna_rov5_bnti,"asv16s_dna_rov5_bnti.csv",quote = F)
write.csv(asv16s_rna_rov1_bnti,"asv16s_rna_rov1_bnti.csv",quote = F)
write.csv(asv16s_rna_rov2_bnti,"asv16s_rna_rov2_bnti.csv",quote = F)
write.csv(asv16s_rna_rov3_bnti,"asv16s_rna_rov3_bnti.csv",quote = F)
write.csv(asv16s_rna_rov4_bnti,"asv16s_rna_rov4_bnti.csv",quote = F)
write.csv(asv16s_rna_rov5_bnti,"asv16s_rna_rov5_bnti.csv",quote = F)

boxplot(asv18s_dna_rov1_bnti,asv18s_dna_rov2_bnti,asv18s_dna_rov3_bnti,asv18s_dna_rov4_bnti,asv18s_dna_rov5_bnti,
        asv18s_rna_rov1_bnti,asv18s_rna_rov2_bnti,asv18s_rna_rov3_bnti,asv18s_rna_rov4_bnti,asv18s_rna_rov5_bnti,
        asv16s_dna_rov1_bnti,asv16s_dna_rov2_bnti,asv16s_dna_rov3_bnti,asv16s_dna_rov4_bnti,asv16s_dna_rov5_bnti),
        asv16s_rna_rov1_bnti,asv16s_rna_rov2_bnti,asv16s_rna_rov3_bnti,asv16s_rna_rov4_bnti,asv16s_rna_rov5_bnti,
        ylab="??NTI")


# betaNTI and ENV ---------------------------------------------------------
b1 =  read.csv("asv16s_dna_rov1_bnti.csv", row.names = 1)
b2 =  read.csv("asv16s_dna_rov2_bnti.csv", row.names = 1)
b3 =  read.csv("asv16s_dna_rov3_bnti.csv", row.names = 1)
b4 =  read.csv("asv16s_dna_rov4_bnti.csv", row.names = 1)
b5 =  read.csv("asv16s_dna_rov5_bnti.csv", row.names = 1)

a1 = 
asv16s_dna_rov1 %>% vegdist() %>% as.matrix() %>% data.frame() %>% 
  rownames_to_column(., var = "Sample1") %>% as_tibble() %>% 
  pivot_longer(!Sample1, names_to = "Sample2", values_to = "BC_distance") %>% 
  filter(Sample1 != Sample2) %>% 
  mutate(Sample1 = Sample1 %>% str_remove("_16s_dna")) %>% 
  mutate(Sample2 = Sample2 %>% str_remove("_16s_dna")) 

a2 = 
  asv16s_dna_rov2 %>% vegdist() %>% as.matrix() %>% data.frame() %>% 
  rownames_to_column(., var = "Sample1") %>% as_tibble() %>% 
  pivot_longer(!Sample1, names_to = "Sample2", values_to = "BC_distance") %>% 
  filter(Sample1 != Sample2) %>% 
  mutate(Sample1 = Sample1 %>% str_remove("_16s_dna")) %>% 
  mutate(Sample2 = Sample2 %>% str_remove("_16s_dna"))

a3 = 
  asv16s_dna_rov3 %>% vegdist() %>% as.matrix() %>% data.frame() %>% 
  rownames_to_column(., var = "Sample1") %>% as_tibble() %>% 
  pivot_longer(!Sample1, names_to = "Sample2", values_to = "BC_distance") %>% 
  filter(Sample1 != Sample2) %>% 
  mutate(Sample1 = Sample1 %>% str_remove("_16s_dna")) %>% 
  mutate(Sample2 = Sample2 %>% str_remove("_16s_dna")) 

a4 = 
  asv16s_dna_rov4 %>% vegdist() %>% as.matrix() %>% data.frame() %>% 
  rownames_to_column(., var = "Sample1") %>% as_tibble() %>% 
  pivot_longer(!Sample1, names_to = "Sample2", values_to = "BC_distance") %>% 
  filter(Sample1 != Sample2) %>% 
  mutate(Sample1 = Sample1 %>% str_remove("_16s_dna")) %>% 
  mutate(Sample2 = Sample2 %>% str_remove("_16s_dna")) 

a5 = 
  asv16s_dna_rov5 %>% vegdist() %>% as.matrix() %>% data.frame() %>% 
  rownames_to_column(., var = "Sample1") %>% as_tibble() %>% 
  pivot_longer(!Sample1, names_to = "Sample2", values_to = "BC_distance") %>% 
  filter(Sample1 != Sample2) %>% 
  mutate(Sample1 = Sample1 %>% str_remove("_16s_dna")) %>% 
  mutate(Sample2 = Sample2 %>% str_remove("_16s_dna")) 

df = rbind(a1,a2,a3,a4,a5)
  
  
env = read.csv("SEM_data_for_16S.csv")
env %>% colnames()
env = env[,c("Sample", "Methane")]

df %>% left_join(env %>% rename(Sample1 = Sample) %>% rename(Methane1 = Methane),
                 by = "Sample1") %>% 
  left_join(env %>% rename(Sample2 = Sample) %>% rename(Methane2 = Methane),
            by = "Sample2") %>% 
  drop_na() %>% 
  mutate(Methane_dif  = abs(Methane1 - Methane2)) %>% 
  mutate(Methane_av = (Methane1 + Methane2)/2) %>% 
  ggplot(aes(x = log(Methane_av,10), y = BC_distance))+
  geom_point()+
  geom_smooth(method = "lm", se=T)+                                    #### set: se = T
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size=15),
        axis.text.y = element_text(color="black", size=15))+
  theme(legend.title=element_text(size=20), legend.text=element_text(size=15),
        axis.title=element_text(size=20))+
  ggpubr::stat_cor(label.y= 0.3 ,method = "spearman", size = 6, color = "red")+
  xlab("Methane concentration (Log10)")+
  ylab("Community dissimilarity")



env = read.csv("SEM_data_for_16S.csv")
env %>% colnames()
env = env[,c("Sample", "Sulfide")]
df %>% left_join(env %>% rename(Sample1 = Sample) %>% rename(Sulfide1 = Sulfide),
                 by = "Sample1") %>% 
  left_join(env %>% rename(Sample2 = Sample) %>% rename(Sulfide2 = Sulfide),
            by = "Sample2") %>% 
  drop_na() %>% 
  mutate(Sulfide_dif  = abs(Sulfide1 - Sulfide2)) %>% 
  mutate(Sulfide_av = (Sulfide1 + Sulfide2)/2) %>% 
  ggplot(aes(x = log(Sulfide_av,10), y = BC_distance))+
  geom_point()+
  geom_smooth(method = "lm", se=T)+                                    #### set: se = T
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size=15),
        axis.text.y = element_text(color="black", size=15))+
  theme(legend.title=element_text(size=20), legend.text=element_text(size=15),
        axis.title=element_text(size=20))+
  ggpubr::stat_cor(label.y= 0.3 ,method = "spearman", size = 6, color = "red")+
  xlab("Sulfide concentration (Log10)")+
  ylab("Community dissimilarity")


env = read.csv("SEM_data_for_16S.csv")
env %>% colnames()
env = env[,c("Sample", "DIC")]

df %>% left_join(env %>% rename(Sample1 = Sample) %>% rename(DIC1 = DIC),
                 by = "Sample1") %>% 
  left_join(env %>% rename(Sample2 = Sample) %>% rename(DIC2 = DIC),
            by = "Sample2") %>% 
  drop_na() %>% 
  mutate(DIC_dif  = abs(DIC1 - DIC2)) %>% 
  mutate(DIC_av = (DIC1 + DIC2)/2) %>% 
  ggplot(aes(x = log(DIC_av,10), y = BC_distance))+
  geom_point()+
  geom_smooth(method = "lm", se=T)+                                    #### set: se = T
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size=15),
        axis.text.y = element_text(color="black", size=15))+
  theme(legend.title=element_text(size=20), legend.text=element_text(size=15),
        axis.title=element_text(size=20))+
  ggpubr::stat_cor(label.y= 0.3 ,method = "spearman", size = 6, color = "red")+
  xlab("DIC concentration (Log10)")+
  ylab("Community dissimilarity")

# NTI ---------------------------------------------------------------------

match.phylo.otu1 = match.phylo.data(tree_18s, t(asv18s_dna_rov1))
match.phylo.otu2 = match.phylo.data(tree_18s, t(asv18s_dna_rov2))
match.phylo.otu3 = match.phylo.data(tree_18s, t(asv18s_dna_rov3))
match.phylo.otu4 = match.phylo.data(tree_18s, t(asv18s_dna_rov4))
match.phylo.otu5 = match.phylo.data(tree_18s, t(asv18s_dna_rov5))
match.phylo.otu6 = match.phylo.data(tree_18s, t(asv18s_rna_rov1))
match.phylo.otu7 = match.phylo.data(tree_18s, t(asv18s_rna_rov2))
match.phylo.otu8 = match.phylo.data(tree_18s, t(asv18s_rna_rov3))
match.phylo.otu9 = match.phylo.data(tree_18s, t(asv18s_rna_rov4))
match.phylo.otu10 = match.phylo.data(tree_18s, t(asv18s_rna_rov5))

match.phylo.otu11 = match.phylo.data(tree_16s, t(asv16s_dna_rov1))
match.phylo.otu12 = match.phylo.data(tree_16s, t(asv16s_dna_rov2))
match.phylo.otu13 = match.phylo.data(tree_16s, t(asv16s_dna_rov3))
match.phylo.otu14 = match.phylo.data(tree_16s, t(asv16s_dna_rov4))
match.phylo.otu15 = match.phylo.data(tree_16s, t(asv16s_dna_rov5))
match.phylo.otu16 = match.phylo.data(tree_16s, t(asv16s_rna_rov1))
match.phylo.otu17 = match.phylo.data(tree_16s, t(asv16s_rna_rov2))
match.phylo.otu18 = match.phylo.data(tree_16s, t(asv16s_rna_rov3))
match.phylo.otu19 = match.phylo.data(tree_16s, t(asv16s_rna_rov4))
match.phylo.otu20 = match.phylo.data(tree_16s, t(asv16s_rna_rov5))

NTI_1 <- ses.mntd(t(match.phylo.otu1$data), cophenetic(match.phylo.otu1$phy),null.model="taxa.labels",runs=999)
NTI_2 <- ses.mntd(t(match.phylo.otu2$data), cophenetic(match.phylo.otu2$phy),null.model="taxa.labels",runs=999)
NTI_3 <- ses.mntd(t(match.phylo.otu3$data), cophenetic(match.phylo.otu3$phy),null.model="taxa.labels",runs=999)
NTI_4 <- ses.mntd(t(match.phylo.otu4$data), cophenetic(match.phylo.otu4$phy),null.model="taxa.labels",runs=999)
NTI_5 <- ses.mntd(t(match.phylo.otu5$data), cophenetic(match.phylo.otu5$phy),null.model="taxa.labels",runs=999)

NTI_6 <- ses.mntd(t(match.phylo.otu6$data), cophenetic(match.phylo.otu6$phy),null.model="taxa.labels",runs=999)
NTI_7 <- ses.mntd(t(match.phylo.otu7$data), cophenetic(match.phylo.otu7$phy),null.model="taxa.labels",runs=999)
NTI_8 <- ses.mntd(t(match.phylo.otu8$data), cophenetic(match.phylo.otu8$phy),null.model="taxa.labels",runs=999)
NTI_9 <- ses.mntd(t(match.phylo.otu9$data), cophenetic(match.phylo.otu9$phy),null.model="taxa.labels",runs=999)
NTI_10 <- ses.mntd(t(match.phylo.otu10$data), cophenetic(match.phylo.otu10$phy),null.model="taxa.labels",runs=999)

NTI_11 <- ses.mntd(t(match.phylo.otu11$data), cophenetic(match.phylo.otu11$phy),null.model="taxa.labels",runs=999)
NTI_12 <- ses.mntd(t(match.phylo.otu12$data), cophenetic(match.phylo.otu12$phy),null.model="taxa.labels",runs=999)
NTI_13 <- ses.mntd(t(match.phylo.otu13$data), cophenetic(match.phylo.otu13$phy),null.model="taxa.labels",runs=999)
NTI_14 <- ses.mntd(t(match.phylo.otu14$data), cophenetic(match.phylo.otu14$phy),null.model="taxa.labels",runs=999)
NTI_15 <- ses.mntd(t(match.phylo.otu15$data), cophenetic(match.phylo.otu15$phy),null.model="taxa.labels",runs=999)

NTI_16 <- ses.mntd(t(match.phylo.otu16$data), cophenetic(match.phylo.otu16$phy),null.model="taxa.labels",runs=999)
NTI_17 <- ses.mntd(t(match.phylo.otu17$data), cophenetic(match.phylo.otu17$phy),null.model="taxa.labels",runs=999)
NTI_18 <- ses.mntd(t(match.phylo.otu18$data), cophenetic(match.phylo.otu18$phy),null.model="taxa.labels",runs=999)
NTI_19 <- ses.mntd(t(match.phylo.otu19$data), cophenetic(match.phylo.otu19$phy),null.model="taxa.labels",runs=999)
NTI_20 <- ses.mntd(t(match.phylo.otu20$data), cophenetic(match.phylo.otu20$phy),null.model="taxa.labels",runs=999)


boxplot(-NTI_1$mntd.obs.z, -NTI_2$mntd.obs.z, -NTI_3$mntd.obs.z, -NTI_4$mntd.obs.z, -NTI_5$mntd.obs.z,
        -NTI_6$mntd.obs.z, -NTI_7$mntd.obs.z, -NTI_8$mntd.obs.z, -NTI_9$mntd.obs.z, -NTI_10$mntd.obs.z,
        -NTI_11$mntd.obs.z, -NTI_12$mntd.obs.z, -NTI_13$mntd.obs.z, -NTI_14$mntd.obs.z, -NTI_15$mntd.obs.z,
        -NTI_16$mntd.obs.z, -NTI_17$mntd.obs.z, -NTI_18$mntd.obs.z, -NTI_19$mntd.obs.z, -NTI_20$mntd.obs.z, 
        ylab = "NTI")

#NRI

NRI_1 <- ses.mpd(t(match.phylo.otu1$data), cophenetic(match.phylo.otu1$phy),null.model="taxa.labels",runs=999)
NRI_2 <- ses.mpd(t(match.phylo.otu2$data), cophenetic(match.phylo.otu2$phy),null.model="taxa.labels",runs=999)
NRI_3 <- ses.mpd(t(match.phylo.otu3$data), cophenetic(match.phylo.otu3$phy),null.model="taxa.labels",runs=999)
NRI_4 <- ses.mpd(t(match.phylo.otu4$data), cophenetic(match.phylo.otu4$phy),null.model="taxa.labels",runs=999)
NRI_5 <- ses.mpd(t(match.phylo.otu5$data), cophenetic(match.phylo.otu5$phy),null.model="taxa.labels",runs=999)

NRI_6 <- ses.mpd(t(match.phylo.otu6$data), cophenetic(match.phylo.otu6$phy),null.model="taxa.labels",runs=999)
NRI_7 <- ses.mpd(t(match.phylo.otu7$data), cophenetic(match.phylo.otu7$phy),null.model="taxa.labels",runs=999)
NRI_8 <- ses.mpd(t(match.phylo.otu8$data), cophenetic(match.phylo.otu8$phy),null.model="taxa.labels",runs=999)
NRI_9 <- ses.mpd(t(match.phylo.otu9$data), cophenetic(match.phylo.otu9$phy),null.model="taxa.labels",runs=999)
NRI_10 <- ses.mpd(t(match.phylo.otu10$data), cophenetic(match.phylo.otu10$phy),null.model="taxa.labels",runs=999)

NRI_11 <- ses.mpd(t(match.phylo.otu11$data), cophenetic(match.phylo.otu11$phy),null.model="taxa.labels",runs=999)
NRI_12 <- ses.mpd(t(match.phylo.otu12$data), cophenetic(match.phylo.otu12$phy),null.model="taxa.labels",runs=999)
NRI_13 <- ses.mpd(t(match.phylo.otu13$data), cophenetic(match.phylo.otu13$phy),null.model="taxa.labels",runs=999)
NRI_14 <- ses.mpd(t(match.phylo.otu14$data), cophenetic(match.phylo.otu14$phy),null.model="taxa.labels",runs=999)
NRI_15 <- ses.mpd(t(match.phylo.otu15$data), cophenetic(match.phylo.otu15$phy),null.model="taxa.labels",runs=999)

NRI_16 <- ses.mpd(t(match.phylo.otu16$data), cophenetic(match.phylo.otu16$phy),null.model="taxa.labels",runs=999)
NRI_17 <- ses.mpd(t(match.phylo.otu17$data), cophenetic(match.phylo.otu17$phy),null.model="taxa.labels",runs=999)
NRI_18 <- ses.mpd(t(match.phylo.otu18$data), cophenetic(match.phylo.otu18$phy),null.model="taxa.labels",runs=999)
NRI_19 <- ses.mpd(t(match.phylo.otu19$data), cophenetic(match.phylo.otu19$phy),null.model="taxa.labels",runs=999)
NRI_20 <- ses.mpd(t(match.phylo.otu20$data), cophenetic(match.phylo.otu20$phy),null.model="taxa.labels",runs=999)


boxplot(-NRI_1$mpd.obs.z, -NRI_2$mpd.obs.z, -NRI_3$mpd.obs.z, -NRI_4$mpd.obs.z, -NRI_5$mpd.obs.z,
        -NRI_6$mpd.obs.z, -NRI_7$mpd.obs.z, -NRI_8$mpd.obs.z, -NRI_9$mpd.obs.z, -NRI_10$mpd.obs.z,
        -NRI_11$mpd.obs.z, -NRI_12$mpd.obs.z, -NRI_13$mpd.obs.z, -NRI_14$mpd.obs.z, -NRI_15$mpd.obs.z,
        -NRI_16$mpd.obs.z, -NRI_17$mpd.obs.z, -NRI_18$mpd.obs.z, -NRI_19$mpd.obs.z, -NRI_20$mpd.obs.z,
        ylab = "NRI")


# gamma diversity ---------------------------------------------------------
#random
#18s
#dna
curve_asv18s_dna_rov1 = specaccum(asv18s_dna_rov1, method = "random")
curve_asv18s_dna_rov2 = specaccum(asv18s_dna_rov2, method = "random")
curve_asv18s_dna_rov3 = specaccum(asv18s_dna_rov3, method = "random")
curve_asv18s_dna_rov4 = specaccum(asv18s_dna_rov4, method = "random")
curve_asv18s_dna_rov5 = specaccum(asv18s_dna_rov5, method = "random")
# m1 <- fitspecaccum(curve_asv18s_dna_rov1, "lomolino")
plot(curve_asv18s_dna_rov1,col = "firebrick1",lwd = 3, main = "18S-DNA",
     xlab = "Number of sites", ylab = "Number of sequences")
plot(curve_asv18s_dna_rov2, add = TRUE, col = "hotpink",lwd = 3)
plot(curve_asv18s_dna_rov3, add = TRUE, col = "darkgoldenrod1",lwd = 3)
plot(curve_asv18s_dna_rov4, add = TRUE, col = "cornflowerblue",lwd = 3)
plot(curve_asv18s_dna_rov5, add = TRUE, col = "darkcyan",lwd = 3)
#18s
#rna
curve_asv18s_rna_rov1 = specaccum(asv18s_rna_rov1, method = "random")
curve_asv18s_rna_rov2 = specaccum(asv18s_rna_rov2, method = "random")
curve_asv18s_rna_rov3 = specaccum(asv18s_rna_rov3, method = "random")
curve_asv18s_rna_rov4 = specaccum(asv18s_rna_rov4, method = "random")
curve_asv18s_rna_rov5 = specaccum(asv18s_rna_rov5, method = "random")
# m1 <- fitspecaccum(curve_asv18s_rna_rov1, "lomolino")
plot(curve_asv18s_rna_rov1,col = "firebrick1",lwd = 3, main = "18S-RNA",
     xlab = "Number of sites", ylab = "Number of sequences")
plot(curve_asv18s_rna_rov2, add = TRUE, col = "hotpink",lwd = 3)
plot(curve_asv18s_rna_rov3, add = TRUE, col = "darkgoldenrod1",lwd = 3)
plot(curve_asv18s_rna_rov4, add = TRUE, col = "cornflowerblue",lwd = 3)
plot(curve_asv18s_rna_rov5, add = TRUE, col = "darkcyan",lwd = 3)

#16s
#dna
curve_asv16s_dna_rov1 = specaccum(asv16s_dna_rov1, method = "random")
curve_asv16s_dna_rov2 = specaccum(asv16s_dna_rov2, method = "random")
curve_asv16s_dna_rov3 = specaccum(asv16s_dna_rov3, method = "random")
curve_asv16s_dna_rov4 = specaccum(asv16s_dna_rov4, method = "random")
curve_asv16s_dna_rov5 = specaccum(asv16s_dna_rov5, method = "random")
# m1 <- fitspecaccum(curve_asv16s_dna_rov1, "lomolino")
plot(curve_asv16s_dna_rov1,col = "firebrick1",lwd = 3,main = "16S-DNA",
     xlab = "Number of sites", ylab = "Number of sequences")
plot(curve_asv16s_dna_rov2, add = TRUE, col = "hotpink",lwd = 3)
plot(curve_asv16s_dna_rov3, add = TRUE, col = "darkgoldenrod1",lwd = 3)
plot(curve_asv16s_dna_rov4, add = TRUE, col = "cornflowerblue",lwd = 3)
plot(curve_asv16s_dna_rov5, add = TRUE, col = "darkcyan",lwd = 3)

#16s
#rna
curve_asv16s_rna_rov1 = specaccum(asv16s_rna_rov1, method = "random")
curve_asv16s_rna_rov2 = specaccum(asv16s_rna_rov2, method = "random")
curve_asv16s_rna_rov3 = specaccum(asv16s_rna_rov3, method = "random")
curve_asv16s_rna_rov4 = specaccum(asv16s_rna_rov4, method = "random")
curve_asv16s_rna_rov5 = specaccum(asv16s_rna_rov5, method = "random")
# m1 <- fitspecaccum(curve_asv16s_rna_rov1, "lomolino")
plot(curve_asv16s_rna_rov1,col = "firebrick1",lwd = 3, main = "16S-RNA",
     xlab = "Number of sites", ylab = "Number of sequences")
plot(curve_asv16s_rna_rov2, add = TRUE, col = "hotpink",lwd = 3)
plot(curve_asv16s_rna_rov3, add = TRUE, col = "darkgoldenrod1",lwd = 3)
plot(curve_asv16s_rna_rov4, add = TRUE, col = "cornflowerblue",lwd = 3)
plot(curve_asv16s_rna_rov5, add = TRUE, col = "darkcyan",lwd = 3)


# ## 13C 15N isotope ------------------------------------------------------

setwd("C:/Users/XZM/Desktop/Coldseep 2nd study/2021 C13 N15 measured by HKU")
cn = read.csv("2021coldseepCN.csv",check.names = F)

# x = 13, y = 15N
library(ggExtra)
p1<-
cn %>% mutate(ROV = str_extract(Sample_name, "ROV\\d")) %>% 
  mutate(Depth_m = str_split(Sample_name,"-") %>% sapply('[', 2)) %>%
  filter(Depth_m %in% c(0,5,10)) %>% 
  ggplot(aes(x = d13C, y = d15N,color = ROV))+
  geom_point(size = 5)+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth = 2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=20),
        axis.text.y = element_text(color="black", size=20))+
  theme(legend.title=element_text(size=20), legend.text=element_text(size=20),
        axis.title=element_text(size=20),
        legend.position = c(0.15,0.85))+
  ylab(expression(paste(delta^15*N)))+
  xlab(expression(paste(delta^13*C)))+
  scale_color_manual(values = c("firebrick1", "hotpink", "darkgoldenrod1", "cornflowerblue","darkcyan"))
  
ggMarginal(p2, groupColour = TRUE, groupFill = TRUE)

# x = %C, y = %N

p2<-
cn %>% mutate(ROV = str_extract(Sample_name, "ROV\\d")) %>% 
  mutate(Depth_m = str_split(Sample_name,"-") %>% sapply('[', 2)) %>%
  filter(Depth_m %in% c(0,5,10)) %>% 
  ggplot(aes(x = `%C`, y = `%N`,color = ROV))+
  geom_point(size = 5)+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth = 2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=20),
        axis.text.y = element_text(color="black", size=20),
        legend.position = c(0.15,0.85))+
  theme(legend.title=element_text(size=20), legend.text=element_text(size=15),
        axis.title=element_text(size=20))+
  scale_color_manual(values = c("firebrick1", "hotpink", "darkgoldenrod1", "cornflowerblue","darkcyan"))

ggarrange(p1,p2,nrow = 1)

# x= C/N, y = 13C
df %>% mutate(ROV = str_extract(Sample_name, "ROV\\d")) %>% 
  mutate(Depth_m = str_split(Sample_name,"-") %>% sapply('[', 2)) %>%
  mutate(`C/N` = `%C` / `%N`) %>% 
  filter(Depth_m %in% c(0,5,10)) %>% 
  ggplot(aes(x = `C/N`, y = d13C,color = ROV))+
  geom_point(size = 5)+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth = 2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=15),
        axis.text.y = element_text(color="black", size=15))+
  theme(legend.title=element_text(size=20), legend.text=element_text(size=15),
        axis.title=element_text(size=20))+
  xlab(expression(paste(C/N)))+
  ylab(expression(paste(delta^13*C)))


# alpha diversity and C N
#Richness
df %>%
  mutate(ROV = str_extract(Sample_name, "ROV\\d")) %>% 
  mutate(Depth_m = str_split(Sample_name,"-") %>% sapply('[', 2) %>% as.numeric()) %>%
  mutate(`C/N` = `%C` / `%N`) %>% 
  mutate(Sample_name = Sample_name %>% str_replace("-","_")) %>% 
  mutate(Richness_18s_dna = 
           asv_18s_apha$Richness[match(df$Sample_name %>% str_replace("-","_"), 
                                       asv_18s_apha %>% rownames_to_column(., var = "sample") %>% as_tibble() %>% 
                                      filter(sample %>% str_detect("dna")) %>% pull(sample) %>% 
                                      str_split("_18s_") %>% sapply('[', 1))],
         Richness_18s_rna = 
           asv_18s_apha$Richness[match(df$Sample_name %>% str_replace("-","_"), 
                                       asv_18s_apha %>% rownames_to_column(., var = "sample") %>% as_tibble() %>% 
                                         filter(sample %>% str_detect("rna")) %>% pull(sample) %>% 
                                         str_split("_18s_") %>% sapply('[', 1))],
         Richness_16s_dna = 
           asv_16s_apha$Richness[match(df$Sample_name %>% str_replace("-","_"), 
                                       asv_16s_apha %>% rownames_to_column(., var = "sample") %>% as_tibble() %>% 
                                         filter(sample %>% str_detect("dna")) %>% pull(sample) %>% 
                                         str_split("_16s_") %>% sapply('[', 1))],
         Richness_16s_rna = 
           asv_16s_apha$Richness[match(df$Sample_name %>% str_replace("-","_"), 
                                       asv_16s_apha %>% rownames_to_column(., var = "sample") %>% as_tibble() %>% 
                                         filter(sample %>% str_detect("rna")) %>% pull(sample) %>% 
                                         str_split("_16s_") %>% sapply('[', 1))]) %>% 
  column_to_rownames(., var = "Sample_name") %>%
  select(-ROV) %>% 
  filter(!is.na(Richness_18s_rna)) %>% 
  filter(!is.na(Richness_16s_rna)) -> df1

M = cor(df1)  
res1 = cor.mtest(df1, conf.level = .95)  
corrplot(M, p.mat = res1$p, sig.level = 0.05)
  
         
# PD
df %>%
  mutate(ROV = str_extract(Sample_name, "ROV\\d")) %>% 
  mutate(Depth_m = str_split(Sample_name,"-") %>% sapply('[', 2) %>% as.numeric()) %>%
  mutate(`C/N` = `%C` / `%N`) %>% 
  mutate(Sample_name = Sample_name %>% str_replace("-","_")) %>% 
  mutate(Pielou_18s_dna = 
           asv_18s_apha$Pielou[match(df$Sample_name %>% str_replace("-","_"), 
                                       asv_18s_apha %>% rownames_to_column(., var = "sample") %>% as_tibble() %>% 
                                         filter(sample %>% str_detect("dna")) %>% pull(sample) %>% 
                                         str_split("_18s_") %>% sapply('[', 1))],
         Pielou_18s_rna = 
           asv_18s_apha$Pielou[match(df$Sample_name %>% str_replace("-","_"), 
                                       asv_18s_apha %>% rownames_to_column(., var = "sample") %>% as_tibble() %>% 
                                         filter(sample %>% str_detect("rna")) %>% pull(sample) %>% 
                                         str_split("_18s_") %>% sapply('[', 1))],
         Pielou_16s_dna = 
           asv_16s_apha$Pielou[match(df$Sample_name %>% str_replace("-","_"), 
                                       asv_16s_apha %>% rownames_to_column(., var = "sample") %>% as_tibble() %>% 
                                         filter(sample %>% str_detect("dna")) %>% pull(sample) %>% 
                                         str_split("_16s_") %>% sapply('[', 1))],
         Pielou_16s_rna = 
           asv_16s_apha$Pielou[match(df$Sample_name %>% str_replace("-","_"), 
                                       asv_16s_apha %>% rownames_to_column(., var = "sample") %>% as_tibble() %>% 
                                         filter(sample %>% str_detect("rna")) %>% pull(sample) %>% 
                                         str_split("_16s_") %>% sapply('[', 1))]) %>% 
  column_to_rownames(., var = "Sample_name") %>%
  select(-ROV) %>% 
  filter(!is.na(Pielou_18s_rna)) %>% 
  filter(!is.na(Pielou_16s_rna)) -> df1

M = cor(df1)  
res1 = cor.mtest(df1, conf.level = .95)  
corrplot(M, p.mat = res1$p, sig.level = 0.05)




# correlation between alpha diversity and all env factors-----------------------------------------

asv_18s = read.table("asv_18s.txt",header = T,check.names = F) %>% t()
tax_18s = read.delim("taxonomy-18s.tsv")
## remove metazoan and land plants
tax_pro = tax_18s[!c(str_detect(tax_18s$Taxon, "Metazoa")),]
tax_pro = tax_pro[!c(str_detect(tax_pro$Taxon, "Streptophyta")),]
asv_18s = asv_18s[,match(tax_pro$Feature.ID, colnames(asv_18s))]
asv_18s = rrarefy(asv_18s,min(rowSums(asv_18s)))

asv_16s = read.table("asv_table_16s.txt",header = T,check.names = F) %>% t()
asv_16s %>% dim
asv16s_dna = asv_16s[asv_16s %>% rownames() %>% str_detect("dna"), ]
asv16s_dna %>% rownames()
asv16s_dna = asv16s_dna[,which(colSums(asv16s_dna)>0)]
asv16s_dna %>% dim
asv16s_dna = rrarefy(asv16s_dna, min(rowSums(asv16s_dna)))
asv16s_dna = asv16s_dna[,which(colSums(asv16s_dna)>0)]

tree_18s = read.tree("rooted_tree_18s.nwk")
tree_16s = read.tree("tree-16s-rooted.nwk")

asv_18s_apha = alpha(asv_18s, tree_18s)
asv_16s_apha = alpha(asv_16s, tree_16s)


asv_16s_apha_meta  <-
asv_16s_apha %>% rownames_to_column(., var = "Sample_16s") %>% as_tibble() %>% 
  mutate(ROV_depth = Sample_16s %>% str_split("_16s") %>% sapply('[', 1)) %>% 
  rename(c(Richness_16s = Richness,
           Shannon_16s = Shannon,
           Simpson_16s = Simpson,
           Pielou_16s = Pielou,
           Chao1_16s = Chao1,
           ACE_16s = ACE,
           goods_coverage_16s = goods_coverage,
           PD_whole_tree_16s = PD_whole_tree)) %>% 
  mutate(nucleid = Sample_16s %>% str_split("_16s_") %>% sapply('[', 2)) %>% 
  mutate(ROV = Sample_16s %>% str_sub(5,8)) %>% 
  mutate(Depth = Sample_16s %>% str_split("_") %>% sapply('[', 2)) 
  

asv_18s_apha_meta <-
asv_18s_apha %>% rownames_to_column(., var = "Sample_18s") %>% as_tibble() %>% 
  mutate(ROV_depth = Sample_18s %>% str_split("_18s") %>% sapply('[', 1)) %>% 
  rename(c(Richness_18s = Richness,
           Shannon_18s = Shannon,
           Simpson_18s = Simpson,
           Pielou_18s = Pielou,
           Chao1_18s = Chao1,
           ACE_18s = ACE,
           goods_coverage_18s = goods_coverage,
           PD_whole_tree_18s = PD_whole_tree)) %>% 
  mutate(nucleid = Sample_18s %>% str_split("_18s_") %>% sapply('[', 2)) %>% 
  mutate(ROV = Sample_18s %>% str_sub(5,8)) %>% 
  mutate(Depth = Sample_18s %>% str_split("_") %>% sapply('[', 2)) %>% 
  filter(Depth %in% c("0","5","10"))


df_dna = merge(asv_16s_apha_meta[which(asv_16s_apha_meta$nucleid == "dna"),],
               asv_18s_apha_meta[which(asv_18s_apha_meta$nucleid == "dna"),], 
               by.x = "ROV_depth", by.y = "ROV_depth")

df %>% as_tibble()

df_dna %>% ggplot(aes(x = Richness_16s, y = Richness_18s))+
  geom_point(aes(color = ROV.x), size = 4)+
  geom_smooth(method = "lm", se=T)+                                    #### set: se = T
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size=15, hjust = 1),
        axis.text.y = element_text(color="black", size=15))+
  theme(legend.title=element_text(size=20), legend.text=element_text(size=15),
        axis.title=element_text(size=20))+
  ggpubr::stat_cor(label.y= 400 ,method = "spearman", size = 5.3)



df_rna = merge(asv_16s_apha_meta[which(asv_16s_apha_meta$nucleid == "rna"),],
               asv_18s_apha_meta[which(asv_18s_apha_meta$nucleid == "rna"),], 
               by.x = "ROV_depth", by.y = "ROV_depth")


df_rna %>% ggplot(aes(x = Richness_16s, y = Richness_18s))+
  geom_point(aes(color = ROV.x), size = 4)+
  geom_smooth(method = "lm")+                                    #### set: se = T
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size=15, hjust = 1),
        axis.text.y = element_text(color="black", size=15))+
  theme(legend.title=element_text(size=20), legend.text=element_text(size=15),
        axis.title=element_text(size=20))+
  ggpubr::stat_cor(label.y= 400 ,method = "spearman", size = 5.3)

## combine with 13C 15N

setwd("C:/Users/XZM/Desktop/Coldseep 2nd study/2021 C13 N15 measured by HKU")
cn = read.csv("2021coldseepCN.csv",check.names = F)

cn<-
cn %>% mutate(ROV_depth = Sample_name %>% str_replace("-","_")) %>% 
  select(!Sample_name)

df_dna = merge(df_dna, cn, by.x = "ROV_depth", by.y = "ROV_depth")

df_dna %>% ggplot(aes(x = ROV.x, y = d13C))+
  geom_point(aes(size = Richness_16s))+
  geom_smooth(method = "lm")+                                    #### set: se = T
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size=15, hjust = 1),
        axis.text.y = element_text(color="black", size=15))+
  theme(legend.title=element_text(size=20), legend.text=element_text(size=15),
        axis.title=element_text(size=20))+
  ggpubr::stat_cor(label.y= 400 ,method = "spearman", size = 5.3)

df_rna = merge(df_rna, cn, by.x = "ROV_depth", by.y = "ROV_depth")

df_rna %>% ggplot(aes(x = -d13C, y = Richness_18s))+
  geom_point(aes(color = ROV.x),size = 4)+
  geom_smooth(method = "lm")+                                    #### set: se = T
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size=15, hjust = 1),
        axis.text.y = element_text(color="black", size=15))+
  theme(legend.title=element_text(size=20), legend.text=element_text(size=15),
        axis.title=element_text(size=20))+
  ggpubr::stat_cor(label.y= 400 ,method = "spearman", size = 5.3)



# correlation with env factors --------------------------------------------

asv_18s_apha = alpha(asv_18s, tree_18s)
asv_16s_apha = alpha(asv_16s, tree_16s)


data_16s = readRDS("cDNA-physeq16S.RDS")
data_18s = readRDS("cDNA-physeq18S.RDS")

env_16s = data_16s@sam_data

env_16s %>% as_tibble() %>% 
  mutate(ROV_depth = Sample.name)

env_18s = data_18s@sam_data


# Mantel test between communities and env --------------------------------

env = read.csv("SEM_data_for_16S.csv")
env %>% head

# for 18s
asv_table_18s = read.table("asv_18s.txt",header = T,check.names = F) %>% t()
tax_18s = read.delim("taxonomy-pr2.tsv")
## remove metazoan
tax_pro = tax_18s[!c(str_detect(tax_18s$Taxon, "Metazoa")),]
tax_pro = tax_pro[!c(str_detect(tax_pro$Taxon, "Streptophyta")),]
asv_table_18s = asv_table_18s[,match(tax_pro$Feature.ID, colnames(asv_table_18s))]

asv18s_dna = asv_table_18s[rownames(asv_table_18s) %>% str_detect("dna"),]
asv18s_rna = asv_table_18s[rownames(asv_table_18s) %>% str_detect("rna"),]

m1 = asv18s_rna %>% rownames() %>% str_remove("_rna")
m2 = asv18s_dna %>% rownames() %>% str_remove("_dna")

asv18s_dna = asv18s_dna[match(m1,m2),]
identical(rownames(asv18s_dna) %>% str_remove("_dna"),
          rownames(asv18s_rna) %>% str_remove("_rna")) # True

asv18s_dna = asv18s_dna[,which(colSums(asv18s_dna)>0)]
asv18s_rna = asv18s_rna[,which(colSums(asv18s_rna)>0)]

## for 16S



###18S
#Breviatea
Breviatea_dna <-
  asv18s_dna %>% t() %>% data.frame() %>% 
  mutate(Tax = tax_18s$Taxon[match(colnames(asv18s_dna),tax_18s$Feature.ID)]) %>%
  filter(str_detect(Tax,"Breviatea")) %>%
  select(-Tax) %>% t() %>% data.frame(check.names = F)

Breviatea_dna = Breviatea_dna[which(rowSums(Breviatea_dna)>0),]
Breviatea_dna = Breviatea_dna[,which(colSums(Breviatea_dna)>0)]

Breviatea_dna = Breviatea_dna %>% data.frame(check.names = F) %>% 
  rownames_to_column(., var = "Sample") %>% as_tibble() %>% 
  mutate(Sample = Sample %>% str_replace("_18s_dna", ""))


df = 
left_join(env, Breviatea_dna, by = "Sample") %>% as_tibble() %>% 
  column_to_rownames(., var = "Sample") %>% 
  drop_na()

df %>% rownames()
df %>% colnames()

mantel(vegdist(df[ ,11]%>% scale, "eu"), vegdist(df[ ,-c(1:37)])) 

#16S
ANME1_dna <-
  asv16s_dna %>% t() %>% data.frame() %>% 
  mutate(Tax = tax_16s$Taxon[match(colnames(asv16s_dna),tax_16s$Feature.ID)]) %>%
  filter(str_detect(Tax,"ANME1")) %>%
  select(-Tax) %>% t() %>% data.frame(check.names = F)

ANME1_dna = ANME1_dna[which(rowSums(ANME1_dna)>0),]
ANME1_dna = ANME1_dna[,which(colSums(ANME1_dna)>0)]

ANME1_dna = ANME1_dna %>% data.frame(check.names = F) %>% 
  rownames_to_column(., var = "Sample") %>% as_tibble() %>% 
  mutate(Sample = Sample %>% str_replace("_16s_dna", ""))


df = 
  left_join(env, ANME1_dna, by = "Sample") %>% as_tibble() %>% 
  column_to_rownames(., var = "Sample") %>% 
  drop_na()

df %>% colnames()

mantel(vegdist(df[ ,11]%>% scale, "eu"), vegdist(df[ ,-c(1:37)])) 

# Relative acitivity: all lines has NA

df <-
rel_act_18s %>% data.frame() %>%
   mutate(Tax = tax_18s$Taxon[match(rownames(rel_act_18s),tax_18s$Feature.ID)]) %>%
   filter(str_detect(Tax,"Breviatea")) %>%
   select(-Tax) %>% rownames_to_column(., var = "ASV") %>% as_tibble() %>% 
  pivot_longer(!ASV, names_to = "Sample", values_to = "RA") %>% 
  na.omit() %>% 
  filter(!RA %>% is.infinite()) %>% 
  mutate(Sample = Sample %>% str_remove("_18s_rna")) %>% 
  pivot_wider(names_from = Sample, values_from = RA,values_fill = 0) %>% 
  column_to_rownames(., var = "ASV") %>% t() %>% data.frame(check.names = F) 

df = df[,which(colSums(df)>0)]
df = df[which(rowSums(df)>0),]

df<-
df %>% rownames_to_column(., var = "Sample") %>% 
  left_join(env2, by = "Sample") %>% na.omit() 

rownames(df) <- df$Sample
df = df[,-1]

df %>% colnames()

mantel(vegdist(df[ ,1:18]), vegdist(df[ ,19:20], "eu")) #r: 0.06419, Significance: 0.238
mantel(vegdist(df[ ,1:18]), vegdist(df[ ,21:23], "eu")) #  r: -0.07932 ,Significance: 0.788
mantel(vegdist(df[ ,1:18]), vegdist(df[ ,24], "eu")) #  r: -0.04649 , Significance: 0.694 
mantel(vegdist(df[ ,1:18]), vegdist(df[ ,29:33], "eu")) # r: 0.1087 ,Significance: 0.116 
mantel(vegdist(df[ ,1:18]), vegdist(df[ ,34:37], "eu"))#  r: -0.03293 ,  Significance: 0.618 
 
  
### use site average RA, liner correlation

rel_act_16s %>% sqrt() %>%  data.frame %>% 
  rownames_to_column(., var = "ASV_ID") %>%
  as_tibble() %>% 
  mutate(Taxonomy = tax_16s$Taxon[match(colnames(asv16s_ra_dna),
                                        tax_16s$Feature.ID)]) %>% 
  # filter(Taxonomy %>% str_detect("ANME-1")) %>%
  select(!Taxonomy) %>% 
  column_to_rownames(., var = "ASV_ID") -> da

da

# ### remove inf and NaN before calculate the mean value of RA at each site
#  remove_na_inf <- function(x) {
#      y = x[!is.na(x)]
#      z = y[!is.infinite(y)]
#      z = mean(z);
#      return(z)
#      }
# 
#  df8 <-
#    data.frame(RA = apply(da, 2, remove_na_inf)) %>%
#    rownames_to_column(., var = "Sample") %>%
#    mutate(Sample = Sample %>% str_remove("_16s_rna")) %>%
#    left_join(env2, by = "Sample") %>% na.omit()
# 
#  rownames(df8) <- df$Sample
# 
#  df8
#  
#  Bre_ra = df1$RA
#  Bre_13C = df1$d13C
#  
#  Api_ra = df2$RA
#  Api_13C = df2$d13C
#  
#  Lab_ra = df3$RA
#  Lab_13C = df3$d13C
#  
#  All_euk_ra = df4$RA
#  All_euk_13C = df4$d13C
#  
#  Anme_ra = df5$RA
#  Anme_13C = df5$d13C
#  
#  Met_ra = df6$RA
#  Met_13C = df6$d13C
#  
#  Gam_ra = df7$RA
#  Gam_13C = df7$d13C
#  
#  All_pro_ra = df8$RA
#  All_pro_13C = df8$d13C
#  
 
#Labyrinthulomycetes

Labyrinthulomycetes_dna <-
  asv_table_18s1 %>% t() %>% data.frame() %>% 
  mutate(Tax = tax_18s$Taxon[match(colnames(asv_table_18s1),tax_18s$Feature.ID)]) %>%
  filter(str_detect(Tax,"Labyrinthulomycetes")) %>% 
  select(-Tax) %>% t() %>% data.frame(check.names = F) %>%  
  rownames_to_column(., var = "Sample") %>% as_tibble() %>% 
  filter(Sample %>% str_detect("dna")) %>% 
  mutate(Sample = Sample %>% str_remove("_18s_dna")) %>% 
  column_to_rownames(., var = "Sample")

Labyrinthulomycetes_dna = Labyrinthulomycetes_dna[,which(colSums(Labyrinthulomycetes_dna)>0)]
Labyrinthulomycetes_dna = Labyrinthulomycetes_dna[which(rowSums(Labyrinthulomycetes_dna)>0),]
Labyrinthulomycetes_dna = Labyrinthulomycetes_dna %>% rownames_to_column(., var = "Sample")

df = 
  left_join(env2, Labyrinthulomycetes_dna, by = "Sample") %>% na.omit() %>% 
  column_to_rownames(., var = "Sample")

df %>% colnames()

mantel(vegdist(df[ ,7]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) # r: 0.1507 ,Significance: 0.021
mantel(vegdist(df[ ,3:5]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #r: -0.05113 , Significance: 0.796 
mantel(vegdist(df[ ,6]%>% scale, "eu"), vegdist(df[ ,-c(1:19)]))   #r: -0.07021  Significance: 0.94
mantel(vegdist(df[ ,11:15]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) # r: 0.09765 Significance: 0.107 
mantel(vegdist(df[ ,16:19]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) # r: 0.1467,Significance: 0.012 



Labyrinthulomycetes_rna <-
  asv_table_18s1 %>% t() %>% data.frame() %>% 
  mutate(Tax = tax_18s$Taxon[match(colnames(asv_table_18s1),tax_18s$Feature.ID)]) %>%
  filter(str_detect(Tax,"Labyrinthulomycetes")) %>% 
  select(-Tax) %>% t() %>% data.frame(check.names = F) %>%  
  rownames_to_column(., var = "Sample") %>% as_tibble() %>% 
  filter(Sample %>% str_detect("rna")) %>% 
  mutate(Sample = Sample %>% str_remove("_18s_rna")) %>% 
  column_to_rownames(., var = "Sample")

Labyrinthulomycetes_rna = Labyrinthulomycetes_rna[,which(colSums(Labyrinthulomycetes_rna)>0)]
Labyrinthulomycetes_rna = Labyrinthulomycetes_rna[which(rowSums(Labyrinthulomycetes_rna)>0),]
Labyrinthulomycetes_rna = Labyrinthulomycetes_rna %>% rownames_to_column(., var = "Sample")

df = 
  left_join(env2, Labyrinthulomycetes_rna, by = "Sample") %>% na.omit() %>% 
  column_to_rownames(., var = "Sample")

df %>% colnames()

mantel(vegdist(df[ ,8]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #  r: 0.1196 ,Significance: 0.034
mantel(vegdist(df[ ,3:5]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #  r: 0.0166 ,Significance: 0.379 
mantel(vegdist(df[ ,6]%>% scale, "eu"), vegdist(df[ ,-c(1:19)]))   #  r: 0.05855 ,Significance: 0.213 
mantel(vegdist(df[ ,11:15]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) # r: 0.1449 ,Significance: 0.039
mantel(vegdist(df[ ,16:19]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #  r: 0.1613 ,Significance: 0.009


# Apicomplexa
Apicomplexa_dna <-
  asv_table_18s1 %>% t() %>% data.frame() %>% 
  mutate(Tax = tax_18s$Taxon[match(colnames(asv_table_18s1),tax_18s$Feature.ID)]) %>%
  filter(str_detect(Tax,"Apicomplexa")) %>% 
  select(-Tax) %>% t() %>% data.frame(check.names = F) %>%  
  rownames_to_column(., var = "Sample") %>% as_tibble() %>% 
  filter(Sample %>% str_detect("dna")) %>% 
  mutate(Sample = Sample %>% str_remove("_18s_dna")) %>% 
  column_to_rownames(., var = "Sample")

Apicomplexa_dna = Apicomplexa_dna[,which(colSums(Apicomplexa_dna)>0)]
Apicomplexa_dna = Apicomplexa_dna[which(rowSums(Apicomplexa_dna)>0),]
Apicomplexa_dna = Apicomplexa_dna %>% rownames_to_column(., var = "Sample")

df = 
  left_join(env2, Apicomplexa_dna, by = "Sample") %>% na.omit() %>% 
  column_to_rownames(., var = "Sample")

df %>% colnames()

mantel(vegdist(df[ ,7]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #  r: 0.01717 ,Significance: 0.365 
mantel(vegdist(df[ ,3:5]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #   r: 0.06604 Significance: 0.056 
mantel(vegdist(df[ ,6]%>% scale, "eu"), vegdist(df[ ,-c(1:19)]))   #  r: 0.0263 Significance: 0.252
mantel(vegdist(df[ ,11:15]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) # r: -0.01736  Significance: 0.626
mantel(vegdist(df[ ,16:19]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #  r: 0.03979  Significance: 0.204 


Apicomplexa_rna <-
  asv_table_18s1 %>% t() %>% data.frame() %>% 
  mutate(Tax = tax_18s$Taxon[match(colnames(asv_table_18s1),tax_18s$Feature.ID)]) %>%
  filter(str_detect(Tax,"Apicomplexa")) %>% 
  select(-Tax) %>% t() %>% data.frame(check.names = F) %>%  
  rownames_to_column(., var = "Sample") %>% as_tibble() %>% 
  filter(Sample %>% str_detect("rna")) %>% 
  mutate(Sample = Sample %>% str_remove("_18s_rna")) %>% 
  column_to_rownames(., var = "Sample")

Apicomplexa_rna = Apicomplexa_rna[,which(colSums(Apicomplexa_rna)>0)]
Apicomplexa_rna = Apicomplexa_rna[which(rowSums(Apicomplexa_rna)>0),]
Apicomplexa_rna = Apicomplexa_rna %>% rownames_to_column(., var = "Sample")

df = 
  left_join(env2, Apicomplexa_rna, by = "Sample") %>% na.omit() %>% 
  column_to_rownames(., var = "Sample")

df %>% colnames()

mantel(vegdist(df[ ,7]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #  r: 0.07812 Significance: 0.093
mantel(vegdist(df[ ,3:5]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #   statistic r: 0.1591  Significance: 0.009 
mantel(vegdist(df[ ,6]%>% scale, "eu"), vegdist(df[ ,-c(1:19)]))   #  r: 0.1389 Significance: 0.028
mantel(vegdist(df[ ,11:15]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) # r: 0.04289 Significance: 0.292
mantel(vegdist(df[ ,16:19]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #  r: 0.02704  Significance: 0.323 

#whole 18S

Meuk_dna <-
  asv_table_18s1 %>%  data.frame(check.names = F) %>%  
  rownames_to_column(., var = "Sample") %>% as_tibble() %>% 
  filter(Sample %>% str_detect("dna")) %>% 
  mutate(Sample = Sample %>% str_remove("_18s_dna")) %>% 
  column_to_rownames(., var = "Sample")

Meuk_dna = Meuk_dna[,which(colSums(Meuk_dna)>0)]
Meuk_dna = Meuk_dna[which(rowSums(Meuk_dna)>0),]
Meuk_dna = Meuk_dna %>% rownames_to_column(., var = "Sample")

df = 
  left_join(env2, Meuk_dna, by = "Sample") %>% na.omit() %>% 
  column_to_rownames(., var = "Sample")

df %>% colnames()


mantel(vegdist(df[ ,7]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #  r: 0.1269   Significance: 0.001 
mantel(vegdist(df[ ,3:5]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #  r: 0.05448  Significance: 0.051
mantel(vegdist(df[ ,6]%>% scale, "eu"), vegdist(df[ ,-c(1:19)]))   #   r: 0.05988 Significance: 0.037 
mantel(vegdist(df[ ,11:15]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #  r: 0.102  Significance: 0.003
mantel(vegdist(df[ ,16:19]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #   r: 0.2167 Significance: 0.001


Meuk_rna <-
  asv_table_18s1 %>%  data.frame(check.names = F) %>%  
  rownames_to_column(., var = "Sample") %>% as_tibble() %>% 
  filter(Sample %>% str_detect("rna")) %>% 
  mutate(Sample = Sample %>% str_remove("_18s_rna")) %>% 
  column_to_rownames(., var = "Sample")

Meuk_rna = Meuk_rna[,which(colSums(Meuk_rna)>0)]
Meuk_rna = Meuk_rna[which(rowSums(Meuk_rna)>0),]
Meuk_rna = Meuk_rna %>% rownames_to_column(., var = "Sample")

df = 
  left_join(env2, Meuk_rna, by = "Sample") %>% na.omit() %>% 
  column_to_rownames(., var = "Sample")

df %>% colnames()


mantel(vegdist(df[ ,7]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #   r: 0.1721 Significance: 0.001
mantel(vegdist(df[ ,3:5]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #   r: 0.07748 Significance: 0.038
mantel(vegdist(df[ ,6]%>% scale, "eu"), vegdist(df[ ,-c(1:19)]))   #    r: 0.07844 Significance: 0.037 
mantel(vegdist(df[ ,11:15]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #   r: 0.1264 Significance: 0.006
mantel(vegdist(df[ ,16:19]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #    r: 0.2109  Significance: 0.001

### 16S
#ANME-1
ANME1_dna <-
  asv_table_16s1 %>% t() %>% data.frame() %>% 
  mutate(Tax = tax_16s$Taxon[match(colnames(asv_table_16s1),tax_16s$Feature.ID)]) %>%
  filter(str_detect(Tax,"ANME-1")) %>% 
  select(-Tax) %>% t() %>% data.frame(check.names = F) %>%  
  rownames_to_column(., var = "Sample") %>% as_tibble() %>% 
  filter(Sample %>% str_detect("dna")) %>% 
  mutate(Sample = Sample %>% str_remove("_16s_dna")) %>% 
  column_to_rownames(., var = "Sample")

ANME1_dna = ANME1_dna[,which(colSums(ANME1_dna)>0)]
ANME1_dna = ANME1_dna[which(rowSums(ANME1_dna)>0),]
ANME1_dna = ANME1_dna %>% rownames_to_column(., var = "Sample")

df = 
  left_join(env2, ANME1_dna, by = "Sample") %>% na.omit() %>% 
  column_to_rownames(., var = "Sample")

df %>% colnames()


mantel(vegdist(df[ ,7]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #  r: -0.04067  Significance: 0.916
mantel(vegdist(df[ ,3:5]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #   statistic r: 0.05722 Significance: 0.054 
mantel(vegdist(df[ ,6]%>% scale, "eu"), vegdist(df[ ,-c(1:19)]))   #   r: 0.0643  Significance: 0.027
mantel(vegdist(df[ ,11:15]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #  r: -0.04983  Significance: 0.977 
mantel(vegdist(df[ ,16:19]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #  r: 0.02613  Significance: 0.245 


ANME1_rna <-
  asv_table_16s1 %>% t() %>% data.frame() %>% 
  mutate(Tax = tax_16s$Taxon[match(colnames(asv_table_16s1),tax_16s$Feature.ID)]) %>%
  filter(str_detect(Tax,"ANME-1")) %>% 
  select(-Tax) %>% t() %>% data.frame(check.names = F) %>%  
  rownames_to_column(., var = "Sample") %>% as_tibble() %>% 
  filter(Sample %>% str_detect("rna")) %>% 
  mutate(Sample = Sample %>% str_remove("_16s_rna")) %>% 
  column_to_rownames(., var = "Sample")

ANME1_rna = ANME1_rna[,which(colSums(ANME1_rna)>0)]
ANME1_rna = ANME1_rna[which(rowSums(ANME1_rna)>0),]
ANME1_rna = ANME1_rna %>% rownames_to_column(., var = "Sample")

df = 
  left_join(env2, ANME1_rna, by = "Sample") %>% na.omit() %>% 
  column_to_rownames(., var = "Sample")

df %>% colnames()


mantel(vegdist(df[ ,7]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #  r: -0.03328 Significance: 0.777 
mantel(vegdist(df[ ,3:5]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #    r: 0.08641  Significance: 0.019 
mantel(vegdist(df[ ,6]%>% scale, "eu"), vegdist(df[ ,-c(1:19)]))   #   r: 0.1265  Significance: 0.001 
mantel(vegdist(df[ ,11:15]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #    r: -0.05259  Significance: 0.891
mantel(vegdist(df[ ,16:19]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #  r: 0.09299 Significance: 0.032 

# Methanosarcinia

Methanosarcinia_dna <-
  asv_table_16s1 %>% t() %>% data.frame() %>% 
  mutate(Tax = tax_16s$Taxon[match(colnames(asv_table_16s1),tax_16s$Feature.ID)]) %>%
  filter(str_detect(Tax,"Methanosarcinia")) %>% 
  select(-Tax) %>% t() %>% data.frame(check.names = F) %>%  
  rownames_to_column(., var = "Sample") %>% as_tibble() %>% 
  filter(Sample %>% str_detect("dna")) %>% 
  mutate(Sample = Sample %>% str_remove("_16s_dna")) %>% 
  column_to_rownames(., var = "Sample")

Methanosarcinia_dna = Methanosarcinia_dna[,which(colSums(Methanosarcinia_dna)>0)]
Methanosarcinia_dna = Methanosarcinia_dna[which(rowSums(Methanosarcinia_dna)>0),]
Methanosarcinia_dna = Methanosarcinia_dna %>% rownames_to_column(., var = "Sample")

df = 
  left_join(env2, Methanosarcinia_dna, by = "Sample") %>% na.omit() %>% 
  column_to_rownames(., var = "Sample")

df %>% colnames()


mantel(vegdist(df[ ,7]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #  r: -0.03209 Significance: 0.765 
mantel(vegdist(df[ ,3:5]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #   r: 0.02908  Significance: 0.25 
mantel(vegdist(df[ ,6]%>% scale, "eu"), vegdist(df[ ,-c(1:19)]))   #   r: 0.02128 Significance: 0.283 
mantel(vegdist(df[ ,11:15]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #   r: -0.07654 Significance: 0.895
mantel(vegdist(df[ ,16:19]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #   r: 0.01964  Significance: 0.344


Methanosarcinia_rna <-
  asv_table_16s1 %>% t() %>% data.frame() %>% 
  mutate(Tax = tax_16s$Taxon[match(colnames(asv_table_16s1),tax_16s$Feature.ID)]) %>%
  filter(str_detect(Tax,"Methanosarcinia")) %>% 
  select(-Tax) %>% t() %>% data.frame(check.names = F) %>%  
  rownames_to_column(., var = "Sample") %>% as_tibble() %>% 
  filter(Sample %>% str_detect("rna")) %>% 
  mutate(Sample = Sample %>% str_remove("_16s_rna")) %>% 
  column_to_rownames(., var = "Sample")

Methanosarcinia_rna = Methanosarcinia_rna[,which(colSums(Methanosarcinia_rna)>0)]
Methanosarcinia_rna = Methanosarcinia_rna[which(rowSums(Methanosarcinia_rna)>0),]
Methanosarcinia_rna = Methanosarcinia_rna %>% rownames_to_column(., var = "Sample")

df = 
  left_join(env2, Methanosarcinia_rna, by = "Sample") %>% na.omit() %>% 
  column_to_rownames(., var = "Sample")

df %>% colnames()


mantel(vegdist(df[ ,8]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #   r: 0.1626 Significance: 0.007 
mantel(vegdist(df[ ,3:5]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #  r: 0.01378 Significance: 0.359 
mantel(vegdist(df[ ,6]%>% scale, "eu"), vegdist(df[ ,-c(1:19)]))   #    r: -0.0223 Significance: 0.63 
mantel(vegdist(df[ ,11:15]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #  r: 0.06142  Significance: 0.204 
mantel(vegdist(df[ ,16:19]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #    r: 0.09305  Significance: 0.064

#Gammaproteobacteria
Gammaproteobacteria_dna <-
  asv_table_16s1 %>% t() %>% data.frame() %>% 
  mutate(Tax = tax_16s$Taxon[match(colnames(asv_table_16s1),tax_16s$Feature.ID)]) %>%
  filter(str_detect(Tax,"Gammaproteobacteria")) %>% 
  select(-Tax) %>% t() %>% data.frame(check.names = F) %>%  
  rownames_to_column(., var = "Sample") %>% as_tibble() %>% 
  filter(Sample %>% str_detect("dna")) %>% 
  mutate(Sample = Sample %>% str_remove("_16s_dna")) %>% 
  column_to_rownames(., var = "Sample")

Gammaproteobacteria_dna = Gammaproteobacteria_dna[,which(colSums(Gammaproteobacteria_dna)>0)]
Gammaproteobacteria_dna = Gammaproteobacteria_dna[which(rowSums(Gammaproteobacteria_dna)>0),]
Gammaproteobacteria_dna = Gammaproteobacteria_dna %>% rownames_to_column(., var = "Sample")

df = 
  left_join(env2, Gammaproteobacteria_dna, by = "Sample") %>% na.omit() %>% 
  column_to_rownames(., var = "Sample")

df %>% colnames()


mantel(vegdist(df[ ,7]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #   r: 0.1681  Significance: 0.004
mantel(vegdist(df[ ,3:5]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #     r: 0.03974 Significance: 0.236 
mantel(vegdist(df[ ,6]%>% scale, "eu"), vegdist(df[ ,-c(1:19)]))   #    r: 0.002996 Significance: 0.433 
mantel(vegdist(df[ ,11:15]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #  r: 0.02776  Significance: 0.303
mantel(vegdist(df[ ,16:19]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #  r: 0.2475  Significance: 0.001 

Gammaproteobacteria_rna <-
  asv_table_16s1 %>% t() %>% data.frame() %>% 
  mutate(Tax = tax_16s$Taxon[match(colnames(asv_table_16s1),tax_16s$Feature.ID)]) %>%
  filter(str_detect(Tax,"Gammaproteobacteria")) %>% 
  select(-Tax) %>% t() %>% data.frame(check.names = F) %>%  
  rownames_to_column(., var = "Sample") %>% as_tibble() %>% 
  filter(Sample %>% str_detect("rna")) %>% 
  mutate(Sample = Sample %>% str_remove("_16s_rna")) %>% 
  column_to_rownames(., var = "Sample")

Gammaproteobacteria_rna = Gammaproteobacteria_rna[,which(colSums(Gammaproteobacteria_rna)>0)]
Gammaproteobacteria_rna = Gammaproteobacteria_rna[which(rowSums(Gammaproteobacteria_rna)>0),]
Gammaproteobacteria_rna = Gammaproteobacteria_rna %>% rownames_to_column(., var = "Sample")

df = 
  left_join(env2, Gammaproteobacteria_rna, by = "Sample") %>% na.omit() %>% 
  column_to_rownames(., var = "Sample")

df %>% colnames()


mantel(vegdist(df[ ,7]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #   r: 0.3099  Significance: 0.001
mantel(vegdist(df[ ,3:5]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #      r: 0.209 Significance: 0.001
mantel(vegdist(df[ ,6]%>% scale, "eu"), vegdist(df[ ,-c(1:19)]))   #     r: 0.2558  Significance: 0.001
mantel(vegdist(df[ ,11:15]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #   r: 0.06065 Significance: 0.179
mantel(vegdist(df[ ,16:19]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #   r: 0.3068 Significance: 0.001



#whole 16s
Prok_dna <-
  asv_table_16s1 %>%  data.frame(check.names = F) %>%  
  rownames_to_column(., var = "Sample") %>% as_tibble() %>% 
  filter(Sample %>% str_detect("dna")) %>% 
  mutate(Sample = Sample %>% str_remove("_16s_dna")) %>% 
  column_to_rownames(., var = "Sample")

Prok_dna = Prok_dna[,which(colSums(Prok_dna)>0)]
Prok_dna = Prok_dna[which(rowSums(Prok_dna)>0),]
Prok_dna = Prok_dna %>% rownames_to_column(., var = "Sample")

df = 
  left_join(env2, Prok_dna, by = "Sample") %>% na.omit() %>% 
  column_to_rownames(., var = "Sample")

colnames(df)

mantel(vegdist(df[ ,7] %>% scale ,"eu"), vegdist(df[ ,-c(1:19)])) #r: 0.02747  Significance: 0.273 
mantel(vegdist(df[ ,3:5]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #     r: 0.03377   Significance: 0.228 
mantel(vegdist(df[ ,6]%>% scale, "eu"), vegdist(df[ ,-c(1:19)]))   #    -0.009   Significance: 0.527 
mantel(vegdist(df[ ,11:15]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #    r: -0.01816  Significance: 0.533
mantel(vegdist(df[ ,16:19]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #    r: 0.1728 Significance: 0.004 

Prok_rna <-
  asv_table_16s1 %>%  data.frame(check.names = F) %>%  
  rownames_to_column(., var = "Sample") %>% as_tibble() %>% 
  filter(Sample %>% str_detect("rna")) %>% 
  mutate(Sample = Sample %>% str_remove("_16s_rna")) %>% 
  column_to_rownames(., var = "Sample")

Prok_rna = Prok_rna[,which(colSums(Prok_rna)>0)]
Prok_rna = Prok_rna[which(rowSums(Prok_rna)>0),]
Prok_rna = Prok_rna %>% rownames_to_column(., var = "Sample")

df = 
  left_join(env2, Prok_rna, by = "Sample") %>% na.omit() %>% 
  column_to_rownames(., var = "Sample")

colnames(df)

mantel(vegdist(df[ ,7] %>% scale ,"eu"), vegdist(df[ ,-c(1:19)])) #
mantel(vegdist(df[ ,3:5]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #     
mantel(vegdist(df[ ,6]%>% scale, "eu"), vegdist(df[ ,-c(1:19)]))   #   
mantel(vegdist(df[ ,11:15]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #  
mantel(vegdist(df[ ,16:19]%>% scale, "eu"), vegdist(df[ ,-c(1:19)])) #   



# Niche breadth -----------------------------------------------------------
library(spaa)
asv18s_dna_rov1.bcom = t(niche.width(asv18s_dna_rov1,method = "levins"))[,1]
asv18s_dna_rov2.bcom = t(niche.width(asv18s_dna_rov2,method = "levins"))[,1]
asv18s_dna_rov3.bcom = t(niche.width(asv18s_dna_rov3,method = "levins"))[,1]
asv18s_dna_rov4.bcom = t(niche.width(asv18s_dna_rov4,method = "levins"))[,1]
asv18s_dna_rov5.bcom = t(niche.width(asv18s_dna_rov5,method = "levins"))[,1]


asv18s_rna_rov1.bcom = t(niche.width(asv18s_rna_rov1,method = "levins"))[,1]
asv18s_rna_rov2.bcom = t(niche.width(asv18s_rna_rov2,method = "levins"))[,1]
asv18s_rna_rov3.bcom = t(niche.width(asv18s_rna_rov3,method = "levins"))[,1]
asv18s_rna_rov4.bcom = t(niche.width(asv18s_rna_rov4,method = "levins"))[,1]
asv18s_rna_rov5.bcom = t(niche.width(asv18s_rna_rov5,method = "levins"))[,1]


asv16s_dna_rov1.bcom = t(niche.width(asv16s_dna_rov1,method = "levins"))[,1]
asv16s_dna_rov2.bcom = t(niche.width(asv16s_dna_rov2,method = "levins"))[,1]
asv16s_dna_rov3.bcom = t(niche.width(asv16s_dna_rov3,method = "levins"))[,1]
asv16s_dna_rov4.bcom = t(niche.width(asv16s_dna_rov4,method = "levins"))[,1]
asv16s_dna_rov5.bcom = t(niche.width(asv16s_dna_rov5,method = "levins"))[,1]


asv16s_rna_rov1.bcom = t(niche.width(asv16s_rna_rov1,method = "levins"))[,1]
asv16s_rna_rov2.bcom = t(niche.width(asv16s_rna_rov2,method = "levins"))[,1]
asv16s_rna_rov3.bcom = t(niche.width(asv16s_rna_rov3,method = "levins"))[,1]
asv16s_rna_rov4.bcom = t(niche.width(asv16s_rna_rov4,method = "levins"))[,1]
asv16s_rna_rov5.bcom = t(niche.width(asv16s_rna_rov5,method = "levins"))[,1]


p1<-
data.frame(bcom = c(asv18s_dna_rov1.bcom, 
                    asv18s_dna_rov2.bcom,
                    asv18s_dna_rov3.bcom,
                    asv18s_dna_rov4.bcom,
                    asv18s_dna_rov5.bcom,
                    asv18s_rna_rov1.bcom, 
                    asv18s_rna_rov2.bcom,
                    asv18s_rna_rov3.bcom,
                    asv18s_rna_rov4.bcom,
                    asv18s_rna_rov5.bcom),
           grp = c(rep("DNA_ROV1", asv18s_dna_rov1.bcom %>% length()),
                   rep("DNA_ROV2", asv18s_dna_rov2.bcom %>% length()),
                   rep("DNA_ROV3", asv18s_dna_rov3.bcom %>% length()),
                   rep("DNA_ROV4", asv18s_dna_rov4.bcom %>% length()),
                   rep("DNA_ROV5", asv18s_dna_rov5.bcom %>% length()),
                   rep("RNA_ROV1", asv18s_rna_rov1.bcom %>% length()),
                   rep("RNA_ROV2", asv18s_rna_rov2.bcom %>% length()),
                   rep("RNA_ROV3", asv18s_rna_rov3.bcom %>% length()),
                   rep("RNA_ROV4", asv18s_rna_rov4.bcom %>% length()),
                   rep("RNA_ROV5", asv18s_rna_rov5.bcom %>% length()))) %>% 
  mutate(ROV = grp %>% str_split("_") %>% sapply('[',2)) %>% 
  ggplot(aes(x = grp, y = bcom))+
  #geom_violin(trim=T,fill="white", linewidth = 1)+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=15),
        axis.text.y = element_text(color="black", size=15))+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), 
             size = 2,aes(color = ROV), alpha = 0.1,show.legend = T)+
  geom_boxplot(color = "grey33", linewidth = 1)+
  theme(legend.title=element_text(size=15), legend.text=element_text(size=15),
        axis.title=element_text(size=15))+
  scale_color_manual(values = c("firebrick1", "hotpink", "darkgoldenrod1", "cornflowerblue","darkcyan"))+
  ylab("Niche breadth (Levins)")+
  ggtitle("Microeukaryotes")

p1

p2<-
data.frame(bcom = c(asv16s_dna_rov1.bcom, 
                    asv16s_dna_rov2.bcom,
                    asv16s_dna_rov3.bcom,
                    asv16s_dna_rov4.bcom,
                    asv16s_dna_rov5.bcom,
                    asv16s_rna_rov1.bcom, 
                    asv16s_rna_rov2.bcom,
                    asv16s_rna_rov3.bcom,
                    asv16s_rna_rov4.bcom,
                    asv16s_rna_rov5.bcom),
           grp = c(rep("DNA_ROV1", asv16s_dna_rov1.bcom %>% length()),
                   rep("DNA_ROV2", asv16s_dna_rov2.bcom %>% length()),
                   rep("DNA_ROV3", asv16s_dna_rov3.bcom %>% length()),
                   rep("DNA_ROV4", asv16s_dna_rov4.bcom %>% length()),
                   rep("DNA_ROV5", asv16s_dna_rov5.bcom %>% length()),
                   rep("RNA_ROV1", asv16s_rna_rov1.bcom %>% length()),
                   rep("RNA_ROV2", asv16s_rna_rov2.bcom %>% length()),
                   rep("RNA_ROV3", asv16s_rna_rov3.bcom %>% length()),
                   rep("RNA_ROV4", asv16s_rna_rov4.bcom %>% length()),
                   rep("RNA_ROV5", asv16s_rna_rov5.bcom %>% length()))) %>% 
  mutate(ROV = grp %>% str_split("_") %>% sapply('[',2)) %>% 
  ggplot(aes(x = grp, y = bcom))+
  # geom_violin(trim=F,fill="white", linewidth = 0.5)+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=15),
        axis.text.y = element_text(color="black", size=15))+
  geom_point(position=position_jitterdodge(jitter.width=0.5, dodge.width = 0), size = 2,aes(color = ROV), alpha = 0.5,show.legend = T)+
  geom_boxplot(color = "grey33", linewidth = 1)+
  theme(legend.title=element_text(size=15), legend.text=element_text(size=15),
        axis.title=element_text(size=15))+
  scale_color_manual(values = c("firebrick1", "hotpink", "darkgoldenrod1", "cornflowerblue","darkcyan"))+
  ylab("Niche breadth (Levins)")+
  ggtitle("Prokaryotes")

ggarrange(p1,p2, nrow = 2, labels = c("(a)", "(b)"))


# data for making SEM -----------------------------------------------------

## 18S (only consider alpha diversity, exclude beta)
#  18S whole community
asv_18s = read.table("asv_18s.txt",header = T,check.names = F) %>% t()
tax_18s = read.delim("taxonomy-pr2.tsv")


## remove metazoan and land plants
tax_pro = tax_18s[!c(str_detect(tax_18s$Taxon, "Metazoa")),]
tax_pro = tax_pro[!c(str_detect(tax_pro$Taxon, "Streptophyta")),]
asv_18s = asv_18s[,match(tax_pro$Feature.ID, colnames(asv_18s))]
asv_18s = rrarefy(asv_18s,min(rowSums(asv_18s)))

asv18s_dna = asv_18s[asv_18s %>% rownames() %>% str_detect("dna"), ]
asv18s_dna %>% rownames()
asv18s_dna = asv18s_dna[,which(colSums(asv18s_dna)>0)]
asv18s_dna %>% dim
asv18s_dna = rrarefy(asv18s_dna, min(rowSums(asv18s_dna)))
asv18s_dna = asv18s_dna[,which(colSums(asv18s_dna)>0)]

asv18s_rna = asv_18s[asv_18s %>% rownames() %>% str_detect("rna"), ]
asv18s_rna %>% rownames()
asv18s_rna = asv18s_rna[,which(colSums(asv18s_rna)>0)]
asv18s_rna %>% dim
asv18s_rna = rrarefy(asv18s_rna, min(rowSums(asv18s_rna)))
asv18s_rna = asv18s_rna[,which(colSums(asv18s_rna)>0)]

# 18S whole
asv18s_dna_alpha = alpha(asv18s_dna)

asv18s_dna_richness =                                   
  asv18s_dna_alpha %>% rownames_to_column(., var = "Sample_name") %>% 
  mutate(Sample = Sample_name %>% str_remove("_18s_dna")) %>% 
  select(Sample, Richness) %>% 
  rename(DNA_18S_Whole_Richness = Richness)

asv18s_rna_alpha = alpha(asv18s_rna)
asv18s_rna_richness = 
  asv18s_rna_alpha %>% rownames_to_column(., var = "Sample_name") %>% 
  mutate(Sample = Sample_name %>% str_remove("_18s_rna")) %>% 
  select(Sample, Richness) %>% 
  rename(RNA_18S_Whole_Richness = Richness)

#18S active groups
#Breviatea
Breviatea_dna_richness <- 
asv18s_dna %>% t() %>% data.frame() %>% rownames_to_column(., var = "ASV") %>% 
  as_tibble() %>% 
  mutate(Tax = tax_18s$Taxon[match(colnames(asv18s_dna),tax_18s$Feature.ID)]) %>%
  filter(str_detect(Tax,"Breviatea")) %>% 
  select(-Tax) %>% column_to_rownames(., var = "ASV") %>% 
  t() %>% alpha() %>% 
  rownames_to_column(., var = "Sample_name") %>%
  mutate(Sample = Sample_name %>% str_remove("_18s_dna")) %>% 
  rename(DNA_Breviatea_Richness = Richness) %>% 
  select(Sample, DNA_Breviatea_Richness)

Breviatea_rna_richness <- 
asv18s_rna %>% t() %>% data.frame() %>% rownames_to_column(., var = "ASV") %>% 
  as_tibble() %>% 
  mutate(Tax = tax_18s$Taxon[match(colnames(asv18s_rna),tax_18s$Feature.ID)]) %>%
  filter(str_detect(Tax,"Breviatea")) %>% 
  select(-Tax) %>% column_to_rownames(., var = "ASV") %>% 
  t() %>% alpha() %>% 
  rownames_to_column(., var = "Sample_name") %>%
  mutate(Sample = Sample_name %>% str_remove("_18s_rna")) %>% 
  rename(RNA_Breviatea_Richness = Richness) %>% 
  select(Sample, RNA_Breviatea_Richness)

# Apicomplexa
Apicomplexa_dna_richness <- 
  asv18s_dna %>% t() %>% data.frame() %>% rownames_to_column(., var = "ASV") %>% 
  as_tibble() %>% 
  mutate(Tax = tax_18s$Taxon[match(colnames(asv18s_dna),tax_18s$Feature.ID)]) %>%
  filter(str_detect(Tax,"Apicomplexa")) %>% 
  select(-Tax) %>% column_to_rownames(., var = "ASV") %>% 
  t() %>% alpha() %>% 
  rownames_to_column(., var = "Sample_name") %>%
  mutate(Sample = Sample_name %>% str_remove("_18s_dna")) %>% 
  rename(DNA_Apicomplexa_Richness = Richness) %>% 
  select(Sample, DNA_Apicomplexa_Richness)

Apicomplexa_rna_richness <- 
  asv18s_rna %>% t() %>% data.frame() %>% rownames_to_column(., var = "ASV") %>% 
  as_tibble() %>% 
  mutate(Tax = tax_18s$Taxon[match(colnames(asv18s_rna),tax_18s$Feature.ID)]) %>%
  filter(str_detect(Tax,"Apicomplexa")) %>% 
  select(-Tax) %>% column_to_rownames(., var = "ASV") %>% 
  t() %>% alpha() %>% 
  rownames_to_column(., var = "Sample_name") %>%
  mutate(Sample = Sample_name %>% str_remove("_18s_rna")) %>% 
  rename(RNA_Apicomplexa_Richness = Richness) %>% 
  select(Sample, RNA_Apicomplexa_Richness)

#Labyrinthulomycetes
Laby_dna_richness <- 
  asv18s_dna %>% t() %>% data.frame() %>% rownames_to_column(., var = "ASV") %>% 
  as_tibble() %>% 
  mutate(Tax = tax_18s$Taxon[match(colnames(asv18s_dna),tax_18s$Feature.ID)]) %>%
  filter(str_detect(Tax,"Laby")) %>% 
  select(-Tax) %>% column_to_rownames(., var = "ASV") %>% 
  t() %>% alpha() %>% 
  rownames_to_column(., var = "Sample_name") %>%
  mutate(Sample = Sample_name %>% str_remove("_18s_dna")) %>% 
  rename(DNA_Laby_Richness = Richness) %>% 
  select(Sample, DNA_Laby_Richness)

Laby_rna_richness <- 
  asv18s_rna %>% t() %>% data.frame() %>% rownames_to_column(., var = "ASV") %>% 
  as_tibble() %>% 
  mutate(Tax = tax_18s$Taxon[match(colnames(asv18s_rna),tax_18s$Feature.ID)]) %>%
  filter(str_detect(Tax,"Laby")) %>% 
  select(-Tax) %>% column_to_rownames(., var = "ASV") %>% 
  t() %>% alpha() %>% 
  rownames_to_column(., var = "Sample_name") %>%
  mutate(Sample = Sample_name %>% str_remove("_18s_rna")) %>% 
  rename(RNA_Laby_Richness = Richness) %>% 
  select(Sample, RNA_Laby_Richness)

SEM_data_for_18S = env2 %>% 
  filter(Sample %>% str_detect("_0|_5|_10")) %>% 
  left_join(asv18s_dna_richness) %>% 
  left_join(asv18s_rna_richness) %>% 
  left_join(Breviatea_dna_richness) %>% 
  left_join(Breviatea_rna_richness) %>% 
  left_join(Laby_dna_richness) %>% 
  left_join(Laby_rna_richness) %>% 
  left_join(Apicomplexa_dna_richness) %>% 
  left_join(Apicomplexa_rna_richness)
  
write.csv(SEM_data_for_18S,"SEM_data_for_18S.csv")

## 16S
asv_table_16s = read.table("asv_table_16s.txt",header = T,check.names = F) %>% t()
asv_table_16s %>% dim
tax_16s = read.delim("taxonomy-16s.tsv") # 108166 taxon
tax_16s = tax_16s[!c(str_detect(tax_16s$Taxon, "Chloroplast")),] # 108061 taxon remained
tax_16s = tax_16s[!c(str_detect(tax_16s$Taxon, "Mitochondria")),] #107915 taxon remained

asv_table_16s = asv_table_16s[,match(tax_16s$Feature.ID, colnames(asv_table_16s))]

asv16s_dna = asv_table_16s[asv_table_16s %>% rownames() %>% str_detect("dna"), ]
asv16s_rna = asv_table_16s[asv_table_16s %>% rownames() %>% str_detect("rna"), ]

asv16s_dna = rrarefy(asv16s_dna, min(rowSums(asv16s_dna)))
asv16s_dna = asv16s_dna[,which(colSums(asv16s_dna)>0)]
asv16s_rna = rrarefy(asv16s_rna, min(rowSums(asv16s_rna)))
asv16s_rna = asv16s_rna[,which(colSums(asv16s_rna)>0)]

# whole 16s
asv16s_dna_alpha = alpha(asv16s_dna)

asv16s_dna_richness =                                   
  asv16s_dna_alpha %>% rownames_to_column(., var = "Sample_name") %>% 
  mutate(Sample = Sample_name %>% str_remove("_16s_dna")) %>% 
  select(Sample, Richness) %>% 
  rename(DNA_16s_Whole_Richness = Richness)

asv16s_rna_alpha = alpha(asv16s_rna)
asv16s_rna_richness = 
  asv16s_rna_alpha %>% rownames_to_column(., var = "Sample_name") %>% 
  mutate(Sample = Sample_name %>% str_remove("_16s_rna")) %>% 
  select(Sample, Richness) %>% 
  rename(RNA_16s_Whole_Richness = Richness)

#16S active groups
#ANME-1
ANME1_dna_richness <- 
  asv16s_dna %>% t() %>% data.frame() %>% rownames_to_column(., var = "ASV") %>% 
  as_tibble() %>% 
  mutate(Tax = tax_16s$Taxon[match(colnames(asv16s_dna),tax_16s$Feature.ID)]) %>%
  filter(str_detect(Tax,"ANME-1")) %>% 
  select(-Tax) %>% column_to_rownames(., var = "ASV") %>% 
  t() %>% alpha() %>% 
  rownames_to_column(., var = "Sample_name") %>%
  mutate(Sample = Sample_name %>% str_remove("_16s_dna")) %>% 
  rename(DNA_ANME1_Richness = Richness) %>% 
  select(Sample, DNA_ANME1_Richness)

ANME1_rna_richness <- 
  asv16s_rna %>% t() %>% data.frame() %>% rownames_to_column(., var = "ASV") %>% 
  as_tibble() %>% 
  mutate(Tax = tax_16s$Taxon[match(colnames(asv16s_rna),tax_16s$Feature.ID)]) %>%
  filter(str_detect(Tax,"ANME-1")) %>% 
  select(-Tax) %>% column_to_rownames(., var = "ASV") %>% 
  t() %>% alpha() %>% 
  rownames_to_column(., var = "Sample_name") %>%
  mutate(Sample = Sample_name %>% str_remove("_16s_rna")) %>% 
  rename(RNA_ANME1_Richness = Richness) %>% 
  select(Sample, RNA_ANME1_Richness)

# Methanosarcinia
Methano_dna_richness <- 
  asv16s_dna %>% t() %>% data.frame() %>% rownames_to_column(., var = "ASV") %>% 
  as_tibble() %>% 
  mutate(Tax = tax_16s$Taxon[match(colnames(asv16s_dna),tax_16s$Feature.ID)]) %>%
  filter(str_detect(Tax,"Methanosarcinia")) %>% 
  select(-Tax) %>% column_to_rownames(., var = "ASV") %>% 
  t() %>% alpha() %>% 
  rownames_to_column(., var = "Sample_name") %>%
  mutate(Sample = Sample_name %>% str_remove("_16s_dna")) %>% 
  rename(DNA_Methano_Richness = Richness) %>% 
  select(Sample, DNA_Methano_Richness)

Methano_rna_richness <- 
  asv16s_rna %>% t() %>% data.frame() %>% rownames_to_column(., var = "ASV") %>% 
  as_tibble() %>% 
  mutate(Tax = tax_16s$Taxon[match(colnames(asv16s_rna),tax_16s$Feature.ID)]) %>%
  filter(str_detect(Tax,"Methanosarcinia")) %>% 
  select(-Tax) %>% column_to_rownames(., var = "ASV") %>% 
  t() %>% alpha() %>% 
  rownames_to_column(., var = "Sample_name") %>%
  mutate(Sample = Sample_name %>% str_remove("_16s_rna")) %>% 
  rename(RNA_Methano_Richness = Richness) %>% 
  select(Sample, RNA_Methano_Richness)

#Gammaproteobacteria
Gamma_dna_richness <- 
  asv16s_dna %>% t() %>% data.frame() %>% rownames_to_column(., var = "ASV") %>% 
  as_tibble() %>% 
  mutate(Tax = tax_16s$Taxon[match(colnames(asv16s_dna),tax_16s$Feature.ID)]) %>%
  filter(str_detect(Tax,"Gammaproteobacteria")) %>% 
  select(-Tax) %>% column_to_rownames(., var = "ASV") %>% 
  t() %>% alpha() %>% 
  rownames_to_column(., var = "Sample_name") %>%
  mutate(Sample = Sample_name %>% str_remove("_16s_dna")) %>% 
  rename(DNA_Gamma_Richness = Richness) %>% 
  select(Sample, DNA_Gamma_Richness)

Gamma_rna_richness <- 
  asv16s_rna %>% t() %>% data.frame() %>% rownames_to_column(., var = "ASV") %>% 
  as_tibble() %>% 
  mutate(Tax = tax_16s$Taxon[match(colnames(asv16s_rna),tax_16s$Feature.ID)]) %>%
  filter(str_detect(Tax,"Gammaproteobacteria")) %>% 
  select(-Tax) %>% column_to_rownames(., var = "ASV") %>% 
  t() %>% alpha() %>% 
  rownames_to_column(., var = "Sample_name") %>%
  mutate(Sample = Sample_name %>% str_remove("_16s_rna")) %>% 
  rename(RNA_Gamma_Richness = Richness) %>% 
  select(Sample, RNA_Gamma_Richness)

SEM_data_for_16S = env2 %>% 
  filter(Sample %>% str_detect("_0|_5|_10")) %>% 
  left_join(asv16s_dna_richness) %>% 
  left_join(asv16s_rna_richness) %>% 
  left_join(Gamma_dna_richness) %>% 
  left_join(Gamma_rna_richness) %>% 
  left_join(ANME1_dna_richness) %>% 
  left_join(ANME1_rna_richness) %>% 
  left_join(Methano_dna_richness) %>% 
  left_join(Methano_rna_richness)


write.csv(SEM_data_for_16S,"SEM_data_for_16S.csv")


# plot of the correlations between richness and and env ---------------

library(corrplot)

M = cor(SEM_data_for_18S %>% na.omit() %>%  column_to_rownames(., var = "Sample"))
testRes = cor.mtest(SEM_data_for_18S %>% na.omit() %>% 
                      column_to_rownames(., var = "Sample"), 
                    conf.level = 0.95)

corrplot(-M, p.mat = testRes$p, method = 'circle',
         diag=FALSE,tl.col = "black",
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.8,
         insig = 'label_sig', pch.col = 'grey30')


M = cor(SEM_data_for_16S %>% na.omit() %>%  column_to_rownames(., var = "Sample"))
testRes = cor.mtest(SEM_data_for_16S %>% na.omit() %>% 
                      column_to_rownames(., var = "Sample"), 
                    conf.level = 0.95)

corrplot(-M, p.mat = testRes$p, method = 'circle',
         diag=FALSE,tl.col = "black",
         sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.8,
         insig = 'label_sig', pch.col = 'grey30')



# correlation between alpha diversity and env -----------------------------

df_18s = read.csv("SEM_data_for_18S.csv")

df_18s %>% as_tibble() %>% 
  mutate(ROV = Sample %>% str_sub(5,8)) %>% 
  ggplot(aes(x = d13C, y = RNA_Breviatea_Richness))+
  geom_point(aes(color = ROV),size = 4)+
  geom_smooth(method = "lm")+                                    #### set: se = T
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size=20),
        axis.text.y = element_text(color="black", size=20))+
  theme(legend.title=element_text(size=20), legend.text=element_text(size=20),
        axis.title=element_text(size=20))+
  ggpubr::stat_cor(label.y= 550 ,method = "pearson", size = 8)+
  xlab("%OC")+
  ylab("18S-RNA Richness")+
  theme(legend.position="none")

df_18s %>% as_tibble() %>% 
  select(2:21)-> df

lm1 <- lm(DNA_18S_Whole_Richness~., data = df)



df_16s = read.csv("SEM_data_for_16S.csv")
df_16s %>% colnames()
df_16s %>% as_tibble() %>% 
  mutate(ROV = Sample %>% str_sub(5,8)) %>% 
  ggplot(aes(x = DIC, y = RNA_16s_Whole_Richness))+
  geom_point(aes(color = ROV),size = 4)+
  geom_smooth(method = "lm")+                                    #### set: se = T
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size=20),
        axis.text.y = element_text(color="black", size=20))+
  theme(legend.position="none")+
  theme(axis.title=element_text(size=20))+
  theme(legend.title=element_text(size=20), legend.text=element_text(size=20),
        axis.title=element_text(size=20))+
  # theme(legend.title=element_text(size=20), legend.text=element_text(size=15),
  #       axis.title=element_text(size=20))+
  ggpubr::stat_cor(label.y= 3000 ,method = "pearson", size = 8)+
  labs(x = expression(DIC))+
  ylab("16S-RNA Richness") 




df_16s = read.csv("SEM_data_for_16S.csv")
df_16s %>% colnames()
df_16s %>% as_tibble() %>% 
  mutate(ROV = Sample %>% str_sub(5,8)) %>% 
  ggplot(aes(x = methane, y = RNA_16s_Whole_Richness))+
  geom_point(aes(color = ROV),size = 4)+
  geom_smooth(method = "lm")+                                    #### set: se = T
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size=15, hjust = 1),
        axis.text.y = element_text(color="black", size=15))+
  theme(legend.title=element_text(size=20), legend.text=element_text(size=15),
        axis.title=element_text(size=20))+
  ggpubr::stat_cor(label.y= 3000 ,method = "spearman", size = 8)+
  xlab("CH4")

df_16s %>% colnames()
df_16s %>% as_tibble() %>% 
  mutate(ROV = Sample %>% str_sub(5,8)) %>% 
  ggplot(aes(x = d15N, y = RNA_16s_Whole_Richness))+
  geom_point(aes(color = ROV),size = 4)+
  geom_smooth(method = "lm")+                                    #### set: se = T
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size=20),
        axis.text.y = element_text(color="black", size=20))+
  theme(legend.title=element_text(size=20), legend.text=element_text(size=20),
        axis.title=element_text(size=20))+
  ggpubr::stat_cor(label.y= 3000 ,method = "pearson", size = 8)+
  xlab(bquote(delta^15*N))+
  ylab("16S-RNA Richness")+
  scale_x_reverse()+
  theme(legend.position="none")


df_16s %>% as_tibble() %>% 
  mutate(ROV = Sample %>% str_sub(5,8)) %>% 
  ggplot(aes(x = DOM_fi, y = RNA_16s_Whole_Richness))+
  geom_point(aes(color = ROV),size = 4)+
  geom_smooth(method = "lm")+                                    #### set: se = T
  theme_classic()+
  theme(axis.text.x = element_text(color="black", size=20),
        axis.text.y = element_text(color="black", size=20))+
  theme(legend.title=element_text(size=20), legend.text=element_text(size=20),
        axis.title=element_text(size=20))+
  ggpubr::stat_cor(label.y= 2950 ,method = "pearson", size = 8)+
  xlab("DOM_FI")+
  ylab("16S-RNA Richness")


d# pielou evenness ---------------------------------------------------------


asv_16s_apha %>% 
  rownames_to_column(., var = "sample_id") %>% 
  as_tibble() %>%  
  mutate(ROV = sample_id %>% str_sub(5,8)) %>%
  mutate(Depth_cmbs = sample_id %>%
           str_split("_") %>% sapply('[', 2) %>% 
           str_replace("^5$","05")) %>%
  mutate(Ribosome = sample_id %>%
           str_split("_") %>% sapply('[', 3)) %>%
  mutate(Nucleic = sample_id %>%
           str_split("_") %>% sapply('[', 4)) %>% 
  filter(Nucleic == "rna" &
           Depth_cmbs < 20) %>%
  mutate(Habitat = ROV %>% 
           str_replace_all(c("ROV1" = "Seep", 
                             "ROV2" = "Seep",
                             "ROV3" = "Seep",
                             "ROV4" = "Non-seep",
                             "ROV5" = "Non-seep")) %>% as.factor())%>% filter(Habitat =="Non-seep") %>% 
  ggplot(aes(x = Habitat, y = Pielou)) +
  geom_boxplot(outlier.colour = "white", outlier.fill = "white", linewidth =2, color = "grey66")+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))+
  theme(axis.text.x = element_text(color="black", size=20),
        axis.text.y = element_text(color="black", size=20))+
  #annotate("text",label="ANOSIM_transplant: R = 0.19, p = 0.027",x=-0.72,y=0.58, size=5,fontface="bold")+
  #annotate("text",label="ANOSIM_tissue: R = 0.38, p = 0.001",x=-0.72,y=0.50, size=5,fontface="bold")+
  theme(legend.title=element_text(size=20), legend.text=element_text(size=20),
        axis.title=element_text(size=20))+
  ggtitle("16S RNA")+
  xlab("")+
  ylab("Pielou evenness")+
  ylim(0,1)-> p4

library(ggpubr)
ggarrange(p1,p3,p2,p4,nrow=2,ncol=2)



# plot of WEOM ------------------------------------------------------------


dom = read.table("WEOM compositions and charateristics.txt",header = T)
# dom %>% as_tibble() %>% 
#   mutate(ROV = Sample %>% str_sub(5,8)) %>% 
#   select(!Sample) %>% 
#   pivot_longer(!ROV, names_to = "DOM_type", values_to = "Value") %>% 
#   ggplot(aes(x = ROV, y = Value))+
#   facet_wrap(~DOM_type) +
#   geom_boxplot()+
#   theme_classic()+
#   theme(axis.text.x = element_text(color="black", size=15),
#         axis.text.y = element_text(color="black", size=15))+
#   theme(axis.text.x = element_text(color="black", size=15),
#         axis.title=element_text(size=20))+
#   theme(legend.title=element_text(size=20), legend.text=element_text(size=20),
#         axis.title=element_text(size=20))

dom1 <-
  dom %>% as_tibble() %>% 
  filter(DOM_a < 4) %>% ## remove 2 extreme values from ROV2
  rename(SR = DOM_sr) %>% 
  rename(BIX = DOM_bix) %>% 
  rename(HIX = DOM_hix) %>% 
  rename(FI = DOM_fi) %>% 
  column_to_rownames(., var = "Sample")
pca_res <- prcomp(dom1, scale = T)  

dom2<-
dom %>% as_tibble() %>% 
  filter(DOM_a < 4) %>%   ## remove 2 extreme values from ROV2
  mutate(ROV = Sample %>% str_sub(5,8)) %>% 
  mutate(DOM_SR = DOM_sr) %>% 
  mutate(DOM_BIX = DOM_bix) %>% 
  mutate(DOM_HIX = DOM_hix) %>% 
  mutate(DOM_FI = DOM_fi) %>% 
  column_to_rownames(., var = "Sample")


library(ggfortify)

autoplot(pca_res, data = dom2, colour = 'ROV',size = 5, alpha = 0.8,
         loadings = TRUE, loadings.colour = 'cornflowerblue', loadings.label.colour = "Black",
         loadings.label = TRUE, loadings.label.size = 5,
         loadings.arrow = grid::arrow(length = grid::unit(5, "points")))+
  geom_hline(yintercept = 0, linetype = 3, linewidth = 1)+
  geom_vline(xintercept = 0, linetype = 3, linewidth = 1)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.title = element_text(size = 20))+ 
  theme(axis.text.x = element_text(color="black", size=15),
        axis.text.y = element_text(color="black", size=15))+
  theme(legend.title=element_text(size=20), legend.text=element_text(size=20),
        axis.title=element_text(size=20))+
  scale_color_manual(values = c("firebrick1", "hotpink", "darkgoldenrod1", "cornflowerblue","darkcyan"))



# Comparison of env factors between ROVs ----------------------------------

df = read.csv("SEM_data_for_16S.csv",check.names = F)

df1<-
df %>% data.frame() %>% as_tibble() %>% 
  select(-Depth) %>% 
  pivot_longer(!c(Sample, ROV), names_to = "Env", values_to = "Value") %>% 
  na.omit() 




df1$Env = factor(df1$Env, levels = c("methane","sulphate","sulfide","ammonium","phosphate",
                                     "DIC", "d13C","d15N","X.C","X.N",
                                     "DOM_a", "DOM_b", "DOM_c","DOM_m","DOM_t",
                                     "DOM_bix", "DOM_fi","DOM_hix" ,"DOM_sr"))
df1 %>% 
  ggplot(aes(x = ROV, y = Value))+
  geom_boxplot(aes(color = ROV))+
  facet_wrap(~Env,scales="free")+
  scale_color_manual(values = c("firebrick1", "hotpink", "darkgoldenrod1", "cornflowerblue","darkcyan"))+
  scale_x_discrete(limits=c("ROV1","ROV2","ROV3","ROV4","ROV5"))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.title = element_text(size = 20))+ 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(color="black", size=15))+
  theme(legend.title=element_text(size=20), legend.text=element_text(size=20),
        axis.title=element_text(size=20))+
  xlab("")+
  ylab("")


df %>% colnames()

df %>% data.frame() %>% as_tibble() %>% 
  select(-Depth) %>% 
  pivot_longer(!c(Sample, ROV), names_to = "Env", values_to = "Value") %>% 
  na.omit() %>% 
  filter(Env == "??15N") %>% 
  filter(ROV %in% c("ROV4","ROV5")) %>% 
  pull(Value) %>% sd()

# relative abundane 

com_18s = read.csv("level-4-18s-pr2.csv",check.names = F)

com_18s<-
com_18s %>% column_to_rownames(., var = "index") %>%
  select(!c(ROV, Depth, Ribosome, Nucleic))
 
com_18s = rrarefy(com_18s, min(rowSums(com_18s)))
com_18s = com_18s/min(rowSums(com_18s))


# fungi: dna:31.72%, rna:19.84%
# Apicom: dna: all 4.7%; seep: 5.38%, nonseep:0.95%; rna: 0.85%; 
# Bre: dna: 0.34% (seep: 0.40%, nonseep: <0.01%);
#      rna: 3.11% (seep: 3.75%, nonseep: 0.15%)
# Cerco: dna:1.76%, rna: 13.31%
# unclassified Opis: dna: 2.62%, rna: 4.87% (seep: 5.42%, non-seep: 2.32%)

com_18s %>% data.frame() %>% 
  rownames_to_column(., var = "Sample") %>% as_tibble() %>% 
  mutate(ROV = Sample %>% str_sub(5,8)) %>% 
  mutate(Nucleic = Sample %>% str_split("_18s_") %>% sapply('[', 2)) %>% 
  select(Sample, ROV, Nucleic, 2:60) %>% 
  filter(Nucleic == "rna") %>% 
  #filter(ROV %in% c("ROV1","ROV2","ROV3")) %>% 
  #filter(ROV %in% c("ROV4","ROV5")) %>% 
  select(contains("Eukaryota.Obazoa.Opisthokonta.__")) %>% 
  pull(Eukaryota.Obazoa.Opisthokonta.__) %>% mean()



  
com_16s = read.csv("level-3-16s.csv",check.names = F)

com_16s<-
  com_16s %>% column_to_rownames(., var = "index") %>%
  select(!c(ROV, Depth, Ribosome, Nucleic))

com_16s = rrarefy(com_16s, min(rowSums(com_16s)))
com_16s = com_16s/min(rowSums(com_16s)) 

#ANME1: dna(rov12: 16.7%, ROV3: 2.55%, rov45: 0.08%);
#       rna(rov12: 28.63%, ROV3: 3.60%, rov45: 0.19%)
#Methano: dna(10.37%, seep:1.22%, non-seep: 0.011%)
#         rna(8.03%, seep:9.53%, 0.22%)
#Gamma: dna(16.94%, seep:15.78%, non-seep:23.32%)
#       rna(27.95%, seep:28.82%, non-seep:23.4%)

## higher in DNA, lower in RNA
#JS1: dna(6.6%, seep:7.69%, non-seep:0.54%)
#       rna(0.3%, seep:0.34%, non-seep:0.11%)

#Bacilli: dna(seep:2.69%, non-seep:14.98%)
#         rna(seep:0.05%, non-seep:0.05%)

#Anaerolineae: dna: 5.07%, rna: 1.12%

com_16s %>% data.frame() %>% 
  rownames_to_column(., var = "Sample") %>% as_tibble() %>% 
  mutate(ROV = Sample %>% str_sub(5,8)) %>% 
  mutate(Nucleic = Sample %>% str_split("_16s_") %>% sapply('[', 2)) %>% 
  select(Sample, ROV, Nucleic, 2:261) %>% 
  filter(Nucleic == "rna") %>% 
  # filter(ROV %in% c("ROV1","ROV2","ROV3")) %>%
  # filter(ROV %in% c("ROV4","ROV5")) %>%
  select(contains("d__Bacteria.p__Firmicutes.c__Bacilli")) %>% 
  pull(d__Bacteria.p__Firmicutes.c__Bacilli) %>% mean()




df = data.frame(a = c(rep(0,15), 1),
                b = c(rep(0,15), 1),
                c = c(rep(0,15), 1),
                d = c(rep(0,16))) %>% t()

df

df1 = df[rowSums(df == 0) < 16,]
df1


data.frame(s1 = c(1:10), 
           s2 = c(2:11), 
           s3 = c(3:12), 
           gp = c("A","A","B","B","C","C","D","D","E","E")) %>% 
  as_tibble() %>% 
  group_by(gp) %>% 
  summarise(av_s1 = mean(s1), 
            av_s2 = mean(s2), 
            av_s3 = mean(s3))
  

# RDA and CCA -------------------------------------------------------------
# asv18s_dna;
# asv18s_rna;
# asv16s_dna;
# asv16s_rna;

#18s
asv_18s = read.table("asv_18s.txt",header = T,check.names = F) %>% t()
tax_18s = read.delim("taxonomy-pr2.tsv")


## remove metazoan and land plants
tax_pro = tax_18s[!c(str_detect(tax_18s$Taxon, "Metazoa")),]
tax_pro = tax_pro[!c(str_detect(tax_pro$Taxon, "Streptophyta")),]
asv_18s = asv_18s[,match(tax_pro$Feature.ID, colnames(asv_18s))]
asv_18s = rrarefy(asv_18s,min(rowSums(asv_18s)))

asv18s_dna = asv_18s[asv_18s %>% rownames() %>% str_detect("dna"), ]
asv18s_dna %>% rownames()
asv18s_dna = asv18s_dna[,which(colSums(asv18s_dna)>0)]
asv18s_dna %>% dim
asv18s_dna = rrarefy(asv18s_dna, min(rowSums(asv18s_dna)))
asv18s_dna = asv18s_dna[,which(colSums(asv18s_dna)>0)]

asv18s_rna = asv_18s[asv_18s %>% rownames() %>% str_detect("rna"), ]
asv18s_rna %>% rownames()
asv18s_rna = asv18s_rna[,which(colSums(asv18s_rna)>0)]
asv18s_rna %>% dim
asv18s_rna = rrarefy(asv18s_rna, min(rowSums(asv18s_rna)))
asv18s_rna = asv18s_rna[,which(colSums(asv18s_rna)>0)]

#16s
asv_table_16s = read.table("asv_table_16s.txt",header = T,check.names = F) %>% t()
asv_table_16s %>% dim
tax_16s = read.delim("taxonomy-16s.tsv") # 108166 taxon
tax_16s = tax_16s[!c(str_detect(tax_16s$Taxon, "Chloroplast")),] # 108061 taxon remained
tax_16s = tax_16s[!c(str_detect(tax_16s$Taxon, "Mitochondria")),] #107915 taxon remained

asv_table_16s = asv_table_16s[,match(tax_16s$Feature.ID, colnames(asv_table_16s))]

asv16s_dna = asv_table_16s[asv_table_16s %>% rownames() %>% str_detect("dna"), ]
asv16s_rna = asv_table_16s[asv_table_16s %>% rownames() %>% str_detect("rna"), ]
asv16s_dna = rrarefy(asv16s_dna, min(rowSums(asv16s_dna)))
asv16s_dna = asv16s_dna[,which(colSums(asv16s_dna)>0)]
asv16s_rna = rrarefy(asv16s_rna, min(rowSums(asv16s_rna)))
asv16s_rna = asv16s_rna[,which(colSums(asv16s_rna)>0)]


env = read.csv("SEM_data_for_16S.csv")
env %>% head

env1 <- 
env %>% as_tibble() %>% 
  filter(Sample %in% c(asv16s_rna %>% rownames() %>% str_split("_16s_") %>% sapply('[', 1))) %>% 
  na.omit() %>% 
  select(1,4:22) %>% 
  rename(`%ON` = `X.ON`) %>%
  rename(`%OC` = `X.OC`) %>%
  column_to_rownames(., var = "Sample")

#16S DNA  
rownames(asv16s_dna) <- 
  asv16s_dna %>% rownames() %>% str_split("_16s_") %>% sapply('[', 1)

asv16s_dna = asv16s_dna[match(rownames(env1), rownames(asv16s_dna)),]

asv16s_dna = asv16s_dna[,which(colSums(asv16s_dna)>0)]

asv16s_dna %>% decorana() # use CCA

spp = decostand(asv16s_dna, method = "hellinger")

env1 = log(abs(env1)+1)

uu = rda(spp, env1)

uu

envfit(uu, env1)

ii = summary(uu)
ii
sp = as.data.frame(ii$species[,1:2])*2
st = as.data.frame(ii$sites[,1:2])
st$ROV = asv16s_dna %>% rownames() %>% str_sub(5,8) %>% str_replace("ROV5","ROV45")

yz = as.data.frame(ii$biplot[,1:2])
spp = as.data.frame(ii$species[,1:2])

library(ggrepel)

P1<-
ggplot()+
  geom_point(data = st, aes(RDA1, RDA2, 
                            color = ROV), size = 5, alpha = 0.7)+
  geom_segment(data = yz, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(angle = 22.5, length = unit(0.35,"cm"), 
                             type = "closed"), 
               linetype = 1, linewidth = 1, color = "deeppink2", alpha = 0.75)+
  geom_text_repel(data = yz, aes(RDA1, RDA2, label = rownames(yz)), size = 4.5)+
  #geom_label_repel(data = spp, aes(CCA1, CCA2, label = rownames(spp)), size = 6)+
  labs(x = paste("RDA1 (", format(100*ii$cont[[1]][2,1], digits = 3), "%)", sep = "" ))+
  labs(y = paste("RDA2 (", format(100*ii$cont[[1]][2,2], digits = 3), "%)", sep = "" ))+
  geom_hline(yintercept = 0, linetype = 3, size = 1)+
  geom_vline(xintercept = 0, linetype = 3, size = 1)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.title = element_text(size = 15))+ 
  theme(axis.text.x = element_text(color="black", size=15),
        axis.text.y = element_text(color="black", size=15))+
  theme(legend.title=element_text(size=15), legend.text=element_text(size=15),
        axis.title=element_text(size=15))+
  ggtitle("16S DNA") 


#16S RNA  
rownames(asv16s_rna) <- 
  asv16s_rna %>% rownames() %>% str_split("_16s_") %>% sapply('[', 1)

asv16s_rna = asv16s_rna[match(rownames(env1), rownames(asv16s_rna)),]

asv16s_rna = asv16s_rna[,which(colSums(asv16s_rna)>0)]

asv16s_rna %>% decorana() # use CCA

spp = decostand(asv16s_rna, method = "hellinger")

env1 = log(abs(env1)+1)

uu = rda(spp, env1)

ii = summary(uu)
ii
sp = as.data.frame(ii$species[,1:2])*2
st = as.data.frame(ii$sites[,1:2])
st$ROV = asv16s_rna %>% rownames() %>% str_sub(5,8) %>% str_replace("ROV5","ROV45")

yz = as.data.frame(ii$biplot[,1:2])
spp = as.data.frame(ii$species[,1:2])


envfit(uu, env1)

library(ggrepel)
P2<-
ggplot()+
  geom_point(data = st, aes(RDA1, RDA2, 
                             color = ROV), size = 5, alpha = 0.7)+
  geom_segment(data = yz, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(angle = 22.5, length = unit(0.35,"cm"), 
                             type = "closed"), 
               linetype = 1, linewidth = 1, color = "deeppink2", alpha = 0.75)+
  geom_text_repel(data = yz, aes(RDA1, RDA2, label = rownames(yz)), size = 4.5)+
  #geom_label_repel(data = spp, aes(CCA1, CCA2, label = rownames(spp)), size = 6)+
  labs(x = paste("RDA1 (", format(100*ii$cont[[1]][2,1], digits = 3), "%)", sep = "" ))+
  labs(y = paste("RDA2 (", format(100*ii$cont[[1]][2,2], digits = 3), "%)", sep = "" ))+
  geom_hline(yintercept = 0, linetype = 3, size = 1)+
  geom_vline(xintercept = 0, linetype = 3, size = 1)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.title = element_text(size = 15))+ 
  theme(axis.text.x = element_text(color="black", size=15),
        axis.text.y = element_text(color="black", size=15))+
  theme(legend.title=element_text(size=15), legend.text=element_text(size=15),
        axis.title=element_text(size=15))+
  ggtitle("16S RNA") 

#18S DNA  
rownames(asv18s_dna) <- 
  asv18s_dna %>% rownames() %>% str_split("_18s_") %>% sapply('[', 1)

asv18s_dna = asv18s_dna[match(rownames(env1), rownames(asv18s_dna)),]

asv18s_dna = asv18s_dna[,which(colSums(asv18s_dna)>0)]

asv18s_dna %>% decorana() # use CCA

spp = decostand(asv18s_dna, method = "hellinger")

env1 = log(abs(env1)+1)

uu = rda(spp, env1)

envfit(uu, env1)

ii = summary(uu)
ii
sp = as.data.frame(ii$species[,1:2])*2
st = as.data.frame(ii$sites[,1:2])
st$ROV = asv18s_dna %>% rownames() %>% str_sub(5,8) %>% str_replace("ROV5","ROV45")

yz = as.data.frame(ii$biplot[,1:2])
spp = as.data.frame(ii$species[,1:2])

library(ggrepel)
P3<-
ggplot()+
  geom_point(data = st, aes(RDA1, RDA2, 
                            color = ROV), size = 5, alpha = 0.7)+
  geom_segment(data = yz, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(angle = 22.5, length = unit(0.35,"cm"), 
                             type = "closed"), 
               linetype = 1, linewidth = 1, color = "royalblue", alpha = 0.75)+
  geom_text_repel(data = yz, aes(RDA1, RDA2, label = rownames(yz)), size = 4.5)+
  #geom_label_repel(data = spp, aes(CCA1, CCA2, label = rownames(spp)), size = 6)+
  labs(x = paste("RDA1 (", format(100*ii$cont[[1]][2,1], digits = 3), "%)", sep = "" ))+
  labs(y = paste("RDA2 (", format(100*ii$cont[[1]][2,2], digits = 3), "%)", sep = "" ))+
  geom_hline(yintercept = 0, linetype = 3, size = 1)+
  geom_vline(xintercept = 0, linetype = 3, size = 1)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.title = element_text(size = 15))+ 
  theme(axis.text.x = element_text(color="black", size=15),
        axis.text.y = element_text(color="black", size=15))+
  theme(legend.title=element_text(size=15), legend.text=element_text(size=15),
        axis.title=element_text(size=15))+
  ggtitle("18S DNA") 

#18S rna  
rownames(asv18s_rna) <- 
  asv18s_rna %>% rownames() %>% str_split("_18s_") %>% sapply('[', 1)

asv18s_rna = asv18s_rna[match(rownames(env1), rownames(asv18s_rna)),]

asv18s_rna = asv18s_rna[,which(colSums(asv18s_rna)>0)]

asv18s_rna %>% decorana() # use CCA

spp = decostand(asv18s_rna, method = "hellinger")

env1 = log(abs(env1)+1)

uu = rda(spp, env1)

ii = summary(uu)
ii
sp = as.data.frame(ii$species[,1:2])*2
st = as.data.frame(ii$sites[,1:2])
st$ROV = asv18s_rna %>% rownames() %>% str_sub(5,8) %>% str_replace("ROV5","ROV45")

yz = as.data.frame(ii$biplot[,1:2])
spp = as.data.frame(ii$species[,1:2])

library(ggrepel)
P4<-
ggplot()+
  geom_point(data = st, aes(RDA1, RDA2, 
                            color = ROV), size = 5, alpha = 0.7)+
  geom_segment(data = yz, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(angle = 22.5, length = unit(0.35,"cm"), 
                             type = "closed"), 
               linetype = 1, linewidth = 1, color = "royalblue", alpha = 0.75)+
  geom_text_repel(data = yz, aes(RDA1, RDA2, label = rownames(yz)), size = 4.5)+
  #geom_label_repel(data = spp, aes(CCA1, CCA2, label = rownames(spp)), size = 6)+
  labs(x = paste("RDA1 (", format(100*ii$cont[[1]][2,1], digits = 3), "%)", sep = "" ))+
  labs(y = paste("RDA2 (", format(100*ii$cont[[1]][2,2], digits = 3), "%)", sep = "" ))+
  geom_hline(yintercept = 0, linetype = 3, size = 1)+
  geom_vline(xintercept = 0, linetype = 3, size = 1)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.title = element_text(size = 15))+ 
  theme(axis.text.x = element_text(color="black", size=15),
        axis.text.y = element_text(color="black", size=15))+
  theme(legend.title=element_text(size=15), legend.text=element_text(size=15),
        axis.title=element_text(size=15))+
  ggtitle("18S RNA") 

library(ggpubr)
ggarrange(P1,P2,P3,P4,nrow = 2,ncol = 2)



# Partial effects between richness and env --------------------------------

library(mgcv)
library(ggplot2)
library(itsadug)

df_18s = read.csv("SEM_data_for_18S.csv")
df_16s = read.csv("SEM_data_for_16S.csv")

#for 18S-dna-richness
mod_18s = lm(DNA_18S_Whole_Richness ~ ammonium + d13C+DOM_sr+DOM_fi, data = df_18s)
library(car)
vif(mod_18s) ## all < 2
gamm_18s_dna<-
  gamm(DNA_18S_Whole_Richness~s(ammonium,bs="cr")+s(d13C,bs="cr")+s(DOM_sr,bs="cr")+s(DOM_fi,bs="cr"),
       data=df_18s,method='REML')
summary(gamm_18s_dna$lme)
summary(gamm_18s_dna$gam)
par(mfrow=c(2,2))  
plot(gamm_18s_dna$gam,shade=TRUE)

#for 18S-rna
mod_18s = lm(RNA_18S_Whole_Richness ~ X.N + X.C, data = df_18s)
library(car)
vif(mod_18s) ## BOTH 3.5
gamm_18s_rna<-
  gamm(RNA_18S_Whole_Richness~s(X.N)+s(X.C),
       data=df_18s,method='REML')
summary(gamm_18s_rna$lme)
summary(gamm_18s_rna$gam)
par(mfrow=c(2,2))  
plot(gamm_18s_rna$gam,shade=TRUE)

#for 16S-dna
mod_16s = lm(DNA_16s_Whole_Richness ~ ??13C+X.OC, data = df_16s)
library(car)
vif(mod_16s) ## <2
gamm_16s_dna<-
  gamm(DNA_16s_Whole_Richness~s(??13C)+s(X.OC),
       data=df_16s,method='REML')
summary(gamm_16s_dna$lme)
summary(gamm_16s_dna$gam)
par(mfrow=c(2,2))  
plot(gamm_16s_dna$gam,shade=TRUE)

#for 16S-rna
mod_16s = lm(RNA_16s_Whole_Richness~Sulfide+Sulphate+Methane+DIC+??15N+??13C+X.ON+
               DOM_b+DOM_SR+DOM_BIX+DOM_FI+DOM_HIX, data = df_16s)
library(car)
vif(mod_16s) ## remove > 5
gamm_16s_rna<-
  gamm(RNA_16s_Whole_Richness~s(Sulfide)+s(??15N)+s(??13C)+s(X.ON)+
         s(DOM_b)+s(DOM_SR)+s(DOM_BIX)+s(DOM_FI)+s(DOM_HIX),
       data=df_16s,method='REML')
summary(gamm_16s_rna$lme)
summary(gamm_16s_rna$gam)
par(mfrow=c(3,3))  
plot(gamm_16s_rna$gam,shade=TRUE)


# STAMP analysis ----------------------------------------------------------

# STAMP -------------------------------------------------------------------
setwd("/Users/xzm/Desktop/cold seep 2nd MS/R codes and files/")
library(vegan)
library(ggplot2)
library(tidyverse)
library(boot)
#18s
com = read.csv("level-4-18s-pr2.csv",check.names = F)
com = com[,c(1:60)] %>% as_tibble() %>% 
  select(!contains("Metazoa")) %>% 
  select(!contains("Streptophyta")) %>% 
  column_to_rownames(., var = "index")

com = rrarefy(com, min(rowSums(com)))
com = com/min(rowSums(com))
# extract the most abundant 15 taxa, with sum of remaining as "Others"
com1 = com %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:11] %>% data.frame() %>% rownames()
m = match(com1, colnames(com))
rownames(com) = gsub("_5_", "_05_",rownames(com))
# barplot of community composition-relative abundance
# library(palettes)
data1 = 
  cbind(com[,m], data.frame(apply(com[,-m],1,sum)))%>% 
  rename(Others = `apply.com....m...1..sum.`) %>% t() %>% 
  data.frame() %>% 
  rownames_to_column(., var = "Microeukaryotes") %>% as_tibble() %>% 
  pivot_longer(!Microeukaryotes, names_to = "Sample", values_to = "Seqs") %>% 
  mutate(ROV = Sample %>% str_sub(5,8)) %>% 
  mutate(Sample = Sample %>% str_replace("_5_", "_05_")) %>% 
  mutate(Nucleic = Sample %>% str_split("_") %>% sapply('[', 4)) %>% 
  mutate(Depth_cmbs = Sample %>% str_split("_") %>% sapply('[', 2)) %>% 
  filter(Nucleic == "dna") %>% 
  select(Sample, Microeukaryotes, Seqs) %>% 
  pivot_wider (names_from = Microeukaryotes, values_from = Seqs) %>% 
  column_to_rownames(., var = "Sample")

grp<-
  data1 %>% rownames() %>% str_sub(5,8) %>% str_replace_all(c("ROV1"="Seep", 
                                                              "ROV2"="Seep",
                                                              "ROV3"="Seep",
                                                              "ROV4"="Non-seep",
                                                              "ROV5"="Non-seep"))

data1$Group = factor(grp, levels = c("Seep","Non-seep"))

diff <- data1 %>%
  select_if(is.numeric) %>%
  map_df(~ broom::tidy(t.test(. ~ Group,data = data1)), .id = 'var')


diff$p.value <- p.adjust(diff$p.value,"bonferroni")
# diff <- diff %>% filter(p.value < 0.05)


abun.bar <- data1[,c(diff$var,"Group")] %>%
  gather(variable,value,-Group) %>%
  group_by(variable,Group) %>%
  summarise(Mean = mean(value))

diff.mean <- diff[,c("var","estimate","conf.low","conf.high","p.value")]
diff.mean$Group <- c(ifelse(diff.mean$estimate >0,levels(data1$Group)[1],
                            levels(data1$Group)[2]))
diff.mean$Group = factor(diff.mean$Group, levels = c("Seep","Non-seep"))
diff.mean <- diff.mean[order(diff.mean$estimate,decreasing = TRUE),]


library(ggplot2)
cbbPalette <- c("#E69F00", "#56B4E9")
abun.bar$variable <- factor(abun.bar$variable,levels = fix_name)

diff.mean$var -> fix_name ## use for next round

p11 <- ggplot(abun.bar,aes(variable,Mean,fill = Group)) +
  scale_x_discrete(limits = levels(diff.mean$var)) +
  coord_flip() +
  xlab("") +
  ylab("Mean proportion (%)") +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.title.x=element_text(colour='black', size=12,face = "bold"),
        axis.text=element_text(colour='black',size=10,face = "bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=12,face = "bold",colour = "black", margin = margin(r = 20)),
        legend.position = c(-1,-0.1),
        legend.direction = "horizontal",
        legend.key.width = unit(0.8,"cm"),
        legend.key.height = unit(0.5,"cm"))


p11 <- p11 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf,
                      fill = ifelse(i %% 2 == 0, 'white', 'gray95'))

p11 <- p11 +
  geom_bar(stat = "identity",position = "dodge",width = 0.7,colour = "black") +
  scale_fill_manual(values=cbbPalette)

p11


diff.mean$var <- factor(diff.mean$var,levels = levels(abun.bar$variable))
diff.mean$p.value <- signif(diff.mean$p.value,3)
diff.mean$p.value <- as.character(diff.mean$p.value)

p22 <- ggplot(diff.mean,aes(var,estimate,fill = Group)) +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.title.x=element_text(colour='black', size=12,face = "bold"),
        axis.text=element_text(colour='black',size=10,face = "bold"),
        axis.text.y = element_blank(),
        legend.position = "none",
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 15,face = "bold",colour = "black",hjust = 0.5)) +
  scale_x_discrete(limits = levels(diff.mean$var)) +
  coord_flip() +
  xlab("") +
  ylab("Difference in mean proportions (%)") +
  labs(title="95% confidence intervals")

for (i in 1:(nrow(diff.mean) - 1))
  p22 <- p22 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf,
                        fill = ifelse(i %% 2 == 0, 'white', 'gray95'))

p22 <- p22 +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                position = position_dodge(0.8), width = 0.5, size = 0.5) +
  geom_point(shape = 21,size = 3) +
  scale_fill_manual(values=cbbPalette) +
  geom_hline(aes(yintercept = 0), linetype = 'dashed', color = 'black')

p22

p3 <- ggplot(diff,aes(var,estimate,fill = Group)) +
  geom_text(aes(y = 0,x = var),label = diff$p.value,
            hjust = 0,fontface = "bold",inherit.aes = FALSE,size = 3) +
  coord_flip() +
  ylim(c(0,1)) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())

p3

library(ggpubr)
ggarrange(p1,p2,p11, p22,nrow = 2,ncol = 2)

## 16S

com = read.csv("level-3-16s.csv",check.names = F)
com = com[,1:261] %>% column_to_rownames(., var = "index")
com = rrarefy(com, min(rowSums(com)))
com = com/46886
# extract the most abundant 15 taxa, with sum of remaining as "Others"
com1 = com %>% colSums() %>% sort(decreasing = TRUE) %>% .[1:11] %>% data.frame() %>% rownames()
m = match(com1, colnames(com))
rownames(com) = gsub("_5_", "_05_",rownames(com))

data1 = 
  cbind(com[,m], data.frame(apply(com[,-m],1,sum)))%>% 
  rename(Others = `apply.com....m...1..sum.`) %>% t() %>% 
  data.frame() %>% 
  rownames_to_column(., var = "Prokaryotes") %>% as_tibble() %>% 
  pivot_longer(!Prokaryotes, names_to = "Sample", values_to = "Seqs") %>% 
  mutate(ROV = Sample %>% str_sub(5,8)) %>% 
  mutate(Sample = Sample %>% str_replace("_5_", "_05_")) %>% 
  mutate(Nucleic = Sample %>% str_split("_") %>% sapply('[', 4)) %>% 
  mutate(Depth_cmbs = Sample %>% str_split("_") %>% sapply('[', 2)) %>% 
  filter(Nucleic == "dna") %>% 
  select(Sample, Prokaryotes, Seqs) %>% 
  pivot_wider (names_from = Prokaryotes, values_from = Seqs) %>% 
  column_to_rownames(., var = "Sample")

grp<-
  data1 %>% rownames() %>% str_sub(5,8) %>% str_replace_all(c("ROV1"="Seep", 
                                                              "ROV2"="Seep",
                                                              "ROV3"="Seep",
                                                              "ROV4"="Non-seep",
                                                              "ROV5"="Non-seep"))

data1$Group = factor(grp, levels = c("Seep","Non-seep"))

diff <- data1 %>%
  select_if(is.numeric) %>%
  map_df(~ broom::tidy(t.test(. ~ Group,data = data1)), .id = 'var')


diff$p.value <- p.adjust(diff$p.value,"bonferroni")
# diff <- diff %>% filter(p.value < 0.05)


abun.bar <- data1[,c(diff$var,"Group")] %>%
  gather(variable,value,-Group) %>%
  group_by(variable,Group) %>%
  summarise(Mean = mean(value))

diff.mean <- diff[,c("var","estimate","conf.low","conf.high","p.value")]
diff.mean$Group <- c(ifelse(diff.mean$estimate >0,levels(data1$Group)[1],
                            levels(data1$Group)[2]))
diff.mean$Group = factor(diff.mean$Group, levels = c("Seep","Non-seep"))
diff.mean <- diff.mean[order(diff.mean$estimate,decreasing = TRUE),]


library(ggplot2)
cbbPalette <- c("#E69F00", "#56B4E9")
abun.bar$variable <- factor(abun.bar$variable,levels = fix_name) ## save "diff.mean$var" and use the same order for rna part

rev(diff.mean$var) -> fix_name ## use for next round

q11 <- ggplot(abun.bar,aes(variable,Mean,fill = Group)) +
  scale_x_discrete(limits = levels(diff.mean$var)) +
  coord_flip() +
  xlab("") +
  ylab("Mean proportion (%)") +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.title.x=element_text(colour='black', size=12,face = "bold"),
        axis.text=element_text(colour='black',size=10,face = "bold"),
        legend.title=element_blank(),
        legend.text=element_text(size=12,face = "bold",colour = "black", margin = margin(r = 20)),
        legend.position = c(-1,-0.1),
        legend.direction = "horizontal",
        legend.key.width = unit(0.8,"cm"),
        legend.key.height = unit(0.5,"cm"))


q11 <- q11 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf,
                      fill = ifelse(i %% 2 == 0, 'white', 'gray95'))

q11 <- q11 +
  geom_bar(stat = "identity",position = "dodge",width = 0.7,colour = "black") +
  scale_fill_manual(values=cbbPalette)

q11


diff.mean$var <- factor(diff.mean$var,levels = levels(abun.bar$variable))
diff.mean$p.value <- signif(diff.mean$p.value,3)
diff.mean$p.value <- as.character(diff.mean$p.value)

q22 <- ggplot(diff.mean,aes(var,estimate,fill = Group)) +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.title.x=element_text(colour='black', size=12,face = "bold"),
        axis.text=element_text(colour='black',size=10,face = "bold"),
        axis.text.y = element_blank(),
        legend.position = "none",
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 15,face = "bold",colour = "black",hjust = 0.5)) +
  scale_x_discrete(limits = levels(diff.mean$var)) +
  coord_flip() +
  xlab("") +
  ylab("Difference in mean proportions (%)") +
  labs(title="95% confidence intervals")

for (i in 1:(nrow(diff.mean) - 1))
  q22 <- q22 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf,
                        fill = ifelse(i %% 2 == 0, 'white', 'gray95'))

q22 <- q22 +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                position = position_dodge(0.8), width = 0.5, size = 0.5) +
  geom_point(shape = 21,size = 3) +
  scale_fill_manual(values=cbbPalette) +
  geom_hline(aes(yintercept = 0), linetype = 'dashed', color = 'black')

q22


ggarrange(q1,q2,q11,q22, nrow=2,ncol=2)

