#
# install.packages(c("tidyverse", "impute", "Rcpp"))
# BiocManager::install("methylclock")
#
library(methylclockData)
library(methylclock)
library(Biobase)
library(tibble)
library(impute)
library(ggplot2)
library(ggpmisc)
library(GEOquery)
library(WriteXLS)
library(ggplot2)
library(RColorBrewer)
library(doParallel)
library(PCAtools)
cl <- makeCluster(6)
registerDoParallel(cl)
# read the input mattrix file
m_matrix <- read.table("series_matrix 2.txt", header = T, sep = "\t", row.names =1)
m_matrix <- as.matrix(m_matrix)
class(m_matrix)
#
checkClocks(m_matrix) # Your data contain the required CpGs for all clocks
#
chrono_age <- DNAmAge(m_matrix)
numeric_cols <- sapply(chrono_age, is.numeric)
chrono_age[, numeric_cols] <- round(chrono_age[, numeric_cols])
chrono_age
#
WriteXLS(c("chrono_age"),"predicted_age.xlsx",row.names=F,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)

###############################################################
###############################################################
###############################################################

chrono_age2 <- chrono_age[,c("id","Horvath")]
chrono_age2 <- as.data.frame(chrono_age2)
#
# chrono_age2$group <- rep(paste0("S", rep(seq(1:10), each = 3)))
# chrono_age2$group <- factor(chrono_age2$group,
#                             levels = paste0("S", seq(1:10)))
chrono_age2$group <- c(paste0("S", rep(1,3)),
                       paste0("S", rep(2,2)),
                       paste0("S", rep(3,3)),
                       paste0("S", rep(4,3)),
                       paste0("S", rep(5,2)),
                       paste0("S", rep(6,2)),
                       paste0("S", rep(7,3)),
                       paste0("S", rep(8,3)),
                       paste0("S", rep(9,3)),
                       paste0("S", rep(10,3)),
                       paste0("S", rep(11,3))
)
                       
chrono_age2$group <- factor(chrono_age2$group,
                           levels = paste0("S", seq(1:11)))

chrono_age2$id <- factor(chrono_age2$id, levels = chrono_age2$id)

#my_colors <- rep(brewer.pal(10, name = "Paired"), each=3)
my_colors <-  c(rep("#F08080", each=3),
                rep("#FFDAB9", each=2),
                rep("#87CEFA", each=3),
                rep("#FFA07A", each=3),
                rep("#FFE4C4", each=2),
                rep("#AFEEEE", each=2),
                rep("#D2B48C", each=3),
                rep("#E6E6FA", each=3),
                rep("#FFC0CB", each=3),
                rep("#BA55D3", each=3),
                rep("#00CED1", each=3)
)

#  bar plot #
pdf("methyl_age_barplot.pdf",width = 9,height = 9)

ggplot(chrono_age2, aes(x = id, y = Horvath, fill=group, colour=group)) +
  geom_bar(stat = "identity" , position = "dodge") +
  labs(x = "", y = "Horvath Age") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_grid(~group, scales = "free_x", space = "free_x")

dev.off()
##
#gesta_age <- DNAmGA(m_matrix)
# #
# chrono_age2 <- DNAmAge(m_matrix, normalize = T)
# numeric_cols <- sapply(chrono_age2, is.numeric)
# chrono_age2[, numeric_cols] <- round(chrono_age2[, numeric_cols])
# chrono_age2$Horvath - chrono_age$Horvath

###################################################################
#######.  PCA plot
###############################################################
#pcaData <- PCAtools::pca(log2(m_matrix[,]), removeVar = 0.9)
pcaData <- PCAtools::pca(m_matrix, removeVar = 0.9)
mylabels <- rownames(pcaData$rotated)
mylabels2 <- factor(mylabels, levels = mylabels)
#my_colors <- rep(brewer.pal(10, name = "Paired"), each=3)

# my_colors2 <- my_colors[1,4,7,10,13,16,19,22,25,28,
#                         2,5,8,11,14,17,20,23,26,29,
#                         3,6,9,12,15,18,21,24,27,30]

pcaData$loadings
pdf("methyl_30_pca.pdf",width = 9,height = 9)
p <-  PCAtools::biplot(pcaData,
                       # lab = NULL,
                       #colby = chrono_age2$group,
                       colkey = my_colors,
                       #shape = '',
                       hline = 0, vline = 0,
                       legendPosition = NULL,
                       title = "",
                       encircle=F, ellipse = F,encircleFill = FALSE,
                       encircleAlpha = 1, encircleLineSize = 5,
                       lab=mylabels2)

p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())
dev.off()
###############################################################
##. top/bottom 10 cpGs , bottom can be changed 
##.  according to variance kept in PCA calculation
#################################################################
load_score <- pcaData$loadings[, 1, drop =F]
load_score_rank <- load_score[order(load_score$PC1, decreasing =T), , drop=F]
#
top10cpGs <- rownames(head(load_score_rank, n=10))
bottom10cpGs <- rownames(tail(load_score_rank, n=10))

top10cpGs %>% as.data.frame()
bottom10cpGs
###
# BiocManager::install(c("IlluminaHumanMethylation450kanno.ilmn12.hg19",
#                        "IlluminaHumanMethylation450kmanifest",
#                        "FDb.InfiniumMethylation.hg19"
#                        ))
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(minfi)
library(FDb.InfiniumMethylation.hg19)

hm450 <- get450k()

###############################################################
## top 10 cpGs
###############################################################

probes <- hm450[top10cpGs]

top10_near_tss <- getNearestTSS(probes)
#getNearestTranscript(probes)

# 
top10_anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19,
              lociNames = top10cpGs)
top10_anno <- as.data.frame(top10_anno)
#getLocations(IlluminaHumanMethylation450kanno.ilmn12.hg19, mergeManifest = FALSE,
#             orderByLocation = FALSE, lociNames = top10cpGs)

WriteXLS(c("top10_near_tss","top10_anno"),"top10_cpGs_anno.xlsx",row.names=F,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)

###############################################################
## bottom 10 cpGs
###################################

probes <- hm450[bottom10cpGs]

bottom10_near_tss <- getNearestTSS(probes)
#getNearestTranscript(probes)
# 
bottom10_anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19,
                               lociNames = bottom10cpGs)
bottom10_anno <- as.data.frame(bottom10_anno)
#getLocations(IlluminaHumanMethylation450kanno.ilmn12.hg19, mergeManifest = FALSE,
#             orderByLocation = FALSE, lociNames = bottom10cpGs)

WriteXLS(c("bottom10_near_tss","bottom10_anno"),"bottom10_cpGs_anno.xlsx",row.names=F,BoldHeaderRow=T,FreezeRow=1,AdjWidth=T)
#######

