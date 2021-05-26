setwd("/Users/wittkopp_member/Code")

# libraries
library(binom)
library(Hmisc)
library(ggplot2)
library(tidyverse)

theme_main <- function() {
  theme_bw() +
  theme(
  #panel.grid.major = element_blank(),
  #panel.grid.minor = element_blank(),
  axis.title = element_text(size = 20),
  strip.text = element_text(size = 20),
  legend.text= element_text(size = 15),
  legend.title = element_text(size = 20),
  plot.title = element_text(size = 13, face = "bold"),
  axis.text.x = element_text(size = 15),
  axis.text.y = element_text(size = 15)
)
}


## do house-keeping genes have fewer enhancers than dev regulated?
# read in single cell data
sc_df <- readRDS(file = "/Users/wittkopp_member/Downloads/Correlation_acc-act.Rds") %>% as.data.frame()
sc_df$score <- as.numeric(sc_df$score)
sc_df <- sc_df[sc_df$score > 0,]
genes <- sc_df[,4] %>% as.data.frame()
colnames(genes)[1] <- "gene"
genes <- ddply(genes,.(gene),nrow)
colnames(genes) <- c("Gene_name","count")

wing_disc_MM <- read.delim("./Fledgling_analysis_ideas/Perc_cis_by_tissue_specificity/Data_tables_generated/Ertl_MM_wing_disc_gene_groups.txt", header = T)
colnames(wing_disc_MM)[2:5] <- c("P1", "P2", "HYB1", "HYB2")
classes <- wing_disc_MM[,c(10,15)]

whole_body_ALL <- right_join(classes, genes, by = "Gene_name") %>% unique()
whole_body_ALL$class[is.na(whole_body_ALL$class)] <- "non_TS_non_HK"

enh_count_plotting <- melt(whole_body_ALL,id.vars = "count", measure.vars = c("class"))
enh_count_plotting[enh_count_plotting$value != "TS",]  %>%
ggplot(aes(x=value, y=count, fill=value)) +
geom_boxplot(notch=TRUE) +
theme_main() +
xlab("") +
ylab("# of enhancers in eye-antennal discs") +
scale_x_discrete(labels=c("Const. active", "Developmentally-\nregulated:\nmany tissues", "Imaginal disc specific")) +
scale_fill_discrete(guide=FALSE)

library(data.table)
count.matrix <- suppressWarnings(data.frame(fread('/Users/wittkopp_member/Downloads/GSM4209233_10X_SCATACSEQ_WT_EA_CTXREGIONS.tsv', sep='\t', verbose=F), row.names=1))
library(cisTopic)
# Initialize object
cisTopicObject <- createcisTopicObject(count.matrix, project.name='EAdisc')
# Run models (updated to cisTopic v3)
cisTopicObject <- runCGSModels(cisTopicObject, topic=c(2, 10, 20, 30:50, 60, 70, 80, 90, 100), seed=987, nCores=10, burnin = 250, iterations = 500)
# Select model
cisTopicObject <- selectModel(cisTopicObject, select=49, keepModels = FALSE)
# Remove data slots
cisTopicObject@count.matrix <- NULL
cisTopicObject@binary.count.matrix <- NULL
# Add run as metadata
run <- as.data.frame(as.factor(sapply(strsplit(cisTopicObject@cell.names, split = "[.]"), "[", 2)))
rownames(run) <- cisTopicObject@cell.names
cisTopicObject <- addCellMetadata(cisTopicObject, run)
# tSNE
cisTopicObject <- runtSNE(cisTopicObject, target='cell', perplexity=100)
# Louvain clustering (with Seurat 2.3.4)
devtools::install_version("Seurat", version = "2.3.4", repos = "http://cran.us.r-project.org")
library(Seurat)
topicCell <- modelMatSelection(cisTopicObject, 'cell', 'Z-score')
seuratObj <- CreateSeuratObject(raw.data = topicCell, min.cells = 0, min.genes = 0, project = "cisTopic_cluster")
seuratObj  <- FindClusters(seuratObj, genes.use=rownames(topicCell), resolution = 1.2)
# Update to Seurat v3
detach("package:Seurat", unload=TRUE)
remove.packages("Seurat", lib="~/R/x86_64-generic-linux-gnu-library/3.4")
install.packages('Seurat')
library(Seurat)
seuratObj <- UpdateSeuratObject(seuratObj)
# Annotate
new.cluster.ids <- c('PMF_Interommatidial', 'AMF_Prog', 'Antenna_A3_Arista', 'Antenna_A1', 'MF_Morphogenetic_Furrow', 'PMF_Interommatidial_Late', 'Antenna_A2a', 'Peripodial_membrane_medial', 'AMF_Prec', 'Antenna_A2b', 'Head_vertex', 'Peripodial_membrane_lateral', 'PMF_PR_Early', 'PMF_PR_Late/CC', 'Unknown_B', 'Glia', 'Unknown_A', 'Brain_A', 'Brain_B', 'twi_cells', 'Hemocytes', 'Unknown_C')
# Set labels on the Seurat object
names(x = new.cluster.ids) <- levels(x = seuratObj)
seuratObj <- RenameIdents(object = seuratObj, new.cluster.ids)
# Add cisTopic coordinates
DimReduc <- setClass(Class = 'DimReduc', slots = c(cell.embeddings = 'matrix', feature.loadings = 'matrix', feature.loadings.projected = 'matrix', assay.used = 'character', global = 'logical', stdev = 'numeric',key = 'character', misc = 'list', jackstraw='ANY'))
tsne_coords <- cisTopicObject@dr$cell$tSNE
colnames(tsne_coords) <- c('tSNE_1', 'tSNE_2')
seuratObj@reductions$tsne <- new('DimReduc', cell.embeddings=tsne_coords, assay.used ='RNA', key='tSNE_')
# Create colVars
colors <- c('#5CC1A4', '#9A52A0', '#77A0D4', '#99BE3F', '#D1992A', '#B4DCB8', '#F68C61', '#42B86E', '#FCC39E', '#A8D05A','#42C6EC', '#63BBD1', '#82D3EA', '#FFFF00', '#E3D022', '#D077AF', '#70CCD8', '#FD79A4', '#B5DEFF', '#FBC8FB', '#3CBC0D', '#FC9CF7')
names(colors) <- c('Peripodial_membrane_medial', 'Peripodial_membrane_lateral', 'Head_vertex', 'Antenna_A1','Antenna_A2a','Antenna_A2b', 'Antenna_A3_Arista', 'AMF_Prec', 'AMF_Prog', 'MF_Morphogenetic_Furrow', 'PMF_PR_Early', 'PMF_PR_Late/CC', 'PMF_Interommatidial', 'PMF_Interommatidial_Late', 'Glia', 'twi_cells', 'Hemocytes', 'Brain_A', 'Brain_B', 'Unknown_A', 'Unknown_B', 'Unknown_C')
# Save objects
saveRDS(colors, file='Figure_2/Processed_data/cisTopic/Cell_Type_ColVars.Rds')
saveRDS(seuratObj, file='Figure_2/Processed_data/cisTopic/LouvainClustering.Rds')
