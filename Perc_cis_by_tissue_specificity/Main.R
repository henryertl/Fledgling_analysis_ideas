# main R code bringing together data from:
### Kondo et al. 2017 G&D for tissue specificity info
### Coolon et al. 2014 for whole body RNA-seq
### own data for wing disc RNA-seq

# set working directory
setwd("/Users/wittkopp_member/Code")

# libraries
library(binom)
library(Hmisc)
library(ggplot2)
library(tidyverse)

##########################################################################################
######### STEP 1 categorize genes by expression profile using Kondo et al info ###########
##########################################################################################

# read in Kondo et al. data
df <- read.delim("./DDivergence/Fly_duplicate_evol/Existing_Datasets/Fly_gene_ages_tau_G_17_TableS1.txt", header = TRUE)

# subset only expression columns to calculate minimum expression for each row across columns (value of gene is lowest expressed tissue)
df_subexpr <- df[,c(8:ncol(df))]

# function to calculate above:
rowMins <- function(df) { do.call(pmin, df) }

# deploy function and add column
df_subexpr$min <- rowMins(df_subexpr)

# add column categorizing as housekeeping or non-housekeeping based on above caluclated min # -- chose mean of min as cutoff: this gives ~estimated # of housekeeping genes
df_subexpr$HK_bin <- NA
for (i in 1:nrow(df_subexpr)) {

if (df_subexpr$min[i] > 0.7943){

	df_subexpr$HK_bin[i] <- "HK"

} else {

  df_subexpr$HK_bin[i] <- "non_HKs"

}
}

# append the housekeeping or non-housekeeping ID back to main df
df$HK_bin <- df_subexpr$HK_bin
df$min <- df_subexpr$min


# subset only relevant columns for downstream analysis becuase done with the all-tissue-data
df_cissub <- df[,c(1,2,4,6,7,(ncol(df) - 1), ncol(df))]

# double check number.. should be around 300
nrow(df_cissub[df_cissub$HK_bin == "HK",])

# now using all info, categorize into tissue specific, housekeeping, etc..
df_cissub$class <- NA
for (i in 1:nrow(df_cissub)) {

if (df_cissub$tau[i] > 0.61){
	df_cissub$class[i] <- "TS"

} else if (df_cissub$tau[i] > 0.61 & df_cissub$HK_bin[i] == "HK"){
  df_cissub$class[i] <- "TS_HK"

} else if (df_cissub$tau[i] < 0.61 & df_cissub$HK_bin[i] == "HK"){
  df_cissub$class[i] <- "non_TS_HK"

} else if (df_cissub$tau[i] < 0.61 & df_cissub$HK_bin[i] == "non_HKs"){
  df_cissub$class[i] <- "non_TS_non_HK"

}
}

write.table(df_cissub, file = "./Fledgling_analysis_ideas/Perc_cis_by_tissue_specificity/Data_tables_generated/Genes_grouped_by_tissue_exp_profile.txt", row.names = F, quote = FALSE, sep = "\t")

######### done with categorizing genes by expression profile ###########

##################################################################################
######## STEP 2: Read in whole body RNA-seq data from Coolon et al 2014 #########
##################################################################################

# read in coolon et al data
df2 <- read.delim("./Fledgling_analysis_ideas/Perc_cis_by_tissue_specificity/Existing_datasets/Supplemental_Dataset2.txt", header = TRUE)

# df2 contains info for 3 spp comparisons - make three separate dfs now to analyze separately **ignoring reciprocal cross info and just choosing one
df2_MM <- df2[,c(1:5)]
df2_SS <- df2[,c(1,8:11)]
df2_MS <- df2[,c(1,14:17)]

# functions to remove < 20 reads and calculate percent cis
remove_less_20 <- function(df) { df[df[,2] >= 20 & df[,3] >= 20 & df[,4] + df[,5] >= 20,] }
parent_ratios <- function(df) { log2_P_ratio <- log2(df[,2]/df[,3]); cbind(df,log2_P_ratio)}
hybrid_ratios <- function(df) { log2_H_ratio <- log2(df[,4]/df[,5]); cbind(df,log2_H_ratio) }
perc_cis_calc <- function(df) { perc_cis <- (abs(df[,7])/(abs(df[,6] - df[,7]) + abs(df[,7]))) * 100; cbind(df,perc_cis) }

# execute functions:
## mel-mel
df2_MM <- df2_MM %>%
remove_less_20() %>%
parent_ratios() %>%
hybrid_ratios() %>%
perc_cis_calc
## sim-sech
df2_SS <- df2_SS %>%
remove_less_20() %>%
parent_ratios() %>%
hybrid_ratios() %>%
perc_cis_calc
## mel-sim
df2_MS <- df2_MS %>%
remove_less_20() %>%
parent_ratios() %>%
hybrid_ratios() %>%
perc_cis_calc

# join above with table from step 1
## first edit column name for gene to CG_ID so it matches
colnames(df2_MM)[1] <- "CG_ID"
colnames(df2_SS)[1] <- "CG_ID"
colnames(df2_MS)[1] <- "CG_ID"

df_final_MM <- left_join(df2_MM, df_cissub, by="CG_ID") %>% unique() %>% na.omit() %>% as.data.frame()
df_final_SS <- left_join(df2_SS, df_cissub, by="CG_ID") %>% unique() %>% na.omit() %>% as.data.frame()
df_final_MS <- left_join(df2_MS, df_cissub, by="CG_ID") %>% unique() %>% na.omit() %>% as.data.frame()

# do binomial tests on read counts to determine if differentialy expressed in parents and/or hybrids
## function didn't work so have to just write it for all three... ugh
pvalsHyb1 <- NULL;
pvalsHyb1_Par <- NULL;
pvalsHyb2 <- NULL;
pvalsHyb2_Par <- NULL;
pvalsPar1 <- NULL;
pvalsPar2 <- NULL;
pvalsHyb1_Hyb2 <- NULL;
pvalsHyb1_Bi <- NULL;
pvalsHyb2_Bi <- NULL;

for (i in 1:nrow(df_final_MM))
{
	Par1_Bi <- binom.test(df_final_MM[[i,2]], (df_final_MM[[i,2]]+df_final_MM[[i,3]]), p = 0.5, alternative = c("t"), conf.level = 0.95);
	Hyb1_Bi <- binom.test(df_final_MM[[i,4]], (df_final_MM[[i,4]]+df_final_MM[[i,5]]), p = 0.5, alternative = c("t"), conf.level = 0.95);
	# collect p-values from binomial tests
	pvalsHyb1 <- rbind(pvalsHyb1,Hyb1_Bi$p.value);
	pvalsPar1 <- rbind(pvalsPar1,Par1_Bi$p.value);
}
  #FDR correct pvalues
  pvalsHyb1.adj <- p.adjust(pvalsHyb1,method="fdr");
  pvalsPar1.adj <- p.adjust(pvalsPar1,method="fdr");

  #append on adjusted pvals
  df_final_MM <- cbind(df_final_MM,pvalsPar1.adj,pvalsHyb1.adj)

pvalsHyb1 <- NULL;
pvalsHyb1_Par <- NULL;
pvalsHyb2 <- NULL;
pvalsHyb2_Par <- NULL;
pvalsPar1 <- NULL;
pvalsPar2 <- NULL;
pvalsHyb1_Hyb2 <- NULL;
pvalsHyb1_Bi <- NULL;
pvalsHyb2_Bi <- NULL;


for (i in 1:nrow(df_final_MS))

{
	Par1_Bi <- binom.test(df_final_MS[[i,2]], (df_final_MS[[i,2]]+df_final_MS[[i,3]]), p = 0.5, alternative = c("t"), conf.level = 0.95);
	Hyb1_Bi <- binom.test(df_final_MS[[i,4]], (df_final_MS[[i,4]]+df_final_MS[[i,5]]), p = 0.5, alternative = c("t"), conf.level = 0.95);

	# collect p-values from binomial tests
	pvalsHyb1 <- rbind(pvalsHyb1,Hyb1_Bi$p.value);
  pvalsPar1 <- rbind(pvalsPar1,Par1_Bi$p.value);
}
  #FDR correct pvalues
  pvalsHyb1.adj <- p.adjust(pvalsHyb1,method="fdr");
  pvalsPar1.adj <- p.adjust(pvalsPar1,method="fdr");
  #append on adjusted pvals
  df_final_MS <- cbind(df_final_MS,pvalsPar1.adj,pvalsHyb1.adj)

pvalsHyb1 <- NULL;
pvalsHyb1_Par <- NULL;
pvalsHyb2 <- NULL;
pvalsHyb2_Par <- NULL;
pvalsPar1 <- NULL;
pvalsPar2 <- NULL;
pvalsHyb1_Hyb2 <- NULL;
pvalsHyb1_Bi <- NULL;
pvalsHyb2_Bi <- NULL;


  for (i in 1:nrow(df_final_SS))
{
  Par1_Bi <- binom.test(df_final_SS[[i,2]], (df_final_SS[[i,2]]+df_final_SS[[i,3]]), p = 0.5, alternative = c("t"), conf.level = 0.95);
  Hyb1_Bi <- binom.test(df_final_SS[[i,4]], (df_final_SS[[i,4]]+df_final_SS[[i,5]]), p = 0.5, alternative = c("t"), conf.level = 0.95);

  # collect p-values from binomial tests
  pvalsHyb1 <- rbind(pvalsHyb1,Hyb1_Bi$p.value);
  pvalsPar1 <- rbind(pvalsPar1,Par1_Bi$p.value);
}
  #FDR correct pvalues
  pvalsHyb1.adj <- p.adjust(pvalsHyb1,method="fdr");
  pvalsPar1.adj <- p.adjust(pvalsPar1,method="fdr");
  #append on adjusted pvals
  df_final_SS <- cbind(df_final_SS,pvalsPar1.adj,pvalsHyb1.adj)


# write tables
write.table(df_final_MM, file = "./Fledgling_analysis_ideas/Perc_cis_by_tissue_specificity/Data_tables_generated/Coolon_MM_gene_groups.txt", row.names = F, quote = FALSE, sep = "\t")
write.table(df_final_SS, file = "./Fledgling_analysis_ideas/Perc_cis_by_tissue_specificity/Data_tables_generated/Coolon_SS_gene_groups.txt", row.names = F, quote = FALSE, sep = "\t")
write.table(df_final_MS, file = "./Fledgling_analysis_ideas/Perc_cis_by_tissue_specificity/Data_tables_generated/Coolon_MS_gene_groups.txt", row.names = F, quote = FALSE, sep = "\t")

##################################################################################
#### STEP 3: Read in my wing disc data and do the same as Coolon whole body ######
##################################################################################
df3 <- read.delim("./Integrative_AS_genomics/AS_ATAC_RNA_2020_10_1/RNA_seq/Data_tables/ZHR_Z30_genic_counts_CPM_final_dm6.txt", header = T)

# subset to only get relevant columns
df3 <- df3[,c(1,2,8,14,20)]

# execute functions from above:
## mel-mel
### round counts up to nearest integer
df3_counts <- df3[,c(2:5)] %>% ceiling()
df3 <- cbind(df3$gene, df3_counts)
df3 <- df3 %>%
remove_less_20() %>%
parent_ratios() %>%
hybrid_ratios() %>%
perc_cis_calc
colnames(df3)[1] <- "CG_ID"

df3_final <- left_join(df3, df_cissub, by="CG_ID") %>% unique() %>% na.omit() %>% as.data.frame()


# do binomial tests on read counts to determine if differentialy expressed in parents and/or hybrids
## function didn't work so have to just write it for all three... ugh
pvalsHyb1 <- NULL;
pvalsHyb1_Par <- NULL;
pvalsHyb2 <- NULL;
pvalsHyb2_Par <- NULL;
pvalsPar1 <- NULL;
pvalsPar2 <- NULL;
pvalsHyb1_Hyb2 <- NULL;
pvalsHyb1_Bi <- NULL;
pvalsHyb2_Bi <- NULL;

for (i in 1:nrow(df3_final))
{
	Par1_Bi <- binom.test(df3_final[[i,2]], (df3_final[[i,2]]+df3_final[[i,3]]), p = 0.5, alternative = c("t"), conf.level = 0.95);
	Hyb1_Bi <- binom.test(df3_final[[i,4]], (df3_final[[i,4]]+df3_final[[i,5]]), p = 0.5, alternative = c("t"), conf.level = 0.95);
	# collect p-values from binomial tests
	pvalsHyb1 <- rbind(pvalsHyb1,Hyb1_Bi$p.value);
	pvalsPar1 <- rbind(pvalsPar1,Par1_Bi$p.value);
}
  #FDR correct pvalues
  pvalsHyb1.adj <- p.adjust(pvalsHyb1,method="fdr");
  pvalsPar1.adj <- p.adjust(pvalsPar1,method="fdr");

  #append on adjusted pvals
  df3_final <- cbind(df3_final,pvalsPar1.adj,pvalsHyb1.adj)

# specify if it is imaginal disc tissue specific
## get max tissue from master df
df <- read.delim("./DDivergence/Fly_duplicate_evol/Existing_Datasets/Fly_gene_ages_tau_G_17_TableS1.txt", header = TRUE)

# calc max tissue and append to our final df for wing disc
df$max_tissue <- colnames(df[,c(8:(ncol(df)))])[apply(df[,c(8:(ncol(df)))],1,which.max)]
df <- df[,c(4,ncol(df))]
df3_final <- left_join(df3_final, df, by = "CG_ID")

# categorize as wing_TS if TS and max tissue == imaginal disc
for (i in 1:nrow(df3_final)) {
if (df3_final$tau[i] > 0.61 & (df3_final$max_tissue[i] == "L3_Imaginal_Discs" | df3_final$max_tissue[i] == "L3_Imaginal_Discs.1" | df3_final$max_tissue[i] == "L3_Imaginal_Discs.2" | df3_final$max_tissue[i] == "L3_Imaginal_Discs.3"| df3_final$max_tissue[i] == "L3_Imaginal_Discs.4")){
	df3_final$class[i] <- "TS v"
}
}


#write table
write.table(df3_final, file = "./Fledgling_analysis_ideas/Perc_cis_by_tissue_specificity/Data_tables_generated/Ertl_MM_wing_disc_gene_groups.txt", row.names = F, quote = FALSE, sep = "\t")

##################################################################################
########################## STEP 4: PLOT PLOT PLOT ################################
##################################################################################

# re-read in datasets generated above and make equivalent column names
whole_body_MM <- read.delim("./Fledgling_analysis_ideas/Perc_cis_by_tissue_specificity/Data_tables_generated/Coolon_MM_gene_groups.txt", header = T)
colnames(whole_body_MM)[2:5] <- c("P1", "P2", "HYB1", "HYB2")
whole_body_SS <- read.delim("./Fledgling_analysis_ideas/Perc_cis_by_tissue_specificity/Data_tables_generated/Coolon_SS_gene_groups.txt", header = T)
colnames(whole_body_SS)[2:5] <- c("P1", "P2", "HYB1", "HYB2")
whole_body_MS <- read.delim("./Fledgling_analysis_ideas/Perc_cis_by_tissue_specificity/Data_tables_generated/Coolon_MS_gene_groups.txt", header = T)
colnames(whole_body_MS)[2:5] <- c("P1", "P2", "HYB1", "HYB2")
wing_disc_MM <- read.delim("./Fledgling_analysis_ideas/Perc_cis_by_tissue_specificity/Data_tables_generated/Ertl_MM_wing_disc_gene_groups.txt", header = T)
colnames(wing_disc_MM)[2:5] <- c("P1", "P2", "HYB1", "HYB2")

## functions to generate same plots for each group
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

# expression level across groups using parent 1
P1_expr_level_by_group <- function(df){
  df[df$Age_classification != "mel-only" & df$Age_classification != "mel-group" & df$Age_classification != "mel-subgroup" & df$Age_classification != "mel-complex",] %>%
    ggplot(aes(x=class, y=log(P1), fill=class)) +
    geom_boxplot(notch=TRUE) +
    scale_x_discrete(labels=c("Const. active", "Developmentally-\nregulated:\nmany tissues", "Tiss. specific")) +
    xlab("") +
    ylab("log10 Parent1 Expression level)") +
    scale_fill_discrete(guide=FALSE)
    theme_main()
}
expr_diff_by_group <- function(df){
  df[df$Age_classification != "mel-only" & df$Age_classification != "mel-group" & df$Age_classification != "mel-subgroup" & df$Age_classification != "mel-complex",] %>%
    ggplot(aes(x=class, y=abs(log2_P_ratio), fill=class)) +
    geom_boxplot(notch=TRUE) +
    scale_x_discrete(labels=c("Const. active", "Developmentally-\nregulated:\nmany tissues", "Tiss. specific")) +
    xlab("") +
    ylab("log2 P1/P2 Expression divergence))") +
    scale_fill_discrete(guide=FALSE) +
    theme_main()
}
perc_cis_by_group <- function(df){
  df[(df$pvalsPar1.adj < 0.05 | df$pvalsHyb1.adj < 0.05) & df$Age_classification != "mel-only" & df$Age_classification != "mel-group" & df$Age_classification != "mel-subgroup" & df$Age_classification != "mel-complex",] %>%
    ggplot(aes(x=class, y=perc_cis, fill=class)) +
    geom_boxplot(notch=TRUE) +
    scale_x_discrete(labels=c("Const. active", "Developmentally-\nregulated:\nmany tissues", "Tiss. specific")) +
    xlab("") +
    ylab("% cis contribution to expression divergence") +
    scale_fill_discrete(guide=FALSE) +
    theme_main()
}

# mel-mel whole body
MM1 <- P1_expr_level_by_group(whole_body_MM) + ggtitle("mel-mel, expression level")
MM2 <- expr_diff_by_group(whole_body_MM) + ggtitle("mel-mel, expression divergence")
MM3 <- perc_cis_by_group(whole_body_MM) + ggtitle("mel-mel, percent cis")

# sim-sech whole body
SS1 <- P1_expr_level_by_group(whole_body_SS) + ggtitle("sim-sech, expression level")
SS2 <- expr_diff_by_group(whole_body_SS) + ggtitle("sim-sech, expression divergence")
SS3 <- perc_cis_by_group(whole_body_SS) + ggtitle("sim-sech, percent cis")

# mel-sim whole body
MS1 <- P1_expr_level_by_group(whole_body_MS) + ggtitle("mel-sim, expression level")
MS2 <- expr_diff_by_group(whole_body_MS) + ggtitle("mel-sim, expression divergence")
MS3 <- perc_cis_by_group(whole_body_MS) + ggtitle("mel-sim, percent cis")

# mel-mel wing disc
P1_expr_level_by_group_wing_discs <- function(df){df[df$class != "TS" & df$Age_classification != "mel-only" & df$Age_classification != "mel-group" & df$Age_classification != "mel-subgroup" & df$Age_classification != "mel-complex",] %>%
  ggplot(aes(x=class, y=log(P1), fill=class)) +
  geom_boxplot(notch=TRUE) +
  scale_x_discrete(labels=c("Const. active", "Developmentally-\nregulated:\nmany tissues", "Wing disc specific")) +
  xlab("") +
  ylab("log10 Parent1 Expression level)") +
  scale_fill_discrete(guide=FALSE) +
  theme_main()
}
expr_diff_by_group_wing_discs <- function(df){ df[df$class != "TS" & df$Age_classification != "mel-only" & df$Age_classification != "mel-group" & df$Age_classification != "mel-subgroup" & df$Age_classification != "mel-complex",] %>%
  ggplot(aes(x=class, y=abs(log2_P_ratio), fill=class)) +
  geom_boxplot(notch=TRUE) +
  scale_x_discrete(labels=c("Const. active", "Developmentally-\nregulated:\nmany tissues", "Wing disc specific")) +
  xlab("") +
  ylab("log2 P1/P2 Expression divergence") +
  scale_fill_discrete(guide=FALSE) +
  theme_main()
}
perc_cis_by_group_wing_discs <- function(df){ df[(df$pvalsPar1.adj < 0.05 | df$pvalsHyb1.adj < 0.05) & df$class != "TS" & df$Age_classification != "mel-only" & df$Age_classification != "mel-group" & df$Age_classification != "mel-subgroup" & df$Age_classification != "mel-complex",] %>%
  ggplot(aes(x=class, y=abs(perc_cis), fill=class)) +
  geom_boxplot(notch=TRUE) +
  scale_x_discrete(labels=c("Const. active", "Developmentally-\nregulated:\nmany tissues", "Wing disc specific")) +
  xlab("") +
  ylab("% cis contribution to expression divergence") +
  scale_fill_discrete(guide=FALSE) +
  theme_main()
}

MM4 <- P1_expr_level_by_group_wing_discs(wing_disc_MM) + ggtitle("mel-mel wing disc, expression level")
MM5 <- expr_diff_by_group_wing_discs(wing_disc_MM) + ggtitle("mel-mel wing disc, expression divergence")
MM6 <- perc_cis_by_group_wing_discs(wing_disc_MM) + ggtitle("mel-mel wing disc, percent cis")

# save figs
ggsave(MM1, file = "./Fledgling_analysis_ideas/Perc_cis_by_tissue_specificity/Figures/Expr_level_by_group_MM.pdf")
ggsave(MM2, file = "./Fledgling_analysis_ideas/Perc_cis_by_tissue_specificity/Figures/Expr_diff_by_group_MM.pdf")
ggsave(MM3, file = "./Fledgling_analysis_ideas/Perc_cis_by_tissue_specificity/Figures/Perc_cis_by_group_MM.pdf")

ggsave(SS1, file = "./Fledgling_analysis_ideas/Perc_cis_by_tissue_specificity/Figures/Expr_level_by_group_SS.pdf")
ggsave(SS2, file = "./Fledgling_analysis_ideas/Perc_cis_by_tissue_specificity/Figures/Expr_diff_by_group_SS.pdf")
ggsave(SS3, file = "./Fledgling_analysis_ideas/Perc_cis_by_tissue_specificity/Figures/Perc_cis_by_group_SS.pdf")

ggsave(MS1, file = "./Fledgling_analysis_ideas/Perc_cis_by_tissue_specificity/Figures/Expr_level_by_group_MS.pdf")
ggsave(MS2, file = "./Fledgling_analysis_ideas/Perc_cis_by_tissue_specificity/Figures/Expr_diff_by_group_MS.pdf")
ggsave(MS3, file = "./Fledgling_analysis_ideas/Perc_cis_by_tissue_specificity/Figures/Perc_cis_by_group_MS.pdf")

ggsave(MM4, file = "./Fledgling_analysis_ideas/Perc_cis_by_tissue_specificity/Figures/Expr_level_by_group_MM_WD.pdf")
ggsave(MM5, file = "./Fledgling_analysis_ideas/Perc_cis_by_tissue_specificity/Figures/Expr_diff_by_group_MM_WD.pdf")
ggsave(MM6, file = "./Fledgling_analysis_ideas/Perc_cis_by_tissue_specificity/Figures/Perc_cis_by_group_MM_WD.pdf")


# now look at metrics acorss comparisons
## need to join all the columns to do that, so will add suffix to each and then join
colnames(whole_body_MM) <- paste(colnames(whole_body_MM), "MM", sep = "_")
colnames(whole_body_MM)[1] <- "CG_ID"

colnames(whole_body_SS) <- paste(colnames(whole_body_SS), "SS", sep = "_")
colnames(whole_body_SS)[1] <- "CG_ID"

colnames(whole_body_MS) <- paste(colnames(whole_body_MS), "MS", sep = "_")
colnames(whole_body_MS)[1] <- "CG_ID"

whole_body_ALL <- join_all(list(whole_body_MM, whole_body_SS, whole_body_MS), by = "CG_ID", type = "full") %>% na.omit() %>% unique()

perc_cis_plotting <- melt(whole_body_ALL, measure.vars = c("perc_cis_MM", "perc_cis_SS", "perc_cis_MS"))

# plot all combined
ALL_per_cis_ALL <- perc_cis_plotting[perc_cis_plotting$pvalsPar1.adj_MM < 0.05 | perc_cis_plotting$pvalsPar1.adj_SS < 0.05 | perc_cis_plotting$pvalsPar1.adj_MS < 0.05,] %>%
ggplot(aes(x=variable, y=value, fill=variable)) +
geom_boxplot(notch=TRUE) +
#facet_wrap(~class_MM, labeller = as_labeller(c('non_TS_HK' = "Housekeeping", 'non_TS_non_HK' = "Developmentally-\nregulated: many tissues", 'TS' = "Tissue specific"))) +
scale_x_discrete(labels=c("Mel-Mel", "Sech-Sim", "Mel-Sim")) +
xlab("Divergence time between comparisons: shortest -> longest") +
ylab("% cis contribution to expression divergence") +
scale_fill_discrete(guide=FALSE) +
theme_main()
ggsave(ALL_per_cis_ALL, file = "./Fledgling_analysis_ideas/Perc_cis_by_tissue_specificity/Figures/perc_cis_ALL_wholebody_ALL.pdf")

# plot each by group
ALL_per_cis <- perc_cis_plotting[perc_cis_plotting$pvalsPar1.adj_MM < 0.05 | perc_cis_plotting$pvalsPar1.adj_SS < 0.05 | perc_cis_plotting$pvalsPar1.adj_MS < 0.05,] %>%
ggplot(aes(x=variable, y=value, fill=class_MM, alpha = variable)) +
geom_boxplot(notch=TRUE) +
facet_wrap(~class_MM, labeller = as_labeller(c('non_TS_HK' = "Const. active", 'non_TS_non_HK' = "Developmentally-\nregulated: many tissues", 'TS' = "Tissue specific"))) +
scale_x_discrete(labels=c("Mel-Mel", "Sech-Sim", "Mel-Sim")) +
xlab("Divergence time between comparisons: shortest -> longest") +
ylab("% cis contribution to expression divergence") +
scale_fill_discrete(guide=FALSE) +
theme_main() +
scale_alpha_discrete(guide=FALSE)
ggsave(ALL_per_cis, file = "./Fledgling_analysis_ideas/Perc_cis_by_tissue_specificity/Figures/perc_cis_ALL_wholebody.pdf", width = 10)
