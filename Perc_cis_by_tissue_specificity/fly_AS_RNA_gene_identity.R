library(binom)
library(Hmisc)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(plyr)


theme_main <- function() {
  theme_bw() +
  theme(
  #panel.grid.major = element_blank(),
  #panel.grid.minor = element_blank(),
  axis.text = element_text(size = 15),
  axis.title = element_text(size = 15),
  strip.text = element_text(size = 15),
  legend.text= element_text(size = 15),
  legend.title = element_text(size = 10),
  plot.title = element_text(size = 15, face = "bold")

)
}



setwd("/Users/wittkopp_member/Code")

# read in dataset
df <- read.delim("./Tissue_environment_specificity_X_perc_cis/Fly_data/Existing_Datasets/Supplemental_dataset2.txt", header = T)

# only select hybrids from one direction
df <- df[,-c(6,7,12,13,18,19)]

# run stats
# Mel-Mel
pvalsHyb1 <- NULL;
pvalsHyb1_Par <- NULL;
pvalsHyb2 <- NULL;
pvalsHyb2_Par <- NULL;
pvalsPar1 <- NULL;
pvalsPar2 <- NULL;
pvalsHyb1_Hyb2 <- NULL;
pvalsHyb1_Bi <- NULL;
pvalsHyb2_Bi <- NULL;
tempH1_P <- NULL

for (i in 1:nrow(df))
{
	Par1_Bi <- binom.test(df[[i,2]], (df[[i,2]]+df[[i,3]]), p = 0.5, alternative = c("t"), conf.level = 0.95);
	Hyb1_Bi <- binom.test(df[[i,4]], (df[[i,4]]+df[[i,5]]), p = 0.5, alternative = c("t"), conf.level = 0.95);
  # collect p-values from binomial tests
	pvalsHyb1 <- rbind(pvalsHyb1,Hyb1_Bi$p.value);
	pvalsPar1 <- rbind(pvalsPar1,Par1_Bi$p.value);

  # FET tests between ratios
  hyb1_par.FET <- fisher.test(matrix(c(df[[i,4]],df[[i,5]],df[[i,2]],df[[i,3]]),nr=2));
  # collect p vals from FETs
  pvalsHyb1_Par <- rbind(pvalsHyb1_Par,hyb1_par.FET$p.value);
}

  #FDR correct pvalues
  pvalsHyb1.adj_MM <- p.adjust(pvalsHyb1,method="fdr");
  pvalsPar1.adj_MM <- p.adjust(pvalsPar1,method="fdr");
  pvalsHyb1_Par.adj_MM <- p.adjust(pvalsHyb1_Par,method="fdr");

  #append on adjusted pvals
  df <- cbind(df,pvalsPar1.adj_MM,pvalsHyb1.adj_MM,pvalsHyb1_Par.adj_MM)


# run stats
# Sim_sech
pvalsHyb1 <- NULL;
pvalsHyb1_Par <- NULL;
pvalsHyb2 <- NULL;
pvalsHyb2_Par <- NULL;
pvalsPar1 <- NULL;
pvalsPar2 <- NULL;
pvalsHyb1_Hyb2 <- NULL;
pvalsHyb1_Bi <- NULL;
pvalsHyb2_Bi <- NULL;

  for (i in 1:nrow(df))
  {
  	Par1_Bi <- binom.test(df[[i,6]], (df[[i,6]]+df[[i,7]]), p = 0.5, alternative = c("t"), conf.level = 0.95);
  	Hyb1_Bi <- binom.test(df[[i,8]], (df[[i,8]]+df[[i,9]]), p = 0.5, alternative = c("t"), conf.level = 0.95);
    # collect p-values from binomial tests
  	pvalsHyb1 <- rbind(pvalsHyb1,Hyb1_Bi$p.value);
  	pvalsPar1 <- rbind(pvalsPar1,Par1_Bi$p.value);

    # FET tests between ratios
    hyb1_par.FET <- fisher.test(matrix(c(df[[i,8]],df[[i,9]],df[[i,6]],df[[i,7]]),nr=2));
    # collect p vals from FETs
    pvalsHyb1_Par <- rbind(pvalsHyb1_Par,hyb1_par.FET$p.value);
  }

  #FDR correct pvalues
  pvalsHyb1.adj_SS <- p.adjust(pvalsHyb1,method="fdr");
  pvalsPar1.adj_SS <- p.adjust(pvalsPar1,method="fdr");
  pvalsHyb1_Par.adj_SS <- p.adjust(pvalsHyb1_Par,method="fdr");

  #append on adjusted pvals
  df <- cbind(df,pvalsPar1.adj_SS,pvalsHyb1.adj_SS,pvalsHyb1_Par.adj_SS)


# run stats
# Mel-Sim
pvalsHyb1 <- NULL;
pvalsHyb1_Par <- NULL;
pvalsHyb2 <- NULL;
pvalsHyb2_Par <- NULL;
pvalsPar1 <- NULL;
pvalsPar2 <- NULL;
pvalsHyb1_Hyb2 <- NULL;
pvalsHyb1_Bi <- NULL;
pvalsHyb2_Bi <- NULL;

  for (i in 1:nrow(df))
  {
  	Par1_Bi <- binom.test(df[[i,10]], (df[[i,10]]+df[[i,11]]), p = 0.5, alternative = c("t"), conf.level = 0.95);
  	Hyb1_Bi <- binom.test(df[[i,12]], (df[[i,12]]+df[[i,13]]), p = 0.5, alternative = c("t"), conf.level = 0.95);
    # collect p-values from binomial tests
  	pvalsHyb1 <- rbind(pvalsHyb1,Hyb1_Bi$p.value);
  	pvalsPar1 <- rbind(pvalsPar1,Par1_Bi$p.value);

    # FET tests between ratios
    hyb1_par.FET <- fisher.test(matrix(c(df[[i,12]],df[[i,13]],df[[i,10]],df[[i,11]]),nr=2));
    # collect p vals from FETs
    pvalsHyb1_Par <- rbind(pvalsHyb1_Par,hyb1_par.FET$p.value);
    }

  #FDR correct pvalues
  pvalsHyb1.adj_MS <- p.adjust(pvalsHyb1,method="fdr");
  pvalsPar1.adj_MS <- p.adjust(pvalsPar1,method="fdr");
  pvalsHyb1_Par.adj_MS <- p.adjust(pvalsHyb1_Par,method="fdr");

  #append on adjusted pvals
  df <- cbind(df,pvalsPar1.adj_MS,pvalsHyb1.adj_MS,pvalsHyb1_Par.adj_MS)

# calculate log2 ratios between parents and hybrid alleles
df$MM_P_ratio <- log2(df[,2]/df[,3])
df$MM_H_ratio <- log2(df[,4]/df[,5])
df$MM_PH_ratio <- df$MM_P_ratio/df$MM_H_ratio
df$SS_P_ratio <- log2(df[,6]/df[,7])
df$SS_H_ratio <- log2(df[,8]/df[,9])
df$SS_PH_ratio <- df$SS_P_ratio/df$SS_H_ratio
df$MS_P_ratio <- log2(df[,10]/df[,11])
df$MS_H_ratio <- log2(df[,12]/df[,13])
df$MS_PH_ratio <- df$MS_P_ratio/df$MS_H_ratio


# categorize cis, trans, etc.
## MM ##
##Set qvalue cut-off
critical_value <- 0.05
​
##Run classifier
df$Regulatory_class_MM <- "Ambiguous"
​
for (i in 1:nrow(df)) {
​
if (df$pvalsPar1.adj_MM[i] > critical_value & df$pvalsHyb1.adj_MM[i] > critical_value & df$pvalsHyb1_Par.adj_MM[i] > critical_value){
​
	df$Regulatory_class_MM[i] <- "Conserved"
​
} else if (df$pvalsPar1.adj_MM[i] < critical_value & df$pvalsHyb1.adj_MM[i] < critical_value & df$pvalsHyb1_Par.adj_MM[i] > critical_value){
​
	df$Regulatory_class_MM[i] <- "Cis"
​
} else if (df$pvalsPar1.adj_MM[i] < critical_value & df$pvalsHyb1.adj_MM[i] > critical_value & df$pvalsHyb1_Par.adj_MM[i] < critical_value){
​
	df$Regulatory_class_MM[i] <- "Trans"
​
} else if (df$pvalsPar1.adj_MM[i] < critical_value & df$pvalsHyb1.adj_MM[i] < critical_value & df$pvalsHyb1_Par.adj_MM[i] < critical_value & sign(df$MM_P_ratio[i]) == sign(df$MM_H_ratio[i]) & df$MM_PH_ratio[i] > 1){
​
	df$Regulatory_class_MM[i] <- "Cis_+_Trans,opposing"
​
} else if (df$pvalsPar1.adj_MM[i] < critical_value & df$pvalsHyb1.adj_MM[i] < critical_value & df$pvalsHyb1_Par.adj_MM[i] < critical_value & sign(df$MM_P_ratio[i]) == sign(df$MM_H_ratio[i]) & df$MM_PH_ratio[i] < 1){
​
  df$Regulatory_class_MM[i] <- "Cis_+_Trans,same"
​
} else if (df$pvalsPar1.adj_MM[i] < critical_value & df$pvalsHyb1.adj_MM[i] < critical_value & df$pvalsHyb1_Par.adj_MM[i] < critical_value & sign(df$MM_P_ratio[i]) != sign(df$MM_H_ratio[i])) {
​
	df$Regulatory_class_MM[i] <- "Cis_*_Trans"
​
} else if (df$pvalsPar1.adj_MM[i] > critical_value & df$pvalsHyb1.adj_MM[i] < critical_value & df$pvalsHyb1_Par.adj_MM[i] < critical_value){
​
	df$Regulatory_class_MM[i] <- "Compensatory"
}
}


##Run classifier for opposing and same
df$Direction_MM <- df$Regulatory_class_MM
​
for (i in 1:nrow(df)) {
​
if (df$Regulatory_class_MM[i] == "Cis_+_Trans,opposing"){
​
	df$Direction_MM[i] <- "Opposing"
​
} else if (df$Regulatory_class_MM[i] == "Cis_*_Trans"){
​
	df$Direction_MM[i] <- "Opposing"
​
} else if (df$Regulatory_class_MM[i] == "Compensatory"){
​
	df$Direction_MM[i] <- "Opposing"
​
} else if (df$Regulatory_class_MM[i] == "Cis_+_Trans,same"){
​
	df$Direction_MM[i] <- "Reinforcing"
}
}

## SS ##
##Run classifier
df$Regulatory_class_SS <- "Ambiguous"
​
for (i in 1:nrow(df)) {
​
if (df$pvalsPar1.adj_SS[i] > critical_value & df$pvalsHyb1.adj_SS[i] > critical_value & df$pvalsHyb1_Par.adj_SS[i] > critical_value){
​
	df$Regulatory_class_SS[i] <- "Conserved"
​
} else if (df$pvalsPar1.adj_SS[i] < critical_value & df$pvalsHyb1.adj_SS[i] < critical_value & df$pvalsHyb1_Par.adj_SS[i] > critical_value){
​
	df$Regulatory_class_SS[i] <- "Cis"
​
} else if (df$pvalsPar1.adj_SS[i] < critical_value & df$pvalsHyb1.adj_SS[i] > critical_value & df$pvalsHyb1_Par.adj_SS[i] < critical_value){
​
	df$Regulatory_class_SS[i] <- "Trans"
​
} else if (df$pvalsPar1.adj_SS[i] < critical_value & df$pvalsHyb1.adj_SS[i] < critical_value & df$pvalsHyb1_Par.adj_SS[i] < critical_value & sign(df$SS_P_ratio[i]) == sign(df$SS_H_ratio[i]) & df$SS_PH_ratio[i] > 1){
​
	df$Regulatory_class_SS[i] <- "Cis_+_Trans,opposing"
​
} else if (df$pvalsPar1.adj_SS[i] < critical_value & df$pvalsHyb1.adj_SS[i] < critical_value & df$pvalsHyb1_Par.adj_SS[i] < critical_value & sign(df$SS_P_ratio[i]) == sign(df$SS_H_ratio[i]) & df$SS_PH_ratio[i] < 1){
​
  df$Regulatory_class_SS[i] <- "Cis_+_Trans,same"
​
} else if (df$pvalsPar1.adj_SS[i] < critical_value & df$pvalsHyb1.adj_SS[i] < critical_value & df$pvalsHyb1_Par.adj_SS[i] < critical_value & sign(df$SS_P_ratio[i]) != sign(df$SS_H_ratio[i])) {
​
	df$Regulatory_class_SS[i] <- "Cis_*_Trans"
​
} else if (df$pvalsPar1.adj_SS[i] > critical_value & df$pvalsHyb1.adj_SS[i] < critical_value & df$pvalsHyb1_Par.adj_SS[i] < critical_value){
​
	df$Regulatory_class_SS[i] <- "Compensatory"
}
}


##Run classifier for opposing and same
df$Direction_SS <- df$Regulatory_class_SS
​
for (i in 1:nrow(df)) {
​
if (df$Regulatory_class_SS[i] == "Cis_+_Trans,opposing"){
​
	df$Direction_SS[i] <- "Opposing"
​
} else if (df$Regulatory_class_SS[i] == "Cis_*_Trans"){
​
	df$Direction_SS[i] <- "Opposing"
​
} else if (df$Regulatory_class_SS[i] == "Compensatory"){
​
	df$Direction_SS[i] <- "Opposing"
​
} else if (df$Regulatory_class_SS[i] == "Cis_+_Trans,same"){
​
	df$Direction_SS[i] <- "Reinforcing"
}
}
​
## MS ##
##Run classifier
df$Regulatory_class_MS <- "Ambiguous"
​
for (i in 1:nrow(df)) {
​
if (df$pvalsPar1.adj_MS[i] > critical_value & df$pvalsHyb1.adj_MS[i] > critical_value & df$pvalsHyb1_Par.adj_MS[i] > critical_value){
​
	df$Regulatory_class_MS[i] <- "Conserved"
​
} else if (df$pvalsPar1.adj_MS[i] < critical_value & df$pvalsHyb1.adj_MS[i] < critical_value & df$pvalsHyb1_Par.adj_MS[i] > critical_value){
​
	df$Regulatory_class_MS[i] <- "Cis"
​
} else if (df$pvalsPar1.adj_MS[i] < critical_value & df$pvalsHyb1.adj_MS[i] > critical_value & df$pvalsHyb1_Par.adj_MS[i] < critical_value){
​
	df$Regulatory_class_MS[i] <- "Trans"
​
} else if (df$pvalsPar1.adj_MS[i] < critical_value & df$pvalsHyb1.adj_MS[i] < critical_value & df$pvalsHyb1_Par.adj_MS[i] < critical_value & sign(df$MS_P_ratio[i]) == sign(df$MS_H_ratio[i]) & df$MS_PH_ratio[i] > 1){
​
	df$Regulatory_class_MS[i] <- "Cis_+_Trans,opposing"
​
} else if (df$pvalsPar1.adj_MS[i] < critical_value & df$pvalsHyb1.adj_MS[i] < critical_value & df$pvalsHyb1_Par.adj_MS[i] < critical_value & sign(df$MS_P_ratio[i]) == sign(df$MS_H_ratio[i]) & df$MS_PH_ratio[i] < 1){
​
  df$Regulatory_class_MS[i] <- "Cis_+_Trans,same"
​
} else if (df$pvalsPar1.adj_MS[i] < critical_value & df$pvalsHyb1.adj_MS[i] < critical_value & df$pvalsHyb1_Par.adj_MS[i] < critical_value & sign(df$MS_P_ratio[i]) != sign(df$MS_H_ratio[i])) {
​
	df$Regulatory_class_MS[i] <- "Cis_*_Trans"
​
} else if (df$pvalsPar1.adj_MS[i] > critical_value & df$pvalsHyb1.adj_MS[i] < critical_value & df$pvalsHyb1_Par.adj_MS[i] < critical_value){
​
	df$Regulatory_class_MS[i] <- "Compensatory"
}
}


##Run classifier for opposing and same
df$Direction_MS <- df$Regulatory_class_MS
​
for (i in 1:nrow(df)) {
​
if (df$Regulatory_class_MS[i] == "Cis_+_Trans,opposing"){
​
	df$Direction_MS[i] <- "Opposing"
​
} else if (df$Regulatory_class_MS[i] == "Cis_*_Trans"){
​
	df$Direction_MS[i] <- "Opposing"
​
} else if (df$Regulatory_class_MS[i] == "Compensatory"){
​
	df$Direction_MS[i] <- "Opposing"
​
} else if (df$Regulatory_class_MS[i] == "Cis_+_Trans,same"){
​
	df$Direction_MS[i] <- "Reinforcing"
}
}
​
df_final <- df[df$Regulatory_class_MM != "Ambiguous" & df$Regulatory_class_SS != "Ambiguous" & df$Regulatory_class_MS != "Ambiguous",]

# append on class suffix for each species comparison
df_final$Direction_MM <- paste(df_final$Direction_MM, "MM", sep = "_")
df_final$Direction_SS <- paste(df_final$Direction_SS, "SS", sep = "_")
df_final$Direction_MS <- paste(df_final$Direction_MS, "MS", sep = "_")



df_final_sankey_MM_SS <- df_final[,c(33, 35)]
df_final_sankey_MM_SS <- ddply(df_final_sankey_MM_SS,.(Direction_MM, Direction_SS),nrow)
colnames(df_final_sankey_MM_SS) <- c("source", "target","value")

df_final_sankey_SS_MS <- df_final[,c(35, 37)]
df_final_sankey_SS_MS <- ddply(df_final_sankey_SS_MS,.(Direction_SS, Direction_MS),nrow)
colnames(df_final_sankey_SS_MS) <- c("source", "target","value")

links <- rbind(df_final_sankey_MM_SS,df_final_sankey_SS_MS)

nodes <- data.frame(
  name=c(as.character(links$source),
  as.character(links$target)) %>% unique()
)

links$IDsource <- match(links$source, nodes$name)-1
links$IDtarget <- match(links$target, nodes$name)-1


p <- sankeyNetwork(Links = links, Nodes = nodes,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name",
              sinksRight=FALSE)
