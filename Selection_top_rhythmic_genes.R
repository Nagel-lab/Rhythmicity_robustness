#################################
# Selection of top rhythmic genes
#################################


# Define directory
setwd("C:/Users/tbonnot/Documents/")

# Load packages
library(tidyverse)
library(ggrepel)
library(gridExtra)
library(eulerr)
library(reshape2)


# Load datasets
#--------------

# Rhythmic genes identified in Romanowski et al. (2020)
Romanowski_cycling <- read.table("./Projects/Side_projects/Rhythmicity_significance/Romanowski_JTK_results.txt", header = T)

# Rhythmic genes identified in Bonnot and Nagel (2021)
Bonnot_cycling <- read.table("./Projects/Side_projects/Rhythmicity_significance/Bonnot_Nagel_rhythmic_transcriptome.txt", header = T)

# BBX genes
BBX <- read.table("./Projects/Side_projects/Rhythmicity_significance/BBX_family.txt", header = T)

# Selected clock genes, according to the Figure 1 in Creux and Harmer (2019)
clock_genes <- read.table("./Projects/Side_projects/Rhythmicity_significance/selected_clock_genes.txt", header = T)

# TFs and the TF family they belong to, according to Pruneda-Paz et al (2014)
Families <- read.table("./Projects/Side_projects/Rhythmicity_significance/Arabidopsis TFs.txt", header = T, sep = "\t")

# ChIP-Seq results (for clock genes, using CAST-R)
ChIP <- read.table("./Projects/Side_projects/Rhythmicity_significance/ChIP_clock_53_cycling_genes.txt", header = T)
names(ChIP)[1] <- "Gene"
ChIP <- left_join(ChIP, clock_genes, by = "Gene")

# DAP-Seq data
dap1 <- read.csv2("./Resources/DAP_Seq_Arabidopsis/Dap_Seq_1.csv", header = T, sep = ",")
dap2 <- read.csv2("./Resources/DAP_Seq_Arabidopsis/Dap_Seq_2.csv", header = T, sep = ",")
dap3 <- read.csv2("./Resources/DAP_Seq_Arabidopsis/Dap_Seq_3.csv", header = T, sep = ",")

# Expression data from Romanowski et al. (2020)
Romanowski_data <- read.table("./Projects/Side_projects/Rhythmicity_significance/Romanowski_expression_data.txt", header = T)

# TAIR annotation
annot <- read.csv2("./Projects/Side_projects/Rhythmicity_significance/TAIR10_gene_description.txt", header = T, sep = "\t")
annot <- annot[!duplicated(annot$AGI),]



# Plot the expression of selected genes in Romanowski data: genes with very high and very low Pvalues 
#---------------------------------------------------------

# Select genes of interest
Gene.plot <- c("AT3G09600", "AT1G72310")

# Filter expression data for these genes
Romanowski_data.selected <- filter(Romanowski_data, Gene_ID %in% Gene.plot)

# Reformat the table
Romanowski_data.selected <- melt(Romanowski_data.selected, id.vars = "Gene_ID")
Romanowski_data.selected <- separate(Romanowski_data.selected, variable, into = c("Genotype", "Light", "Time", "Rep"))

# Calculate mean and sd
str(Romanowski_data.selected)
Romanowski_data.selected$value <- as.numeric(Romanowski_data.selected$value)
Romanowski_data.selected$Time <- as.numeric(Romanowski_data.selected$Time)
Romanowski.selected.mean <- Romanowski_data.selected %>% 
  group_by(Gene_ID, Time) %>%
  summarise(
    mean = mean(value, na.rm = T),
    SD = sd(value, na.rm = T)
  )

# Plot the selected genes
ggplot(Romanowski.selected.mean, aes(x = Time, y = mean, color = Gene_ID)) +
  geom_rect(aes(xmin=36, xmax=48, ymin=-Inf, ymax=Inf), 
            fill="grey", alpha=0.06, color = NA)+
  geom_rect(aes(xmin=60, xmax=72, ymin=-Inf, ymax=Inf), 
            fill="grey", alpha=0.06, color = NA)+
  geom_point(size = 2.5)+
  geom_line(linewidth = 1)+
  geom_errorbar(aes(ymin=mean-SD, ymax=mean+SD),width=1.2, size=0.6,
                position=position_dodge(width=0.1),show.legend=F)+
  scale_color_manual(values = c("#12A78D", "#086757"))+
  facet_wrap(~Gene_ID, scales = "free_y")+
  theme_bw()+
  xlab("Time of day (ZT)")+
  ylab("Normalized transcript level")+
  theme(panel.grid=element_blank(), axis.ticks.length=unit(-0.15, "cm"),
        axis.text=element_text(size=12, colour = "black"),axis.text.x=element_text(margin = unit(c(1, 0, 1, 0), "mm")),
        axis.text.y=element_text(margin=unit(c(0,1,0,0),"mm")),
        axis.title = element_text(size = 12),
        strip.background = element_blank(),
        legend.position = "none")



# Transform adjusted P-values of the rhythmicity significance to get positive values: -log10(P-value)
#---------------------------------------------------------------------------------
Romanowski_cycling$log_BH.Q <- -log10(Romanowski_cycling$BH.Q)
Bonnot_cycling <- Bonnot_cycling[,c("AGI","meta2d_BH.Q")]
Bonnot_cycling$log_BH.Q <- -log10(Bonnot_cycling$meta2d_BH.Q)

# In Bonnot data, replace Inf values by 13 (WARNING: mention this in the figure legend!!)
#Bonnot_cycling[sapply(Bonnot_cycling, is.infinite)] <- 13.0


# Select TF families of interest from the two datasets
#-----------------------------------------------------
# All rhythmic genes
Romanowski_all <- Romanowski_cycling[,c(1,3)]
Romanowski_all$Group <- "All rhythmic genes"
Bonnot_all <- Bonnot_cycling[,c(1,3)]
Bonnot_all$Group <- "All rhythmic genes"

# TFs
Romanowski_TF <- filter(Romanowski_all, AGI %in% Families$AGI)
Romanowski_TF$Group <- "All rhythmic TFs"
Bonnot_TF <- filter(Bonnot_all, AGI %in% Families$AGI)
Bonnot_TF$Group <- "All rhythmic TFs"

# MYB
Romanowski_MYB <- filter(Romanowski_all, AGI %in% Families[Families$Family == "MYB",1])
Romanowski_MYB$Group <- "MYB"
Bonnot_MYB <- filter(Bonnot_all, AGI %in% Families[Families$Family == "MYB",1])
Bonnot_MYB$Group <- "MYB"

# NAC
Romanowski_NAC <- filter(Romanowski_all, AGI %in% Families[Families$Family == "NAC",1])
Romanowski_NAC$Group <- "NAC"
Bonnot_NAC <- filter(Bonnot_all, AGI %in% Families[Families$Family == "NAC",1])
Bonnot_NAC$Group <- "NAC"

# AP2-EREBP
Romanowski_AP2_EREBP <- filter(Romanowski_all, AGI %in% Families[Families$Family == "AP2-EREBP",1])
Romanowski_AP2_EREBP$Group <- "AP2-EREBP"
Bonnot_AP2_EREBP <- filter(Bonnot_all, AGI %in% Families[Families$Family == "AP2-EREBP",1])
Bonnot_AP2_EREBP$Group <- "AP2-EREBP"

# bZIP
Romanowski_bZIP <- filter(Romanowski_all, AGI %in% Families[Families$Family == "bZIP",1])
Romanowski_bZIP$Group <- "BZIP"
Bonnot_bZIP <- filter(Bonnot_all, AGI %in% Families[Families$Family == "bZIP",1])
Bonnot_bZIP$Group <- "BZIP"

# bHLH
Romanowski_BHLH <- filter(Romanowski_all, AGI %in% Families[Families$Family == "bHLH",1])
Romanowski_BHLH$Group <- "BHLH"
Bonnot_BHLH <- filter(Bonnot_all, AGI %in% Families[Families$Family == "bHLH",1])
Bonnot_BHLH$Group <- "BHLH"

# BBX
Romanowski_BBX <- filter(Romanowski_all, AGI %in% BBX$AGI)
Romanowski_BBX$Group <- "BBX"
Bonnot_BBX <- filter(Bonnot_all, AGI %in% BBX$AGI)
Bonnot_BBX$Group <- "BBX"

# Clock genes
Romanowski_clock <- filter(Romanowski_all, AGI %in% clock_genes$AGI)
Romanowski_clock$Group <- "Clock genes"
Bonnot_clock <- filter(Bonnot_all, AGI %in% clock_genes$AGI)
Bonnot_clock$Group <- "Clock genes"

# Group dataframes
Romanowski_selected <- rbind(Romanowski_clock, Romanowski_TF, Romanowski_MYB,
                             Romanowski_NAC, Romanowski_AP2_EREBP, Romanowski_bZIP,
                             Romanowski_BHLH, Romanowski_BBX)
Bonnot_selected <- rbind(Bonnot_clock, Bonnot_TF, Bonnot_MYB,
                         Bonnot_NAC, Bonnot_AP2_EREBP, Bonnot_bZIP,
                         Bonnot_BHLH, Bonnot_BBX)

Romanowski_selected$Study <- "Romanowski"
Bonnot_selected$Study <- "Bonnot"
Data_boxplot <- rbind(Romanowski_selected, Bonnot_selected)


# Plot the distribution of transformed P-values for the different sets of genes
#------------------------------------------------------------------------------

Data_boxplot$Group <- factor(Data_boxplot$Group, levels = c("All rhythmic TFs", 
                                                            "Clock genes", 
                                                            "AP2-EREBP",
                                                            "BZIP",
                                                            "MYB",
                                                            "NAC",
                                                            "BHLH",
                                                            "BBX"))
Data_boxplot$Study <- factor(Data_boxplot$Study, levels = c("Romanowski", "Bonnot"))

# Make the plot (here we have removed outliers)
g2 <- ggplot(Data_boxplot, aes(x = Group, y = log_BH.Q))+
  geom_jitter(color='#7c8ce5ff',alpha=0.5, width = 0.3, size = 2)+
  geom_boxplot(fill="#7c8ce5ff",color="#2237aaff",alpha=0.2, width = 0.8, size = 0.6, outlier.alpha = 0)+ 
  facet_wrap(~Study, nrow = 2, scales = "free_y")+
  guides(color='none')+
  theme_bw()+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.length = unit(-0.1, "cm"))+
  xlab("")+
  ylab("-log10(Adjusted P-value)")



# Select the top 100 rhythmic TFs in the two studies and get the overlap
#-----------------------------------------------------------------------

# Select the top 100 rhythmic TFs in each study
Romanowski_top100 <- Romanowski_TF %>% slice_max(log_BH.Q, n = 100)
Bonnot_top100 <- Bonnot_TF %>% slice_max(log_BH.Q, n = 100)

# Venn diagram to see the overlapping genes
top100 <- list(A = Bonnot_top100$AGI, B = Romanowski_top100$AGI)
g3 <- plot(euler(top100, shape = "ellipse"), quantities = TRUE)

# Summarize the info of the top rhythmic TFs
Top_rhythmic_TFs <- full_join(Romanowski_top100[,1:2], Bonnot_top100[,1:2], by = "AGI")

# Add gene annotation
Top_rhythmic_TFs <- left_join(Top_rhythmic_TFs, annot[,c("AGI","Description")], by = "AGI")

# Add the TF family info
Top_rhythmic_TFs <- left_join(Top_rhythmic_TFs, Families, by = "AGI")

# Add the BBX gene info
Top_rhythmic_TFs <- left_join(Top_rhythmic_TFs, BBX, by = "AGI")

# Add the clock gene info
Top_rhythmic_TFs <- left_join(Top_rhythmic_TFs, clock_genes, by = "AGI")

# Export the table
#write.table(Top_rhythmic_TFs, "./Projects/Side_projects/Rhythmicity_significance/Top_100_rhythmic_TFs.txt", row.names = F, sep = "\t")


# Select the top 100 rhythmic genes (that are not TFs) in the two studies and get the overlap
#--------------------------------------------------------------------------------------------

# Exclude TFs from rhythmic genes in each study
Romanowski_nonTF <- anti_join(Romanowski_all, Families, by = "AGI")
Bonnot_nonTF <- anti_join(Bonnot_all, Families, by = "AGI")

# Select the top 100 rhythmic non-TF genes in each study
Romanowski_nonTF_top100 <- Romanowski_nonTF[,c("AGI","log_BH.Q")] %>% slice_max(log_BH.Q, n = 100)
Bonnot_nonTF_top100 <- Bonnot_nonTF[,c("AGI","log_BH.Q")] %>% slice_max(log_BH.Q, n = 100)

# Venn diagram to see the overlapping genes
top100_nonTF <- list(A = Bonnot_nonTF_top100$AGI, B = Romanowski_nonTF_top100$AGI)
plot(euler(top100_nonTF, shape = "ellipse"), quantities = TRUE)

# Summarize the info of the top rhythmic non-TF genes
Top_rhythmic_nonTF <- full_join(Romanowski_nonTF_top100, Bonnot_nonTF_top100, by = "AGI")

# Add gene annotation
Top_rhythmic_nonTF <- left_join(Top_rhythmic_nonTF, annot[,c("AGI","Description")], by = "AGI")

# Add the clock gene info
Top_rhythmic_nonTF <- left_join(Top_rhythmic_nonTF, clock_genes, by = "AGI")

# Export the table
#write.table(Top_rhythmic_nonTF, "./Projects/Side_projects/Rhythmicity_significance/Top_100_rhythmic_nonTFs.txt", row.names = F, sep = "\t")



# Investigate interactions between the 53 top rhythmic TFs (using ChIP-Seq and DAP-Seq data)
#-------------------------------------------------------------------------------------------

# Define the list of 53 
Overlap_genes <- intersect(Romanowski_top100$AGI, Bonnot_top100$AGI)

# Filter for the 53 genes in the ChIP data
names(ChIP) <- c("Clock","to","from")
ChIP <- filter(ChIP, from%in% Overlap_genes)
ChIP <- filter(ChIP, to%in% Overlap_genes)

# Filter DAP-Seq data
dap1 <- filter(dap1, TF_ID%in% Overlap_genes)
dap1 <- filter(dap1, Target_ID%in% Overlap_genes)
dap2 <- filter(dap2, TF_ID%in% Overlap_genes)
dap2 <- filter(dap2, Target_ID%in% Overlap_genes)
dap3 <- filter(dap3, TF_ID%in% Overlap_genes)
dap3 <- filter(dap3, Target_ID%in% Overlap_genes)

dap <- rbind(dap1, dap2, dap3)
unique(c(dap$from, dap$to))
names(dap) <- c("from","to")

# Combine ChIP-Seq and DAP-Seq data
ChIP <- ChIP[,c(3,2)]
ChIP$Method <- "ChIP"
dap$Method <- "DAP-Seq"
interaction_53_genes <- rbind(ChIP,dap)

# Export results to build a network in cytoscape
#write.table(interaction_53_genes, "./Projects/Side_projects/Rhythmicity_significance/interaction_53_cycling_genes.txt", row.names = F, sep = "\t")



