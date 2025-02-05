#############################################################
# Robustness of rhythmicity vs amplitude and expression level
#############################################################

setwd("C:/Users/tbonnot/Documents/")

# Load packages
library(reshape2)
library(tidyverse)
library(ggalign)
library(ggpointdensity)
library(viridis)
library(ggrepel)


# Load datasets
#--------------

# Rhythmic genes identified in Romanowski et al. (2020)
Romanowski_cycling <- read.table("./Projects/Side_projects/Rhythmicity_significance/Romanowski_BHQ_Amp.txt", header = T)

# Rythmic genes identified in Bonnot and Nagel (2021)
Bonnot_cycling <- read.table("./Projects/Side_projects/Rhythmicity_significance/Bonnot_Nagel_rhythmic_transcriptome.txt", header = T)
Bonnot_cycling <- Bonnot_cycling[,c("AGI","meta2d_BH.Q","JTK_amplitude")]
names(Bonnot_cycling) <- c("AGI","BH.Q","AMP")

# Romanowski's expression data
data_romanowski <- read.table("./Projects/Side_projects/Rhythmicity_significance/Romanowski_expression_data.txt", header = T)

# Bonnot' expression data
data_bonnot <- read.table("./Projects/Side_projects/Rhythmicity_significance/rlog_values_100.txt", header = T)
data_bonnot <- data_bonnot[,grepl("AGI", colnames(data_bonnot))|grepl("TOT",colnames(data_bonnot))] # keep transcriptome data
data_bonnot <- data_bonnot[,!grepl("H", colnames(data_bonnot))] # keep control condition

# TFs and the TF family they belong to, according to Pruneda-Paz et al (2014)
Families <- read.table("./Projects/Side_projects/Rhythmicity_significance/Arabidopsis TFs.txt", header = T, sep = "\t")
TFs <- Families[!duplicated(Families$AGI),1]

# Clock genes
clock_genes <- read.table("./Projects/Side_projects/Rhythmicity_significance/selected_clock_genes.txt", header = T)
names(clock_genes)[2] <- "Gene_ID"


# Select rhythmic TFs in datasets
#--------------------------------
Romanowski_cycling_TF <- filter(Romanowski_cycling, Gene_ID %in% TFs)
Bonnot_cycling_TF <- filter(Bonnot_cycling, AGI %in% TFs)

data_romanowski_cycling_TF <- filter(data_romanowski, Gene_ID %in% Romanowski_cycling_TF$Gene_ID)
data_bonnot_cycling_TF <- filter(data_bonnot, AGI %in% Bonnot_cycling_TF$AGI)
  

# Calculate average expression by time of day, then by gene
#----------------------------------------------------------
data_romanowski_cycling_TF <- melt(data_romanowski_cycling_TF, id.vars = "Gene_ID")
data_romanowski_cycling_TF <- separate(data_romanowski_cycling_TF, variable, into = c("Genotype","Condition","Time","Rep"))
data_bonnot_cycling_TF <- melt(data_bonnot_cycling_TF, id.vars = "AGI")
data_bonnot_cycling_TF <- separate(data_bonnot_cycling_TF, variable, into = c("RNA","Time","Rep"))

str(data_romanowski_cycling_TF)
data_romanowski_cycling_TF$value <- as.numeric(data_romanowski_cycling_TF$value)
romanowski.cycTF.mean <- data_romanowski_cycling_TF[,c("Gene_ID","Time","value")] %>%
  group_by(Gene_ID, Time) %>% summarise(Mean = mean(value))
romanowski.cycTF.mean <- romanowski.cycTF.mean %>%
  group_by(Gene_ID) %>% summarise(Mean = mean(Mean))

str(data_bonnot_cycling_TF)
bonnot.cycTF.mean <- data_bonnot_cycling_TF[,c("AGI","Time","value")] %>%
  group_by(AGI, Time) %>% summarise(Mean = mean(value))
bonnot.cycTF.mean <- bonnot.cycTF.mean %>%
  group_by(AGI) %>% summarise(Mean = mean(Mean))


# Calculate -log10(BH.Q)
#-----------------------

# Merge BH.Q, amplitude and mean expression
Romanowski_cycling_TF <- left_join(Romanowski_cycling_TF, romanowski.cycTF.mean, by = "Gene_ID")
Bonnot_cycling_TF <- left_join(Bonnot_cycling_TF, bonnot.cycTF.mean, by = "AGI")

# -log10(BH.Q)
Romanowski_cycling_TF$log10BHQ <- -log10(Romanowski_cycling_TF$BH.Q)
Bonnot_cycling_TF$log10BHQ <- -log10(Bonnot_cycling_TF$BH.Q)

# In Bonnot and Nagel dataset, replace Inf values by 13 
Bonnot_cycling_TF[sapply(Bonnot_cycling_TF, is.infinite)] <- 13.0


# Identify the -log10(BH.Q) cutoff value corresponding to the top 100 rhythmic TFs
#---------------------------------------------------------------------------
Romanowski_top100 <- Romanowski_cycling_TF %>% slice_max(log10BHQ, n = 100) 
Bonnot_top100 <- Bonnot_cycling_TF %>% slice_max(log10BHQ, n = 100)

Romanowski_min_value_top100 <- min(Romanowski_top100$log10BHQ)
Bonnot_min_value_top100 <- min(Bonnot_top100$log10BHQ)
log_BH.Q_cutoff <- data.frame(Dataset = c("Romanowski","Bonnot"),
                              Cutoff = c(Romanowski_min_value_top100, Bonnot_min_value_top100))


# Group datasets
#---------------
names(Bonnot_cycling_TF)[1] <- "Gene_ID"
Bonnot_cycling_TF$Dataset <- "Bonnot"
Romanowski_cycling_TF$Dataset <- "Romanowski"
cycling_TF <- rbind(Romanowski_cycling_TF, Bonnot_cycling_TF)


# Add the clock gene info
#------------------------
cycling_TF <- left_join(cycling_TF, clock_genes, by = "Gene_ID")


# Make the plots
#---------------

cycling_TF$Gene <- as.factor(cycling_TF$Gene)
cycling_TF$Dataset <- factor(cycling_TF$Dataset, levels = c("Romanowski","Bonnot"))
log_BH.Q_cutoff$Dataset <- factor(log_BH.Q_cutoff$Dataset, levels = c("Romanowski","Bonnot"))

labels <- c("Bonnot" = "Bonnot and Nagel (2021)",
            "Romanowski" = "Romanowski et al. (2020)")


amp <- ggplot(cycling_TF, aes(x = log10BHQ, y = AMP, label = Gene))+
  geom_pointdensity()+
  scale_color_viridis(name = "Point density")+
  geom_smooth(alpha=0.25, color="black", fill="black")+
  facet_wrap(~Dataset, scales = "free",
             labeller = as_labeller(labels))+
  geom_vline(log_BH.Q_cutoff, mapping = aes(xintercept = Cutoff), 
             color = "purple", linetype = "dashed", linewidth = 0.5)+
  theme_bw()+
  theme(axis.text = element_text(color = "black", size = 12),
        panel.grid = element_blank(),
        axis.ticks.length=unit(-0.15, "cm"),
        axis.text.x=element_text(margin = unit(c(1, 0, 0, 0), "mm")),
        axis.text.y=element_text(margin=unit(c(0,1,0,0),"mm")),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12))+
  labs(x = "-log10(Adjusted P-value)",
       y = "Amplitude of oscillation")+
  geom_label_repel(data = cycling_TF[!is.na(cycling_TF$Gene),],
                   size = 2,
                   box.padding = 0.3,
                   max.overlaps = 50,
                   segment.size = 0.5, 
                   force = 50,
                   min.segment.length = 0,
                   point.padding = unit(0.1, "lines"))


exp <- ggplot(cycling_TF, aes(x = log10BHQ, y = Mean, label = Gene))+
  geom_pointdensity()+
  scale_color_viridis(name = "Point Density")+
  geom_smooth(alpha=0.25, color="black", fill="black")+
  facet_wrap(~Dataset, scales = "free",
             labeller = as_labeller(labels))+
  geom_vline(log_BH.Q_cutoff, mapping = aes(xintercept = Cutoff), 
             color = "purple", linetype = "dashed", linewidth = 0.5)+
  theme_bw()+
  theme(axis.text = element_text(color = "black", size = 12),
        panel.grid = element_blank(),
        axis.ticks.length=unit(-0.15, "cm"),
        axis.text.x=element_text(margin = unit(c(1, 0, 0, 0), "mm")),
        axis.text.y=element_text(margin=unit(c(0,1,0,0),"mm")),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12))+
  labs(x = "-log10(Adjusted P-value)",
       y = "Average expression level")+
  geom_label_repel(data = cycling_TF[!is.na(cycling_TF$Gene),],
                  size = 2,
                  box.padding = 0.3,
                  max.overlaps = 50,
                  segment.size = 0.5, 
                  force = 50,
                  min.segment.length = 0,
                  point.padding = unit(0.1, "lines"))

align_plots(amp, exp, ncol = 1)


# Calculate rank of the BH.Q values within each Dataset
#------------------------------------------------------
cycling_TF <- cycling_TF %>%
  group_by(Dataset) %>%
  mutate(Rank = dense_rank(desc(log10BHQ))) %>%
  ungroup()


