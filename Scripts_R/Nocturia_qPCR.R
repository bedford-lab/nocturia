require(tidyverse)
require(lme4)
require(lmerTest)

# Set working directory and plot directory
setwd("/Volumes/bedfordlab/Nocturia/Data_Analysis")
plotdir = paste(getwd(), "/Plots/", sep = "") 

# Read in the qPCR data:
df<- read.csv("/Volumes/bedfordlab/Nocturia/Data_Analysis/Nocturia Lookup Tables - qPCR.csv")
df<- read.csv("/Volumes/bedfordlab/Nocturia/Data_Analysis/qPCR_Bruns/FoldChange_data_long_format.csv")
str(df)
range(df$FoldChange)

# Make new grouping variables:
df<- df %>% mutate(AgeTime = paste(Age, Time),
                   TissueGene = paste(Tissue, Gene))

df$AgeTime<- factor(df$AgeTime, levels = c("Young AM", "Young PM",
                                           "Old AM", "Old PM"))

# Spot check fold change values
hist(df$FoldChange, breaks = 20)
hist(log10(df$FoldChange), breaks = 20)
hist(log1p(df$FoldChange), breaks = 20)

# Plot fold change for each gene
tissueGene<- unique(df$TissueGene)

for (i in tissueGene){
  fo <- df %>%
    filter(TissueGene == i)
  
  ggplot(data = fo, aes(x = AgeTime, y = log10(FoldChange))) +
    facet_grid(.~ Sex) +
    geom_boxplot() +
    #scale_y_continuous(trans = 'log10') +
    geom_point(colour="black", pch=21, position=position_jitter(0.1), size=1) +
    ggtitle(i) +
    theme_classic()
  ggsave(paste0(plotdir, "qPCR/", i, ".pdf"), width = 20, height = 10, units = "cm", useDingbats = FALSE)
  
}

# Plot gene specific histograms:
for (i in tissueGene){
  fo <- df %>%
    filter(TissueGene == i)
  
  p = ggplot(data = fo, aes(FoldChange)) +
    geom_histogram() +
    ggtitle(i)
  
  print(p)
  
}

# Remove fold change values higher than 15
df = df %>% filter(FoldChange < 20) %>% droplevels()

# Run stats:
library(broom)

stats<- df %>% filter(!TissueGene == "kidney Cry") %>% group_by(TissueGene) %>%
  do(broom::tidy(lm(log10(FoldChange) ~ Age*Time + Sex, .))) %>% 
  filter(!term == "(Intercept)") %>%
  mutate(P_round = round(p.value, digits = 4)) %>%
  ungroup()
write.csv(stats, "Nocturia_qPCR_stats_log10.csv")

# Sample size: 
ss = df %>% select(Age, Sex, Mouse_ID) %>% unique()
length(unique(ss$Mouse_ID))
table(ss$Sex, ss$Age)

## BLADDER MODELS:

# Bladder Bmal1:
summary(lm(log10(FoldChange) ~ Age + Time, data = df %>% filter(Sex == "F", TissueGene == "Bladder Bmal1"))) # NS
summary(lm(log10(FoldChange) ~ Time, data = df %>% filter(Sex == "F", Age == "Young", TissueGene == "Bladder Bmal1"))) # NS
summary(lm(log10(FoldChange) ~ Time, data = df %>% filter(Sex == "F", Age == "Old", TissueGene == "Bladder Bmal1"))) # NS

summary(lm(log10(FoldChange) ~ Age*Time, data = df %>% filter(Sex == "M", TissueGene == "Bladder Bmal1"))) # Age*Time: 0.0358 *
summary(lm(log10(FoldChange) ~ Time, data = df %>% filter(Sex == "M", Age == "Young", TissueGene == "Bladder Bmal1"))) # 0.0119 *
summary(lm(log10(FoldChange) ~ Time, data = df %>% filter(Sex == "M", Age == "Old", TissueGene == "Bladder Bmal1"))) # NS

# Bladder Per2:
summary(lm(log10(FoldChange) ~ Age + Time, data = df %>% filter(Sex == "F", TissueGene == "Bladder Per2"))) # NS
summary(lm(log10(FoldChange) ~ Time, data = df %>% filter(Sex == "F", Age == "Young", TissueGene == "Bladder Per2"))) # NS
summary(lm(log10(FoldChange) ~ Time, data = df %>% filter(Sex == "F", Age == "Old", TissueGene == "Bladder Per2"))) # NS

summary(lm(log10(FoldChange) ~ Age + Time, data = df %>% filter(Sex == "M", TissueGene == "Bladder Per2"))) # NS
summary(lm(log10(FoldChange) ~ Time, data = df %>% filter(Sex == "M", Age == "Young", TissueGene == "Bladder Per2"))) # NS
summary(lm(log10(FoldChange) ~ Time, data = df %>% filter(Sex == "M", Age == "Old", TissueGene == "Bladder Per2"))) # NS

# Bladder Piezo1:
summary(lm(log10(FoldChange) ~ Age*Time, data = df %>% filter(Sex == "F", TissueGene == "Bladder Piezo1"))) # NS
summary(lm(log10(FoldChange) ~ Time, data = df %>% filter(Sex == "F", Age == "Young", TissueGene == "Bladder Piezo1"))) # 0.0188 *
summary(lm(log10(FoldChange) ~ Time, data = df %>% filter(Sex == "F", Age == "Old", TissueGene == "Bladder Piezo1"))) # NS

summary(lm(log10(FoldChange) ~ Age + Time, data = df %>% filter(Sex == "M", TissueGene == "Bladder Piezo1"))) # Time: NS
summary(lm(log10(FoldChange) ~ Time, data = df %>% filter(Sex == "M", Age == "Young", TissueGene == "Bladder Piezo1"))) # NS
summary(lm(log10(FoldChange) ~ Time, data = df %>% filter(Sex == "M", Age == "Old", TissueGene == "Bladder Piezo1"))) # NS

## KIDNEY MODELS:

# kidney Bmal1:
summary(lm(log10(FoldChange) ~ Age*Time, data = df %>% filter(Sex == "F", TissueGene == "kidney Bmal1"))) # NS
summary(lm(log10(FoldChange) ~ Time, data = df %>% filter(Sex == "F", Age == "Young", TissueGene == "kidney Bmal1"))) # 0.00245 **
summary(lm(log10(FoldChange) ~ Time, data = df %>% filter(Sex == "F", Age == "Old", TissueGene == "kidney Bmal1"))) # 0.0773 .

summary(lm(log10(FoldChange) ~ Age*Time, data = df %>% filter(Sex == "M", TissueGene == "kidney Bmal1"))) # NS
summary(lm(log10(FoldChange) ~ Time, data = df %>% filter(Sex == "M", Age == "Young", TissueGene == "kidney Bmal1"))) # 0.0026 **
summary(lm(log10(FoldChange) ~ Time, data = df %>% filter(Sex == "M", Age == "Old", TissueGene == "kidney Bmal1"))) # NS

# kidney Per2:
summary(lm(log10(FoldChange) ~ Age*Time, data = df %>% filter(Sex == "F", TissueGene == "kidney Per2"))) # NS
summary(lm(log10(FoldChange) ~ Time, data = df %>% filter(Sex == "F", Age == "Young", TissueGene == "kidney Per2"))) # 0.0319 *
summary(lm(log10(FoldChange) ~ Time, data = df %>% filter(Sex == "F", Age == "Old", TissueGene == "kidney Per2"))) # 0.0533 .

summary(lm(log10(FoldChange) ~ Age*Time, data = df %>% filter(Sex == "M", TissueGene == "kidney Per2"))) # NS
summary(lm(log10(FoldChange) ~ Time, data = df %>% filter(Sex == "M", Age == "Young", TissueGene == "kidney Per2"))) # 6.98e-05 ***
summary(lm(log10(FoldChange) ~ Time, data = df %>% filter(Sex == "M", Age == "Old", TissueGene == "kidney Per2"))) # 0.0665 .

df<- df %>%
  mutate(Age = fct_relevel(Age, c("Young", "Old"))) %>%
  mutate(FoldChange_unTrans = FoldChange) %>%
  mutate(FoldChange = log1p(FoldChange)) # Express as log1p just for making plots

# What is the distribution of fold changes? Are there outliers?
hist(df$FoldChange_unTrans, breaks = 50)
hist(df$FoldChange, breaks = 50)

## BLADDER PLOTS:
library(Rmisc)

## Effect size table for Piezo1 bladder:
ES <- summarySE(data = df %>% filter(TissueGene == "Bladder Piezo1"),
                measurevar = "FoldChange_unTrans", groupvars = c("Sex", "Age", "Time"), na.rm=T)

# (dark - light) / light
round((ES[6,5] - ES[5,5]) / ES[5,5], digits = 3) 
round((ES[8,5] - ES[7,5]) / ES[7,5], digits = 3) 

## Bladder Bmal1:

df %>% filter(TissueGene == "Bladder Bmal1") %>% 
  dplyr::summarise(min = min(FoldChange, na.rm = TRUE),
                   max = max(FoldChange, na.rm = TRUE))

summaryF <- summarySE(data = df %>% filter(TissueGene == "Bladder Bmal1", Sex == "F"),
                      measurevar = "FoldChange", groupvars = c("Age", "Time"), na.rm=T)
ggplot(data=df %>% filter(TissueGene == "Bladder Bmal1", Sex == "F"), 
       aes(x=Time, y=FoldChange, fill = Time)) +
  facet_wrap(~Age) +
  geom_bar(data=summaryF, aes(fill=Time), stat="identity", width=0.85) +
  geom_errorbar(data=summaryF, aes(ymin=FoldChange-se, ymax=FoldChange+se, x=Time), size=0.5, width=0.25) + 
  geom_point(aes(fill=Time), colour="black", pch=21, position=position_jitter(width = 0.1), size=1) +
  scale_fill_manual(values=c("#f7941d", "#262262")) +
  #scale_y_continuous(trans = 'log10') +
  coord_cartesian(ylim = c(0, 1.689639)) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig4.Bladder Bmal1_FoldChange_F.pdf", width=5.5, height=6, units="cm", dpi=1500, useDingbats=FALSE)

summaryM <- summarySE(data = df %>% filter(TissueGene == "Bladder Bmal1", Sex == "M"),
                      measurevar = "FoldChange", groupvars = c("Age", "Time"), na.rm=T)
ggplot(data=df %>% filter(TissueGene == "Bladder Bmal1", Sex == "M"), 
       aes(x=Time, y=FoldChange, fill = Time)) +
  facet_wrap(~Age) +
  geom_bar(data=summaryM, aes(fill=Time), stat="identity", width=0.85) +
  geom_errorbar(data=summaryM, aes(ymin=FoldChange-se, ymax=FoldChange+se, x=Time), size=0.5, width=0.25) + 
  geom_point(aes(fill=Time), colour="black", pch=21, position=position_jitter(width = 0.1), size=1) +
  scale_fill_manual(values=c("#f7941d", "#262262")) +
  #scale_y_continuous(trans = 'log10') +
  coord_cartesian(ylim = c(0, 1.689639)) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig4.Bladder Bmal1_FoldChange_M.pdf", width=5.5, height=6, units="cm", dpi=1500, useDingbats=FALSE)

## Bladder Per2:

df %>% filter(TissueGene == "Bladder Per2") %>% 
  dplyr::summarise(min = min(FoldChange, na.rm = TRUE),
                   max = max(FoldChange, na.rm = TRUE))

summaryF <- summarySE(data = df %>% filter(TissueGene == "Bladder Per2", Sex == "F"),
                      measurevar = "FoldChange", groupvars = c("Age", "Time"), na.rm=T)
ggplot(data=df %>% filter(TissueGene == "Bladder Per2", Sex == "F"), 
       aes(x=Time, y=FoldChange, fill = Time)) +
  facet_wrap(~Age) +
  geom_bar(data=summaryF, aes(fill=Time), stat="identity", width=0.85) +
  geom_errorbar(data=summaryF, aes(ymin=FoldChange-se, ymax=FoldChange+se, x=Time), size=0.5, width=0.25) + 
  geom_point(aes(fill=Time), colour="black", pch=21, position=position_jitter(width = 0.1), size=1) +
  scale_fill_manual(values=c("#f7941d", "#262262")) +
  #scale_y_continuous(trans = 'log10') +
  coord_cartesian(ylim = c(0, 2.502843)) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig4.Bladder Per2_FoldChange_F.pdf", width=5.5, height=6, units="cm", dpi=1500, useDingbats=FALSE)

summaryM <- summarySE(data = df %>% filter(TissueGene == "Bladder Per2", Sex == "M"),
                      measurevar = "FoldChange", groupvars = c("Age", "Time"), na.rm=T)
ggplot(data=df %>% filter(TissueGene == "Bladder Per2", Sex == "M"), 
       aes(x=Time, y=FoldChange, fill = Time)) +
  facet_wrap(~Age) +
  geom_bar(data=summaryM, aes(fill=Time), stat="identity", width=0.85) +
  geom_errorbar(data=summaryM, aes(ymin=FoldChange-se, ymax=FoldChange+se, x=Time), size=0.5, width=0.25) + 
  geom_point(aes(fill=Time), colour="black", pch=21, position=position_jitter(width = 0.1), size=1) +
  scale_fill_manual(values=c("#f7941d", "#262262")) +
  #scale_y_continuous(trans = 'log10') +
  coord_cartesian(ylim = c(0, 2.502843)) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig4.Bladder Per2_FoldChange_M.pdf", width=5.5, height=6, units="cm", dpi=1500, useDingbats=FALSE)

## Bladder Piezo1:

df %>% filter(TissueGene == "Bladder Piezo1") %>% 
  dplyr::summarise(min = min(FoldChange, na.rm = TRUE),
                   max = max(FoldChange, na.rm = TRUE))

summaryF <- summarySE(data = df %>% filter(TissueGene == "Bladder Piezo1", Sex == "F"),
                      measurevar = "FoldChange", groupvars = c("Age", "Time"), na.rm=T)
ggplot(data=df %>% filter(TissueGene == "Bladder Piezo1", Sex == "F"), 
       aes(x=Time, y=FoldChange, fill = Time)) +
  facet_wrap(~Age) +
  geom_bar(data=summaryF, aes(fill=Time), stat="identity", width=0.85) +
  geom_errorbar(data=summaryF, aes(ymin=FoldChange-se, ymax=FoldChange+se, x=Time), size=0.5, width=0.25) + 
  geom_point(aes(fill=Time), colour="black", pch=21, position=position_jitter(width = 0.1), size=1) +
  scale_fill_manual(values=c("#f7941d", "#262262")) +
  #scale_y_continuous(trans = 'log10') +
  coord_cartesian(ylim = c(0, 2.725656)) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig4.Bladder Piezo1_FoldChange_F.pdf", width=5.5, height=6, units="cm", dpi=1500, useDingbats=FALSE)

summaryM <- summarySE(data = df %>% filter(TissueGene == "Bladder Piezo1", Sex == "M"),
                      measurevar = "FoldChange", groupvars = c("Age", "Time"), na.rm=T)
ggplot(data=df %>% filter(TissueGene == "Bladder Piezo1", Sex == "M"), 
       aes(x=Time, y=FoldChange, fill = Time)) +
  facet_wrap(~Age) +
  geom_bar(data=summaryM, aes(fill=Time), stat="identity", width=0.85) +
  geom_errorbar(data=summaryM, aes(ymin=FoldChange-se, ymax=FoldChange+se, x=Time), size=0.5, width=0.25) + 
  geom_point(aes(fill=Time), colour="black", pch=21, position=position_jitter(width = 0.1), size=1) +
  scale_fill_manual(values=c("#f7941d", "#262262")) +
  #scale_y_continuous(trans = 'log10') +
  coord_cartesian(ylim = c(0, 2.725656)) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig4.Bladder Piezo1_FoldChange_M.pdf", width=5.5, height=6, units="cm", dpi=1500, useDingbats=FALSE)

## KIDNEY PLOTS:
## kidney Bmal1:

df %>% filter(TissueGene == "kidney Bmal1") %>% 
  dplyr::summarise(min = min(FoldChange, na.rm = TRUE),
                   max = max(FoldChange, na.rm = TRUE))

summaryF <- summarySE(data = df %>% filter(TissueGene == "kidney Bmal1", Sex == "F"),
                      measurevar = "FoldChange", groupvars = c("Age", "Time"), na.rm=T)
ggplot(data=df %>% filter(TissueGene == "kidney Bmal1", Sex == "F"), 
       aes(x=Time, y=FoldChange, fill = Time)) +
  facet_wrap(~Age) +
  geom_bar(data=summaryF, aes(fill=Time), stat="identity", width=0.85) +
  geom_errorbar(data=summaryF, aes(ymin=FoldChange-se, ymax=FoldChange+se, x=Time), size=0.5, width=0.25) + 
  geom_point(aes(fill=Time), colour="black", pch=21, position=position_jitter(width = 0.1), size=1) +
  scale_fill_manual(values=c("#f7941d", "#262262")) +
  #scale_y_continuous(trans = 'log10') +
  coord_cartesian(ylim = c(0, 2.439539)) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig4.kidney Bmal1_FoldChange_F.pdf", width=5.5, height=6, units="cm", dpi=1500, useDingbats=FALSE)

summaryM <- summarySE(data = df %>% filter(TissueGene == "kidney Bmal1", Sex == "M"),
                      measurevar = "FoldChange", groupvars = c("Age", "Time"), na.rm=T)
ggplot(data=df %>% filter(TissueGene == "kidney Bmal1", Sex == "M"), 
       aes(x=Time, y=FoldChange, fill = Time)) +
  facet_wrap(~Age) +
  geom_bar(data=summaryM, aes(fill=Time), stat="identity", width=0.85) +
  geom_errorbar(data=summaryM, aes(ymin=FoldChange-se, ymax=FoldChange+se, x=Time), size=0.5, width=0.25) + 
  geom_point(aes(fill=Time), colour="black", pch=21, position=position_jitter(width = 0.1), size=1) +
  scale_fill_manual(values=c("#f7941d", "#262262")) +
  #scale_y_continuous(trans = 'log10') +
  coord_cartesian(ylim = c(0, 2.439539)) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig4.kidney Bmal1_FoldChange_M.pdf", width=5.5, height=6, units="cm", dpi=1500, useDingbats=FALSE)

## kidney Per2:

df %>% filter(TissueGene == "kidney Per2") %>% 
  dplyr::summarise(min = min(FoldChange, na.rm = TRUE),
                   max = max(FoldChange, na.rm = TRUE))

summaryF <- summarySE(data = df %>% filter(TissueGene == "kidney Per2", Sex == "F"),
                      measurevar = "FoldChange", groupvars = c("Age", "Time"), na.rm=T)
ggplot(data=df %>% filter(TissueGene == "kidney Per2", Sex == "F"), 
       aes(x=Time, y=FoldChange, fill = Time)) +
  facet_wrap(~Age) +
  geom_bar(data=summaryF, aes(fill=Time), stat="identity", width=0.85) +
  geom_errorbar(data=summaryF, aes(ymin=FoldChange-se, ymax=FoldChange+se, x=Time), size=0.5, width=0.25) + 
  geom_point(aes(fill=Time), colour="black", pch=21, position=position_jitter(width = 0.1), size=1) +
  scale_fill_manual(values=c("#f7941d", "#262262")) +
  #scale_y_continuous(trans = 'log10') +
  coord_cartesian(ylim = c(0, 2.945766)) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig4.kidney Per2_FoldChange_F.pdf", width=5.5, height=6, units="cm", dpi=1500, useDingbats=FALSE)

summaryM <- summarySE(data = df %>% filter(TissueGene == "kidney Per2", Sex == "M"),
                      measurevar = "FoldChange", groupvars = c("Age", "Time"), na.rm=T)
ggplot(data=df %>% filter(TissueGene == "kidney Per2", Sex == "M"), 
       aes(x=Time, y=FoldChange, fill = Time)) +
  facet_wrap(~Age) +
  geom_bar(data=summaryM, aes(fill=Time), stat="identity", width=0.85) +
  geom_errorbar(data=summaryM, aes(ymin=FoldChange-se, ymax=FoldChange+se, x=Time), size=0.5, width=0.25) + 
  geom_point(aes(fill=Time), colour="black", pch=21, position=position_jitter(width = 0.1), size=1) +
  scale_fill_manual(values=c("#f7941d", "#262262")) +
  #scale_y_continuous(trans = 'log10') +
  coord_cartesian(ylim = c(0, 2.945766)) +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(legend.position = "none", axis.text.y=element_blank(), axis.text.x=element_blank(), panel.spacing = unit(0, "lines"),
        strip.background = element_blank(), strip.text.x = element_blank())
ggsave("Plots/Fig4.kidney Per2_FoldChange_M.pdf", width=5.5, height=6, units="cm", dpi=1500, useDingbats=FALSE)