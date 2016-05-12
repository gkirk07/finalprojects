setwd("~/Desktop/")
library(readr)
library(tidyr)
library(ggplot2)
library(cowplot)

ETV6.data <- read_csv("condensed_ETV6data.csv")

#BTLA Analysis

#ETV6/BTLA analysis
ETV6.deleted <- subset(ETV6.data, select = c("ETV6_alteration",
                                             "B_cell_pathway_lesion_description",
                                             "Other_key_lesions")) %>%
  #for some reason || does not work here and returns all 2013 rows of the spreadsheet but | works
  filter(ETV6_alteration == "Yes" |
           ETV6_alteration == "Mutation" |
           ETV6_alteration == "Deletion / Rearrangement" |
           ETV6_alteration == "Rearrangement")

ETV6.wt <- subset(ETV6.data, select = c("ETV6_alteration",
                                        "B_cell_pathway_lesion_description",
                                        "Other_key_lesions")) %>%
  filter(ETV6_alteration == "No")

BTLA.deleted.ETV6.deleted <- ETV6.deleted[grep("BTLA", ETV6.deleted$Other_key_lesions), ]
BTLA.deleted.ETV6.wt <- ETV6.wt[grep("BTLA", ETV6.wt$Other_key_lesions), ]

Percent.BTLA.deleted.ETV6.deleted <- (nrow(BTLA.deleted.ETV6.deleted) / nrow(ETV6.deleted))*100
Percent.BTLA.deleted.ETV6.wt <- (nrow(BTLA.deleted.ETV6.wt) / nrow(ETV6.wt))*100

#PAX5/BTLA analysis
PAX5.deleted <- subset(ETV6.data, select = c("PAX5_alteration", "B_cell_pathway_lesion_description", "Other_key_lesions")) %>%
  filter(PAX5_alteration == "Yes")
PAX5.wt <- subset(ETV6.data, select = c("PAX5_alteration", "B_cell_pathway_lesion_description", "Other_key_lesions")) %>%
  filter(PAX5_alteration == "No")

BTLA.deleted.PAX5.deleted <- PAX5.deleted[grep("BTLA", PAX5.deleted$Other_key_lesions), ]
BTLA.deleted.PAX5.wt <- PAX5.wt[grep("BTLA", PAX5.wt$Other_key_lesions), ]

Percent.BTLA.deleted.PAX5.deleted <- (nrow(BTLA.deleted.PAX5.deleted) / nrow(PAX5.deleted))*100
Percent.BTLA.deleted.PAX5.wt <- (nrow(BTLA.deleted.PAX5.wt) / nrow(PAX5.wt))*100

#EBF1/BTLA analysis
EBF1.deleted <- subset(ETV6.data, select = c("EBF1_deletion", "B_cell_pathway_lesion_description", "Other_key_lesions")) %>%
  filter(EBF1_deletion == "Yes")
EBF1.wt <- subset(ETV6.data, select = c("EBF1_deletion", "B_cell_pathway_lesion_description", "Other_key_lesions")) %>%
  filter(EBF1_deletion == "No")

BTLA.deleted.EBF1.deleted <- EBF1.deleted[grep("BTLA", EBF1.deleted$Other_key_lesions), ]
BTLA.deleted.EBF1.wt <- EBF1.wt[grep("BTLA", EBF1.wt$Other_key_lesions), ]

Percent.BTLA.deleted.EBF1.deleted <- (nrow(BTLA.deleted.EBF1.deleted) / nrow(EBF1.deleted))*100
Percent.BTLA.deleted.EBF1.wt <- (nrow(BTLA.deleted.EBF1.wt) / nrow(EBF1.wt))*100

#Bar graph in ggplot
Del <- c(Percent.BTLA.deleted.ETV6.deleted, Percent.BTLA.deleted.PAX5.deleted, Percent.BTLA.deleted.EBF1.deleted)
WT <- c(Percent.BTLA.deleted.ETV6.wt, Percent.BTLA.deleted.PAX5.wt, Percent.BTLA.deleted.EBF1.wt)

untidy = data.frame(WT, Del)
tidy <- gather(untidy, status, Percent_BTLA_Deleted)
gene <- c("ETV6", "PAX5", "EBF1")
df <- data.frame(gene, tidy)
ggplot(df, aes(gene, Percent_BTLA_Deleted, fill=status)) +
  geom_bar(color="black", width=0.8, stat="identity", position = "dodge") +
  scale_x_discrete(limits=c("ETV6", "PAX5", "EBF1")) +
  scale_fill_brewer(palette = "Greens", name="Gene\nStatus", labels=c("Altered", "WT")) +
  theme_cowplot() + ylab("Percent BTLA-CD200 Deleted") + xlab("") +
  theme(text=element_text(size = 20), axis.text.x=element_text(size = 20))


#CDKN Analysis

#ETV6/CDKN analysis
CDKN.deleted.ETV6.deleted <- ETV6.deleted[grep("CDKN", ETV6.deleted$Other_key_lesions), ]
CDKN.deleted.ETV6.wt <- ETV6.wt[grep("CDKN", ETV6.wt$Other_key_lesions), ]

Percent.CDKN.deleted.ETV6.deleted <- (nrow(CDKN.deleted.ETV6.deleted) / nrow(ETV6.deleted))*100
Percent.CDKN.deleted.ETV6.wt <- (nrow(CDKN.deleted.ETV6.wt) / nrow(ETV6.wt))*100

#PAX5/CDKN analysis
CDKN.deleted.PAX5.deleted <- PAX5.deleted[grep("CDKN", PAX5.deleted$Other_key_lesions), ]
CDKN.deleted.PAX5.wt <- PAX5.wt[grep("CDKN", PAX5.wt$Other_key_lesions), ]

Percent.CDKN.deleted.PAX5.deleted <- (nrow(CDKN.deleted.PAX5.deleted) / nrow(PAX5.deleted))*100
Percent.CDKN.deleted.PAX5.wt <- (nrow(CDKN.deleted.PAX5.wt) / nrow(PAX5.wt))*100

#EBF1/CDKN analysis
CDKN.deleted.EBF1.deleted <- EBF1.deleted[grep("CDKN", EBF1.deleted$Other_key_lesions), ]
CDKN.deleted.EBF1.wt <- EBF1.wt[grep("CDKN", EBF1.wt$Other_key_lesions), ]

Percent.CDKN.deleted.EBF1.deleted <- (nrow(CDKN.deleted.EBF1.deleted) / nrow(EBF1.deleted))*100
Percent.CDKN.deleted.EBF1.wt <- (nrow(CDKN.deleted.EBF1.wt) / nrow(EBF1.wt))*100

#CDKN Bar graph in ggplot
Del2 <- c(Percent.CDKN.deleted.ETV6.deleted, Percent.CDKN.deleted.PAX5.deleted, Percent.CDKN.deleted.EBF1.deleted)
WT2 <- c(Percent.CDKN.deleted.ETV6.wt, Percent.CDKN.deleted.PAX5.wt, Percent.CDKN.deleted.EBF1.wt)

untidy2 = data.frame(WT2, Del2)
tidy2 <- gather(untidy2, status, Percent_CDKN_Deleted)
df2 <- data.frame(gene, tidy2)
ggplot(df2, aes(gene, Percent_CDKN_Deleted, fill=status)) +
  geom_bar(color="black", width=0.8, stat="identity", position = "dodge") +
  scale_x_discrete(limits=c("ETV6", "PAX5", "EBF1")) +
  scale_fill_brewer(palette ="Paired", name="Gene\nStatus", labels=c("Altered", "WT")) +
  theme_cowplot() + ylab("Percent CDKN2A/B Deleted") + xlab("") +
  theme(text=element_text(size = 20), axis.text.x=element_text(size = 20))

#BTG1 Analysis

#ETV6/BTG1 analysis
BTG1.deleted.ETV6.deleted <- ETV6.deleted[grep("BTG1", ETV6.deleted$Other_key_lesions), ]
BTG1.deleted.ETV6.wt <- ETV6.wt[grep("BTG1", ETV6.wt$Other_key_lesions), ]

Percent.BTG1.deleted.ETV6.deleted <- (nrow(BTG1.deleted.ETV6.deleted) / nrow(ETV6.deleted))*100
Percent.BTG1.deleted.ETV6.wt <- (nrow(BTG1.deleted.ETV6.wt) / nrow(ETV6.wt))*100

#PAX5/BTG1 analysis
BTG1.deleted.PAX5.deleted <- PAX5.deleted[grep("BTG1", PAX5.deleted$Other_key_lesions), ]
BTG1.deleted.PAX5.wt <- PAX5.wt[grep("BTG1", PAX5.wt$Other_key_lesions), ]

Percent.BTG1.deleted.PAX5.deleted <- (nrow(BTG1.deleted.PAX5.deleted) / nrow(PAX5.deleted))*100
Percent.BTG1.deleted.PAX5.wt <- (nrow(BTG1.deleted.PAX5.wt) / nrow(PAX5.wt))*100

#EBF1/BTG1 analysis
BTG1.deleted.EBF1.deleted <- EBF1.deleted[grep("BTG1", EBF1.deleted$Other_key_lesions), ]
BTG1.deleted.EBF1.wt <- EBF1.wt[grep("BTG1", EBF1.wt$Other_key_lesions), ]

Percent.BTG1.deleted.EBF1.deleted <- (nrow(BTG1.deleted.EBF1.deleted) / nrow(EBF1.deleted))*100
Percent.BTG1.deleted.EBF1.wt <- (nrow(BTG1.deleted.EBF1.wt) / nrow(EBF1.wt))*100

#BTG1 Bar graph in ggplot
Del3 <- c(Percent.BTG1.deleted.ETV6.deleted, Percent.BTG1.deleted.PAX5.deleted, Percent.BTG1.deleted.EBF1.deleted)
WT3 <- c(Percent.BTG1.deleted.ETV6.wt, Percent.BTG1.deleted.PAX5.wt, Percent.BTG1.deleted.EBF1.wt)

untidy3 = data.frame(WT3, Del3)
tidy3 <- gather(untidy3, status, Percent_BTG1_Deleted)
df3 <- data.frame(gene, tidy3)
ggplot(df3, aes(gene, Percent_BTG1_Deleted, fill=status)) +
  geom_bar(color="black", width=0.8, stat="identity", position = "dodge") +
  scale_x_discrete(limits=c("ETV6", "PAX5", "EBF1")) +
  scale_fill_brewer(palette = "Reds", name="Gene\nStatus", labels=c("Altered", "WT")) +
  theme_cowplot() + ylab("Percent BTG1 Deleted") + xlab("") +
  theme(text=element_text(size = 20), axis.text.x=element_text(size = 20))
