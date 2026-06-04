########################################
## Packages - don't package shame people
########################################

# .libPaths(c(
# 	"/data/san/data2/users/david/miniforge/envs/R/lib/R/library",
# 	"/data/san/data0/users/david/rstudio/packages", .libPaths()))
library(data.table)
library(dplyr)
library(ggside)
library(ggsci)
library(ggpubr)
setwd("/data/san/data2/users/david/2024_snp_investigation/plot_ellen_GFP")

df <- fread("results.tsv", col.names = c("sample", "au"))

# set factor levels for sample
df$sample <- factor(df$sample, levels = c("WT", "ER3", "ER6"))

# boxplot of au by sample
plot <- ggplot(df, aes(x = sample, y = au)) +
	geom_boxplot(aes(fill = sample)) +
	scale_fill_manual(values = c("#3C5488", "#8899bc", "#b4cbf9")) +
	theme_bw() +
	xlab(element_blank()) +
	ylab("Mean GFP fluorescence (a.u.)") +
	stat_compare_means(method = "t.test", 
		comparisons = list(c("WT", "ER3"), 
		c("WT", "ER6")), label = "p.format", size = 3) +
	theme(legend.position = "none") +
	coord_cartesian(clip = "off", ylim = c(0, 65))

ggsave("GFP_fluorescence_boxplot.png", 
	plot, unit="cm", width = 6, height = 8)

# top 10 samples by 