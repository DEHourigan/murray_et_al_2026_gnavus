########################################
## Packages - don't package shame people
########################################

.libPaths(c(
	# "/data/san/data2/users/david/miniforge/envs/R/lib/R/library",
	"/data/san/data0/users/david/rstudio/packages", 
	.libPaths()))
library(data.table)
library(dplyr)
library(gggenomes)
library(thacklr)
library(janitor)
library(svglite)
library(stringr)
library(tidyr)
library(svglite)
sqlFile <- "/data/san/data0/users/david/taxonomizer/accessionTaxa.sql"
texlib <- "/data/san/data0/users/david/rstudio/tinytex"
setwd("/data/san/data1/users/david/ellen/r_gnavus/circular_map")

blue <- rgb(148, 158, 186, maxColorValue = 255)
green <- rgb(116, 199, 186, maxColorValue = 255)


# read in all gbff files use read_gbk
gbff_data <- read_gbk('/data/san/data1/users/david/ellen/r_gnavus/snp_analysis/JCM6515T_spades_scaffold_circular.current.gb')

gbff_df <- gbff_data %>%
	filter(start != 1 & end != 3552639) %>%
	filter(type %in% c("CDS","tRNA")) %>%
	mutate(file_id = "gbff_df") %>% 
	select(where(~ !all(is.na(.)))) %>%
	mutate(color = ifelse(grepl("FXV78_00545", feat_id), "#3C5488FF", "grey")) %>%
	# remove "cds-" from feat_id
	mutate(feat_id = gsub("cds-", "", feat_id))
 
# gbff_df %>% filter(grepl("synthetase",product)) %>% View

focus_locus <- c("FXV78_00545")

p1 <- gggenomes(genes = gbff_df) %>%
	focus(feat_id %in% focus_locus, .expand = c(10000), .overhang = c("trim")) +
	geom_seq() +
	geom_gene(aes(fill = color), color = "black", size=5) +
	scale_fill_identity() +
	geom_bin_label(size = 4, expand_left = 0.2)	+
	geom_gene_label(aes(label = feat_id), size = 2.5) +
	theme(plot.background = element_rect(fill = "transparent", color = NA),
		  panel.background = element_rect(fill = "transparent", color = NA))
	


ggsave("plot_1_ppGpp.png", 
	plot = p1, width = 22, height = 10, dpi = 300, units = "cm")

# save as SVG
detach("package:svglite", unload = TRUE, character.only = TRUE)
.libPaths(c(
	"/data/san/data2/users/david/miniforge/envs/R/lib/R/library"))
library(svglite)
ggsave("plot_1_ppGpp.svg",
	plot = p1, width = 20, height = 4, dpi = 300, units = "cm")

p1_no_labels <- gggenomes(genes = gbff_df) %>%
	focus(feat_id %in% focus_locus, .expand = c(10000), .overhang = c("trim")) +
	geom_seq() +
	geom_gene(aes(fill = color), color = "black", size=5) +
	scale_fill_identity() +
	geom_bin_label(size = 3, expand_left = 0.2) +
	theme(plot.background = element_rect(fill = "transparent", color = NA),
		  panel.background = element_rect(fill = "transparent", color = NA))

ggsave("plot_1_ppGpp_no_labels.png",
	plot = p1_no_labels, width = 20, height = 4, dpi = 300, units = "cm", bg = "transparent")

ggsave("plot_1_ppGpp_no_labels.svg",
	plot = p1_no_labels, width = 20, height = 4, dpi = 300, units = "cm")


# read in katia


# read in all gbff files use read_gbk
gbff_data_ER1_katia <- read_gbk('/data/san/data2/users/david/2024_snp_investigation/data_katia/bakta/ER3-1_1/scaffolds.gbff')

gbff_ER1_katia <- gbff_data_ER1_katia %>%
	filter(type %in% c("CDS","tRNA")) %>%
	mutate(file_id = "gbff_df") %>% 
	select(where(~ !all(is.na(.)))) %>%
	mutate(color = ifelse(grepl("OKFGEK_00005", feat_id), "#3C5488FF", "grey")) %>%
	# remove "cds-" from feat_id
	mutate(feat_id = gsub("cds-", "", feat_id))
 
# gbff_df %>% filter(grepl("synthetase",product)) %>% View

focus_locus_katia <- c("OKFGEK_00005")

p1_ER1_katia <- gggenomes(genes = gbff_ER1_katia) %>%
	flip("NODE_1_length_343981_cov_99.100464") %>%
	focus(feat_id %in% focus_locus_katia, .expand = c(10000), .overhang = c("trim")) +
	geom_seq() +
	geom_gene(aes(fill = color), color = "black", size=5) +
	scale_fill_identity() +
	geom_bin_label(size = 4, expand_left = 0.2)	+
	geom_gene_label(aes(label = feat_id), size = 2.5) +
	theme(plot.background = element_rect(fill = "transparent", color = NA),
		  panel.background = element_rect(fill = "transparent", color = NA))

ggsave("plot_1_ppGpp_katia.png",
	plot = p1_ER1_katia, width = 22, height = 10, dpi = 300, units = "cm")



ava <- read_links("/data/san/data1/users/david/ellen/r_gnavus/circular_map/ava/ava.o6")

# combine dataframes and plot together
common_cols <- intersect(names(gbff_df), names(gbff_ER1_katia))

joined_df <- bind_rows(
	dplyr::select(gbff_df, all_of(common_cols)),
	dplyr::select(gbff_ER1_katia, all_of(common_cols))
) %>%
	filter(seq_id %in% c("NODE_1_length_343981_cov_99.100464","CP043051"))

joined_plot <- gggenomes(genes = joined_df) %>%
	add_links(links = ava) +
	# focus(feat_id %in% c(focus_locus, focus_locus_katia), .expand = c(10000), .overhang = c("trim")) +
	geom_seq() +
	geom_gene(aes(fill = color), color = "black", size = 5) +
	geom_link() +
	scale_fill_identity() +
	geom_bin_label(size = 4, expand_left = 0.2) +
	geom_gene_label(aes(label = feat_id), size = 2.5) +
	theme(
		plot.background = element_rect(fill = "transparent", color = NA),
		panel.background = element_rect(fill = "transparent", color = NA)
	)

ggsave("joined_plot.png",
	plot = joined_plot, width = 22, height = 10, dpi = 300, units = "cm")

joined_plot_no_links <- gggenomes(genes = joined_df) %>%
	flip(2) %>%
	focus(feat_id %in% c(focus_locus, focus_locus_katia), .expand = c(3000), .overhang = c("trim")) +
	geom_seq() +
	geom_gene(aes(fill = color), color = "black", size = 5) +
	scale_fill_identity() +
	geom_bin_label(size = 4, expand_left = 0.2) +
	geom_gene_label(aes(label = feat_id), size = 2.5) +
	theme(
		plot.background = element_rect(fill = "transparent", color = NA),
		panel.background = element_rect(fill = "transparent", color = NA)
	)
	
ggsave("joined_plot_no_links.png",
	plot = joined_plot_no_links, width = 22, height = 10, dpi = 300, units = "cm", bg = "transparent")

# save as SVG
# save as SVG
detach("package:svglite", unload = TRUE, character.only = TRUE)
.libPaths(c(
	"/data/san/data2/users/david/miniforge/envs/R/lib/R/library"))
library(svglite)
ggsave("joined_plot_no_links.svg",
	plot = joined_plot_no_links, width = 22, height = 10, dpi = 300, units = "cm", bg = "transparent")
