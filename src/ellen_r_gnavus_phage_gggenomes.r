#r_gnavus phage gggenomes


.libPaths()
.libPaths(c("/data/san/data0/users/david/rstudio/packages", .libPaths()))
newlib <- "/data/san/data0/users/david/rstudio/packages"
.libPaths()
library(data.table)
library(dplyr)
library(stringr)
library(gggenomes)
library(thacklr)
library(janitor)

awk '/\t/' *.gff > master_tsv_file.tsv  


.libPaths(c(
	"/data/san/data2/users/david/miniforge/envs/R/lib/R/library",
	"/data/san/data0/users/david/rstudio/packages", 
	.libPaths()))
library(svglite)


# read_bl
setwd("/data/san/data1/users/david/ellen/r_gnavus/phage_genomes_2024")
# Load the data
table <- fread("master_tsv_file.tsv", col.names = c("seq_id", "analysis", "type", "start", "end", "score", "strand", "unknown", "dbref"))

extract_value <- function(dbref_entry, field_name) {
  pattern <- paste0(field_name, "=([^;]+)")
  matches <- regmatches(dbref_entry, regexpr(pattern, dbref_entry))
  if (length(matches) > 0) {
    sub(paste0(field_name, "="), "", matches)
  } else {
    NA
  }
}

table$locus_tag <- sapply(table$dbref, extract_value, "locus_tag")
table$fun <- sapply(table$dbref, extract_value, "function")
table$product <- sapply(table$dbref, extract_value, "product")
table = table %>% mutate(feat_id = locus_tag)

table$strand <- case_when(
  table$strand == "+" ~ "+",
  table$strand == "-" ~ "-",
  table$strand == "?" ~ NA,
  TRUE ~ table$strand)

table = table %>% 
	mutate(color = case_when(product == "endolysin" ~ "#34528B",
		TRUE ~ as.character("grey"))) %>%
	mutate(seq_id = case_when(
		seq_id == "MT980836.1" ~ "phiRg507T2/2",
		seq_id == "MT980837.1" ~ "phiRg507T2/3",
		seq_id == "MT980838.1" ~ "phiRg519T2",
		seq_id == "MT980839.1" ~ "phiRgPS6",
		seq_id == "MT980840.1" ~ "phiRgIBDN1",
		TRUE ~ as.character(seq_id))) %>%
  filter(!is.na(strand)) %>%
  filter(strand != "") %>%
  mutate()

ids = unique(table$seq_id)
table[75,strand]

#read in blast
ava = read_sublinks("/data/san/data1/users/david/ellen/r_gnavus/genome_analyis/all_proteins.blast", format = "blast") %>% 
  filter(pident > 35 & evalue < 0.001) %>%
  group_by(feat_id, feat_id2) %>%
  arrange(desc(length), evalue, desc(pident)) %>%
  filter(row_number() == 1) %>%
  ungroup()

plot = gggenomes(genes = table) %>%
	pick_seqs(ids[c(2,5,1,4,3)]) %>%
	add_sublinks(ava) +
    geom_seq() +
	geom_link(aes(fill=pident),
          offset = 0.05) +
	scale_color_gradient(low="grey88", high="grey50") +
  scale_fill_gradient(low="grey88", high="grey50") +
    ggnewscale::new_scale("fill") +
    geom_gene(aes(fill = color),
        size = 10) +
     scale_fill_identity() +
     coord_cartesian(clip = 'off') +
     scale_x_continuous(expand = c(0.2,0)) +
	 geom_bin_label(size = 5)

ggsave(plot, filename = "/data/san/data1/users/david/ellen/r_gnavus/figures/r_gnavus_phages_gggenomes_ellen.png",
    width = 16, height = 6)

# save as svg
ggsave(plot, filename = "/data/san/data1/users/david/ellen/r_gnavus/figures/r_gnavus_phages_gggenomes_ellen.svg",
    width = 16, height = 6)

