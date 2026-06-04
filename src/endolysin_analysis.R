########################################
## Packages - don't package shame people
########################################

.libPaths(c(
	"/data/san/data2/users/david/miniforge/envs/R/lib/R/library",
	"/data/san/data0/users/david/rstudio/packages", 
	.libPaths()))
library(data.table)
library(ape)
library(dplyr)
library(gggenomes)
library(thacklr)
library(tools)
library(msa)
library(Biostrings)
library(tinytex)
library(remotes)
library(data.table)
library(dplyr)
library(stringr)
library(ggtree)
library(ggsci)
library(taxonomizr)
library(Biostrings)
library(stringr)
library(svglite)

sqlFile <- "/data/san/data0/users/david/taxonomizer/accessionTaxa.sql"
texlib <- "/data/san/data0/users/david/rstudio/tinytex"
setwd("/data/san/data1/users/david/ellen/r_gnavus")


faa <- readAAStringSet("endolysin/Rg_endolysins_Renamed.fasta", format = "fasta")


# msa on all sequences
alignment_all = msaMuscle(faa, order = "aligned")


# RAXML was run on the alignment file externally via snakemake
msaPrettyPrint(alignment_all, file="endolysin/Rg_endolysins_Renamed.tex", output="tex",
			   showNames="left", showNumbering="none", showLogo="top",
			   showConsensus="bottom",
			   consensusThreshold= c(50,100),
			   shadingMode = "identical",
			   logoColors="rasmol",
			   verbose=FALSE, askForOverwrite=FALSE, paperWidth = 11, 
			   	paperHeight = 10, 
			   alFile =  "endolysin/Rg_endolysins_Renamed_aln.fasta")

texfile <- "endolysin/Rg_endolysins_Renamed.tex"
tinytex::pdflatex(texfile)



alignment_all = msa(faa, order = "aligned")


# RAXML was run on the alignment file externally via snakemake
msaPrettyPrint(alignment_all, file="endolysin/Rg_endolysins_Renamed.tex", output="tex",
			   showNames="left", showNumbering="none", showLogo="top",
			   showConsensus="bottom",
			   consensusThreshold= c(50,100),
			   shadingMode = "identical",
			   logoColors="rasmol",
			   verbose=FALSE, askForOverwrite=FALSE, paperWidth = 11, 
			   	paperHeight = 10, 
			   alFile =  "endolysin/Rg_endolysins_Renamed_aln.fasta")

texfile <- "endolysin/Rg_endolysins_Renamed.tex"
tinytex::pdflatex(texfile)


########################################
## load tree
########################################
tree <- read.tree("endolysin/Rg_endolysins_Renamed_aln.fasta.raxml.supportFBP")
tree$tip.label
rooting <- "WP_269335750.1"

tree_lysin <- root(tree, outgroup = rooting, resolve.root = TRUE)


ggtree <- ggtree(tree_lysin, size = 0.5, layout = "rectangular")  

gg_tree_plot <- ggtree +
	geom_tiplab(
		size = 3, align = FALSE, linesize = 0.25, offset = 0.3,
		aes(color = label == "LysIBDN1")) +
	scale_color_manual(values = c("TRUE" = "#34528B", "FALSE" = "black"), guide = "none") +
	expand_limits(x = 4, y = 0.1) +
	ggtitle("") +
	theme(
		legend.title = element_text(size = 8, family = "Arial"),
		legend.text  = element_text(size = 6, family = "Arial"),
		legend.key.size = grid::unit(0.3, "cm"),
		text = element_text(family = "Arial"),
		plot.background = element_rect(fill = "transparent", color = NA),
		panel.background = element_rect(fill = "transparent", color = NA)
	) 

ggsave(gg_tree_plot, 
	file ="endolysin/endolysin_tree.png",
	height = 10, width = 6, dpi = 300, unit="cm")

# save as svg
ggsave(gg_tree_plot, 
	file ="endolysin/endolysin_tree.svg",
	height = 10, width = 6, dpi = 300, unit="cm")
