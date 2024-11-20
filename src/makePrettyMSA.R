# Create an MSA of a limited region and add the mammalian ZFX/ZFY phylogeny along side
# Just a pretty colourful image for outreach

#### Imports #####

source("src/functions.R")
load.packages()
ALIGNMENTS <- read.alignments()



# Extract all nt and aa sequences from the MSA for the given region, and calculate the
# consensus sequence
extract.alignment.region <- function(nt.start=NULL, nt.end=NULL){

  nt.aln <- as.data.frame(as.matrix(ALIGNMENTS$nt.mammal.biostrings)[,nt.start:nt.end]) %>%
    dplyr::mutate(Sequence = rownames(.),
                  Sequence = factor(Sequence, levels = mammal.taxa.name.order)) %>%
    dplyr::arrange(as.integer(Sequence)) %>%
    dplyr::rename_with(.cols=starts_with("V"), .fn=\(x) paste0("Site_",as.integer(gsub("V", "", x))+nt.start-1) )

  # find the most common nucleotide in the consensus matrix
  nt.consensus <- paste(unlist(
    apply(consensusMatrix(ALIGNMENTS$nt.mammal.biostrings)[,(nt.start):(nt.end)], 
          2, \(x) names( which(x==max(x)) )[1]  ) # Take only the first consensus in case of ties
  ), 
  collapse = "")
  
  list(
    nt.aln   = nt.aln,
    nt.start = nt.start,
    nt.end   = nt.end,
    nt.length = nt.end - nt.start + 1,
    nt.consensus = nt.consensus
  )
}

# Take a region from the ZF domain that is well conserved
aligned.region <- extract.alignment.region(2124, 2220)$nt.aln %>%
  dplyr::mutate(Sequence = rownames(.)) %>%
  tidyr::pivot_longer(-Sequence, values_to = "Base", names_to = "Site") %>%
  dplyr::mutate(Site = as.numeric(str_extract(Site, "\\d+")))



# Plot an MSA
plot.msa <- ggplot(data=aligned.region,  aes(x = Site, y = Sequence, fill=Base))+
  geom_tile()+
  geom_text(aes(label=Base), col="black", size=1)+
  scale_fill_manual(values = c("A"="lightgreen", "C"="#FFFF66", "G"="salmon", "T"="lightblue"))+
  theme_void()+
  theme(legend.position = "none")

# Read the mammal-only phylogentic tree
nt.aln.tree <- ape::read.tree(FILES$mammal.nt.aln.treefile)

# Root the tree on platypus
nt.aln.tree <- phytools::reroot(nt.aln.tree, which(nt.aln.tree$tip.label=="Platypus_ZFX"), position = 0.015)

# Find the nodes that are ZFY vs ZFX and add to tree
mammal.gene.groups <- split(METADATA$mammal$common.name, METADATA$mammal$group)
nt.aln.tree <- tidytree::groupOTU(nt.aln.tree, mammal.gene.groups, group_name = "group")
nt.aln.tree$tip.label <- rep("", length(nt.aln.tree$tip.label)) # remove tip labels


plot.zfx.zfy <-  ggtree(nt.aln.tree) + 
  geom_tree(col="black", linewidth=0.5) +
  geom_tiplab(align=TRUE, linetype="dotted", linesize=.3) + # use tiplab to get lines
  scale_color_manual(values = c(OUT.TREE.COLOUR, ZFX.TREE.COLOUR, ZFY.TREE.COLOUR))+
  theme_tree() +
  theme(legend.position = "none")+
  coord_cartesian(clip="off", xlim = c(0, 0.25))

# Combine the plots
plot.zfx.zfy + plot.msa + patchwork::plot_layout(widths = c(0.2, 0.8))

ggsave(filename = "figure/ZFX_MSA_Tree.png", last_plot(), dpi=600, units="mm", width=105, height=90)

