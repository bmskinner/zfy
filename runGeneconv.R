# Run GENECONV to test for gene conversion
# GENECONV is expected on the path
source("functions.R")
load.packages()

fa.files <- list.files(path = "fasta/nt", pattern = "*.fa$", 
                       include.dirs = T, full.names = T)

fa.read  <- lapply(fa.files, read.fasta)
metadata.mammal <- read.metadata(fa.read)

# We want to look at gene conversion within species lineages. We also want to
# compare the ancestral ZFXs and ZFYs at each node To do this, create a separate
# tree for each of ZFX and ZFY, confirm that they have equivalent branches, then
# add the ancestral reconstruction to the geneconv config file.
nt.aln.file <- "aln/zfxy.nt.aln"
msa.nt.aln <- Biostrings::readDNAMultipleAlignment(nt.aln.file, format="fasta")
ape.nt.aln <- ape::read.FASTA(nt.aln.file)

# Write ZFX sequences to file
zfx.nt.aln <- msa.nt.aln@unmasked[names(msa.nt.aln@unmasked) %in% c(metadata.mammal$common.name[metadata.mammal$group=="ZFX"], 
                                                                    metadata.mammal$common.name[metadata.mammal$group=="Outgroup"])]
Biostrings::writeXStringSet(zfx.nt.aln,  file = "aln/zfx_only/zfx.aln", format = "fasta")

# Write ZFY sequences to file
zfy.nt.aln <- msa.nt.aln@unmasked[names(msa.nt.aln@unmasked) %in% c(metadata.mammal$common.name[metadata.mammal$group=="ZFY"], 
                                                                    metadata.mammal$common.name[metadata.mammal$group=="Outgroup"])]
Biostrings::writeXStringSet(zfy.nt.aln,  file = "aln/zfy_only/zfy.aln", format = "fasta")

# Align the sequences with IQTREE
system2("iqtree", paste("-s ", "aln/zfx_only/zfx.aln", 
                        "-bb 1000", # number of bootstrap replicates
                        "-alrt 1000", # number of replicates to perform SH-like approximate likelihood ratio test (SH-aLRT) 
                        "-nt AUTO", # number of threads
                        "-asr")) # ancestral sequence reconstruction
system2("iqtree", paste("-s ", "aln/zfy_only/zfy.aln", 
                        "-bb 1000", # number of bootstrap replicates
                        "-alrt 1000", # number of replicates to perform SH-like approximate likelihood ratio test (SH-aLRT) 
                        "-nt AUTO", # number of threads
                        "-asr")) # ancestral sequence reconstruction

# Read the ZFX and ZFY trees
zfx.nt.aln.tree <- ape::read.tree("aln/zfx_only/zfx.aln.treefile")
zfy.nt.aln.tree <- ape::read.tree("aln/zfy_only/zfy.aln.treefile")

# Drop the second ZFYs in mouse and rat
zfy.nt.aln.tree <- tidytree::drop.tip(zfy.nt.aln.tree, "Mouse_Zfy2") 
zfy.nt.aln.tree <- tidytree::drop.tip(zfy.nt.aln.tree, "African_Grass_Rat_ZFY2-like_1") 

# Root the trees on platypus
zfx.nt.aln.tree <- phytools::reroot(zfx.nt.aln.tree, which(zfx.nt.aln.tree$tip.label=="Platypus_ZFX"), position = 0.015)
zfy.nt.aln.tree <- phytools::reroot(zfy.nt.aln.tree, which(zfy.nt.aln.tree$tip.label=="Platypus_ZFX"), position = 0.015)

# Remove gene names so tip labels are comparable
zfx.nt.aln.tree$tip.label <- str_replace(zfx.nt.aln.tree$tip.label, "_Z[F|f][X|x].*", "")
zfy.nt.aln.tree$tip.label <- str_replace(zfy.nt.aln.tree$tip.label, "(_putative)?(_|-)Z[F|f][X|x|Y|y].*", "")

# Export comparison of the trees
png(filename = "figure/zfx.zfy.nt.ancestral.treediff.png")
treespace::plotTreeDiff(zfx.nt.aln.tree, zfy.nt.aln.tree, treesFacing=TRUE)
dev.off()

# Plot the two trees with node labels
zfx.nt.aln.tree.plot <- plot.tree(zfx.nt.aln.tree)
zfy.nt.aln.tree.plot <- plot.tree(zfy.nt.aln.tree)
zfx.zfy.aln.tree.plot <- zfx.nt.aln.tree.plot + zfy.nt.aln.tree.plot + patchwork::plot_annotation(tag_levels = list(c("ZFX", "ZFY")))
save.double.width("figure/zfx.zfy.ancestral.comparison.tree.png", zfx.zfy.aln.tree.plot)

# Find the matching nodes
zfy.zfx.common.nodes <- ape::comparePhylo(zfx.nt.aln.tree, zfy.nt.aln.tree)$NODES %>%
  tidyr::separate_wider_delim(cols = zfx.nt.aln.tree, delim = " ", names = c("ZFX_node", "zfx_nnodes") ) %>%
  tidyr::separate_wider_delim(cols = zfy.nt.aln.tree, delim = " ", names = c("ZFY_node", "zfy_nnodes") ) %>%
  # Remove booststrap values from node names
  dplyr::mutate(ZFX_node = str_replace(ZFX_node, "\\/.*", ""),
                ZFY_node = str_replace(ZFY_node, "\\/.*", "")) %>% 
  dplyr::filter(ZFX_node!="Root" & ZFX_node!="") 

# Read the ancestral states for the nodes
ancestral.zfx.seqs <- read.table("aln/zfx_only/zfx.aln.state", header=TRUE)
ancestral.zfy.seqs <- read.table("aln/zfy_only/zfy.aln.state", header=TRUE)

get.node.sequence <- function(node.name, ancestral.seq.table, type){
  ancestral.seq.table %>%
    dplyr::filter(Node==node.name) %>%
    dplyr::arrange(Site) %>%
    dplyr::select(State) %>%
    dplyr::summarise(Seq =  paste0(">", type, "_", node.name, "\n", paste(State, collapse = "")))
}

# Get the ancestral sequences for these node pairs
zfy.zfx.common.nodes$anc.zfx <- sapply(zfy.zfx.common.nodes$ZFX_node, 
                                       get.node.sequence, 
                                       ancestral.seq.table=ancestral.zfx.seqs, type="ZFX")
zfy.zfx.common.nodes$anc.zfy <- sapply(zfy.zfx.common.nodes$ZFY_node, 
                                       get.node.sequence, 
                                       ancestral.seq.table=ancestral.zfy.seqs, type="ZFY")

# Write the node sequences to fasta
anc.node.seqs <- paste0(c(zfy.zfx.common.nodes$anc.zfx, zfy.zfx.common.nodes$anc.zfy, "\n"), collapse = "\n")
write_file(anc.node.seqs, "aln/ancestral.zfx.zfy.nodes.fa")

# Combine with existing ZFX/Y sequence file
ape::write.dna(ape.nt.aln, "aln/ancestral.zfx.zfy.nodes.fa", format="fasta", 
               append = TRUE, colsep = "", nbcol=-1)

# Write geneconv configuration file specifying which sequences are in the same
# group

node.string.names <- zfy.zfx.common.nodes %>% 
  dplyr::mutate(GroupString= paste0("-group ZFX_", ZFX_node, "_ZFY_",ZFY_node, " ZFX_", ZFX_node, " ZFY_",ZFY_node ))
node.string <- paste(node.string.names$GroupString, collapse = "\n")

split.gene.names <- metadata.mammal %>%
  # Since we're about to split on 'Z', remove secondary instances of the string
  # keeping whatever - or _ was in the original name
  dplyr::mutate(common.name = str_replace(common.name, "putative-Zfy", "putative-Y"),
                common.name = str_replace(common.name, "putative_Zfy", "putative_Y")) %>%
  separate_wider_delim(common.name, 
                       delim= "Z", names = c("Species", "Gene"), too_few = "debug") %>%
  dplyr::mutate(Species = str_replace(Species, "_$", ""),
                Gene = paste0("Z", Gene)) %>%
  # Repair gene names following the split earlier
  dplyr::mutate( Gene = str_replace(Gene, fixed("Zfx-like_putative_Y"), "Zfx-like_putative_Zfy"),
                 Gene = str_replace(Gene, fixed("ZFX-like_putative_Y"), "ZFX-like_putative_ZFY"),
                 Gene = str_replace(Gene, fixed("Zfx-like_putative-Y"), "Zfx-like_putative-Zfy"),
                 Gene = str_replace(Gene, fixed("ZFX-like_putative-Y"), "ZFX-like_putative-ZFY"),
                 Seq = paste0(Species, "_", Gene)) %>%
  dplyr::arrange(Species, Gene) %>%
  dplyr::group_by(Species) %>%
  dplyr::reframe(GroupString =  paste("-group", Species, paste(Seq, collapse = " "))) %>%
  dplyr::distinct()


group.string <- paste(split.gene.names$GroupString, collapse = "\n")

# Make the geneconv control file
# e.g 
#GCONV_CONFIG
# -Startseed=123   -MaxSimGlobalPval=0.05
# -group GRPI   S1  S2
# -group GRPII  H1 C1 I2
# -group GRPII  O1  O2
# -group GRPIII G1 G2 R1 R2 
anc.gene.conv.control <- paste0("#GCONV_CONFIG\n",
                                "-Startseed=123   -MaxSimGlobalPval=0.05\n",
                                group.string, "\n", node.string)

write_file(anc.gene.conv.control, "aln/anc.zfx.zfy.geneconv.cfg")

#  Only on PATH in Windows, not running on cluster
system2("geneconv", paste("aln/ancestral.zfx.zfy.nodes.fa",
                          "aln/anc.zfx.zfy.geneconv.cfg",
                          "aln/anc.zfx.zfy.geneconv.out /lp"),
        stdout = "aln/ancestral.zfx.zfy.nodes.geneconv.log", 
        stderr = "aln/ancestral.zfx.zfy.nodes.geneconv.log")


# Read the geneconv output file
aln.geneconv.data <- read_table("aln/anc.zfx.zfy.geneconv.frags", comment = "#", 
                                col_names = c("Type", "Pair", "Sim_Pvalue", "KA_pvalue", "begin", "end", 
                                              "len", "num_poly", "num_dif", "tot_difs", "mism_pen")) %>%
  na.omit %>%
  dplyr::mutate(Species = str_replace_all(Pair, "_Z[F|f][Y|y|X|x].*", ""),
                Type = case_when(Type == "GO" ~ "Global outer",
                                 Type == "GI" ~ "Global inner",
                                 Type == "PO" ~ "Pairwise outer",
                                 Type == "PI" ~ "Pairwise inner",)) %>%
  dplyr::filter(Type == "Pairwise inner" | Type == "Pairwise outer") %>%
  dplyr::mutate(y = row_number())

# What are the species descending from each node?
# zfy.zfx.common.nodes
# zfx.nt.aln.tree
# offspring(zfx.nt.aln.tree, zfy.zfx.common.nodes$ZFX_node[1])




# types:
# Global - looks across the whole dataset
# Pairwise - looks at the paired list given
# Inner -  evidence of a possible gene conversion event between ancestors of two sequences in the alignment. 
# Outer - evidence of past gene conversion events that may have originated from outside of the alignment,
# or else from within the alignment but such that evidence of the source has been destroyed by later mutation or gene conversion
# GI: Global inner (Inner fragments are runs of matching sites)
# GO: Global outer (Outer-sequence fragments are runs of sites that are unique in that group)
# PI: Pairwise inner
# PO: Pairwise outer-sequence

# We only really care about pairwise inner
# Plot the geneconv fragements in the ancestral nodes

aln.geneconv.plot <- ggplot(aln.geneconv.data)+
  # geom_rect(data=mouse.exons, aes( xmin = start-0.5, xmax = end+0.5, ymin=0, ymax=0.08, fill=exon), alpha=0.5)+
  # geom_text(data=mouse.exons, aes(x=(start+end)/2,y = 0.075, label = exon), size = 3)+
  # geom_rect(data = locations.zf, aes(xmin=start_nt, xmax=end_nt, ymin=i-0.5, ymax=i+0.5), fill="grey", alpha=0.5)+
  # geom_rect(data = locations.9aaTAD, aes(xmin=start_nt, xmax=end_nt, ymin=i-0.5, ymax=i+0.5), fill="blue", alpha=0.5)+
  # geom_rect(data = locations.NLS, aes(xmin=start_nt, xmax=end_nt, ymin=i-0.5, ymax=i+0.5), fill="green", alpha=0.5)+
  geom_segment(aes(x=begin, y = y, col=Species, xend = end, yend = y), linewidth = 2)+
  # coord_cartesian(xlim = c(1, 2600))+
  # scale_fill_manual(values = c("grey", "white", "grey", "white", "grey", "white", "grey"))+
  # scale_y_continuous(labels = rev(taxa.name.order), breaks = 1:length(taxa.name.order))+
  # geom_text(aes(x=(begin+end)/2,y = KA_pvalue+0.002, label = Pair), size = 2)+
  labs(x = "Position", y = "Bonferroni corrected p-value")+
  guides(fill = "none")+
  theme_bw()+
  theme(axis.title = element_blank())

save.double.width("figure/aln.ancestral.zfxy_geneconv.png", aln.geneconv.plot)
