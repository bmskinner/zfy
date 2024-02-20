# Identify rodent species with WGS data for extra ZFX and ZFY sequences beyond
# annotations

library(rentrez)
library(tidyverse)
library(ape)
# List search fields
# entrez_db_searchable("sra")

# search.1 <- entrez_search(db="sra", term = "male[word] AND dna[word] AND WGS[word] AND paired[lay] AND rodent[orgn:__txid9989]")
# data <- entrez_fetch(db="sra", id=search.1$ids, rettype = "xml", parsed = T)
# data.list <- XML::xmlToList(data)

# Align existing species
nt.raw <- do.call(c, lapply(list.files(path = "fasta/nt", include.dirs = T, full.names = T), ape::read.FASTA))

fasta <- ape::read.FASTA("Jerboa_Zfx_Zfy.fa")
existing <- ape::read.FASTA("ZFY shared evolution folder/Nucleotide sequences/ZFY+ZFX_All nucleotide sequences.fa")
fasta <- c(fasta, existing)
# Write the combined fasta to file
ape::write.FASTA(fasta, file = "Jerboa_combined.fas")

# Align with muscle
# aln <- ape::muscle5(fasta, mc.cores = 10, MoreArgs = "")
# ape::image.DNAbin(aln)

# Run iqtree to find the best model and create a tree
# system2("iqtree","-s Jerboa_combined.fas -st CODON -bb 1000 -alrt 1000 -nt AUTO", stdout = TRUE, stderr = TRUE)

# Root the tree on platypus
# tree <- ape::read.tree("Jerboa_combined.fas.treefile")
# tree <- ape::root(tree, "XM_029079877.2_PREDICTED_Ornithorhynchus_1-2445")
# plot(tree, font = 1)


# Try a codon aware alignment
# java -jar c:/portable/macse_v2.07.jar
system2("java","-jar c:/portable/macse_v2.07.jar -prog alignSequences -seq Jerboa_combined.fas", stdout = TRUE, stderr = TRUE)
macse.aln <- ape::read.FASTA("Jerboa_combined_NT.fas")
ape::image.DNAbin(macse.aln)
# reconstruct ancestral sequences
system2("iqtree","-s Jerboa_combined_NT.fas -bb 1000 -alrt 1000 -nt AUTO -asr", stdout = TRUE, stderr = TRUE)
macse.tree <- ape::read.tree("Jerboa_combined_NT.fas.treefile")
macse.tree <- ape::root(macse.tree, "XM_029079877.2_PREDICTED_Ornithorhynchus_1-2445")
plot(macse.tree, show.node.label = T, align.tip.label = T)

# Ancestral reconstruction
tab <- read.table('Jerboa_combined_NT.fas.state',header=TRUE)
#
# ggplot(tab, aes(x=Site))+
#   geom_line(aes(y = p_A), alpha = 0.2)+
#   facet_wrap(~Node)
