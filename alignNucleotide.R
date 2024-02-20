# Create nucleotide alignments for ZFX ad ZFY combined. 
# Make subtrees for just ZFX or ZFY from the combined tree

#### Setup #####
library(tidyverse)
library(ape) # for some plotting
library(msa) # for pretty plotting
library(filesstrings)
library(ggmsa) # alternate pretty plots
library(ggtree)

# Clear previous analyses
filesstrings::dir.remove("aln")
filesstrings::dir.remove("figure")

# Create missing dirs if needed
if(!dir.exists("aln")) dir.create("aln")
if(!dir.exists("bin")) dir.create("bin")
if(!dir.exists("figure")) dir.create("figure")

# Putative Zfy sequences in rodents detected with NCBI gene search:
# rodent[orgn:__txid9989] AND zinc finger X-chromosomal protein-like 

# Read all unaligned sequence files with .fa extension
nt.raw <- do.call(c, lapply(list.files(path = "fasta/nt", pattern = "*.fa$", 
                                       include.dirs = T, full.names = T), 
                            ape::read.FASTA))

# Save the original FASTA header, extract the common name
original.names <- names(nt.raw)
common.names <- str_extract(original.names, "\\[.*\\]") %>%
  str_replace_all("\\[", "")  %>%
  str_replace_all("\\]", "")
names(nt.raw) <- common.names

# Write the combined fasta to file with .fas extension
nt.combined.file <- "fasta/zfxy.nt.fas"
ape::write.FASTA(nt.raw, file = nt.combined.file)

#### Align #####

# Run a codon aware alignment with MACSE
# Expect java on the PATH. Macse download from https://www.agap-ge2pop.org/macsee-pipelines/

# Direct download link:
# https://www.agap-ge2pop.org/wp-content/uploads/macse/releases/macse_v2.07.jar
if(!file.exists("bin/macse_v2.07.jar")) stop("MACSE not present in bin directory")

# Run the alignment
nt.aln.file <- "aln/zfxy.nt.aln"
aa.aln.file <- "aln/zfxy.aa.aln"
system2("java", paste("-jar bin/macse_v2.07.jar -prog alignSequences",
                      "-seq", nt.combined.file, # input
                      "-out_NT", nt.aln.file,  # output nt alignment
                      "-out_AA", aa.aln.file), # output aa alignment
        stdout = paste0(nt.aln.file, ".macse.log"),  # logs
        stderr = paste0(nt.aln.file, ".macse.log"))  # error logs
macse.aln <- ape::read.FASTA(nt.aln.file)
# ape::image.DNAbin(macse.aln, "-", base.col = "white")


#### Tree #####

# Make ML tree and reconstruct ancestral sequences
# Expect iqtree on the PATH.
# Note model testing is automatically performed in v1.5.4 onwards
# Note: we can use a partition model if we specify exon coordinates
system2("iqtree", paste("-s ", nt.aln.file, 
                        "-bb 1000", # number of bootstrap replicates
                        "-alrt 1000", # number of replicates to perform SH-like approximate likelihood ratio test (SH-aLRT) 
                        "-nt AUTO", # number of threads
                        "-asr"), # ancestral sequence reconstruction
        stdout = paste0(nt.aln.file, ".iqtree.log"), 
        stderr = paste0(nt.aln.file, ".iqtree.log"))

#### Plot tree #####

macse.tree <- ape::read.tree(paste0(nt.aln.file, ".treefile"))
# Root the tree on platypus
macse.tree <- ape::root(macse.tree, "Platypus_ZFX")

# Find the nodes that are ZFY vs ZFX and add to tree
zfy.nodes <- gsub(" ", "_", common.names[str_detect(common.names, "Z[F|f][Y|y|a]")])
zfy.nodes <- gsub("\\?", "_", zfy.nodes)
macse.tree <- groupOTU(macse.tree, zfy.nodes, group_name = "ZFY")

plot.zfx.zfy <- ggtree(macse.tree) + 
  geom_tree() +
  theme_tree() +
  geom_tiplab(aes(col = ZFY)) +
  geom_treescale()+
  coord_cartesian(clip="off")+
  xlim(0, 0.6) +
  theme(legend.position = "none")

ggsave(paste0("figure/zfx.zfy.tree.png"),  
       plot.zfx.zfy, dpi = 300, units = "mm", width = 170, height = 170)


# Drop the ZFY and Zfa sequences and just look at the ZFX nodes in the tree
tree.zfx <- ape::drop.tip(macse.tree, zfy.nodes)
ape::write.tree(tree.zfx, file = paste0(nt.aln.file, ".zfx.treefile"))
plot.zfx <- ggplot(tree.zfx, aes(x, y)) + 
  geom_tree() + 
  theme_tree() +
  geom_tiplab() +
  geom_treescale()+
  coord_cartesian(clip="off")+
  xlim(0, 0.6)
ggsave(paste0("figure/zfx.tree.png"),  
       plot.zfx, dpi = 300, units = "mm", width = 170, height = 170)

# Keep Zfy and Zfa, drop Zfx nodes. Keep outgroups
tree.zfy <- ape::keep.tip(macse.tree, c(zfy.nodes, "Platypus_ZFX", "Opossum_ZFX"))
ape::write.tree(tree.zfy, file = paste0(nt.aln.file, ".zfy.treefile"))

plot.zfy <- ggplot(tree.zfy, aes(x, y)) + 
  geom_tree() + 
  theme_tree() +
  geom_tiplab() +
  geom_treescale()+
  coord_cartesian(clip="off")+
  xlim(0, 0.6)
ggsave(paste0("figure/zfy.tree.png"),  
       plot.zfy, dpi = 300, units = "mm", width = 170, height = 170)


#### Ancestral nodes #####

# Ancestral reconstruction
ancestral.seqs <- read.table(paste0(nt.aln.file, ".state") ,header=TRUE)

#### Exon by exon MSA #####

msa.aln = readDNAMultipleAlignment(nt.aln.file, format="fasta")

# Check for manual exon boundary calling
# tidy.msa <- tidy_msa(msa.aln) %>%
#   dplyr::filter(name == "Mouse Zfy1") 
# paste(tidy.msa$character, collapse = "")

# Exon by exon coordinates of the alignment will be needed for clear
# testing of selection. Match these from the final alignment.
exons <- data.frame("exon" = c("3", "4", "5",  "6", "7", "8",  "9"), 
                    "start" = c( 1,  62, 689,  839, 1016, 1142, 1373),
                    "end"   = c(61, 688, 838, 1015, 1141, 1372, 2580))

for(i in seq(1, nrow(exons))){
  start <- exons$start[i]
  end <- exons$end[i]
  exon <- exons$exon[i]
  msa.plot <- ggmsa(msa.aln, start, end, seq_name = T, font=NULL, 
                    border=NA, color="Chemistry_NT") +
    theme(axis.text = element_text(size=2))+
    facet_msa(field = 100)
  
  ggsave(paste0("figure/msa_exon_", exon, ".png"),  
         msa.plot, dpi = 300, units = "mm", width = 170)
}

#### Plot full MSA #####

msa.plot <- ggmsa(msa.aln, start = 1, end=2600, seq_name = F, 
                  border = NA, font=NULL, color="Chemistry_NT") + 
  facet_msa(field = 400)
ggsave(paste0("figure/msa.complete.png"),  
       msa.plot, dpi = 300, units = "mm", width = 170, height = 170)

#### Create partition model for IQ-TREE #####

