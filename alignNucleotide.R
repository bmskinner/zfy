# Create nucleotide alignments for ZFX ad ZFY combined. 
# Make subtrees for just ZFX or ZFY from the combined tree

#### Setup #####
library(tidyverse)
library(ape) # for some plotting
library(msa) # for pretty plotting
library(filesstrings)
library(ggmsa) # alternate pretty plots
library(ggtree)
library(seqinr)
library(phangorn)
library(installr)
library(treeio)
library(httr)
library(seqLogo)
library(seqvisr) # remotes::install_github("vragh/seqvisr")

#### Read files #####

# Clear previous analyses
filesstrings::dir.remove("aln")
filesstrings::dir.remove("figure")

# Create missing dirs if needed
if(!dir.exists("aln")) dir.create("aln")
if(!dir.exists("aln/outgroup")) dir.create("aln/outgroup")
if(!dir.exists("aln/exons")) dir.create("aln/exons")
if(!dir.exists("bin")) dir.create("bin")
if(!dir.exists("figure")) dir.create("figure")
if(!dir.exists("paml/site-specific")) dir.create("paml/site-specific")
if(!dir.exists("paml/site-branch")) dir.create("paml/site-branch")

# Putative Zfy sequences in rodents detected with NCBI gene search:
# rodent[orgn:__txid9989] AND zinc finger X-chromosomal protein-like 

# Read FASTA file and extract metadata
read.fasta <- function(f){
  fa.data <- ape::read.FASTA(f)
  original.name <- names(fa.data)
  common.names <-  str_extract(original.name, "\\[.*\\]") %>%
    str_replace_all("\\[", "")  %>%
    str_replace_all("\\]", "") %>%
    str_replace_all(" ", "_")
  names(fa.data) <- common.names
  
  species.name <-gsub("_", " ",  gsub(".fa$", "", gsub("fasta/(nt|aa)/", "",  f)))
  
  list("fa" = fa.data,
       metadata = data.frame(
        "accession" = str_extract(original.name, "[^:]+"),
       "original.name" = original.name,
       "common.name" = common.names,
       "species" = rep(species.name, length(common.names))))
}

# Read all unaligned sequence files with .fa extension
fa.files <- list.files(path = "fasta/nt", pattern = "*.fa$", 
                       include.dirs = T, full.names = T)

fa.read <- lapply(fa.files, read.fasta)

nt.raw <- do.call(c, lapply(fa.read, function(x) x$fa))
metadata <- do.call(rbind, lapply(fa.read, function(x) x$metadata))

# Grouping for trees
metadata %<>%
  dplyr::mutate(group = case_when(grepl("ZFY", common.name, ignore.case=T) ~ "ZFY",
                                  common.name %in% c("Platypus_ZFX", "Opossum_ZFX") ~ "Outgroup",
                                  T ~ "ZFX"))

# Write the combined fasta to file with .fas extension
nt.combined.file <- "fasta/zfxy.nt.fas"
ape::write.FASTA(nt.raw, file = nt.combined.file)

#### Create NT and AA initial alignment #####

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

# Convert to interleaved for display
# write.dna(macse.aln, file="aln/nt.interleaved.fa", format = "interleaved")

# nt.phylip.phy.file <- "aln/nt.phylip.phy"
# write.phyDat(macse.aln, file=nt.phylip.phy.file, format = "phylip")

# Exon by exon coordinates of the alignment will be needed for clear
# testing of selection. Match these from the final alignment via mouse Zfy1
find.exons <- function(){
  macse.nt.aln <- seqinr::read.alignment(nt.aln.file, format="fasta")
  mouse.zfy1 <- toupper(macse.nt.aln$seq[macse.nt.aln$nam=="Mouse_Zfy1"])
  
  mouse.exons <- data.frame("exon" = c("1", "2", "3",  "4", "5", "6",  "7"),
                            "start" = c("ATGGATGAA", "GAGCTGATGCA", "TGGATGAACC", "GAGAAACTAT", "AAGTAATTGT", "ATAATAATTCT", "CAATATTTGTT"),
                            "end" = c("TGGAATAG", "ATGATGTCTT", "GGATGAATTAG", "---------ACTG", "GACAGCAGCTTATG", "CAGTACCAGTCAG", "CCTGCCCTAA"))
  
  find.exon.bounds <- function(exon, start.seq, end.seq){
    start <- str_locate(mouse.zfy1, start.seq)
    end <- str_locate(mouse.zfy1, end.seq)
    data.frame("exon" = exon, "start" = start[1], "end"= end[2])
  }
  
  do.call(rbind, lapply(1:nrow(mouse.exons), function(i) find.exon.bounds(exon = mouse.exons$exon[i], start.seq = mouse.exons$start[i], end.seq = mouse.exons$end[i])))
}

mouse.exons <- find.exons()

# Find the acidic domain and ZFs in the alignment from mouse Zfy1
find.domains <- function(){
  
}

mouse.domains <- find.domains()


# Convert FASTA format to clustal style format
# alignment - seqinr alignment format
# names - optional character vector of names (if null, alignment names are used)
# names.length - override the default width of the names column (if na, default is used)
# chunksize - number of letters per row
printMultipleAlignment <- function(alignment, names=NULL, names.length=NA, chunksize=60){
  # this function requires the Biostrings package
  # find the number of sequences in the alignment
  numseqs <- alignment$nb
  
  if(is.null(names)){
    names <- alignment$nam
  }
  
  if(is.na(names.length)){
    names.length <- max(nchar(names))
  }
  
  # find the length of the alignment
  alignmentlen <- nchar(alignment$seq[[1]])
  
  # Calculate the start position of each line
  line.starts <- seq(1, alignmentlen, by=chunksize)
  
  # How many blocks are needed
  n.blocks <- length(line.starts)
  
  # get the alignment for each  sequences
  aln <- unlist(alignment$seq)
  lettersprinted <- rep(0, numseqs)
  
  create.block <- function(start){
    block.lines <- rep("", numseqs+1)
    block.lines[numseqs+1] <- "\n"
    for (j in 1:numseqs){
      alnj <- aln[j]
      chunkseq <- toupper(substring(alnj, start, start+chunksize-1))
      
      # Calculate how many residues of the sequence we have printed so far in the alignment
      # Total minus gaps
      lettersprinted[j] <<- lettersprinted[j] + chunksize - Biostrings::countPattern("-",chunkseq)
      block.lines[j] <- paste0(sprintf( paste0("%", names.length,"s"), names[j]), "\t", chunkseq, " ", lettersprinted[j])
    }

    paste0(block.lines, "\n")
  }
  
  result <- unlist(lapply(line.starts, create.block))
  cat(paste0(result, "\n"))
  result
}

# Check the AA alignment
macse.aa.aln <- seqinr::read.alignment(aa.aln.file, format="fasta")
# aa.clustal <- paste(printMultipleAlignment(macse.aa.aln, names = common.names), collapse = "")
# write_file(aa.clustal, file="aln/aa.clustal.aln")

# Calculate conservation at each site
msa.nt.aln <- Biostrings::readDNAMultipleAlignment(nt.aln.file, format="fasta")
tidy.msa <- ggmsa::tidy_msa(msa.nt.aln) %>%
  dplyr::filter(character != "-") %>%
  dplyr::group_by(position, character) %>%
  dplyr::summarise(n = n(), fraction = n/nrow(msa.nt.aln)) %>%
  dplyr::group_by(position) %>%
  dplyr::arrange(desc(fraction)) %>%
  dplyr::slice_head(n=1) %>%
  dplyr::arrange(position)
        
# Plot exon by exon conservation           
conservation.plot <- ggplot(tidy.msa)+
  geom_rect(data=mouse.exons, aes( xmin = start-0.5, xmax = end+0.5, ymin=0, ymax=1, fill=exon), alpha=1)+
  geom_line( aes(x=position, y=fraction))+
  scale_fill_manual(values = rep(c("grey", "white"), 4))+
  labs(x = "Position in alignment", y = "Fraction of species with consensus nucleotide")+
  theme_bw()


# macse.aa.plot <- ggmsa(macse.aa.aln, start = 1, end=860, seq_name = T, 
#       border = NA, font=NULL, color="Chemistry_AA") + 
#   facet_msa(field = 100)+
#   theme(axis.text.y = element_text(size = 2))
# ggsave(paste0("figure/msa.aa.complete.png"),  
#        macse.aa.plot, dpi = 300, units = "mm", width = 170, height = 370)

#### Make tree from whole CDS #####

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

#### Plot whole CDS tree #####

save.double.width <- function(filename, plot) ggsave(filename, plot, dpi = 300, 
                                                     units = "mm", width = 170, 
                                                     height = 170)

plot.tree <- function(tree.data){
  p <- ggtree(tree.data) + 
    geom_tree() +
    geom_tiplab(aes(col = group), size=2)+
    geom_nodelab(size=2, nudge_x = -0.003, nudge_y = 0.5, hjust=1,  node = "internal")+
    geom_treescale(fontsize =2, y = -1)+
    coord_cartesian(clip="off")+
    theme_tree() +
    theme(legend.position = "none")
  p
}

macse.tree <- ape::read.tree(paste0(nt.aln.file, ".treefile"))

# Root the tree on platypus and opossum and resave
macse.tree <- ape::root(macse.tree, outgroup = c("Platypus_ZFX"), resolve.root = TRUE)
ape::write.tree(macse.tree, file = paste0(nt.aln.file, ".rooted.treefile"))

# Remove node names, leaving just bootstrap values
macse.tree$node.label <- gsub("(Node\\d+\\/|Root)", "", macse.tree$node.label)

# Find the nodes that are ZFY vs ZFX and add to tree
group_info <- split(metadata$common.name, metadata$group)
macse.tree <- groupOTU(macse.tree, group_info, group_name = "group")

plot.zfx.zfy <- plot.tree(macse.tree) + coord_cartesian(clip="off", xlim = c(0, 0.5))
save.double.width("figure/zfx.zfy.tree.png", plot.zfx.zfy)


# Drop the ZFY sequences and just look at the ZFX nodes in the tree
tree.zfx <- ape::drop.tip(macse.tree, group_info$ZFY)
tree.zfx <- groupOTU(tree.zfx, group_info, group_name = "group")
ape::write.tree(tree.zfx, file = paste0(nt.aln.file, ".zfx.treefile"))
plot.zfx <- plot.tree(tree.zfx)
save.double.width("figure/zfx.tree.png", plot.zfx)

# Keep Zfy and Zfa, drop Zfx nodes. Keep outgroups
tree.zfy <- ape::keep.tip(macse.tree, c(group_info$ZFY, group_info$Outgroup))
tree.zfy <- groupOTU(tree.zfy, group_info, group_name = "group")
ape::write.tree(tree.zfy, file = paste0(nt.aln.file, ".zfy.treefile"))
plot.zfy <- plot.tree(tree.zfy)
save.double.width("figure/zfy.tree.png", plot.zfy)

#### Make tree from each exon separately ####

# Aim here is to look for signs of gene conversion in the final exon as per other 
# papers.

for(i in 1:nrow(mouse.exons)){

  exon.aln <- as.matrix(macse.aln)[,mouse.exons$start[i]:mouse.exons$end[i]]
  exon.aln.file <- paste0("aln/exons/exon_", mouse.exons$exon[i], ".aln")
  ape::write.FASTA(exon.aln, file = exon.aln.file)

  
  system2("iqtree", paste("-s ", exon.aln.file,
                          "-bb 1000", # number of bootstrap replicates
                          "-alrt 1000", # number of replicates to perform SH-like approximate likelihood ratio test (SH-aLRT)
                          "-nt AUTO"), # number of threads
          stdout = paste0(exon.aln.file, ".iqtree.log"),
          stderr = paste0(exon.aln.file, ".iqtree.log"))
  
  # Some exons will fail - too many gaps
  if(!file.exists(paste0(exon.aln.file, ".treefile"))) next
  
  exon.tree <- ape::read.tree(paste0(exon.aln.file, ".treefile"))
  # Root the tree on platypus and opossum
  exon.tree <- ape::root(exon.tree, c("Platypus_ZFX"), resolve.root = TRUE)
  
  # Find the nodes that are ZFY vs ZFX and add to tree
  group_info <- split(metadata$common.name, metadata$group)
  exon.tree <- groupOTU(exon.tree, group_info, group_name = "group")
  
  plot.exon.tree <- plot.tree(exon.tree)
  exon.fig.file <- paste0("figure/exon_", mouse.exons$exon[i], ".zfx.zfy.tree.png")
  save.double.width(exon.fig.file, plot.exon.tree)
  
}

# We also want to look at all except exon 9

exon.aln <- as.matrix(macse.aln)[,mouse.exons$start[1]:mouse.exons$end[6]]
exon.aln.file <- paste0("aln/exons/exon_1-7.aln")
ape::write.FASTA(exon.aln, file = exon.aln.file)

system2("iqtree", paste("-s ", exon.aln.file,
                        "-bb 1000", # number of bootstrap replicates
                        "-alrt 1000", # number of replicates to perform SH-like approximate likelihood ratio test (SH-aLRT)
                        "-nt AUTO"), # number of threads
        stdout = paste0(exon.aln.file, ".iqtree.log"),
        stderr = paste0(exon.aln.file, ".iqtree.log"))

exon.tree <- ape::read.tree(paste0(exon.aln.file, ".treefile"))
# Root the tree on platypus and opossum
exon.tree <- ape::root(exon.tree, c("Platypus_ZFX"), resolve.root = TRUE)

# Find the nodes that are ZFY vs ZFX and add to tree
# Find the nodes that are ZFY vs ZFX and add to tree
group_info <- split(metadata$common.name, metadata$group)
exon.tree <- groupOTU(exon.tree, group_info, group_name = "group")

plot.exon.tree <- plot.tree(exon.tree)
exon.fig.file <- paste0("figure/exon_1-7.zfx.zfy.tree.png")
save.double.width(exon.fig.file, plot.exon.tree)

#### Run GENECONV to test for gene conversion ####

# We want to look at gene conversion within species lineages only
# Write a geneconv configuration file specifying which sequences are in the same
# group

# e.g 
#GCONV_CONFIG
# -Startseed=123   -MaxSimGlobalPval=0.05
# -group GRPI   S1  S2
# -group GRPII  H1 C1 I2
# -group GRPII  O1  O2
# -group GRPIII G1 G2 R1 R2 

split.names <- metadata %>%
  # Since we're about to split on 'Z', remove secondary instances of the string
  # keeping whatever - or _ was in the original name
  dplyr::mutate(common.name = str_replace(common.name, "putative-Zfy", "putative-Y"),
                common.name = str_replace(common.name, "putative_Zfy", "putative_Y")) %>%
  separate_wider_delim(common.name, 
                       delim= "Z", names = c("Species", "Gene")) %>%
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

group.string <- paste(split.names$GroupString, collapse = "\n")

# Make the geneconv control file
gene.conv.control <- paste("#GCONV_CONFIG\n",
                           "-Startseed=123   -MaxSimGlobalPval=0.05\n",
                           group.string)

write_file(gene.conv.control, "aln/zfy.geneconv.cfg")

#  Only on PATH in Windows, not on cluster yet
if(installr::is.windows()){
  system2("geneconv", paste(nt.aln.file,
                            "aln/zfy.geneconv.cfg",
                            "aln/geneconv.out /lp"),
          stdout = paste0(nt.aln.file, ".geneconv.log"), 
          stderr = paste0(nt.aln.file, ".geneconv.log"))
  
  
  # Read the geneconv output file
  geneconv.data <- read_table("aln/geneconv.frags", comment = "#", 
                              col_names = c("Type", "Pair", "Sim_Pvalue", "KA_pvalue", "begin", "end", 
                                            "len", "num_poly", "num_dif", "tot_difs", "mism_pen")) %>%
    na.omit %>%
    dplyr::mutate(Species = str_replace_all(Pair, "_Z[F|f][Y|y|X|x].*", ""),
                  Type = case_when(Type == "GO" ~ "Global outer",
                                   Type == "GI" ~ "Global inner",
                                   Type == "PO" ~ "Pairwise outer",
                                   Type == "PI" ~ "Pairwise inner",)) %>%
    dplyr::filter(Type == "Pairwise inner")
  
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
  
  geneconv.plot <- ggplot(geneconv.data)+
    geom_rect(data=mouse.exons, aes( xmin = start-0.5, xmax = end+0.5, ymin=0, ymax=0.08, fill=exon), alpha=0.5)+
    geom_text(data=mouse.exons, aes(x=(start+end)/2,y = 0.075, label = exon), size = 3)+
    geom_segment(aes(x=begin, y = KA_pvalue, col=Species, xend = end, yend = KA_pvalue), linewidth = 2)+
    coord_cartesian(xlim = c(1, 2600))+
    scale_fill_manual(values = c("grey", "white", "grey", "white", "grey", "white", "grey"))+
    geom_text(aes(x=(begin+end)/2,y = KA_pvalue+0.002, label = Pair), size = 2)+
    labs(x = "Position", y = "Bonferroni corrected p-value")+
    guides(fill = "none")+
    theme_bw()
  
  save.double.width("figure/zfxy_geneconv.png", geneconv.plot)
  
}

#### Ancestral sequence reconstruction #####

# Ancestral reconstruction
# ancestral.seqs <- read.table(paste0(nt.aln.file, ".state") ,header=TRUE)

#### Exon by exon MSA #####

nt.aln.tidy <- tidy_msa(msa.nt.aln)
for(i in 1:nrow(mouse.exons)){
  start <- mouse.exons$start[i]
  end <- mouse.exons$end[i]
  exon <- mouse.exons$exon[i]
  
  msa.plot <- ggplot()+
    geom_msa(data = nt.aln.tidy, seq_name = T, font=NULL, 
             border=NA, color="Chemistry_NT", consensus_views = T, ref = "Platypus_ZFX", )+
    coord_cartesian(xlim = c(start, end))+
    theme_minimal()+
    theme(axis.text = element_text(size=6),
          axis.title.y = element_blank(),
          panel.grid = element_blank())
  
  ggsave(paste0("figure/msa_exon_", exon, ".png"),  
         msa.plot, dpi = 300, units = "mm", width = 170)
}

#### Plot full MSA #####

msa.plot <- ggplot()+
  geom_msa(data = nt.aln.tidy, seq_name = T, font=NULL, 
           border=NA, color="Chemistry_NT", consensus_views = T, ref = "Platypus_ZFX", )+
  facet_msa(field = 400)+
  theme_minimal()+
  theme(axis.text = element_text(size=2),
        axis.title.y = element_blank(),
        panel.grid = element_blank())

ggsave(paste0("figure/msa.complete.png"),
       msa.plot, dpi = 300, units = "mm", width = 170, height = 170)

#### Test selection globally ####

# ape::dnds(macse.aln) # errors
seqin.aln <- seqinr::read.alignment(nt.aln.file, format = "fasta")
kaks.data <- seqinr::kaks(seqin.aln)
kaks.ratio <- kaks.data$ka / kaks.data$ks
# plot(kaks.ratio)
# Strong purifying selection in all pairs

#### Run codeml to check for site-specific selection ####

# Adapted from Beginner's Guide on the Use of PAML to Detect Positive Selection
# https://academic.oup.com/mbe/article/40/4/msad041/7140562 for details

# Two files need to be uploaded to the cluster for running codeml:
# - alignment
# - tree file (unrooted)
# The treefile created earlier needs node names and branch lengths removing: 
# cat aln/zfxy.nt.aln.treefile | sed -e 's/Node[0-9]\+\/[0-9\.]\+\/[0-9\.]\+:[0-9\.]\+//g'  > paml/site-specific/zfxy.nt.aln.paml.treefile

# cat aln/zfxy.nt.aln.unrooted.treefile | sed -e 's/Node[0-9]\+\/[0-9\.]\+\/[0-9\.]\+:[0-9\.]\+//g' | sed -e 's/:[0-9\.]\+//g'  > paml/site-specific/zfxy.nt.aln.paml.treefile

# sed -e 's/:[0-9\.]\+//g'

labelled.tree <- macse.tree
labelled.tree$node.label <- rep("", length(labelled.tree$node.label))
labelled.tree$edge.length <- rep(1, length(labelled.tree$edge.length))
ape::write.tree(labelled.tree, file = "paml/site-specific/zfxy.nt.aln.paml.treefile")

# Then remove the branch lengths, separate the labels from taxon names and resave
newick.test <- read_file("paml/site-specific/zfxy.nt.aln.paml.treefile")
newick.test <- gsub(":1", "", newick.test)
write_file(newick.test, "paml/site-specific/zfxy.nt.aln.paml.treefile")

# This control file tests site models With heterogeneous ω Across Sites
paml.site.file <- paste0("seqfile   = ", nt.aln.file, " * alignment file\n",
                       "treefile  = paml/site-specific/zfxy.nt.aln.paml.treefile * tree in Newick format without nodes\n", # 
                       "outfile   = paml/site-specific/site.specific.paml.out.txt\n",
                       "\n",
                       "noisy     = 3 * on screen logging\n",
                       "verbose   = 1 * detailed output in file\n",
                       "\n",
                       "seqtype   = 1 * codon data\n",
                       "ndata     = 1 * one gene alignment\n",
                       "icode     = 0 * universal genetic code\n",
                       "cleandata = 0 * keep sites with ambiguity data\n",
                       "\n",
                       "model     = 0 * ω consistent across branches\n",
                       "NSsites   = 0 1 2 7 8 * ω variation across sites\n",
                       "CodonFreq = 7 * use mutation selection model\n",
                       "estFreq   = 0 * use observed frequencies to calc fitness/freq pairs\n",
                       "clock     = 0 * assume no clock\n",
                       "fix_omega = 0 * enables option to estimate omega\n",
                       "omega     = 0.5 * initial omega value\n")
write_file(paml.site.file, "paml/site-specific/zfy.site-specific.paml.ctl")

if(!installr::is.windows()){
  # Run codeml
  system2("codeml", "paml/site-specific/zfy.site-specific.paml.ctl",
          stdout = "paml/site-specific/zfy.site-specific.paml.log", 
          stderr = "paml/site-specific/zfy.site-specific.paml.log")
  
  # Extract lnl
  system2("cat", "zfy.out.txt | grep --before-context=5 'lnL' | grep -e 'lnL' -e 'Model'| paste -d ' '  - - > zfy.out.lnl.txt")
}

# Process the output file to find log likelihood values to calculate LRT
# (likelihood ratio test): twice the difference in log-likelihood ℓ between the
# null and alternative hypotheses, 2Δℓ = 2(ℓ1 − ℓ0), where ℓ0 is the
# log-likelihood score for the null model, whereas ℓ1 is the log-likelihood
# under the alternative model.

# ℓ is in the output file at lines starting lnL
# Grep the lnL and previous 5 lines (which has model name).
# Get just the lnL and model lines from these 5
# Paste alternate lines together with a space
# cat zfy.out.txt | grep --before-context=5 'lnL' | grep -e 'lnL' -e 'Model'| paste -d " "  - -
# Note that there are other lines with 'model' in the file so we still need the first grep

# e.g. values from testing
# Model 0: one-ratio lnL(ntime: 95  np:160): -20797.229748      +0.000000
# Model 1: NearlyNeutral (2 categories) lnL(ntime: 95  np:161): -20712.485759      +0.000000
# Model 2: PositiveSelection (3 categories) lnL(ntime: 95  np:163): -20712.488036      +0.000000
# Model 7: beta (10 categories) lnL(ntime: 95  np:161): -20614.279484      +0.000000
# Model 8: Model 8: beta&w>1 (11 categories) lnL(ntime: 95  np:163): -20614.287541      +0.000000

# Calculate the liklihood ratio test for two models
# lnl - the log likelihoods
# np - the number of free parameters
# This calculates the LRT and tests it against the chi-distribution where the
# degrees of freedom are the difference in the number of free parameters between
# the models. 
calc.LRT <- function(lnl0, lnl1, np0, np1){
  lrt <- 2 * (lnl1-lnl0)
  df <- abs(np1-np0)
  crit.value =  qchisq(p=0.05, df=df, lower.tail = FALSE)
  list("crit.value" = crit.value,
       "p.value" = pchisq(lrt, df, lower.tail = FALSE),
       "lrt" = lrt)
}

# M0 vs. M1a (one-ratio vs. nearly neutral)

# This is a test for variability of selective pressure among amino acid sites
# rather than a test of positive selection. M1a fits the data much better than
# M0,indicating that the selective pressure reflected by ω
# varies hugely among sites.
m0m1a <- calc.LRT(-20797.229748, -20712.485759, 160, 161)

# Compared with M1a, M2a adds a class of sites under positive selection with ω2
# > 1 (in proportion p2). This does not improve the fit of the model
# significantly
m1am2a <- calc.LRT(-20712.485759, -20712.488036, 161, 163) # (nearly neutral vs. positive selection)

# Additional test for positive selection by comparing M7 (beta, null model)
# against M8 (beta&ω, alternative model).
m7m8 <- calc.LRT(-20614.279484, -20614.287541, 161, 163) # (positive selection vs null model)
# No evidence for sites under positive selection 

#### Run codeml to look for selection specifically in rodents after beaver ####

# To look at the rodent clade, we need a rooted tree

# Find the MRCA of the rodents with rapid ZFY evolution
rodent.node <- ape::getMRCA(macse.tree, c("Mouse_Zfy2", "Desert_hamster_Zfx-like_putative-Zfy"))

# Remove existing node labels, label the nodes and tips of the tree with #1
# for foreground branches
labelled.tree <- macse.tree
labelled.tree$node.label <- rep("", length(labelled.tree$node.label))
labelled.tree$edge.length <- rep(1, length(labelled.tree$edge.length))
labelled.tree <- treeio::label_branch_paml(labelled.tree, rodent.node, "#1")
ape::write.tree(labelled.tree, file = "paml/site-branch/zfxy.nt.aln.paml.treefile")

# Then remove the branch lengths, separate the labels from taxon names and resave
newick.test <- read_file("paml/site-branch/zfxy.nt.aln.paml.treefile")
newick.test <- gsub(":1", "", newick.test) # branch lengths we set to 1
newick.test <- gsub("_#1", " #1", newick.test) # fg labels
write_file(newick.test, "paml/site-branch/zfxy.nt.aln.paml.fg.treefile")

# This control file tests site models With heterogeneous ω Across Sites
paml.site.branch.file <- paste0("seqfile   = ", nt.aln.file, " * alignment file\n",
                                "treefile  = paml/site-branch/zfxy.nt.aln.paml.fg.treefile * tree in Newick format without nodes\n", # 
                                "outfile   = paml/site-branch/site.branch.paml.out.txt\n",
                                "\n",
                                "noisy     = 3 * on screen logging\n",
                                "verbose   = 1 * detailed output in file\n",
                                "\n",
                                "seqtype   = 1 * codon data\n",
                                "ndata     = 1 * one gene alignment\n",
                                "icode     = 0 * universal genetic code\n",
                                "cleandata = 0 * keep sites with ambiguity data\n",
                                "\n",
                                "model     = 2 * 2 ω values across branches\n",
                                "NSsites   = 2 * Model M2a\n",
                                "CodonFreq = 7 * use mutation selection model\n",
                                "estFreq   = 0 * use observed frequencies to calc fitness/freq pairs\n",
                                "clock     = 0 * assume no clock\n",
                                "fix_omega = 0 * enables option to estimate omega\n",
                                "omega     = 0.5 * initial omega value\n")
write_file(paml.site.branch.file, "paml/site-branch/zfy.site-branch.paml.ctl")


#### Extend the AA alignments with non-mammalian outgroups ####

# Read all unaligned sequence files with .fa extension
outgroup.aa.files <- list.files(path = "fasta/aa", pattern = "*.fa$", 
                       include.dirs = T, full.names = T)

outgroup.fa.read <- lapply(outgroup.aa.files, read.fasta)

outgroup.nt.raw <- do.call(c, lapply(outgroup.fa.read, function(x) x$fa))
outgroup.metadata <- do.call(rbind, lapply(outgroup.fa.read, function(x) x$metadata))

# Find the nodes that are ZFY vs ZFX to add to tree
outgroup.metadata %<>%
  dplyr::mutate(group = case_when(grepl("ZFY", common.name, ignore.case=T) ~ "ZFY",
                                  common.name %in% c(outgroup.metadata$common.name, "Platypus_ZFX", "Opossum_ZFX") ~ "Outgroup",
                                  T ~ "ZFX"))

combined.metadata <- rbind(metadata, outgroup.metadata)

# Write the combined fasta to file with .fas extension
outgroup.nt.file <- "fasta/outgroup.nt.fas"
ape::write.FASTA(outgroup.nt.raw, file = outgroup.nt.file)

# Use macse to align the outgroup NT sequences
system2("java", paste("-jar bin/macse_v2.07.jar -prog alignSequences",
                      "-seq ", outgroup.nt.file,
                      "-out_NT", "aln/outgroup/outgroup.nt.aln",
                      "-out_AA", "aln/outgroup/outgroup.aa.aln"), 
        stdout = paste0("aln/outgroup/outgroup.macse.log"),  # logs
        stderr = paste0("aln/outgroup/outgroup.macse.log"))  # error logs

# Use macse to align the outgroup NT sequences with the existing NT alignment
system2("java", paste("-jar bin/macse_v2.07.jar -prog alignTwoProfiles ",
                      "-p1", nt.aln.file,
                      "-p2", "aln/outgroup/outgroup.nt.aln",
                      "-out_NT", "aln/outgroup/combined.nt.aln",
                      "-out_AA", "aln/outgroup/combined.aa.aln"),
        stdout = paste0("aln/combined.macse.log"),  # logs
        stderr = paste0("aln/combined.macse.log"))  # error logs


#### View the extended outgroup AA MSA ####

msa.outgroup.aln <- Biostrings::readAAMultipleAlignment("aln/outgroup/combined.aa.aln", format="fasta")
msa.outgroup.aln.tidy <- tidy_msa(msa.outgroup.aln)

ggplot()+
  geom_msa(data = msa.outgroup.aln.tidy, seq_name = T, font=NULL, border=NA,
           consensus_views = T, ref = "Platypus_ZFX", alpha = 0.5
  )+
  labs(x = "Amino acid")+
  theme_minimal()+
  theme(axis.text = element_text(size=6),
        axis.title.y = element_blank(),
        panel.grid = element_blank())

#### Make a tree with the outgroups MSA using AA only ####

system2("iqtree", paste("-s ", "aln/outgroup/combined.aa.aln", 
                        "-bb 1000", # number of bootstrap replicates
                        "-alrt 1000", # number of replicates to perform SH-like approximate likelihood ratio test (SH-aLRT) 
                        "-nt AUTO"), # number of threads
        stdout = paste0("aln/outgroup/combined.aa.iqtree.log"), 
        stderr = paste0("aln/outgroup/combined.aa.iqtree.log"))


outgroup.tree <- ape::read.tree(paste0("aln/outgroup/combined.aa.aln.treefile"))

# Root the tree on Mouse_Zfp711 and resave
outgroup.tree <- ape::root(outgroup.tree, outgroup = c("Mouse_Zfp711"), resolve.root = TRUE)
ape::write.tree(outgroup.tree, file = paste0("aln/outgroup/combined.aa.aln.rooted.treefile"))

group_info <- split(combined.metadata$common.name, combined.metadata$group)
outgroup.tree <- groupOTU(outgroup.tree, group_info, group_name = "group")

outgroup.plot <- plot.tree(outgroup.tree)+
  coord_cartesian(clip="off", xlim = c(0, 0.9), ylim= c(-2, 62))

# outgroup.plot <- ggtree(outgroup.tree) + 
#   geom_tree() +
#   theme_tree() +
#   geom_tiplab(size=3, aes(col = group))+
#   geom_treescale(y=-1)+
#   coord_cartesian(clip="off", xlim = c(0, 0.9), ylim= c(-1, 60))+
#   theme(legend.position = "none")

save.double.width("figure/outgroup.tree.png", outgroup.plot)



#### What are the binding motifs of the ZFs in each species? ####

# predict the ZF targets in each sequence via pwm_predict (http://zf.princeton.edu/download2.php)

# pwm_predict zfy.aa.aln

# Remove header and split the output PWMs to separate files
# cat zfxy.aa.pwm | grep -v ';' | split -l 5 - zf_

# Read each file, get headers and PWMs
pwm.files <- list.files("pwm", pattern = "zf_", full.names = T)

read.pwm <- function(f){
  pwm <- read.table(f, skip=1)
  header <- read.table(f, nrows=1)
  p <- seqLogo::makePWM(pwm)
  consensus <- seqLogo::consensus(p)
  data.frame("species"= header$V1,
       "zf" = header$V2,
       "consensus" = consensus)
}

pwm.predictions <- do.call(rbind, lapply(pwm.files, read.pwm)) %>%
  tidyr::pivot_wider(id_cols = species, names_from = zf, values_from = consensus)

#### Annotate the locations of the ZFs in the AA MSA ####

# We can use the regex .C.{2,4}C.{12}H.{3,5}H to find ZFs as per https://pubmed.ncbi.nlm.nih.gov/21572177/

msa.aa.aln <- Biostrings::readAAMultipleAlignment(aa.aln.file, format="fasta")

zf.regex <- ".C(.{2,4}?)C.{12}H(.{3,5}?)H"

find.zf <- function(aa, sequence.name){
  hits <- as.data.frame(str_locate_all(as.character(aa), zf.regex))
  hits$species <- sequence.name
  hits
}

zf.locations <- do.call(rbind, mapply(find.zf, aa=msa.aa.aln@unmasked, sequence.name = names(msa.aa.aln@unmasked), SIMPLIFY = FALSE))

msa.aa.aln.tidy <- tidy_msa(msa.aa.aln)
# Reverse the order of rows to match the y order of the MSA
msa.aa.aln.tidy.order <- msa.aa.aln.tidy %>% dplyr::filter(position==1) %>% dplyr::mutate(i = n()+1-row_number())
zf.locations <- merge(zf.locations, msa.aa.aln.tidy.order[,c(1, 4)], by.x="species", by.y="name")

# Calculate the midpoint of each ZF for labelling, smooth out species with
# different start points
zf.labels <- zf.locations %>% 
  dplyr::mutate(mid = round( (start+end)/2, digits = 0),
                mid = ifelse(mid%%2==0, mid, mid+1)
                ) %>%
  dplyr::group_by(mid) %>%
  summarise() %>%
  dplyr::mutate(label = row_number())

make.aa.msa <- function(start, end){
  ggplot()+
    geom_msa(data = msa.aa.aln.tidy, seq_name = T, font=NULL, border=NA,
             consensus_views = T, ref = "Platypus_ZFX", alpha = 0.5
    )+
    geom_rect(data = zf.locations, aes(xmin=start, xmax=end, ymin=i-0.5, ymax=i+0.5), fill="grey", alpha=0.5)+
    geom_text(data = zf.labels, aes(x = mid, y = 57, label = label))+
    coord_cartesian(xlim = c(start, end),
                    ylim = c(0, 57))+
    labs(x = "Amino acid")+
    theme_minimal()+
    theme(axis.text = element_text(size=6),
          axis.title.y = element_blank(),
          panel.grid = element_blank())
}

# too wide to display in one figure, split in ~middle
aa.msa.plot.1 <- make.aa.msa(min(zf.locations$start)-10, 650)
aa.msa.plot.2 <- make.aa.msa(640, max(zf.locations$end)+10)

save.double.width("figure/aa.msa.1.png", aa.msa.plot.1)
save.double.width("figure/aa.msa.2.png", aa.msa.plot.2)

#### Get divergence times to highlight the rapid evolution in the rodents ####

# Get the NCBI taxon ids for each species
taxon.data <- lapply(metadata$species, function(x) httr::GET( paste0("http://timetree.temple.edu/api/taxon/",curl::curl_escape(x)))  ) 
taxon.ids <- sapply(lapply(taxon.data, httr::content), function(x) x$taxon_id)
metadata$taxon_id <- taxon.ids

# Find the pairwise distances between each species
pairwise.species <- expand.grid(unique(metadata$taxon_id), unique(metadata$taxon_id)) %>%
  dplyr::filter(Var1!=Var2, Var1<Var2)

get.time.tree <- function(tax.a, tax.b){
  tryCatch({
    tt <- httr::GET(paste0("http://timetree.temple.edu/api/pairwise/",tax.a, "/", tax.b))
    ct <- content(tt)
    Sys.sleep(0.5)
    
    
    data <- as.data.frame( str_split(ct, "\r\n"), col.names = c("V1")) %>% 
      dplyr::slice_tail(n=1) %>%
      tidyr::separate_wider_delim( cols = V1, delim = ",", names = c("taxon_a_id","taxon_b_id","scientific_name_a","scientific_name_b","all_total","precomputed_age","precomputed_ci_low","precomputed_ci_high","adjusted_age"))
    return(data)
  }, error = function(e) { 
    print(e)
    return(data.frame("taxon_a_id" = c(),"taxon_b_id" = c(),"scientific_name_a" = c(),"scientific_name_b" = c(),"all_total" = c(),"precomputed_age" = c(),"precomputed_ci_low" = c(),"precomputed_ci_high" = c(),"adjusted_age" = c()))
  }  ) 
}

# Save to avoid repeated API calls
if(!file.exists( "time.tree.data.tsv")){
  pairwise.times <- do.call(rbind, mapply(get.time.tree, pairwise.species$Var1, pairwise.species$Var2))
  write_tsv(pairwise.times, file = "time.tree.data.tsv", col_names = T, quote = "none")
} else{
  pairwise.times <- read_tsv("time.tree.data.tsv")
}

#### ####


# testing of side-by-side plots
p <- ggtree(macse.tree) + geom_tiplab(size=3, aes(col = group))
n <- msaplot(p, aa.aln.file, offset=0.4, width=2)+theme(legend.position = "none")

