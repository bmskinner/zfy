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

# Clear previous analyses
filesstrings::dir.remove("aln")
filesstrings::dir.remove("figure")

# Create missing dirs if needed
if(!dir.exists("aln")) dir.create("aln")
if(!dir.exists("aln/exons")) dir.create("aln/exons")
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
  str_replace_all("\\]", "") %>%
  str_replace_all(" ", "_")
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

# Convert to interleaved for display
write.dna(macse.aln, file="aln/nt.interleaved.fa", format = "interleaved")

nt.phylip.phy.file <- "aln/nt.phylip.phy"
write.phyDat(macse.aln, file=nt.phylip.phy.file, format = "phylip")

# Exon by exon coordinates of the alignment will be needed for clear
# testing of selection. Match these from the final alignment via mouse Zfy1
find.exons <- function(){
  macse.nt.aln <- seqinr::read.alignment(nt.aln.file, format="fasta")
  mouse.zfy1 <- toupper(macse.nt.aln$seq[macse.nt.aln$nam=="Mouse_Zfy1"])
  
  mouse.exons <- data.frame("exon" = c("3", "4", "5",  "6", "7", "8",  "9"),
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
aa.clustal <- paste(printMultipleAlignment(macse.aa.aln, names = common.names), collapse = "")
write_file(aa.clustal, file="aln/aa.clustal.aln")

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
ggplot(tidy.msa)+
  geom_rect(data=mouse.exons, aes( xmin = start-0.5, xmax = end+0.5, ymin=0, ymax=1, fill=exon), alpha=0.5)+
  geom_line( aes(x=position, y=fraction))+
  labs(x = "Position in alignement", y = "Fraction of species with consensus nucleotide")+
  theme_bw()


macse.aa.plot <- ggmsa(macse.aa.aln, start = 1, end=860, seq_name = T, 
      border = NA, font=NULL, color="Chemistry_AA") + 
  facet_msa(field = 100)+
  theme(axis.text.y = element_text(size = 2))
ggsave(paste0("figure/msa.aa.complete.png"),  
       macse.aa.plot, dpi = 300, units = "mm", width = 170, height = 370)

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
  ggtree(tree.data) + 
    geom_tree() +
    theme_tree() +
    geom_tiplab(aes(col = ZFY), size=3) +
    # geom_nodelab()+
    geom_treescale()+
    coord_cartesian(clip="off")+
    xlim(0, 0.6) +
    theme(legend.position = "none")
}

macse.tree <- ape::read.tree(paste0(nt.aln.file, ".treefile"))
# Root the tree on platypus
macse.tree <- ape::root(macse.tree, "Platypus_ZFX")

# Find the nodes that are ZFY vs ZFX and add to tree
zfy.nodes <- gsub(" ", "_", common.names[str_detect(common.names, "Z[F|f][Y|y|a]")])
zfy.nodes <- gsub("\\?", "_", zfy.nodes)
macse.tree <- groupOTU(macse.tree, zfy.nodes, group_name = "ZFY")

plot.zfx.zfy <- plot.tree(macse.tree)
save.double.width("figure/zfx.zfy.tree.png", plot.zfx.zfy)


# Drop the ZFY and Zfa sequences and just look at the ZFX nodes in the tree
tree.zfx <- ape::drop.tip(macse.tree, zfy.nodes)
ape::write.tree(tree.zfx, file = paste0(nt.aln.file, ".zfx.treefile"))
plot.zfx <- plot.tree(tree.zfx)
save.double.width("figure/zfx.tree.png", plot.zfx)

# Keep Zfy and Zfa, drop Zfx nodes. Keep outgroups
tree.zfy <- ape::keep.tip(macse.tree, c(zfy.nodes, "Platypus_ZFX", "Opossum_ZFX"))
ape::write.tree(tree.zfy, file = paste0(nt.aln.file, ".zfy.treefile"))

plot.zfy <- plot.tree(tree.zfy)
save.double.width("figure/zfy.tree.png", plot.zfy)

#### Make tree from each exon separately ####

# Aim here is to look for signs of gene conversion as per other papers

for(i in 1:nrow(mouse.exons)){

  exon.aln <- as.matrix(macse.aln)[,mouse.exons$start[i]:mouse.exons$end[i]]
  exon.aln.file <- paste0("aln/exon/exon_", mouse.exons$exon[i], ".aln")
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
  # Root the tree on platypus
  exon.tree <- ape::root(exon.tree, "Platypus_ZFX")
  
  # Find the nodes that are ZFY vs ZFX and add to tree
  zfy.nodes <- gsub(" ", "_", common.names[str_detect(common.names, "Z[F|f][Y|y|a]")])
  zfy.nodes <- gsub("\\?", "_", zfy.nodes)
  exon.tree <- groupOTU(exon.tree, zfy.nodes, group_name = "ZFY")
  
  plot.exon.tree <- plot.tree(exon.tree)
  exon.fig.file <- paste0("figure/exon_", mouse.exons$exon[i], ".zfx.zfy.tree.png")
  save.double.width(exon.fig.file, plot.exon.tree)
  
}

#### Ancestral sequence reconstruction #####

# Ancestral reconstruction
# ancestral.seqs <- read.table(paste0(nt.aln.file, ".state") ,header=TRUE)

#### Exon by exon MSA #####

for(i in 1:nrow(mouse.exons)){
  start <- mouse.exons$start[i]
  end <- mouse.exons$end[i]
  exon <- mouse.exons$exon[i]
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

#### Test selection exon by exon ####

# ape::dnds(macse.aln) # errors
seqin.aln <- seqinr::read.alignment(nt.aln.file, format = "fasta")
kaks.data <- seqinr::kaks(seqin.aln)
kaks.ratio <- kaks.data$ka / kaks.data$ks
plot(kaks.ratio)
# Strong purifying selection in all pairs



#### Create partition model for IQ-TREE #####

