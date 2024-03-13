# Analysis pipeline for ZFX/ZFY evolutionary analysis

# The following binaries must be supplied in ./bin:
# macse_v2.07.jar
# muscle5.1.win64.exe / muscle5.1.linux_intelx64
# nlstradamul.pl
# pwm_predict/pwm_predict

# Other binaries are expected on the PATH:
# PAML
# IQ-TREE
# GENECONV

# HyPhy is a python package installed into a conda environment:
# conda create -n hyphy # create conda environment
# conda activate hyphy # enter the environment
# conda install -c bioconda hyphy # install hyphy from bioconda
# To activate conda from a shell script:
# source activate hyphy

RUN_PAML = as.logical(commandArgs(trailingOnly = T)[1])

if(is.na(RUN_PAML)) RUN_PAML <- FALSE

cat("Run PAML is", RUN_PAML, "\n")
#### Imports #####

install.cran <- function(package){
  if(!require(package, character.only = TRUE, quietly = TRUE)){
    install.packages(package, repos = "https://cran.ma.imperial.ac.uk")
    library(package, character.only = TRUE, quietly = TRUE)
  }
}

install.bioconductor <- function(package){
  if(!require(package, character.only = TRUE, quietly = TRUE)){
    BiocManager::install(package, update = FALSE)
    library(package, character.only = TRUE, quietly = TRUE)
  }
}

install.github <- function(package){
  pkg.name <- gsub("^.*\\/", "", package)
  if(!require(pkg.name, character.only = TRUE, quietly = TRUE)){
    remotes::install_github(package)
    library(pkg.name, character.only = TRUE, quietly = TRUE)
  }
}

cran.packages <- c("tidyverse", "ape", "filesstrings", "seqinr", "phangorn",
                   "installr","treespace", "httr", "seqLogo", "assertthat", "aplot",
                   "paletteer", "ggnewscale", "slider", "BiocManager",
                   "remotes", "patchwork")

sapply(cran.packages, install.cran)

github.packages <- c('YuLab-SMU/ggtree', # since ggtree cannot install in Bioconductor 3.15 on cluster
                     "vragh/seqvisr", "vmikk/metagMisc")
sapply(github.packages, install.github)

bioconductor.packages <- c("msa", "ggmsa", "treeio")
sapply(bioconductor.packages, install.bioconductor)

# library(tidyverse) # CRAN
# library(ape) # CRAN
# library(msa) #  BiocManager::install("msa")
# library(filesstrings)
# library(ggmsa) # BiocManager
# library(ggtree) # BiocManager
# library(seqinr)
# library(phangorn)
# library(installr)
# library(treeio) # BiocManager
# library(httr)
# library(seqLogo)
# library(seqvisr) # remotes::install_github("vragh/seqvisr")
# library(metagMisc) # remotes::install_github("vmikk/metagMisc")
# library(assertthat)
# library(aplot)
# library(treespace)
# library(paletteer)
# library(ggnewscale) 
# library(slider)

source("find9aaTADs.R")
source("findZF.R")
source("calcCharge.R")
source("calcHydrophobicity.R")

cat("Packages loaded\n")

#### Common functions ####

save.double.width <- function(filename, plot, height=170) ggsave(filename, plot, dpi = 300, 
                                                                 units = "mm", width = 170, 
                                                                 height = height)

# Translate ungapped coordinates back to gapped
# site.no.gap - the integer site in an ungapped sequence to convert
# gapped.seq - the sequence with gaps from an alignment
convert.to.gapped.coordinate <- function(site.no.gap, gapped.seq){

  if(is.na(site.no.gap)) return(NA)
  gapped.seq.char <- as.character(gapped.seq)
  # find gaps and stop codons
  gaps <- str_locate_all(gapped.seq.char, "-|\\*")[[1]][,1]
  n <- site.no.gap
  for(i in gaps){
    if(i<n) n <- n + 1
  }
  n
}

# Exon by exon coordinates of the alignment will be needed for clear
# testing of selection. Match these in the final alignment via mouse Zfy1
# biostrings.nt.alignment - an MSA from Biostrings::readDNAMultipleAlignment
# biostrings.aa.alignment - an MSA from Biostrings::readAAMultipleAlignment
find.exons <- function(biostrings.nt.alignment, biostrings.aa.alignment){
  
  mouse.zfy1.nt <- as.character(biostrings.nt.alignment@unmasked$Mouse_Zfy1)
  mouse.zfy1.nt.ungapped <- str_remove_all(mouse.zfy1.nt, "-|\\*")
  
  mouse.zfy1.aa <- as.character(biostrings.aa.alignment@unmasked$Mouse_Zfy1)
  mouse.zfy1.aa.ungapped <- str_remove_all(mouse.zfy1.aa, "-|\\*")
  
  mouse.exons <- data.frame("exon"     = c("1", "2", "3",  "4", "5", "6",  "7"),
                            "start_nt" = c("ATGGATGAA", "GAGCTGATGCA", "TGGATGAACC", 
                                           "GAGAAACTAT", "AAGTAATTGT", "ATAATAATTCT", "CAATATTTGTT"),
                            "end_nt"   = c("TGGAATAG", "ATGATGTCTT", "GGATGAATTAG", 
                                         "GAAGAAGATACTG", "GACAGCAGCTTATG", "CAGTACCAGTCAG", "CCTGCCC"),
                            "start_aa" = c("MDEDEIEL", "GADAVHMD", "LDEPSKADL", "LGETIHAVE",
                                           "EVIVGDED", "DNNSDEIE", "AIFVAPDGQ"),
                            "end_aa"   = c("KSFFDGIG", "INCEDYLMMSL","ADSEVDEL", "SQKEEEDTE",
                                         "PIAWTAAYD", "PESKQYQSA", "RHHKVGLP"))
  
  start.nt <- sapply(mouse.exons$start_nt, str_locate, string=mouse.zfy1.nt.ungapped)[1,]
  end.nt   <- sapply(mouse.exons$end_nt,   str_locate, string=mouse.zfy1.nt.ungapped)[2,]
  
  start.aa <- sapply(mouse.exons$start_aa, str_locate, string=mouse.zfy1.aa.ungapped)[1,]
  end.aa   <- sapply(mouse.exons$end_aa,   str_locate, string=mouse.zfy1.aa.ungapped)[2,]
  
  # Note that the aa and nt positions are not from equivalent alignments!
  # NT is from mammals, AA is from mammals + outgroups
  data.frame("exon"     = mouse.exons$exon,
             "start"    = sapply(start.nt, convert.to.gapped.coordinate, mouse.zfy1.nt),
             "end"      = sapply(end.nt,   convert.to.gapped.coordinate, mouse.zfy1.nt),
             "start_aa" = sapply(start.aa, convert.to.gapped.coordinate, mouse.zfy1.aa),
             "end_aa"   = sapply(end.aa,   convert.to.gapped.coordinate, mouse.zfy1.aa),
             "is_even"  = sapply(mouse.exons$exon, function(i) as.numeric(i)%%2==0))%>%
    dplyr::mutate(length_nt = end - start,
                   # fix the offsets to get ORF
                  start_nt_codon_offset = case_when(exon==1 ~ start, # hardcode the codon offsets for subsetting
                                                    exon>=2 ~ start-1),
                  end_nt_codon_offset   = case_when(exon<7 ~ end-1,
                                                    exon==7 ~ end),
                  corrected_offset_length = (end_nt_codon_offset - start_nt_codon_offset)%%3) # check divisible by 3 
}

plot.tree <- function(tree.data, ...){
  # Get the complete node labels
  # Separate out bootstrap info
  # Numbers in parentheses are SH-aLRT support (%) / ultrafast bootstrap support (%)
  node.label.values <- data.frame("label" = tree.data$node.label) %>%
    tidyr::separate_wider_delim(label, delim = "/", names = c("name","SHaLRT", "UFBoot"),
                                too_few = "align_end") %>%
    dplyr::mutate(UFBoot = as.numeric(UFBoot),
                  SHaLRT = as.numeric(SHaLRT),
                  isSupportedUFBoot = UFBoot>=95 & !is.na(UFBoot),
                  isSupportedSHalRT = SHaLRT>=80 & !is.na(SHaLRT),
                  colour = case_when(isSupportedUFBoot & isSupportedSHalRT ~ "black",
                                     isSupportedUFBoot |isSupportedSHalRT ~ "grey",
                                     .default = "white"))
  ggtree(tree.data) + 
    geom_tree() +
    geom_tiplab(size=2, aes_string(...))+
    # geom_nodelab(size=2, nudge_x = -0.003, nudge_y = 0.5, hjust=1,  node = "internal")+
    geom_nodepoint(size=1.5,  col="black")+
    geom_nodepoint(size=0.75,  col=node.label.values$colour)+
    geom_treescale(fontsize =2, y = -1) +
    coord_cartesian(clip="off")+
    theme_tree() +
    theme(legend.position = "none")
}

# Create a line-only tree scaled to the given max y. This allows the tree
# to be combined with other plots with larger max y values.
make.outgroup.mini.tree <-function(combined.outgroup.tree, text.labels){
  combined.outgroup.tree.mini <- combined.outgroup.tree
  
  # How many labels do we need above the tree?
  # 2 rows per label, plus one spacer
  max.y <- length(combined.outgroup.tree.mini$tip.label) + (length(text.labels)*3)

  
  # Replace the tip labels with empty string so no labels are plotted
  combined.outgroup.tree.mini$tip.label <- rep("", length(combined.outgroup.tree.mini$tip.label))
  result <- ggtree(combined.outgroup.tree.mini, aes(color=group)) + 
    geom_tiplab(align=TRUE, linetype='dashed', linesize=.3) + # use tiplab to get lines
    coord_cartesian(ylim = c(1, max.y))
  
  curr.y <- length(combined.outgroup.tree.mini$tip.label) + 2.5
  for(i in 1:length(text.labels)){
    result <- result + annotate(geom="text", x=0.6, y=curr.y, 
                                label=text.labels[i], size=2, hjust=1, vjust=0.5)
    curr.y <- curr.y + 3
  }
    # 
    # 
    # annotate(geom="text", x=-20, y=n.taxa+4, label="Chicken", size=1)+
    # annotate(geom="text", x=-20, y=n.taxa+7, label="Opossum", size=1)+
    
  result + theme(legend.position = "none",
          plot.margin = margin(r=0))
}

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

# Add a rectangle track to a plot
add.track <- function(ranges, y.start, y.end, ...){
  geom_rect(data=ranges, 
            aes(xmin = start, xmax = end, ymin = y.start, ymax=y.end),
            ...)
}
# Add a rectangle track labels to a plot
add.track.labels <- function(ranges, y.start, y.end, ...){
  geom_text(data=ranges,
            aes(x =(start+end)/2, y=(y.start+y.end)/2, label=label), size=1.8, col="white")
}
# Add a conservation track to a plot
add.conservation.track <- function(ranges,  y.start, y.end, ...){
  geom_rect(data=ranges,  aes(xmin=position-0.45, 
                              xmax=position+0.45, 
                              ymin=y.start, 
                              ymax=y.end, 
                              fill=smoothed9))
}

# Add an exon track to a plot
add.exon.track <- function(y.start, y.end, ...){
  geom_rect(data = mouse.exons, aes(xmin = start_aa, xmax = end_aa, 
                                    ymin = y.start, ymax = y.end, 
                                    fill = exon),...)
}
# Add exon labels track to a plot
add.exon.labels <- function(y.start, y.end, ...){
  geom_text(data=mouse.exons,
            aes(x =(start_aa+end_aa)/2, y=(y.start+y.end)/2, label=exon), size=1.8, col="black")
}

# Annotate Zfy structures to a plot
# plot - the plot to annotate
# n.taxa - the number of rows within data; annotation tracks are above this
annotate.structure.plot <- function(plot, n.taxa){
  
  plot <- plot+
    # Draw the conservation with Xenopus, chicken and opossum
    new_scale_fill()+
    scale_fill_viridis_c(limits = c(0, 1))+
    labs(fill="Conservation (9-window smooth)")+
    add.conservation.track(msa.aa.aln.tidy.frog.conservation,    n.taxa,   n.taxa+2)+
    add.conservation.track(msa.aa.aln.tidy.chicken.conservation, n.taxa+3, n.taxa+5)+
    add.conservation.track(msa.aa.aln.tidy.opossum.conservation, n.taxa+6, n.taxa+8)+
    
    # Draw the structures
    add.track(ranges.ZF.common,     n.taxa+9, n.taxa+11, fill="lightgrey")+
    add.track(ranges.NLS.common,    n.taxa+9, n.taxa+11, fill="green", alpha = 0.5)+
    add.track(ranges.9aaTAD.common, n.taxa+9, n.taxa+11, fill="blue",  alpha = 0.5)+
    add.track.labels(ranges.9aaTAD.common, n.taxa+9, n.taxa+11)+   # Label the 9aaTADs
    
    new_scale_fill()+
    scale_fill_manual(values=c("grey", "salmon", "white", "grey", "white", "grey", "white"))+
    guides(fill = "none")+
    add.exon.track(n.taxa+12, n.taxa+14, col = "black")+
    add.exon.labels(n.taxa+12, n.taxa+14)+
    
    scale_x_continuous(expand = c(0, 0))+
    coord_cartesian(xlim = c(0, max(msa.aa.aln.tidy.frog.conservation$position)))+
    theme_bw()+
    theme(axis.text.y = element_blank(),
          axis.title = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size=6),
          legend.position = "top",
          legend.title = element_text(size = 6, vjust = 0.7),
          legend.text = element_text(size = 6),
          panel.grid = element_blank())
  
  # Add the tree with the outgroups
  make.outgroup.mini.tree(combined.outgroup.tree,
                          text.labels = c("Xenopus", "Chicken", "Opossum", "Features", "Exons"))+
    plot +
    patchwork::plot_layout(widths = c(0.1, 0.9))
}


#### Create output directory structure #####

# Clear previous analyses
filesstrings::dir.remove("aln")
filesstrings::dir.remove("figure")
filesstrings::dir.remove("nls")
filesstrings::dir.remove("pwm")

# Create missing dirs if needed
filesstrings::create_dir("aln")
filesstrings::create_dir("aln/outgroup")
filesstrings::create_dir("aln/exons")
filesstrings::create_dir("aln/zfx_only")
filesstrings::create_dir("aln/zfy_only")
filesstrings::create_dir("bin")
filesstrings::create_dir("figure")
filesstrings::create_dir("hyphy")
filesstrings::create_dir("nls")
filesstrings::create_dir("paml")
filesstrings::create_dir("paml/site-specific")
filesstrings::create_dir("paml/branch-site")
filesstrings::create_dir("paml/exon_1_3-6")
filesstrings::create_dir("paml/exon_2")
filesstrings::create_dir("paml/exon_7")
filesstrings::create_dir("pwm")

writeLines(capture.output(sessionInfo()), "figure/session_info.txt")

#### Read mammal NT FA files #####

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
mammal.nt.file <- "fasta/mammal.nt.fas"
ape::write.FASTA(nt.raw, file = mammal.nt.file)

#### Read outgroup NT FA files #####

# Read all unaligned sequence files with .fa extension
outgroup.nt.files <- list.files(path = "fasta/aa", pattern = "*.fa$", 
                                include.dirs = T, full.names = T)

outgroup.fa.read <- lapply(outgroup.nt.files, read.fasta)

outgroup.nt.raw <- do.call(c, lapply(outgroup.fa.read, function(x) x$fa))
outgroup.metadata <- do.call(rbind, lapply(outgroup.fa.read, function(x) x$metadata))

# Find the nodes that are ZFY, ZFX or outgroup
outgroup.metadata %<>%
  dplyr::mutate(group = case_when(grepl("ZFY", common.name, ignore.case=T) ~ "ZFY",
                                  common.name %in% c(outgroup.metadata$common.name, 
                                                     "Platypus_ZFX", 
                                                     "Opossum_ZFX") ~ "Outgroup",
                                  T ~ "ZFX"))

# Combine the outgroups with the mammals
combined.metadata <- rbind(metadata, outgroup.metadata)

# Write the unaligned combined fasta to file with .fas extension
combined.nt.raw <- c(outgroup.nt.raw, nt.raw) # all sequences
combined.nt.file <- "fasta/combined.nt.fas"
ape::write.FASTA(combined.nt.raw, file = combined.nt.file)

combined.aa.file <- "fasta/combined.aa.fas"
combined.aa.raw <- ape::trans(combined.nt.raw)
ape::write.FASTA(combined.aa.raw, file = combined.aa.file)

# Create supplementary table with all accessions and sequence info
supplementary.accessions.table <- combined.metadata %>%
  dplyr::rename(Accession = accession, Description = original.name,
                Name_in_figures = common.name, Species = species, Group = group)
write_tsv(supplementary.accessions.table, "figure/accessions.supplement.tsv")

#### Run mammal NT alignment guided by AA #####

# Run a codon aware alignment with MACSE
# Expect java on the PATH. Macse download from https://www.agap-ge2pop.org/macsee-pipelines/

# Direct download link:
# https://www.agap-ge2pop.org/wp-content/uploads/macse/releases/macse_v2.07.jar
if(!file.exists("bin/macse_v2.07.jar")) stop("MACSE not present in bin directory")

# Run the alignment
nt.aln.file <- "aln/zfxy.nt.aln"
aa.aln.file <- "aln/zfxy.aa.aln"
system2("java", paste("-jar bin/macse_v2.07.jar -prog alignSequences",
                      "-seq", mammal.nt.file, # input
                      "-out_NT", nt.aln.file,  # output nt alignment
                      "-out_AA", aa.aln.file), # output aa alignment
        stdout = paste0(nt.aln.file, ".macse.log"),  # logs
        stderr = paste0(nt.aln.file, ".macse.err"))  # error logs
ape.nt.aln <- ape::read.FASTA(nt.aln.file)

# Check the AA alignment
# ape.aa.aln <- seqinr::read.alignment(aa.aln.file, format="fasta")
# aa.clustal <- paste(printMultipleAlignment(ape.aa.aln, names = common.names), collapse = "")
# write_file(aa.clustal, file="aln/aa.clustal.aln")

# Calculate conservation at each site
msa.nt.aln <- Biostrings::readDNAMultipleAlignment(nt.aln.file, format="fasta")

#### Run extended outgroup AA alignment ####

# Use muscle to align the aa files and save out for iqtree
combined.aa.aln.file <- "aln/outgroup/combined.aa.aln"

# Cluster has muscle 3.8.31, desktop has muscle 5
if(installr::is.windows()){
  system2("bin/muscle5.1.win64.exe", paste("-align", combined.aa.file,
                                            "-output", combined.aa.aln.file))

}else {
  system2("bin/muscle5.1.linux_intel64", paste("-align", combined.aa.file,
                                            "-output", combined.aa.aln.file))
}

combined.aa.aln <- ape::read.FASTA(combined.aa.aln.file)
# Read in a format for exon detection also
msa.aa.aln <- Biostrings::readAAMultipleAlignment(combined.aa.aln.file, format="fasta")
mouse.exons <- find.exons(msa.nt.aln, msa.aa.aln)

#### Plot extended outgroup AA tree ####

system2("iqtree", paste("-s ", "aln/outgroup/combined.aa.aln", 
                        "-bb 1000", # number of bootstrap replicates
                        "-alrt 1000", # number of replicates to perform SH-like approximate likelihood ratio test (SH-aLRT) 
                        "-nt AUTO"), # number of threads
        stdout = paste0("aln/outgroup/combined.aa.iqtree.log"), 
        stderr = paste0("aln/outgroup/combined.aa.iqtree.log"))


read.combined.outgroup.tree <- function(file){
  combined.outgroup.tree <- ape::read.tree(paste0(file))
  
  # Root the tree in the edge between Xenopus nodes and chicken
  xenopus.node <- ape::getMRCA(combined.outgroup.tree, c("Xenopus_ZFX.S","Xenopus_ZFX.L"))
  combined.outgroup.tree <- phytools::reroot(combined.outgroup.tree, xenopus.node, position = 0.1)
  ape::write.tree(combined.outgroup.tree, file = paste0("aln/outgroup/combined.aa.aln.rooted.treefile"))
  
  group_info <- split(combined.metadata$common.name, combined.metadata$group)
  combined.outgroup.tree <- groupOTU(combined.outgroup.tree, group_info, group_name = "group")
  combined.outgroup.tree
}

combined.outgroup.tree <- read.combined.outgroup.tree("aln/outgroup/combined.aa.aln.treefile")

outgroup.plot <- plot.tree(combined.outgroup.tree, col = "group")+
  coord_cartesian(clip="off", xlim = c(0, 0.9), ylim= c(-2, 62))

save.double.width("figure/combined.aa.tree.png", outgroup.plot)

# Also save the order of taxa in the outgroup tree to use later
outgroup.taxa.name.order <- get_taxa_name(outgroup.plot) 

#### Make mammal CDS NT tree #####

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

#### Plot mammal CDS NT tree #####

nt.aln.tree <- ape::read.tree(paste0(nt.aln.file, ".treefile"))

# Root the tree on platypus and resave
nt.aln.tree <- phytools::reroot(nt.aln.tree, which(nt.aln.tree$tip.label=="Platypus_ZFX"), position = 0.015)
ape::write.tree(nt.aln.tree, file = paste0(nt.aln.file, ".rooted.treefile"))

# Find the nodes that are ZFY vs ZFX and add to tree
group_info <- split(metadata$common.name, metadata$group)
nt.aln.tree <- groupOTU(nt.aln.tree, group_info, group_name = "group")

plot.zfx.zfy <- plot.tree(nt.aln.tree, col= "group") + coord_cartesian(clip="off", xlim = c(0, 0.5))

# Emphasise the ZFX / ZFY splits more by rotating the Laurasiatheria node
laurasiatheria.node <- ape::getMRCA(nt.aln.tree, c("Cat_ZFY", "Cat_ZFX"))
plot.zfx.zfy <- rotate(plot.zfx.zfy, laurasiatheria.node)

save.double.width("figure/zfx.zfy.tree.png", plot.zfx.zfy)

# Get the order of taxa names for reordering other plots later
taxa.name.order <- get_taxa_name(plot.zfx.zfy) 


# Drop the ZFY sequences and just look at the ZFX nodes in the tree
tree.zfx <- ape::drop.tip(nt.aln.tree, group_info$ZFY)
tree.zfx <- groupOTU(tree.zfx, group_info, group_name = "group")
ape::write.tree(tree.zfx, file = paste0(nt.aln.file, ".zfx.treefile"))
plot.zfx <- plot.tree(tree.zfx)
save.double.width("figure/zfx.tree.png", plot.zfx)

# Keep Zfy and Zfa, drop Zfx nodes. Keep outgroups
tree.zfy <- ape::keep.tip(nt.aln.tree, c(group_info$ZFY, group_info$Outgroup))
tree.zfy <- groupOTU(tree.zfy, group_info, group_name = "group")
ape::write.tree(tree.zfy, file = paste0(nt.aln.file, ".zfy.treefile"))
plot.zfy <- plot.tree(tree.zfy)
save.double.width("figure/zfy.tree.png", plot.zfy)

#### Make and plot mammal exon NT trees ####

# Aim here is to look for signs of gene conversion in the final exon as per other 
# papers.

create.exon.plot <- function(i){
  exon.aln <- as.matrix(ape.nt.aln)[,mouse.exons$start[i]:mouse.exons$end[i]]
  exon.aln <- ape::del.colgapsonly(exon.aln, threshold = 0.2) # remove columns with >20% gaps
  exon.aln.file <- paste0("aln/exons/exon_", mouse.exons$exon[i], ".aln")
  
  ape::write.FASTA(exon.aln, file = exon.aln.file)
  
  system2("iqtree", paste("-s ", exon.aln.file,
                          "-bb 1000", # number of bootstrap replicates
                          "-alrt 1000", # number of replicates to perform SH-like approximate likelihood ratio test (SH-aLRT)
                          "-nt AUTO"), # number of threads
          stdout = paste0(exon.aln.file, ".iqtree.log"),
          stderr = paste0(exon.aln.file, ".iqtree.log"))
  
  # Some exons will fail - too many gaps
  if(!file.exists(paste0(exon.aln.file, ".treefile"))) return()
  
  exon.tree <- ape::read.tree(paste0(exon.aln.file, ".treefile"))
  # Root the tree on platypus
  exon.tree <- phytools::reroot(exon.tree, which(exon.tree$tip.label=="Platypus_ZFX"), position = 0.015)
  
  # Find the nodes that are ZFY vs ZFX and add to tree
  group_info <- split(metadata$common.name, metadata$group)
  exon.tree <- groupOTU(exon.tree, group_info, group_name = "group")
  
  plot.exon.tree <- plot.tree(exon.tree, col="group")  + coord_cartesian(clip="off", xlim = c(0, 0.8))
  exon.fig.file <- paste0("figure/exon_", mouse.exons$exon[i], ".zfx.zfy.tree.png")
  save.double.width(exon.fig.file, plot.exon.tree)
  
  # Return for playing
  plot.exon.tree
}

exon.1.7.plots <- lapply(1:nrow(mouse.exons),create.exon.plot)

# We also want to look at all except exon 7

exon.1.6.aln <- as.matrix(ape.nt.aln)[,mouse.exons$start[1]:mouse.exons$end[6]]
exon.1.6.aln.file <- paste0("aln/exons/exon_1-6.aln")
ape::write.FASTA(exon.1.6.aln, file = exon.1.6.aln.file)

system2("iqtree", paste("-s ", exon.1.6.aln.file,
                        "-bb 1000", # number of bootstrap replicates
                        "-alrt 1000", # number of replicates to perform SH-like approximate likelihood ratio test (SH-aLRT)
                        "-nt AUTO"), # number of threads
        stdout = paste0(exon.1.6.aln.file, ".iqtree.log"),
        stderr = paste0(exon.1.6.aln.file, ".iqtree.log"))

exon.1.6.tree <- ape::read.tree(paste0(exon.1.6.aln.file, ".treefile"))
# Root the tree on platypus
exon.1.6.tree <- phytools::reroot(exon.1.6.tree, which(exon.1.6.tree$tip.label=="Platypus_ZFX"), position = 0.015)
# Find the nodes that are ZFY vs ZFX and add to tree
# Find the nodes that are ZFY vs ZFX and add to tree
group_info <- split(metadata$common.name, metadata$group)
exon.1.6.tree <- groupOTU(exon.1.6.tree, group_info, group_name = "group")

plot.exon.1.6.tree <- plot.tree(exon.1.6.tree, col="group")  + coord_cartesian(clip="off", xlim = c(0, 0.8))
exon.1.6.fig.file <- paste0("figure/exon_1-6.zfx.zfy.tree.png")
save.double.width(exon.1.6.fig.file, plot.exon.1.6.tree)

# Create a joint figure of exons 1-6 and exon 7

exon.joint.tree <- plot.exon.1.6.tree + exon.1.7.plots[[7]] + 
  patchwork::plot_annotation(tag_levels = list(c("Exons 1-6", "Exon 7"))) &
  theme(plot.tag = element_text(size = 6))
save.double.width("figure/exon.joint.tree.png", exon.joint.tree, height=120)
#### Plot exon by exon NT MSA #####

nt.aln.tidy <- tidy_msa(msa.nt.aln) %>%
  dplyr::mutate(name = factor(name, 
                                 levels = rev(taxa.name.order))) # sort reverse to match tree
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
          axis.title.x = element_blank(),
          panel.grid = element_blank())
  
  ggsave(paste0("figure/msa_exon_", exon, ".png"),  
         msa.plot, dpi = 300, units = "mm", width = 170)
}

#### Plot mammal CDS NT MSA #####

msa.nt.plot <- ggplot()+
  geom_msa(data = nt.aln.tidy, seq_name = T, font=NULL, 
           border=NA, color="Chemistry_NT", consensus_views = T, ref = "Platypus_ZFX", )+
  facet_msa(field = 400)+
  theme_minimal()+
  theme(axis.text = element_text(size=2),
        axis.title.y = element_blank(),
        panel.grid = element_blank())

ggsave(paste0("figure/msa.nt.mammal.png"),
       msa.nt.plot, dpi = 300, units = "mm", width = 170, height = 170)

#### Test selection globally in mammals ####

# ape::dnds(ape.nt.aln) # errors
seqin.aln <- seqinr::read.alignment(nt.aln.file, format = "fasta")
kaks.data <- seqinr::kaks(seqin.aln)

kaks.ratio <- kaks.data$ka / kaks.data$ks
kaks.pairwise <- dist2list(kaks.ratio, tri = F)
kaks.pairwise$col <- factor(kaks.pairwise$col, 
                            levels = taxa.name.order)
kaks.pairwise$row <- factor(kaks.pairwise$row, 
                            levels = taxa.name.order)

kaks.pairwise %<>% 
  dplyr::rowwise() %>%
  dplyr::mutate(rownum = which(row==taxa.name.order),
                colnum = which(col==taxa.name.order)) %>%
  dplyr::filter(rownum > colnum)

kaks.pairwise.plot <- ggplot(kaks.pairwise, aes(x = col, y = row))+
  geom_tile(aes(fill=value))+
  scale_fill_viridis_c(limits = c(0, 1), direction = -1)+
  labs(fill="dNdS")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 6),
        axis.title = element_blank(),
        legend.position = c(0.9, 0.15),
        legend.background = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8))
save.double.width("figure/dnds.png", kaks.pairwise.plot)

# Strong purifying selection in all pairs, but weaker in rodents

#### Test selection globally partition by partition in mammals ####

# Look at the final exon versus exons 1&3-6; is purifying selection more 
# pronounced in the ZFs versus exon 2 and exon 7?
exon1.3_6.locs <- c(mouse.exons$start_nt_codon_offset[1]:mouse.exons$end_nt_codon_offset[1], 
                    mouse.exons$start_nt_codon_offset[3]:mouse.exons$end_nt_codon_offset[6])
exon1.3_6.aln <- as.matrix(ape.nt.aln)[,exon1.3_6.locs]

ape::write.FASTA(exon1.3_6.aln, file = "aln/exons/exon_1_3-6.kaks.aln")
seqin.aln.exon.1.3_6 <- seqinr::read.alignment("aln/exons/exon_1_3-6.kaks.aln", format = "fasta")

exon2.aln <- as.matrix(ape.nt.aln)[,mouse.exons$start_nt_codon_offset[2]:(mouse.exons$end_nt_codon_offset[2])]
ape::write.FASTA(exon2.aln, file = "aln/exons/exon_2.kaks.aln")
seqin.aln.exon.2 <- seqinr::read.alignment("aln/exons/exon_2.kaks.aln", format = "fasta")
exon2.aln.msa <- ape::read.FASTA("aln/exons/exon_2.kaks.aln")

exon.7.aln <- as.matrix(ape.nt.aln)[,(mouse.exons$start_nt_codon_offset[7]):mouse.exons$end_nt_codon_offset[7]] 
ape::write.FASTA(exon.7.aln, file = "aln/exons/exon_7.kaks.aln")
seqin.aln.exon.7 <- seqinr::read.alignment("aln/exons/exon_7.kaks.aln", format = "fasta")

create.pairwise.kaks.data <- function(seqinr.aln){
  kaks.data <- seqinr::kaks(seqinr.aln, rmgap = FALSE)
  kaks.ratio <- kaks.data$ka / kaks.data$ks
  kaks.pairwise <- dist2list(kaks.ratio, tri = F)
  kaks.pairwise$col <- factor(kaks.pairwise$col, 
                                       levels = taxa.name.order)
  kaks.pairwise$row <- factor(kaks.pairwise$row, 
                                       levels = taxa.name.order)
  
  kaks.pairwise %>% 
    dplyr::rowwise() %>%
    dplyr::mutate(rownum = which(row==taxa.name.order),
                  colnum = which(col==taxa.name.order),
                  value  = ifelse(value==1, NA, value)) %>%  # values of exactly 1 are from missing data
    dplyr::filter(rownum > colnum)
}

kaks.exon.1.3_6 <- create.pairwise.kaks.data(seqin.aln.exon.1.3_6)
kaks.exon.2 <- create.pairwise.kaks.data(seqin.aln.exon.2)
kaks.exon.7 <- create.pairwise.kaks.data(seqin.aln.exon.7)

plot.kaks <- function(kaks.pairwise){
  ggplot(kaks.pairwise, aes(x = col, y = row))+
    geom_tile(aes(fill=value))+
    scale_fill_viridis_c(limits = c(0, 1), direction = -1, na.value="white")+
    labs(fill="dNdS")+
    theme_bw()+
    theme(axis.text.x = element_text(size = 4, angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(size = 4),
          axis.title = element_blank(),
          legend.position = c(0.8, 0.3),
          legend.background = element_blank(),
          legend.title = element_text(size = 5),
          legend.text = element_text(size = 5))
}

exon.1.3_6.kaks.pairwise.plot <- plot.kaks(kaks.exon.1.3_6)
exon.2.kaks.pairwise.plot <- plot.kaks(kaks.exon.2)
exon.7.kaks.pairwise.plot <- plot.kaks(kaks.exon.7)

exon.kaks.plot <- exon.1.3_6.kaks.pairwise.plot + exon.2.kaks.pairwise.plot + exon.7.kaks.pairwise.plot +
  patchwork::plot_annotation(tag_levels = c("A")) + plot_layout(axes="collect")

save.double.width("figure/exon.1_3-6.2.7.dnds.png", exon.kaks.plot, height=110)

#### Plot extended outgroup AA MSA ####

msa.outgroup.aln <- Biostrings::readAAMultipleAlignment("aln/outgroup/combined.aa.aln", format="fasta")
msa.outgroup.aln.tidy <- tidy_msa(msa.outgroup.aln)

msa.combined.aa.plot <- ggplot()+
  geom_msa(data = msa.outgroup.aln.tidy, seq_name = T, font=NULL, border=NA,
           consensus_views = T, ref = "Platypus_ZFX", alpha = 0.5
  )+
  labs(x = "Amino acid")+
  theme_minimal()+
  theme(axis.text = element_text(size=6),
        axis.title.y = element_blank(),
        panel.grid = element_blank())

#### Identify the locations of the ZFs in the extended outgroup AA MSA ####

combined.aa.aln <- Biostrings::readAAMultipleAlignment(combined.aa.aln.file, format="fasta")

locations.zf <- do.call(rbind, mapply(find.zf, aa=combined.aa.aln@unmasked, 
                                      sequence.name = names(combined.aa.aln@unmasked), 
                                      SIMPLIFY = FALSE)) %>%
  dplyr::mutate(sequence = factor(sequence, 
                                  levels = rev(outgroup.taxa.name.order))) %>% # sort reverse to match tree
  dplyr::rowwise() %>%
  dplyr::mutate(i = as.integer(sequence)) %>%  # Set the row indexes for plotting
  
  # Add the gapped nt alignment coordinates for nt sequences
  dplyr::mutate(start_nt_gapped = ifelse( sequence %in% names(msa.nt.aln@unmasked), # we have the nt alignment
                                          convert.to.gapped.coordinate(start_nt_ungapped,  msa.nt.aln@unmasked[[sequence]]),
                                          NA),
                end_nt_gapped = ifelse( sequence %in% names(msa.nt.aln@unmasked), # we have the nt alignment
                                          convert.to.gapped.coordinate(end_nt_ungapped,  msa.nt.aln@unmasked[[sequence]]),
                                          NA))

#### Identify the locations of the 9aaTAD in the extended outgroup AA MSA ####

# Only report the highest confidence 9aaTADs
locations.9aaTAD <- do.call(rbind, mapply(find.9aaTAD, aa=combined.aa.aln@unmasked, 
                                               sequence.name = names(combined.aa.aln@unmasked),
                                               rc.threshold=0, SIMPLIFY = FALSE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(sequence = factor(sequence, levels = rev(outgroup.taxa.name.order))) %>% # sort reverse to match tree
  dplyr::rowwise() %>%
  dplyr::mutate(i = as.integer(sequence)) %>%  # Set the row indexes for plotting
  
  # Add the gapped nt alignment coordinates for nt sequences
  dplyr::mutate(start_nt_gapped = ifelse( sequence %in% names(msa.nt.aln@unmasked), # we have the nt alignment
                                          convert.to.gapped.coordinate(start_nt_ungapped,  msa.nt.aln@unmasked[[sequence]]),
                                          NA),
                end_nt_gapped = ifelse( sequence %in% names(msa.nt.aln@unmasked), # we have the nt alignment
                                        convert.to.gapped.coordinate(end_nt_ungapped,  msa.nt.aln@unmasked[[sequence]]),
                                        NA))

#### Identify the locations of the NLS in the extended outgroup AA MSA ####
 
# Nuclear localisation sequence
 # using NLStradamus
 # Nguyen Ba AN, Pogoutse A, Provart N, Moses AM. NLStradamus: a simple Hidden Markov Model for nuclear localization signal prediction. BMC Bioinformatics. 2009 Jun 29;10(1):202. 
 
# Use aa translated sequence, no gaps
 
# aa.raw <- ape::trans(nt.raw)
write.FASTA(combined.aa.raw, file = "fasta/combined.aa.fas")

if(!installr::is.windows()){
  # ensure relatively lax threshold for broad detection
  # look for bipartite NLS
  # perl bin/nlstradamus.pl -i fasta/combined.aa.fas -t 0.5 -m 2 > nls/combined.aa.nls.out
  system2("perl", paste(" bin/nlstradamus.pl -i fasta/combined.aa.fas -t 0.5 > nls/combined.aa.nls.out"))
  
  # Remove the non-table output
  system2("cat", "nls/combined.aa.nls.out | grep -v 'Finished' | grep -v '=' | grep -v 'Analyzed' | grep -v 'sites' | grep -v 'Input' | grep -v 'Threshold' > nls/combined.aa.nls.filt.out")
}

# Read in the NLS prections
locations.NLS <- read_table("nls/combined.aa.nls.filt.out", 
                       col_names = c("sequence", "type", "posterior_prob", "start_ungapped", "end_ungapped", "aa"))

# Adjust the raw sequences to their positions in aa msa
locations.NLS$start_gapped <- sapply(1:nrow(locations.NLS), function(i) convert.to.gapped.coordinate(locations.NLS$start_ungapped[i], gapped.seq = combined.aa.aln@unmasked[[locations.NLS$sequence[i]]]))
locations.NLS$end_gapped <- sapply(1:nrow(locations.NLS), function(i) convert.to.gapped.coordinate(locations.NLS$end_ungapped[i], gapped.seq = combined.aa.aln@unmasked[[locations.NLS$sequence[i]]]))

locations.NLS %<>%
  dplyr::mutate(sequence = factor(sequence, levels = rev(outgroup.taxa.name.order)), # sort reverse to match tree
              start_nt_ungapped = start_ungapped * 3,
              end_nt_ungapped = end_ungapped * 3) %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(i = as.integer(sequence)) %>%  # Set the row indexes for plotting
  
  # Add the gapped nt alignment coordinates for nt sequences
  dplyr::mutate(start_nt_gapped = ifelse( sequence %in% names(msa.nt.aln@unmasked), # we have the nt alignment
                                          convert.to.gapped.coordinate(start_nt_ungapped,  msa.nt.aln@unmasked[[sequence]]),
                                          NA),
                end_nt_gapped = ifelse( sequence %in% names(msa.nt.aln@unmasked), # we have the nt alignment
                                        convert.to.gapped.coordinate(end_nt_ungapped,  msa.nt.aln@unmasked[[sequence]]),
                                        NA))

# Export the table of NLS sequences
write_tsv(locations.NLS, file = "figure/nls.output.tsv")
 
#### Export the locations of the  ZFs,  9aaTADs, and NLS ####

write_tsv(locations.zf %>% 
            dplyr::select(sequence, start_ungapped, end_ungapped, 
                          start_gapped, end_gapped, start_nt_ungapped, end_nt_ungapped, 
                          start_nt_gapped, end_nt_gapped),
          "figure/locations.zf.tsv")
write_tsv(locations.9aaTAD %>% 
            dplyr::select(sequence, hit, rc_score, start_ungapped, end_ungapped, 
                          start_gapped, end_gapped, start_nt_ungapped, end_nt_ungapped, 
                          start_nt_gapped, end_nt_gapped),
          "figure/locations.9aaTAD.tsv")
write_tsv(locations.NLS %>% 
            dplyr::select(sequence, aa, type, posterior_prob, start_ungapped, end_ungapped, 
                          start_gapped, end_gapped, start_nt_ungapped, end_nt_ungapped, 
                          start_nt_gapped, end_nt_gapped),
          "figure/locations.NLS.tsv")

#### Plot ZFs,  9aaTADs, and NLS in the extended outgroup AA MSA ####

# We need to combine the full set of structure locations into an overlapping set
# to be plotted in a single row. Keep those that overlap in >=5 species

find.common.overlaps <- function(locations.data){
  ranges.9aaTAD <- IRanges(start=locations.data$start_gapped, end = locations.data$end_gapped, names = locations.data$sequence)
  ranges.9aaTAD.reduce <- IRanges::reduce(ranges.9aaTAD)
  findOverlaps(ranges.9aaTAD, ranges.9aaTAD.reduce)
  n.sequences.in.range <- sapply(1:length(ranges.9aaTAD.reduce), function(i) length(subsetByOverlaps(ranges.9aaTAD, ranges.9aaTAD.reduce[i,])))
  as.data.frame(ranges.9aaTAD.reduce[n.sequences.in.range>4,])
}
ranges.ZF.common <- find.common.overlaps(locations.zf)
ranges.9aaTAD.common <- find.common.overlaps(locations.9aaTAD)
ranges.9aaTAD.common$label <- LETTERS[1:nrow(ranges.9aaTAD.common)]
ranges.NLS.common <- find.common.overlaps(locations.NLS)

#### Identify binding motifs of the ZFs in each species ####

# PWM prediction http://zf.princeton.edu/logoMain.php
# predict the ZF targets in each sequence via pwm_predict (http://zf.princeton.edu/download2.php)

if(!installr::is.windows()){
  old.wd <- getwd()
  setwd("./bin/pwm_predict")
  system2("./pwm_predict", "../../fasta/combined.aa.fas")
  setwd(old.wd)
  filesstrings::move_files(files = c("fasta/combined.aa.pwm"),
                           destinations = c("pwm"),
                           overwrite = TRUE)
  # Remove header and split the output PWMs to separate files
  system2("cat", "pwm/combined.aa.pwm | grep -v ';' | split -l 5 - pwm/zf_")
}

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
  dplyr::group_by(species) %>%
  dplyr::mutate(zf = paste0("zf_target_",row_number()),
                species = str_replace_all(species, ">", "")) %>%
  dplyr::rename(sequence = species) %>%
  tidyr::pivot_wider(id_cols = sequence, names_from = zf, values_from = consensus) %>%
  dplyr::group_by(zf_target_1, zf_target_2) %>%
  dplyr::summarise(total = n(), sequences = paste(sequence, collapse = ", ")) %>%
  dplyr::mutate(sequences = ifelse(total>40, "All others", sequences))


write_tsv(pwm.predictions, "figure/zf_targets.tsv")


#### Fetch divergence times to highlight the rapid evolution in the rodents ####

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
    Sys.sleep(1)
    
    
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

#### Run codeml to check for site-specific selection ####

# Adapted from Beginner's Guide on the Use of PAML to Detect Positive Selection
# https://academic.oup.com/mbe/article/40/4/msad041/7140562 for details

# Two files need to be uploaded to the cluster for running codeml:
# - alignment
# - tree file (unrooted)

# The treefile created earlier needs node names and branch lengths removing.
# Set node names to the empty string, and set branch lengths to 1
labelled.tree <- ape::read.tree(paste0(nt.aln.file, ".treefile"))
labelled.tree$node.label <- rep("", length(labelled.tree$node.label))
labelled.tree$edge.length <- rep(1, length(labelled.tree$edge.length))
ape::write.tree(labelled.tree, file = "paml/site-specific/zfxy.nt.aln.paml.treefile")

# We now read the tree as a raw string to remove the branch lengths entirely,
# separate the labels from taxon names and resave
labelled.tree <- readr::read_file("paml/site-specific/zfxy.nt.aln.paml.treefile")
labelled.tree <- gsub(":1", "", labelled.tree)
readr::write_file(labelled.tree, "paml/site-specific/zfxy.nt.aln.paml.treefile")

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

if(RUN_PAML & !installr::is.windows()){
  # Run codeml
  old.wd <- getwd()
  setwd("paml/site-specific")
  system2("codeml", "zfy.site-specific.paml.ctl",
          stdout = "zfy.site-specific.paml.log", 
          stderr = "zfy.site-specific.paml.log")
  
  # Extract lnl
  system2("cat", "zfy.out.txt | grep --before-context=5 'lnL' | grep -e 'lnL' -e 'Model'| paste -d ' '  - - > zfy.out.lnl.txt")
  
  setwd(old.wd)
  
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
}



#### Run codeml to look for selection specifically in Muroidea ####

# To look at the rodent clade, we need a rooted tree

# Find the MRCA of the rodents with rapid ZFY evolution
rodent.node <- ape::getMRCA(nt.aln.tree, c("Mouse_Zfy2", "Desert_hamster_Zfx-like_putative-Zfy"))

# Remove existing node labels, label the nodes and tips of the tree with #1
# for foreground branches
labelled.tree <- nt.aln.tree
labelled.tree$node.label <- rep("", length(labelled.tree$node.label))
labelled.tree$edge.length <- rep(1, length(labelled.tree$edge.length))
labelled.tree <- treeio::label_branch_paml(labelled.tree, rodent.node, "#1")
ape::write.tree(labelled.tree, file = "paml/branch-site/zfxy.nt.aln.paml.treefile")

# Then remove the branch lengths, separate the labels from taxon names and resave
newick.test <- read_file("paml/branch-site/zfxy.nt.aln.paml.treefile")
newick.test <- gsub(":1", "", newick.test) # branch lengths we set to 1
newick.test <- gsub("_#1", " #1", newick.test) # fg labels
write_file(newick.test, "paml/branch-site/zfxy.nt.aln.paml.fg.treefile")

# This control file tests site models With heterogeneous ω Across Sites
paml.site.branch.file <- paste0("seqfile   = ", nt.aln.file, " * alignment file\n",
                                "treefile  = paml/branch-site/zfxy.nt.aln.paml.fg.treefile * tree in Newick format without nodes\n", # 
                                "outfile   = paml/branch-site/site.branch.paml.out.txt\n",
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
write_file(paml.site.branch.file, "paml/branch-site/zfy.site-branch.paml.ctl")


paml.branch.site.output <- "paml/branch-site/zfy.site-branch.positive.sites.txt"
if(RUN_PAML & !installr::is.windows() & !file.exists(paml.branch.site.output)){
  
  old.wd <- getwd()
  setwd("paml/branch-site")
  # Run codeml
  system2("codeml", "zfy.branch-site.paml.ctl",
          stdout = "zfy.branch-site.paml.log", 
          stderr = "zfy.branch-site.paml.log")
  
  system2("cat", 'branch-site.paml.out.txt | grep "^ \\{2,5\\} [0-9]\\{1,3\\} [A-Z\\-]" > zfy.branch-site.positive.sites.txt' )
  setwd(old.wd)
  
  
  # Bayes Empirical Bayes (BEB) analysis (Yang, Wong & Nielsen 2005. Mol. Biol. Evol. 22:1107-1118)
  # Extract sites under positive selection
  # cat paml/site-branch/site.branch.paml.out.txt | grep "^ \{2,5\} [0-9]\{1,3\} [A-Z\-]" > zfy.site-branch.positive.sites.txt
  
  #system2("cat", 'paml/site-branch/site.branch.paml.out.txt | grep "^ \{2,5\} [0-9]\{1,3\} [A-Z\-]" > paml/site-branch/zfy.site-branch.positive.sites.txt' )
}

#### Create codeml files to look for selection in Muroidea across exon 1,3-6 ####

# We need to extract just the exon 1, 3-6 coordinates; however we must be careful
# not to disrupt codon boundaries. Coordinates have been adjusted in the mouse
# exon detection function
exon1.3_6.locs <- c(mouse.exons$start_nt_codon_offset[1]:mouse.exons$end_nt_codon_offset[1], 
                    mouse.exons$start_nt_codon_offset[3]:mouse.exons$end_nt_codon_offset[6])

exon1.3_6.aln <- as.matrix(ape.nt.aln)[,exon1.3_6.locs]
# exon.aln <- ape::del.colgapsonly(exon.aln, threshold = 0.2) # remove columns with >20% gaps
ape::write.FASTA(exon1.3_6.aln, file = "paml/exon_1_3-6/exon_1_3-6.aln")

# This control file tests branch site models With heterogeneous ω Across Sites for exons 1 and 3-6
paml.exon1.3_6.file <- paste0("seqfile     = exon_1_3-6.aln  * alignment file\n",
                                "treefile  = zfxy.nt.aln.paml.fg.treefile * tree in Newick format without nodes\n", # 
                                "outfile   = exon_1_3-6.paml.out.txt\n",
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
write_file(paml.exon1.3_6.file, "paml/exon_1_3-6/exon_1_3-6.paml.ctl")


exon1.3_6.output.file <- "paml/exon_1_3-6/exon_1_3-6.positive.sites.txt"
if(file.exists(exon1.3_6.output.file)){
  # Read the positive sites file 
  positive.sites <- read_table(exon1.3_6.output.file, col_names = c("site", "aa", "p"))
  positive.sites$p <- as.numeric(gsub("\\*+", "", positive.sites$p))
  
  positive.sites.y <- max(locations.zf$i) + 1.5
  
  positive.sites.plot <- ggplot()+
    # geom_rect(data = locations.zf,     aes(xmin=start_gapped, xmax=end_gapped, ymin=i-0.5, ymax=i+0.5), fill="grey", alpha=0.5)+
    # geom_rect(data = locations.9aaTAD, aes(xmin=start_gapped, xmax=end_gapped, ymin=i-0.5, ymax=i+0.5), fill="blue", alpha=0.5)+
    # geom_rect(data = locations.NLS,    aes(xmin=start_gapped, xmax=end_gapped, ymin=i-0.5, ymax=i+0.5), fill="green", alpha=0.5)+
    geom_rect(data = positive.sites,   aes(xmin=site-0.5, xmax=site+0.5, ymin=positive.sites.y, ymax=positive.sites.y+2, fill=p))+
    labs(x = "Site", fill = "p(ω>1)")+
    scale_fill_viridis_c()+
    theme_bw()
  save.double.width("figure/exon_1_3-6.positive.sites.png", positive.sites.plot, height = 85)
}


#### Create codeml files to look for selection in Muroidea across exon 2 ####

exon2.aln <- as.matrix(ape.nt.aln)[,mouse.exons$start_nt_codon_offset[2]:(mouse.exons$end_nt_codon_offset[2]-1)]
# exon.aln <- ape::del.colgapsonly(exon.aln, threshold = 0.2) # remove columns with >20% gaps
ape::write.FASTA(exon2.aln, file = "paml/exon_2/exon_2.aln")

# This control file tests branch site models With heterogeneous ω Across Sites for exons 1 and 3-6
paml.exon2.file <- paste0("seqfile   = exon_2.aln  * alignment file\n",
                              "treefile  = zfxy.nt.aln.paml.fg.treefile * tree in Newick format without nodes\n", # 
                              "outfile   = exon_2.paml.out.txt\n",
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
write_file(paml.exon2.file, "paml/exon_2/exon_2.paml.ctl")

exon.2.output.file <- "paml/exon_2/exon_2.positive.sites.txt"
if(file.exists(exon.2.output.file)){
  # Read the positive sites file 
  positive.sites <- read_table(exon.2.output.file, col_names = c("site", "aa", "p"))
  positive.sites$p <- as.numeric(gsub("\\*+", "", positive.sites$p))
  
  positive.sites.y <- max(locations.zf$i) + 1.5
  
  positive.sites.plot <- ggplot()+
    # geom_rect(data = locations.zf,     aes(xmin=start_gapped, xmax=end_gapped, ymin=i-0.5, ymax=i+0.5), fill="grey", alpha=0.5)+
    # geom_rect(data = locations.9aaTAD, aes(xmin=start_gapped, xmax=end_gapped, ymin=i-0.5, ymax=i+0.5), fill="blue", alpha=0.5)+
    # geom_rect(data = locations.NLS,    aes(xmin=start_gapped, xmax=end_gapped, ymin=i-0.5, ymax=i+0.5), fill="green", alpha=0.5)+
    geom_rect(data = positive.sites,   aes(xmin=site-0.5, xmax=site+0.5, ymin=positive.sites.y, ymax=positive.sites.y+2, fill=p))+
    labs(x = "Site", fill = "p(ω>1)")+
    scale_fill_viridis_c()+
    theme_bw()
  save.double.width("figure/exon_2.positive.sites.png", positive.sites.plot, height = 85)
}


#### Create codeml files to look for selection in Muroidea across exon 7 ####
exon7.aln <- as.matrix(ape.nt.aln)[,mouse.exons$start_nt_codon_offset[7]:mouse.exons$end_nt_codon_offset[7]]
# exon.aln <- ape::del.colgapsonly(exon.aln, threshold = 0.2) # remove columns with >20% gaps
ape::write.FASTA(exon7.aln, file = "paml/exon_7/exon_7.aln")

# This control file tests branch site models With heterogeneous ω Across Sites for exons 1 and 3-6
paml.exon7.file <- paste0("seqfile   = exon_7.aln  * alignment file\n",
                              "treefile  = zfxy.nt.aln.paml.fg.treefile * tree in Newick format without nodes\n", # 
                              "outfile   = exon_7.paml.out.txt\n",
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
write_file(paml.exon7.file, "paml/exon_7/exon_7.paml.ctl")

exon.7.output.file <- "paml/exon_7/exon_7.positive.sites.txt"
if(file.exists(exon.2.output.file)){
  # Read the positive sites file 
  positive.sites <- read_table(exon.2.output.file, col_names = c("site", "aa", "p"))
  positive.sites$p <- as.numeric(gsub("\\*+", "", positive.sites$p))
  
  positive.sites.y <- max(locations.zf$i) + 1.5
  
  positive.sites.plot <- ggplot()+
    # geom_rect(data = locations.zf,     aes(xmin=start_gapped, xmax=end_gapped, ymin=i-0.5, ymax=i+0.5), fill="grey", alpha=0.5)+
    # geom_rect(data = locations.9aaTAD, aes(xmin=start_gapped, xmax=end_gapped, ymin=i-0.5, ymax=i+0.5), fill="blue", alpha=0.5)+
    # geom_rect(data = locations.NLS,    aes(xmin=start_gapped, xmax=end_gapped, ymin=i-0.5, ymax=i+0.5), fill="green", alpha=0.5)+
    geom_rect(data = positive.sites,   aes(xmin=site-0.5, xmax=site+0.5, ymin=positive.sites.y, ymax=positive.sites.y+2, fill=p))+
    labs(x = "Site", fill = "p(ω>1)")+
    scale_fill_viridis_c()+
    theme_bw()
  save.double.width("figure/exon_7.positive.sites.png", positive.sites.plot, height = 85)
}


#### Run HyPhy RELAX to test for relaxed selection in Muroidea ####

# Newick tree tips and nodes need to be tagged with {Test} and {Reference}
# Set all branches and nodes to reference
hyphy.tree <- nt.aln.tree
hyphy.tree$node.label <- rep("", length(hyphy.tree$node.label)) # remove all node names
hyphy.tree$node.label <- paste0(hyphy.tree$node.label, "{Reference}") # everything is reference
hyphy.tree$tip.label <- paste0(hyphy.tree$tip.label, "{Reference}") # everything is reference
hyphy.tree <- treeio::label_branch_paml(hyphy.tree, rodent.node, "{Test}") # now rodents are test
hyphy.tree$node.label <- gsub("\\{Reference\\} \\{Test\\}", "{Test}", hyphy.tree$node.label) # remove reference from anything in test
hyphy.tree$tip.label <- gsub("\\{Reference\\} \\{Test\\}", "{Test}", hyphy.tree$tip.label) # remove reference from anything in test
ape::write.tree(hyphy.tree, file = "hyphy/mammal.nt.aln.hyphy.treefile")

# HyPhy should be installed in a conda environment
# This script should be invoked within the conda environment
if(!installr::is.windows()){
  
  # hyphy relax --alignment ../aln/zfxy.nt.aln --tree mammal.nt.aln.hyphy.treefile --reference "Reference" --test "Test" --output mammal.relax.json
  
  # Test relaxation of selection across the entire sequence in Muroidea
  system2("hyphy", " relax --alignment aln/zfxy.nt.aln --tree hyphy/mammal.nt.aln.hyphy.treefile --reference 'Reference' --test 'Test' --output hyphy/mammal.relax.json")
  
  # also run separately on the 3 exon partitions
  system2("hyphy", " relax --alignment paml/exon_1_3-6/exon_1_3-6.aln --tree hyphy/mammal.nt.aln.hyphy.treefile --reference 'Reference' --test 'Test' --output hyphy/exon_1_3-6.relax.json")
  system2("hyphy", " relax --alignment paml/exon_2/exon_2.aln --tree hyphy/mammal.nt.aln.hyphy.treefile --reference 'Reference' --test 'Test' --output hyphy/exon_2.relax.json")
  system2("hyphy", " relax --alignment paml/exon_7/exon_7.aln --tree hyphy/mammal.nt.aln.hyphy.treefile --reference 'Reference' --test 'Test' --output hyphy/exon_7.relax.json")
  
  # create an sh file to submit these manually
  hyphy.sh.data <- paste0("#!/bin/bash\n",
                          "source activate hyphy\n",
                          "hyphy relax --alignment paml/exon_1_3-6/exon_1_3-6.aln --tree hyphy/mammal.nt.aln.hyphy.treefile --reference 'Reference' --test 'Test' --output hyphy/exon_1_3-6.relax.json\n",
                          "hyphy relax --alignment paml/exon_2/exon_2.aln --tree hyphy/mammal.nt.aln.hyphy.treefile --reference 'Reference' --test 'Test' --output hyphy/exon_2.relax.json\n",
                          "hyphy relax --alignment paml/exon_2/exon_7.aln --tree hyphy/mammal.nt.aln.hyphy.treefile --reference 'Reference' --test 'Test' --output hyphy/exon_7.relax.json\n",
                          "hyphy relax --alignment aln/zfxy.nt.aln --tree hyphy/mammal.nt.aln.hyphy.treefile --reference 'Reference' --test 'Test' --output hyphy/mammal.relax.json\n"
  )
  write_file(hyphy.sh.data, "hyphy.sh")
  
  
  create.relax.k.tree <- function(json.file){
    # Read the json file and parse results
    hyphy.data <- jsonlite::read_json(json.file)
    # Note that the tree needs to terminate with ; otherwise read.tree returns NULL
    hyphy.input.tree <- ape::read.tree(text = paste0(hyphy.data$input$trees[["0"]], ";"))
    
    # Coloration of tree by k based on https://observablehq.com/@spond/plotting-relax-k-values-on-branches-of-the-tree
    
    # Make a dataframe with the k values and node numbers
    k.vals <- data.frame("k" = sapply(hyphy.data$`branch attributes`[["0"]], function(x) x$`k (general descriptive)`))
    k.vals$NodeLab <- rownames(k.vals)
    k.vals$node <- sapply(k.vals$NodeLab,  treeio::nodeid, tree = hyphy.input.tree)
    # Rescale values above 1 to the range 1-2 so we get a clean diverging scale
    k.vals$adj.k <- ifelse(k.vals$k <= 1, k.vals$k, (k.vals$k/max(k.vals$k)+1))
    
    # Add a new row for the root node
    k.vals[nrow(k.vals)+1,] <- list(0, "", length(hyphy.input.tree$tip.label)+1, 1)
    
    # Get the branch lengths from the HyPhy output
    k.vals$branch.length <-  sapply(k.vals$NodeLab, function(x) ifelse(x=="", 0, hyphy.data$`branch attributes`[["0"]][[x]]$`MG94xREV with separate rates for branch sets`))
    
    # Reorder the branches to match the node/tip order of the tree
    # Ordering in hyphy.input.tree$edge[,2] (the destination node, lengths are for incoming branches)
    branch.lengths <- unlist(sapply(hyphy.input.tree$edge[,2], function(x)  k.vals[k.vals$node==x,"branch.length"]))
    hyphy.input.tree$edge.length <- branch.lengths
    
    p <- ggtree(hyphy.input.tree) + geom_tiplab()
    p <- p %<+% k.vals + aes(colour=adj.k) + 
      scale_color_paletteer_c("ggthemes::Classic Red-Blue", 
                              direction = -1, 
                              limits = c(0, 2),
                              labels = c("0.0", "0.5", "1.0", round(max(k.vals$k)/2, digits = 0), round(max(k.vals$k), digits = 0)))+
      labs(color = "K")+
      coord_cartesian(xlim = c(0, 0.7))+
      annotate(geom="text", x = 0.3, y = 62, 
               label = paste("K(Muroidea) =", round(hyphy.data$`test results`$`relaxation or intensification parameter`, digits = 2)))+
      theme(legend.position = c(0.8, 0.2))
    
    p
  }
  
  
  if(file.exists("hyphy/mammal.relax.json")){
    relax.tree.all <- create.relax.k.tree("hyphy/mammal.relax.json")
    save.double.width("figure/mammal.relax.K.png", relax.tree.all)
  }
  if(file.exists("hyphy/exon_1_3-6.relax.json")){
    relax.tree.e1_3_6 <- create.relax.k.tree("hyphy/exon_1_3-6.relax.json")
    save.double.width("figure/exon_1_3-6.relax.K.png", relax.tree.e1_3_6)
  }
  if(file.exists("hyphy/exon_2.relax.json")){
    relax.tree.e2 <- create.relax.k.tree("hyphy/exon_2.relax.json")
    save.double.width("figure/exon_2.relax.K.png", relax.tree.e2)
  }
  if(file.exists("hyphy/exon_7.relax.json")){
    relax.tree.e7 <- create.relax.k.tree("hyphy/exon_7.relax.json")
    save.double.width("figure/exon_7.relax.K.png", relax.tree.e7)
  }
}

#### Run HyPhy GARD to test for recombination ####

if(!installr::is.windows()){
  
  system2("hyphy", " gard --alignment aln/zfxy.nt.aln --type codon --output hyphy/mammal.gard.json")
  
  # Read the json file and parse results
  hyphy.data <- jsonlite::read_json("hyphy/mammal.gard.json")
  
  siteBreakPointSupport <- data.frame("site" = as.numeric(names(hyphy.data$siteBreakPointSupport)),
                                         "support" = unlist(hyphy.data$siteBreakPointSupport))
  
  ggplot() +
    geom_rect(data = mouse.exons,          aes(xmin = start_aa, xmax = end_aa, ymin = 2, ymax = 3, fill = is_even))+
    scale_fill_manual(values=c("grey", "white"))+
    geom_rect(data = ranges.ZF.common,     aes(xmin=start, xmax=end, ymin=0, ymax=1), fill="grey", alpha=0.5)+
    geom_rect(data = ranges.9aaTAD.common, aes(xmin=start, xmax=end, ymin=0, ymax=1), fill="blue", alpha=0.5)+
    geom_rect(data = ranges.NLS.common,    aes(xmin=start, xmax=end, ymin=0, ymax=1), fill="green", alpha=0.5)+
    geom_vline(xintercept = 488)+
    geom_vline(xintercept = 884)+
    scale_x_continuous(expand = c(0,0))+
    guides(fill = "none")+
    # geom_point(data =siteBreakPointSupport, aes(x = site, y = log10(support) ))+
    theme_bw()+
    theme(axis.text.y = element_blank())
  
  tree.1.488 <- hyphy.data$breakpointData[[1]]$tree
  readr::write_file(tree.1.488, file = "hyphy/tree.1.488.treefile")
  tree.1.488 <- read.tree("hyphy/tree.1.488.treefile")
  
  tree.489.884 <- hyphy.data$breakpointData[[2]]$tree
  readr::write_file(tree.489.884, file = "hyphy/tree.489.884.treefile")
  tree.489.884 <- read.tree("hyphy/tree.489.884.treefile")
    
  # Trees are entirely unconvincing, leave this
  
}


#### Plot codeml branch-site model output  ####

if(file.exists(paml.branch.site.output)){
  # Read the positive sites file 
  positive.sites <- read_table(paml.branch.site.output, col_names = c("site", "aa", "p"))
  positive.sites$p <- as.numeric(gsub("\\*+", "", positive.sites$p))
  
  positive.sites.y <- max(locations.zf$i) + 1.5
  
  positive.sites.plot <- ggplot()+
    geom_rect(data = locations.zf,     aes(xmin=start_gapped, xmax=end_gapped, ymin=i-0.5, ymax=i+0.5), fill="grey", alpha=0.5)+
    geom_rect(data = locations.9aaTAD, aes(xmin=start_gapped, xmax=end_gapped, ymin=i-0.5, ymax=i+0.5), fill="blue", alpha=0.5)+
    geom_rect(data = locations.NLS,    aes(xmin=start_gapped, xmax=end_gapped, ymin=i-0.5, ymax=i+0.5), fill="green", alpha=0.5)+
    geom_rect(data = positive.sites,   aes(xmin=site-0.5, xmax=site+0.5, ymin=positive.sites.y, ymax=positive.sites.y+2, fill=p))+
    labs(x = "Site", fill = "p(ω>1)")+
    scale_fill_viridis_c()+
    theme_bw()
  save.double.width("figure/positive.sites.png", positive.sites.plot, height = 85)
}
#### Ancestral sequence reconstruction #####

# Read the ancestral reconstruction
ancestral.nt.seqs <- read.table(paste0(nt.aln.file, ".state") ,header=TRUE)

# We care about the eutherian common ancestor and the rodent ancestor
# Find these nodes

# Read the tree back to keep full node names
nt.aln.tree.nodes <- ape::read.tree(paste0(nt.aln.file, ".treefile"))
nt.aln.tree.nodes <- ape::root(nt.aln.tree.nodes, "Platypus_ZFX")

ggtree(nt.aln.tree.nodes) + 
  geom_tree() +
  geom_tiplab(size=2)+
  geom_nodelab(size=2, nudge_x = -0.003, nudge_y = 0.5, hjust=1,  node = "internal")+
  geom_treescale(fontsize =2, y = -1)+
  coord_cartesian(clip="off")+
  theme_tree() +
  theme(legend.position = "none")

# Node number returned includes number of tip labels; subtract to get node
rodent.node <- ape::getMRCA(nt.aln.tree.nodes, c("Mouse_Zfy1", "Desert_hamster_Zfx-like_putative-Zfy")) - length(nt.aln.tree.nodes$tip.label)
eutheria.node <- getMRCA(nt.aln.tree.nodes, c("Mouse_Zfy1", "African_bush_elephant_ZFY")) - length(nt.aln.tree.nodes$tip.label)

# Get ancestral rodent sequence
rodent.anc.nt <- ancestral.nt.seqs %>%
  dplyr::filter(Node == paste0("Node", rodent.node)) %>%
  dplyr::select(State)
rodent.anc.nt <- ape::as.DNAbin(rodent.anc.nt$State)

# Bind together as a matrix
rodent.plus.anc.nt.aln <- rbind(as.matrix(ape.nt.aln), rodent.anc.nt)

# Convert back to list and add names
rodent.ancestor.label <- "Ancestral_post-beaver_rodent"
rodent.plus.anc.nt.aln <- as.list.DNAbin(rodent.plus.anc.nt.aln)
names(rodent.plus.anc.nt.aln) <- c(names(ape.nt.aln), rodent.ancestor.label)

# Plot the MSA

# Find the tips under the ancestral node and get the labels
rodent.tip.labels <- tidytree::offspring( as_tibble(nt.aln.tree.nodes), .node = rodent.node+length(nt.aln.tree.nodes$tip.label), tiponly = T)


rodent.plus.anc.nt.aln.tidy <- tidy_msa(rodent.plus.anc.nt.aln) %>%
  dplyr::filter(name %in% c(rodent.ancestor.label, rodent.tip.labels$label))
rodent.plus.anc.nt.aln.tidy$name <- factor(rodent.plus.anc.nt.aln.tidy$name, 
                               levels = rev( c(taxa.name.order, rodent.ancestor.label))) # sort reverse to match tree

# Filter the zf locations and correct y locations
rodent.plus.anc.zf <- locations.zf %>%
  dplyr::filter(sequence %in% rodent.plus.anc.nt.aln.tidy$name) %>%
  dplyr::rowwise() %>%
  # Correct for different number of taxa in the y axis
  dplyr::mutate(i = which( levels(rodent.plus.anc.nt.aln.tidy$name) == sequence ) - (1+length(unique(taxa.name.order))- length(unique(rodent.plus.anc.nt.aln.tidy$name))))

# Filter 9aaTAD locations and correct the y locations
rodent.plus.anc.9aaTAD <- locations.9aaTAD %>%
  dplyr::filter(sequence %in% rodent.plus.anc.nt.aln.tidy$name) %>%
  dplyr::rowwise() %>%
  # Correct for different number of taxa in the y axis
  dplyr::mutate(i = which( levels(rodent.plus.anc.nt.aln.tidy$name) == sequence ) - (1+length(unique(taxa.name.order))- length(unique(rodent.plus.anc.nt.aln.tidy$name))))

rodent.plus.anc.NLS <- locations.NLS %>%
  dplyr::filter(sequence %in% rodent.plus.anc.nt.aln.tidy$name) %>%
  dplyr::rowwise() %>%
  # Correct for different number of taxa in the y axis
  dplyr::mutate(i = which( levels(rodent.plus.anc.nt.aln.tidy$name) == sequence ) - (1+length(unique(taxa.name.order))- length(unique(rodent.plus.anc.nt.aln.tidy$name))))



rodent.plus.anc.nt.msa.plot <- ggplot()+
  geom_msa(data = rodent.plus.anc.nt.aln.tidy, seq_name = T, font=NULL, 
           border=NA, color="Chemistry_NT", consensus_views = T, ref = rodent.ancestor.label, )+
  geom_rect(data = rodent.plus.anc.zf,     aes(xmin=start_nt_gapped, xmax=end_nt_gapped, ymin=i-0.5, ymax=i+0.5), fill="grey", alpha=0.5)+
  geom_rect(data = rodent.plus.anc.9aaTAD, aes(xmin=start_nt_gapped, xmax=end_nt_gapped, ymin=i-0.5, ymax=i+0.5), fill="blue", alpha=0.5)+
  geom_rect(data = rodent.plus.anc.NLS,    aes(xmin=start_nt_gapped, xmax=end_nt_gapped, ymin=i-0.5, ymax=i+0.5), fill="green", alpha=0.5)
if(file.exists(paml.branch.site.output)){
  rodent.plus.anc.nt.msa.plot <- rodent.plus.anc.nt.msa.plot +
    geom_point(data=positive.sites[positive.sites$p>0.9,], aes(x = site*3, y = p+6.5) , size=0.25)
}

rodent.plus.anc.nt.msa.plot <- rodent.plus.anc.nt.msa.plot +
  theme_minimal()+
  theme(axis.text = element_text(size=5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid = element_blank())

save.double.width("figure/rodent.ancestral.nt.msa.png", rodent.plus.anc.nt.msa.plot, height = 30)

#### Plot the conservation across the NT domains ####

# Plot conservation versus structural features

# Use Platypus ZFX as the outgroup?
platypus.nt.aln <-  ggmsa::tidy_msa(msa.nt.aln) %>%
  dplyr::filter(name=="Platypus_ZFX") %>%
  dplyr::select(-name, position, ref_char = character)

msa.nt.aln.tidy.conservation <- ggmsa::tidy_msa(msa.nt.aln) %>%
  merge( platypus.nt.aln, by = c("position")) %>%
  dplyr::filter(ref_char!="-") %>%
  dplyr::mutate(matchesRef = character==ref_char) %>%
  dplyr::group_by(position, matchesRef) %>%
  dplyr::summarise(n = n(), fraction = n/nrow(msa.nt.aln)) %>%
  dplyr::filter(matchesRef)

plot.conservation <- function(start, end){
  
  conservation.y <- max(locations.zf$i) + 1.5
  
  ggplot(msa.nt.aln.tidy.conservation)+
    # geom_rect(data=mouse.exons, aes( xmin = start-0.5, xmax = end+0.5, ymin=0, ymax=1, fill=exon), alpha=1)+
    geom_rect(data = locations.zf,     aes(xmin=start_nt_gapped, xmax=end_nt_gapped, ymin=i-0.5, ymax=i+0.5), fill="grey", alpha=0.5)+
    geom_rect(data = locations.9aaTAD, aes(xmin=start_nt_gapped, xmax=end_nt_gapped, ymin=i-0.5, ymax=i+0.5), fill="blue", alpha=0.5)+
    geom_rect(data = locations.NLS,    aes(xmin=start_nt_gapped, xmax=end_nt_gapped, ymin=i-0.5, ymax=i+0.5), fill="green", alpha=0.5)+
    
    geom_rect(data=msa.nt.aln.tidy.conservation,  aes(xmin=position-0.5, xmax=position+0.5,  ymin=conservation.y, ymax=conservation.y+2, fill=fraction))+
    scale_fill_viridis_c(direction = -1)+
    # geom_line( aes(x=position, y=fraction*nrow(msa.nt.aln)))+
    # scale_fill_manual(values = rep(c("grey", "white"), 4))+
    coord_cartesian(xlim = c(start, end))+
    labs(x = "Position in alignment", y = "Fraction of species with consensus nucleotide")+
    theme_bw()+
    theme(axis.text.y = element_blank())
}

conservation.plot.1 <- plot.conservation(1, ncol(msa.nt.aln)/3)
conservation.plot.2 <- plot.conservation(ncol(msa.nt.aln)/3+1, (ncol(msa.nt.aln)/3)*2)
conservation.plot.3 <- plot.conservation((ncol(msa.nt.aln)/3)*2+1, ncol(msa.nt.aln))

conservation.plot <- conservation.plot.1 / conservation.plot.2 / conservation.plot.3 + patchwork::plot_layout(axis_titles = "collect_y", guides = "collect")
save.double.width("figure/conservation_nt.png", conservation.plot, height = 150)

#### Calculate the conservation across the AA domains for outgroup levels ####

# Mammalian outgroup - opossum

# Calculate the fraction of sequences conserved with the given outgroup
calculate.conservation <-function(aa.aln, outgroup.name){
  
  # Find the characters in the reference sequence
  ref.aa.aln <- ggmsa::tidy_msa(aa.aln) %>%
    dplyr::filter(name==outgroup.name) %>%
    dplyr::select(-name, position, ref_char = character)
  
  # Filter the alignment to only mammalia after opossum
  # We can use the sequence level order for this
  outgroup.level <- which(outgroup.taxa.name.order=="Opossum_ZFX")
  species.to.calc <- outgroup.taxa.name.order[1:outgroup.level-1]
  
  ggmsa::tidy_msa(aa.aln) %>%
    merge( ref.aa.aln, by = c("position")) %>%
    dplyr::filter(ref_char!="-", 
                  name %in% species.to.calc) %>%
    dplyr::mutate(matchesRef = character==ref_char) %>%
    dplyr::group_by(position, matchesRef) %>%
    dplyr::summarise(n = n(), fraction = n/length(species.to.calc)) %>%
    dplyr::filter(matchesRef) %>%
    # Perform smoothing over aa moving windows
    dplyr::mutate(smoothed9 = slider::slide_dbl(fraction, mean, .before=4, .after = 4),
                  smoothed5 = slider::slide_dbl(fraction, mean, .before=2, .after = 2))
}

combined.aa.aln <- readAAMultipleAlignment(combined.aa.aln.file, format = "fasta")

# Use Xenopus ZFX.S as the comparison group
msa.aa.aln.tidy.frog.conservation    <- calculate.conservation(combined.aa.aln,"Xenopus_ZFX.S" )
msa.aa.aln.tidy.chicken.conservation <- calculate.conservation(combined.aa.aln,"Chicken_ZFX" )
msa.aa.aln.tidy.opossum.conservation <- calculate.conservation(combined.aa.aln,"Opossum_ZFX" )

# Combine the structural conservation plot with the aa tree
aa.structure.plot <- ggplot()+
  geom_tile(data = locations.zf,     aes(x=(start_gapped+end_gapped)/2,
                                         width = (end_gapped-start_gapped),
                                         y=sequence),
            fill="grey", alpha=0.5)+
  geom_tile(data = locations.9aaTAD,     aes(x=(start_gapped+end_gapped)/2,
                                         width = (end_gapped-start_gapped),
                                         y=sequence,
                                         fill=as.factor(round(rc_score))),
            alpha=0.9)+
  geom_tile(data = locations.NLS,     aes(x=(start_gapped+end_gapped)/2,
                                         width = (end_gapped-start_gapped),
                                         y=sequence),
            fill="green", alpha=0.5)+
  labs(fill = "RC score (%)")
aa.structure.plot <- annotate.structure.plot(aa.structure.plot, length(outgroup.taxa.name.order) + 1.5)

save.double.width("figure/aa.structure.png", aa.structure.plot, height = 120)


#### Plot the conservation of hydrophobicity across AA MSA ####

msa.aa.aln.tidy.hydrophobicity <- do.call(rbind, mapply(calc.hydrophobicity, aa=combined.aa.aln@unmasked, 
                                                sequence.name = names(combined.aa.aln@unmasked), 
                                                window.size = 9,
                                                SIMPLIFY = FALSE))  %>%
  dplyr::mutate(sequence = factor(sequence, levels = rev(outgroup.taxa.name.order))) # sort reverse to match tree

n.taxa <- length(outgroup.taxa.name.order) + 1.5

hydrophobicity.plot <- ggplot()+
  geom_tile(data=msa.aa.aln.tidy.hydrophobicity,  aes(x = position_gapped, y = sequence, fill=hydrophobicity_smoothed))+
  scale_fill_paletteer_c("ggthemes::Classic Red-Blue", direction = -1, limits = c(0, 1))+
  labs(fill="Hydrophobicity (9-window smooth)")
hydrophobicity.plot <- annotate.structure.plot(hydrophobicity.plot, n.taxa)
save.double.width("figure/hydrophobicity.convervation.tree.png", hydrophobicity.plot, height = 120)

#### Plot conservation of charge across AA MSA ####

msa.aa.aln.tidy.charge <- do.call(rbind, mapply(calc.charge, aa=combined.aa.aln@unmasked, 
                                                sequence.name = names(combined.aa.aln@unmasked), 
                                                window.size = 9,
                                                SIMPLIFY = FALSE))  %>%
  dplyr::mutate(sequence = factor(sequence, levels = rev(outgroup.taxa.name.order))) # sort reverse to match tree

charge.plot <- ggplot()+
  # Draw the charges per sequence
  geom_tile(data=msa.aa.aln.tidy.charge,  aes(x = position_gapped, y = sequence, fill=charge_smoothed))+
  scale_fill_paletteer_c("ggthemes::Classic Red-Blue", direction = 1, limits = c(-1, 1))+
  labs(fill="Charge (9-window smooth)")
charge.plot <- annotate.structure.plot(charge.plot, n.taxa)
save.double.width("figure/charge.convervation.tree.png", charge.plot, height = 120)

#### Combine all structure plots ####

structure.plot <- (aa.structure.plot) / (hydrophobicity.plot) / (charge.plot) + plot_layout(ncol = 1)
save.double.width("figure/structure.plot.png", structure.plot, height = 230)

#### Run GENECONV to test for gene conversion ####

# We want to look at gene conversion within species lineages. We also want to
# compare the ancestral ZFXs and ZFYs at each node To do this, create a separate
# tree for each of ZFX and ZFY, confirm that they have equivalent branches, then
# add the ancestral reconstruction to the geneconv config file.


# Write ZFX sequences to file
zfx.nt.aln <- msa.nt.aln@unmasked[names(msa.nt.aln@unmasked) %in% c(metadata$common.name[metadata$group=="ZFX"], 
                                                      metadata$common.name[metadata$group=="Outgroup"])]
Biostrings::writeXStringSet(zfx.nt.aln,  file = "aln/zfx_only/zfx.aln", format = "fasta")

# Write ZFY sequences to file
zfy.nt.aln <- msa.nt.aln@unmasked[names(msa.nt.aln@unmasked) %in% c(metadata$common.name[metadata$group=="ZFY"], 
                                                                    metadata$common.name[metadata$group=="Outgroup"])]
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

# Combine with existing ZFX/Y sequences
ape::write.dna(ape.nt.aln, "aln/ancestral.zfx.zfy.nodes.fa", format="fasta", 
               append = TRUE, colsep = "", nbcol=-1)

# Write geneconv configuration file specifying which sequences are in the same
# group

node.string.names <- zfy.zfx.common.nodes %>% 
  dplyr::mutate(GroupString= paste0("-group ZFX_", ZFX_node, "_ZFY_",ZFY_node, " ZFX_", ZFX_node, " ZFY_",ZFY_node ))
node.string <- paste(node.string.names$GroupString, collapse = "\n")

split.gene.names <- metadata %>%
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

#  Only on PATH in Windows, not on cluster yet
if(installr::is.windows()){
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
  
}




#### ####


# testing of side-by-side plots
# p <- ggtree(combined.outgroup.tree) + geom_tiplab(size=3, aes(col = group), align=F)
 # msaplot(p, combined.aa.aln.file, offset=0.4, width=2)+theme(legend.position = "none")

# p %>% aplot::insert_right(msa.combined.aa.plot)




# ggmsa::treeMSA_plot(ggtree(combined.outgroup.tree)+geom_tiplab(size=3, aes(col = group)), msa.outgroup.aln.tidy, color = "Chemistry_AA",
#                     border=NA, font=NULL)
#  