# Common functions across scripts

#### Global variables #####
ZFY.TREE.COLOUR <- "#860086"  #3B4992
ZFX.TREE.COLOUR <- "#3CB22D" #EE0000
OUT.TREE.COLOUR <- "#303030"
TAD.COLOUR      <- "#00366C" # fill color max from "grDevices::Blues 3"
NLS.COLOUR      <- "#3CB22D"  #60CC52 #2CA02C
ZF.COLOUR       <- "grey"
  
# Load all R packages installing if needed
load.packages <- function(){
  
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
                     "remotes", "patchwork", "ggpattern", "xlsx")
  
  sapply(cran.packages, install.cran)
  
  github.packages <- c('YuLab-SMU/ggtree', # since ggtree cannot install in Bioconductor 3.15 on cluster
                       "vmikk/metagMisc"  # for converting distance matrices to data frames
                        )#"vragh/seqvisr"
  sapply(github.packages, install.github)
  
  bioconductor.packages <- c("msa", "ggmsa", "treeio")
  sapply(bioconductor.packages, install.bioconductor)
  return()
}
load.packages() # load on source

# Set common file paths
set.file.paths <-function(){
  files <- list() # common file paths
  files$mammal.nt.fas <- "fasta/mammal.nt.fas"
  files$combined.nt.fas <- "fasta/combined.nt.fas"
  files$combined.aa.fas <- "fasta/combined.aa.fas"
  files$combined.aa.aln <- "aln/combined/combined.aa.aln"
  files$mammal.nt.aln <- "aln/mammal/mammal.nt.aln"
  files$mammal.aa.aln <- "aln/mammal/mammal.aa.aln"
  files$mammal.nt.aln.treefile <- paste0(files$mammal.nt.aln, ".treefile")
  files$combined.aa.aln.treefile <- paste0(files$combined.aa.aln, ".treefile")
  files$paml.branch.site.output <- "paml/branch-site/zfy.branch-site.positive.sites.txt"
  files
}

files <- set.file.paths()

prepare.fas.files <- function(){
  
  # Putative Zfy sequences in rodents detected with NCBI gene search:
  # rodent[orgn:__txid9989] AND zinc finger X-chromosomal protein-like 
  
  # Read all unaligned sequence files with .fa extension
  fa.files <- list.files(path = "fasta/nt", pattern = "*.fa$", 
                         include.dirs = T, full.names = T)
  
  fa.read  <- lapply(fa.files, read.fasta)
  nt.raw   <- read.sequences(fa.read)
  metadata.mammal <<- read.metadata(fa.read) %>%
    dplyr::mutate(Species_common_name = str_replace(common.name, "_Z[F|f][X|x].*", ""),
                  Species_common_name = str_replace(Species_common_name, "(_putative)?(_|-)Z[F|f][X|x|Y|y].*", ""))
  
  # Write the combined fasta to file with .fas extension
  ape::write.FASTA(nt.raw, file = files$mammal.nt.fas)
  
  
  # Read outgroup NT FA files
  # Read all unaligned sequence files with .fa extension
  outgroup.nt.files <- list.files(path = "fasta/aa", pattern = "*.fa$", 
                                  include.dirs = T, full.names = T)
  
  outgroup.fa.read  <- lapply(outgroup.nt.files, read.fasta)
  outgroup.nt.raw   <- read.sequences(outgroup.fa.read)
  metadata.outgroup <<-  read.metadata(outgroup.fa.read) %>%
    dplyr::mutate(Species_common_name = str_replace(common.name, "_Z[F|f][X|x].*", ""),
                  Species_common_name = str_replace(Species_common_name, "(_putative)?(_|-)Z[F|f][X|x|Y|y].*", ""))
  
  # Combine the outgroups with the mammals
  metadata.combined <<- rbind(metadata.mammal, metadata.outgroup)
  
  # Write the unaligned combined fasta to file with .fas extension
  combined.nt.raw <- c(outgroup.nt.raw, nt.raw) # all sequences
  combined.aa.raw <- ape::trans(combined.nt.raw)
  
  ape::write.FASTA(combined.nt.raw, file = files$combined.nt.fas)
  ape::write.FASTA(combined.aa.raw, file = files$combined.aa.fas)
  
  # Create supplementary table with all accessions and sequence info
  metadata.combined %>%
    dplyr::rename(Accession = accession, Species = species, Group = group,
                  Name_in_figures = common.name,Description = original.name  ) %>%
    dplyr::select(Accession,Species,Group,  Name_in_figures, Description ) %>%
    create.xlsx(., "figure/accessions.supplement.xlsx")
}

read.alignments <- function(){
  alignments <- list()
  alignments$aa.combined.ape <- ape::read.FASTA(files$combined.aa.aln, type="AA")
  # Read in Biostrings format for exon detection also
  alignments$aa.combined.biostrings <- Biostrings::readAAMultipleAlignment(files$combined.aa.aln, format="fasta")
  alignments$nt.mammal.ape <- ape::read.FASTA(files$mammal.nt.aln)
  alignments$nt.mammal.biostrings <- Biostrings::readDNAMultipleAlignment(files$mammal.nt.aln, format="fasta")
  alignments
}

# Read the given FASTA file and extract metadata
# Returns as a list containing FA sequence and metadata dataframe
# The FASTA file headers have common names manually added in [] brackets
# The file name should be the binomial species name
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

prepare.fas.files <- function(){
  
  # Putative Zfy sequences in rodents detected with NCBI gene search:
  # rodent[orgn:__txid9989] AND zinc finger X-chromosomal protein-like 
  
  # Read all unaligned sequence files with .fa extension
  fa.files <- list.files(path = "fasta/nt", pattern = "*.fa$", 
                         include.dirs = T, full.names = T)
  
  fa.read  <- lapply(fa.files, read.fasta)
  nt.raw   <- read.sequences(fa.read)
  metadata.mammal <<- read.metadata(fa.read)
  
  # Write the combined fasta to file with .fas extension
  ape::write.FASTA(nt.raw, file = files$mammal.nt.fas)
  
  
  # Read outgroup NT FA files
  # Read all unaligned sequence files with .fa extension
  outgroup.nt.files <- list.files(path = "fasta/aa", pattern = "*.fa$", 
                                  include.dirs = T, full.names = T)
  
  outgroup.fa.read  <- lapply(outgroup.nt.files, read.fasta)
  outgroup.nt.raw   <- read.sequences(outgroup.fa.read)
  metadata.outgroup <<-  read.metadata(outgroup.fa.read)
  
  # Combine the outgroups with the mammals
  metadata.combined <<- rbind(metadata.mammal, metadata.outgroup)
  
  # Write the unaligned combined fasta to file with .fas extension
  combined.nt.raw <- c(outgroup.nt.raw, nt.raw) # all sequences
  combined.aa.raw <- ape::trans(combined.nt.raw)
  
  ape::write.FASTA(combined.nt.raw, file = files$combined.nt.fas)
  ape::write.FASTA(combined.aa.raw, file = files$combined.aa.fas)
  
  # Create supplementary table with all accessions and sequence info
  metadata.combined %>%
    dplyr::rename(Accession = accession, Species = species, Group = group,
                  Name_in_figures = common.name,Description = original.name  ) %>%
    dplyr::select(Accession,Species,Group,  Name_in_figures, Description ) %>%
    create.xlsx(., "figure/accessions.supplement.xlsx")
}

read.sequences <- function(read.fasta.output){
  do.call(c, lapply(read.fasta.output, function(x) x$fa))
}

# Extract and combine metadata from a list of outputs of read.fasta
read.metadata <- function(read.fasta.output){
  metadata <- do.call(rbind, lapply(read.fasta.output, function(x) x$metadata))
  
  # Grouping for trees
  metadata %<>%
    dplyr::mutate(group = case_when(grepl("ZFY", common.name, ignore.case=T) ~ "ZFY",
                                    common.name %in% c("Platypus_ZFX", "Opossum_ZFX", 
                                                       "Xenopus_ZFX.S", "Xenopus_ZFX.L", "Chicken_ZFX") ~ "Outgroup",
                                    T ~ "ZFX"))
  metadata
}


save.double.width <- function(filename, plot, height=170){ 
  ggsave(filename, plot, dpi = 600, units = "mm", width = 170, height = height)
  ggsave(str_replace(filename, ".png$", ".svg"), plot,dpi = 300, 
         units = "mm", width = 170, height = height)
}

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
    if(i<=n) n <- n + 1
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
  data <- data.frame("exon"     = mouse.exons$exon,
             "start"    = sapply(start.nt, convert.to.gapped.coordinate, mouse.zfy1.nt),
             "end"      = sapply(end.nt,   convert.to.gapped.coordinate, mouse.zfy1.nt),
             "start_aa" = sapply(start.aa, convert.to.gapped.coordinate, mouse.zfy1.aa),
             "end_aa"   = sapply(end.aa,   convert.to.gapped.coordinate, mouse.zfy1.aa),
             "is_even"  = sapply(mouse.exons$exon, function(i) as.numeric(i)%%2==0)) %>%
    dplyr::mutate(
                  # Correct for gaps in the alignment affecting exon boundaries
                  # Extend the next exon to start at the end of the current
                  start = ifelse(lag(end)+1!=start & !is.na(lag(end)), lag(end)+1, start),
      
                  length_nt = end - start + 1,
                  length_aa = end_aa - start_aa + 1,
                  # fix the offsets to get ORF of each exon
                  start_nt_codon_offset = case_when(exon==1 ~ start, # hardcode the codon offsets for subsetting
                                                    exon>=2 ~ start-1),
                  end_nt_codon_offset   = case_when(exon<7 ~ end-1,
                                                    exon==7 ~ end),
                  corrected_offset_length = (end_nt_codon_offset - start_nt_codon_offset +1)%%3) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(exon_orf = subset.sequence(biostrings.nt.alignment, "Mouse_Zfy1", start_nt_codon_offset, end_nt_codon_offset),
                  exon_orf_length = nchar(exon_orf),
                  exon_orf_triplet = exon_orf_length/3) # check divisible by 3
  
  
  data
  
}

plot.tree <- function(tree.data, ...){
  
  # Remove underscores for pretty printing
  tree.data$tip.label <- str_replace_all(tree.data$tip.label, "_"," ")
  
  # Get the complete node labels
  # Separate out bootstrap info
  # Numbers in parentheses are SH-aLRT support (%) / ultrafast bootstrap support (%)
  node.label.values <- data.frame("label" = tree.data$node.label) %>%
    tidyr::separate_wider_delim(label, delim = "/", names = c("name","SHaLRT", "UFBoot"),
                                too_few = "align_end", too_many = "merge") %>%
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
    scale_color_manual(values = c(OUT.TREE.COLOUR, ZFX.TREE.COLOUR, ZFY.TREE.COLOUR))+
    # geom_nodelab(size=2, nudge_x = -0.003, nudge_y = 0.5, hjust=1,  node = "internal")+
    geom_nodepoint(size=1.5,  col="black")+
    geom_nodepoint(size=0.75,  col=node.label.values$colour)+
    geom_treescale(fontsize =2, y = -1, width = 0.05) +
    coord_cartesian(clip="off", ylim = c(-2, length(tree.data$tip.label)+1))+
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
    geom_tiplab(align=TRUE, linetype="dotted", linesize=.3) + # use tiplab to get lines
    scale_color_manual(values = c(OUT.TREE.COLOUR, ZFX.TREE.COLOUR, ZFY.TREE.COLOUR))+
    coord_cartesian(ylim = c(0.5, max.y+1), expand = FALSE)
  
  # Draw marker lines for ZFX and ZFX sequences
  result <- result + annotate(geom="text", x=0.05, y = 49, size = 2,
                              angle = 90, label="ZFY", col = ZFY.TREE.COLOUR)
  result <- result + annotate(geom="text", x=0.05, y = 20, size = 2,
                              angle = 90, label="ZFX", col = ZFX.TREE.COLOUR)
  
  result <- result + annotate(geom="segment", x=0.05, y = 6, size = 0.5,
                              xend = 0.05, yend=18, col = ZFX.TREE.COLOUR)
  result <- result + annotate(geom="segment", x=0.05, y = 22, size = 0.5,
                              xend = 0.05, yend=33, col = ZFX.TREE.COLOUR)
  result <- result + annotate(geom="segment", x=0.05, y = 34, size = 0.5,
                              xend = 0.05, yend=47, col = ZFY.TREE.COLOUR)
  result <- result + annotate(geom="segment", x=0.05, y = 51, size = 0.5,
                              xend = 0.05, yend=63, col = ZFY.TREE.COLOUR)
  
  # Whare is the top of the tree?
  curr.y <- length(combined.outgroup.tree.mini$tip.label) + 3
  for(i in 1:length(text.labels)){
    result <- result + annotate(geom="text", x=0.6, y=curr.y, 
                                label=text.labels[i], size=2, hjust=1, vjust=0.5)
    curr.y <- curr.y + 3
  }

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
add.track <- function(ranges, y.start, y.end, start_col = "start", end_col = "end", ...){
  geom_rect(data=ranges, 
            aes(xmin = .data[[start_col]], xmax = .data[[end_col]], ymin = y.start, ymax=y.end),
            ...)
}
# Add a rectangle track labels to a plot
add.track.labels <- function(ranges, y.start, y.end, start_col = "start", end_col = "end", label_col ="label", ...){
  geom_text(data=ranges,
            aes(x =(.data[[start_col]]+.data[[end_col]])/2, y=(y.start+y.end)/2, 
                label=.data[[label_col]]), size=2, ...)
}
# Add a conservation track to a plot
add.conservation.track <- function(ranges,  y.start, y.end, ...){
  geom_rect(data=ranges,  aes(xmin=position-0.55,
                              xmax=position+0.55,
                              ymin=y.start,
                              ymax=y.end,
                              fill=smoothed5))
}

# Add an exon track to a plot
add.exon.track <- function(y.start, y.end, start_col = "start_aa", end_col = "end_aa", ...){
  # The second exon is alternatively spliced, so mark it with a striped fill
  # via ggpattern::geom_rect_pattern
  geom_rect_pattern(data = mouse.exons, aes(xmin = .data[[start_col]], xmax = .data[[end_col]], 
                                            ymin = y.start, ymax = y.end,
                                            pattern = exon==2, fill=exon), 
                    pattern_xoffset = 0.04, # so the stripes don't obscure labels
                    pattern_yoffset = 0.0375, 
                    pattern_angle = 45, 
                    pattern_density = 0.5, # equal stripe widths
                    pattern_spacing = 0.05, # enough space for text label
                    pattern_fill = "grey", # match background color of exon
                     ...)
}
# Add exon labels track to a plot
add.exon.labels <- function(y.start, y.end, start_col = "start_aa", end_col = "end_aa", ...){
  geom_text(data=mouse.exons,
            aes(x =(.data[[start_col]]+.data[[end_col]])/2, y=(y.start+y.end)/2, label=exon), size=1.8, col="black")
}

# Annotate Zfy structures to a plot
# plot - the plot to annotate
# n.taxa - the number of rows within data; annotation tracks are above this
annotate.structure.plot <- function(plot, n.taxa){
  
  plot <- plot+
    # Draw the conservation with Xenopus, chicken and opossum
    
    
    # Draw the structures
    add.track(ranges.NLS.common,    n.taxa+8.5, n.taxa+11.5, fill=NLS.COLOUR, alpha = 1)+ # +8.5
    add.track(ranges.ZF.common,     n.taxa+9, n.taxa+11, fill=ZF.COLOUR)+ # 9 - 11
    add.track.labels(ranges.ZF.common, n.taxa+9, n.taxa+11, col="white", label_col = "motif_number")+   # Label the ZFs
    add.track(ranges.9aaTAD.common, n.taxa+9, n.taxa+11, fill=TAD.COLOUR,  alpha = 0.9)+  #9
    add.track.labels(ranges.9aaTAD.common, n.taxa+9, n.taxa+11, col="white")+   # Label the 9aaTADs
    
    new_scale_fill()+ 
    scale_fill_paletteer_c("grDevices::Cividis", direction = 1, limits = c(0, 1))+
    labs(fill="Conservation (5 site average)")+
    add.conservation.track(msa.aa.aln.tidy.frog.conservation,    n.taxa,   n.taxa+2)+
    add.conservation.track(msa.aa.aln.tidy.chicken.conservation, n.taxa+3, n.taxa+5)+
    add.conservation.track(msa.aa.aln.tidy.opossum.conservation, n.taxa+6, n.taxa+8)+
    
    new_scale_fill()+
    scale_fill_manual(values=c("white", "grey", "white", "grey", "white", "grey", "white"))+
    scale_pattern_color_manual(values=c("white", "white"))+

    scale_pattern_manual(values = c("none", "stripe")) + # which exons are patterned
    guides(fill = "none", pattern="none")+
    add.exon.track(n.taxa+12, n.taxa+14, col = "black")+ # color of border
    add.exon.labels(n.taxa+12, n.taxa+14)+
    
    scale_x_continuous(expand = c(0, 0), breaks = seq(0, 950, 50))+
    coord_cartesian(xlim = c(0, max(msa.aa.aln.tidy.frog.conservation$position)))+
    theme_bw()+
    theme(axis.text.y = element_blank(),
          axis.title = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size=6),
          legend.position = "top",
          legend.title = element_text(size = 6, vjust = 0.85),
          legend.text = element_text(size = 6),
          legend.key.height = unit(3, "mm"),
          legend.spacing.y = unit(2, "mm"),
          legend.box.spacing = unit(2, "mm"),
          panel.border = element_blank(),
          axis.line.x.bottom = element_line(),
          panel.grid = element_blank(),
          plot.margin = margin(l=0))
  
  # Add the tree with the outgroups
  make.outgroup.mini.tree(combined.outgroup.tree,
                          text.labels = c("Xenopus", "Chicken", "Opossum", "Features", "Exons"))+
    plot +
    patchwork::plot_layout(widths = c(0.1, 0.9))
}

locate.zfs.in.alignment <- function(taxa.order){
  
  taxa.order <- gsub(" ", "_", taxa.order) # in case _ had previously been replaced
  
  # Run the prediction from pwm_predict
  system2("hmmsearch", "--domtblout aln/pwm/combined.aa.hmm.dom.txt  bin/pwm_predict/zf_C2H2.ls.hmm fasta/combined.aa.fas")
  
  # Parse the resulting table
  if(!file.exists("aln/pwm/combined.aa.hmm.dom.txt")){
    stop("Missing hmmsearch output, cannot find ZFs")
  }
  zf.data <- read_table("aln/pwm/combined.aa.hmm.dom.txt", comment = "#",
             col_names = FALSE)
  zf.data <- zf.data[,c(1, 20, 21)]
  colnames(zf.data) <- c("sequence", "start_ungapped", "end_ungapped")

  # Find ZFs in each aa sequence, then find the gapped coordinates in the nt alignment
  do.call(rbind, mapply(find.zf, 
                        sequence.name = zf.data$sequence, 
                        start = zf.data$start_ungapped,
                        end = zf.data$end_ungapped,
                        MoreArgs = list(aa=alignments$aa.combined.biostrings),
                        SIMPLIFY = FALSE)) %>%
    dplyr::mutate(sequence = factor(sequence, 
                                    levels = rev(taxa.order))) %>% # sort reverse to match tree
    dplyr::rowwise() %>%
    dplyr::mutate(i = as.integer(sequence)) %>%  # Set the row indexes for plotting
    
    # Add the gapped nt alignment coordinates for nt sequences
    dplyr::mutate(start_nt_gapped = ifelse( sequence %in% names(alignments$nt.mammal.biostrings@unmasked), # we have the nt alignment
                                            convert.to.gapped.coordinate(start_nt_ungapped,  alignments$nt.mammal.biostrings@unmasked[[sequence]]),
                                            NA),
                  end_nt_gapped = ifelse( sequence %in% names(alignments$nt.mammal.biostrings@unmasked), # we have the nt alignment
                                          convert.to.gapped.coordinate(end_nt_ungapped,  alignments$nt.mammal.biostrings@unmasked[[sequence]]),
                                          NA)) %>%
    # Add the AA sequence covered by the ZF and motifs
    dplyr::mutate(aa_motif = gsub("-", "", alignments$aa.combined.biostrings@unmasked[[sequence]][start_gapped:end_gapped]),
                                                                             # 6    32  -1
                  # find the contact motif within the ZF (if available) .*H.{3}H(.)..(..).(.).{5}C..C.*
                  # via https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1287-y
                  # but reverse to match our sequences : C..C.{5}(.).(..)..(.)H.{3,4}.H
                  contact_bases = paste(str_match(aa_motif, "C..C.....(.).(..)..(.)H.{3,4}")[,2:4], collapse = "" ),
                  contact_bases = ifelse(contact_bases=="NANANA", NA, contact_bases) # clean up NAs
                  )
}

# Identify 9aaTADs and group overlapping TADS above a given rc threshold into 'superTADs'
locate.9aaTADs.in.alignment <- function(aa.alignment.file, nt.alignment.file, taxa.order){
  
  taxa.order <- gsub(" ", "_", taxa.order) # in case _ had previously been replaced
  
  aa.aln <- Biostrings::readAAMultipleAlignment(aa.alignment.file, format="fasta")
  nt.aln <- Biostrings::readDNAMultipleAlignment(nt.alignment.file, format="fasta")
  
  # Find ZFs in each aa sequence, then find the gapped coordinates in the nt alignment
 do.call(rbind, mapply(find.9aaTAD, aa=aa.aln@unmasked, 
                        sequence.name = names(aa.aln@unmasked),
                        rc.threshold=0, SIMPLIFY = FALSE)) %>%
    dplyr::rename(aa_motif = hit) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sequence = factor(sequence, levels = rev(taxa.order))) %>% # sort reverse to match tree
    dplyr::rowwise() %>%
    dplyr::mutate(i = as.integer(sequence)) %>%  # Set the row indexes for plotting
    
    # Add the gapped nt alignment coordinates for nt sequences
    dplyr::mutate(start_nt_gapped = ifelse( sequence %in% names(nt.aln@unmasked), # we have the nt alignment
                                            convert.to.gapped.coordinate(start_nt_ungapped,  nt.aln@unmasked[[sequence]]),
                                            NA),
                  end_nt_gapped = ifelse( sequence %in% names(nt.aln@unmasked), # we have the nt alignment
                                          convert.to.gapped.coordinate(end_nt_ungapped,  nt.aln@unmasked[[sequence]]),
                                          NA)) %>%
    dplyr::group_by(sequence)
}

locate.NLS.in.alignment <- function(aa.alignment.file, nt.alignment.file, taxa.order){
  
  taxa.order <- gsub(" ", "_", taxa.order) # in case _ had previously been replaced
  
  aa.aln <- Biostrings::readAAMultipleAlignment(aa.alignment.file, format="fasta")
  nt.aln <- Biostrings::readDNAMultipleAlignment(nt.alignment.file, format="fasta")
  
  if(!file.exists("nls/combined.aa.nls.filt.out")){
    # Nuclear localisation sequence
    # using NLStradamus
    # Nguyen Ba AN, Pogoutse A, Provart N, Moses AM. NLStradamus: a simple Hidden Markov Model for nuclear localization signal prediction. BMC Bioinformatics. 2009 Jun 29;10(1):202. 
    
    # Use aa translated sequence, no gaps
    # ensure relatively lax threshold for broad detection
    # look for bipartite NLS
    # perl bin/nlstradamus.pl -i fasta/combined.aa.fas -t 0.5 -m 2 > nls/combined.aa.nls.out
    system2("perl", paste(" bin/nlstradamus.pl -i fasta/combined.aa.fas -t 0.5 > nls/combined.aa.nls.out"))
    
    # Remove the non-table output
    system2("cat", "nls/combined.aa.nls.out | grep -v 'Finished' | grep -v '=' | grep -v 'Analyzed' | grep -v 'sites' | grep -v 'Input' | grep -v 'Threshold' > nls/combined.aa.nls.filt.out")
    
  }
  # Read in the NLS prections
  locations.NLS <- read_table("nls/combined.aa.nls.filt.out", 
                              col_names = c("sequence", "type", "posterior_prob", "start_ungapped", "end_ungapped", "aa_motif"))
  
  # Adjust the raw sequences to their positions in aa msa
  locations.NLS$start_gapped <- sapply(1:nrow(locations.NLS),
                                       function(i) convert.to.gapped.coordinate(locations.NLS$start_ungapped[i], 
                                                                                gapped.seq = aa.aln@unmasked[[locations.NLS$sequence[i]]]))
  locations.NLS$end_gapped <- sapply(1:nrow(locations.NLS), 
                                     function(i) convert.to.gapped.coordinate(locations.NLS$end_ungapped[i], 
                                                                              gapped.seq = aa.aln@unmasked[[locations.NLS$sequence[i]]]))
  
  locations.NLS %>%
    dplyr::mutate(sequence = factor(sequence, levels = rev(taxa.order)), # sort reverse to match tree
                  start_nt_ungapped = start_ungapped * 3,
                  end_nt_ungapped = end_ungapped * 3) %>% 
    dplyr::rowwise() %>%
    dplyr::mutate(i = as.integer(sequence)) %>%  # Set the row indexes for plotting
    
    # Add the gapped nt alignment coordinates for nt sequences
    dplyr::mutate(start_nt_gapped = ifelse( sequence %in% names(nt.aln@unmasked), # we have the nt alignment
                                            convert.to.gapped.coordinate(start_nt_ungapped,  nt.aln@unmasked[[sequence]]),
                                            NA),
                  end_nt_gapped = ifelse( sequence %in% names(nt.aln@unmasked), # we have the nt alignment
                                          convert.to.gapped.coordinate(end_nt_ungapped,  nt.aln@unmasked[[sequence]]),
                                          NA))
}

# Given a Biostrings alignment, extract the sequence string at the given coordinates
# aln - msa in Biostrings format
# sequence.name - the name of the sequence from the alignment to extact
# start, end - coordinates in the gapped alignment
subset.sequence <- function(aln, sequence.name, start, end){
  as.character(aln@unmasked[[sequence.name]][start:end])
}

# Create an Excel file with column filtering
create.xlsx = function(data, file.name, cols.to.fixed.size.font = NULL, cols.to.rich.text = NULL){
  
  # Set rich text formatting and highlight VV motifs in the given column index
  set.rich.text.on.vv <- function(wb, sh, col.index){
    
    normal.font.ref <-  xlsx::Font(wb, heightInPoints = 10, isBold=FALSE, name = "Courier New")$ref
    highlight.font.ref <- xlsx::Font(wb, heightInPoints = 10, color="red", isBold=TRUE,  name = "Courier New")$ref
    
    for(cell in xlsx::getCells(xlsx::getRows(sh), colIndex=col.index)){
      
      oldval <- xlsx::getCellValue(cell)
      vv.locs <- str_locate_all(oldval, "VV")
      
      if( nrow(vv.locs[[1]]) > 0 ){
        
        # Create a rich text string
        new.value <- rJava::.jnew("org/apache/poi/xssf/usermodel/XSSFRichTextString",
                                  oldval )
        
        # Set entire cell to normal style
        rJava::.jcall(obj=new.value, returnSig = "V",  # void return
                      method="applyFont", normal.font.ref)
        
        for(r in 1:nrow(vv.locs[[1]])){
          # Apply the new font to the correct indexes (0-indexed inclusive)
          rJava::.jcall(obj=new.value,returnSig = "V",  # void return
                        method="applyFont", 
                        as.integer(vv.locs[[1]][r,1]-1), # start index
                        as.integer(vv.locs[[1]][r,2]), # end index
                        highlight.font.ref)
        }
        
        
        # Set the new cell value and cast to a rich text string
        rJava::.jcall(cell, "V", "setCellValue",
                      rJava::.jcast(new.value, "org/apache/poi/ss/usermodel/RichTextString"))
      }
    }
  }
  
  oldOpt = options()
  options(xlsx.date.format="yyyy-mm-dd") # change date format
  wb = xlsx::createWorkbook(type = "xlsx")
  sh = xlsx::createSheet(wb)
  xlsx::addDataFrame(data, sh, row.names = F)
  xlsx::createFreezePane(sh, 2, 2, 2, 2) # freeze top row and first column
  cs <- xlsx::CellStyle(wb) + 
    xlsx::Font(wb,heightInPoints = 10, isBold = FALSE, name="Courier New")
  if(!is.null(cols.to.fixed.size.font)){
    for(i in xlsx::getCells(xlsx::getRows(sh), colIndex=cols.to.fixed.size.font)){
      xlsx::setCellStyle(i, cs)
    }
  }
  
  if(!is.null(cols.to.rich.text)) {
    for(col in cols.to.rich.text) set.rich.text.on.vv(wb, sh, col)
    
  }
  xlsx::autoSizeColumn(sh, 1:ncol(data))
  xlsx::saveWorkbook(wb, file=file.name)
  options(oldOpt)
}

read.time.tree.data <- function(){
  pairwise.times <- read_tsv("time.tree.data.tsv")
  # Swap the reciprocals and use adjusted age where available to keep tree ultrametric
  rbind(pairwise.times, pairwise.times %>% dplyr::rename(taxon_a_id = 2, 
                                                                           taxon_b_id = 1,
                                                                           scientific_name_a = 4,
                                                                           scientific_name_b = 3)) %>%
    dplyr::mutate(Mya = ifelse(adjusted_age > 0, adjusted_age, precomputed_age))
}


# Calculate the fraction of AA sequences conserved with the given outgroup
calculate.conservation <-function(aa.aln, outgroup.name){
  
  # Find the characters in the reference sequence
  ref.aa.aln <- ggmsa::tidy_msa(aa.aln) %>%
    dplyr::filter(name==outgroup.name) %>%
    dplyr::select(-name, position, ref_char = character)
  
  # Filter the alignment to only mammalia after opossum
  # We can use the sequence level order for this
  outgroup.level <- which(combined.taxa.name.order=="Opossum_ZFX")
  species.to.calc <- combined.taxa.name.order[1:outgroup.level-1]
  
  ggmsa::tidy_msa(aa.aln) %>%
    merge(ref.aa.aln, by = c("position")) %>%
    dplyr::filter(ref_char!="-", 
                  name %in% species.to.calc) %>%
    dplyr::mutate(matchesRef = character==ref_char) %>%
    dplyr::group_by(position, matchesRef) %>%
    dplyr::summarise(fraction = n()/length(species.to.calc)) %>%
    tidyr::pivot_wider(names_from = matchesRef, values_from = fraction, values_fill=0) %>%
    dplyr::rename(fraction = `TRUE`) %>%
    dplyr::select(-`FALSE`) %>%
    dplyr::ungroup() %>%
    # Perform smoothing over aa moving windows
    dplyr::mutate(smoothed9 = slider::slide_dbl(fraction, mean, .before=4, .after = 4),
                  smoothed5 = slider::slide_dbl(fraction, mean, .before=2, .after = 2))
}