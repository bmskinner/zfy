# Common functions across scripts

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
                     "remotes", "patchwork", "ggpattern")
  
  sapply(cran.packages, install.cran)
  
  github.packages <- c('YuLab-SMU/ggtree', # since ggtree cannot install in Bioconductor 3.15 on cluster
                       "vragh/seqvisr", "vmikk/metagMisc")
  sapply(github.packages, install.github)
  
  bioconductor.packages <- c("msa", "ggmsa", "treeio")
  sapply(bioconductor.packages, install.bioconductor)
  return()
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
  # The second exon is alternatively spliced, so mark it with a striped fill
  # via ggpattern::geom_rect_pattern
  geom_rect_pattern(data = mouse.exons, aes(xmin = start_aa, xmax = end_aa, 
                                            ymin = y.start, ymax = y.end,
                                            pattern = exon==2, fill=exon), 
                    pattern_angle = 45, 
                    pattern_density = 0.5, # equal stripe widths
                    pattern_spacing = 0.05, # enough space for text label
                    pattern_fill = "grey", # match background color of exon
                     ...)
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
    labs(fill="Conservation (9 site average)")+
    add.conservation.track(msa.aa.aln.tidy.frog.conservation,    n.taxa,   n.taxa+2)+
    add.conservation.track(msa.aa.aln.tidy.chicken.conservation, n.taxa+3, n.taxa+5)+
    add.conservation.track(msa.aa.aln.tidy.opossum.conservation, n.taxa+6, n.taxa+8)+
    
    # Draw the structures
    add.track(ranges.ZF.common,     n.taxa+9, n.taxa+11, fill="lightgrey")+
    add.track(ranges.NLS.common,    n.taxa+9, n.taxa+11, fill="green", alpha = 0.5)+
    add.track(ranges.9aaTAD.common, n.taxa+9, n.taxa+11, fill="#19506F",  alpha = 0.9)+
    add.track.labels(ranges.9aaTAD.common, n.taxa+9, n.taxa+11)+   # Label the 9aaTADs
    
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
          panel.grid = element_blank())
  
  # Add the tree with the outgroups
  make.outgroup.mini.tree(combined.outgroup.tree,
                          text.labels = c("Xenopus", "Chicken", "Opossum", "Features", "Exons"))+
    plot +
    patchwork::plot_layout(widths = c(0.1, 0.9))
}

locate.zfs.in.alignment <- function(aa.alignment.file, nt.alignment.file, taxa.order){
  aa.aln <- Biostrings::readAAMultipleAlignment(aa.alignment.file, format="fasta")
  nt.aln <- Biostrings::readDNAMultipleAlignment(nt.alignment.file, format="fasta")
  
  # Find ZFs in each aa sequence, then find the gapped coordinates in the nt alignment
  do.call(rbind, mapply(find.zf, aa=aa.aln@unmasked, 
                                        sequence.name = names(aa.aln@unmasked), 
                                        SIMPLIFY = FALSE)) %>%
    dplyr::mutate(sequence = factor(sequence, 
                                    levels = rev(taxa.order))) %>% # sort reverse to match tree
    dplyr::rowwise() %>%
    dplyr::mutate(i = as.integer(sequence)) %>%  # Set the row indexes for plotting
    
    # Add the gapped nt alignment coordinates for nt sequences
    dplyr::mutate(start_nt_gapped = ifelse( sequence %in% names(nt.aln@unmasked), # we have the nt alignment
                                            convert.to.gapped.coordinate(start_nt_ungapped,  nt.aln@unmasked[[sequence]]),
                                            NA),
                  end_nt_gapped = ifelse( sequence %in% names(msa.nt.aln@unmasked), # we have the nt alignment
                                          convert.to.gapped.coordinate(end_nt_ungapped,  nt.aln@unmasked[[sequence]]),
                                          NA))
}

locate.9aaTADs.in.alignment <- function(aa.alignment.file, nt.alignment.file, taxa.order){
  aa.aln <- Biostrings::readAAMultipleAlignment(aa.alignment.file, format="fasta")
  nt.aln <- Biostrings::readDNAMultipleAlignment(nt.alignment.file, format="fasta")
  
  # Find ZFs in each aa sequence, then find the gapped coordinates in the nt alignment
  do.call(rbind, mapply(find.9aaTAD, aa=aa.aln@unmasked, 
                        sequence.name = names(aa.aln@unmasked),
                        rc.threshold=0, SIMPLIFY = FALSE)) %>%
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
                                          NA))
}

locate.NLS.in.alignment <- function(aa.alignment.file, nt.alignment.file, taxa.order){
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
                              col_names = c("sequence", "type", "posterior_prob", "start_ungapped", "end_ungapped", "aa"))
  
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