# Common functions used across scripts

#### Load all R packages installing if needed #####

load.packages <- function(){
  
  install.cran <- function(package){
    if(!require(package, character.only = TRUE, quietly = TRUE)){
      install.packages(package, lib = Sys.getenv("R_LIBS_USER"), repos = "https://cran.ma.imperial.ac.uk")
      library(package, character.only = TRUE, quietly = TRUE)
    }
  }
  
  install.bioconductor <- function(package){
    if(!require(package, character.only = TRUE, quietly = TRUE)){
      BiocManager::install(package, lib = Sys.getenv("R_LIBS_USER"), update = FALSE)
      library(package, character.only = TRUE, quietly = TRUE)
    }
  }
  
  install.github <- function(package){
    pkg.name <- gsub("^.*\\/", "", package)
    if(!require(pkg.name, character.only = TRUE, quietly = TRUE)){
      remotes::install_github(package, lib = Sys.getenv("R_LIBS_USER"))
      library(pkg.name, character.only = TRUE, quietly = TRUE)
    }
  }
  
  cran.packages <- c("tidyverse", "ape", "filesstrings", "seqinr", "phangorn",
                     "installr","treespace", "httr", "assertthat", "aplot",
                     "paletteer", "ggnewscale", "slider", "BiocManager",
                     "remotes", "patchwork", "ggpattern", "xlsx", "svglite")
  
  suppressPackageStartupMessages(sapply(cran.packages, install.cran))
  
  github.packages <- c('YuLab-SMU/ggtree', # since ggtree cannot install in Bioconductor 3.15 on cluster
                       "vmikk/metagMisc",  # for converting distance matrices to data frames
                       "fmichonneau/chopper" # functions for converting FASTA to NEXUS format
  )#"vragh/seqvisr"
  suppressPackageStartupMessages(sapply(github.packages, install.github))
  
  bioconductor.packages <- c("msa", "ggmsa", "treeio", "seqLogo")
  suppressPackageStartupMessages(sapply(bioconductor.packages, install.bioconductor))
  return()
}
load.packages() # Must be run on source so that FASTA files prepare automatically

#### Global variables #####
ZFY.TREE.COLOUR <- "#860086"  #3B4992
ZFX.TREE.COLOUR <- "#3CB22D" #EE0000
OUT.TREE.COLOUR <- "#303030"
TAD.COLOUR      <- "#00366C" # fill color max from "grDevices::Blues 3"
NLS.COLOUR      <- "#3CB22D"  #60CC52 #2CA02C
ZF.COLOUR       <- "grey"

MUSCLE.PATH   <- ifelse(installr::is.windows(), "bin/muscle5.1.win64.exe", "bin/muscle5.1.linux_intel64")
MACSE.PATH    <- "bin/macse_v2.07.jar"
CLUSTALO.PATH <- ifelse(installr::is.windows(), "bin/clustal-omega-1.2.2-win64/clustalo.exe", "clustalo") 

# Common file paths
FILES <- list(
  mammal.nt.fas = "fasta/mammal.nt.fas",
  combined.nt.fas = "fasta/combined.nt.fas",
  combined.nt.aln = "aln/combined/combined.nt.aln",
  combined.nt.nexus = "aln/combined/combined.nt.nex",
  combined.nt.partition = "aln/combined/combined.nt.partition",
  combined.aa.fas = "fasta/combined.aa.fas",
  combined.aa.aln = "aln/combined/combined.aa.aln",
  mammal.nt.aln = "aln/mammal/mammal.nt.aln",
  mammal.nt.nexus = "aln/mammal/mammal.nt.nex",
  mammal.nt.partition = "aln/mammal/mammal.nt.partition",
  mammal.aa.aln = "aln/mammal/mammal.aa.aln",
  mammal.nt.aln.treefile = "aln/mammal/mammal.nt.aln.treefile",
  combined.nt.aln.treefile = "aln/combined/combined.nt.aln.treefile",
  combined.aa.aln.treefile = "aln/combined/combined.aa.aln.treefile",
  paml.branch.site.output = "paml/branch-site/zfy.branch-site.positive.sites.txt",
  
  cox1.nt.fas = "fasta/cox1.nt.fas",
  cox1.nt.aln = "aln/cox1/cox1.nt.aln",
  cox1.aa.aln = "aln/cox1/cox1.aa.aln",
  
  final.intron.nt.fas     = "fasta/final.intron.nt.fas", # combined Zfx and Zfy
  final.intron.nt.aln     = "aln/final.intron/final.intron.nt.aln", # combined Zfx and Zfy
  final.intron.zfy.nt.fas = "fasta/final.intron.zfy.nt.fas",
  final.intron.zfx.nt.fas = "fasta/final.intron.zfx.nt.fas",
  final.intron.zfy.nt.aln = "aln/final.intron.zfy/final.intron.zfy.nt.aln",
  final.intron.zfx.nt.aln = "aln/final.intron.zfx/final.intron.zfx.nt.aln",
  
  final.intron.zfy.nt.aln.treefile = "aln/final.intron.zfy/final.intron.zfy.nt.aln.treefile",
  final.intron.zfx.nt.aln.treefile = "aln/final.intron.zfx/final.intron.zfx.nt.aln.treefile",
  
  final.intron.zfy.nt.aln.divvy.aln.treefile = "aln/final.intron.zfy/final.intron.zfy.nt.aln.divvy.aln.treefile",
  final.intron.zfx.nt.aln.divvy.aln.treefile = "aln/final.intron.zfx/final.intron.zfx.nt.aln.divvy.aln.treefile"
  
  
  
  
)

#### Invokations of external binaries #####

# Run muscle on the given input FASTA and output to the given alignment file base name
run.macse <- function(fa.file, aln.file, ...){
  if(!file.exists(MACSE.PATH)) stop("MACSE not present in bin directory")
  
  aa.out <- paste0(aln.file, ".aa.aln")
  nt.out <- paste0(aln.file, ".nt.aln")
  # Run a codon aware alignment with MACSE
  system2("java", paste("-jar ", MACSE.PATH, 
                        " -prog alignSequences",
                        "-seq",    fa.file, # input
                        "-out_NT", nt.out,  # output nt alignment
                        "-out_AA", aa.out,  # output aa alignment
                        ...), # any other options
          stdout = paste0(aln.file, ".macse.log"),  # logs
          stderr = paste0(aln.file, ".macse.log"))  # error logs
  
  # Return output file names
  list("aa" = aa.out, "nt" = nt.out)
}

# Create a alignment of two existing alignments (or an alignment and a single sequence FASTA)
# This preserves the structure of the existing alignments
run.clustalo.dual.profile <- function(aln.file.1, aln.file.2, out.file){
  if(!file.exists(CLUSTALO.PATH)) stop("ClustalO not found")
  system2(CLUSTALO.PATH, paste("--profile1", aln.file.1,
                               "--profile2", aln.file.2,
                               "-o", out.file,
                               "--force")) # overwrite existing alignments
}

# Run muscle on the given input FASTA and output to the given alignment file
run.muscle <- function(fa.file, aln.file){
  if(!file.exists(MUSCLE.PATH)) stop("Muscle not present in bin directory")
  cat(aln.file, "\n")
  system2(MUSCLE.PATH, paste("-align",  fa.file,
                             "-output", aln.file))
}

# Run IQTREE on the given alignment file (should end .aln)
run.iqtree <- function(aln.file, ...){
  system2("iqtree", paste("-s ", aln.file, 
                          "-nt 6",  # number of threads
                          ...), # any other arguments to iqtree
          stdout = gsub(".aln$", ".iqtree.log", aln.file), 
          stderr = gsub(".aln$", ".iqtree.log", aln.file))
  
  # Return file with tree
  paste0(aln.file, ".treefile")
}

# Run HyPhy MEME on the given alignment and tree. Returns the output file path
run.hyphy.meme <- function(nex.file, tree.file){
  
  json.file <- paste0( tree.file, ".json")
  bash.file <- paste0( tree.file, ".sh")
  # Create control script for MEME
  write_file(paste0("#!/bin/bash\n",
                    "source activate hyphy\n\n",
                    
                    "# MEME test for positive selection at individual sites\n",
                    "hyphy meme --alignment  ", nex.file, " --tree ", tree.file, " --branches Test --output ", json.file, "\n"
                    
  ),
  bash.file)
  system2("bash", bash.file)
  json.file
}

# Run divvier from within a conda environment (named divvier).Creates a divvied
# alignment.
# aln.file - the alignment to divvy
# Returns - the divvied alignment file path
run.divvier <- function(aln.file){
  bash.file <- paste0(aln.file, ".sh")
  write_file(paste0("#!/bin/bash\n",
                    "source activate divvier\n\n",
                    
                    "# Identify phylogenetically informatve sites with indels\n",
                    "divvier -divvygap ", aln.file,  "\n",
                    
                    "# Ensure file extension is aln for use in IQTREE\n",
                    "mv ", paste0(aln.file, ".divvy.fas"), " ", paste0(aln.file, ".divvy.aln"), "\n"
                    
  ),
  bash.file)
  system2("bash", bash.file)
  # Return the divvied file
  paste0(aln.file, ".divvy.aln")
}

# Run pwm_predict on the combined aa fasta file
run.pwm.predict <- function(){
  old.wd <- getwd()
  setwd("./bin/pwm_predict")
  system2("./pwm_predict", "-l 20 ../../fasta/combined.aa.fas") # ensure all ZFs linked
  setwd(old.wd)
  filesstrings::move_files(files = c("fasta/combined.aa.pwm"),
                           destinations = c("aln/pwm"),
                           overwrite = TRUE)
  # Remove header and split the output PWMs to separate files
  system2("cat", "aln/pwm/combined.aa.pwm | grep -v ';' | split -l 5 - aln/pwm/zf_")
}


#### Functions needed to load data files #####

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
  
  outgroup.sequences <-  c("Platypus_ZFX", "Opossum_ZFX", 
                           "Xenopus_ZFX.S", "Xenopus_ZFX.L", 
                           "Chicken_ZFX", "Zebra_finch_ZFX")
  
  # Grouping for trees
  metadata %<>%
    dplyr::mutate(group = case_when(grepl("ZFY", common.name, ignore.case=T) ~ "ZFY",
                                    common.name %in% outgroup.sequences ~ "Outgroup",
                                    T ~ "ZFX"))
  metadata
}

# Read all alignments and return as a list
read.alignments <- function(){
  list(
    aa.combined.ape = ape::read.FASTA(FILES$combined.aa.aln, type="AA"),
    # Read in Biostrings format for exon detection also
    aa.combined.biostrings = Biostrings::readAAMultipleAlignment(FILES$combined.aa.aln, format="fasta"),
    
    nt.combined.ape = ape::read.FASTA(FILES$combined.nt.aln, type="DNA"),
    # Read in Biostrings format for exon detection also
    nt.combined.biostrings = Biostrings::readAAMultipleAlignment(FILES$combined.nt.aln, format="fasta"),
    
    nt.mammal.ape = ape::read.FASTA(FILES$mammal.nt.aln),
    nt.mammal.biostrings = Biostrings::readDNAMultipleAlignment(FILES$mammal.nt.aln, format="fasta"),
    
    aa.mammal.ape <- ape::read.FASTA(FILES$mammal.aa.aln, type="AA"),
    aa.mammal.biostrings = Biostrings::readAAMultipleAlignment(FILES$mammal.aa.aln, format="fasta")
  )
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

#### Load input data #####

# Read source FASTA files, combine to single .fas and translate
# to protein. Read metadata
prepare.fas.files <- function(){
  
  # Putative Zfy sequences in rodents detected with NCBI gene search:
  # rodent[orgn:__txid9989] AND zinc finger X-chromosomal protein-like 
  
  # Read all unaligned sequence files with .fa extension
  fa.files <- list.files(path = "fasta/nt", pattern = "*.fa$", 
                         include.dirs = T, full.names = T)
  
  fa.read  <- lapply(fa.files, read.fasta)
  nt.raw   <- read.sequences(fa.read)
  
  metadata <- list()
  metadata$mammal <- read.metadata(fa.read) %>%
    dplyr::mutate(Species_common_name = str_replace(common.name, "_Z[F|f][X|x].*", ""),
                  Species_common_name = str_replace(Species_common_name, "(_putative)?(_|-)Z[F|f][X|x|Y|y].*", ""))
  
  # Write the combined fasta to file with .fas extension
  ape::write.FASTA(nt.raw, file = FILES$mammal.nt.fas)
  
  
  # Read outgroup NT FA files
  # Read all unaligned sequence files with .fa extension
  outgroup.nt.files <- list.files(path = "fasta/aa", pattern = "*.fa$", 
                                  include.dirs = T, full.names = T)
  
  outgroup.fa.read  <- lapply(outgroup.nt.files, read.fasta)
  outgroup.nt.raw   <- read.sequences(outgroup.fa.read)
  metadata$outgroup <-  read.metadata(outgroup.fa.read) %>%
    dplyr::mutate(Species_common_name = str_replace(common.name, "_Z[F|f][X|x].*", ""),
                  Species_common_name = str_replace(Species_common_name, "(_putative)?(_|-)Z[F|f][X|x|Y|y].*", ""))
  
  # Combine the outgroups with the mammals
  metadata$combined <- rbind(metadata$mammal, metadata$outgroup)
  
  # Write the unaligned combined fasta to file with .fas extension
  combined.nt.raw <- c(outgroup.nt.raw, nt.raw) # all sequences
  combined.aa.raw <- ape::trans(combined.nt.raw)
  
  ape::write.FASTA(combined.nt.raw, file = FILES$combined.nt.fas)
  ape::write.FASTA(combined.aa.raw, file = FILES$combined.aa.fas)
  
  # Create supplementary table with all accessions and sequence info
  cat("Writing metadata in Excel format\n")
  metadata$combined %>%
    dplyr::rename(Accession = accession, Species = species, Group = group,
                  Name_in_figures = common.name,Description = original.name  ) %>%
    dplyr::select(Accession,Species,Group,  Name_in_figures, Description ) %>%
    create.xlsx(., "figure/accessions.supplement.xlsx")
  
  metadata
}



#### Plotting functions ####

save.double.width <- function(filename, plot, height=170){ 
  ggsave(filename, plot, dpi = 600, units = "mm", width = 170, height = height)
  ggsave(str_replace(filename, ".png$", ".svg"), plot,dpi = 300, 
         units = "mm", width = 170, height = height)
}

plot.tree <- function(tree.data, tiplab.font.size = 2, ...){
  
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
    geom_tiplab(size=tiplab.font.size, aes_string(...))+
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
  
  result <- result + annotate(geom="segment", x=0.05, y = 6, linewidth = 0.5,
                              xend = 0.05, yend=18, col = ZFX.TREE.COLOUR)
  result <- result + annotate(geom="segment", x=0.05, y = 22, linewidth = 0.5,
                              xend = 0.05, yend=33, col = ZFX.TREE.COLOUR)
  result <- result + annotate(geom="segment", x=0.05, y = 34, linewidth = 0.5,
                              xend = 0.05, yend=47, col = ZFY.TREE.COLOUR)
  result <- result + annotate(geom="segment", x=0.05, y = 51, linewidth = 0.5,
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

# Given an alignment, calculate KaKs and return in tidy format
calc.kaks <- function(nt.aln.file){
  seqin.aln <- seqinr::read.alignment(nt.aln.file, format = "fasta")
  kaks.data <- seqinr::kaks(seqin.aln)
  
  kaks.ratio <- kaks.data$ka / kaks.data$ks
  
  # Convert to long format and remove pairwise diagonal
  metagMisc::dist2list(kaks.ratio, tri = F)
}

# Given an alignment and order of species, calculate KaKs and make a pairwise plot
# nt.aln.file - the nucleotide alignment
# species.order - a vector with the plotting order for species in the file
plot.kaks <- function(nt.aln.file, species.order, kaks.limits=c(0, 1)){
  
  # Convert to long format and remove pairwise diagonal
  kaks.pairwise <- calc.kaks(nt.aln.file)
  
  # Ensure spaces are converted to underscores in the species names
  # species.order <- gsub(" ", "_", species.order)
  
  kaks.pairwise %<>%
    dplyr::mutate(col = gsub("_", " ", col),
                  row = gsub("_", " ", row),
                  col = factor(col, levels = species.order, ordered = T),
                  row = factor(row, levels = species.order, ordered = T),
                  colnum = as.integer(col),
                  rownum = as.integer(row)) %>%
    dplyr::filter(rownum < colnum)
  
  max.y <- max(kaks.pairwise$value, na.rm = TRUE)
  
  y.limit = ifelse(is.infinite(max.y), 1, max.y)
  
  palette.choice <- scale_fill_viridis_c(limits = kaks.limits, direction = -1)

  ggplot(kaks.pairwise, aes(x = col, y = row))+
    geom_tile(aes(fill=value))+
    palette.choice+
    labs(fill="dNdS")+
    scale_x_discrete(limits=rev)+
    theme_bw()+
    theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 6),
          axis.title = element_blank(),
          legend.position = c(0.9, 0.8),
          legend.background = element_blank(),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8))
}


#### Finding structural features #####
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
                        MoreArgs = list(aa=ALIGNMENTS$aa.combined.biostrings),
                        SIMPLIFY = FALSE)) %>%
    dplyr::mutate(sequence = factor(sequence, 
                                    levels = rev(taxa.order))) %>% # sort reverse to match tree
    dplyr::rowwise() %>%
    dplyr::mutate(i = as.integer(sequence)) %>%  # Set the row indexes for plotting
    
    # Add the gapped nt alignment coordinates for nt sequences
    dplyr::mutate(start_nt_gapped = ifelse( sequence %in% names(ALIGNMENTS$nt.combined.biostrings@unmasked), # we have the nt alignment
                                            convert.to.gapped.coordinate(start_nt_ungapped,  ALIGNMENTS$nt.combined.biostrings@unmasked[[sequence]]),
                                            NA),
                  end_nt_gapped = ifelse( sequence %in% names(ALIGNMENTS$nt.combined.biostrings@unmasked), # we have the nt alignment
                                          convert.to.gapped.coordinate(end_nt_ungapped,  ALIGNMENTS$nt.combined.biostrings@unmasked[[sequence]]),
                                          NA)) %>%
    # Add the AA sequence covered by the ZF and motifs
    dplyr::mutate(aa_motif_gapped = as.character(ALIGNMENTS$aa.combined.biostrings@unmasked[[sequence]][start_gapped:end_gapped]),
                  aa_motif = gsub("-", "", aa_motif_gapped),
                  # 6    32  -1
                  # find the contact motif within the ZF (if available) .*H.{3}H(.)..(..).(.).{5}C..C.*
                  # via https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1287-y
                  # but reverse to match our sequences : C..C.{5}(.).(..)..(.)H.{3,4}.H
                  contact_bases = paste(str_match(aa_motif, "C..C.....(.).(..)..(.)H.{3,4}")[,2:4], collapse = "" ),
                  contact_bases = ifelse(contact_bases=="NANANA", NA, contact_bases) # clean up NAs
                  
                  # Find where the contact bases motif starts with respect to the start position of the ZF already identified
                  # contact_bases_start = str_locate(aa_motif, "C..C.....(.).(..)..(.)H.{3,4}")[,1],
                  # contact_base_ungapped_start = start_nt_ungapped-contact_bases_start+1,
                  # 
                  # # Get the nt motif in the ungapped alignment
                  # nt_seq = gsub("-", "", ALIGNMENTS$nt.combined.biostrings@unmasked[[sequence]]),
                  # nt_motif = substr(nt_seq, contact_base_ungapped_start, end_nt_ungapped),
                  # nt_motif_vec = ifelse(is.na(contact_bases), NA, 
                  #                       list(s2c(nt_motif))),
                  # nt_motif_tr = ifelse(is.na(contact_bases), NA, 
                  #                      paste(seqinr::translate(nt_motif_vec), collapse = "")),
                  # motif_ok = nt_motif_tr==aa_motif
                  # 
                  # # Find the coordinate of the contact bases in the ungapped nt alignment
                  # contact_base_1 = ifelse(is.na(contact_bases), 
                  #                         NA, contact_base_start+31), # 12aa from start
                  # contact_base_1_l = ifelse(is.na(contact_bases), 
                  #                            NA,  substr(nt_seq, contact_base_1, contact_base_1+2)),
                  # contact_base_23 = ifelse(is.na(contact_bases), 
                  #                         NA, contact_base_start+37), # 14+15aa from start
                  # contact_base_23_l = ifelse(is.na(contact_bases), 
                  #                           NA,  substr(nt_seq, contact_base_23, contact_base_23+5)),
                  # contact_base_4 = ifelse(is.na(contact_bases), 
                  #                         NA, contact_base_start+49), # 18aa from start
                  # contact_base_4_l = ifelse(is.na(contact_bases), 
                  #                           NA,  substr(nt_seq, contact_base_4, contact_base_4+2)),
                  # contact_base_nt = paste0(contact_base_1_l, contact_base_23_l, contact_base_4_l),
                  # contact_nt_tr = paste(seqinr::translate(s2c(contact_base_nt)), collapse = ""),
                  # contact_motif_ok = contact_nt_tr==contact_bases
                  
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
  
  if(!file.exists("aln/nls/combined.aa.nls.filt.out")){
    # Nuclear localisation sequence
    # using NLStradamus
    # Nguyen Ba AN, Pogoutse A, Provart N, Moses AM. NLStradamus: a simple Hidden Markov Model for nuclear localization signal prediction. BMC Bioinformatics. 2009 Jun 29;10(1):202. 
    
    # Use aa translated sequence, no gaps
    # ensure relatively lax threshold for broad detection
    # look for bipartite NLS
    # perl bin/nlstradamus.pl -i fasta/combined.aa.fas -t 0.5 -m 2 > nls/combined.aa.nls.out
    system2("perl", paste(" bin/nlstradamus.pl -i fasta/combined.aa.fas -t 0.5 > aln/nls/combined.aa.nls.out"))
    
    # Remove the non-table output
    system2("cat", "aln/nls/combined.aa.nls.out | grep -v 'Finished' | grep -v '=' | grep -v 'Analyzed' | grep -v 'sites' | grep -v 'Input' | grep -v 'Threshold' > aln/nls/combined.aa.nls.filt.out")
    
  }
  # Read in the NLS prections
  locations.NLS <- read_table("aln/nls/combined.aa.nls.filt.out", 
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

# Read the given treefile for mammals and outgroup species
read.combined.outgroup.tree <- function(file){
  combined.outgroup.tree <- ape::read.tree(paste0(file))
  
  rooted.file <- gsub("treefile", "rooted.treefile", file)
  
  # Root the tree in the edge between Xenopus nodes and chicken
  xenopus.node <- ape::getMRCA(combined.outgroup.tree, c("Xenopus_ZFX.S","Xenopus_ZFX.L"))
  combined.outgroup.tree <- phytools::reroot(combined.outgroup.tree, xenopus.node, position = 0.1)
  ape::write.tree(combined.outgroup.tree, file = rooted.file)
  
  mammal.gene.groups <- split(METADATA$combined$common.name, METADATA$combined$group)
  combined.outgroup.tree <- groupOTU(combined.outgroup.tree, mammal.gene.groups, group_name = "group")
  combined.outgroup.tree
}

#### Other functions #####

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
find.exons <- function(){
  
  # Mammal only NT alignment
  mouse.zfy1.nt <- as.character(ALIGNMENTS$nt.mammal.biostrings@unmasked$Mouse_Zfy1)
  mouse.zfy1.nt.ungapped <- str_remove_all(mouse.zfy1.nt, "-|\\*")
  
  # Combined NT alignment
  mouse.zfy1.nt.combined <- as.character(ALIGNMENTS$nt.combined.biostrings@unmasked$Mouse_Zfy1)
  mouse.zfy1.nt.combined.ungapped <- str_remove_all(mouse.zfy1.nt.combined, "-|\\*")
  
  # Combined AA alignment
  mouse.zfy1.aa <- as.character(ALIGNMENTS$aa.combined.biostrings@unmasked$Mouse_Zfy1)
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
  
  start.nt.combined <- sapply(mouse.exons$start_nt, str_locate, string=mouse.zfy1.nt.combined.ungapped)[1,]
  end.nt.combined   <- sapply(mouse.exons$end_nt,   str_locate, string=mouse.zfy1.nt.combined.ungapped)[2,]
  
  start.aa <- sapply(mouse.exons$start_aa, str_locate, string=mouse.zfy1.aa.ungapped)[1,]
  end.aa   <- sapply(mouse.exons$end_aa,   str_locate, string=mouse.zfy1.aa.ungapped)[2,]
  
  # Note that the aa and nt positions are not from equivalent alignments!
  # NT is from mammals, AA is from mammals + outgroups
  data <- data.frame("exon"     = mouse.exons$exon,
                     "start_nt"    = sapply(start.nt, convert.to.gapped.coordinate, mouse.zfy1.nt),
                     "end_nt"      = sapply(end.nt,   convert.to.gapped.coordinate, mouse.zfy1.nt),
                     "start_nt_combined"    = sapply(start.nt.combined, convert.to.gapped.coordinate, mouse.zfy1.nt),
                     "end_nt_combined"      = sapply(end.nt.combined,   convert.to.gapped.coordinate, mouse.zfy1.nt),
                     "start_aa" = sapply(start.aa, convert.to.gapped.coordinate, mouse.zfy1.aa),
                     "end_aa"   = sapply(end.aa,   convert.to.gapped.coordinate, mouse.zfy1.aa),
                     "is_even"  = sapply(mouse.exons$exon, function(i) as.numeric(i)%%2==0)) %>%
    dplyr::mutate(
      # Correct for gaps in the alignment affecting exon boundaries
      # Extend the next exon to start at the end of the current
      start_nt = ifelse(lag(end_nt)+1!=start_nt & !is.na(lag(end_nt)), lag(end_nt)+1, start_nt),
      
      length_nt = end_nt - start_nt + 1,
      length_aa = end_aa - start_aa + 1,
      # fix the offsets to get ORF of each exon
      start_nt_codon_offset = case_when(exon==1 ~ start_nt, # hardcode the codon offsets for subsetting
                                        exon>=2 ~ start_nt-1),
      end_nt_codon_offset   = case_when(exon<7 ~ end_nt-1,
                                        exon==7 ~ end_nt),
      corrected_offset_length = (end_nt_codon_offset - start_nt_codon_offset +1)%%3) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(exon_orf = subset.sequence(ALIGNMENTS$nt.combined.biostrings, "Mouse_Zfy1", start_nt_codon_offset, end_nt_codon_offset),
                  exon_orf_length = nchar(exon_orf),
                  exon_orf_triplet = exon_orf_length/3) # check divisible by 3
  
  
  data
  
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

# Given a Biostrings alignment, extract the sequence string at the given coordinates
# aln - msa in Biostrings format
# sequence.name - the name of the sequence from the alignment to extact
# start, end - coordinates in the gapped alignment
subset.sequence <- function(aln, sequence.name, start, end){
  as.character(aln@unmasked[[sequence.name]][start:end])
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

find.common.aa.overlaps <- function(locations.data){
  ranges.9aaTAD <- IRanges(start=locations.data$start_gapped, end = locations.data$end_gapped, names = locations.data$sequence)
  ranges.9aaTAD.reduce <- IRanges::reduce(ranges.9aaTAD)
  
  n.sequences.in.range <- sapply(1:length(ranges.9aaTAD.reduce), \(i) length(subsetByOverlaps(ranges.9aaTAD, ranges.9aaTAD.reduce[i,])))
  as.data.frame(ranges.9aaTAD.reduce[n.sequences.in.range>4,]) %>%
    dplyr::mutate(motif_number = row_number())
}

find.common.nt.overlaps <- function(locations.data){
  locations.data %<>%
    dplyr::filter(!is.na(start_nt_gapped) & !is.na(end_nt_gapped))
  
  ranges.9aaTAD <- IRanges(start=locations.data$start_nt_gapped, end = locations.data$end_nt_gapped, names = locations.data$sequence)
  ranges.9aaTAD.reduce <- IRanges::reduce(ranges.9aaTAD)
  
  n.sequences.in.range <- sapply(1:length(ranges.9aaTAD.reduce), \(i) length(subsetByOverlaps(ranges.9aaTAD, ranges.9aaTAD.reduce[i,])))
  as.data.frame(ranges.9aaTAD.reduce[n.sequences.in.range>4,]) %>%
    dplyr::mutate(motif_number = row_number()) %>%
    dplyr::rename(start_nt = start, end_nt = end, width_nt = width)
}

find.matching.range <-function(start.col, end.col, start, end, ranges){
  hit <- ranges[ranges[[start.col]]<=start & ranges[[end.col]]>=end,]$motif_number
  if(length(hit)==0) 0 else hit
}

# Reroot the given tree to the given nodes. If node labels contains more than one node,
# their MRCA will be used as the root.
reroot.tree <- function(tree, node.labels, position=0.01){
  if(!all(node.labels %in% tree$tip.label)){
    cat("Not all nodes (", node.labels, ") are in the tree, cannot reroot, not changing tree\n")
    return(tree)
  }
  if(length(node.labels)==1){
    
    return(phytools::reroot(tree, which(tree$tip.label==node.labels), position = position))
  }
  root.node <- ape::getMRCA(tree, node.labels)
  return(phytools::reroot(tree, root.node, position = position))
}


# Extract all nt and aa sequences from a combined-outgroup MSA for the given region, and calculate the
# consensus sequence. Ranges are inclusive.
extract.combined.alignment.region <- function(nt.start=NULL, nt.end=NULL, aa.start=NULL, aa.end=NULL){

  if(is.null(nt.start) & is.null(nt.end)){
    nt.start <- floor(max(1, aa.start*3)-2)
    nt.end <- ceiling(max(1, aa.end*3))
  }
  
  if(is.null(aa.start) & is.null(aa.end)){
    aa.start <- floor(max(1, nt.start/3))
    aa.end <- ceiling(max(1, nt.end/3))
  }
  
  nt.aln <- as.data.frame(as.matrix(ALIGNMENTS$nt.combined.biostrings)[,nt.start:nt.end]) %>%
    dplyr::mutate(Sequence = rownames(.),
                  Sequence = str_replace_all(Sequence, "_", " "), # remove underscores for pretty printing
                  Sequence = factor(Sequence, levels = combined.taxa.name.order)) %>%
    dplyr::arrange(as.integer(Sequence)) %>%
    dplyr::rename_with(.cols=starts_with("V"), .fn=\(x) paste0("Site_",as.integer(gsub("V", "", x))+nt.start-1) )
  
  aa.aln <- as.data.frame(as.matrix(ALIGNMENTS$aa.combined.biostrings)[,aa.start:aa.end]) %>%
    dplyr::mutate(Sequence = rownames(.),
                  Sequence = str_replace_all(Sequence, "_", " "), # remove underscores for pretty printing
                  Sequence = factor(Sequence, levels = combined.taxa.name.order)) %>%
    dplyr::arrange(as.integer(Sequence)) %>%
    dplyr::rename_with(.cols=starts_with("V"), .fn=\(x) paste0("Site_",as.integer(gsub("V", "", x))+aa.start-1) )
  
  
  # find the most common nucleotide in the consensus matrix
  nt.consensus <- paste(unlist(
    apply(consensusMatrix(ALIGNMENTS$nt.combined.biostrings)[,(nt.start):(nt.end)], 
          2, \(x) names( which(x==max(x)) )[1]  ) # Take only the first consensus in case of ties
  ), 
  collapse = "")
  
  aa.consensus <- paste(unlist(
    apply(consensusMatrix(ALIGNMENTS$aa.combined.biostrings)[,(aa.start):(aa.end)], 
          2, \(x) names( which(x==max(x)) )[1] ) # Take only the first consensus in case of ties
  ), 
  collapse = "")
  
  list(
    nt.aln   = nt.aln,
    aa.aln   = aa.aln,
    aa.start = aa.start,
    aa.end   = aa.end,
    aa.length = aa.end - aa.start + 1,
    nt.start = nt.start,
    nt.end   = nt.end,
    nt.length = nt.end - nt.start + 1,
    nt.consensus = nt.consensus,
    aa.consensus = aa.consensus
  )
}


# Extract all nt and aa sequences from a mammal-only MSA for the given region, and calculate the
# consensus sequence. Ranges are inclusive.
extract.mammal.alignment.region <- function(nt.start=NULL, nt.end=NULL, aa.start=NULL, aa.end=NULL){

  
  if(is.null(nt.start) & is.null(nt.end)){
    nt.start <- floor(max(1, aa.start*3)-2)
    nt.end <- ceiling(max(1, aa.end*3))
  }
  
  if(is.null(aa.start) & is.null(aa.end)){
    aa.start <- floor(max(1, nt.start/3))
    aa.end <- ceiling(max(1, nt.end/3))
  }
  
  nt.aln <- as.data.frame(as.matrix(ALIGNMENTS$nt.mammal.biostrings)[,nt.start:nt.end]) %>%
    dplyr::mutate(Sequence = rownames(.),
                  Sequence = str_replace_all(Sequence, "_", " "), # remove underscores for pretty printing
                  Sequence = factor(Sequence, levels = mammal.taxa.name.order)) %>%
    dplyr::arrange(as.integer(Sequence)) %>%
    dplyr::rename_with(.cols=starts_with("V"), .fn=\(x) paste0("Site_",as.integer(gsub("V", "", x))+nt.start-1) )
  
  aa.aln <- as.data.frame(as.matrix(ALIGNMENTS$aa.mammal.biostrings)[,aa.start:aa.end]) %>%
    dplyr::mutate(Sequence = rownames(.),
                  Sequence = str_replace_all(Sequence, "_", " "), # remove underscores for pretty printing
                  Sequence = factor(Sequence, levels = mammal.taxa.name.order)) %>%
    dplyr::arrange(as.integer(Sequence)) %>%
    dplyr::rename_with(.cols=starts_with("V"), .fn=\(x) paste0("Site_",as.integer(gsub("V", "", x))+aa.start-1) )
  
  
  # find the most common nucleotide in the consensus matrix
  nt.consensus <- paste(unlist(
    apply(consensusMatrix(ALIGNMENTS$nt.mammal.biostrings)[,(nt.start):(nt.end)], 
          2, \(x) names( which(x==max(x)) )[1]  ) # Take only the first consensus in case of ties
  ), 
  collapse = "")
  
  aa.consensus <- paste(unlist(
    apply(consensusMatrix(ALIGNMENTS$aa.mammal.biostrings)[,(aa.start):(aa.end)], 
          2, \(x) names( which(x==max(x)) )[1] ) # Take only the first consensus in case of ties
  ), 
  collapse = "")
  
  list(
    nt.aln   = nt.aln,
    aa.aln   = aa.aln,
    aa.start = aa.start,
    aa.end   = aa.end,
    aa.length = aa.end - aa.start + 1,
    nt.start = nt.start,
    nt.end   = nt.end,
    nt.length = nt.end - nt.start + 1,
    nt.consensus = nt.consensus,
    aa.consensus = aa.consensus
  )
}

#### ####