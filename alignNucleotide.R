# Analysis pipeline for ZFX/ZFY evolutionary analysis

# The following binaries must be supplied in ./bin:
# macse_v2.07.jar
# muscle5.1.win64.exe / muscle5.1.linux_intelx64
# nlstradamul.pl
# pwm_predict/pwm_predict

# Other binaries are expected on the PATH:
# PAML
# IQ-TREE
# hmmsearch (HMMR v3.3)

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

source("functions.R")
load.packages()

source("find9aaTADs.R")
source("findZF.R")
source("calcCharge.R")
source("calcHydrophobicity.R")

cat("Packages loaded\n")

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

# Read all unaligned sequence files with .fa extension
fa.files <- list.files(path = "fasta/nt", pattern = "*.fa$", 
                       include.dirs = T, full.names = T)

fa.read  <- lapply(fa.files, read.fasta)
nt.raw   <- read.sequences(fa.read)
metadata.mammal <- read.metadata(fa.read)

# Write the combined fasta to file with .fas extension
mammal.nt.file <- "fasta/mammal.nt.fas"
ape::write.FASTA(nt.raw, file = mammal.nt.file)

#### Read outgroup NT FA files #####

# Read all unaligned sequence files with .fa extension
outgroup.nt.files <- list.files(path = "fasta/aa", pattern = "*.fa$", 
                                include.dirs = T, full.names = T)

outgroup.fa.read  <- lapply(outgroup.nt.files, read.fasta)
outgroup.nt.raw   <- read.sequences(outgroup.fa.read)
metadata.outgroup <-  read.metadata(outgroup.fa.read)

# Combine the outgroups with the mammals
metadata.combined <- rbind(metadata.mammal, metadata.outgroup)

# Write the unaligned combined fasta to file with .fas extension
combined.nt.raw <- c(outgroup.nt.raw, nt.raw) # all sequences
combined.nt.file <- "fasta/combined.nt.fas"
ape::write.FASTA(combined.nt.raw, file = combined.nt.file)

combined.aa.file <- "fasta/combined.aa.fas"
combined.aa.raw <- ape::trans(combined.nt.raw)
ape::write.FASTA(combined.aa.raw, file = combined.aa.file)

# Create supplementary table with all accessions and sequence info
supplementary.accessions.table <- metadata.combined %>%
  dplyr::rename(Accession = accession, Description = original.name,
                Name_in_figures = common.name, Species = species, Group = group)
write_tsv(supplementary.accessions.table, "figure/accessions.supplement.tsv")

#### Run combined mammal/outgroup AA alignment ####

# Use muscle to align the aa files and save out for iqtree
combined.aa.aln.file <- "aln/outgroup/combined.aa.aln"

# Cluster has muscle 3.8.31 on path., so specify the binary directly
if(!file.exists("bin/muscle5.1.linux_intel64")) stop("Muscle5.1 not present in bin directory")
system2("bin/muscle5.1.linux_intel64", paste("-align", combined.aa.file,
                                             "-output", combined.aa.aln.file))

combined.aa.aln <- ape::read.FASTA(combined.aa.aln.file)
# Read in Biostrings format for exon detection also
msa.aa.aln <- Biostrings::readAAMultipleAlignment(combined.aa.aln.file, format="fasta")

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
msa.nt.aln <- Biostrings::readDNAMultipleAlignment(nt.aln.file, format="fasta")

# Identify the coordinates of the exon boundaries in the gapped alignments
# Based on the mouse Zfy1 sequence
mouse.exons <- find.exons(msa.nt.aln, msa.aa.aln)

#### Create combined mammal/outgroup AA tree ####

system2("iqtree", paste("-s ", "aln/outgroup/combined.aa.aln", 
                        "-bb 1000", # number of bootstrap replicates
                        "-alrt 1000", # number of replicates to perform SH-like approximate likelihood ratio test (SH-aLRT) 
                        "-nt AUTO"), # number of threads
        stdout = paste0("aln/outgroup/combined.aa.iqtree.log"), 
        stderr = paste0("aln/outgroup/combined.aa.iqtree.log"))


#### Plot combined mammal/outgroup AA tree ####
read.combined.outgroup.tree <- function(file){
  combined.outgroup.tree <- ape::read.tree(paste0(file))
  
  # Root the tree in the edge between Xenopus nodes and chicken
  xenopus.node <- ape::getMRCA(combined.outgroup.tree, c("Xenopus_ZFX.S","Xenopus_ZFX.L"))
  combined.outgroup.tree <- phytools::reroot(combined.outgroup.tree, xenopus.node, position = 0.1)
  ape::write.tree(combined.outgroup.tree, file = paste0("aln/outgroup/combined.aa.aln.rooted.treefile"))
  
  mammal.gene.groups <- split(metadata.combined$common.name, metadata.combined$group)
  combined.outgroup.tree <- groupOTU(combined.outgroup.tree, mammal.gene.groups, group_name = "group")
  combined.outgroup.tree
}

combined.outgroup.tree <- read.combined.outgroup.tree("aln/outgroup/combined.aa.aln.treefile")

combined.aa.tree <- plot.tree(combined.outgroup.tree, col = "group")+
  coord_cartesian(clip="off", xlim = c(0, 0.9), ylim= c(-2, 62))

save.double.width("figure/combined.aa.tree.png", combined.aa.tree)

# Also save the order of taxa in the outgroup tree to use later
combined.taxa.name.order <- ggtree::get_taxa_name(combined.aa.tree) 

#### Create mammal CDS NT tree #####

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
# The root is arbitrarily placed in the platypus branch to fit neatly
nt.aln.tree <- phytools::reroot(nt.aln.tree, which(nt.aln.tree$tip.label=="Platypus_ZFX"), position = 0.015)
ape::write.tree(nt.aln.tree, file = paste0(nt.aln.file, ".rooted.treefile"))

# Find the nodes that are ZFY vs ZFX and add to tree
mammal.gene.groups <- split(metadata.mammal$common.name, metadata.mammal$group)
nt.aln.tree <- tidytree::groupOTU(nt.aln.tree, mammal.gene.groups, group_name = "group")

plot.zfx.zfy <- plot.tree(nt.aln.tree, col= "group") + coord_cartesian(clip="off", xlim = c(0, 0.5))

# Emphasise the ZFX / ZFY splits more by rotating the Laurasiatheria node
laurasiatheria.node <- ape::getMRCA(nt.aln.tree, c("Cat_ZFY", "Cat_ZFX"))
plot.zfx.zfy <- rotate(plot.zfx.zfy, laurasiatheria.node)

save.double.width("figure/zfx.zfy.tree.png", plot.zfx.zfy)

# Get the order of taxa names for reordering other plots later
mammal.taxa.name.order <- get_taxa_name(plot.zfx.zfy) 

# Now we want to view separate trees for Zfx and Zfy based on this combined tree

# Drop the ZFY sequences and just look at the ZFX nodes in the tree
tree.zfx <- ape::drop.tip(nt.aln.tree, mammal.gene.groups$ZFY)
tree.zfx <- groupOTU(tree.zfx, mammal.gene.groups, group_name = "group")
ape::write.tree(tree.zfx, file = paste0(nt.aln.file, ".zfx.treefile"))
plot.zfx <- plot.tree(tree.zfx)
save.double.width("figure/zfx.tree.png", plot.zfx)

# Keep Zfy and outgroups, drop other tips
tree.zfy <- ape::keep.tip(nt.aln.tree, c(mammal.gene.groups$ZFY, mammal.gene.groups$Outgroup))
tree.zfy <- groupOTU(tree.zfy, mammal.gene.groups, group_name = "group")
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
  mammal.gene.groups <- split(metadata.mammal$common.name, metadata.mammal$group)
  exon.tree <- groupOTU(exon.tree, mammal.gene.groups, group_name = "group")
  
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
mammal.gene.groups <- split(metadata.mammal$common.name, metadata.mammal$group)
exon.1.6.tree <- groupOTU(exon.1.6.tree, mammal.gene.groups, group_name = "group")

plot.exon.1.6.tree <- plot.tree(exon.1.6.tree, col="group")  + coord_cartesian(clip="off", xlim = c(0, 0.8))
exon.1.6.fig.file <- paste0("figure/exon_1-6.zfx.zfy.tree.png")
save.double.width(exon.1.6.fig.file, plot.exon.1.6.tree)

# Create a joint figure of exons 1-6 and exon 7

exon.joint.tree <- plot.exon.1.6.tree + exon.1.7.plots[[2]] + exon.1.7.plots[[7]] + 
  patchwork::plot_annotation(tag_levels = list(c("Exons 1-6", "Exon 2", "Exon 7"))) &
  theme(plot.tag = element_text(size = 6))
save.double.width("figure/exon.joint.tree.png", exon.joint.tree, height=120)
#### Test selection globally in mammals ####

# ape::dnds(ape.nt.aln) # errors
seqin.aln <- seqinr::read.alignment(nt.aln.file, format = "fasta")
kaks.data <- seqinr::kaks(seqin.aln)

kaks.ratio <- kaks.data$ka / kaks.data$ks
kaks.pairwise <- dist2list(kaks.ratio, tri = F)
kaks.pairwise$col <- factor(kaks.pairwise$col, 
                            levels = mammal.taxa.name.order)
kaks.pairwise$row <- factor(kaks.pairwise$row, 
                            levels = mammal.taxa.name.order)

kaks.pairwise %<>% 
  dplyr::rowwise() %>%
  dplyr::mutate(rownum = which(row==mammal.taxa.name.order),
                colnum = which(col==mammal.taxa.name.order)) %>%
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
                                       levels = mammal.taxa.name.order)
  kaks.pairwise$row <- factor(kaks.pairwise$row, 
                                       levels = mammal.taxa.name.order)
  
  kaks.pairwise %>% 
    dplyr::rowwise() %>%
    dplyr::mutate(rownum = which(row==mammal.taxa.name.order),
                  colnum = which(col==mammal.taxa.name.order),
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

#### Identify the locations of the ZFs, 9aaTADs and NLS in the AA & NT MSAs ####

locations.zf <- locate.zfs.in.alignment(combined.aa.aln.file, nt.aln.file, combined.taxa.name.order)
locations.9aaTAD <- locate.9aaTADs.in.alignment(combined.aa.aln.file, nt.aln.file, combined.taxa.name.order)
locations.NLS <- locate.NLS.in.alignment(combined.aa.aln.file, nt.aln.file, combined.taxa.name.order)

# Export the locations of the  ZFs,  9aaTADs, and NLS
write_tsv(locations.zf %>% 
            dplyr::select(sequence, aa_motif, start_ungapped, end_ungapped, 
                          start_gapped, end_gapped, start_nt_ungapped, end_nt_ungapped, 
                          start_nt_gapped, end_nt_gapped),
          "figure/locations.zf.tsv")
write_tsv(locations.9aaTAD %>% 
            dplyr::select(sequence, aa_motif, rc_score, start_ungapped, end_ungapped, 
                          start_gapped, end_gapped, start_nt_ungapped, end_nt_ungapped, 
                          start_nt_gapped, end_nt_gapped),
          "figure/locations.9aaTAD.tsv")
write_tsv(locations.NLS %>% 
            dplyr::select(sequence, aa_motif, type, posterior_prob, start_ungapped, end_ungapped, 
                          start_gapped, end_gapped, start_nt_ungapped, end_nt_ungapped, 
                          start_nt_gapped, end_nt_gapped),
          "figure/locations.NLS.tsv")

# We need to combine the full set of structure locations into an overlapping set
# to be plotted in a single row. Keep those that overlap in >=5 species

find.common.overlaps <- function(locations.data){
  ranges.9aaTAD <- IRanges(start=locations.data$start_gapped, end = locations.data$end_gapped, names = locations.data$sequence)
  ranges.9aaTAD.reduce <- IRanges::reduce(ranges.9aaTAD)

  n.sequences.in.range <- sapply(1:length(ranges.9aaTAD.reduce), function(i) length(subsetByOverlaps(ranges.9aaTAD, ranges.9aaTAD.reduce[i,])))
  as.data.frame(ranges.9aaTAD.reduce[n.sequences.in.range>4,]) %>%
    dplyr::mutate(motif_number = row_number())
}
ranges.ZF.common <- find.common.overlaps(locations.zf)
# Only keep the high confidence 9aaTADs for the track
ranges.9aaTAD.common <- find.common.overlaps(locations.9aaTAD[locations.9aaTAD$rc_score==100,])
ranges.9aaTAD.common$label <- LETTERS[1:nrow(ranges.9aaTAD.common)]
ranges.NLS.common <- find.common.overlaps(locations.NLS)

find.matching.range <-function(start, end, ranges){
  hit <- ranges[ranges$start<=start & ranges$end>=end,]$motif_number
  if(length(hit)==0) 0 else hit
}

# Annotate the individual ZFs, 9aaTADs and NLS with which consensus motif they belong to (0 if none)
locations.zf %<>% 
  dplyr::rowwise() %>%
  dplyr::mutate(motif_number = find.matching.range(start_gapped, end_gapped, ranges.ZF.common )) %>%
  dplyr::ungroup()

locations.9aaTAD %<>% 
  dplyr::rowwise() %>%
  dplyr::mutate(motif_number = find.matching.range(start_gapped, end_gapped, ranges.9aaTAD.common )) %>%
  dplyr::ungroup()

locations.NLS %<>% 
  dplyr::rowwise() %>%
  dplyr::mutate(motif_number = find.matching.range(start_gapped, end_gapped, ranges.NLS.common )) %>%
  dplyr::ungroup()

#### Export the residues within each ZF and 9aaTAD

# Create an Excel file with column filtering
create.xlsx = function(data, file.name, cols.to.fixed.size.font, cols.to.rich.text = NULL){
  oldOpt = options()
  options(xlsx.date.format="yyyy-mm-dd") # change date format
  wb = xlsx::createWorkbook(type = "xlsx")
  sh = xlsx::createSheet(wb)
  xlsx::addDataFrame(data, sh, row.names = F)
  cols.to.filter = paste0("A1:", LETTERS[ncol(data)], "1")
  xlsx::addAutoFilter(sh, cols.to.filter)
  xlsx::createFreezePane(sh, 2, 2, 2, 2) # freeze top row and first column
  cs <- xlsx::CellStyle(wb) + 
    xlsx::Font(wb,heightInPoints = 10, isBold = FALSE, name="Courier New")
  for(i in xlsx::getCells(xlsx::getRows(sh), colIndex=cols.to.fixed.size.font)){
    xlsx::setCellStyle(i, cs)
  }
  
  if(!is.null(cols.to.rich.text)) {
   for(col in cols.to.rich.text) set.rich.text.on.vv(wb, sh, col)

  }
  xlsx::autoSizeColumn(sh, 1:ncol(data))
  xlsx::saveWorkbook(wb, file=file.name)
  options(oldOpt)
}

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

# Full ZF motifs
locations.zf %>%
  dplyr::select(Sequence = sequence, motif_number, aa_motif) %>%
  dplyr::arrange(desc(motif_number)) %>% # so we display 13 - 1
  dplyr::mutate(motif_number = paste0("ZF_", motif_number)) %>%
  tidyr::pivot_wider(id_cols = Sequence, names_from = motif_number, values_from = aa_motif) %>%
  as.data.frame %>%
  dplyr::arrange(as.integer(Sequence)) %>%
  create.xlsx(., "figure/locations.zf.xlsx", cols.to.fixed.size.font = 2:14)

# Just the ZF contact bases
locations.zf %>%
  dplyr::select(Sequence = sequence, motif_number, contact_bases) %>%
  dplyr::arrange(desc(motif_number)) %>% # so we display 13 - 1
  dplyr::mutate(motif_number = paste0("ZF_", motif_number)) %>%
  tidyr::pivot_wider(id_cols = Sequence, names_from = motif_number, values_from = contact_bases) %>%
  as.data.frame %>%
  dplyr::arrange(as.integer(Sequence)) %>%
  create.xlsx(., "figure/locations.zf.contact_bases.xlsx", cols.to.fixed.size.font = 2:14)

locations.9aaTAD %>%
  dplyr::select(Sequence = sequence, motif_number, aa_motif) %>%
  dplyr::filter(motif_number >0) %>%
  dplyr::mutate(motif_number = paste0("9aaTAD_", LETTERS[motif_number])) %>%
  # summarise unique motifs
  dplyr::group_by(motif_number, aa_motif) %>%
  dplyr::summarise(Count = n(),
                   Sequences =case_when(Count > 6 ~ "Others",
                   .default = paste(Sequence, collapse = ", ")))%>%
  as.data.frame %>%
  create.xlsx(., "figure/locations.9aaTAD.unique.xlsx", cols.to.fixed.size.font = 2, cols.to.rich.text = 2)

locations.9aaTAD %>%
  dplyr::select(Sequence = sequence, motif_number, aa_motif) %>%
  dplyr::filter(motif_number >0) %>%
  dplyr::mutate(motif_number = paste0("9aaTAD_", LETTERS[motif_number])) %>%
  tidyr::pivot_wider(id_cols = c(Sequence), names_from = motif_number, 
                     values_from = aa_motif, values_fn = ~paste(.x, collapse = ", ")) %>%
  as.data.frame %>%
  create.xlsx(., "figure/locations.9aaTAD.xlsx", cols.to.fixed.size.font = 2:5, 
              cols.to.rich.text = 2:5)

locations.NLS %>%
  dplyr::select(Sequence = sequence, motif_number, aa_motif) %>%
  dplyr::mutate(motif_number = paste0("NLS_", motif_number)) %>%
  tidyr::pivot_wider(id_cols = Sequence, names_from = motif_number, values_from = aa_motif,
                     values_fn = ~paste(.x, collapse = ", ")) %>%
  as.data.frame %>%
  dplyr::arrange(as.integer(Sequence)) %>%
  dplyr::
  create.xlsx(., "figure/locations.NLS.xlsx", cols.to.fixed.size.font = 2:4)

#### Identify binding motifs of the ZFs in each species ####

# PWM prediction http://zf.princeton.edu/logoMain.php
# predict the ZF targets in each sequence via pwm_predict (http://zf.princeton.edu/download2.php)
old.wd <- getwd()
setwd("./bin/pwm_predict")
system2("./pwm_predict", "-l 20 ../../fasta/combined.aa.fas") # ensure all ZFs linked
setwd(old.wd)
filesstrings::move_files(files = c("fasta/combined.aa.pwm"),
                         destinations = c("pwm"),
                         overwrite = TRUE)
# Remove header and split the output PWMs to separate files
system2("cat", "pwm/combined.aa.pwm | grep -v ';' | split -l 5 - pwm/zf_")


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
taxon.data <- lapply(metadata.mammal$species, function(x) httr::GET( paste0("http://timetree.temple.edu/api/taxon/",curl::curl_escape(x)))  ) 
taxon.ids <- sapply(lapply(taxon.data, httr::content), function(x) x$taxon_id)
metadata.mammal$taxon_id <- taxon.ids

# Find the pairwise distances between each species
pairwise.species <- expand.grid(unique(metadata.mammal$taxon_id), unique(metadata.mammal$taxon_id)) %>%
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

#### Run codeml site model to check for site-specific selection ####

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



#### Run codeml branch-site model to look for selection specifically in Muroidea ####

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
rodent.ancestor.label <- "Muroidea ancestral Zfy"
rodent.plus.anc.nt.aln <- as.list.DNAbin(rodent.plus.anc.nt.aln)
names(rodent.plus.anc.nt.aln) <- c(names(ape.nt.aln), rodent.ancestor.label)

# Plot the MSA

# Find the tips under the ancestral node and get the labels
rodent.tip.labels <- tidytree::offspring( as_tibble(nt.aln.tree.nodes), .node = rodent.node+length(nt.aln.tree.nodes$tip.label), tiponly = T)


rodent.plus.anc.nt.aln.tidy <- tidy_msa(rodent.plus.anc.nt.aln) %>%
  dplyr::filter(name %in% c(rodent.ancestor.label, rodent.tip.labels$label))
rodent.plus.anc.nt.aln.tidy$name <- factor(rodent.plus.anc.nt.aln.tidy$name, 
                               levels = rev( c(mammal.taxa.name.order, rodent.ancestor.label))) # sort reverse to match tree

# Filter the zf locations and correct y locations
rodent.plus.anc.zf <- locations.zf %>%
  dplyr::filter(sequence %in% rodent.plus.anc.nt.aln.tidy$name) %>%
  dplyr::rowwise() %>%
  # Correct for different number of taxa in the y axis
  dplyr::mutate(i = which( levels(rodent.plus.anc.nt.aln.tidy$name) == sequence ) - (1+length(unique(mammal.taxa.name.order))- length(unique(rodent.plus.anc.nt.aln.tidy$name))))

# Filter 9aaTAD locations and correct the y locations
rodent.plus.anc.9aaTAD <- locations.9aaTAD %>%
  dplyr::filter(sequence %in% rodent.plus.anc.nt.aln.tidy$name) %>%
  dplyr::rowwise() %>%
  # Correct for different number of taxa in the y axis
  dplyr::mutate(i = which( levels(rodent.plus.anc.nt.aln.tidy$name) == sequence ) - (1+length(unique(mammal.taxa.name.order))- length(unique(rodent.plus.anc.nt.aln.tidy$name))))

rodent.plus.anc.NLS <- locations.NLS %>%
  dplyr::filter(sequence %in% rodent.plus.anc.nt.aln.tidy$name) %>%
  dplyr::rowwise() %>%
  # Correct for different number of taxa in the y axis
  dplyr::mutate(i = which( levels(rodent.plus.anc.nt.aln.tidy$name) == sequence ) - (1+length(unique(mammal.taxa.name.order))- length(unique(rodent.plus.anc.nt.aln.tidy$name))))



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

#### Calculate the conservation across the mammal AA domains for outgroup levels ####

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
# This should show all the 9aaTADs
aa.structure.plot <- ggplot()+
  geom_tile(data = locations.zf,     aes(x=(start_gapped+end_gapped)/2,
                                         width = (end_gapped-start_gapped),
                                         y=sequence),
            fill="grey", alpha=0.5)+
  geom_tile(data = locations.9aaTAD,     aes(x=(start_gapped+end_gapped)/2,
                                         width = (end_gapped-start_gapped),
                                         y=sequence,
                                         fill=round(rc_score)),
            alpha=0.9)+
  paletteer::scale_fill_paletteer_c("grDevices::Blues 3", "9aaTAD RC score (%)", direction = -1, limits = c(50, 100)) +
  geom_tile(data = locations.NLS,     aes(x=(start_gapped+end_gapped)/2,
                                         width = (end_gapped-start_gapped),
                                         y=sequence),
            fill="green", alpha=0.5)
aa.structure.plot <- annotate.structure.plot(aa.structure.plot, length(combined.taxa.name.order) + 1.5)

save.double.width("figure/aa.structure.png", aa.structure.plot, height = 120)

# Also create a trimmed down version that has only the 83, 92, 100% confidence 9aaTADs
aa.structure.confident.plot <- ggplot()+
  geom_tile(data = locations.zf,     aes(x=(start_gapped+end_gapped)/2,
                                         width = (end_gapped-start_gapped),
                                         y=sequence),
            fill="grey", alpha=0.5)+
  geom_tile(data = locations.9aaTAD[locations.9aaTAD$rc_score>80,],     aes(x=(start_gapped+end_gapped)/2,
                                             width = (end_gapped-start_gapped),
                                             y=sequence,
                                             fill=round(rc_score)),
            alpha=0.9)+
  paletteer::scale_fill_paletteer_c("grDevices::Blues 3", "9aaTAD RC score (%)", direction = -1, limits = c(50, 100)) +
  geom_tile(data = locations.NLS,     aes(x=(start_gapped+end_gapped)/2,
                                          width = (end_gapped-start_gapped),
                                          y=sequence),
            fill="green", alpha=0.5)
aa.structure.confident.plot <- annotate.structure.plot(aa.structure.confident.plot, length(combined.taxa.name.order) + 1.5)

save.double.width("figure/aa.structure.confident.png", aa.structure.confident.plot, height = 120)


#### Plot the conservation of hydrophobicity across mammal/outgroup AA MSA ####

msa.aa.aln.tidy.hydrophobicity <- do.call(rbind, mapply(calc.hydrophobicity, aa=combined.aa.aln@unmasked, 
                                                sequence.name = names(combined.aa.aln@unmasked), 
                                                window.size = 9,
                                                SIMPLIFY = FALSE))  %>%
  dplyr::mutate(sequence = factor(sequence, levels = rev(combined.taxa.name.order))) # sort reverse to match tree

n.taxa <- length(combined.taxa.name.order) + 1.5

hydrophobicity.plot <- ggplot()+
  geom_tile(data=msa.aa.aln.tidy.hydrophobicity,  aes(x = position_gapped, y = sequence, fill=hydrophobicity_smoothed))+
  scale_fill_paletteer_c("ggthemes::Classic Red-Blue", direction = -1, limits = c(0, 1))+
  labs(fill="Hydrophobicity (9 site average)")
hydrophobicity.plot <- annotate.structure.plot(hydrophobicity.plot, n.taxa)
save.double.width("figure/hydrophobicity.convervation.tree.png", hydrophobicity.plot, height = 120)

#### Plot conservation of charge across mammal/outgroup AA MSA ####

msa.aa.aln.tidy.charge <- do.call(rbind, mapply(calc.charge, aa=combined.aa.aln@unmasked, 
                                                sequence.name = names(combined.aa.aln@unmasked), 
                                                window.size = 9,
                                                SIMPLIFY = FALSE))  %>%
  dplyr::mutate(sequence = factor(sequence, levels = rev(combined.taxa.name.order))) # sort reverse to match tree

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

#### What are the hydrophobic patches in exons 3 and 5 that are not 9aaTADs?

exon3.patch <- msa.aa.aln.tidy.hydrophobicity %>%
  dplyr::filter(position_gapped>mouse.exons$start_aa[3]+42 & position_gapped<mouse.exons$end_aa[3]-2)

exon3.hydro.plot <- ggplot()+
  # Draw the charges per sequence
  geom_tile(data=exon3.patch,  aes(x = position_gapped, y = sequence, fill=hydrophobicity_smoothed))+
  scale_fill_paletteer_c("ggthemes::Classic Red-Blue", direction = -1, limits = c(0, 1))+
  labs(fill="Hydrophobicity (9 site average)")

# Get the region from the msa
exon3.patch.table <- do.call(rbind, lapply(metadata.combined$common.name,  function(i) list("Sequence" = i, "exon_3_285-305"= as.character(msa.aa.aln@unmasked[[i]][285:305])))) %>%
  as.data.frame %>%
  dplyr::mutate(Sequence = factor(Sequence, levels = combined.taxa.name.order)) %>%
  dplyr::arrange(as.integer(Sequence))

exon5.patch <- msa.aa.aln.tidy.hydrophobicity %>%
  dplyr::filter(position_gapped>mouse.exons$start_aa[5]+38 & position_gapped<mouse.exons$end_aa[5]+3)
exon5.hydro.plot <- ggplot()+
  # Draw the charges per sequence
  geom_tile(data=exon5.patch,  aes(x = position_gapped, y = sequence, fill=hydrophobicity_smoothed))+
  scale_fill_paletteer_c("ggthemes::Classic Red-Blue", direction = -1, limits = c(0, 1))+
  labs(fill="Hydrophobicity (9 site average)")

exon5.patch.table <- do.call(rbind, lapply(metadata.combined$common.name,  function(i) list("Sequence" = i, "exon_5_399-418"= as.character(msa.aa.aln@unmasked[[i]][399:418])))) %>%
  as.data.frame %>%
  dplyr::mutate(Sequence = factor(Sequence, levels = combined.taxa.name.order)) %>%
  dplyr::arrange(as.integer(Sequence))

exon.patch.table <- merge(exon3.patch.table, exon5.patch.table, by=c("Sequence")) %>%
  create.xlsx(., "figure/hydrophobic_patches.xlsx", cols.to.fixed.size.font = 2:3)

# Get consensus strings for the two patches
paste(s2c(consensusString(msa.aa.aln))[285:305], collapse = "")
paste(s2c(consensusString(msa.aa.aln))[399:418], collapse = "")

#### ####


# testing of side-by-side plots
# p <- ggtree(combined.outgroup.tree) + geom_tiplab(size=3, aes(col = group), align=F)
 # msaplot(p, combined.aa.aln.file, offset=0.4, width=2)+theme(legend.position = "none")

# p %>% aplot::insert_right(msa.combined.aa.plot)




# ggmsa::treeMSA_plot(ggtree(combined.outgroup.tree)+geom_tiplab(size=3, aes(col = group)), msa.outgroup.aln.tidy, color = "Chemistry_AA",
#                     border=NA, font=NULL)
#  