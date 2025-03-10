# Analyse COX1
# We want to compare the evolutionary rate of the Zfy final intron versus exons

# This requires divvier 1.01 in a conda environment
# RepeatMasker 4.1.4

#### Imports #####

source("src/functions.R")
load.packages()

METADATA <- prepare.fas.files() # load FASTA files and write metadata table

# Clear previous analyses
filesstrings::dir.remove(c("aln/final.intron.zfy", "aln/final.intron.zfx", "aln/final.intron"))
filesstrings::create_dir(c("aln/final.intron.zfy", "aln/final.intron.zfx", "aln/final.intron"))

#### Read sequences #####

# Read all unaligned sequence files with .fa extension
zfy.files <- list.files(path = "fasta/final.intron.zfy", pattern = "*.fa$", 
                       include.dirs = T, full.names = T)

# Combine to single FASTA file
zfy.read  <- lapply(zfy.files, ape::read.FASTA)
# # Remove existing names from sequences
zfy.read <- do.call(c, lapply(zfy.read,\(x)x))

# Reduce to simple names
names(zfy.read) <- gsub(" ", "_", stringr::str_extract(names(zfy.read), "\\[(.*)\\]", group = 1))
ape::write.FASTA(zfy.read, file = FILES$final.intron.zfy.nt.fas)

# Read all unaligned sequence files with .fa extension
zfx.files <- list.files(path = "fasta/final.intron.zfx", pattern = "*.fa$", 
                        include.dirs = T, full.names = T)

# Combine to single FASTA file
zfx.read  <- lapply(zfx.files, ape::read.FASTA)
zfx.read <- do.call(c, lapply(zfx.read,\(x)x))

# Reduce to simple names
names(zfx.read) <- gsub(" ", "_", stringr::str_extract(names(zfx.read), "\\[(.*)\\]", group = 1))
ape::write.FASTA(zfx.read, file = FILES$final.intron.zfx.nt.fas)

#### Write metadata to supplement ####

list.files(path = c("fasta/final.intron.zfx","fasta/final.intron.zfy" ), 
           pattern = "*.fa$", include.dirs = T, full.names = T) %>%
  lapply(., read.fasta) %>%
  read.metadata(.) %>%
  dplyr::mutate(Species_common_name = str_replace(common.name, "_Z[F|f][X|x].*", ""),
                Species_common_name = str_replace(Species_common_name, "(_putative)?(_|-)Z[F|f][X|x|Y|y].*", "")) %>%
  dplyr::rename(Accession = accession, Species = species, Group = group,
                Name_in_figures = common.name, Description = original.name  ) %>%
  dplyr::distinct() %>%
  dplyr::select(Accession, Species, Group, Name_in_figures, Description ) %>%
  dplyr::arrange(Group, Species) %>%
  create.xlsx(., "figure/accessions.final.intron.supplement.xlsx")

#### Remove repeats via RepeatMasker ####

cat("Running RepeatMasker\n")
# Repeats and low complexity regions make the alignments harder. Remove as much
# as possible before we align
system2("RepeatMasker", paste("-s -species mammals", 
                              "-dir aln/final.intron.zfy/zfy.masked", 
                              "-libdir /usr/local/RepeatMasker/Libraries/ ", 
                              FILES$final.intron.zfy.nt.fas))
# Remove masked repeats and blank lines
system2("sed", paste("'s/NN\\+//g' aln/final.intron.zfy/zfy.masked/final.intron.zfy.nt.fas.masked | sed -e '/^$/d' > aln/final.intron.zfy/final.intron.zfy.nt.trim.fas"))

system2("RepeatMasker", paste("-s -species mammals", 
                              "-dir aln/final.intron.zfx/zfx.masked", 
                              "-libdir /usr/local/RepeatMasker/Libraries/ ", 
                              FILES$final.intron.zfx.nt.fas))
system2("sed", paste("'s/NN\\+//g' aln/final.intron.zfx/zfx.masked/final.intron.zfx.nt.fas.masked | sed -e '/^$/d' > aln/final.intron.zfx/final.intron.zfx.nt.trim.fas"))


# Create a single combined ZFX/Y file for overall alignment
zf.all <- c(zfx.read, zfy.read)
ape::write.FASTA(zf.all, file = FILES$final.intron.nt.fas)

zfy.trim <- ape::read.FASTA("aln/final.intron.zfy/final.intron.zfy.nt.trim.fas")
zfx.trim <- ape::read.FASTA("aln/final.intron.zfx/final.intron.zfx.nt.trim.fas")
zfx.trim$Platypus_ZFX <- NULL # remove duplicate outgroups before combining
zfx.trim$Koala_ZFX <- NULL
zfx.trim$Opossum_ZFX <- NULL
zfx.trim$Australian_echidna_ZFX <- NULL

zf.trim  <- c(zfx.trim, zfy.trim)
ape::write.FASTA(zf.trim, file = "aln/final.intron/final.intron.nt.trim.fas")

#### Align sequences #####

# Use MUSCLE to align  - non coding
cat("Aligning sequences\n")
run.muscle("aln/final.intron.zfy/final.intron.zfy.nt.trim.fas", FILES$final.intron.zfy.nt.aln)
run.muscle("aln/final.intron.zfx/final.intron.zfx.nt.trim.fas", FILES$final.intron.zfx.nt.aln)
run.muscle("aln/final.intron/final.intron.nt.trim.fas", FILES$final.intron.nt.aln)

# Filter a Biostrings multiple sequence alignment to mask sites with more than a
# given number of gaps. Returns the masked Biostrings alignment
# aln - the Biostrings alignment
# max.gaps - the maximum number of sequences with a gap for a site to be kept
# min.sequences - the minimum number of sequences with a nucleotide for a site 
#                 to be kept. Overrides max.gaps if specified
filter.alignment.on.gaps <- function(aln, max.gaps = 4, min.sequences = NA){
  
  n.sequences <- length(aln@unmasked)
  min.sequences <- ifelse(is.na(min.sequences), 0, min.sequences)
  
  # Find and count sites with gaps
  matches <- do.call(c, sapply(aln@unmasked, \(x) Biostrings::start( Biostrings::matchPattern("-",x))))
  # Mask gapped sites in the alignment
  max.gaps.mask <- as.integer(names(table(matches)[table(matches) > max.gaps]))
  min.seqs.mask <- as.integer(names(table(matches)[table(matches) > n.sequences-min.sequences]))
  
  # Use min.sequences if specified
  if(min.sequences>0){
    colmask(aln) <- IRanges(start=min.seqs.mask, width = 1)
  } else {
    colmask(aln) <- IRanges(start=max.gaps.mask, width = 1)
  }

  aln
}

# Convert a Biostrings MSA to a character list format that can be exported as FASTA
biostrings.aln.to.list <- function(aln) as.list(apply(as.matrix(aln), 1, paste, collapse=""))

# Filter ZFY alignment
zfy.aln <- Biostrings::readDNAMultipleAlignment(FILES$final.intron.zfy.nt.aln, format="fasta")
zfy.aln <- filter.alignment.on.gaps(zfy.aln, min.sequences = 8)
zfy.aln.vec <- biostrings.aln.to.list(zfy.aln)
seqinr::write.fasta(zfy.aln.vec, names(zfy.aln.vec), "aln/final.intron.zfy/final.intron.zfy.nt.filt.aln")

# Filter ZFX alignment
zfx.aln <- Biostrings::readDNAMultipleAlignment(FILES$final.intron.zfx.nt.aln, format="fasta")
zfx.aln <- filter.alignment.on.gaps(zfx.aln, min.sequences = 8)
zfx.aln.vec <- biostrings.aln.to.list(zfx.aln)
seqinr::write.fasta(zfx.aln.vec, names(zfx.aln.vec), "aln/final.intron.zfx/final.intron.zfx.nt.filt.aln")

# Filter combined ZFX/Y alignment
zf.aln <- Biostrings::readDNAMultipleAlignment(FILES$final.intron.nt.aln, format="fasta")
zf.aln <- filter.alignment.on.gaps(zf.aln, min.sequences = 8)
zf.aln.vec <- biostrings.aln.to.list(zf.aln)
seqinr::write.fasta(zf.aln.vec, names(zf.aln.vec), "aln/final.intron/final.intron.nt.filt.aln")


#### Run divvier to improve high-confidence homologies in alignments ####

final.intron.zfy.nt.divvy.aln <- run.divvier(FILES$final.intron.zfy.nt.aln, "-mincol 4")
final.intron.zfx.nt.divvy.aln <- run.divvier(FILES$final.intron.zfx.nt.aln, "-mincol 4")
final.intron.nt.divvy.aln <- run.divvier(FILES$final.intron.nt.aln, "-mincol 4")

#### Make species phylogenies for Zfx and Zfy #####

# Specify species trees for each of Zfx and Zfy
# Note that not all species have intron info available - make a new species tree
# The alignment does not include marsupials etc so drop these
# ape::read.tree("aln/zfy_only/zfy.nt.species.nwk") %>%
#   ape::write.tree(., "aln/final.intron.zfy/zfy.nt.species.nwk")
# 
# ape::read.tree("aln/zfx_only/zfx.nt.species.nwk") %>%
#   ape::write.tree(., "aln/final.intron.zfx/zfx.nt.species.nwk")


#### Make ML tree based on raw alignments using species phylogienies #####

# cat("Creating raw alignment speicies trees\n")
# # Make the tree using the species phylogeny
# final.intron.zfy.nt.aln.treefile <- run.iqtree(FILES$final.intron.zfy.nt.aln, 
#                                                "-nt AUTO", # number of threads
#                                                "-keep-ident",
#                                                "-te aln/final.intron.zfy/zfy.nt.species.nwk" # user tree guide
# ) 
# 
# 
# final.intron.zfy.nt.aln.tree <- ape::read.tree(FILES$final.intron.zfy.nt.aln.treefile) %>%
#   reroot.tree(., c("Platypus_ZFX", "Autralian_echidna_ZFX"), position = 0.015)
# 
# final.intron.zfx.nt.aln.treefile <- run.iqtree(FILES$final.intron.zfx.nt.aln, 
#                                                "-nt AUTO", # number of threads
#                                                "-keep-ident",
#                                                "-te aln/final.intron.zfx/zfx.nt.species.nwk" # user tree guide
# ) 
# final.intron.zfx.nt.aln.tree <- ape::read.tree(FILES$final.intron.zfx.nt.aln.treefile) %>%
#   reroot.tree(., c("Platypus_ZFX", "Autralian_echidna_ZFX"), position = 0.015)

#### Make ML tree based on raw alignments #####

cat("Creating raw alignment trees\n")
# Make ZFY tree without species phylogeny
final.intron.zfy.nt.aln.treefile <- run.iqtree("aln/final.intron.zfy/final.intron.zfy.nt.filt.aln", 
                                               "-nt AUTO", # number of threads
                                               "-bb 1000",
                                               "-alrt 1000",
                                               "-keep-ident"
) 

final.intron.zfy.nt.aln.tree <- ape::read.tree(final.intron.zfy.nt.aln.treefile) %>%
        reroot.tree(., c("Platypus_ZFX", "Australian_echidna_ZFX"), position = 0.015)

# Make ZFX tree without species phylogeny
final.intron.zfx.nt.aln.treefile <- run.iqtree("aln/final.intron.zfx/final.intron.zfx.nt.filt.aln", 
                                               "-nt AUTO", # number of threads
                                               "-bb 1000",
                                               "-alrt 1000",
                                               "-keep-ident"
)

final.intron.zfx.nt.aln.tree <- ape::read.tree(final.intron.zfx.nt.aln.treefile) %>%
        reroot.tree(., c("Platypus_ZFX", "Australian_echidna_ZFX"), position = 0.015)

final.intron.nt.aln.treefile <- run.iqtree("aln/final.intron/final.intron.nt.filt.aln", 
                                           "-nt AUTO",
                                           "-bb 1000",
                                           "-alrt 1000",
                                           "-keep-ident"
) 
final.intron.nt.aln.tree <- ape::read.tree(final.intron.nt.aln.treefile) %>%
  reroot.tree(., c("Platypus_ZFX", "Australian_echidna_ZFX"), position = 0.015)
mammal.gene.groups <- split(METADATA$combined$common.name, METADATA$combined$group)
final.intron.nt.aln.tree <- tidytree::groupOTU(final.intron.nt.aln.tree, mammal.gene.groups, group_name = "group")

#### Make ML tree based on divvied alignments ####
cat("Creating divvied alignment trees\n")

# ZFY
final.intron.zfy.nt.divvy.aln.treefile <- run.iqtree(final.intron.zfy.nt.divvy.aln,
                                                     "-nt AUTO", # number of threads
                                                     "-bb 1000",
                                                     "-alrt 1000",
                                                     "-keep-ident"
)

final.intron.zfy.nt.divvy.aln.tree <- ape::read.tree(FILES$final.intron.zfy.nt.aln.divvy.aln.treefile) %>%
  reroot.tree(., c("Platypus_ZFX", "Australian_echidna_ZFX"), position = 0.015)

# ZFX
final.intron.zfx.nt.divvy.aln.treefile <- run.iqtree(final.intron.zfx.nt.divvy.aln,
                                                     "-nt AUTO", # number of threads
                                                     "-bb 1000",
                                                     "-alrt 1000",
                                                     "-keep-ident"
)
final.intron.zfx.nt.divvy.aln.tree <- ape::read.tree(FILES$final.intron.zfx.nt.aln.divvy.aln.treefile) %>%
  reroot.tree(.,c("Platypus_ZFX", "Australian_echidna_ZFX"), position = 0.015)

# ZFX/Y
final.intron.zfx.nt.divvy.aln.treefile <- run.iqtree(final.intron.nt.divvy.aln,
                                                     "-nt AUTO", # number of threads
                                                     "-bb 1000",
                                                     "-alrt 1000",
                                                     "-keep-ident"
)
final.intron.nt.divvy.aln.tree <- ape::read.tree(FILES$final.intron.nt.aln.divvy.aln.treefile) %>%
  reroot.tree(.,c("Platypus_ZFX", "Australian_echidna_ZFX"), position = 0.015)


#### Plot the trees ####

# Raw
final.intron.zfy.nt.aln.tree.plot <- plot.tree(final.intron.zfy.nt.aln.tree)  + xlim(-0.1, 2) + labs(title = "ZFY")
final.intron.zfx.nt.aln.tree.plot <- plot.tree(final.intron.zfx.nt.aln.tree) + xlim(0, 2)+ labs(title = "ZFX")
save.double.width("figure/final.intron.tree.png", final.intron.zfx.nt.aln.tree.plot/final.intron.zfy.nt.aln.tree.plot)

# Divvied
final.intron.zfy.nt.divvy.aln.tree.plot <- plot.tree(final.intron.zfy.nt.divvy.aln.tree)  + xlim(0, 2) + labs(title = "ZFY (divvied)")
final.intron.zfx.nt.divvy.aln.tree.plot <- plot.tree(final.intron.zfx.nt.divvy.aln.tree) + xlim(0, 2)+ labs(title = "ZFX (divvied)")
save.double.width("figure/final.intron.divvy.tree.png", final.intron.zfy.nt.divvy.aln.tree.plot/final.intron.zfx.nt.divvy.aln.tree.plot)

# Raw ZFX/Y
final.intron.nt.aln.tree.plot <- plot.tree(final.intron.nt.aln.tree, col= "group")  + xlim(-0.1, 1.5) + labs(title = "ZFX/Y final intron")
save.double.width("figure/final.intron.zfx.zfy.tree.png", final.intron.nt.aln.tree.plot)

# Divvied ZFX/Y
final.intron.nt.divvy.aln.tree.plot <- plot.tree(final.intron.nt.divvy.aln.tree)  + xlim(0, 2) + labs(title = "ZFX/Y (divvied)")
save.double.width("figure/final.intron.zfx.zfy.divvy.tree.png", final.intron.nt.divvy.aln.tree.plot)

# Panel figure
# combined.plot <- (final.intron.zfy.nt.aln.tree.plot + final.intron.zfy.nt.divvy.aln.tree.plot) / (final.intron.zfx.nt.aln.tree.plot + final.intron.zfx.nt.divvy.aln.tree.plot)
# save.double.width("figure/final.intron.combined.tree.png", combined.plot)

#### End ####
cat("Done!")
