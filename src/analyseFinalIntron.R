# Analyse COX1
# We want to compare the evolutionary rate of the Zfy final intron versus exons

#### Imports #####

source("src/functions.R")
load.packages()

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
names(zfy.read) <- gsub(" ", "_", str_extract(names(zfy.read), "\\[(.*)\\]", group = 1))
ape::write.FASTA(zfy.read, file = FILES$final.intron.zfy.nt.fas)

# Read all unaligned sequence files with .fa extension
zfx.files <- list.files(path = "fasta/final.intron.zfx", pattern = "*.fa$", 
                        include.dirs = T, full.names = T)

# Combine to single FASTA file
zfx.read  <- lapply(zfx.files, ape::read.FASTA)
zfx.read <- do.call(c, lapply(zfx.read,\(x)x))

# Reduce to simple names
names(zfx.read) <- gsub(" ", "_", str_extract(names(zfx.read), "\\[(.*)\\]", group = 1))
ape::write.FASTA(zfx.read, file = FILES$final.intron.zfx.nt.fas)

# Add to a single combined file for checking overall alignment
zf.all <- c(zfx.read, zfy.read)
ape::write.FASTA(zf.all, file = FILES$final.intron.nt.fas)

#### Align sequences #####

# Use MUSCLE to align  - non coding
run.muscle(FILES$final.intron.zfy.nt.fas, FILES$final.intron.zfy.nt.aln)
run.muscle(FILES$final.intron.zfx.nt.fas, FILES$final.intron.zfx.nt.aln)
run.muscle(FILES$final.intron.nt.fas,     FILES$final.intron.nt.aln)

#### Make ML tree with distances #####

final.intron.zfy.nt.aln.tree <- run.iqtree(FILES$final.intron.zfy.nt.aln) %>%
        ape::read.tree(.) %>%
        reroot.tree(.,"Platypus", position = 0.015)

final.intron.zfx.nt.aln.tree <- run.iqtree(FILES$final.intron.zfx.nt.aln) %>%
        ape::read.tree(.) %>%
        reroot.tree(.,"Platypus", position = 0.015)

final.intron.nt.aln.tree <- run.iqtree(FILES$final.intron.nt.aln) %>%
        ape::read.tree(.) %>%
        reroot.tree(., "Platypus", position = 0.015)

#### Plot the trees ####

final.intron.zfy.nt.aln.tree.plot <- plot.tree(final.intron.zfy.nt.aln.tree)  + xlim(0, 1.5) + labs(title = "ZFY")

final.intron.zfx.nt.aln.tree.plot <- plot.tree(final.intron.zfx.nt.aln.tree) + xlim(0, 1.5)+ labs(title = "ZFX")
save.double.width("figure/final.intron.tree.png", final.intron.zfx.nt.aln.tree.plot/final.intron.zfy.nt.aln.tree.plot)

final.intron.nt.aln.tree.plot <- plot.tree(final.intron.nt.aln.tree)  + xlim(0, 1.5) + labs(title = "Combined")
save.double.width("figure/final.intron.combined.tree.png", final.intron.nt.aln.tree.plot)


