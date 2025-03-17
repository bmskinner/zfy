# Analysis pipeline for ZFX/ZFY evolutionary analysis

# The following binaries must be supplied in ./bin:
# macse_v2.07.jar
# muscle5.1.win64.exe / muscle5.1.linux_intelx64
# nlstradamus.pl
# pwm_predict/pwm_predict

if(!file.exists("./bin/macse_v2.07.jar")) stop("MACSE not found in ./bin/\nwget https://www.agap-ge2pop.org/wp-content/uploads/macse/releases/macse_v2.07.jar")
if(!file.exists("./bin/nlstradamus.pl")) stop("NLStradamus not found in ./bin/\nwget http://www.moseslab.csb.utoronto.ca/NLStradamus/NLStradamus/NLStradamus.1.8.tar.gz")
if(!file.exists("./bin/pwm_predict/pwm_predict")) stop("pwm_predict not found in ./bin/\nwget https://zf.princeton.edu/downloads/pwm_predict.1.0.tar.gz")
if(!file.exists("./bin/muscle5.1.linux_intel64")) stop("Muscle not found in ./bin/\nwget https://github.com/rcedgar/muscle/releases/download/v5.1/muscle5.1.linux_intel64")

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
#### Imports #####

source("src/functions.R")
# load.packages() - automatically loaded by functions.R

source("src/find9aaTADs.R")
source("src/findZF.R")
source("src/calcCharge.R")
source("src/calcHydrophobicity.R")

cat("Packages loaded\n")

#### Create output directory structure #####

# Create directories that remain constant
filesstrings::create_dir("paml")
filesstrings::create_dir("paml/site-specific")
filesstrings::create_dir("paml/exon_1_3-6")
filesstrings::create_dir("paml/exon_2")
filesstrings::create_dir("paml/exon_7")

# Clear previous analyses if present
cat(timestamp(), "Removing existing output directories\n")
filesstrings::remove_dir(c("aln", "figure"))

# Create missing dirs if needed
cat(timestamp(), "Creating output directories\n")
create.output.dir <- function(dir.path){
  if(!filesstrings::create_dir(dir.path)) stop(paste0("Could not create output directory '", dir.path, "'"))
}

create.output.dir("aln")
create.output.dir("aln/mammal")
create.output.dir("aln/combined")
create.output.dir("aln/exons")
create.output.dir("aln/hyphy")
create.output.dir("aln/nls")
create.output.dir("aln/pwm")
create.output.dir("aln/zfx_only")
create.output.dir("aln/zfy_only")
create.output.dir("aln/final.intron.zfx")
create.output.dir("aln/final.intron.zfy")
create.output.dir("aln/final.intron")
create.output.dir("figure")

# Output the packages and versions used
writeLines(capture.output(sessionInfo()), "figure/session_info.txt")

# Read species-specific FASTA files, combine into files for analysis.
# Make metadata from file names and FASTA headers, export in Excel format.
METADATA <- prepare.fas.files()


#### Run combined mammal/outgroup AA alignment ####
cat(timestamp(), "Creating mammal plus outgroup alignments\n")
# Expect java on the PATH. Macse download from https://www.agap-ge2pop.org/macsee-pipelines/
# Direct download link:
# https://www.agap-ge2pop.org/wp-content/uploads/macse/releases/macse_v2.07.jar
# Use macse to align the nt files for combined species
run.macse(FILES$combined.nt.fas, "aln/combined/combined")

#### Run mammal NT alignment guided by AA #####
cat(timestamp(), "Creating mammal alignments\n")
# Run a codon aware alignment with MACSE
run.macse(FILES$mammal.nt.fas, "aln/mammal/mammal")

# Read all alignments in ape and Biostrings formats
ALIGNMENTS <- read.alignments()

# Identify the coordinates of the exon boundaries in the gapped alignments
# Based on the mouse Zfy1 sequence
mouse.exons <- find.exons()

# Create NEXUS format alignments with partition information for the mammal and
# combined alignments 
# We consider three partitions: Exons 1 & 3-6 Exon 2 Exon 7
write_file(paste0("BEGIN SETS;\n",
                  "\tcharset exon13456 = ", mouse.exons$start_nt_codon_offset_mammal[1], "-", mouse.exons$end_nt_codon_offset_mammal[1],
                  " ", mouse.exons$start_nt_codon_offset_mammal[3], "-", mouse.exons$end_nt_codon_offset_mammal[6],";\n",
                  "\tcharset exon2 = ", mouse.exons$start_nt_codon_offset_mammal[2], "-", mouse.exons$end_nt_codon_offset_mammal[2],";\n",
                  "\tcharset exon7 = ", mouse.exons$start_nt_codon_offset_mammal[7], "-", mouse.exons$end_nt_codon_offset_mammal[7],";\n",
                  "END;\n"),
           FILES$mammal.nt.partition)

write_file(paste0("BEGIN SETS;\n",
                  "\tcharset exon13456 = ", mouse.exons$start_nt_codon_offset_combined[1], "-", mouse.exons$end_nt_codon_offset_combined[1],
                  " ", mouse.exons$start_nt_codon_offset_combined[3], "-", mouse.exons$end_nt_codon_offset_combined[6],";\n",
                  "\tcharset exon2 = ", mouse.exons$start_nt_codon_offset_combined[2], "-", mouse.exons$end_nt_codon_offset_combined[2],";\n",
                  "\tcharset exon7 = ", mouse.exons$start_nt_codon_offset_combined[7], "-", mouse.exons$end_nt_codon_offset_combined[7],";\n",
                  "END;\n"),
           FILES$combined.nt.partition)

# Resave the alignment in NEXUS format with the partitions
alg2nex(FILES$mammal.nt.aln, format = "fasta", interleaved = FALSE, gap = "-",
        missing = "?", partition.file = FILES$mammal.nt.partition)

alg2nex(FILES$combined.nt.aln, format = "fasta", interleaved = FALSE, gap = "-",
        missing = "?", partition.file = FILES$combined.nt.partition)


#### Create combined mammal/outgroup AA trees ####
cat(timestamp(), "Creating trees\n")
# All outgroups
FILES$combined.aa.aln.treefile <- run.iqtree(FILES$combined.aa.aln, 
                                             "-bb 1000", # number of bootstrap replicates
                                             "-alrt 1000") # number of replicates to perform SH-like approximate likelihood ratio test (SH-aLRT) 
# Mammal only
FILES$mammal.aa.aln.treefile <- run.iqtree(FILES$mammal.aa.aln, 
                                             "-bb 1000", # number of bootstrap replicates
                                             "-alrt 1000") # number of replicates to perform SH-like approximate likelihood ratio test (SH-aLRT) 


#### Create CDS NT trees #####

# Make ML tree and reconstruct ancestral sequences
# Expect iqtree on the PATH.
# Note model testing is automatically performed in v1.5.4 onwards

# All outgroups
FILES$combined.nt.aln.treefile <- run.iqtree(FILES$combined.nt.aln, 
                                             "-bb 1000", # number of bootstrap replicates
                                             "-alrt 1000", # number of replicates to perform SH-like approximate likelihood ratio test (SH-aLRT) 
                                             "-asr" # Ancestral sequence resconstruction
                                             )

# Mammal only
FILES$mammal.nt.aln.treefile <- run.iqtree(FILES$mammal.nt.aln,
                                           "-bb 1000", # number of bootstrap replicates
                                           "-alrt 1000", # number of replicates to perform SH-like approximate likelihood ratio test (SH-aLRT) 
                                           "-asr" # Ancestral sequence resconstruction
)

#### Make individual mammal exon NT trees ####

# Aim here is to look for signs of gene conversion in the final exon versus others
create.exon.alignment <- function(exon.aln, name){
  exon.aln <- ape::del.colgapsonly(exon.aln, threshold = 0.2) # remove columns with >20% gaps
  exon.aln.file <- paste0("aln/exons/", name)
  
  ape::write.FASTA(exon.aln, file = exon.aln.file)
  
  system2("iqtree", paste("-s ", exon.aln.file,
                          "-bb 1000", # number of bootstrap replicates
                          "-alrt 1000", # number of replicates to perform SH-like approximate likelihood ratio test (SH-aLRT)
                          "-nt AUTO"), # number of threads
          stdout = paste0(exon.aln.file, ".iqtree.log"),
          stderr = paste0(exon.aln.file, ".iqtree.log"))
}

exon.2.aln <- as.matrix(ALIGNMENTS$nt.mammal.ape)[,mouse.exons$start_nt_mammal[2]:mouse.exons$end_nt_mammal[2]]
create.exon.alignment(exon.2.aln, "exon_2.aln")

exon.7.aln <- as.matrix(ALIGNMENTS$nt.mammal.ape)[,mouse.exons$start_nt_mammal[7]:mouse.exons$end_nt_mammal[7]]
create.exon.alignment(exon.7.aln, "exon_7.aln")

# All except exon 7
exon.1.6.aln <- as.matrix(ALIGNMENTS$nt.mammal.ape)[,mouse.exons$start_nt_mammal[1]:mouse.exons$end_nt_mammal[6]]
create.exon.alignment(exon.1.6.aln, "exon_1-6.aln")

# All except exon 2 and 7
exon1.3_6.locs <- c(mouse.exons$start_nt_codon_offset_mammal[1]:mouse.exons$end_nt_codon_offset_mammal[1],
                    mouse.exons$start_nt_codon_offset_mammal[3]:mouse.exons$end_nt_codon_offset_mammal[6])
exon1.3_6.aln <- as.matrix(ALIGNMENTS$nt.mammal.ape)[,exon1.3_6.locs]
create.exon.alignment(exon1.3_6.aln, "exon_1.3-6.aln")

#### Identify binding motifs of the ZFs in each species ####

run.pwm.predict()

#### Find and export structural features ####
cat(timestamp(), "Identifying structural features\n")
# Read the combined tree and find the taxa order
combined.outgroup.tree <- read.combined.outgroup.tree(FILES$combined.aa.aln.treefile)
combined.aa.tree <- plot.tree(combined.outgroup.tree, col = "group")
combined.taxa.name.order <- ggtree::get_taxa_name(combined.aa.tree) 

# Find and export ZFs
locations.zf <- locate.zfs.in.alignment(FILES$combined.aa.aln, FILES$combined.nt.aln, combined.taxa.name.order)
readr::write_tsv(locations.zf, "aln/locations.zf.combined.tsv")

# Find and export 9aaTADS
locations.9aaTAD <- locate.9aaTADs.in.alignment(FILES$combined.aa.aln, FILES$combined.nt.aln, combined.taxa.name.order)
readr::write_tsv(locations.9aaTAD, "aln/locations.9aaTAD.combined.tsv")

# Find and export NLS
locations.NLS <- locate.NLS.in.alignment(FILES$combined.aa.aln, FILES$combined.nt.aln, combined.taxa.name.order)
readr::write_tsv(locations.NLS, "aln/locations.NLS.combined.tsv")

# Also find the locations of these features in the mammal-only alignments

mammal.outgroup.tree <- read.mammal.outgroup.tree(FILES$mammal.aa.aln.treefile)
mammal.aa.tree <- plot.tree(mammal.outgroup.tree, col = "group")
mammal.taxa.name.order <- ggtree::get_taxa_name(mammal.aa.tree) 

locations.zf <- locate.zfs.in.alignment(FILES$mammal.aa.aln, FILES$mammal.nt.aln, mammal.taxa.name.order)
readr::write_tsv(locations.zf, "aln/locations.zf.mammal.tsv")

# Find and export 9aaTADS
locations.9aaTAD <- locate.9aaTADs.in.alignment(FILES$mammal.aa.aln, FILES$mammal.nt.aln, mammal.taxa.name.order)
readr::write_tsv(locations.9aaTAD, "aln/locations.9aaTAD.mammal.tsv")

# Find and export NLS
locations.NLS <- locate.NLS.in.alignment(FILES$mammal.aa.aln, FILES$mammal.nt.aln, mammal.taxa.name.order)
readr::write_tsv(locations.NLS, "aln/locations.NLS.mammal.tsv")

#### Create species trees for ZFX and ZFY using TimeTree Newick tree ####

# Create a species list for use in bulk TimeTree website. Note that this may not 
# retrieve the full species tree - check manually. The manual output is saved
# to ./species_names.nwk
write.table(unique(METADATA$mammal$species), file = "figure/species_names.tsv", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# Read the manually created species tree
# Replace the Latin names with common names, and ZFX/ZFY suffix to make
# gene specific species trees
species.tree <- ape::read.tree("species_names.nwk")
species.tree$tip.label <- gsub("_", " ", species.tree$tip.label)

# Set a new node label for the common ancestral node of the given species
# tree - the tree to update
# new.name - the new name for a node
# ... the tip labels defining the node to update; MRCA of given tips.
update.node.label <- function(tree, new.name, ...){
  test.node <- ape::getMRCA(tree, c(...))
  tree$node.label[test.node-length(tree$tip.label)] <- new.name
  tree
}

# Set node names for species tree
species.tree <- update.node.label(species.tree, "Mus-Arvicanthis", "Mus musculus", "Arvicanthis niloticus")
species.tree <- update.node.label(species.tree, "Rattus", "Rattus norvegicus", "Rattus rattus")
species.tree <- update.node.label(species.tree, "Murinae", "Mus musculus", "Rattus rattus")
species.tree <- update.node.label(species.tree, "Muridae", "Mus musculus", "Meriones unguiculatus")
species.tree <- update.node.label(species.tree, "Cricetidae", "Peromyscus maniculatus", "Phodopus roborovskii")
species.tree <- update.node.label(species.tree, "Eumuroida", "Mus musculus", "Phodopus roborovskii")
species.tree <- update.node.label(species.tree, "Muroidea", "Mus musculus", "Castor canadensis")
species.tree <- update.node.label(species.tree, "Muroidea-Fukomys", "Mus musculus", "Fukomys damarensis")
species.tree <- update.node.label(species.tree, "Xerinae", "Marmota marmota", "Urocitellus parryii")
species.tree <- update.node.label(species.tree, "Sciuridae", "Marmota marmota", "Sciurus carolinensis")
species.tree <- update.node.label(species.tree, "Rodentia", "Mus musculus", "Sciurus carolinensis")
species.tree <- update.node.label(species.tree, "Hominini", "Pan troglodytes", "Homo sapiens")
species.tree <- update.node.label(species.tree, "Homininae", "Gorilla gorilla", "Homo sapiens")
species.tree <- update.node.label(species.tree, "Hominidae", "Pongo pygmaeus", "Homo sapiens")
species.tree <- update.node.label(species.tree, "Cercopithecinae", "Macaca mulatta", "Papio anubis")
species.tree <- update.node.label(species.tree, "Cercopithecidae", "Rhinopithecus roxellana", "Papio anubis")
species.tree <- update.node.label(species.tree, "Catarrhini", "Homo sapiens", "Papio anubis")
species.tree <- update.node.label(species.tree, "Simiiformes", "Callithrix jacchus", "Papio anubis")
species.tree <- update.node.label(species.tree, "Euarchonoglires", "Mus musculus", "Papio anubis")
species.tree <- update.node.label(species.tree, "Otarioidea", "Zalophus californianus", "Odobenus rosmarus")
species.tree <- update.node.label(species.tree, "Otarioidea-Mustela", "Mustela erminea", "Odobenus rosmarus")
species.tree <- update.node.label(species.tree, "Ursus", "Ursus arctos", "Ursus maritimus")
species.tree <- update.node.label(species.tree, "Arctoidea", "Mustela erminea", "Ursus maritimus")
species.tree <- update.node.label(species.tree, "Canidae", "Canis lupus", "Vulpes vulpes")
species.tree <- update.node.label(species.tree, "Caniformia", "Canis lupus", "Odobenus rosmarus")
species.tree <- update.node.label(species.tree, "Carnivora", "Canis lupus", "Felis catus")
species.tree <- update.node.label(species.tree, "Felinae", "Lynx canadensis", "Felis catus")
species.tree <- update.node.label(species.tree, "Ferae", "Manis pentadactyla", "Felis catus")
species.tree <- update.node.label(species.tree, "Panperissodactyla", "Equus caballus", "Felis catus")
species.tree <- update.node.label(species.tree, "Bovidae", "Capra hircus", "Bos taurus")
species.tree <- update.node.label(species.tree, "Pecora", "Odocoileus virginianus", "Bos taurus")
species.tree <- update.node.label(species.tree, "Cetacea", "Balaenoptera musculus", "Monodon monoceros")
species.tree <- update.node.label(species.tree, "Cetruminantia", "Balaenoptera musculus", "Bos taurus")
species.tree <- update.node.label(species.tree, "Artiodactyla", "Sus scrofa", "Bos taurus")
species.tree <- update.node.label(species.tree, "Euungulata", "Sus scrofa", "Zalophus californianus")
species.tree <- update.node.label(species.tree, "Scrotifera", "Phyllostomus discolor", "Bos taurus")
species.tree <- update.node.label(species.tree, "Boreoeutheria", "Mus musculus", "Bos taurus")
species.tree <- update.node.label(species.tree, "Atlantogenata", "Choloepus didactylus", "Loxodonta africana")
species.tree <- update.node.label(species.tree, "Eutheria", "Mus musculus", "Loxodonta africana")
species.tree <- update.node.label(species.tree, "Marsupialia", "Monodelphis domestica", "Phascolarctos cinereus")
species.tree <- update.node.label(species.tree, "Theria", "Mus musculus", "Phascolarctos cinereus")
species.tree <- update.node.label(species.tree, "Monotremata", "Ornithorhynchus anatinus", "Tachyglossus aculeatus")
species.tree <- update.node.label(species.tree, "Mammalia", "Mus musculus", "Tachyglossus aculeatus")
ape::write.tree(species.tree, "aln/node.labeled.species.tree.nwk")

# Now remove time info
species.tree$edge.length <- NULL # remove times

# Make the ZFX tree
zfx.phylogeny <- species.tree
zfx.phylogeny$tip.label <- sapply(zfx.phylogeny$tip.label, \(x) unique(METADATA$mammal[METADATA$mammal$species==x & (METADATA$mammal$group=="ZFX" | METADATA$mammal$group=="Outgroup"),]$common.name), simplify = TRUE)
zfx.phylogeny$tip.label <- gsub(" ", "_", zfx.phylogeny$tip.label)
ape::write.tree(zfx.phylogeny, "aln/zfx_only/zfx.nt.species.nwk")

# For ZFY, we need to account for the two copies in mouse and Nile Rat. Split
# these into two tips each
zfy.phylogeny <- species.tree
mus.node <- which(zfy.phylogeny$tip.label=="Mus musculus")
zfy.phylogeny <- phytools::bind.tip(zfy.phylogeny, "Mus musculus", where=mus.node)
arvicanthis.node <- which(zfy.phylogeny$tip.label=="Arvicanthis niloticus")
zfy.phylogeny <- phytools::bind.tip(zfy.phylogeny, "Arvicanthis niloticus", where=arvicanthis.node)

# Replace node numbers for the newly added sequences
zfy.phylogeny <- update.node.label(zfy.phylogeny, "Mus", "Mus musculus", "Mus musculus")
zfy.phylogeny <- update.node.label(zfy.phylogeny, "Arvicanthis", "Arvicanthis niloticus", "Arvicanthis niloticus")

# Get the ordered set of tip labels from metadata based on the tree order. 
# Search only unique species names, gets two hits for Mus+Arvicanthis, flatten the list result
zfy.phylogeny$tip.label <- unlist( lapply(unique(zfy.phylogeny$tip.label),
                                  \(x) unique(METADATA$mammal[METADATA$mammal$species==x & 
                                                                (METADATA$mammal$group=="ZFY" | 
                                                                   METADATA$mammal$group=="Outgroup"),]$common.name) ))
zfy.phylogeny$tip.label <- gsub(" ", "_", zfy.phylogeny$tip.label)
ape::write.tree(zfy.phylogeny, "aln/zfy_only/zfy.nt.species.nwk")

#### Fetch TimeTree divergence times to highlight the rapid evolution in the rodents ####

get.time.tree <- function(tax.a, tax.b){
  tryCatch({
    tt <- httr::GET(paste0("http://timetree.temple.edu/api/pairwise/",tax.a, "/", tax.b))
    Sys.sleep(0.5) # rate limit
    if(as.numeric(tt$headers$`content-length`)>1000){
      message(paste("Cannot get data for taxa", tax.a, "and", tax.b))
      return(data.frame("taxon_a_id" = tax.a,"taxon_b_id" = tax.b,"scientific_name_a" = NA,
                        "scientific_name_b" = NA,"all_total" = NA,"precomputed_age" = NA,
                        "precomputed_ci_low" = NA,"precomputed_ci_high" = NA,"adjusted_age" = NA))
    }
    
    ct <- content(tt)
    # Parse the response
    data <- as.data.frame( str_split(ct, "\r\n"), col.names = c("V1")) %>% 
      dplyr::slice_tail(n=1) %>%
      tidyr::separate_wider_delim( cols = V1, delim = ",", 
                                   names = c("taxon_a_id","taxon_b_id","scientific_name_a",
                                             "scientific_name_b","all_total","precomputed_age",
                                             "precomputed_ci_low","precomputed_ci_high","adjusted_age"))
    return(data)
  }, error = function(e) { 
    message(paste("Cannot get data for taxa", tax.a, "and", tax.b))
    message("Original error message:")
    message(conditionMessage(e))
    
    return(data.frame("taxon_a_id" = tax.a,"taxon_b_id" = tax.b,"scientific_name_a" = NA,
                      "scientific_name_b" = NA,"all_total" = NA,"precomputed_age" = NA,
                      "precomputed_ci_low" = NA,"precomputed_ci_high" = NA,"adjusted_age" = NA))
  }  ) 
}

# Save to avoid repeated API calls
if(!file.exists( "time.tree.data.tsv")){
  
  # Get the NCBI taxon ids for each species and add to metadata
  taxon.data <- lapply( METADATA$mammal$species, \(x) httr::GET( paste0("http://timetree.temple.edu/api/taxon/",curl::curl_escape(x)))  ) 
  taxon.ids <- sapply(lapply(taxon.data, httr::content), \(x) x$taxon_id)
  METADATA$mammal$taxon_id <- taxon.ids
  
  # Find the pairwise distances between each species
  pairwise.species <- expand.grid(unique(METADATA$mammal$taxon_id), unique(METADATA$mammal$taxon_id)) %>%
    dplyr::filter(Var1!=Var2, Var1<Var2) %>% # only call each pair once
    dplyr::arrange(Var1, Var2)
  
  pairwise.times <- do.call(rbind, mapply(get.time.tree, pairwise.species$Var1, pairwise.species$Var2, SIMPLIFY = FALSE))
  write_tsv(pairwise.times, file = "time.tree.data.tsv", col_names = T, quote = "none")
}

#### Export ZFX and ZFY separately ####

# Write ZFX sequences to file
zfx.nt.aln <- ALIGNMENTS$nt.mammal.biostrings@unmasked[names(ALIGNMENTS$nt.mammal.biostrings@unmasked) %in% c(METADATA$mammal$common.name[METADATA$mammal$group=="ZFX"], 
                                                                                                              METADATA$mammal$common.name[METADATA$mammal$group=="Outgroup"])]
Biostrings::writeXStringSet(zfx.nt.aln,  file = "aln/zfx_only/zfx.aln", format = "fasta")

# Write ZFY sequences to file
zfy.nt.aln <- ALIGNMENTS$nt.mammal.biostrings@unmasked[names(ALIGNMENTS$nt.mammal.biostrings@unmasked) %in% c(METADATA$mammal$common.name[METADATA$mammal$group=="ZFY"], 
                                                                                                              METADATA$mammal$common.name[METADATA$mammal$group=="Outgroup"])]
Biostrings::writeXStringSet(zfy.nt.aln,  file = "aln/zfy_only/zfy.aln", format = "fasta")

#### Create independent trees for ZFX and ZFY sequences with species trees ####

cat(timestamp(), "Running ancestral sequence reconstruction with species tree\n")

# Run the ancestral reconstructions
run.iqtree("aln/zfx_only/zfx.aln", 
           "-nt AUTO", # number of threads
           "-te aln/zfx_only/zfx.nt.species.nwk", # user tree guide
           "-asr") # ancestral sequence reconstruction

run.iqtree("aln/zfy_only/zfy.aln", 
           "-nt AUTO", # number of threads
           "-te aln/zfy_only/zfy.nt.species.nwk", # user tree guide
           "-asr") # ancestral sequence reconstruction

#### Remove duplicate species nodes from the trees  ####

# For a ~species tree, we only need one sequence per species. Either can be dropped,
# since they give the same branching order

# Read the ML ZFX and ZFY trees 
zfx.nt.aln.tree <- ape::read.tree("aln/zfx_only/zfx.aln.treefile")
zfy.nt.aln.tree <- ape::read.tree("aln/zfy_only/zfy.aln.treefile")

# Drop the second ZFYs in mouse and rat
zfy.nt.aln.tree <- tidytree::drop.tip(zfy.nt.aln.tree, "Mouse_Zfy2") 
zfy.nt.aln.tree <- tidytree::drop.tip(zfy.nt.aln.tree, "African_Grass_Rat_ZFY2-like_1") 

# Root the trees on monotremes
zfx.nt.aln.tree <- reroot.tree(zfx.nt.aln.tree, c("Platypus_ZFX", "Autralian_echidna_ZFX"), position = 0.015)
zfy.nt.aln.tree <- reroot.tree(zfy.nt.aln.tree, c("Platypus_ZFX", "Autralian_echidna_ZFX"), position = 0.015)

# Remove gene names so tip labels are comparable
zfx.nt.aln.tree$tip.label <- str_replace(zfx.nt.aln.tree$tip.label, "_Z[F|f][X|x].*", "")
zfy.nt.aln.tree$tip.label <- str_replace(zfy.nt.aln.tree$tip.label, "(_putative)?(_|-)Z[F|f][X|x|Y|y].*", "")

#### Plot ZFX / ZFY tree comparisons  ####

# Plot the two trees with node labels
zfx.nt.aln.tree.plot <- plot.tree(zfx.nt.aln.tree) + geom_nodelab(size=2, nudge_x = -0.003, nudge_y = 0.5, hjust=1,  node = "internal")
zfy.nt.aln.tree.plot <- plot.tree(zfy.nt.aln.tree) + geom_nodelab(size=2, nudge_x = -0.003, nudge_y = 0.5, hjust=1,  node = "internal")

save.double.width("figure/species.tree.zfx.png", zfx.nt.aln.tree.plot)
save.double.width("figure/species.tree.zfy.png", zfy.nt.aln.tree.plot)

#### Export ancestral ZFX/Y reconstructions for each node in the species tree ####

# # Read the ancestral states for the nodes
ancestral.zfx.seqs <- read.table("aln/zfx_only/zfx.aln.state", header=TRUE)
ancestral.zfy.seqs <- read.table("aln/zfy_only/zfy.aln.state", header=TRUE)

get.node.sequence <- function(node.name, ancestral.seq.table, type){
  ancestral.seq.table %>%
    dplyr::filter(Node==node.name) %>%
    dplyr::arrange(Site) %>%
    dplyr::select(State) %>%
    dplyr::summarise(Seq =  paste0(">", type, "_", node.name, "\n", paste(State, collapse = "")))
}

# 
# # Get the ancestral sequences for these node pairs
anc.zfx <- sapply(species.tree$node.label[-1], # skip Mammalia
                                       get.node.sequence,
                                       ancestral.seq.table=ancestral.zfx.seqs, type="ZFX")
anc.zfy <- sapply(species.tree$node.label[-1], # skip Mammalia
                                       get.node.sequence,
                                       ancestral.seq.table=ancestral.zfy.seqs, type="ZFY")


# Write the node sequences to fasta
anc.node.seqs <- paste0(c(anc.zfx, anc.zfy, "\n"), collapse = "\n")
write_file(anc.node.seqs, "aln/ancestral.zfx.zfy.nodes.fa")

# Combine with existing ZFX/Y sequence file
ape::write.dna(ALIGNMENTS$nt.mammal.ape, "aln/ancestral.zfx.zfy.nodes.fa", format="fasta", 
               append = TRUE, colsep = "", nbcol=-1)

#### Write geneconv configuration file ####

# Put Zfx and Zfy for each ancestral node in a group string
ancestral.node.string <- paste( paste0("-group ", species.tree$node.label[-1], " ZFX_", species.tree$node.label[-1], " ZFY_", species.tree$node.label[-1]), collapse = "\n")

# Now do the same for the species sequences
split.gene.names <- METADATA$mammal %>%
  # Since we're about to split on 'Z', remove secondary instances of the string
  # keeping whatever - or _ was in the original name
  dplyr::mutate(common.name = str_replace(common.name, "(?i)putative-Zfy", "putative-Y"),
                common.name = str_replace(common.name, "(?i)putative_Zfy", "putative_Y")) %>%
  separate_wider_delim(common.name, 
                       delim= "Z", names = c("Species", "Gene"), too_few = "debug", too_many = "debug") %>%
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


species.string <- paste(split.gene.names$GroupString, collapse = "\n")

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
                                species.string, "\n", ancestral.node.string, "\n")

write_file(anc.gene.conv.control, "aln/anc.zfx.zfy.geneconv.cfg")

#### Read final intron sequences #####
cat(timestamp(), "Processing final intron sequences\n")
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

#### Align final intron sequences #####

# Use MUSCLE to align  - non coding
cat(timestamp(), "Aligning final inton sequences\n")
run.muscle(FILES$final.intron.zfy.nt.fas, FILES$final.intron.zfy.nt.aln)
run.muscle(FILES$final.intron.zfx.nt.fas, FILES$final.intron.zfx.nt.aln)
run.muscle(FILES$final.intron.nt.fas,     FILES$final.intron.nt.aln)

#### Run divvier to improve high-confidence homologies in final intron alignments ####

final.intron.zfy.nt.divvy.aln <- run.divvier(FILES$final.intron.zfy.nt.aln)
final.intron.zfx.nt.divvy.aln <- run.divvier(FILES$final.intron.zfx.nt.aln)

#### Make species phylogenies for Zfx and Zfy final intron  #####
cat(timestamp(), "Creating final intron species trees\n")

# We want to specify actual species tree for each of Zfx and Zfy
# Note that not all species have intron info available - make a new species tree
zfx.phylogeny <- paste0("", # Outgroups
                        "(", # Eutheria
                        "(Southern_two-toed_sloth_ZFX, African_bush_elephant_ZFX)Atlantogenata, ", # Atlantogenata  
                        "(", # Boreoeutheria
                        "(",  # Euarchonoglires
                        "( ", # Simiiformes
                        "Common_marmoset_ZFX,", # New world monkeys
                        "(", # Catarrhini (Old world monkeys & apes)
                        "(", #Cercopithecidae (Old world monkeys)
                        "Golden_snub-nosed_monkey_ZFX,",   #Colobinae
                        "(Olive_baboon_ZFX, Macaque_ZFX)Cercopithecinae", # Cercopithecinae
                        ")Cercopithecidae,", # /Cercopithecidae (Old world monkeys)
                        "(", # Hominidae
                        "Gorilla_ZFX, (Chimpanzee_ZFX, Human_ZFX)Hominini",
                        ")Hominidae",  # /Hominidae
                        ")Catarrhini", # /Catarrhini
                        ")Simiiformes,",  # /Simiiformes
                        "(", # Rodentia
                        "(Gray_squirrel_Zfx,(Arctic_ground_squirrel_Zfx, Alpine_marmot_ZFX)Xerinae)Sciuridae,", # Sciuridae 
                        "(Damara_mole-rat_Zfx,", # Muroidea-Fukomys
                        "(Beaver_Zfx, (", # Muroidea
                        "(North_American_deer_mouse_Zfx, Desert_hamster_Zfx)Cricetidae,", # Cricetidae 
                        "(Mongolian_gerbil_Zfx, (Rat_Zfx, (Mouse_Zfx, African_Grass_Rat_Zfx)Mus-Arvicanthis)Murinae)Muridae", # Muridae 
                        ")Eumuroida)Muroidea",  # /Muroidea
                        ")Muroidea-Fukomys", # /Muroidea-Fukomys
                        ")Rodentia", # /Rodentia
                        ")Euarchonoglires,", # /Euarchonoglires
                        "(", # Laurasiatheria
                        "(",  # Carnivora
                        "Cat_ZFX, (Dog_ZFX, (Stoat_ZFX, Polar_bear_ZFX)Arctoidea)Caniformia",
                        ")Carnivora,",  # /Carnivora
                        "(",  # Euungulata 
                        "Horse_ZFX, (Pig_ZFX, (White_tailed_deer_ZFX, (Cattle_ZFX, Goat_ZFX)Bovidae)Pecora)Artiodactyla",
                        ")Euungulata", # /Euungulata 
                        ")Laurasiatheria", # /Laurasiatheria
                        ")Boreoeutheria", # /Boreoeutheria
                        ")Eutheria;" # /Eutheria
) # /Outgroups

write_file(zfx.phylogeny, "aln/final.intron.zfx/zfx.nt.species.nwk")

zfy.phylogeny <- paste0("", # Outgroups
                        "(", # Eutheria
                        "(Southern_two-toed_sloth_ZFY, African_bush_elephant_ZFY)Atlantogenata, ", # Afrotheria & Xenarthra
                        "(", # Boreoeutheria
                        "(",  # Euarchonoglires
                        "( ", # Simiiformes
                        "Common_marmoset_ZFY,", # New world monkeys
                        "(", # Catarrhini (Old world monkeys & apes)
                        "(", #Cercopithecidae (Old world monkeys)
                        "Golden_snub-nosed_monkey_ZFY,",   #Colobinae
                        "(Olive_baboon_ZFY, Macaque_ZFY)Cercopithecinae", # Cercopithecinae
                        ")Cercopithecidae,", # /Cercopithecidae (Old world monkeys)
                        "(", # Hominidae
                        "Gorilla_ZFY, (Chimpanzee_ZFY, Human_ZFY)Hominini",
                        ")Hominidae",  # /Hominidae
                        ")Catarrhini", # /Catarrhini
                        ")Simiiformes,",  # /Simiiformes
                        "(", # Rodentia
                        "(Gray_squirrel_Zfy,(Arctic_ground_squirrel_Zfx-like_putative-Zfy, Alpine_marmot_ZFY)Xerinae)Sciuridae,", # Sciuridae 
                        "(Damara_mole-rat_Zfy,", 
                        "(Beaver_Zfx-like_putative-Zfy, (", # Muroidea
                        "(North_American_deer_mouse_Zfx-like_putative-Zfy, Desert_hamster_Zfx-like_putative-Zfy)Cricetidae,", # Cricetidae 
                        "(Mongolian_gerbil_Zfx-like_putative-Zfy, (Rat_Zfy2, ((Mouse_Zfy1, Mouse_Zfy2), (African_Grass_Rat_ZFY2-like_1, African_Grass_Rat_ZFY2-like_2))Mus-Arvicanthis)Murinae)Muridae", # Muridae 
                        ")Eumuroida)Muroidea)Muroidea-Fukomys", # /Muroidea
                        ")Rodentia", # /Rodentia
                        ")Euarchonoglires,", # /Euarchonoglires
                        "(", # Laurasiatheria
                        "(",  # Carnivora
                        "Dog_ZFY, (Stoat_ZFY, Polar_bear_ZFY)Arctoidea",
                        ")Carnivora,",  # /Carnivora
                        "(",  # Euungulata 
                        "Horse_ZFY, (Pig_ZFY, (White_tailed_deer_ZFY, (Cattle_ZFY, Goat_ZFY)Bovidae)Pecora)Artiodactyla",
                        ")Euungulata", # /Euungulata 
                        ")Laurasiatheria", # /Laurasiatheria
                        ")Boreoeutheria", # /Boreoeutheria
                        ")Eutheria;" # /Eutheria
) # /Outgroups

write_file(zfy.phylogeny, "aln/final.intron.zfy/zfy.nt.species.nwk")

#### Make ML tree based on final intron raw alignments #####

cat(timestamp(), "Creating final intron ML trees using species phylogenies\n")
# Make the tree using the species phylogeny
final.intron.zfy.nt.aln.treefile <- run.iqtree(FILES$final.intron.zfy.nt.aln, 
                                               "-nt AUTO", # number of threads
                                               "-keep-ident",
                                               "-te aln/final.intron.zfy/zfy.nt.species.nwk" # user tree guide
) 

final.intron.zfx.nt.aln.treefile <- run.iqtree(FILES$final.intron.zfx.nt.aln, 
                                               "-nt AUTO", # number of threads
                                               "-keep-ident",
                                               "-te aln/final.intron.zfx/zfx.nt.species.nwk" # user tree guide
) 

#### Make ML tree based on final intron divvied alignments ####
cat(timestamp(), "Creating divvied final intron ML trees using species phylogenies\n")

final.intron.zfy.nt.divvy.aln.treefile <- run.iqtree(final.intron.zfy.nt.divvy.aln, 
                                                     "-nt AUTO", # number of threads
                                                     "-keep-ident",
                                                     "-te aln/final.intron.zfy/zfy.nt.species.nwk" # user tree guide
)

final.intron.zfx.nt.divvy.aln.treefile <- run.iqtree(final.intron.zfx.nt.divvy.aln, 
                                                     "-nt AUTO", # number of threads
                                                     "-keep-ident",
                                                     "-te aln/final.intron.zfx/zfx.nt.species.nwk" # user tree guide
)

#### Prepare codeml site model to check for site-specific selection ####

cat(timestamp(), "Creating CODEML site model control files\n")
# Adapted from Beginner's Guide on the Use of PAML to Detect Positive Selection
# https://academic.oup.com/mbe/article/40/4/msad041/7140562 for details

# Two files need to be uploaded to the cluster for running codeml:
# - alignment
# - tree file (unrooted)

# The treefile created earlier needs node names and branch lengths removing.
# Set node names to the empty string, and set branch lengths to 1
labelled.tree <- ape::read.tree(paste0(FILES$mammal.nt.aln, ".treefile"))
labelled.tree$node.label <- rep("", length(labelled.tree$node.label))
labelled.tree$edge.length <- rep(1, length(labelled.tree$edge.length))
ape::write.tree(labelled.tree, file = "paml/site-specific/zfxy.nt.aln.paml.treefile")

# We now read the tree as a raw string to remove the branch lengths entirely,
# separate the labels from taxon names and resave
labelled.tree <- readr::read_file("paml/site-specific/zfxy.nt.aln.paml.treefile")
labelled.tree <- gsub(":1", "", labelled.tree)
readr::write_file(labelled.tree, "paml/site-specific/zfxy.nt.aln.paml.treefile")

# This control file tests site models With heterogeneous ω Across Sites
paml.site.file <- paste0("seqfile   = ../../", FILES$mammal.nt.aln, " * alignment file\n",
                         "treefile  = zfxy.nt.aln.paml.treefile * tree in Newick format without nodes\n", # 
                         "outfile   = site.specific.paml.out.txt\n",
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

#### Prepare codeml branch-site model to look for selection specifically in Eumuroida ####

cat(timestamp(), "Creating CODEML branch-site model control files\n")
# To look at the rodent clade, we need a rooted tree
nt.aln.tree <- ape::read.tree(FILES$mammal.nt.aln.treefile)
# Root the tree on platypus and resave
# The root is arbitrarily placed in the platypus branch to fit neatly
nt.aln.tree <- reroot.tree(nt.aln.tree, c("Platypus_ZFX", "Australian_echidna_ZFX"), position = 0.015)
ape::write.tree(nt.aln.tree, file = paste0(FILES$mammal.nt.aln, ".rooted.treefile"))

# Find the nodes that are ZFY vs ZFX and add to tree
mammal.gene.groups <- split(METADATA$mammal$common.name, METADATA$mammal$group)
nt.aln.tree <- tidytree::groupOTU(nt.aln.tree, mammal.gene.groups, group_name = "group")



# Label the given node and its children and save to the given files
# Creates null model and test model directories
write.paml.fg.tree <- function(fg.node, dir.base){
  filesstrings::create_dir(dir.base)
  filesstrings::create_dir(paste0(dir.base, "-null"))
  # Remove existing node labels, label the nodes and tips of the tree with #1
  # for foreground branches
  labelled.tree <- nt.aln.tree
  labelled.tree$node.label <- rep("", length(labelled.tree$node.label))
  labelled.tree$edge.length <- rep(1, length(labelled.tree$edge.length))
  labelled.tree <- treeio::label_branch_paml(labelled.tree, fg.node, "#1")
  
  tree.file <- paste0(dir.base, "/tree.treefile")
  tree.file.null <- paste0(dir.base, "-null/tree.treefile")
  fg.tree.file <- paste0(dir.base, "/tree.fg.treefile")
  fg.tree.file.null <- paste0(dir.base, "-null/tree.fg.treefile")
  control.file <- paste0(dir.base, "/paml.ctl")
  control.file.null <- paste0(dir.base, "-null/paml.ctl")
  
  ape::write.tree(labelled.tree, file = tree.file)
  ape::write.tree(labelled.tree, file = tree.file.null)
  
  # Then remove the branch lengths, separate the labels from taxon names and resave
  newick.test <- read_file(tree.file)
  newick.test <- gsub(":1", "", newick.test) # branch lengths we set to 1
  newick.test <- gsub("_#1", " #1", newick.test) # fg labels
  write_file(newick.test, fg.tree.file)
  write_file(newick.test, fg.tree.file.null)
  
  # Create control files for the test and null codeml
  # This control file tests site models With heterogeneous ω across Sites
  control.file.data <- paste0("seqfile   = ../../", FILES$mammal.nt.aln, " * alignment file\n",
                              "treefile  = tree.fg.treefile * tree in Newick format without nodes\n", # 
                              "outfile   = paml.out.txt\n",
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
  write_file(control.file.data, control.file)
  
  # Prepare codeml branch-site null model to look for selection specifically in Muroidea
  # This control file tests null branch site model With fixed ω across sites
  control.file.data.null <- paste0("seqfile   = ../../", FILES$mammal.nt.aln, " * alignment file\n",
                                   "treefile  = tree.fg.treefile * tree in Newick format without nodes\n", # 
                                   "outfile   = paml.out.txt\n",
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
                                   "fix_omega = 1 * The value of ω2 for foreground branches will be fixed\n",
                                   "omega     = 1 * The fixed value will be ω2 = 1\n")
  write_file(control.file.data.null, control.file.null)
  
}

rodentia.node <- ape::getMRCA(nt.aln.tree, c("Mouse_Zfy2", "Gray_squirrel_Zfy"))
write.paml.fg.tree(rodentia.node, 
                   "paml/branch-site-rodentia")

# Find the MRCA of the rodents with rapid ZFY evolution
eumuroida.node <- ape::getMRCA(nt.aln.tree, c("Mouse_Zfy2", "Desert_hamster_Zfx-like_putative-Zfy"))
# Create labelled tree for rodent node
write.paml.fg.tree(eumuroida.node, 
                   "paml/branch-site-eumuroida")

cat(timestamp(), "Creating CODEML control script\n")
# Add a script to submit these jobs to the cluster
# Run the commands manually - this takes a long time to run, so only invoke when needed
paml.shell.script <- paste0("#!/bin/bash\n\n",
                            "# qsubme is an alias to submit the job to the cluster\n",
                            "shopt -s expand_aliases\n",
                            "source ~/.bashrc\n\n",
                            "cd paml/site-specific\n",
                            "# codeml zfy.site-specific.paml.ctl\n",
                            "qsubme codeml zfy.site-specific.paml.ctl\n\n",

                            "cd ../branch-site-eumuroida\n",
                            "qsubme codeml paml.ctl\n",
                            "cd ../branch-site-eumuroida-null\n",
                            "qsubme codeml paml.ctl\n\n",

)
write_file(paml.shell.script, "run_paml.sh")

#### Run HyPhy RELAX to test for relaxed selection and MEME for diversifying selection ####

cat(timestamp(), "Creating HyPhy RELAX and MEME control files\n")
# Create and save a tree in HyPhy format. Set test branches to the given node, all others
# as reference. 
create.mammal.hyphy.relax.tree.file <- function(fg.node, node.name){
  # Newick tree tips and nodes need to be tagged with {Test} and {Reference}
  # Set all branches and nodes to reference
  hyphy.tree <- ape::read.tree(FILES$mammal.nt.aln.treefile)
  
  # Root the tree on platypus
  # The root is arbitrarily placed in the platypus branch to fit neatly
  hyphy.tree <- reroot.tree(hyphy.tree, c("Platypus_ZFX", "Australian_echidna_ZFX"), position = 0.015)
  
  # Find the nodes that are ZFY vs ZFX and add to tree
  mammal.gene.groups <- split(METADATA$mammal$common.name, METADATA$mammal$group)
  hyphy.tree <- tidytree::groupOTU(hyphy.tree, mammal.gene.groups, group_name = "group")
  
  hyphy.tree$node.label <- rep("", length(hyphy.tree$node.label)) # remove all node names
  hyphy.tree$node.label <- paste0(hyphy.tree$node.label, "{Reference}") # everything is reference
  hyphy.tree$tip.label <- paste0(hyphy.tree$tip.label, "{Reference}") # everything is reference
  hyphy.tree <- treeio::label_branch_paml(hyphy.tree, fg.node, "{Test}") # now rodents are test
  hyphy.tree$node.label <- gsub("\\{Reference\\} \\{Test\\}", "{Test}", hyphy.tree$node.label) # remove reference from anything in test
  hyphy.tree$tip.label <- gsub("\\{Reference\\} \\{Test\\}", "{Test}", hyphy.tree$tip.label) # remove reference from anything in test
  
  tree.file <- paste0("aln/hyphy/mammal.", node.name, ".hyphy.relax.treefile")
  
  ape::write.tree(hyphy.tree, file = tree.file)
}

create.combined.hyphy.relax.tree.file <- function(fg.node, node.name){
  # Newick tree tips and nodes need to be tagged with {Test} and {Reference}
  # Set all branches and nodes to reference
  hyphy.tree <- ape::read.tree(FILES$combined.nt.aln.treefile)
  
  # Root on Xenopus
  xenopus.node <- ape::getMRCA(hyphy.tree, c("Xenopus_ZFX.S","Xenopus_ZFX.L"))
  hyphy.tree <- phytools::reroot(hyphy.tree, xenopus.node, position = 0.01)
  
  # Find the nodes that are ZFY vs ZFX and add to tree
  gene.groups <- split(METADATA$combined$common.name, METADATA$combined$group)
  hyphy.tree <- tidytree::groupOTU(hyphy.tree, gene.groups, group_name = "group")
  
  hyphy.tree$node.label <- rep("", length(hyphy.tree$node.label)) # remove all node names
  hyphy.tree$node.label <- paste0(hyphy.tree$node.label, "{Reference}") # everything is reference
  hyphy.tree$tip.label <- paste0(hyphy.tree$tip.label, "{Reference}") # everything is reference
  hyphy.tree <- treeio::label_branch_paml(hyphy.tree, fg.node, "{Test}") # now rodents are test
  hyphy.tree$node.label <- gsub("\\{Reference\\} \\{Test\\}", "{Test}", hyphy.tree$node.label) # remove reference from anything in test
  hyphy.tree$tip.label <- gsub("\\{Reference\\} \\{Test\\}", "{Test}", hyphy.tree$tip.label) # remove reference from anything in test
  
  tree.file <- paste0("aln/hyphy/combined.", node.name, ".hyphy.relax.treefile")
  
  ape::write.tree(hyphy.tree, file = tree.file)
}

# MEME cannot handle either / or . in node names. Replace both
create.mammal.hyphy.meme.tree.file <- function(fg.node, node.name){
  hyphy.tree <- ape::read.tree(FILES$mammal.nt.aln.treefile)
  
  # Root the tree on platypus
  # The root is arbitrarily placed in the platypus branch to fit neatly
  hyphy.tree <- reroot.tree(hyphy.tree, c("Platypus_ZFX", "Australian_echidna_ZFX"), position = 0.015)
  
  # Find the nodes that are ZFY vs ZFX and add to tree
  mammal.gene.groups <- split(METADATA$mammal$common.name, METADATA$mammal$group)
  hyphy.tree <- tidytree::groupOTU(hyphy.tree, mammal.gene.groups, group_name = "group")
  
  hyphy.tree$node.label <- gsub("\\.", "_", hyphy.tree$node.label)
  hyphy.tree$node.label <- gsub("/", "_", hyphy.tree$node.label)
  
  hyphy.tree$node.label <- paste0(hyphy.tree$node.label, "{Reference}") # everything is reference
  hyphy.tree$tip.label <- paste0(hyphy.tree$tip.label, "{Reference}") # everything is reference
  hyphy.tree <- treeio::label_branch_paml(hyphy.tree, fg.node, "{Test}") # now rodents are test
  hyphy.tree$node.label <- gsub("\\{Reference\\} \\{Test\\}", "{Test}", hyphy.tree$node.label) # remove reference from anything in test
  hyphy.tree$tip.label <- gsub("\\{Reference\\} \\{Test\\}", "{Test}", hyphy.tree$tip.label) # remove reference from anything in test
  
  tree.file <- paste0("aln/hyphy/mammal.", node.name,".meme.hyphy.treefile")
  ape::write.tree(hyphy.tree, file = tree.file)
}

create.combined.hyphy.meme.tree.file <- function(fg.node, node.name){
  hyphy.tree <- ape::read.tree(FILES$combined.nt.aln.treefile)
  
  # Root on Xenopus
  hyphy.tree <- reroot.tree(hyphy.tree, c("Xenopus_ZFX.S","Xenopus_ZFX.L"), position = 0.01)
  
  # Find the nodes that are ZFY vs ZFX and add to tree
  gene.groups <- split(METADATA$combined$common.name, METADATA$combined$group)
  hyphy.tree <- tidytree::groupOTU(hyphy.tree, gene.groups, group_name = "group")
  
  hyphy.tree$node.label <- gsub("\\.", "_", hyphy.tree$node.label)
  hyphy.tree$node.label <- gsub("/", "_", hyphy.tree$node.label)
  
  hyphy.tree$node.label <- paste0(hyphy.tree$node.label, "{Reference}") # everything is reference
  hyphy.tree$tip.label <- paste0(hyphy.tree$tip.label, "{Reference}") # everything is reference
  hyphy.tree <- treeio::label_branch_paml(hyphy.tree, fg.node, "{Test}") # now rodents are test
  hyphy.tree$node.label <- gsub("\\{Reference\\} \\{Test\\}", "{Test}", hyphy.tree$node.label) # remove reference from anything in test
  hyphy.tree$tip.label <- gsub("\\{Reference\\} \\{Test\\}", "{Test}", hyphy.tree$tip.label) # remove reference from anything in test
  
  tree.file <- paste0("aln/hyphy/combined.", node.name,".meme.hyphy.treefile")
  ape::write.tree(hyphy.tree, file = tree.file)
}

nodes <- c(rodentia.node, eumuroida.node, muridae.node, murinae.node)
node.names <- c("rodentia", "eumuroida") #, "muridae", "murinae"

mapply(create.combined.hyphy.relax.tree.file, nodes, node.names)
mapply(create.mammal.hyphy.relax.tree.file,   nodes, node.names)
mapply(create.combined.hyphy.meme.tree.file,  nodes, node.names)
mapply(create.mammal.hyphy.meme.tree.file,    nodes, node.names)

# Create Hyphy RELAX invokations for shell scripting for the given node names
create.hyphy.relax.invokation <- function(node.name){
  paste0(
    "# Mammal alignment for ", node.name, "\n",
    "hyphy relax --alignment  ", FILES$mammal.nt.aln, " --tree aln/hyphy/mammal.", node.name , ".hyphy.relax.treefile --reference 'Reference' --test 'Test' --output aln/hyphy/mammal.", node.name, ".relax.json\n",
    "# Combined alignment for ", node.name, "\n",
    "hyphy relax --alignment  ", FILES$combined.nt.aln, " --tree aln/hyphy/combined.", node.name , ".hyphy.relax.treefile --reference 'Reference' --test 'Test' --output aln/hyphy/combined.", node.name, ".relax.json\n"
  )
}

# Create Hyphy MEME invokations for shell scripting for the given node names
# Use a partition model to allow three different classes of evolution:
# Exons 1, 3-6, exon 2 and exon 7
create.hyphy.meme.invokation <- function(node.name){
  paste0(
    "# Mammal partitioned alignment for ", node.name, "\n",
    "hyphy meme --alignment  ", FILES$mammal.nt.nexus, " --tree aln/hyphy/mammal.", node.name , ".meme.hyphy.treefile --branches Test --output aln/hyphy/mammal.", node.name, ".meme.json\n",
    "# Combined partitioned alignment for ", node.name, "\n",
    "hyphy meme --alignment  ", FILES$combined.nt.nexus, " --tree aln/hyphy/combined.", node.name , ".meme.hyphy.treefile --branches Test --output aln/hyphy/combined.", node.name, ".meme.json\n")
}



# HyPhy is installed in a conda environment (hence we can't use 'system2' to
# directly call 'hyphy'). This creates a script to activate the conda
# environment and run RELAX.
write_file(paste0("#!/bin/bash\n",
                  "source activate hyphy\n\n",
                  "# RELAX test for realxed purifying selection in test branches\n",
                  paste(sapply(node.names, create.hyphy.relax.invokation), 
                        collapse = ""), "\n",
                  "# MEME test for positive selection at individual sites in all branches\n",
                  paste(sapply(node.names, create.hyphy.meme.invokation), 
                        collapse = "")
),
"run_hyphy.sh")

cat(timestamp(), "Running HyPhy RELAX and MEME across nodes", paste(node.names, collapse=";"), "\n")
system2("bash", "run_hyphy.sh")

#### Tar the outputs ####
system2("tar", "czf aln.tar.gz aln")
cat(timestamp(), "Done!\n")