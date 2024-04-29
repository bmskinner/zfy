# Analysis pipeline for ZFX/ZFY evolutionary analysis

# The following binaries must be supplied in ./bin:
# macse_v2.07.jar
# muscle5.1.win64.exe / muscle5.1.linux_intelx64
# nlstradamus.pl
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
#### Imports #####

source("src/functions.R")
load.packages()

source("src/find9aaTADs.R")
source("src/findZF.R")
source("src/calcCharge.R")
source("src/calcHydrophobicity.R")

cat("Packages loaded\n")

#### Create output directory structure #####

# Clear previous analyses
filesstrings::dir.remove("aln")
filesstrings::dir.remove("figure")
filesstrings::dir.remove("nls")
filesstrings::dir.remove("pwm")

# Create missing dirs if needed
filesstrings::create_dir("aln")
filesstrings::create_dir("aln/mammal")
filesstrings::create_dir("aln/combined")
filesstrings::create_dir("aln/exons")
filesstrings::create_dir("aln/hyphy")
filesstrings::create_dir("aln/nls")
filesstrings::create_dir("aln/pwm")
filesstrings::create_dir("aln/zfx_only")
filesstrings::create_dir("aln/zfy_only")
filesstrings::create_dir("bin")
filesstrings::create_dir("figure")
filesstrings::create_dir("paml")
filesstrings::create_dir("paml/site-specific")
filesstrings::create_dir("paml/branch-site")
filesstrings::create_dir("paml/exon_1_3-6")
filesstrings::create_dir("paml/exon_2")
filesstrings::create_dir("paml/exon_7")

writeLines(capture.output(sessionInfo()), "figure/session_info.txt")

#### Read NT FA files #####

prepare.fas.files()

#### Run combined mammal/outgroup AA alignment ####

# Use muscle to align the aa files and save out for iqtree
# Cluster has muscle 3.8.31 on path., so specify the binary directly
muscle.path <- ifelse(installr::is.windows(), "bin/muscle5.1.win64.exe", "bin/muscle5.1.linux_intel64")
if(!file.exists(muscle.path)) stop("Muscle5.1 not present in bin directory")
system2(muscle.path, paste("-align",  files$combined.aa.fas,
                           "-output", files$combined.aa.aln))

#### Run mammal NT alignment guided by AA #####

# Run a codon aware alignment with MACSE
# Expect java on the PATH. Macse download from https://www.agap-ge2pop.org/macsee-pipelines/

# Direct download link:
# https://www.agap-ge2pop.org/wp-content/uploads/macse/releases/macse_v2.07.jar
if(!file.exists("bin/macse_v2.07.jar")) stop("MACSE not present in bin directory")

# Run the alignment
system2("java", paste("-jar bin/macse_v2.07.jar -prog alignSequences",
                      "-seq",    files$mammal.nt.fas, # input
                      "-out_NT", files$mammal.nt.aln,  # output nt alignment
                      "-out_AA", files$mammal.aa.aln), # output aa alignment
        stdout = paste0(files$mammal.nt.aln, ".macse.log"),  # logs
        stderr = paste0(files$mammal.nt.aln, ".macse.err"))  # error logs

# Read all alignments in ape and Biostrings formats
alignments <- read.alignments()

# Identify the coordinates of the exon boundaries in the gapped alignments
# Based on the mouse Zfy1 sequence
mouse.exons <- find.exons(alignments$nt.mammal.biostrings, alignments$aa.combined.biostrings)

#### Create combined mammal/outgroup AA tree ####

system2("iqtree", paste("-s ", files$combined.aa.aln, 
                        "-bb 1000", # number of bootstrap replicates
                        "-alrt 1000", # number of replicates to perform SH-like approximate likelihood ratio test (SH-aLRT) 
                        "-nt AUTO"), # number of threads
        stdout = gsub(".aln$", ".iqtree.log", files$combined.aa.aln), 
        stderr = gsub(".aln$", ".iqtree.log", files$combined.aa.aln))

#### Create mammal CDS NT tree #####

# Make ML tree and reconstruct ancestral sequences
# Expect iqtree on the PATH.
# Note model testing is automatically performed in v1.5.4 onwards
# Note: we can use a partition model if we specify exon coordinates
system2("iqtree", paste("-s ", files$mammal.nt.aln, 
                        "-bb 1000", # number of bootstrap replicates
                        "-alrt 1000", # number of replicates to perform SH-like approximate likelihood ratio test (SH-aLRT) 
                        "-nt AUTO", # number of threads
                        "-asr"), # ancestral sequence reconstruction
        stdout = gsub(".aln$", ".iqtree.log", files$mammal.nt.aln), 
        stderr = gsub(".aln$", ".iqtree.log", files$mammal.nt.aln))

#### Make mammal exon NT trees ####

# Aim here is to look for signs of gene conversion in the final exon as per other 
# papers.

create.exon.alignment <- function(i){
  exon.aln <- as.matrix(alignments$nt.mammal.ape)[,mouse.exons$start[i]:mouse.exons$end[i]]
  exon.aln <- ape::del.colgapsonly(exon.aln, threshold = 0.2) # remove columns with >20% gaps
  exon.aln.file <- paste0("aln/exons/exon_", mouse.exons$exon[i], ".aln")
  
  ape::write.FASTA(exon.aln, file = exon.aln.file)
  
  system2("iqtree", paste("-s ", exon.aln.file,
                          "-bb 1000", # number of bootstrap replicates
                          "-alrt 1000", # number of replicates to perform SH-like approximate likelihood ratio test (SH-aLRT)
                          "-nt AUTO"), # number of threads
          stdout = paste0(exon.aln.file, ".iqtree.log"),
          stderr = paste0(exon.aln.file, ".iqtree.log"))
}

exon.1.7.plots <- lapply(1:nrow(mouse.exons),create.exon.alignment)

# We also want to look at all except exon 7
exon.1.6.aln <- as.matrix(alignments$nt.mammal.ape)[,mouse.exons$start[1]:mouse.exons$end[6]]
exon.1.6.aln.file <- paste0("aln/exons/exon_1-6.aln")
ape::write.FASTA(exon.1.6.aln, file = exon.1.6.aln.file)

system2("iqtree", paste("-s ", exon.1.6.aln.file,
                        "-bb 1000", # number of bootstrap replicates
                        "-alrt 1000", # number of replicates to perform SH-like approximate likelihood ratio test (SH-aLRT)
                        "-nt AUTO"), # number of threads
        stdout = paste0(exon.1.6.aln.file, ".iqtree.log"),
        stderr = paste0(exon.1.6.aln.file, ".iqtree.log"))

#### Identify binding motifs of the ZFs in each species ####
old.wd <- getwd()
setwd("./bin/pwm_predict")
system2("./pwm_predict", "-l 20 ../../fasta/combined.aa.fas") # ensure all ZFs linked
setwd(old.wd)
filesstrings::move_files(files = c("fasta/combined.aa.pwm"),
                         destinations = c("aln/pwm"),
                         overwrite = TRUE)
# Remove header and split the output PWMs to separate files
system2("cat", "aln/pwm/combined.aa.pwm | grep -v ';' | split -l 5 - aln/pwm/zf_")
#### Fetch divergence times to highlight the rapid evolution in the rodents ####

# Get the NCBI taxon ids for each species
taxon.data <- lapply(metadata.mammal$species, \(x) httr::GET( paste0("http://timetree.temple.edu/api/taxon/",curl::curl_escape(x)))  ) 
taxon.ids <- sapply(lapply(taxon.data, httr::content), \(x) x$taxon_id)
metadata.mammal$taxon_id <- taxon.ids

# Find the pairwise distances between each species
pairwise.species <- expand.grid(unique(metadata.mammal$taxon_id), unique(metadata.mammal$taxon_id)) %>%
  dplyr::filter(Var1!=Var2, Var1<Var2) # only call each pair once

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
}

# check the tree is sensible
# pairwise.times <- read.time.tree.data()
# Cluster a matrix of times
# pairwise.time.matrix <- pairwise.times %>%
#   dplyr::select(scientific_name_a, scientific_name_b, age) %>%
#   tidyr::pivot_wider(names_from = scientific_name_b, values_from = age,
#                      names_sort = TRUE
#                      ) %>%
#   dplyr::arrange(scientific_name_a) %>%
#   dplyr::select(-scientific_name_a) %>%
#   as.matrix()
# rownames(pairwise.time.matrix) <- colnames(pairwise.time.matrix)
# plot(hclust(as.dist(pairwise.time.matrix)))

#### Prepare codeml site model to check for site-specific selection ####

# Adapted from Beginner's Guide on the Use of PAML to Detect Positive Selection
# https://academic.oup.com/mbe/article/40/4/msad041/7140562 for details

# Two files need to be uploaded to the cluster for running codeml:
# - alignment
# - tree file (unrooted)

# The treefile created earlier needs node names and branch lengths removing.
# Set node names to the empty string, and set branch lengths to 1
labelled.tree <- ape::read.tree(paste0(files$mammal.nt.aln, ".treefile"))
labelled.tree$node.label <- rep("", length(labelled.tree$node.label))
labelled.tree$edge.length <- rep(1, length(labelled.tree$edge.length))
ape::write.tree(labelled.tree, file = "paml/site-specific/zfxy.nt.aln.paml.treefile")

# We now read the tree as a raw string to remove the branch lengths entirely,
# separate the labels from taxon names and resave
labelled.tree <- readr::read_file("paml/site-specific/zfxy.nt.aln.paml.treefile")
labelled.tree <- gsub(":1", "", labelled.tree)
readr::write_file(labelled.tree, "paml/site-specific/zfxy.nt.aln.paml.treefile")

# This control file tests site models With heterogeneous ω Across Sites
paml.site.file <- paste0("seqfile   = ../../", files$mammal.nt.aln, " * alignment file\n",
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

#### Prepare codeml branch-site model to look for selection specifically in Muroidea ####

# To look at the rodent clade, we need a rooted tree

nt.aln.tree <- ape::read.tree(files$mammal.nt.aln.treefile)
# Root the tree on platypus and resave
# The root is arbitrarily placed in the platypus branch to fit neatly
nt.aln.tree <- phytools::reroot(nt.aln.tree, which(nt.aln.tree$tip.label=="Platypus_ZFX"), position = 0.015)
ape::write.tree(nt.aln.tree, file = paste0(files$mammal.nt.aln, ".rooted.treefile"))

# Find the nodes that are ZFY vs ZFX and add to tree
mammal.gene.groups <- split(metadata.mammal$common.name, metadata.mammal$group)
nt.aln.tree <- tidytree::groupOTU(nt.aln.tree, mammal.gene.groups, group_name = "group")
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

# This control file tests site models With heterogeneous ω across Sites
paml.branch.site.file <- paste0("seqfile   = ../../", files$mammal.nt.aln, " * alignment file\n",
                                "treefile  = zfxy.nt.aln.paml.fg.treefile * tree in Newick format without nodes\n", # 
                                "outfile   = branch-site.paml.out.txt\n",
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
write_file(paml.branch.site.file, "paml/branch-site/zfy.branch-site.paml.ctl")

# Add a script to submit these jobs to the cluster
# Run the commands manually
paml.shell.script <- paste0("#!/bin/bash\n\n",
                            "# qsubme is an alias to submit the job to the cluster\n",
                            "shopt -s expand_aliases\n",
                            "source ~/.bashrc\n\n",
                            "cd paml/site-specific\n",
                            "# codeml zfy.site-specific.paml.ctl\n",
                            "qsubme codeml zfy.site-specific.paml.ctl\n\n",
                            "cd ../branch-site\n",
                            "# codeml zfy.branch-site.paml.ctl\n",
                            "qsubme codeml zfy.branch-site.paml.ctl\n"
                            
)
write_file(paml.shell.script, "run_paml.sh")

paml.branch.site.output <- "paml/branch-site/zfy.branch-site.positive.sites.txt"

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
ape::write.tree(hyphy.tree, file = "aln/hyphy/mammal.nt.aln.hyphy.treefile")

# create an sh file to submit these
# HyPhy should be installed in a conda environment
# This script should be invoked within the conda environment
hyphy.sh.data <- paste0("#!/bin/bash\n",
                        "source activate hyphy\n",
                        "hyphy relax --alignment paml/exon_1_3-6/exon_1_3-6.aln --tree aln/hyphy/mammal.nt.aln.hyphy.treefile --reference 'Reference' --test 'Test' --output aln/hyphy/exon_1_3-6.relax.json\n",
                        "hyphy relax --alignment paml/exon_2/exon_2.aln --tree aln/hyphy/mammal.nt.aln.hyphy.treefile --reference 'Reference' --test 'Test' --output aln/hyphy/exon_2.relax.json\n",
                        "hyphy relax --alignment paml/exon_7/exon_7.aln --tree aln/hyphy/mammal.nt.aln.hyphy.treefile --reference 'Reference' --test 'Test' --output aln/hyphy/exon_7.relax.json\n",
                        "hyphy relax --alignment  ", files$mammal.nt.aln, " --tree aln/hyphy/mammal.nt.aln.hyphy.treefile --reference 'Reference' --test 'Test' --output aln/hyphy/mammal.relax.json\n"
)
write_file(hyphy.sh.data, "hyphy.sh")

system2("bash", "hyphy.sh")

#### Run HyPhy GARD to test for recombination ####

# Output was entirely unconvincing - drop this from the analysis
# if(!installr::is.windows()){
#   
#   system2("hyphy", " gard --alignment aln/zfxy.nt.aln --type codon --output hyphy/mammal.gard.json")
#   
#   # Read the json file and parse results
#   hyphy.data <- jsonlite::read_json("hyphy/mammal.gard.json")
#   
#   siteBreakPointSupport <- data.frame("site" = as.numeric(names(hyphy.data$siteBreakPointSupport)),
#                                          "support" = unlist(hyphy.data$siteBreakPointSupport))
#   
#   ggplot() +
#     geom_rect(data = mouse.exons,          aes(xmin = start_aa, xmax = end_aa, ymin = 2, ymax = 3, fill = is_even))+
#     scale_fill_manual(values=c("grey", "white"))+
#     geom_rect(data = ranges.ZF.common,     aes(xmin=start, xmax=end, ymin=0, ymax=1), fill="grey", alpha=0.5)+
#     geom_rect(data = ranges.9aaTAD.common, aes(xmin=start, xmax=end, ymin=0, ymax=1), fill="blue", alpha=0.5)+
#     geom_rect(data = ranges.NLS.common,    aes(xmin=start, xmax=end, ymin=0, ymax=1), fill="green", alpha=0.5)+
#     geom_vline(xintercept = 488)+
#     geom_vline(xintercept = 884)+
#     scale_x_continuous(expand = c(0,0))+
#     guides(fill = "none")+
#     # geom_point(data =siteBreakPointSupport, aes(x = site, y = log10(support) ))+
#     theme_bw()+
#     theme(axis.text.y = element_blank())
#   
#   tree.1.488 <- hyphy.data$breakpointData[[1]]$tree
#   readr::write_file(tree.1.488, file = "hyphy/tree.1.488.treefile")
#   tree.1.488 <- read.tree("hyphy/tree.1.488.treefile")
#   
#   tree.489.884 <- hyphy.data$breakpointData[[2]]$tree
#   readr::write_file(tree.489.884, file = "hyphy/tree.489.884.treefile")
#   tree.489.884 <- read.tree("hyphy/tree.489.884.treefile")
#     
#   # Trees are entirely unconvincing, leave this
#   
# }


#### Tar the outputs ####
system2("tar", "czf aln.tar.gz aln")
cat("Done!\n")