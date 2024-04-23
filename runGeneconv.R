# Run GENECONV to test for gene conversion
# GENECONV is expected on the path
source("functions.R")

filesstrings::dir.remove("aln/zfx_only")
filesstrings::dir.remove("aln/zfy_only")
filesstrings::create_dir("aln/zfx_only")
filesstrings::create_dir("aln/zfy_only")

alignments <- list()
alignments$nt.mammal.ape <- ape::read.FASTA(files$mammal.nt.aln)

#### Read the aligned FASTA files and metadata, export ZFX and ZFY separately ####
fa.files <- list.files(path = "fasta/nt", pattern = "*.fa$", 
                       include.dirs = T, full.names = T)

fa.read  <- lapply(fa.files, read.fasta)
metadata.mammal <- read.metadata(fa.read)

# We want to look at gene conversion within species lineages. We also want to
# compare the ancestral ZFXs and ZFYs at each node To do this, create a separate
# tree for each of ZFX and ZFY, confirm that they have equivalent branches, then
# add the ancestral reconstruction to the geneconv config file.
nt.aln.file <- "aln/mammal/mammal.nt.aln"
msa.nt.aln <- Biostrings::readDNAMultipleAlignment(nt.aln.file, format="fasta")
ape.nt.aln <- ape::read.FASTA(nt.aln.file)

# Write ZFX sequences to file
zfx.nt.aln <- msa.nt.aln@unmasked[names(msa.nt.aln@unmasked) %in% c(metadata.mammal$common.name[metadata.mammal$group=="ZFX"], 
                                                                    metadata.mammal$common.name[metadata.mammal$group=="Outgroup"])]
Biostrings::writeXStringSet(zfx.nt.aln,  file = "aln/zfx_only/zfx.aln", format = "fasta")

# Write ZFY sequences to file
zfy.nt.aln <- msa.nt.aln@unmasked[names(msa.nt.aln@unmasked) %in% c(metadata.mammal$common.name[metadata.mammal$group=="ZFY"], 
                                                                    metadata.mammal$common.name[metadata.mammal$group=="Outgroup"])]
Biostrings::writeXStringSet(zfy.nt.aln,  file = "aln/zfy_only/zfy.aln", format = "fasta")

#### Create independent trees for ZFX and ZFY sequences ####

# We can specify the true phylogeny of the sequences for ancestral reconstruction
zfx.phylogeny <- paste0("(Platypus_ZFX, (Opossum_ZFX, ", # Outgroups
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
                              "(Beaver_Zfx, (", # Muroidea
                                "(North_American_deer_mouse_Zfx, Desert_hamster_Zfx)Cricetidae,", # Cricetidae 
                                "(Mongolian_gerbil_Zfx, (Rat_Zfx, (Mouse_Zfx, African_Grass_Rat_Zfx)Mus-Arvicanthis)Murinae)Muridae", # Muridae 
                              ")Eumuroida)Muroidea", # /Muroidea
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
                        ")Eutheria", # /Eutheria
                        ")Theria)Mammalia;") # /Outgroups

write_file(zfx.phylogeny, "aln/zfx_only/zfx.nt.species.nwk")

zfy.phylogeny <- paste0("(Platypus_ZFX, (Opossum_ZFX, ", # Outgroups
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
                        "(Beaver_Zfx-like_putative-Zfy, (", # Muroidea
                        "(North_American_deer_mouse_Zfx-like_putative-Zfy, Desert_hamster_Zfx-like_putative-Zfy)Cricetidae,", # Cricetidae 
                        "(Mongolian_gerbil_Zfx-like_putative-Zfy, (Rat_Zfy2, ((Mouse_Zfy1, Mouse_Zfy2), (African_Grass_Rat_ZFY2-like_1, African_Grass_Rat_ZFY2-like_2))Mus-Arvicanthis)Murinae)Muridae", # Muridae 
                        ")Eumuroida)Muroidea", # /Muroidea
                        ")Rodentia", # /Rodentia
                        ")Euarchonoglires,", # /Euarchonoglires
                        "(", # Laurasiatheria
                        "(",  # Carnivora
                        "Cat_ZFY, (Dog_ZFY, (Stoat_ZFY, Polar_bear_ZFY)Arctoidea)Caniformia",
                        ")Carnivora,",  # /Carnivora
                        "(",  # Euungulata 
                        "Horse_ZFY, (Pig_ZFY, (White_tailed_deer_ZFY, (Cattle_ZFY, Goat_ZFY)Bovidae)Pecora)Artiodactyla",
                        ")Euungulata", # /Euungulata 
                        ")Laurasiatheria", # /Laurasiatheria
                        ")Boreoeutheria", # /Boreoeutheria
                        ")Eutheria", # /Eutheria
                        ")Theria)Mammalia;") # /Outgroups

write_file(zfy.phylogeny, "aln/zfy_only/zfy.nt.species.nwk")

# Run the ancestral reconstructions
system2("iqtree", paste("-s ", "aln/zfx_only/zfx.aln", 
                        # "-bb 1000 -alrt 1000", # bootstrapping
                        "-nt AUTO", # number of threads
                        "-te aln/zfx_only/zfx.nt.species.nwk", # user tree guide
                        "-asr")) # ancestral sequence reconstruction

system2("iqtree", paste("-s ", "aln/zfy_only/zfy.aln", 
                        # "-bb 1000 -alrt 1000", # bootstrapping
                        "-nt AUTO", # number of threads
                        "-te aln/zfy_only/zfy.nt.species.nwk", # user tree guide
                        "-asr")) # ancestral sequence reconstruction

#### Remove duplicate species nodes from the trees  ####

# For a ~species tree, we only need one sequence per species. Either can be dropped,
# since they give the same branching order

# Read the ML ZFX and ZFY trees 
zfx.nt.aln.tree <- ape::read.tree("aln/zfx_only/zfx.aln.treefile")
zfy.nt.aln.tree <- ape::read.tree("aln/zfy_only/zfy.aln.treefile")

# zfx.nt.aln.tree <- ape::read.tree(text = zfx.phylogeny)
# zfy.nt.aln.tree <- ape::read.tree(text = zfy.phylogeny)

# Drop the second ZFYs in mouse and rat
zfy.nt.aln.tree <- tidytree::drop.tip(zfy.nt.aln.tree, "Mouse_Zfy2") 
zfy.nt.aln.tree <- tidytree::drop.tip(zfy.nt.aln.tree, "African_Grass_Rat_ZFY2-like_1") 

# Root the trees on platypus
zfx.nt.aln.tree <- phytools::reroot(zfx.nt.aln.tree, which(zfx.nt.aln.tree$tip.label=="Platypus_ZFX"), position = 0.015)
zfy.nt.aln.tree <- phytools::reroot(zfy.nt.aln.tree, which(zfy.nt.aln.tree$tip.label=="Platypus_ZFX"), position = 0.015)

# Remove gene names so tip labels are comparable
zfx.nt.aln.tree$tip.label <- str_replace(zfx.nt.aln.tree$tip.label, "_Z[F|f][X|x].*", "")
zfy.nt.aln.tree$tip.label <- str_replace(zfy.nt.aln.tree$tip.label, "(_putative)?(_|-)Z[F|f][X|x|Y|y].*", "")


#### Plot ZFX / ZFY tree comparisons  ####

# Export comparison of the ML trees
png(filename = "figure/zfx.zfy.nt.ancestral.treediff.png")
treespace::plotTreeDiff(zfx.nt.aln.tree, zfy.nt.aln.tree, treesFacing=TRUE)
dev.off()

# Plot the two trees with node labels
zfx.nt.aln.tree.plot <- plot.tree(zfx.nt.aln.tree) + geom_nodelab(size=2, nudge_x = -0.003, nudge_y = 0.5, hjust=1,  node = "internal")
zfy.nt.aln.tree.plot <- plot.tree(zfy.nt.aln.tree) + geom_nodelab(size=2, nudge_x = -0.003, nudge_y = 0.5, hjust=1,  node = "internal")
# zfx.zfy.aln.tree.plot <- zfx.nt.aln.tree.plot + zfy.nt.aln.tree.plot + patchwork::plot_annotation(tag_levels = list(c("ZFX", "ZFY")))
save.double.width("figure/ancestral.zfx.png", zfx.nt.aln.tree.plot)
save.double.width("figure/ancestral.zfy.png", zfy.nt.aln.tree.plot)

#### Find matching nodes in the two trees and get ancestral reconstructions ####

# Nodes with the same branching order will be viable for comparison of ancestral
# reconstructions

# Find the matching nodes
zfy.zfx.common.nodes <- ape::comparePhylo(zfx.nt.aln.tree, zfy.nt.aln.tree)$NODES %>%
  tidyr::separate_wider_delim(cols = zfx.nt.aln.tree, delim = " ", names = c("ZFX_node", "zfx_nnodes") ) %>%
  tidyr::separate_wider_delim(cols = zfy.nt.aln.tree, delim = " ", names = c("ZFY_node", "zfy_nnodes") ) %>%
  # Remove booststrap values from node names
  dplyr::mutate(ZFX_node = str_replace(ZFX_node, "\\/\\d.*", ""),
                ZFY_node = str_replace(ZFY_node, "\\/\\d.*", "")) %>% 
  dplyr::filter(ZFX_node!="Mammalia" & ZFX_node!="Root" & ZFX_node!="") 

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

# Combine with existing ZFX/Y sequence file
ape::write.dna(alignments$nt.mammal.ape, "aln/ancestral.zfx.zfy.nodes.fa", format="fasta", 
               append = TRUE, colsep = "", nbcol=-1)

#### Run geneconv with the ancestral nodes and real sequences ####

# Write geneconv configuration file specifying which sequences are in the same
# group
node.string.names <- zfy.zfx.common.nodes %>% 
  dplyr::mutate(GroupString= paste0("-group ZFX_", ZFX_node, "_ZFY_",ZFY_node, " ZFX_", ZFX_node, " ZFY_",ZFY_node ))
node.string <- paste(node.string.names$GroupString, collapse = "\n")

split.gene.names <- metadata.mammal %>%
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

#  Only on PATH in Windows, not running on cluster
system2("geneconv", paste("aln/ancestral.zfx.zfy.nodes.fa",
                          "aln/anc.zfx.zfy.geneconv.cfg",
                          "aln/anc.zfx.zfy.geneconv.out /lp"),
        stdout = "aln/ancestral.zfx.zfy.nodes.geneconv.log", 
        stderr = "aln/ancestral.zfx.zfy.nodes.geneconv.log")

#### Read geneconv output and plot potential gene conversion fragments ####


get.offspring.of.node <- \(i){ 
  tips <- offspring(zfx.nt.aln.tree, length(zfx.nt.aln.tree$tip.label)+i, type="tips")
  paste(zfx.nt.aln.tree$tip.label[tips], collapse = ", ")
}

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
  dplyr::mutate(
    # y = row_number(),
                # Get just the ZFX node names for ancestral nodes
                SingleSpecies = ifelse(str_detect(Species, ";"), str_extract(Species, "^[^;]+"), Species),
                SingleSpecies = str_replace(SingleSpecies, "ZFX_", ""),
                isNode = str_detect(Species, "ZFX"),
                nodeType = ifelse(isNode, "Ancestral node", "Species")) %>%
  # What are the species descending from each of the ancestral nodes? Make labels clearer
  dplyr::rowwise() %>%
  dplyr::group_by(SingleSpecies) %>%
  dplyr::mutate(label = cur_group_id()) # same y axis value for each species




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
# Plot the geneconv fragments in the ancestral nodes

aln.geneconv.plot <- ggplot(aln.geneconv.data)+
  
  labs(x = "Position", y = "Species")+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 2800, 200))+
  scale_y_discrete(labels = function(x) str_wrap( str_replace_all(x, "_", " "), width = 30))+
  
  # Draw the structures
  add.track(ranges.ZF.common,     -Inf, Inf, start_col = "start_nt", end_col = "end_nt", fill="lightgrey")+
  add.track(ranges.NLS.common,    -Inf, Inf,  start_col = "start_nt", end_col = "end_nt",fill="green", alpha = 0.5)+
  add.track(ranges.9aaTAD.common,  -Inf, Inf,   start_col = "start_nt", end_col = "end_nt",fill="#00366C",  alpha = 0.9)+ # fill color max from "grDevices::Blues 3"
  add.track.labels(ranges.9aaTAD.common,  0.5, 1, start_col = "start_nt", end_col = "end_nt", col="white")+   # Label the 9aaTADs
  add.track.labels(ranges.ZF.common,  0.5, 1, start_col = "start_nt", end_col = "end_nt",
                   label_col = "motif_number", col="black")+   # Label the ZFs
  
  new_scale_fill()+
  scale_fill_manual(values=c("white", "grey", "white", "grey", "white", "grey", "white"))+
  scale_pattern_color_manual(values=c("white", "white"))+
  scale_pattern_manual(values = c("none", "stripe")) + # which exons are patterned
  guides(fill = "none", pattern="none")+
  add.exon.track(-0.5, 0.5, start_col = "start", end_col = "end", col = "black")+
  add.exon.labels(-0.5, 0.5, start_col = "start", end_col = "end")+
  
  geom_segment(aes(x=begin, y = SingleSpecies, xend = end, yend = SingleSpecies), linewidth = 2)+
  coord_cartesian(xlim = c(0, max(mouse.exons$end)))+
  facet_wrap(~nodeType, scale="free", ncol = 1)+
  theme_bw()+
  theme(axis.title = element_blank())

save.double.width("figure/aln.ancestral.zfxy_geneconv.png", aln.geneconv.plot)


#### How do substitution rates compare to the divergence times? ####
# We want to calculate the number of substitutions per million years
# Combine the branch lengths with the TimeTree dates

pairwise.times <- read.time.tree.data()
mammal.nt.tree.data <- tidytree::as_tibble(zfy.nt.aln.tree)

# Manually match the divergence times from the TimeTree data
divergence.point.data <- c(
  "Mammalia" = 180.06610,
  "Theria" = 160,
  "Eutheria" = 99.18870,
  "Atlantogenata" = 97,
  "Boreoeutheria" = 94.00000,
  "Laurasiatheria" = 76.00000,
  "Carnivora" = 55.36164,
  "Caniformia" = 45.10000,
  "Arctoidea" = 40.12000,
  "Euungulata" = 72.7,
  "Artiodactyla" = 61.84265,
  "Pecora" = 23.7,
  "Bovidae" = 21.62563,
  "Euarchonoglires" = 87.20000,
  "Simiiformes" = 42.90000,
  "Catarrhini" = 28.82000,
  "Hominidae" = 8.60000,
  "Hominini" = 6.40000,
  "Cercopithecidae" = 17.75500,
  "Cercopithecinae" = 10.45400,
  "Rodentia" = 70.20250,
  "Sciuridae" = 34.46259,
  "Xerinae" = 11.38264,
  "Muroidea" = 68.31756,
  "Eumuroida" = 26.2,
  "Cricetidae" = 18.6,
  "Muridae" = 12.44458,
  "Murinae" = 11.64917,
  "Mus-Arvicanthis" = 10.05224
)

divergence.point.species <- rep(0, length(zfy.nt.aln.tree$tip.label))
names(divergence.point.species) <- zfy.nt.aln.tree$tip.label
divergence.point.data <- c(divergence.point.data, divergence.point.species)

# TODO - Use these times to create weights for each edge - plot like the relax K trees
zfy.nt.aln.tree$node.label[1] <- "Mammalia" # replace default 'Root' label
zfy.nt.aln.tree$tip.label <- str_replace_all(zfy.nt.aln.tree$tip.label, "_", " ")

get.edge.time <- function(node, tree){
  nodelabel <- treeio::nodelab(tree, node)
  parentnode <- tree$edge[tree$edge[,2]==node,1] # parent of edge leading to node
  parentlabel <- treeio::nodelab(tree, parentnode) # label of the parent

  v1 = max(divergence.point.data[parentlabel], divergence.point.data[nodelabel])
  v2 = min(divergence.point.data[parentlabel], divergence.point.data[nodelabel])
  time = v1 - v2
  edge.row <-  ifelse(length(parentnode)==0,0, which(tree$edge[,2]==node))
  
  edgelength <- ifelse(length(parentnode)==0,0, tree$edge.length[edge.row])
  subsPerMyr <- ifelse(length(parentnode)==0,0, edgelength / time)

  return(data.frame(parentnode = ifelse(length(parentnode)==0, NA, parentnode),
                    parentlabel =ifelse(length(parentnode)==0, NA, parentlabel),
                    node = node,
                    nodelabel = nodelabel,
                    name = paste0(parentlabel, "-", nodelabel),
                    parent.time = v1,
                    node.time = v2,
                    time = time,
                    edge.row,
                    subsPerSite = edgelength,
                    subsPerMyr = subsPerMyr))
}

time.vals <- do.call(rbind, lapply(1:59, get.edge.time, tree=zfy.nt.aln.tree))

subs.site.mya.plot <- ggtree(zfy.nt.aln.tree, size = 1) %<+%
  time.vals +
  aes(colour = log(subsPerMyr)) +
  scale_color_paletteer_c("ggthemes::Classic Red-Blue", 
                          direction = -1, limits =c(-9, -3))+
  labs(color = "Log substitutions per site\nper million years")+
  geom_nodelab(size=2, nudge_x = -0.005, nudge_y = 0.5, hjust = 1, color = "black")+
  geom_tiplab(size=2, color = "black")+
  geom_treescale(fontsize =2, y = -1) +
  coord_cartesian(xlim = c(-0.05, 0.4))+
  theme_tree() +
  theme(legend.position = c(0.2, 0.8),
        legend.background = element_blank(),
        legend.text = element_text(size=6),
        legend.title = element_text(size=6))
save.double.width("figure/subs.per.site.per.Myr.png", subs.site.mya.plot)

# Redraw the tree with branch lengths from the actual dates. Confirm it adds up
# to about the same length per species.

zfy.nt.aln.tree.time <- zfy.nt.aln.tree
get.time.for.node <- function(node) time.vals[time.vals$node==node,"time"]
zfy.nt.aln.tree.time$edge.length <- sapply(zfy.nt.aln.tree.time$edge[,2], get.time.for.node)

time.plot <- ggtree(zfy.nt.aln.tree.time, size = 1) %<+%
  time.vals +
  aes(colour = log(subsPerMyr)) +
  scale_color_paletteer_c("ggthemes::Classic Red-Blue", 
                          direction = -1, limits =c(-9, -3))+
  labs(color = "Log substitutions per site\nper million years")+
  geom_nodelab(size=2, nudge_x = -3, nudge_y = 0.5, hjust = 1, color = "black")+
  geom_tiplab(size=2, color = "black")+
  geom_treescale(fontsize =2, y = -1) +
  coord_cartesian(xlim = c(-5, 210))+
  theme_tree() +
  theme(legend.position = c(0.2, 0.8),
        legend.background = element_blank(),
        legend.text = element_text(size=6),
        legend.title = element_text(size=6))
save.double.width("figure/subs.per.site.time.png", time.plot)

