# Run GENECONV to test for gene conversion
# GENECONV is expected on the path
source("src/functions.R")

filesstrings::dir.remove("aln/zfx_only")
filesstrings::dir.remove("aln/zfy_only")
filesstrings::create_dir("aln/zfx_only")
filesstrings::create_dir("aln/zfy_only")

ALIGNMENTS <- read.alignments()

#### Run geneconv with the ancestral nodes and real sequences ####
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
  annotate("rect", xmin=0.12, ymin=25.8, xmax=0.19, ymax=27.5, fill="darkgreen", alpha=0.4)+
  annotate("text", x=0.13, y=27, label="Ssty appears", size=2, hjust=0)+
  annotate("text", x=0.13, y=26.2, label="Zfy testis specific", size=2, hjust=0)+
  annotate("rect", xmin=0.29, ymin=29.5, xmax=0.34, ymax=31, fill="darkgreen", alpha=0.4)+
  annotate("text", x=0.295, y=30.5, label="Sly amplifies", size=2, hjust=0)+
  annotate("rect", xmin=0.24, ymin=27.4, xmax=0.27, ymax=29, fill="darkgreen", alpha=0.4)+
  annotate("text", x=0.242, y=28, label="Slxl1\nacquired?", size=2, hjust=0)+
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
  geom_treescale(fontsize =2, y = -1, width = 10) +
  coord_cartesian(xlim = c(-5, 210))+
  theme_tree() +
  theme(legend.position = c(0.2, 0.8),
        legend.background = element_blank(),
        legend.text = element_text(size=6),
        legend.title = element_text(size=6))
save.double.width("figure/subs.per.site.time.png", time.plot)

