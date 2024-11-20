# Run GENECONV to test for gene conversion
# GENECONV is expected on the path
source("src/functions.R")

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
  # dplyr::filter(str_detect(Type, "Global" )) %>%
  dplyr::mutate(
                # Get just the ZFX node names for ancestral nodes
                SingleSpecies = ifelse(str_detect(Species, ";"), str_extract(Species, "^[^;]+"), Species),
                SingleSpecies = str_replace(SingleSpecies, "ZFX_", ""),
                isNode = str_detect(Species, "ZFX"),
                nodeType = ifelse(isNode, "Ancestral node", "Species"),
                ManualPosition = as.factor(SingleSpecies),
                ManualPosition = fct_relevel(ManualPosition,
                                             "African_Grass_Rat",
                                             "Mouse",
                                             "Mus-Arvicanthis",
                                             "Rat",
                                             "Murinae",
                                             "Mongolian_gerbil",
                                             "Muridae",
                                             "North_American_deer_mouse",
                                             "Cricetidae",
                                             "Eumuroida",
                                             "Beaver",
                                             "Arctic_ground_squirrel",
                                             "Xerinae",
                                             "Gray_squirrel",
                                             "Sciuridae",
                                             "Simiiformes",
                                             "Cattle",
                                             "Polar_bear",
                                             "Stoat",
                                             "Arctoidea",
                                             "Caniformia",
                                             "Cat",
                                             "Carnivora",
                                             "Laurasiatheria",
                                             "Boreoeutheria",
                                             "African_bush_elephant",
                                             "Atlantogenata",
                                             "Eutheria"
                                             
                )) %>%
  # What are the species descending from each of the ancestral nodes? Make labels clearer
  dplyr::rowwise() %>%
  dplyr::group_by(ManualPosition) %>%
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
  
  # Add shaded boxes to show the phylogeny
  annotate("rect", xmin = -880, xmax = 0, ymin = -1, ymax = 29, fill = "grey", col="black", alpha=0.5)+ # Eutheria
  annotate("rect", xmin = -850, xmax = 0, ymin = 25.6, ymax = 27.5, fill = "grey", col="black",  alpha=0.5)+ # Atlantogenata
  annotate("rect", xmin = -850, xmax = 0, ymin = 0, ymax = 25.4, fill = "grey", col="black",  alpha=0.5)+ # Boreoeutheria
  annotate("rect", xmin = -750, xmax = 0, ymin = 16.5, ymax = 24.5, fill = "grey",col="black",  alpha=0.5)+ # Laurasiatheria
  
  annotate("rect", xmin = -700, xmax = 0, ymin = 17.5, ymax = 23.5, fill = "grey",col="black",  alpha=0.5)+ # Carnivora
  annotate("rect", xmin = -700, xmax = 0, ymin = 11.5, ymax = 15.5, fill = "grey",col="black",  alpha=0.5)+ # Sciuridae
  annotate("rect", xmin = -780, xmax = 0, ymin = 0.25, ymax = 10.5, fill = "grey",col="black",  alpha=0.5)+ # Eumuroida
  annotate("rect", xmin = -700, xmax = 0, ymin = 0.5, ymax = 7.5, fill = "grey",col="black",  alpha=0.5)+ # Muridae
  
  labs(x = "Position", y = "Species")+
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 2800, 200))+
  scale_y_discrete(labels = function(x) str_wrap( str_replace_all(x, "_", " "), width = 30))+
  
  # Draw the structures
  add.track(ranges.ZF.common,     -Inf, Inf, start_col = "start_nt", end_col = "end_nt", fill=ZF.COLOUR)+
  add.track(ranges.NLS.common,    -Inf, Inf,  start_col = "start_nt", end_col = "end_nt",fill=NLS.COLOUR, alpha = 0.5)+
  add.track(ranges.9aaTAD.common,  -Inf, Inf,   start_col = "start_nt", end_col = "end_nt",fill=TAD.COLOUR,  alpha = 0.9)+ # fill color max from "grDevices::Blues 3"
  add.track.labels(ranges.9aaTAD.common,  0.5, 1, start_col = "start_nt", end_col = "end_nt", col="white")+   # Label the 9aaTADs
  add.track.labels(ranges.ZF.common,  0.5, 1, start_col = "start_nt", end_col = "end_nt",
                   label_col = "motif_number", col="black")+   # Label the ZFs
  
  new_scale_fill()+
  scale_fill_manual(values=c("white", "grey", "white", "grey", "white", "grey", "white"))+
  scale_pattern_color_manual(values=c("white", "white"))+
  scale_pattern_manual(values = c("none", "stripe")) + # which exons are patterned
  guides(fill = "none", pattern="none")+
  add.exon.track(-0.5, 0.5, start_col = "start_nt", end_col = "end_nt", col = "black")+
  add.exon.labels(-0.5, 0.5, start_col = "start_nt", end_col = "end_nt")+
  
  # Draw the fragments
  geom_segment(aes(x=begin, y = ManualPosition, xend = end, yend = ManualPosition, col=Type), linewidth = 2, alpha = 0.5)+
  scale_colour_manual(values = c("Pairwise inner"="orange", 
                                 "Pairwise outer"="green", 
                                 "Global inner"="blue", 
                                 "Global outer"="red"))+
  coord_cartesian(xlim = c(0, max(mouse.exons$end_nt)), clip = "off")+
  
  # Internal phylogeny boxes
  annotate("rect", xmin = 0, xmax = Inf, ymin = 17.5, ymax = 23.5, col="black", fill=NA)+ # Carnivora
  annotate("rect", xmin = 0, xmax = Inf, ymin = 11.5, ymax = 15.5, fill = NA ,col="black",  alpha=0.5)+ # Sciuridae
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0.5, ymax = 7.5, fill = NA,col="black",  alpha=0.5)+ # Muridae
 
  # facet_wrap(~Type, ncol = 1)+

  theme_bw()+
  theme(axis.title = element_blank(),
        legend.position = "top",
        axis.text.y.left = element_text(hjust = 1))

save.double.width("figure/aln.ancestral.zfxy_geneconv.png", aln.geneconv.plot)

