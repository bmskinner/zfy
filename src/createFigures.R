# Create figures

#### Imports #####

source("src/functions.R")
load.packages()

source("src/find9aaTADs.R")
source("src/findZF.R")
source("src/calcCharge.R")
source("src/calcHydrophobicity.R")

cat("Packages loaded\n")
METADATA <- prepare.fas.files() # load FASTA files and write metadata table
ALIGNMENTS <- read.alignments()

# Identify the coordinates of the exon boundaries in the gapped alignments
# Based on the mouse Zfy1 sequence
mouse.exons <- find.exons()

#### Plot combined mammal/outgroup AA tree ####

combined.outgroup.tree <- read.combined.outgroup.tree(FILES$combined.aa.aln.treefile)

combined.aa.tree <- plot.tree(combined.outgroup.tree, col = "group")+
  coord_cartesian(clip="off", xlim = c(0, 0.85), ylim= c(-2, length(combined.outgroup.tree$tip.label)))
  
# Find the location of the Eumuroida clade in the plot
eumuroida.node <- ape::getMRCA(combined.outgroup.tree, c("Mouse_Zfy1", "North_American_deer_mouse_Zfx-like_putative-Zfy"))
eumuroida.node.position <- ggtree::get_clade_position(combined.aa.tree, eumuroida.node)

# Add text annotation just below the y midpoint of the Eumuroida clade
combined.aa.tree <- combined.aa.tree + annotate("text", x = 0.3, 
                                                y = (eumuroida.node.position$ymax+eumuroida.node.position$ymin)/2-1, 
                                                label = "Eumuroida", size=2)

save.double.width("figure/Figure_1_aa_tree.png", combined.aa.tree)

# Also save the order of taxa in the outgroup tree to use later
combined.taxa.name.order <- ggtree::get_taxa_name(combined.aa.tree) 


#### Plot mammal CDS NT tree #####

nt.aln.tree <- ape::read.tree(FILES$combined.nt.aln.treefile)

# Root the tree and resave
# The root is arbitrarily placed to fit neatly
xenopus.node <- ape::getMRCA(nt.aln.tree, c("Xenopus_ZFX.S","Xenopus_ZFX.L"))
nt.aln.tree <- phytools::reroot(nt.aln.tree, xenopus.node, position = 0.01)
ape::write.tree(nt.aln.tree, file = paste0(FILES$combined.nt.aln, ".rooted.treefile"))

# Find the nodes that are ZFY vs ZFX and add to tree
mammal.gene.groups <- split(METADATA$combined$common.name, METADATA$combined$group)
nt.aln.tree <- tidytree::groupOTU(nt.aln.tree, mammal.gene.groups, group_name = "group")

plot.zfx.zfy <- plot.tree(nt.aln.tree, col= "group") + coord_cartesian(clip="off", xlim = c(0, 0.65)) 

# Find the location of the Eumuroida clade in the plot
eumuroida.node <- ape::getMRCA(nt.aln.tree, c("Mouse_Zfy1", "North_American_deer_mouse_Zfx-like_putative-Zfy"))
eumuroida.node.position <- ggtree::get_clade_position(plot.zfx.zfy, eumuroida.node)


plot.zfx.zfy <- plot.zfx.zfy + annotate("text", x = 0.3, 
                                        y = (eumuroida.node.position$ymax+eumuroida.node.position$ymin)/2-1, 
                                        label = "Eumuroida", size=2)

# Emphasise the ZFX / ZFY splits more by rotating the Laurasiatheria node
# laurasiatheria.node <- ape::getMRCA(nt.aln.tree, c("Cat_ZFY", "Cat_ZFX"))
# plot.zfx.zfy <- rotate(plot.zfx.zfy, laurasiatheria.node)

save.double.width("figure/Figure_xxxx_nucleotide_tree.png", plot.zfx.zfy)

# Get the order of taxa names for reordering other plots later
mammal.taxa.name.order <- get_taxa_name(plot.zfx.zfy) 

# Now we want to view separate trees for Zfx and Zfy based on this combined tree

# Drop the ZFY sequences and just look at the ZFX nodes in the tree
tree.zfx <- ape::drop.tip(nt.aln.tree, mammal.gene.groups$ZFY)
tree.zfx <- groupOTU(tree.zfx, mammal.gene.groups, group_name = "group")
ape::write.tree(tree.zfx, file = paste0(FILES$mammal.nt.aln, ".zfx.treefile"))
plot.zfx <- plot.tree(tree.zfx) + coord_cartesian(clip="off", xlim = c(0, 0.65))
save.double.width("figure/mammal.zfx.tree.png", plot.zfx)

# Keep Zfy and outgroups, drop other tips
tree.zfy <- ape::keep.tip(nt.aln.tree, c(mammal.gene.groups$ZFY, mammal.gene.groups$Outgroup))
tree.zfy <- groupOTU(tree.zfy, mammal.gene.groups, group_name = "group")
ape::write.tree(tree.zfy, file = paste0(FILES$mammal.nt.aln, ".zfy.treefile"))
plot.zfy <- plot.tree(tree.zfy) + coord_cartesian(clip="off", xlim = c(0, 0.65))
save.double.width("figure/mammal.zfy.tree.png", plot.zfy)

#### Plot mammal exon NT trees ####

create.exon.plot <- function(name){
  exon.aln.file <- paste0("aln/exons/", name)

  # Some exons will fail - too many gaps
  if(!file.exists(paste0(exon.aln.file, ".treefile"))) return()
  
  exon.tree <- ape::read.tree(paste0(exon.aln.file, ".treefile"))
  # Root the tree on platypus
  exon.tree <- reroot.tree(exon.tree, c("Platypus_ZFX", "Australian_echidna_ZFX"), position = 0.015 )
  
  # Find the nodes that are ZFY vs ZFX and add to tree
  mammal.gene.groups <- split(METADATA$mammal$common.name, METADATA$mammal$group)
  exon.tree <- groupOTU(exon.tree, mammal.gene.groups, group_name = "group")
  
  plot.exon.tree <- plot.tree(exon.tree, tiplab.font.size = 1.5,  col="group")  + 
    coord_cartesian(clip="off", xlim = c(0.16, 0.83))
  # Return for playing
  plot.exon.tree
}


exon.plots <- lapply(c("exon_1-6.aln", "exon_2.aln", 
                       "exon_7.aln", "exon_1.3-6.aln"), 
                     create.exon.plot)

# Create a joint figure of exons 1-6, exon 2, and exon 7
exon.joint.tree <- exon.plots[[4]] + exon.plots[[2]] + exon.plots[[3]] + 
  patchwork::plot_annotation(tag_levels = list(c("Exons 1,3-6", "Exon 2", "Exon 7"))) &
  theme(plot.tag = element_text(size = 6),
        plot.margin = margin(t=0, l=0, r=0, b=0))
save.plot("figure/Figure_Sxxxx_exons_tree.png", exon.joint.tree, width=270, height=170)

#### Plot dN/dS globally in mammals ####

kaks.pairwise.plot <- plot.kaks(FILES$mammal.nt.aln, 
                                species.order = combined.taxa.name.order[1:which(combined.taxa.name.order=="Platypus ZFX")], # use the aa tree order so ZFX and ZFY are grouped
                                kaks.limits = c(0, 1.5))


# Annotate plot with lines to delineate eumuroida Zfx
# 
# eumuroida.zfx.ymin <- which(levels(kaks.pairwise.plot$data$row) =="Black rat Zfx")-0.5
# eumuroida.zfx.ymax <- which(levels(kaks.pairwise.plot$data$row) =="Mongolian gerbil Zfx")+0.5
# 
# eumuroida.zfx.xmin <- which(rev(levels(kaks.pairwise.plot$data$col)) =="Black rat Zfx")+0.5
# eumuroida.zfx.xmax <- which(rev(levels(kaks.pairwise.plot$data$col)) =="Mongolian gerbil Zfx")-0.5
# 
# kaks.pairwise.plot <-kaks.pairwise.plot +
#   annotate(geom = "segment", x = 0, xend = eumuroida.zfx.xmax, y = eumuroida.zfx.ymin, yend =  eumuroida.zfx.ymin)+
#   annotate(geom = "segment", x = 0, xend = eumuroida.zfx.xmin, y = eumuroida.zfx.ymax, yend =  eumuroida.zfx.ymax)+
#   annotate(geom = "segment", x = eumuroida.zfx.xmin, xend = eumuroida.zfx.xmin, y = 0, yend =  eumuroida.zfx.ymax)+
#   annotate(geom = "segment", x = eumuroida.zfx.xmax, xend = eumuroida.zfx.xmax, y = 0, yend =  eumuroida.zfx.ymin)

save.plot("figure/dnds.png", kaks.pairwise.plot, width = 270, height = 170)

# Strong purifying selection in all pairs, but weaker in rodents
#### Plot dN/dS globally partition by partition in mammals ####

# Look at the final exon versus exons 1&3-6; is purifying selection more 
# pronounced in the ZFs versus exon 2 and exon 7?
exon1.3_6.locs <- c(mouse.exons$start_nt_codon_offset_mammal[1]:mouse.exons$end_nt_codon_offset_mammal[1], 
                    mouse.exons$start_nt_codon_offset_mammal[3]:mouse.exons$end_nt_codon_offset_mammal[6])
exon1.3_6.aln <- as.matrix(ALIGNMENTS$nt.mammal.ape)[,exon1.3_6.locs]

ape::write.FASTA(exon1.3_6.aln, file = "aln/exons/exon_1_3-6.kaks.aln")
seqin.aln.exon.1.3_6 <- seqinr::read.alignment("aln/exons/exon_1_3-6.kaks.aln", format = "fasta")

exon2.aln <- as.matrix(ALIGNMENTS$nt.mammal.ape)[,mouse.exons$start_nt_codon_offset_mammal[2]:(mouse.exons$end_nt_codon_offset_mammal[2])]
ape::write.FASTA(exon2.aln, file = "aln/exons/exon_2.kaks.aln")
seqin.aln.exon.2 <- seqinr::read.alignment("aln/exons/exon_2.kaks.aln", format = "fasta")
exon2.aln.msa <- ape::read.FASTA("aln/exons/exon_2.kaks.aln")

exon.7.aln <- as.matrix(ALIGNMENTS$nt.mammal.ape)[,(mouse.exons$start_nt_codon_offset_mammal[7]):mouse.exons$end_nt_codon_offset_mammal[7]] 
ape::write.FASTA(exon.7.aln, file = "aln/exons/exon_7.kaks.aln")
seqin.aln.exon.7 <- seqinr::read.alignment("aln/exons/exon_7.kaks.aln", format = "fasta")

create.pairwise.kaks.data <- function(seqinr.aln){
  kaks.data <- seqinr::kaks(seqinr.aln, rmgap = FALSE)
  kaks.ratio <- kaks.data$ka / kaks.data$ks
  
  kaks.pairwise <- metagMisc::dist2list(kaks.ratio, tri = F) %>%
    dplyr::mutate(col = str_replace_all(col, "_", " "),
                  row = str_replace_all(row, "_", " "),
                  col = factor(col, levels = combined.taxa.name.order[1:which(combined.taxa.name.order=="Platypus ZFX")],),
                  row = factor(row, levels = combined.taxa.name.order[1:which(combined.taxa.name.order=="Platypus ZFX")],),
                  colnum = as.integer(col),
                  rownum = as.integer(row),
                  kaks  = ifelse(value==1, NA, value)) %>%  # values of exactly 1 are from missing data)
    dplyr::filter(rownum < colnum) %>%
    dplyr::select(-value)
}

kaks.exon.1.3_6 <- create.pairwise.kaks.data(seqin.aln.exon.1.3_6)
kaks.exon.2 <- create.pairwise.kaks.data(seqin.aln.exon.2)
kaks.exon.7 <- create.pairwise.kaks.data(seqin.aln.exon.7)

plot.pairwise.kaks <- function(kaks.pairwise){
  ggplot(kaks.pairwise, aes(x = col, y = row))+
    geom_tile(aes(fill=kaks))+
    scale_fill_viridis_c(limits = c(0, 1.5), direction = -1, na.value="white")+
    labs(fill="dNdS")+
    scale_x_discrete(limits=rev)+
    theme_bw()+
    theme(axis.text.x = element_text(size = 3, angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(size = 3.5),
          axis.title = element_blank(),
          legend.position = c(0.8, 0.7),
          legend.background = element_blank(),
          legend.title = element_text(size = 5),
          legend.text = element_text(size = 5))
}

exon.1.3_6.kaks.pairwise.plot <- plot.pairwise.kaks(kaks.exon.1.3_6)
exon.2.kaks.pairwise.plot <- plot.pairwise.kaks(kaks.exon.2)
exon.7.kaks.pairwise.plot <- plot.pairwise.kaks(kaks.exon.7)

exon.kaks.plot <- exon.1.3_6.kaks.pairwise.plot + exon.2.kaks.pairwise.plot + exon.7.kaks.pairwise.plot +
  patchwork::plot_annotation(tag_levels = c("A")) + plot_layout(axes="collect")

save.plot("figure/exon.1_3-6.2.7.dnds.png", exon.kaks.plot,  width = 270, height = 170)

# Create supplementary data tables with the pairwise values
create.xlsx(kaks.exon.1.3_6, "figure/kaks.exon.1.3-6.xlsx")
create.xlsx(kaks.exon.2, "figure/kaks.exon.2.xlsx")
create.xlsx(kaks.exon.7, "figure/kaks.exon.7.xlsx")

#### Identify the locations of the ZFs in the AA & NT MSAs ####

LOCATIONS <- list()
RANGES <- list()

LOCATIONS$combined.zf <- readr::read_tsv("aln/locations.zf.combined.tsv")
LOCATIONS$mammal.zf <- readr::read_tsv("aln/locations.zf.mammal.tsv")

write_tsv(LOCATIONS$combined.zf %>% 
            dplyr::select(sequence, aa_motif, start_ungapped, end_ungapped, 
                          start_gapped, end_gapped, start_nt_ungapped, end_nt_ungapped, 
                          start_nt_gapped, end_nt_gapped),
          "figure/locations.zf.tsv")

RANGES$combined.zf <- merge(find.common.aa.overlaps(LOCATIONS$combined.zf), find.common.nt.overlaps(LOCATIONS$combined.zf), by = "motif_number")
RANGES$mammal.zf <- merge(find.common.aa.overlaps(LOCATIONS$mammal.zf), find.common.nt.overlaps(LOCATIONS$mammal.zf), by = "motif_number")

# Annotate the individual ZFs with which consensus motif they belong to (0 if none)
LOCATIONS$combined.zf %<>% 
  dplyr::rowwise() %>%
  dplyr::mutate(motif_number = find.matching.range("start", "end", start_gapped, end_gapped, RANGES$combined.zf )) %>%
  dplyr::ungroup()

LOCATIONS$mammal.zf %<>% 
  dplyr::rowwise() %>%
  dplyr::mutate(motif_number = find.matching.range("start", "end", start_gapped, end_gapped, RANGES$mammal.zf )) %>%
  dplyr::ungroup()


#### Identify the locations of the 9aaTADs in the AA & NT MSAs ####

LOCATIONS$combined.9aaTAD <- readr::read_tsv("aln/locations.9aaTAD.combined.tsv")
LOCATIONS$mammal.9aaTAD <- readr::read_tsv("aln/locations.9aaTAD.mammal.tsv")

write_tsv(LOCATIONS$combined.9aaTAD %>% 
            dplyr::select(sequence, aa_motif, rc_score, start_ungapped, end_ungapped, 
                          start_gapped, end_gapped, start_nt_ungapped, end_nt_ungapped, 
                          start_nt_gapped, end_nt_gapped),
          "figure/locations.9aaTAD.tsv")

# Only keep the high confidence 9aaTADs for the track
RANGES$combined.9aaTAD <-  merge(find.common.aa.9aaTADs(LOCATIONS$combined.9aaTAD, 
                                                      rc.threshold = 80,
                                                      coverage.threshold = 21),
                               find.common.nt.9aaTADs(LOCATIONS$combined.9aaTAD, 
                                                      rc.threshold = 80,
                                                      coverage.threshold = 21), 
                               by = c("motif_number", "label"), all.x = T)

RANGES$mammal.9aaTAD <- merge(find.common.aa.9aaTADs(LOCATIONS$mammal.9aaTAD, 
                                                     rc.threshold = 80,
                                                     coverage.threshold = 21),
                              find.common.nt.9aaTADs(LOCATIONS$mammal.9aaTAD, 
                                                     rc.threshold = 80,
                                                     coverage.threshold = 21), 
                              by = c("motif_number", "label"), all.x = T)

LOCATIONS$combined.9aaTAD %<>%
  dplyr::rowwise() %>%
  dplyr::mutate(motif_number = find.matching.range("start", "end", 
                                                   start_gapped, end_gapped, RANGES$combined.9aaTAD )) %>%
  dplyr::ungroup() %>%
  merge(., RANGES$combined.9aaTAD, by = "motif_number", all.x = TRUE)

LOCATIONS$mammal.9aaTAD %<>%
  dplyr::rowwise() %>%
  dplyr::mutate(motif_number = find.matching.range("start", "end", 
                                                   start_gapped, end_gapped, RANGES$mammal.9aaTAD )) %>%
  dplyr::ungroup() %>%
  merge(., RANGES$mammal.9aaTAD, by = "motif_number", all.x = TRUE)

#### Identify the locations of the NLS in the AA & NT MSAs ####

LOCATIONS$combined.NLS <- readr::read_tsv("aln/locations.NLS.combined.tsv")
LOCATIONS$mammal.NLS <- readr::read_tsv("aln/locations.NLS.mammal.tsv")

# Export the locations of the NLS
write_tsv(LOCATIONS$combined.NLS %>% 
            dplyr::select(sequence, aa_motif, type, posterior_prob, start_ungapped, end_ungapped, 
                          start_gapped, end_gapped, start_nt_ungapped, end_nt_ungapped, 
                          start_nt_gapped, end_nt_gapped),
          "figure/locations.NLS.tsv")

# We need to combine the full set of structure locations into an overlapping set
# to be plotted in a single row. Keep those that overlap in >=5 species
RANGES$combined.nls <- merge(find.common.aa.overlaps(LOCATIONS$combined.NLS), 
                            find.common.nt.overlaps(LOCATIONS$combined.NLS), 
                            by = "motif_number")
RANGES$mammal.nls <-  merge(find.common.aa.overlaps(LOCATIONS$mammal.NLS), 
                            find.common.nt.overlaps(LOCATIONS$mammal.NLS), 
                            by = "motif_number")

# Annotate the individual NLS with which consensus motif they belong to (0 if none)
LOCATIONS$combined.NLS %<>% 
  dplyr::rowwise() %>%
  dplyr::mutate(motif_number = find.matching.range("start", "end", start_gapped, end_gapped, RANGES$combined.nls )) %>%
  dplyr::ungroup()

LOCATIONS$mammal.NLS %<>% 
  dplyr::rowwise() %>%
  dplyr::mutate(motif_number = find.matching.range("start", "end", start_gapped, end_gapped, RANGES$mammal.nls )) %>%
  dplyr::ungroup()

#### Export the residues within each ZF ####

# Full ZF motifs
LOCATIONS$combined.zf %>%
  dplyr::select(Sequence = sequence, motif_number, aa_motif) %>%
  dplyr::arrange(desc(motif_number)) %>% # so we display 13 - 1
  dplyr::mutate(motif_number = paste0("ZF_", motif_number)) %>%
  tidyr::pivot_wider(id_cols = Sequence, names_from = motif_number, values_from = aa_motif) %>%
  as.data.frame %>%
  dplyr::arrange(as.integer(Sequence)) %>%
  create.xlsx(., "figure/locations.zf.xlsx", cols.to.fixed.size.font = 2:14)

#### Export the contact bases  within each ZF ####

zf.contact.bases <- LOCATIONS$combined.zf %>%
  dplyr::select(Sequence = sequence, motif_number, contact_bases) %>%
  merge(., METADATA$combined, by.x="Sequence", by.y = "common.name") %>%
  dplyr::rename_with(., str_to_title) %>%
  dplyr::arrange(desc(Motif_number)) %>% # so we display 13 - 1
  dplyr::mutate(Motif_number = paste0("ZF_", sprintf("%02d", Motif_number))) %>%
  dplyr::select(-c(Species, Species_common_name, Accession, Original.name)) %>%
  tidyr::pivot_wider(id_cols = c(Sequence, Group), names_from = Motif_number, values_from = Contact_bases) %>%
  as.data.frame %>%
  dplyr::mutate(Sequence = factor(Sequence, levels = combined.taxa.name.order)) %>%
  dplyr::arrange(Group, as.integer(Sequence)) %>%
  dplyr::mutate(Sequence = str_replace_all(Sequence, "_", " "))

# What are the most common ZF contact motifs?
zf.contact.bases.conserved <- LOCATIONS$combined.zf %>%
  dplyr::select(Sequence = sequence, motif_number, contact_bases) %>%
  dplyr::arrange(desc(motif_number)) %>% # so we display 13 - 1
  dplyr::mutate(motif_number = paste0("ZF_", motif_number)) %>%
  dplyr::group_by(motif_number, contact_bases) %>%
  dplyr::summarise(count = n()) %>%
  dplyr::group_by(motif_number) %>%
  dplyr::arrange(desc(count)) %>%
  dplyr::slice_head(n=1)

# Create excel file, colouring the contact bases by whether they are conserved
create.contact.base.xlsx <- function(contact.bases){
  wb = xlsx::createWorkbook(type = "xlsx")
  sh = xlsx::createSheet(wb)
  xlsx::addDataFrame(contact.bases, sh, row.names = F)
  xlsx::createFreezePane(sh, 2, 2, 2, 2) # freeze top row and first column
  
  fixed.style <- xlsx::CellStyle(wb) + 
    xlsx::Font(wb, heightInPoints = 10, isBold = FALSE, name="Courier New")
  
  fixed.style.low <- xlsx::CellStyle(wb) + 
    xlsx::Font(wb, heightInPoints = 10, isBold = FALSE, name="Courier New")+
    xlsx::Fill(backgroundColor = "orange",
               foregroundColor = "orange")
  
  fixed.style.missing <- xlsx::CellStyle(wb) + 
    xlsx::Font(wb, heightInPoints = 10, isBold = FALSE, name="Courier New")+
    xlsx::Fill(backgroundColor = "salmon",
               foregroundColor = "salmon")
  
  
  rows <- getRows(sh) 
  for(i in 2:length(xlsx::getRows(sh))){
    cells <- getCells(rows[i]) 
    cell.names <- names(cells)
    for(col in 3:15){
      cell.name <- cell.names[col]
      cell.ref <- cells[[cell.name]]
      zf.motif <- paste0("ZF_",(16-col)) # zfs are in descending order
      cell.val <-  contact.bases[i-1, col]
      
      
      if( is.na(cell.val)){
        xlsx::setCellStyle(cell.ref, fixed.style.missing)
      } else {
        
        common.zf.base <- zf.contact.bases.conserved[zf.contact.bases.conserved$motif_number==zf.motif, "contact_bases"]
        # print(paste("Col", col, "Row", i, "Value", cell.val, "Motif", zf.motif, "Common", common.zf.base))
        if(cell.val != common.zf.base){
          xlsx::setCellStyle(cell.ref, fixed.style.low)
        } else {
          xlsx::setCellStyle(cell.ref, fixed.style)
        }
      }
    }
  }
  
  xlsx::autoSizeColumn(sh, 1:ncol(contact.bases))
  xlsx::saveWorkbook(wb, file="figure/locations.zf.contact_bases.xlsx")
}
create.contact.base.xlsx(zf.contact.bases)

#### Export the residues within each NLS ####

LOCATIONS$combined.NLS %>%
  dplyr::select(Sequence = sequence, motif_number, aa_motif) %>%
  dplyr::mutate(motif_number = paste0("NLS_", motif_number)) %>%
  tidyr::pivot_wider(id_cols = Sequence, names_from = motif_number, values_from = aa_motif,
                     values_fn = ~paste(.x, collapse = ", ")) %>%
  as.data.frame %>%
  dplyr::arrange(as.integer(Sequence)) %>%
  create.xlsx(., "figure/locations.NLS.xlsx", cols.to.fixed.size.font = 2:4)

#### Export the residues within each 9aaTAD ####

# All unique 9aaTADs identified
LOCATIONS$combined.9aaTAD %>%
  dplyr::select(Sequence = sequence, label, aa_motif) %>%
  # summarise unique motifs
  dplyr::group_by(label, aa_motif) %>%
  dplyr::summarise(Count = n(),
                   Sequences = case_when(Count > 6 ~ "Others",
                                         .default = paste(Sequence, collapse = ", ")))%>%
  as.data.frame %>%
  create.xlsx(., "figure/locations.9aaTAD.unique.xlsx", cols.to.fixed.size.font = 2, cols.to.rich.text = 2)


# For each high-confidence 9aaTAD region:
# Extract the aa_motif if available, whatever the rc score
# If there is no detected 9aaTAD, extract the superTAD region sequence
extract.superTAD.motifs <- function(locations.9aaTAD, ranges.9aaTAD){
  
  # Look at all superTAD locations in all sequences
  sequences <- unique(locations.9aaTAD$sequence)
  tad.labels <- na.omit(unique(locations.9aaTAD$label))
  combos <- expand.grid("sequence" = sequences, "tad.label"= tad.labels)
  combos <- merge(combos, ranges.9aaTAD, by.x = "tad.label", by.y = "label",
                  all.x = TRUE)
  
  # Combine the superTAD locations for each row
  locations.superTAD <- merge(combos, locations.9aaTAD, 
                              by.x = c("sequence", "tad.label", "start",
                                       "end","width", "motif_number"),
                              by.y = c("sequence", "label", "start",
                                       "end", "width", "motif_number"), all.x = TRUE) %>%
    merge(., METADATA$combined, by.x=c("sequence"), by.y = c("common.name")) %>%
    
    # Reduce overlapping TAD ranges for each sequence
    dplyr::group_by(sequence, tad.label) %>%
    dplyr::mutate(tad.start = min(start_gapped),
                  tad.end = max(end_gapped),
                  tad.start = ifelse(is.na(tad.start), start, tad.start),
                  tad.end = ifelse(is.na(tad.end), end, tad.end)) %>%
    # EXtract the sequence for the tad range from the aa alignment
    dplyr::rowwise() %>%
    dplyr::mutate(tad.sequence = subset.sequence(ALIGNMENTS$aa.combined.biostrings,
                                                 sequence,
                                                 tad.start,
                                                 tad.end)) %>%
    # Drop unneeded columns
    dplyr::select(Sequence = sequence, tad.label, rc_score, tad.sequence, Group = group) %>%
    dplyr::mutate(Sequence = factor(Sequence, levels = str_replace_all(combined.taxa.name.order, " ", "_"))) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Sequence, tad.label) %>%
    dplyr::arrange(Sequence, tad.label, desc(rc_score)) %>%
    dplyr::slice_head(n=1) %>% # if there are multiple rc score per supertad, take only the highest
    dplyr::distinct() %>%
    
    # Pivot to each 9aaTAD as a separate column
    tidyr::pivot_wider(id_cols = c(Sequence, Group), names_from = tad.label, 
                       values_from = c(tad.sequence, rc_score), values_fn = ~paste(.x, collapse = ", "),
                       names_glue = "9aaTAD_{tad.label}_{.value}") %>%
    as.data.frame %>%
    dplyr::mutate(across(`9aaTAD_A_rc_score`:`9aaTAD_G_rc_score`, as.numeric)) %>%
    dplyr::rename_with(.fn = function(x) gsub("_tad.sequence", "", x)) %>%
    dplyr::arrange(Group, Sequence)
  
  # Create a custom xlsx colouring background of cells by rc score
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
      } else {
        # Set entire cell to normal style
        new.value <- rJava::.jnew("org/apache/poi/xssf/usermodel/XSSFRichTextString",
                                  oldval )
        rJava::.jcall(obj=new.value, returnSig = "V",  # void return
                      method="applyFont", normal.font.ref)
      }
    }
  }
  
  wb = xlsx::createWorkbook(type = "xlsx")
  sh = xlsx::createSheet(wb)
  xlsx::addDataFrame(locations.superTAD[,1:9], sh, row.names = F)
  xlsx::createFreezePane(sh, 2, 2, 2, 2) # freeze top row and first column
  
  fixed.style <- xlsx::CellStyle(wb) + 
    xlsx::Font(wb, heightInPoints = 10, isBold = FALSE, name="Courier New")
  
  fixed.style.low.rc <- xlsx::CellStyle(wb) + 
    xlsx::Font(wb, heightInPoints = 10, isBold = FALSE, name="Courier New")+
    xlsx::Fill(backgroundColor = "orange",
               foregroundColor = "orange")
  
  fixed.style.missing <- xlsx::CellStyle(wb) + 
    xlsx::Font(wb, heightInPoints = 10, isBold = FALSE, name="Courier New")+
    xlsx::Fill(backgroundColor = "salmon",
               foregroundColor = "salmon")
  
  
  rows <- getRows(sh) 
  
  set.cell.formatting <- function(sh, data, row.indexes){
    for(i in row.indexes){
      cells <- getCells(rows[i]) 
      cell.names <- names(cells)
      for(col in 3:9){
        cell.name <- cell.names[col]
        cell.ref <- cells[[cell.name]]
        rc.val <-  data[i-1, col+7]
        
        if( is.na(rc.val)){
          xlsx::setCellStyle(cell.ref, fixed.style.missing)
        } else {
          if(rc.val < 80){
            xlsx::setCellStyle(cell.ref, fixed.style.low.rc)
          } else {
            xlsx::setCellStyle(cell.ref, fixed.style)
          }
        }
      }
    }
  }
  
  set.cell.formatting(sh, locations.superTAD, 2:length(xlsx::getRows(sh)))
  
  
  for(col in 3:9) set.rich.text.on.vv(wb, sh, col)
  
  # Save the full table
  xlsx::autoSizeColumn(sh, 1:ncol(locations.superTAD))
  xlsx::saveWorkbook(wb, file="figure/locations.9aaTAD.xlsx")
  
  # Save only the rows needed for an in-text table of rodents and human
  table.values <- locations.superTAD %>% 
    dplyr::filter(Sequence %in% c("Human_ZFY", "Mouse_Zfy1", "Mouse_Zfy2",
                                  "African_Grass_Rat_ZFY2-like_1", "African_Grass_Rat_ZFY2-like_2",
                                  "Black_rat_Zfy2",
                                  "Norwegian_Rat_Zfy2",
                                  "Mongolian_gerbil_Zfx-like_putative-Zfy",
                                  "Desert_hamster_Zfx-like_putative-Zfy",
                                  "North_American_deer_mouse_Zfx-like_putative-Zfy",
                                  "Alpine_marmot_ZFY",
                                  "Arctic_ground_squirrel_Zfx-like_putative-Zfy",
                                  "Gray_squirrel_Zfy",
                                  "Beaver_Zfx-like_putative-Zfy",
                                  "Damara_mole-rat_Zfy"),
                  Group == "ZFY")
  
  wb = xlsx::createWorkbook(type = "xlsx")
  sh = xlsx::createSheet(wb)
  xlsx::addDataFrame(table.values[,1:9], sh, row.names = F)
  xlsx::createFreezePane(sh, 2, 2, 2, 2) # freeze top row and first column
  xlsx::autoSizeColumn(sh, 1:ncol(table.values))
  set.cell.formatting(sh, table.values, 2:length(xlsx::getRows(sh)))
  for(col in 3:9) set.rich.text.on.vv(wb, sh, col)
  xlsx::saveWorkbook(wb, file="figure/locations.9aaTAD.table.xlsx")
  
}

extract.superTAD.motifs(LOCATIONS$combined.9aaTAD, RANGES$combined.9aaTAD)

#### Export binding motifs of the ZFs in each species ####

# Read each file, get headers and PWMs
pwm.files <- list.files("aln/pwm", pattern = "zf_", full.names = T)

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
  tidyr::replace_na(list(zf_target_1="", zf_target_2="")) %>%
  dplyr::mutate(zf_target = paste0(zf_target_1, zf_target_2)) %>%
  dplyr::group_by(zf_target) %>%
  dplyr::summarise(total = n(), sequences = paste(sequence, collapse = ", ")) %>%
  dplyr::mutate(sequences = ifelse(total>40, "All others", sequences)) %>%
  as.data.frame


create.xlsx(pwm.predictions, "figure/zf_binding_targets.xlsx", cols.to.fixed.size.font = 1)

#### Ancestral sequence reconstruction #####

# Read the ancestral reconstruction
ancestral.nt.seqs <- read.table(paste0(FILES$mammal.nt.aln, ".state") ,header=TRUE)

# We care about the eutherian common ancestor and the rodent ancestor
# Find these nodes

# Read the tree back to keep full node names
nt.aln.tree.nodes <- ape::read.tree(paste0(FILES$mammal.nt.aln, ".treefile"))
nt.aln.tree.nodes.root <- ape::getMRCA(nt.aln.tree.nodes, c("Platypus_ZFX", "Australian_echidna_ZFX"))
nt.aln.tree.nodes <- ape::root(nt.aln.tree.nodes, node=nt.aln.tree.nodes.root)

# Node number returned includes number of tip labels; subtract to get node
rodent.node <- ape::getMRCA(nt.aln.tree.nodes, c("Mouse_Zfy1", "Desert_hamster_Zfx-like_putative-Zfy")) - length(nt.aln.tree.nodes$tip.label)
eutheria.node <- getMRCA(nt.aln.tree.nodes, c("Mouse_Zfy1", "African_bush_elephant_ZFY")) - length(nt.aln.tree.nodes$tip.label)

# Get ancestral rodent sequence
rodent.anc.nt <- ancestral.nt.seqs %>%
  dplyr::filter(Node == paste0("Node", rodent.node)) %>%
  dplyr::select(State)
rodent.anc.nt <- ape::as.DNAbin(rodent.anc.nt$State)

# Bind together as a matrix
rodent.plus.anc.nt.aln <- rbind(as.matrix(ALIGNMENTS$nt.mammal.ape), rodent.anc.nt)

# Convert back to list and add names
rodent.ancestor.label <- "Muroidea ancestral Zfy"
rodent.plus.anc.nt.aln <- as.list.DNAbin(rodent.plus.anc.nt.aln)
names(rodent.plus.anc.nt.aln) <- c(names(ALIGNMENTS$nt.mammal.ape), rodent.ancestor.label)

# Plot the MSA


mammal.taxa.name.order <- gsub(" ", "_", mammal.taxa.name.order)
combined.taxa.name.order <- gsub(" ", "_", combined.taxa.name.order)
# Find the tips under the ancestral node and get the labels
rodent.tip.labels <- tidytree::offspring( as_tibble(nt.aln.tree.nodes), .node = rodent.node+length(nt.aln.tree.nodes$tip.label), tiponly = T)


rodent.plus.anc.nt.aln.tidy <- tidy_msa(rodent.plus.anc.nt.aln) %>%
  dplyr::filter(name %in% c(rodent.ancestor.label, rodent.tip.labels$label))
rodent.plus.anc.nt.aln.tidy$name <- factor(rodent.plus.anc.nt.aln.tidy$name, 
                                           levels = rev( c(mammal.taxa.name.order, rodent.ancestor.label))) # sort reverse to match tree

# Filter the zf locations and correct y locations
rodent.plus.anc.zf <- LOCATIONS$combined.zf %>%
  dplyr::filter(sequence %in% rodent.plus.anc.nt.aln.tidy$name) %>%
  dplyr::rowwise() %>%
  # Correct for different number of taxa in the y axis
  dplyr::mutate(i = which( levels(rodent.plus.anc.nt.aln.tidy$name) == sequence ) - (1+length(unique(mammal.taxa.name.order))- length(unique(rodent.plus.anc.nt.aln.tidy$name))))

# Filter 9aaTAD locations and correct the y locations
rodent.plus.anc.9aaTAD <- LOCATIONS$combined.9aaTAD %>%
  dplyr::filter(sequence %in% rodent.plus.anc.nt.aln.tidy$name) %>%
  dplyr::rowwise() %>%
  # Correct for different number of taxa in the y axis
  dplyr::mutate(i = which( levels(rodent.plus.anc.nt.aln.tidy$name) == sequence ) - (1+length(unique(mammal.taxa.name.order))- length(unique(rodent.plus.anc.nt.aln.tidy$name))))

rodent.plus.anc.NLS <- LOCATIONS$combined.NLS %>%
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

if(file.exists(FILES$paml.branch.site.output)){
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

# Use Xenopus ZFX.S as the comparison group
msa.aa.aln.tidy.frog.conservation    <- calculate.conservation(ALIGNMENTS$aa.combined.biostrings,"Xenopus_ZFX.S" )
msa.aa.aln.tidy.chicken.conservation <- calculate.conservation(ALIGNMENTS$aa.combined.biostrings,"Chicken_ZFX" )
msa.aa.aln.tidy.opossum.conservation <- calculate.conservation(ALIGNMENTS$aa.combined.biostrings,"Opossum_ZFX" )

#### Create overall structure plot ####

# Combine the structural conservation plot with the aa tree
# This should show all the 9aaTADs
aa.structure.plot <- ggplot()+
  geom_tile(data = LOCATIONS$combined.zf,     aes(x=(start_gapped+end_gapped)/2,
                                         width = (end_gapped-start_gapped),
                                         y=sequence),
            fill="grey", alpha=0.5)+
  geom_tile(data = LOCATIONS$combined.9aaTAD,     aes(x=(start_gapped+end_gapped)/2,
                                             width = (end_gapped-start_gapped),
                                             y=sequence,
                                             fill=round(rc_score)),
            alpha=0.9)+
  paletteer::scale_fill_paletteer_c("grDevices::Blues 3", "9aaTAD RC score (%)", direction = -1, limits = c(50, 100)) +
  geom_tile(data = LOCATIONS$combined.NLS,     aes(x=(start_gapped+end_gapped)/2,
                                          width = (end_gapped-start_gapped),
                                          y=sequence),
            fill=NLS.COLOUR, alpha=0.5)
aa.structure.plot <- annotate.structure.plot(aa.structure.plot, length(combined.taxa.name.order) + 1.5)

save.double.width("figure/Figure_Sxxxx_aa.structure.png", aa.structure.plot, height = 120)

# Also create a trimmed down version that has only the 83, 92, 100% confidence 9aaTADs
aa.structure.confident.plot <- ggplot()+
  geom_tile(data = LOCATIONS$combined.zf,     aes(x=(start_gapped+end_gapped)/2,
                                         width = (end_gapped-start_gapped),
                                         y=sequence),
            fill="grey", alpha=0.5)+
  geom_tile(data = LOCATIONS$combined.9aaTAD[LOCATIONS$combined.9aaTAD$rc_score>80,],     aes(x=(start_gapped+end_gapped)/2,
                                                                            width = (end_gapped-start_gapped),
                                                                            y=sequence,
                                                                            fill=round(rc_score)),
            alpha=0.9)+
  paletteer::scale_fill_paletteer_c("grDevices::Blues 3", "9aaTAD RC score (%)", direction = -1, limits = c(50, 100)) +
  geom_tile(data = LOCATIONS$combined.NLS,     aes(x=(start_gapped+end_gapped)/2,
                                          width = (end_gapped-start_gapped),
                                          y=sequence),
            fill="green", alpha=0.5)
aa.structure.confident.plot <- annotate.structure.plot(aa.structure.confident.plot, length(combined.taxa.name.order) + 1.5)

save.double.width("figure/Figure_2_aa.structure.confident.png", aa.structure.confident.plot, height = 120)


#### Calculate conservation between X and Y gametologues for each species ####

# Given species names, find the X and Y sequences from the aa alignment.
# Generate a vector of 0 (diferent aa at a site) or 1 (same aa at a site).
calculate.gametologue.conservation <- function(zfx, zfy, common.name){
  
  # cat(zfx, "\n")
  tryCatch({
    # Find the characters in the reference sequence
    aa.aln <- ggmsa::tidy_msa(ALIGNMENTS$aa.combined.biostrings)  %>%
      dplyr::filter( name==zfx | name==zfy) %>%
      dplyr::mutate(type = case_when(str_detect(name, "(?i)zfy") ~ "ZFY",
                                     .default = "ZFX")) %>%
      tidyr::pivot_wider(id_cols = position, names_from = type, values_from = character) %>%
      dplyr::mutate(conservation = as.numeric(ZFX==ZFY))
    
    return(data.frame("zfx"=zfx, "zfy"=zfy, "common.name" = common.name,
                "conservation" = aa.aln$conservation, site = 1:length( aa.aln$conservation)))
  },
  error = function(e){
    message(conditionMessage(e))
    return(data.frame("zfx"=zfx, "zfy"=zfy, "common.name" = common.name, "conservation"=NA, site=NA ))
  }
  )
  
}

# Get only one ZFX/Y pair per species
valid.pairs <- METADATA$mammal %>%
  dplyr::filter(group=="ZFX" | group=="ZFY") %>%
  dplyr::group_by(species) %>%
  dplyr::arrange(species, group) %>%
  dplyr::slice_head(n=2) %>% # skip second zfy if present
  tidyr::pivot_wider(id_cols = Species_common_name, names_from = group, values_from = common.name)

# Calculate homologue conservation across the alignments
gametologue.conservation <- do.call(rbind, mapply(calculate.gametologue.conservation, 
                                                  valid.pairs$ZFX,   
                                                  valid.pairs$ZFY, 
                                                  valid.pairs$Species_common_name, 
                                                  SIMPLIFY = FALSE)) %>%
  dplyr::mutate(zfx = factor(zfx, levels = combined.taxa.name.order),
                zfy = factor(zfy, levels = combined.taxa.name.order),
                y.val = as.integer(zfy)) %>%
  dplyr::group_by(y.val) %>%
  dplyr::mutate( cum.diff =  cumsum(conservation==0),
                 group = case_when(common.name=="Mouse" ~ "Eumuroida",
                                   common.name=="African_Grass_Rat" ~ "Eumuroida",
                                   common.name=="Black_rat" ~ "Eumuroida",
                                   common.name=="Norwegian_Rat" ~ "Eumuroida",
                                   common.name=="Mongolian_gerbil" ~ "Eumuroida",
                                   common.name=="Desert_hamster" ~ "Eumuroida",
                                   common.name=="North_American_deer_mouse" ~ "Eumuroida",
                                   common.name=="Damara_mole-rat" ~ "Damara mole-rat & beaver",
                                   common.name=="Beaver" ~ "Damara mole-rat & beaver",
                                   common.name=="Grey_squirrel" ~ "Other Rodentia",
                                   common.name=="Alpine_marmot" ~ "Other Rodentia",
                                   common.name=="Arctic_ground_squirrel" ~ "Other Rodentia",
                                   .default = "Other Mammalia"),
                 group = forcats::as_factor(group),
                 group = forcats::fct_relevel(group, "Eumuroida", "Damara mole-rat & beaver", 
                                              "Other Rodentia", "Other Mammalia"))


# Create a cumulative difference plot
gametologue.cumdiff.plot <-  ggplot()
gametologue.cumdiff.plot <- add.structures(gametologue.cumdiff.plot, y.start= 0, y.end = 270, alpha=0.5)
gametologue.cumdiff.plot <- gametologue.cumdiff.plot+
  geom_line(data = gametologue.conservation, aes(x=site, y=cum.diff, col=group, group = common.name ))+
  geom_line(data = gametologue.conservation[gametologue.conservation$group=="Other Rodentia",], 
            aes(x=site, y=cum.diff, col=group, group = common.name ))+
  labs(y = "Cumulative X-Y differences", x = "Site in alignment")+
  scale_color_manual(values = c("Eumuroida"="#002AFFFF", "Other Mammalia"="#bbbbbbFF", 
                                "Other Rodentia"="#FF6619FF", "Damara mole-rat & beaver"="darkgreen"))+
  scale_x_continuous(breaks = seq(0, 900, 100))+
  scale_y_continuous(breaks = seq(0, 300, 50))+
  theme_bw()+
  theme(legend.position = c(0.2, 0.7),
        legend.title = element_blank())

gametologue.cumdiff.plot <- gametologue.cumdiff.plot +
new_scale_fill()+
  scale_fill_manual(values=c("white", "grey", "white", "grey", "white", "grey", "white"))+
  scale_pattern_color_manual(values=c("white", "white"))+
  scale_pattern_manual(values = c("none", "stripe")) + # which exons are patterned
  guides(fill = "none", pattern="none")+
  add.exon.track(270, 285, col = "black")+ # color of border
  add.exon.labels(270, 285)


  
save.double.width(filename = "figure/gametologue.cumdiff.png", gametologue.cumdiff.plot, height = 100)

#### Calculate conservation between X and Y gametologues for ancestral nodes species ####
# # Read the ancestral states for the nodes
ancestral.zfx.seqs <- read.table("aln/zfx_only/zfx.indel.filtered.state", header=TRUE) %>%
  dplyr::mutate(Type = "ZFX") %>%
  dplyr::select(Node, Type, Site, State = State.x)

ancestral.zfy.seqs <- read.table("aln/zfy_only/zfy.indel.filtered.state", header=TRUE) %>%
  dplyr::mutate(Type = "ZFY") %>%
  dplyr::select(Node, Type, Site, State = State.x)


ancestral.seqs <- rbind(ancestral.zfx.seqs, ancestral.zfy.seqs) %>%
  dplyr::filter(Node != "Arvicanthis", Node != "Mus", Node !="Monotremata", 
                Node != "Mammalia", Node !="Marsupialia")

calculate.ancestral.conservation <- function(node.name){
  ancestral.seqs %>% 
    dplyr::filter(Node==node.name) %>%
    tidyr::pivot_wider(id_cols = c(Node, Site), names_from = Type, values_from = State) %>%
    dplyr::arrange(Node, Site) %>%
    dplyr::mutate(Conservation = as.numeric(ZFX==ZFY),
                  CumSum = cumsum(Conservation==0))
}

ancestral.conservation <- do.call(rbind, lapply(unique(ancestral.seqs$Node), calculate.ancestral.conservation))

# Reorder the nodes that will be highlighted so legend ordering is sensible
ancestral.conservation$Node <- as.factor(ancestral.conservation$Node)
ancestral.conservation$Node <- forcats::fct_relevel(ancestral.conservation$Node, "Murinae","Eumuroida", "Muroidea","Eutheria","Theria", after = 0)

highlight.colours <- RColorBrewer::brewer.pal(5, "Dark2")

ancestral.cumdiff.plot <-  ggplot()
ancestral.cumdiff.plot <- add.structures(ancestral.cumdiff.plot, y.start= 0, y.end = 250, alpha=0.5, 
                                         start_col="start_aa_mammal", end_col = "end_aa_mammal")
ancestral.cumdiff.plot <-  ancestral.cumdiff.plot+
  geom_line(data = ancestral.conservation, aes(x=Site, y=CumSum, group=Node), col="#bbbbbbFF")+
  geom_line(data = ancestral.conservation[ancestral.conservation$Node %in% c("Muroidea","Eumuroida","Theria","Murinae","Eutheria" ),], aes(x=Site, y=CumSum, col=Node))+
  labs(y = "Cumulative X-Y differences", x = "Site in alignment")+
  scale_color_manual(values = c("Theria"=highlight.colours[1],"Eutheria"=highlight.colours[2], 
                                "Muroidea"=highlight.colours[3],"Eumuroida"=highlight.colours[4],
                                "Murinae"=highlight.colours[5]))+
  # scale_color_manual(values = c("Muroidea"="purple","Eumuroida"="#002AFFFF","Theria"="darkgreen", "Murinae"="orange", "Eutheria"="red"))+
  scale_x_continuous(breaks = seq(0, 900, 100))+
  scale_y_continuous(breaks = seq(0, 900, 50))+
  theme_bw()+
  theme(legend.position = c(0.15, 0.68),
        legend.title = element_blank())

ancestral.cumdiff.plot <- ancestral.cumdiff.plot +
  new_scale_fill()+
  scale_fill_manual(values=c("white", "grey", "white", "grey", "white", "grey", "white"))+
  scale_pattern_color_manual(values=c("white", "white"))+
  scale_pattern_manual(values = c("none", "stripe")) + # which exons are patterned
  guides(fill = "none", pattern="none")+
  add.exon.track(250, 270, col = "black", start_col = "start_aa_mammal", end_col = "end_aa_mammal")+ # color of border
  add.exon.labels(250, 270, start_col = "start_aa_mammal", end_col = "end_aa_mammal")

save.double.width(filename = "figure/ancestral.cumdiff.png", ancestral.cumdiff.plot, height = 100)

#### Plot the conservation of hydrophobicity across mammal/outgroup AA MSA ####
cat("Plotting conservation of hydrophobicity\n")
msa.aa.aln.hydrophobicity <- do.call(rbind, mapply(calc.hydrophobicity, aa=ALIGNMENTS$aa.combined.biostrings@unmasked, 
                                                   sequence.name = names(ALIGNMENTS$aa.combined.biostrings@unmasked), 
                                                   window.size = 9,
                                                   SIMPLIFY = FALSE))  %>%
  dplyr::mutate(sequence = factor(sequence, levels = rev(combined.taxa.name.order))) # sort reverse to match tree

n.taxa <- length(combined.taxa.name.order) + 1.5


# TODO: guides does nto work with new_scale(), so we may need to specify the order of the colour bar manually
# scale_fill_gradient2(low = "#f0f0f0", mid = "#969696", high = "#252525", midpoint = 21, n.breaks = 5, 
#                      guide = guide_colourbar(title = "Salinity", title.theme = element_text(size = 16, face = "bold"),
#                                              label.theme = element_text(size = 14, face = "bold"), label.position = "top",
#                                              barwidth = 10.7, barheight = 1, order = 1, frame.linetype = 1, frame.colour = "#000000",
#                                              ticks.colour = "#000000", direction = "horizontal")))

hydrophobicity.plot <- ggplot()+
  geom_raster(data=msa.aa.aln.hydrophobicity,  aes(x = position_gapped,
                                                   y = sequence,
                                                   fill=hydrophobicity_smoothed),
                                                   hjust = 0, vjust = 0)+ 
  scale_fill_paletteer_c("ggthemes::Classic Red-Black", direction = -1, limits = c(0, 1))+
  labs(fill="Hydrophobicity (9 site average)")
hydrophobicity.plot <- annotate.structure.plot(hydrophobicity.plot, n.taxa)
save.double.width("figure/hydrophobicity.convervation.tree.png", hydrophobicity.plot, height = 120)

#### Plot conservation of charge across mammal/outgroup AA MSA ####
cat("Plotting conservation of charge\n")
msa.aa.aln.charge <- do.call(rbind, mapply(calc.charge, aa=ALIGNMENTS$aa.combined.biostrings@unmasked, 
                                           sequence.name = names(ALIGNMENTS$aa.combined.biostrings@unmasked), 
                                           window.size = 9,
                                           SIMPLIFY = FALSE))  %>%
  dplyr::mutate(sequence = factor(sequence, levels = rev(combined.taxa.name.order))) # sort reverse to match tree

n.taxa <- length(combined.taxa.name.order) + 1.5
charge.plot <- ggplot()+
  geom_raster(data=msa.aa.aln.charge,  aes(x = position_gapped,
                                                   y = sequence,
                                                   fill=charge_smoothed),
              hjust = 0, vjust = 0)+ 
  scale_fill_paletteer_c("ggthemes::Red-Blue Diverging", direction = 1, limits = c(-1, 1))+
  labs(fill="Charge (9 site average)")
charge.plot <- annotate.structure.plot(charge.plot, n.taxa)
save.double.width("figure/charge.convervation.tree.png", charge.plot, height = 120)

#### Combine charge and hydorphobicity structure plots ####
cat("Plotting combined hydrophobicity and charge\n")
# To test spacing and balance
structure.plot <- (charge.plot / hydrophobicity.plot) +
  plot_layout(ncol = 1)+ 
  plot_annotation(tag_levels = list(c("A", "", "B", "")))
save.double.width("figure/Figure_xxxx_combined_structure.plot.png", structure.plot, height = 230)


#### Plot the hydrophobic patches in exons 2, 3 and 5 that are not 9aaTADs ####
cat("Identifying hydrophobic patches\n")

plot.hydrophobic.patch <- function(patch.data, patch.start, patch.end, patch.consensus){
  
  y.max <- length(unique(patch.data$sequence)) + 2 # scale plot height to number of sequences
  y.consensus.label <- y.max - 0.8 # the consensus logo should be at the top
  
  ggplot()+
    geom_tile(data=patch.data,  aes(x = position_gapped, y = sequence, fill=hydrophobicity))+
    scale_fill_paletteer_c("ggthemes::Classic Red-Black", direction = -1, limits = c(0, 1))+
    labs(fill="Hydrophobicity (per residue)")+
    geom_text(data = patch.data, aes(x = position_gapped, y = sequence, label=character), size=2, family="mono", col="white")+
    scale_y_discrete(labels = gsub("_", " ", rev(combined.taxa.name.order)))+
    coord_cartesian(xlim = c(patch.start,patch.end),  ylim = c(0, y.max), clip = 'on')+
    annotate("text", label=s2c(patch.consensus), size=2, x=(patch.start):(patch.start+nchar(patch.consensus)-1), 
             y=y.consensus.label, hjust=0.5, family="mono", fontface="bold")
}

# TODO - coordinates are hard coded from the alignments. Will need updating if more
# sequences are added
exon2.patch <- extract.combined.alignment.region(aa.start = mouse.exons$start_aa_combined[2]+93,
                                                 aa.end   = mouse.exons$start_aa_combined[2]+108)

exon3.patch <- extract.combined.alignment.region(aa.start = mouse.exons$start_aa_combined[3]+28,
                                        aa.end   = mouse.exons$end_aa_combined[3]-2)

exon5.patch <- extract.combined.alignment.region(aa.start = mouse.exons$start_aa_combined[5]+27,
                                        aa.end   = mouse.exons$end_aa_combined[5]+1)

exon7.patch <- extract.combined.alignment.region(aa.start = mouse.exons$start_aa_combined[7]-2,
                                                 aa.end   = mouse.exons$start_aa_combined[7]+22)

# Export patches to file
merge(exon2.patch$aa.aln, exon3.patch$aa.aln, by=c("Sequence")) %>%
  merge(., exon5.patch$aa.aln, by=c("Sequence")) %>%
  merge(., exon7.patch$aa.aln, by=c("Sequence")) %>%
  create.xlsx(., "figure/hydrophobic_patches.xlsx")

# Plot the patch MSAs
exon2.hydro.plot <- plot.hydrophobic.patch(msa.aa.aln.hydrophobicity, exon2.patch$aa.start, exon2.patch$aa.end, exon2.patch$aa.consensus)
exon3.hydro.plot <- plot.hydrophobic.patch(msa.aa.aln.hydrophobicity, exon3.patch$aa.start, exon3.patch$aa.end, exon3.patch$aa.consensus)
exon5.hydro.plot <- plot.hydrophobic.patch(msa.aa.aln.hydrophobicity, exon5.patch$aa.start, exon5.patch$aa.end, exon5.patch$aa.consensus)
exon7.hydro.plot <- plot.hydrophobic.patch(msa.aa.aln.hydrophobicity, exon7.patch$aa.start, exon7.patch$aa.end, exon7.patch$aa.consensus)

# Collect the plots

patch.plot.complete <- exon2.hydro.plot + exon3.hydro.plot + exon5.hydro.plot +
  plot_layout(nrow=1, guides = "collect", axes = "collect") &
  theme_bw()+
  theme(legend.position = "top",
        axis.text = element_text(size=6),
        axis.title = element_blank(),
        legend.text = element_text(size=6),
        legend.title = element_text(size=8))
save.double.width(filename = "figure/hydrophobic.patch.all.png", patch.plot.complete)

#### Plot HyPhy RELAX test for relaxed selection ####
cat("Plotting RELAX result\n")
create.relax.k.tree <- function(json.file){
  # Read the json file and parse results
  hyphy.data <- jsonlite::read_json(json.file)
  # Note that the tree needs to terminate with ; otherwise read.tree returns NULL
  hyphy.input.tree <- ape::read.tree(text = paste0(hyphy.data$input$trees[["0"]], ";"))
  
  # Coloration of tree by k based on https://observablehq.com/@spond/plotting-relax-k-values-on-branches-of-the-tree
  
  # Make a dataframe with the k values and node numbers
  k.vals <- data.frame("k" = sapply(hyphy.data$`branch attributes`[["0"]], \(x) x$`k (general descriptive)`))
  k.vals$NodeLab <- rownames(k.vals)
  k.vals$node <- sapply(k.vals$NodeLab,  treeio::nodeid, tree = hyphy.input.tree)
  # Rescale values above 1 to the range 1-2 so we get a clean diverging scale
  # k.vals$adj.k <- ifelse(k.vals$k <= 1, k.vals$k, (k.vals$k/50)+1)
  k.vals$adj.k <- log(k.vals$k)
  
  # Add a new row for the root node
  k.vals[nrow(k.vals)+1,] <- list(0, "", length(hyphy.input.tree$tip.label)+1, 1)
  
  # Get the branch lengths from the HyPhy output
  k.vals$branch.length <-  sapply(k.vals$NodeLab, \(x) ifelse(x=="", 0, hyphy.data$`branch attributes`[["0"]][[x]]$`MG94xREV with separate rates for branch sets`))
  
  # Reorder the branches to match the node/tip order of the tree
  # Ordering in hyphy.input.tree$edge[,2] (the destination node, lengths are for incoming branches)
  branch.lengths <- unlist(sapply(hyphy.input.tree$edge[,2], \(x)  k.vals[k.vals$node==x,"branch.length"]))
  hyphy.input.tree$edge.length <- branch.lengths
  
  hyphy.input.tree <- reroot.tree(hyphy.input.tree, node.labels=c("Platypus_ZFX", "Australian_echidna_ZFX"),  position=0.015)
  hyphy.input.tree <- reroot.tree(hyphy.input.tree, node.labels= c("Xenopus_ZFX_L", "Xenopus_ZFX_S"),  position = 0.015) # will fail for mammal-only

  cat(json.file,  paste("K =", round(hyphy.data$`test results`$`relaxation or intensification parameter`, digits = 2), 
      "p =", round(hyphy.data$`test results`$`p-value`, digits = 2), collapse = " "), "\n", sep = "\t")
  
  p <- ggtree(hyphy.input.tree, size=1) + 
    geom_tiplab(size=2)
  p <- p %<+% k.vals + aes(colour=adj.k) + 
    scale_color_paletteer_c("ggthemes::Classic Red-Blue",
                            direction = 1,
                            limits = c(-3, 3))+
    labs(color = "log(K)")+
    coord_cartesian(xlim = c(0, 0.6), clip = "off")+
    geom_treescale(fontsize =2, y = -1, width = 0.05) +
    theme(legend.position = c(0.8, 0.5),
          legend.background = element_blank(),
          legend.text = element_text(size=6),
          legend.title = element_text(size=6),
          plot.margin = margin(r=0, l=0))
  
  p
}

if(file.exists("aln/hyphy/mammal.eumuroida.relax.json")){
  
  # node.names <- c("rodentia", "eumuroida", "muridae", "murinae")
  
  make.tree <- \(name){
    relax.combined.tree <- create.relax.k.tree(paste0("aln/hyphy/combined.", name, ".relax.json"))
    save.double.width(paste0("figure/RELAX.combined.", name, ".png"), relax.combined.tree)
    
    relax.mammal.tree <- create.relax.k.tree(paste0("aln/hyphy/mammal.", name, ".relax.json"))
    save.double.width(paste0("figure/RELAX.mammal.", name, ".png"), relax.mammal.tree)
  } 
  
  # sapply(node.names, make.tree)
  make.tree("eumuroida")
}

#### Plot HyPhy MEME test for directional selection ####

cat("Plotting MEME result\n")

# Read the json file of meme results for a given node
# Filter to those of interest at p<0.05
# outgroup.type: combined or mammal
read.meme.results <- function(node.name="rodentia", outgroup.type="combined"){
  
  json.file <- paste0("aln/hyphy/", outgroup.type, ".", node.name, ".meme.json")
  meme.data <- jsonlite::read_json(json.file)
  
  # Read a given partition
  read.partition <- function(partition){
    
    # Extract the values for the given partition
    meme.cols <- as_tibble(do.call(rbind, meme.data$MLE$content[[as.name(partition)]])) %>% 
      unnest_longer(everything())

    # Extract the sites covered by the partition. Correct for sites are 0-indexed
    meme.cols$site <- unlist(meme.data$`data partitions`[[as.name(partition)]]$coverage)+1
    meme.cols$partition <- partition
    
    # Extract column names
    colnames(meme.cols) <- c(unlist(lapply(meme.data$MLE$headers, \(i) i[[1]])), "site", "partition")
    meme.cols
  }
  
  # Read all partitions and combine
  meme.data <-  do.call(rbind, lapply(0:2, read.partition))
  
  meme.data$node <- node.name
  meme.data$outgroup.type <- outgroup.type
  
  # Plot the MEME results - window around each site
  meme.sites <- meme.data %>%
    dplyr::filter(`p-value`<0.0025) %>%
    dplyr::mutate(start = site-1,
                  end = site+1)
  
  list(meme.data = meme.data,
       meme.sites = meme.sites)
}

combined.meme.results <- lapply( c("eumuroida") , read.meme.results, outgroup.type="combined") # node.names

mammal.meme.results <- lapply(c("eumuroida") , read.meme.results, outgroup.type="mammal") #node.names

# Summarise the MEME results to give an overview of the significant sites at each test branch
calculate.meme.overview.data <- function(meme.results){
  # Find the unique set of sites in the meme data
  do.call(rbind, lapply(meme.results, \(x) x$meme.data)) %>%
    tidyr::pivot_wider(id_cols = c(site, partition), names_from = node, values_from = `p-value`) %>%
    dplyr::filter(if_any(.cols= eumuroida, ~.<0.01)) # filter to sites significant in any node
  
}

# Plot MEME overview data
# meme.overview - the output of calculate.meme.overview.data
# outgroup.type - combined or mammal, so features are plotted wrt the correct alignment
create.meme.overview.plot <- function(meme.overview, outgroup.type){
  
  if(outgroup.type=="mammal"){
    ranges.NLS.common <- RANGES$mammal.nls
    ranges.ZF.common <- RANGES$mammal.zf
    ranges.9aaTAD.common <- RANGES$mammal.9aaTAD
    aa.aln <- ALIGNMENTS$aa.mammal.ape
  } else {
    ranges.NLS.common <- RANGES$combined.nls
    ranges.ZF.common <- RANGES$combined.zf
    ranges.9aaTAD.common <- RANGES$combined.9aaTAD
    aa.aln <- ALIGNMENTS$aa.combined.ape
  }
  
  p.value.threshold <- 0.01
  
  # Annotate MEME sites over the exons
  meme.location.plot <- ggplot()+
    # Draw the structures
    add.track(ranges.NLS.common,    1.95, 2.55, fill=NLS.COLOUR, alpha = 1)+ # +8.5
    add.track(ranges.ZF.common,     2, 2.5, fill=ZF.COLOUR)+ # 9 - 11
    add.track.labels(ranges.ZF.common, 2, 2.5, col="white", label_col = "motif_number")+   # Label the ZFs
    add.track(ranges.9aaTAD.common, 2, 2.5, fill=TAD.COLOUR,  alpha = 0.9)+  #9
    add.track.labels(ranges.9aaTAD.common, 2, 2.5, col="white")+   # Label the 9aaTADs
    
    # annotate("rect", xmin = exon2.patch$aa.start, xmax = exon2.patch$aa.end, 
    #          ymin = 2, ymax = 2.5, fill = "pink", alpha=0.5)+
    # annotate("rect", xmin = exon3.patch$aa.start, xmax = exon3.patch$aa.end, 
    #          ymin = 2, ymax = 2.5, fill = "pink", alpha=0.5)+
    # annotate("rect", xmin = exon5.patch$aa.start, xmax = exon5.patch$aa.end, 
    #          ymin = 2, ymax = 2.5, fill = "pink", alpha=0.5)+
    
    new_scale_fill()+
    scale_fill_manual(values=c("white", "grey", "white", "grey", "white", "grey", "white"))+
    scale_pattern_color_manual(values=c("white", "white"))+
    
    scale_pattern_manual(values = c("none", "stripe")) + # which exons are patterned
    guides(fill = "none", pattern="none")+
    add.exon.track(y.start = 2.6, y.end = 3, col="black", start_col = "start_aa_mammal", end_col = "end_aa_mammal")+
    add.exon.labels(y.start = 2.6, y.end = 3, start_col = "start_aa_mammal", end_col = "end_aa_mammal")+
    
    # Show which test clades the sites show up in
    # geom_point(data = meme.overview, aes(x=site, y = 3.1, alpha = rodentia<p.value.threshold), col="black",size = 0.2)+
    geom_point(data = meme.overview, aes(x=site, y = 3.2, col = eumuroida), size = 0.2)+
    # geom_point(data = meme.overview, aes(x=site, y = 3.3, alpha = muridae<p.value.threshold), col="black",size = 0.2)+
    # geom_point(data = meme.overview, aes(x=site, y = 3.4, alpha = murinae<p.value.threshold), col="black", size = 0.2)+
    # scale_alpha_manual(values=c(0, 1))+
    # scale_color_viridis_c()+
    scale_color_paletteer_c("ggthemes::Classic Red-Black", direction = 1, limits = c(0, p.value.threshold))+
    
    geom_tile(data = meme.overview, aes(x=site, y = 2.5, height = 1, width=1), fill = "brown")+
    # annotate("text", x= 800, y=3.1, label = "rodentia", hjust=0, size = 2)+
    annotate("text", x= 790, y=3.2, label = "Eumuroida p-value", hjust=0, size = 2)+
    # annotate("text", x= 800, y=3.3, label = "muridae", hjust=0, size = 2)+
    # annotate("text", x= 800, y=3.4, label = "murinae", hjust=0, size = 2)+
    
    
    scale_x_continuous(breaks = seq(0, length(aa.aln[[1]]), 50), expand = c(0, 0.05))+
    coord_cartesian(xlim = c(0, length(aa.aln[[1]])))+
    theme_bw()+
    theme(axis.text.y = element_blank(),
          axis.title = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size=6),
          legend.position = "none",
          legend.title = element_text(size = 6, vjust = 0.85),
          legend.text = element_text(size = 6),
          legend.key.height = unit(3, "mm"),
          legend.spacing.y = unit(2, "mm"),
          legend.box.spacing = unit(2, "mm"),
          panel.border = element_blank(),
          axis.line.x.bottom = element_line(),
          panel.grid = element_blank())
  
  save.double.width(filename = paste0("figure/meme.",outgroup.type,".site.locations.png"), meme.location.plot, height = 35)
  meme.location.plot
}

mammal.meme.overview <- calculate.meme.overview.data(mammal.meme.results)
mammal.meme.overview.plot <- create.meme.overview.plot(mammal.meme.overview, outgroup.type="mammal")

# combined.meme.overview <- calculate.meme.overview.data(combined.meme.results)
# combined.meme.overview.plot <-create.meme.overview.plot(combined.meme.overview,"combined")

# both.meme.overview.plot <- combined.meme.overview.plot / mammal.meme.overview.plot
# save.double.width(filename = paste0("figure/meme.all.site.locations.png"), both.meme.overview.plot, height = 90)

# Given an alignment region produced by `extract.alignment.region`, plot it
plot.alignment.region <- function(region.data, meme.overview){
  
  # Convert aa site to nt
  meme.overview$site.nt <- meme.overview$site*3-1
  
  tidy.region.nt.data <- region.data$nt.aln %>%
    tidyr::pivot_longer( -Sequence, names_to = "Site", values_to = "Base") %>%
    dplyr::mutate(Site = as.integer(gsub("Site_", "", Site)))
  
  tidy.region.aa.data <- region.data$aa.aln %>%
    tidyr::pivot_longer( -Sequence, names_to = "Site", values_to = "Residue") %>%
    dplyr::mutate(Site = as.integer(gsub("Site_", "", Site))*3-1) %>%
    dplyr::mutate(isSignificant = Site %in% (meme.overview$site*3-1)) %>%
    merge(., meme.overview, by.x = "Site", by.y = "site.nt", all.x=T) %>%
    dplyr::mutate(sequenceInt = as.integer(Sequence),
                  clade = case_when(str_detect(Sequence, "Grass_Rat_ZFY") ~ "rodentia eumuroida muridae murinae",
                                    str_detect(Sequence, "Mouse_Zfy") ~ "rodentia eumuroida muridae murinae",
                                    str_detect(Sequence, "Black_rat_Zfy") ~ "rodentia eumuroida muridae murinae",
                                    str_detect(Sequence, "Norwegian_Rat_Zfy") ~ "rodentia eumuroida muridae murinae",
                                    str_detect(Sequence, "Mongolian_gerbil_Zfx-like") ~ "rodentia eumuroida muridae",
                                    str_detect(Sequence, "deer_mouse_Zfx-like") ~ "rodentia eumuroida",
                                    str_detect(Sequence, "hamster_Zfx-like") ~ "rodentia eumuroida",
                                    str_detect(Sequence, "marmot_ZFY") ~ "rodentia",
                                    str_detect(Sequence, "squirrel_Zfx-like") ~ "rodentia",
                                    str_detect(Sequence, "squirrel_Zfy") ~ "rodentia",
                                    str_detect(Sequence, "Damara_mole_rat_Zfy") ~ "rodentia",
                                    str_detect(Sequence, "Beaver_Zfx-like") ~ "rodentia",
                                    .default = "other"
                  ),
                  isSignificantSite = case_when( str_detect(clade, "eumuroida") & eumuroida<0.01 ~ T,
                                                 
                                                 .default = F)
                  
    )
  
  
  
  nt.consensus.y <- max(as.integer(tidy.region.nt.data$Sequence)) + 1
  aa.consensus.y <- max(as.integer(tidy.region.nt.data$Sequence)) + 2
  
  p <- ggplot()+
    geom_tile(data=tidy.region.nt.data,  aes(x = Site, y = Sequence, fill=Base))+
    scale_fill_manual(values = c("G"="#EB413C", "C"="#FFB340", "T"="#3C88EE", "A"="#64F73F", "-"="grey"))+ # Jalview colours
    labs(fill="Base")+
    
    # Label bases in each tile
    new_scale_fill()+
    geom_tile(data=tidy.region.aa.data,  aes(x = Site, y = Sequence, fill=Residue, width=0.5, height=1))+
    geom_text(data = tidy.region.aa.data, aes(x = Site, y = Sequence, label=Residue, col=isSignificantSite), 
              size=2, family="mono", )+
    scale_color_manual(values = c(`FALSE`="white", `TRUE`="black"))+
    guides(col = FALSE)+
    guides(fill = FALSE)+
    
    # Overlay diversifying selection from MEME
    new_scale_fill()+
    labs(fill="Diversifying selection")+
    geom_tile(data = tidy.region.aa.data, aes(x = Site, y = aa.consensus.y, width=3, height=2, fill = isSignificant))+
    scale_fill_manual(values = c(`TRUE`="lightblue", `FALSE`="grey"))+
    guides(fill = FALSE)+
    
    # Draw a line between each codon
    geom_segment(aes(
      x = seq(region.data$nt.start-0.5, region.data$nt.end+0.5, 3),
      xend = seq(region.data$nt.start-0.5, region.data$nt.end+0.5, 3),
      y = -Inf, yend=Inf
    ), linewidth = 0.75, col = "black")+
    
    # Set chart scales
    scale_x_continuous(breaks = seq(region.data$nt.start+1, region.data$nt.end+1, 3),
                       name = "NT site",
                       sec.axis = sec_axis(~., breaks = derive(), 
                                           labels=seq(region.data$aa.start, region.data$aa.end, 1),
                                           name = "AA site"
                                            ))+
    coord_cartesian(xlim = c(region.data$nt.start,region.data$nt.end),  ylim = c(0, aa.consensus.y+1.5), clip = 'on')+
    
     # Add NT consensus
    annotate("text", label=s2c(region.data$nt.consensus), size=2, 
             x=(region.data$nt.start):(region.data$nt.start+nchar(region.data$nt.consensus)-1),
             y= nt.consensus.y+0.4, hjust=0.5, family="mono", fontface="bold")+
    # Add AA consensus
    annotate("text", label=s2c(region.data$aa.consensus), size=2,
             x=seq((region.data$nt.start)+1,
                   (region.data$nt.start+nchar(region.data$nt.consensus)),
                   3),
             y=aa.consensus.y+0.4, hjust=0.5, family="mono", fontface="bold")+
    theme_bw()
  
  list(plot = p,
       region = region.data)
}


# Create plots of the MSA around the sites of interest
plot.individual.meme.sites <- function(meme.overview, outgroup.type){
  
  extract.function <- ifelse(outgroup.type=="combined", extract.combined.alignment.region, extract.mammal.alignment.region)
  
  site.ranges <- with(meme.overview, IRanges(start=site-2, end=site+2))
  collapsed.ranges <- as.data.frame(IRanges::reduce(site.ranges))
  unique.regions <- mapply(extract.function, aa.start = collapsed.ranges$start+1, 
                           aa.end =collapsed.ranges$end-1, SIMPLIFY = FALSE)
  
  site.plots <- lapply(unique.regions, plot.alignment.region, meme.overview=meme.overview)
  
  
  lapply(site.plots, \(x){ x$plot+theme(legend.position = "top",
                                        axis.text = element_text(size=6),
                                        axis.title.x = element_text(size=6),
                                        axis.title.y = element_blank(),
                                        legend.text = element_text(size=6),
                                        legend.title = element_text(size=8))
    save.double.width(filename = paste0("figure/meme.", outgroup.type, ".", x$region$aa.start,"-",x$region$aa.end,".png"), plot = last_plot())
  } )
}

plot.individual.meme.sites(mammal.meme.overview, "mammal")
# plot.individual.meme.sites(combined.meme.overview, "combined")

#### codeml site models to check for site-specific and branch-site selection ####
cat("Reading codeml results\n")
# Read the test and null codeml outputs
# These were created with a foreground node name
read.branch.site.codeml.output <- function(fg.node.name){
  
  calc.LRT <- function(lnl0, lnl1, np0, np1){
    lrt <- 2 * (lnl1-lnl0)
    df <- abs(np1-np0)
    crit.value =  qchisq(p=0.05, df=df, lower.tail = FALSE)
    list("crit.value" = crit.value,
         "p.value" = pchisq(lrt, df, lower.tail = FALSE),
         "lrt" = lrt)
  }
  
  test.output <- paste0("paml/branch-site-", fg.node.name, "/paml.out.txt")
  null.output <- paste0("paml/branch-site-", fg.node.name, "-null/paml.out.txt")
  
  if(!file.exists(test.output) || !file.exists(null.output)){
    return()
  }
  
  test.lines <- read_lines(test.output) %>%
    as.data.frame %>%
    dplyr::rename_with(.fn = function(i) "line") %>%
    dplyr::filter(grepl("lnL", lead(line, n=5)) | grepl("lnL", line))
  
  null.lines <- read_lines(null.output) %>%
    as.data.frame %>%
    dplyr::rename_with(.fn = function(i) "line") %>%
    dplyr::filter(grepl("lnL", lead(line, n=5)) | grepl("lnL", line))
  
  
  test.data <- data.frame(model = test.lines[seq(1, 9, 2),],
                     value = test.lines[seq(2, 10, 2),]) %>%
    dplyr::mutate(lnL = str_extract(value, "-\\d+\\.\\d+"),
                  np = str_replace(str_extract(value, "np:\\d+"), "np:", ""),
                  across(lnL:np, as.numeric)) %>%
    na.omit()
  
  null.data <- data.frame(model = null.lines[seq(1, 9, 2),],
                          value = null.lines[seq(2, 10, 2),]) %>%
    dplyr::mutate(lnL = str_extract(value, "-\\d+\\.\\d+"),
                  np = str_replace(str_extract(value, "np:\\d+"), "np:", ""),
                  across(lnL:np, as.numeric)) %>%
    na.omit()
  
  lrt.data <- calc.LRT(test.data$lnL, null.data$lnL, test.data$np, null.data$np)
  
  cat(fg.node.name, ": p=", lrt.data$p.value, "; LRT=", lrt.data$lrt, "\n")
}

# read.branch.site.codeml.output("rodentia")
read.branch.site.codeml.output("eumuroida")
# read.branch.site.codeml.output("muridae")
# read.branch.site.codeml.output("murinae")

# Test if any sites are evolving in a non-neutral manner
read.site.specific.codeml.output <- function(){
  
  # Only run this section if the codeml analysis has completed
  if(!file.exists("paml/site-specific/site.specific.paml.out.txt")){
    return()
  }
  
  # Extract lnl from output
  lines <- read_lines("paml/site-specific/site.specific.paml.out.txt") %>%
    as.data.frame %>%
    dplyr::rename_with(.fn = function(i) "line") %>%
    dplyr::filter(grepl("lnL", lead(line, n=5)) | grepl("lnL", line))
  
  data <- data.frame(model = lines[seq(1, 9, 2),],
                     value = lines[seq(2, 10, 2),]) %>%
    dplyr::mutate(lnL = str_extract(value, "-\\d+\\.\\d+"),
                  np = str_replace(str_extract(value, "np:\\d+"), "np:", ""),
                  across(lnL:np, as.numeric))
  
  
  # Process the output file to find log likelihood values to calculate LRT
  # (likelihood ratio test): twice the difference in log-likelihood ℓ between the
  # null and alternative hypotheses, 2Δℓ = 2(ℓ1 − ℓ0), where ℓ0 is the
  # log-likelihood score for the null model, whereas ℓ1 is the log-likelihood
  # under the alternative model.
  
  # ℓ is in the output file at lines starting lnL
  # Grep the lnL and previous 5 line (which has model name).
  
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
  
  pairs <- list( c(1, 2), c(2, 3), c(4, 5) )
  
  calc.values <- function(v) calc.LRT(data$lnL[v[1]], data$lnL[v[2]], data$np[v[1]], data$np[v[2]])
  
  lrts <- lapply(pairs, calc.values)
  
  # M0 vs. M1a (one-ratio vs. nearly neutral)
  
  # This is a test for variability of selective pressure among amino acid sites
  # rather than a test of positive selection. M1a fits the data much better than
  # M0,indicating that the selective pressure reflected by ω
  # varies hugely among sites.
  cat("M0 vs. M1a (is there variability of selective pressure among sites?)\n")
  cat("p=", lrts[[1]]$p.value, "lrt=", lrts[[1]]$lrt, "\n")
  
  # Compared with M1a, M2a adds a class of sites under positive selection with ω2
  # > 1 (in proportion p2). This does not improve the fit of the model
  # significantly
  # (nearly neutral vs. positive selection)
  cat("M1a, M2a(is there positive selection versus nearly neutral selection?)\n")
  cat("p=", lrts[[2]]$p.value, "lrt=", lrts[[2]]$lrt, "\n")
  
  # Additional test for positive selection by comparing M7 (beta, null model)
  # against M8 (beta&ω, alternative model).
  # (positive selection vs null model)
  # Evidence for sites under positive selection 
  cat("M7, M8 (is there positive selection at specific sites?)\n")
  cat("p=", lrts[[3]]$p.value, "lrt=", lrts[[3]]$lrt, "\n")
}

read.site.specific.codeml.output()





#### How do substitution rates compare to the divergence times? ####
# We want to calculate the number of substitutions per million years
# Combine the branch lengths with the TimeTree dates
cat("Comparing branch lengths with divergence times\n")

# Read the ML ZFX and ZFY trees 
zfx.nt.aln.tree <- ape::read.tree("aln/zfx_only/zfx.aln.treefile")
zfy.nt.aln.tree <- ape::read.tree("aln/zfy_only/zfy.aln.treefile")

# Drop the second ZFYs in mouse and rat
zfy.nt.aln.tree <- tidytree::drop.tip(zfy.nt.aln.tree, "Mouse_Zfy2") 
zfy.nt.aln.tree <- tidytree::drop.tip(zfy.nt.aln.tree, "African_Grass_Rat_ZFY2-like_1") 

# Root the trees on monotremes
zfx.nt.aln.tree <- reroot.tree(zfx.nt.aln.tree, c("Platypus_ZFX", "Australian_echidna_ZFX"), position = 0.015)
zfy.nt.aln.tree <- reroot.tree(zfy.nt.aln.tree, c("Platypus_ZFX", "Australian_echidna_ZFX"), position = 0.015)

# Remove gene names so tip labels are comparable
zfx.nt.aln.tree$tip.label <- str_replace(zfx.nt.aln.tree$tip.label, "_Z[F|f][X|x].*", "")
zfy.nt.aln.tree$tip.label <- str_replace(zfy.nt.aln.tree$tip.label, "(_putative)?(_|-)Z[F|f][X|x|Y|y].*", "")

# Replace underscores for tidy names
zfy.nt.aln.tree$tip.label <- str_replace_all(zfy.nt.aln.tree$tip.label, "_", " ")

mammal.nt.tree.data <- tidytree::as_tibble(zfy.nt.aln.tree)

# Read the species tree with divergence times to get the node divergences
species.tree <- ape::read.tree("aln/node.labeled.species.tree.nwk")
# Note node positions are reversed (0 is root). Convert to times.
node.times <- max(node.depth.edgelength(species.tree)) - node.depth.edgelength(species.tree)
# Get only the node times, ignore tips
divergence.point.data <- node.times[(length(species.tree$tip.label)+1):(length(species.tree$edge.length)+1)]
names(divergence.point.data) <- species.tree$node.label

# Create an zero vector for each tip (time for extant species is 0Mya)
divergence.point.species <- rep(0, length(zfy.nt.aln.tree$tip.label))
names(divergence.point.species) <- zfy.nt.aln.tree$tip.label
# Combine edge lengths with tips for complete weighted edges
divergence.point.data <- c(divergence.point.data, divergence.point.species)

# Use these times to create weights for each edge
zfy.nt.aln.tree$node.label[1] <- "Mammalia" # replace default 'Root' label

# Find the length of the edge from a node to its parent in Mya
get.edge.time <- function(node, tree){
  nodelabel <- treeio::nodelab(tree, node) # label of the node
  parentnode <- tree$edge[tree$edge[,2]==node,1] # parent of edge leading to node
  parentlabel <- treeio::nodelab(tree, parentnode) # label of the parent
  
  # Find the divergence time of the node and it's parent
  v1 = max(divergence.point.data[parentlabel], divergence.point.data[nodelabel])
  v2 = min(divergence.point.data[parentlabel], divergence.point.data[nodelabel])
  time = v1 - v2
  
  # If no parent node was found, zero. Otherwise, the row index in the tree containing this edge
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

# Calculate the length of each branch in Myr
time.vals <- do.call(rbind, lapply(zfy.nt.aln.tree$edge[,2], get.edge.time, tree=zfy.nt.aln.tree))


# Plot the Zfy tree, with edge lengths replaced by Myr
subs.site.mya.plot <- ggtree(zfy.nt.aln.tree, size = 1) %<+%
  time.vals +
  aes(colour = log(subsPerMyr)) +
  scale_color_paletteer_c("ggthemes::Classic Red-Blue", 
                          direction = -1, limits =c(-9, -3))+
  labs(color = "Log substitutions per site\nper million years")+

  geom_tiplab(size=2, color = "black")+
  geom_treescale(fontsize =2, y = -1) +
  coord_cartesian(xlim = c(-0.05, 0.4))+
  theme_tree() +
  theme(legend.position = c(0.2, 0.8),
        legend.background = element_blank(),
        legend.text = element_text(size=6),
        legend.title = element_text(size=6))

# Find the coordinates for the labels in the plot
zf.to.xy.y <- subs.site.mya.plot$data[subs.site.mya.plot$data$label=="Eutheria","y"]$y
eumuroida.y <- subs.site.mya.plot$data[subs.site.mya.plot$data$label=="Eumuroida","y"]$y
slxl1.y <- subs.site.mya.plot$data[subs.site.mya.plot$data$label=="Murinae","y"]$y
sly.amplifies.y <- subs.site.mya.plot$data[subs.site.mya.plot$data$label=="Mouse","y"]$y

zf.to.xy.x <- subs.site.mya.plot$data[subs.site.mya.plot$data$label=="Eutheria","x"]$x
eumuroida.x <- subs.site.mya.plot$data[subs.site.mya.plot$data$label=="Eumuroida","x"]$x
slxl1.x <- subs.site.mya.plot$data[subs.site.mya.plot$data$label=="Murinae","x"]$x
sly.amplifies.x <- subs.site.mya.plot$data[subs.site.mya.plot$data$label=="Mouse","x"]$x

# Annotate plot with labels
subs.site.mya.plot <- subs.site.mya.plot +

  # ZF* moves to sex chromosomes
  annotate("text", x=zf.to.xy.x-0.04, y=zf.to.xy.y, label="ZF* to\nX/Y", size=2, hjust=0)+
  annotate("rect", xmin=zf.to.xy.x-0.05, xmax=zf.to.xy.x-0.01, 
           ymin=zf.to.xy.y - 1, ymax=zf.to.xy.y+1, fill="darkgreen", alpha=0.4)+
  
  # Ssty box
  annotate("text", x=eumuroida.x-0.08, y=eumuroida.y, label="Ssty appears\nZfy testis specific", size=2, hjust=0)+
  annotate("rect", xmin=eumuroida.x-0.09,  xmax=eumuroida.x-0.01,
           ymin=eumuroida.y - 1, ymax=eumuroida.y + 1, fill="darkgreen", alpha=0.4)+

  # Slxl1 acquired box
  annotate("text", x=slxl1.x-0.0325, y=slxl1.y, label="Slxl1\nacquired? ", size=2, hjust=0)+
  annotate("rect", xmin=slxl1.x-0.035,  xmax=slxl1.x-0.005, 
           ymin=slxl1.y - 1,ymax=slxl1.y + 1, fill="darkgreen", alpha=0.4)+

  # Sly amplifies box
  annotate("text", x=sly.amplifies.x-0.05, y=sly.amplifies.y, label="Sly amplifies\n ", size=2, hjust=0)+
  annotate("rect", xmin=sly.amplifies.x-0.06, xmax=sly.amplifies.x-0.005, 
           ymin=sly.amplifies.y-0.75, ymax=sly.amplifies.y+0.75, fill="darkgreen", alpha=0.4)


save.double.width("figure/subs.per.site.per.Myr.png", subs.site.mya.plot)

# Also create a copy with node labels for reference
subs.site.mya.plot <- subs.site.mya.plot+ geom_nodelab(size=2, nudge_x = -0.005, nudge_y = 0.5, hjust = 1, color = "black")
save.double.width("figure/subs.per.site.per.Myr.node.labels.png", subs.site.mya.plot)


# Redraw the tree with branch lengths from the actual dates. Confirm it adds up
# to about the same length per species.

zfy.nt.aln.tree.time <- zfy.nt.aln.tree # copy the original tree

# For a given child node, find the branch time to its parent
get.time.for.node <- function(node) time.vals[time.vals$node==node,"time"]
# Get the branch times from every child node in the tree to replace edge lengths
zfy.nt.aln.tree.time$edge.length <- sapply(zfy.nt.aln.tree.time$edge[,2], get.time.for.node)

# Make fully annotated time tree
time.plot <- ggtree(zfy.nt.aln.tree.time, size = 1) %<+%
  time.vals +
  aes(colour = log(subsPerMyr)) +
  scale_color_paletteer_c("ggthemes::Classic Red-Blue", 
                          direction = -1, limits =c(-9, -3))+
  labs(color = "Log substitutions per site\nper million years")+
  # geom_nodelab(size=2, nudge_x = -3, nudge_y = 0.5, hjust = 1, color = "black")+
  geom_tiplab(size=2, color = "black")+
  theme_tree2() +
  theme(legend.position = c(0.2, 0.8),
        legend.background = element_blank(),
        legend.text = element_text(size=6),
        legend.title = element_text(size=6))
time.plot <- revts(time.plot)+scale_x_continuous(breaks=seq(0, -180, -20), labels=abs, limits = c(-185, 30))
save.double.width("figure/subs.per.site.time.png", time.plot)

# Make figures for presentations. These build up the tree and colors

eumuroida.x <- time.plot$data[time.plot$data$label=="Eumuroida","x"]$x
zf.to.xy.x <- time.plot$data[time.plot$data$label=="Eutheria","x"]$x
slxl1.x <- time.plot$data[time.plot$data$label=="Murinae","x"]$x
sly.amplifies.x <- time.plot$data[time.plot$data$label=="Mouse","x"]$x

# Scaled timetree view, no colours
time.plot.base <- ggtree(zfy.nt.aln.tree.time, size = 1, color="black") %<+%
  time.vals +
  # geom_nodelab(size=2, nudge_x = -3, nudge_y = 0.5, hjust = 1, color = "black")+
  geom_tiplab(size=2, color = "black")+
  theme_tree2() +
  theme(legend.position = "none",
        legend.background = element_blank(),
        legend.text = element_text(size=6),
        legend.title = element_text(size=6))
time.plot.base <- revts(time.plot.base)+scale_x_continuous(breaks=seq(0, -180, -20), labels=abs, limits = c(-185, 30))
save.double.width("figure/subs.per.site.time.base.png", time.plot.base)

# Scaled timetree view, with colours
time.plot.colours <- ggtree(zfy.nt.aln.tree.time, size = 1) %<+%
  time.vals +
  aes(colour = log(subsPerMyr)) +
  scale_color_paletteer_c("ggthemes::Classic Red-Blue", 
                          direction = -1, limits =c(-9, -3))+
  labs(color = "Log substitutions per site\nper million years")+
  # geom_nodelab(size=2, nudge_x = -3, nudge_y = 0.5, hjust = 1, color = "black")+
  geom_tiplab(size=2, color = "black")+
  geom_treescale(fontsize =2, y = -1, width = 10) +
  theme_tree2() +
  theme(legend.position = c(0.2, 0.8),
        legend.background = element_blank(),
        legend.text = element_text(size=6),
        legend.title = element_text(size=6))
time.plot.colours <- revts(time.plot.colours)+scale_x_continuous(breaks=seq(0, -180, -20), labels=abs, limits = c(-185, 30))
save.double.width("figure/subs.per.site.time.colours.png", time.plot.colours)

# Annotate with labels
time.plot.annotated <- ggtree(zfy.nt.aln.tree.time, size = 1) %<+%
  time.vals +
  # ZF* moves to sex chromosomes
  annotate("rect", xmin=zf.to.xy.x-50, xmax=zf.to.xy.x-20,
           ymin=zf.to.xy.y - 0.75 , ymax=zf.to.xy.y+0.75, fill="darkgreen", alpha=0.4)+
  annotate("text", x=zf.to.xy.x-45, y=zf.to.xy.y,label="ZFX moves to\nsex chromosomes", size=2, hjust=0)+
  
  
  # Ssty box
  annotate("rect", xmin=eumuroida.x-40, xmax=eumuroida.x-10, 
           ymin=eumuroida.y - 0.75, ymax=eumuroida.y +0.75, fill="darkgreen", alpha=0.4)+
  annotate("text", x=eumuroida.x-35, y=eumuroida.y, label="Ssty appears\nZfy testis specific", size=2, hjust=0)+
  
  
  # Slxl1 acquired box
  annotate("rect", xmin=slxl1.x-27, xmax=slxl1.x-4,
           ymin=slxl1.y + 0.25, ymax=slxl1.y+1.75, fill="darkgreen", alpha=0.4)+
  annotate("text", x=slxl1.x-5, y=slxl1.y+1, label="Slxl1 acquired? ", size=2, hjust=1)+
  annotate("segment", slxl1.x-5, xend = slxl1.x, y = slxl1.y+1, yend = slxl1.y, col ="darkgreen")+
  
  
  # Sly amplifies box
  annotate("rect", 
           xmin=sly.amplifies.x-10, xmax=sly.amplifies.x, 
           ymin=sly.amplifies.y-0.75, ymax=sly.amplifies.y+0.75, fill="darkgreen", alpha=0.4)+
  annotate("text", x=sly.amplifies.x-10, y=sly.amplifies.y-0.5, label="Sly\namplifies\n ", size=2, hjust=0)+

  aes(colour = log(subsPerMyr)) +
  scale_color_paletteer_c("ggthemes::Classic Red-Blue", 
                          direction = -1, limits =c(-9, -3))+
  labs(color = "Log substitutions per site\nper million years")+
  # geom_nodelab(size=2, nudge_x = -3, nudge_y = 0.5, hjust = 1, color = "black")+
  geom_tiplab(size=2, color = "black")+
  theme_tree2() +
  theme(legend.position = c(0.2, 0.8),
        legend.background = element_blank(),
        legend.text = element_text(size=6),
        legend.title = element_text(size=6))

time.plot.annotated <- revts(time.plot.annotated)+scale_x_continuous(breaks=seq(0, -180, -20), labels=abs, limits = c(-185, 30))
save.double.width("figure/subs.per.site.time.annotated.png", time.plot.annotated)

#### Read final intron ML trees ####
cat("Reading final intron trees\n")

final.intron.zfy.nt.aln.tree <- ape::read.tree(FILES$final.intron.zfy.nt.filt.aln.treefile) %>%
  reroot.tree(., c("Opossum_ZFX", "Koala_ZFX"), position = 0.015)

final.intron.zfx.nt.aln.tree <- ape::read.tree(FILES$final.intron.zfx.nt.filt.aln.treefile) %>%
  reroot.tree(., c("Opossum_ZFX", "Koala_ZFX"), position = 0.015)

final.intron.zfy.nt.divvy.aln.tree <- ape::read.tree(FILES$final.intron.zfy.nt.aln.divvy.aln.treefile) %>%
  reroot.tree(., c("Opossum_ZFX", "Koala_ZFX"), position = 0.015)

final.intron.zfx.nt.divvy.aln.tree <- ape::read.tree(FILES$final.intron.zfx.nt.aln.divvy.aln.treefile) %>%
  reroot.tree(., c("Opossum_ZFX", "Koala_ZFX"), position = 0.015)

mammal.gene.groups <- split(METADATA$combined$common.name, METADATA$combined$group)

final.intron.nt.aln.tree <- ape::read.tree(FILES$final.intron.nt.filt.aln.treefile) %>%
  reroot.tree(., c("Opossum_ZFX", "Koala_ZFX"), position = 0.015) %>%
  tidytree::groupOTU(., mammal.gene.groups, group_name = "group")

final.intron.nt.divvy.aln.tree <- ape::read.tree(FILES$final.intron.nt.aln.divvy.aln.treefile) %>%
  reroot.tree(., c("Opossum_ZFX", "Koala_ZFX"), position = 0.015) %>%
  tidytree::groupOTU(., mammal.gene.groups, group_name = "group")

#### Plot the final intron trees ####

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
final.intron.nt.divvy.aln.tree.plot <- plot.tree(final.intron.nt.divvy.aln.tree, col= "group")  + xlim(0, 2) + labs(title = "ZFX/Y (divvied)")
save.double.width("figure/final.intron.zfx.zfy.divvy.tree.png", final.intron.nt.divvy.aln.tree.plot)

# Panel figure
# combined.plot <- (final.intron.zfy.nt.aln.tree.plot + final.intron.zfy.nt.divvy.aln.tree.plot) / (final.intron.zfx.nt.aln.tree.plot + final.intron.zfx.nt.divvy.aln.tree.plot)
# save.double.width("figure/final.intron.combined.tree.png", combined.plot)

#### Plot the final intron MSAs ####

plot.msa <- function(tidy.aln){
  ggplot()+
    geom_msa(data = tidy.aln, seq_name = T, font=NULL, 
             border=NA, color="Chemistry_NT", consensus_views = F )+
    coord_cartesian(expand = FALSE)+
    scale_x_continuous(breaks = seq(0, max(tidy.aln$position), 100))+
    scale_y_discrete( labels = \(x) str_replace_all(x, "_", " ") )+
    theme_bw()+
    theme(axis.title = element_blank())
}

# Read alignment, set order to match phylogeny

zfy.raw.aln <- tidy_msa(Biostrings::readDNAMultipleAlignment("aln/final.intron.zfy/final.intron.zfy.nt.aln", format="fasta"))
zfy.raw.aln$name <- factor(zfy.raw.aln$name, levels = mammal.taxa.name.order)

zfy.divvy.aln <- tidy_msa(Biostrings::readDNAMultipleAlignment("aln/final.intron.zfy/final.intron.zfy.nt.aln.divvy.aln", format="fasta"))
zfy.divvy.aln$name <- factor(zfy.divvy.aln$name, levels = mammal.taxa.name.order)

zfx.raw.aln <- tidy_msa(Biostrings::readDNAMultipleAlignment("aln/final.intron.zfx/final.intron.zfx.nt.aln", format="fasta"))
zfx.raw.aln$name <- factor(zfx.raw.aln$name, levels = mammal.taxa.name.order)

zfx.divvy.aln <- tidy_msa(Biostrings::readDNAMultipleAlignment("aln/final.intron.zfx/final.intron.zfx.nt.aln.divvy.aln", format="fasta"))
zfx.divvy.aln$name <- factor(zfx.divvy.aln$name, levels = mammal.taxa.name.order)

zf.raw.aln <- tidy_msa(Biostrings::readDNAMultipleAlignment("aln/final.intron/final.intron.nt.aln", format="fasta"))
zf.raw.aln$name <- factor(zf.raw.aln$name, levels = mammal.taxa.name.order)

zf.divvy.aln <- tidy_msa(Biostrings::readDNAMultipleAlignment("aln/final.intron/final.intron.nt.aln.divvy.aln", format="fasta"))
zf.divvy.aln$name <- factor(zf.divvy.aln$name, levels = mammal.taxa.name.order)

zfy.msa.raw.plot <- plot.msa(zfy.raw.aln)
zfy.msa.plot <- plot.msa(zfy.divvy.aln)

zfx.msa.raw.plot <- plot.msa(zfx.raw.aln)
zfx.msa.plot <- plot.msa(zfx.divvy.aln)

zf.msa.raw.plot <- plot.msa(zf.raw.aln)
zf.msa.plot <- plot.msa(zf.divvy.aln)

save.plot("figure/final.intron.msa.zfy.raw.png", zfy.msa.raw.plot, width=270, height = 170)
save.plot("figure/final.intron.msa.zfy.divvy.png", zfy.msa.plot, width=270, height = 170)
save.plot("figure/final.intron.msa.zfx.raw.png", zfx.msa.raw.plot, width=270, height = 170)
save.plot("figure/final.intron.msa.zfx.divvy.png", zfx.msa.plot, width=270, height = 170)
save.plot("figure/final.intron.msa.zf.raw.png", zf.msa.raw.plot, width=270, height = 170)
save.plot("figure/final.intron.msa.zf.divvy.png", zf.msa.plot, width=270, height = 170)


#### Plot full protein MSA ####
full.aa.msa <- tidy_msa(Biostrings::readAAMultipleAlignment("aln/combined/combined.aa.aln", format="fasta")) %>%
  dplyr::mutate(name = factor(name, levels = rev(combined.taxa.name.order) ))
full.aa.msa.plot <- ggplot()+
  geom_msa(data = full.aa.msa, seq_name = T, font=NULL, 
           border=NA, color="Chemistry_AA", consensus_views = F )+
  coord_cartesian(expand = FALSE)+
  scale_x_continuous(breaks = seq(0, max(full.aa.msa$position), 100))+
  scale_y_discrete( labels = \(x) str_replace_all(x, "_", " ") )+
  theme_bw()+
  theme(axis.title = element_blank(),
        axis.text.y = element_text(size = 6))

save.plot("figure/full.msa.png", full.aa.msa.plot, width=270, height = 170)

#### Plot Therian and Eutherian ancestral ZFX/ZFY MSAs ####

# How similar are the reconstructed sequences?

ancestral.msa <- ancestral.seqs %>% 
  dplyr::filter(Node %in% c("Theria", "Eutheria")) %>%
  dplyr::mutate(name = paste0(Node,"_", Type),
                name = factor(name, levels = c("Theria_ZFY", "Eutheria_ZFY", "Theria_ZFX", "Eutheria_ZFX"))) %>%
  dplyr::select(name, position = Site, character = State)

ancestral.msa.plot <- ggplot()+
  geom_msa(data = ancestral.msa, seq_name = T, font=NULL, 
           border=NA, color="Chemistry_AA", consensus_views = T )+
  coord_cartesian(expand = FALSE)+
  scale_x_continuous(breaks = seq(0, max(ancestral.msa$position), 100))+
  scale_y_discrete( labels = \(x) str_replace_all(x, "_", " ") )+
  theme_bw()+
  theme(axis.title = element_blank())

save.plot("figure/ancestral.msa.png", ancestral.msa.plot, width=270, height = 50)


#### Tar the figure outputs ####
system2("tar", "czf figure.tar.gz figure")
cat("Done!\n")

