# Create figures

#### Imports #####

source("src/functions.R")
load.packages()

source("src/find9aaTADs.R")
source("src/findZF.R")
source("src/calcCharge.R")
source("src/calcHydrophobicity.R")

cat("Packages loaded\n")
filesstrings::create_dir("figure")
ALIGNMENTS <- read.alignments()

# Identify the coordinates of the exon boundaries in the gapped alignments
# Based on the mouse Zfy1 sequence
mouse.exons <- find.exons()

#### Plot combined mammal/outgroup AA tree ####
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

combined.outgroup.tree <- read.combined.outgroup.tree(FILES$combined.aa.aln.treefile)

combined.aa.tree <- plot.tree(combined.outgroup.tree, col = "group")+
  coord_cartesian(clip="off", xlim = c(0, 0.7), ylim= c(-2, length(combined.outgroup.tree$tip.label)))+
  annotate("text", x = 0.32, y = 58.5, label = "Eumuroida", size=2)

save.double.width("figure/Figure_1_aa_tree.png", combined.aa.tree)

# Also save the order of taxa in the outgroup tree to use later
combined.taxa.name.order <- ggtree::get_taxa_name(combined.aa.tree) 

#### Plot mammal CDS NT tree #####

nt.aln.tree <- ape::read.tree(FILES$combined.nt.aln.treefile)

# Root the tree on platypus and resave
# The root is arbitrarily placed in the platypus branch to fit neatly
xenopus.node <- ape::getMRCA(nt.aln.tree, c("Xenopus_ZFX.S","Xenopus_ZFX.L"))
nt.aln.tree <- phytools::reroot(nt.aln.tree, xenopus.node, position = 0.01)
# nt.aln.tree <- phytools::reroot(nt.aln.tree, which(nt.aln.tree$tip.label=="Platypus_ZFX"), position = 0.015)
ape::write.tree(nt.aln.tree, file = paste0(FILES$combined.nt.aln, ".rooted.treefile"))

# Find the nodes that are ZFY vs ZFX and add to tree
mammal.gene.groups <- split(METADATA$combined$common.name, METADATA$combined$group)
nt.aln.tree <- tidytree::groupOTU(nt.aln.tree, mammal.gene.groups, group_name = "group")

plot.zfx.zfy <- plot.tree(nt.aln.tree, col= "group") + coord_cartesian(clip="off", xlim = c(0, 0.5))

# Emphasise the ZFX / ZFY splits more by rotating the Laurasiatheria node
laurasiatheria.node <- ape::getMRCA(nt.aln.tree, c("Cat_ZFY", "Cat_ZFX"))
plot.zfx.zfy <- rotate(plot.zfx.zfy, laurasiatheria.node)

save.double.width("figure/Figure_xxxx_nucleotide_tree.png", plot.zfx.zfy)

# Get the order of taxa names for reordering other plots later
mammal.taxa.name.order <- get_taxa_name(plot.zfx.zfy) 

# Now we want to view separate trees for Zfx and Zfy based on this combined tree

# Drop the ZFY sequences and just look at the ZFX nodes in the tree
tree.zfx <- ape::drop.tip(nt.aln.tree, mammal.gene.groups$ZFY)
tree.zfx <- groupOTU(tree.zfx, mammal.gene.groups, group_name = "group")
ape::write.tree(tree.zfx, file = paste0(FILES$mammal.nt.aln, ".zfx.treefile"))
plot.zfx <- plot.tree(tree.zfx)
save.double.width("figure/mammal.zfx.tree.png", plot.zfx)

# Keep Zfy and outgroups, drop other tips
tree.zfy <- ape::keep.tip(nt.aln.tree, c(mammal.gene.groups$ZFY, mammal.gene.groups$Outgroup))
tree.zfy <- groupOTU(tree.zfy, mammal.gene.groups, group_name = "group")
ape::write.tree(tree.zfy, file = paste0(FILES$mammal.nt.aln, ".zfy.treefile"))
plot.zfy <- plot.tree(tree.zfy)
save.double.width("figure/mammal.zfy.tree.png", plot.zfy)

#### Plot mammal exon NT trees ####

create.exon.plot <- function(i){
  exon.aln.file <- paste0("aln/exons/exon_", mouse.exons$exon[i], ".aln")

  # Some exons will fail - too many gaps
  if(!file.exists(paste0(exon.aln.file, ".treefile"))) return()
  
  exon.tree <- ape::read.tree(paste0(exon.aln.file, ".treefile"))
  # Root the tree on platypus
  exon.tree <- phytools::reroot(exon.tree, which(exon.tree$tip.label=="Platypus_ZFX"), position = 0.015)
  
  # Find the nodes that are ZFY vs ZFX and add to tree
  mammal.gene.groups <- split(METADATA$mammal$common.name, METADATA$mammal$group)
  exon.tree <- groupOTU(exon.tree, mammal.gene.groups, group_name = "group")
  
  plot.exon.tree <- plot.tree(exon.tree, tiplab.font.size = 1.5,  col="group")  + 
    coord_cartesian(clip="off", xlim = c(0.16, 0.83))
  # exon.fig.file <- paste0("figure/exon_", mouse.exons$exon[i], ".zfx.zfy.tree.png")
  # save.double.width(exon.fig.file, plot.exon.tree)
  
  # Return for playing
  plot.exon.tree
}

exon.1.7.plots <- lapply(1:nrow(mouse.exons),create.exon.plot)

# We also want to look at all except exon 7
exon.1.6.aln.file <- paste0("aln/exons/exon_1-6.aln")
exon.1.6.tree <- ape::read.tree(paste0(exon.1.6.aln.file, ".treefile"))
# Root the tree on platypus
exon.1.6.tree <- phytools::reroot(exon.1.6.tree, which(exon.1.6.tree$tip.label=="Platypus_ZFX"), position = 0.015)
# Find the nodes that are ZFY vs ZFX and add to tree
# Find the nodes that are ZFY vs ZFX and add to tree
mammal.gene.groups <- split(METADATA$mammal$common.name, METADATA$mammal$group)
exon.1.6.tree <- groupOTU(exon.1.6.tree, mammal.gene.groups, group_name = "group")

plot.exon.1.6.tree <- plot.tree(exon.1.6.tree, tiplab.font.size = 1.5,  col="group")  + coord_cartesian(clip="off", xlim = c(0.16, 0.83))
# exon.1.6.fig.file <- paste0("figure/exon_1-6.zfx.zfy.tree.png")
# save.double.width(exon.1.6.fig.file, plot.exon.1.6.tree)

# Create a joint figure of exons 1-6, exon 2, and exon 7

exon.joint.tree <- plot.exon.1.6.tree + exon.1.7.plots[[2]] + exon.1.7.plots[[7]] + 
  patchwork::plot_annotation(tag_levels = list(c("Exons 1-6", "Exon 2", "Exon 7"))) &
  theme(plot.tag = element_text(size = 6),
        plot.margin = margin(t=0, l=0, r=0, b=0))
save.double.width("figure/Figure_Sxxxx_exons_tree.png", exon.joint.tree, height=120)

#### Test selection globally in mammals ####

# ape::dnds(ALIGNMENTS$nt.mammal.ape) # errors
seqin.aln <- seqinr::read.alignment(FILES$mammal.nt.aln, format = "fasta")
kaks.data <- seqinr::kaks(seqin.aln)

kaks.ratio <- kaks.data$ka / kaks.data$ks
kaks.pairwise <- metagMisc::dist2list(kaks.ratio, tri = F) %>%
  dplyr::mutate(col = str_replace_all(col, "_", " "),
                row = str_replace_all(row, "_", " "),
                col = factor(col, levels = mammal.taxa.name.order),
                row = factor(row, levels = mammal.taxa.name.order),
                colnum = as.integer(col),
                rownum = as.integer(row)) %>%
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
exon1.3_6.aln <- as.matrix(ALIGNMENTS$nt.mammal.ape)[,exon1.3_6.locs]

ape::write.FASTA(exon1.3_6.aln, file = "aln/exons/exon_1_3-6.kaks.aln")
seqin.aln.exon.1.3_6 <- seqinr::read.alignment("aln/exons/exon_1_3-6.kaks.aln", format = "fasta")

exon2.aln <- as.matrix(ALIGNMENTS$nt.mammal.ape)[,mouse.exons$start_nt_codon_offset[2]:(mouse.exons$end_nt_codon_offset[2])]
ape::write.FASTA(exon2.aln, file = "aln/exons/exon_2.kaks.aln")
seqin.aln.exon.2 <- seqinr::read.alignment("aln/exons/exon_2.kaks.aln", format = "fasta")
exon2.aln.msa <- ape::read.FASTA("aln/exons/exon_2.kaks.aln")

exon.7.aln <- as.matrix(ALIGNMENTS$nt.mammal.ape)[,(mouse.exons$start_nt_codon_offset[7]):mouse.exons$end_nt_codon_offset[7]] 
ape::write.FASTA(exon.7.aln, file = "aln/exons/exon_7.kaks.aln")
seqin.aln.exon.7 <- seqinr::read.alignment("aln/exons/exon_7.kaks.aln", format = "fasta")

create.pairwise.kaks.data <- function(seqinr.aln){
  kaks.data <- seqinr::kaks(seqinr.aln, rmgap = FALSE)
  kaks.ratio <- kaks.data$ka / kaks.data$ks
  
  kaks.pairwise <- metagMisc::dist2list(kaks.ratio, tri = F) %>%
    dplyr::mutate(col = str_replace_all(col, "_", " "),
                  row = str_replace_all(row, "_", " "),
                  col = factor(col, levels = mammal.taxa.name.order),
                  row = factor(row, levels = mammal.taxa.name.order),
                  colnum = as.integer(col),
                  rownum = as.integer(row),
                  value  = ifelse(value==1, NA, value)) %>%  # values of exactly 1 are from missing data)
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

#### Identify the locations of the ZFs in the AA & NT MSAs ####

locations.zf <- locate.zfs.in.alignment(combined.taxa.name.order)

write_tsv(locations.zf %>% 
            dplyr::select(sequence, aa_motif, start_ungapped, end_ungapped, 
                          start_gapped, end_gapped, start_nt_ungapped, end_nt_ungapped, 
                          start_nt_gapped, end_nt_gapped),
          "figure/locations.zf.tsv")

ranges.ZF.common <- merge(find.common.aa.overlaps(locations.zf), find.common.nt.overlaps(locations.zf), by = "motif_number")

# Annotate the individual ZFs with which consensus motif they belong to (0 if none)
locations.zf %<>% 
  dplyr::rowwise() %>%
  dplyr::mutate(motif_number = find.matching.range("start", "end", start_gapped, end_gapped, ranges.ZF.common )) %>%
  dplyr::ungroup()

#### Identify the locations of the 9aaTADs in the AA & NT MSAs ####

locations.9aaTAD <- locate.9aaTADs.in.alignment(FILES$combined.aa.aln, FILES$combined.nt.aln, combined.taxa.name.order)

write_tsv(locations.9aaTAD %>% 
            dplyr::select(sequence, aa_motif, rc_score, start_ungapped, end_ungapped, 
                          start_gapped, end_gapped, start_nt_ungapped, end_nt_ungapped, 
                          start_nt_gapped, end_nt_gapped),
          "figure/locations.9aaTAD.tsv")

# Only keep the high confidence 9aaTADs for the track
ranges.9aaTAD.common <-  merge(find.common.aa.9aaTADs(locations.9aaTAD, 
                                                      rc.threshold = 80,
                                                      coverage.threshold = 21),
                               find.common.nt.9aaTADs(locations.9aaTAD, 
                                                      rc.threshold = 80,
                                                      coverage.threshold = 21), 
                               by = c("motif_number", "label"), all.x = T)


locations.9aaTAD %<>%
  dplyr::rowwise() %>%
  dplyr::mutate(motif_number = find.matching.range("start", "end", 
                                                   start_gapped, end_gapped, ranges.9aaTAD.common )) %>%
  dplyr::ungroup() %>%
  merge(., ranges.9aaTAD.common, by = "motif_number", all.x = TRUE)

#### Identify the locations of the NLS in the AA & NT MSAs ####

locations.NLS <- locate.NLS.in.alignment(FILES$combined.aa.aln, FILES$combined.nt.aln, combined.taxa.name.order)

# Export the locations of the NLS
write_tsv(locations.NLS %>% 
            dplyr::select(sequence, aa_motif, type, posterior_prob, start_ungapped, end_ungapped, 
                          start_gapped, end_gapped, start_nt_ungapped, end_nt_ungapped, 
                          start_nt_gapped, end_nt_gapped),
          "figure/locations.NLS.tsv")

# We need to combine the full set of structure locations into an overlapping set
# to be plotted in a single row. Keep those that overlap in >=5 species
ranges.NLS.common <- merge(find.common.aa.overlaps(locations.NLS), 
                           find.common.nt.overlaps(locations.NLS), 
                           by = "motif_number")


# Annotate the individual NLS with which consensus motif they belong to (0 if none)
locations.NLS %<>% 
  dplyr::rowwise() %>%
  dplyr::mutate(motif_number = find.matching.range("start", "end", start_gapped, end_gapped, ranges.NLS.common )) %>%
  dplyr::ungroup()

#### Export the residues within each ZF ####

# Full ZF motifs
locations.zf %>%
  dplyr::select(Sequence = sequence, motif_number, aa_motif) %>%
  dplyr::arrange(desc(motif_number)) %>% # so we display 13 - 1
  dplyr::mutate(motif_number = paste0("ZF_", motif_number)) %>%
  tidyr::pivot_wider(id_cols = Sequence, names_from = motif_number, values_from = aa_motif) %>%
  as.data.frame %>%
  dplyr::arrange(as.integer(Sequence)) %>%
  create.xlsx(., "figure/locations.zf.xlsx", cols.to.fixed.size.font = 2:14)

#### Export the contact bases  within each ZF ####

zf.contact.bases <- locations.zf %>%
  dplyr::select(Sequence = sequence, motif_number, contact_bases) %>%
  dplyr::arrange(desc(motif_number)) %>% # so we display 13 - 1
  dplyr::mutate(motif_number = paste0("ZF_", motif_number)) %>%
  tidyr::pivot_wider(id_cols = Sequence, names_from = motif_number, values_from = contact_bases) %>%
  as.data.frame %>%
  dplyr::arrange(as.integer(Sequence))

# What are the most common ZF contact motifs?
zf.contact.bases.conserved <- locations.zf %>%
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
    for(col in 2:14){
      cell.name <- cell.names[col]
      cell.ref <- cells[[cell.name]]
      zf.motif <- paste0("ZF_",(15-col)) # zfs are in descending order
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

locations.NLS %>%
  dplyr::select(Sequence = sequence, motif_number, aa_motif) %>%
  dplyr::mutate(motif_number = paste0("NLS_", motif_number)) %>%
  tidyr::pivot_wider(id_cols = Sequence, names_from = motif_number, values_from = aa_motif,
                     values_fn = ~paste(.x, collapse = ", ")) %>%
  as.data.frame %>%
  dplyr::arrange(as.integer(Sequence)) %>%
  create.xlsx(., "figure/locations.NLS.xlsx", cols.to.fixed.size.font = 2:4)

#### Export the residues within each 9aaTAD ####

# All unique 9aaTADs identified
locations.9aaTAD %>%
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
extract.superTAD.motifs <- function(locations.9aaTAD){
  
  # Look at all superTAD locations in all sequences
  sequences <- unique(locations.9aaTAD$sequence)
  tad.labels <- na.omit(unique(locations.9aaTAD$label))
  combos <- expand.grid("sequence" = sequences, "tad.label"= tad.labels)
  combos <- merge(combos, ranges.9aaTAD.common, by.x = "tad.label", by.y = "label",
                  all.x = TRUE)
  
  # Combine the superTAD locations for each row
  locations.superTAD <- merge(combos, locations.9aaTAD, 
                              by.x = c("sequence", "tad.label", "start",
                                       "end","width", "motif_number"),
                              by.y = c("sequence", "label", "start",
                                       "end", "width", "motif_number"), all.x = TRUE) %>%
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
    dplyr::select(Sequence = sequence, tad.label, rc_score, tad.sequence) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Sequence, tad.label) %>%
    dplyr::arrange(Sequence, tad.label, desc(rc_score)) %>%
    dplyr::slice_head(n=1) %>% # if there are multiple rc score per supertad, take only the highest
    dplyr::distinct() %>%
    
    # Pivot to each 9aaTAD as a separate column
    tidyr::pivot_wider(id_cols = c(Sequence), names_from = tad.label, 
                       values_from = c(tad.sequence, rc_score), values_fn = ~paste(.x, collapse = ", "),
                       names_glue = "9aaTAD_{tad.label}_{.value}") %>%
    as.data.frame %>%
    dplyr::mutate(across(`9aaTAD_A_rc_score`:`9aaTAD_G_rc_score`, as.numeric)) %>%
    dplyr::rename_with(.fn = function(x) gsub("_tad.sequence", "", x)) %>%
    dplyr::arrange(as.integer(Sequence))
  
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
      }
    }
  }
  
  wb = xlsx::createWorkbook(type = "xlsx")
  sh = xlsx::createSheet(wb)
  xlsx::addDataFrame(locations.superTAD[,1:8], sh, row.names = F)
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
  for(i in 2:length(xlsx::getRows(sh))){
    cells <- getCells(rows[i]) 
    cell.names <- names(cells)
    for(col in 2:8){
      cell.name <- cell.names[col]
      cell.ref <- cells[[cell.name]]
      rc.val <-  locations.superTAD[i-1, col+7]
      
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
  
  for(col in 2:8) set.rich.text.on.vv(wb, sh, col)
  
  
  xlsx::autoSizeColumn(sh, 1:ncol(locations.superTAD))
  xlsx::saveWorkbook(wb, file="figure/locations.9aaTAD.xlsx")
}

extract.superTAD.motifs(locations.9aaTAD)

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
  dplyr::group_by(zf_target_1, zf_target_2) %>%
  dplyr::summarise(total = n(), sequences = paste(sequence, collapse = ", ")) %>%
  dplyr::mutate(sequences = ifelse(total>40, "All others", sequences))
write_tsv(pwm.predictions, "figure/zf_targets.tsv")

#### Ancestral sequence reconstruction #####

# Read the ancestral reconstruction
ancestral.nt.seqs <- read.table(paste0(FILES$mammal.nt.aln, ".state") ,header=TRUE)

# We care about the eutherian common ancestor and the rodent ancestor
# Find these nodes

# Read the tree back to keep full node names
nt.aln.tree.nodes <- ape::read.tree(paste0(FILES$mammal.nt.aln, ".treefile"))
nt.aln.tree.nodes <- ape::root(nt.aln.tree.nodes, "Platypus_ZFX")

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
            fill=NLS.COLOUR, alpha=0.5)
aa.structure.plot <- annotate.structure.plot(aa.structure.plot, length(combined.taxa.name.order) + 1.5)

save.double.width("figure/Figure_Sxxxx_aa.structure.png", aa.structure.plot, height = 120)

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

save.double.width("figure/Figure_2_aa.structure.confident.png", aa.structure.confident.plot, height = 120)


#### Plot the conservation of hydrophobicity across mammal/outgroup AA MSA ####

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
  scale_fill_paletteer_c("ggthemes::Classic Red-Blue", direction = -1, limits = c(0, 1))+
  labs(fill="Hydrophobicity (9 site average)")
hydrophobicity.plot <- annotate.structure.plot(hydrophobicity.plot, n.taxa)
save.double.width("figure/hydrophobicity.convervation.tree.png", hydrophobicity.plot, height = 120)

#### Plot conservation of charge across mammal/outgroup AA MSA ####

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
  scale_fill_paletteer_c("ggthemes::Classic Red-Black", direction = -1, limits = c(-1, 1))+
  labs(fill="Charge (9 site average)")
charge.plot <- annotate.structure.plot(charge.plot, n.taxa)
save.double.width("figure/charge.convervation.tree.png", charge.plot, height = 120)

#### Combine all structure plots ####

# To test spacing and balance
structure.plot <- (hydrophobicity.plot) / (charge.plot) + plot_layout(ncol = 1)
save.double.width("figure/Figure_xxxx_combined_structure.plot.png", structure.plot, height = 230)

#### Plot HyPhy RELAX test for relaxed selection in Muroidea ####

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
  
  hyphy.input.tree <- phytools::reroot(hyphy.input.tree, which(hyphy.input.tree$tip.label=="Platypus_ZFX"), position = 0.015)
  
  p <- ggtree(hyphy.input.tree) + 
    geom_tiplab(size=2)
  p <- p %<+% k.vals + aes(colour=adj.k) + 
    scale_color_paletteer_c("ggthemes::Classic Red-Blue",
                            direction = 1,
                            limits = c(-3, 3))+
                            # 
                            # labels = c("-2", "0", "2", round(25, digits = 0), round(50, digits = 0)))+
    labs(color = "log(K)")+
    coord_cartesian(xlim = c(0, 0.7), clip = "off")+
    geom_treescale(fontsize =2, y = -1, width = 0.05) +
    # annotate(geom="text", x = 0.3, y = 62, size = 2,
    #          label = paste("K(Muroidea) =", round(hyphy.data$`test results`$`relaxation or intensification parameter`, digits = 2)))+
    theme(legend.position = c(0.8, 0.5),
          legend.background = element_blank(),
          legend.text = element_text(size=6),
          legend.title = element_text(size=6),
          plot.margin = margin(r=0, l=0))
  
  p
}

if(file.exists("hyphy/mammal.relax.json")){
  relax.tree.all <- create.relax.k.tree("hyphy/mammal.relax.json")
  save.double.width("figure/Figure_xxxx_RELAX_mammal.png", relax.tree.all)
  
  relax.tree.e1_3_6 <- create.relax.k.tree("hyphy/exon_1_3-6.relax.json")
  relax.tree.e2 <- create.relax.k.tree("hyphy/exon_2.relax.json")
  relax.tree.e7 <- create.relax.k.tree("hyphy/exon_7.relax.json")

  # Combine the exon trees into one figure
  exon.relax.tree <- relax.tree.e1_3_6 + relax.tree.e2 + relax.tree.e7 + 
    plot_annotation(tag_levels = list(c("Exons 1,3-6", "Exon 2", "Exon 7"))) +
    plot_layout(guides = "collect") &
    theme(legend.position='bottom', 
          plot.tag = element_text(size=6),
          plot.tag.position = c(0.2, 1))
  save.double.width("figure/exons.relax.K.png", exon.relax.tree)
}




#### codeml site model to check for site-specific selection ####

# Only run this section if the codeml analysis has completed
if(file.exists("paml/site-specific/site.specific.paml.out.txt")){
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
  cat("M0 vs. M1a\n")
  cat(paste(lrts[[1]], collapse = " | "), "\n")
  
  # Compared with M1a, M2a adds a class of sites under positive selection with ω2
  # > 1 (in proportion p2). This does not improve the fit of the model
  # significantly
  # (nearly neutral vs. positive selection)
  cat("M1a, M2a\n")
  cat(paste(lrts[[2]], collapse = " | "), "\n")
  
  # Additional test for positive selection by comparing M7 (beta, null model)
  # against M8 (beta&ω, alternative model).
  # (positive selection vs null model)
  # Evidence for sites under positive selection 
  cat("M7, M8\n")
  cat(paste(lrts[[3]], collapse = " | "), "\n")
  
  
  #### codeml branch-site model to look for selection specifically in Muroidea ####
  
  # Bayes Empirical Bayes (BEB) analysis (Yang, Wong & Nielsen 2005. Mol. Biol. Evol. 22:1107-1118)
  # Extract sites under positive selection
  positive.sites <- read_lines("paml/branch-site/branch-site.paml.out.txt") %>%
    as.data.frame %>%
    dplyr::rename_with(.fn = function(i) "line") %>%
    dplyr::filter(grepl("^ {2,5} \\d{1,3} [A-Z\\-]", line)  ) %>%
    dplyr::mutate(line = str_trim(line)) %>%
    tidyr::separate_wider_delim(line, delim = " ", names = c("site", "aa", "p")) %>%
    dplyr::mutate(p = str_replace_all(p, "\\*", ""),
                  site = as.numeric(site),
                  p = as.numeric(p))
  
  positive.sites.y <- 1
  
  positive.sites.plot <- ggplot()+
    geom_rect(data = positive.sites,   aes(xmin=site-0.5, xmax=site+0.5, ymin=positive.sites.y, ymax=positive.sites.y+2, fill=p))+
    geom_rect(data = positive.sites[positive.sites$p>0.9,],   aes(xmin=site-0.5, xmax=site+0.5, ymin=positive.sites.y+3, ymax=positive.sites.y+5, fill=p))+
    labs(x = "Site", fill = "p(ω>1)")+
    scale_fill_viridis_c()+
    theme_bw()
  
  # msa.aa.aln.tidy.frog.conservation    <- calculate.conservation(ALIGNMENTS$aa.combined.biostrings,"Xenopus_ZFX.S" )
  # msa.aa.aln.tidy.chicken.conservation <- calculate.conservation(ALIGNMENTS$aa.combined.biostrings,"Chicken_ZFX" )
  # msa.aa.aln.tidy.opossum.conservation <- calculate.conservation(ALIGNMENTS$aa.combined.biostrings,"Opossum_ZFX" )
  
  n.taxa <- 2
  positive.sites.plot <-  positive.sites.plot+
    # Draw the conservation with Xenopus, chicken and opossum
    new_scale_fill()+
    scale_fill_viridis_c(limits = c(0, 1))+
    labs(fill="Conservation (5 site average)")+
    # add.conservation.track(msa.aa.aln.tidy.frog.conservation,    n.taxa,   n.taxa+2)+
    # add.conservation.track(msa.aa.aln.tidy.chicken.conservation, n.taxa+3, n.taxa+5)+
    # add.conservation.track(msa.aa.aln.tidy.opossum.conservation, n.taxa+6, n.taxa+8)+
    add.track(ranges.NLS.common,    n.taxa+4.5, n.taxa+7.5, fill=NLS.COLOUR, alpha = 1)+
    # Draw the structures
    add.track(ranges.ZF.common,     n.taxa+5, n.taxa+7, fill=ZF.COLOUR)+
    add.track.labels(ranges.ZF.common, n.taxa+5, n.taxa+7, col="white", label_col = "motif_number")+
    
    add.track(ranges.9aaTAD.common, n.taxa+5, n.taxa+7, fill=TAD.COLOUR,  alpha = 0.9)+ # fill color max from "grDevices::Blues 3"
    add.track.labels(ranges.9aaTAD.common, n.taxa+3, n.taxa+7, col="white")+   # Label the 9aaTADs
    
    new_scale_fill()+
    scale_fill_manual(values=c("white", "grey", "white", "grey", "white", "grey", "white"))+
    scale_pattern_color_manual(values=c("white", "white"))+
    
    scale_pattern_manual(values = c("none", "stripe")) + # which exons are patterned
    guides(fill = "none", pattern="none")+
    add.exon.track(n.taxa+8, n.taxa+10, col = "black")+ # color of border
    add.exon.labels(n.taxa+8, n.taxa+10)+
    
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
          panel.grid = element_blank())
  
  save.double.width("figure/positive.sites.png", positive.sites.plot, height = 85)
  
  
  # Plot the fraction of sites under possible selection in a 9-base window
  sites <- data.frame(site = 1:max(mouse.exons$end_aa))
  sites$isSelected <- sapply(sites$site, \(i) nrow(positive.sites[positive.sites$site==i,"p"])>0)
  sites$isSelected.0.9 <- sapply(sites$site, \(i) ifelse(nrow(positive.sites[positive.sites$site==i,"p"])>0, positive.sites[positive.sites$site==i,"p"]>0.9, FALSE))
  
  
  sites %<>% dplyr::mutate(
    smoothSelected = slider::slide_dbl(isSelected, sum, .before=4, .after = 4),
    smoothSelected0.9 = slider::slide_dbl(isSelected.0.9, sum, .before=4, .after = 4),
  )
  
  sites.plot <- ggplot()+
    
    annotate("segment", x=0, xend=0, y=0, yend=8, linewidth=1)+
    new_scale_fill()+
    scale_fill_manual(values=c("white", "grey", "white", "grey", "white", "grey", "white"))+
    scale_pattern_color_manual(values=c("white", "white"))+
    
    scale_pattern_manual(values = c("none", "stripe")) + # which exons are patterned
    guides(fill = "none", pattern="none")+
    add.exon.track(-2, -1, col = "black")+ # color of border
    add.exon.labels(-2, -1)+
    
    add.track(ranges.NLS.common, -0.75, 7.25, fill=NLS.COLOUR, alpha = 1)+
    # Draw the structures
    add.track(ranges.ZF.common,    -0.5, 7,, fill=ZF.COLOUR)+
    add.track.labels(ranges.ZF.common, -0.5, 0, col="white", label_col = "motif_number")+
    
    add.track(ranges.9aaTAD.common, -0.5, 7, fill=TAD.COLOUR,  alpha = 0.9)+ # fill color max from "grDevices::Blues 3"
    add.track.labels(ranges.9aaTAD.common, -0.5, 0, col="white")+  # Label the 9aaTADs
    geom_line(data=sites,aes(x=site, y=smoothSelected0.9), linewidth=0.5 )+
    coord_cartesian(xlim = c(0, max(sites$site)))+
    scale_x_continuous(expand = c(0, 0), breaks = seq(0, 950, 50))+
    scale_y_continuous(breaks= c(0, 2, 4, 6, 8), labels = round(c(0, 2/9, 4/9, 6/9, 8/9), digits=2))+
    labs(y = "Fraction of bases under selection")+
    theme_bw()+
    theme(
      axis.text.y = element_text(size=7),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=6),
      # axis.ticks.y = element_blank(),
      axis.text.x = element_text(size=6),
      legend.position = "top",
      legend.title = element_text(size = 6, vjust = 0.85),
      legend.text = element_text(size = 6),
      legend.key.height = unit(3, "mm"),
      legend.spacing.y = unit(2, "mm"),
      legend.box.spacing = unit(2, "mm"),
      panel.border = element_blank(),
      axis.line.x.bottom = element_line(),
      panel.grid = element_blank())
  
  save.double.width("figure/positive.sites.filt.png", sites.plot, height = 35)
}

#### What are the hydrophobic patches in exons 2, 3 and 5 that are not 9aaTADs? ####

exon2.patch.start <- mouse.exons$start_aa[2]+95
exon2.patch.end   <- mouse.exons$start_aa[2]+118
exon2.patch <- msa.aa.aln.hydrophobicity %>%
  dplyr::filter(position_gapped>exon2.patch.start & position_gapped<exon2.patch.end)

exon2.patch.table <- do.call(rbind, 
                             lapply(METADATA$combined$common.name,  
                                    \(i) list("Sequence" = i, 
                                              as.character(ALIGNMENTS$aa.combined.biostrings@unmasked[[i]][exon2.patch.start:exon2.patch.end])))) %>%
  as.data.frame %>%
  dplyr::mutate(Sequence = factor(Sequence, levels = combined.taxa.name.order)) %>%
  dplyr::arrange(as.integer(Sequence)) 
colnames(exon2.patch.table) <- c("Sequence", paste0("exon_2_",exon2.patch.start,"-", 
                                                    exon2.patch.end))


exon3.patch.start <- mouse.exons$start_aa[3]+25
exon3.patch.end   <- mouse.exons$end_aa[3]
exon3.patch <- msa.aa.aln.hydrophobicity %>%
  dplyr::filter(position_gapped>exon3.patch.start & position_gapped<exon3.patch.end)

# Get the region from the msa

exon3.patch.table <- do.call(rbind, 
                             lapply(METADATA$combined$common.name,  
                                    \(i) list("Sequence" = i, 
                                              as.character(ALIGNMENTS$aa.combined.biostrings@unmasked[[i]][exon3.patch.start:exon3.patch.end])))) %>%
  as.data.frame %>%
  dplyr::mutate(Sequence = factor(Sequence, levels = combined.taxa.name.order)) %>%
  dplyr::arrange(as.integer(Sequence)) 
colnames(exon3.patch.table) <- c("Sequence", paste0("exon_3_",exon3.patch.start,"-", 
                                                    exon3.patch.end))

exon5.patch.start <- mouse.exons$start_aa[5]+28
exon5.patch.end   <- mouse.exons$end_aa[5]+1
exon5.patch <- msa.aa.aln.hydrophobicity %>%
  dplyr::filter(position_gapped>exon5.patch.start & position_gapped<exon5.patch.end)

exon5.patch.table <- do.call(rbind, lapply(METADATA$combined$common.name,  
                                           \(i) list("Sequence" = i,
                                                     as.character(ALIGNMENTS$aa.combined.biostrings@unmasked[[i]][exon5.patch.start:exon5.patch.end])))) %>%
  as.data.frame %>%
  dplyr::mutate(Sequence = factor(Sequence, levels = combined.taxa.name.order)) %>%
  dplyr::arrange(as.integer(Sequence))
colnames(exon5.patch.table)<-c("Sequence", paste0("exon_5_",exon5.patch.start,"-", 
                                                  exon5.patch.end))

exon.patch.table <- merge(exon2.patch.table, exon3.patch.table, by=c("Sequence")) %>%
  merge(., exon5.patch.table, by=c("Sequence")) %>%
  create.xlsx(., "figure/hydrophobic_patches.xlsx", cols.to.fixed.size.font = 2:4)





# Get consensus strings for the patches
# Note Biostrings::consensusString will fail if there are non-standard / ambiguity characters
try({
  # Get the consensus matrix and find the most frequent value per site
  exon.2.patch.consensus <- paste(unlist(apply(consensusMatrix(ALIGNMENTS$aa.combined.biostrings)[,(exon2.patch.start+1):(exon2.patch.end-1)], 
                                               2, 
                                               \(x) names(which(x==max(x))))), 
                                  collapse = "")
  
  exon.3.patch.consensus <- paste(unlist(apply(consensusMatrix(ALIGNMENTS$aa.combined.biostrings)[,(exon3.patch.start+1):(exon3.patch.end-1)], 
                                               2, 
                                               \(x) names(which(x==max(x))))), 
                                  collapse = "")
  
  exon.5.patch.consensus <- paste(unlist(apply(consensusMatrix(ALIGNMENTS$aa.combined.biostrings)[,(exon5.patch.start+1):(exon5.patch.end-1)], 
                                               2, 
                                               \(x) names(which(x==max(x))))), 
                                  collapse = "")
  
  cat("Exon 2 patch consensus: ", exon.2.patch.consensus, "\n")
  cat("Exon 3 patch consensus: ", exon.3.patch.consensus, "\n")
  cat("Exon 5 patch consensus: ", exon.5.patch.consensus, "\n")
})

exon2.hydro.plot <- ggplot()+
  geom_tile(data=exon2.patch,  aes(x = position_gapped, y = sequence, fill=hydrophobicity))+
  geom_text(data = exon2.patch, aes(x = position_gapped, y = sequence, label=character), size=2, family="mono", col="white")+
  geom_tile(data = sites, aes(x = site, y = 65, fill = as.integer(isSelected.0.9)))+
  scale_fill_paletteer_c("ggthemes::Classic Red-Blue", direction = -1, limits = c(0, 1))+
  scale_y_discrete(labels = gsub("_", " ", rev(combined.taxa.name.order)))+
  labs(fill="Hydrophobicity (per residue)")+
  coord_cartesian(xlim = c(exon2.patch.start,exon2.patch.end),  ylim = c(0, 65), clip = 'on')+
  annotate("text", label=s2c(exon.2.patch.consensus), size=2, x=(exon2.patch.start+1):(exon2.patch.end-1), 
           y=64.2, hjust=0.5, family="mono", fontface="bold")

exon3.hydro.plot <- ggplot()+
  geom_tile(data=exon3.patch,  aes(x = position_gapped, y = sequence, fill=hydrophobicity))+
  geom_text(data = exon3.patch, aes(x = position_gapped, y = sequence, label=character), size=2, family="mono", col="white")+
  scale_fill_paletteer_c("ggthemes::Classic Red-Blue", direction = -1, limits = c(0, 1))+
  scale_y_discrete(labels = gsub("_", " ", rev(combined.taxa.name.order)))+
  labs(fill="Hydrophobicity (per residue)")+
  coord_cartesian(xlim = c(exon3.patch.start,exon3.patch.end),  ylim = c(0, 63), clip = 'off')+
  annotate("text", label=s2c(exon.3.patch.consensus), size=2, x=(exon3.patch.start+1):(exon3.patch.end-1), 
           y=64.2, hjust=0.5, family="mono", fontface="bold")

exon5.hydro.plot <- ggplot()+
  geom_tile(data=exon5.patch,  aes(x = position_gapped, y = sequence, fill=hydrophobicity))+
  geom_text(data = exon5.patch, aes(x = position_gapped, y = sequence, label=character), size=2, family="mono", col="white")+
  scale_fill_paletteer_c("ggthemes::Classic Red-Blue", direction = -1, limits = c(0, 1))+
  scale_y_discrete(labels = gsub("_", " ", rev(combined.taxa.name.order)))+
  labs(fill="Hydrophobicity (per residue)")+
  coord_cartesian(xlim = c(exon5.patch.start,exon5.patch.end),  ylim = c(0, 63), clip = 'off')+
  annotate("text", label=s2c(exon.5.patch.consensus), size=2, x=(exon5.patch.start+1):(exon5.patch.end-1), 
           y=64.2, hjust=0.5, family="mono", fontface="bold")

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



#### Tar the outputs ####
system2("tar", "czf figure.tar.gz figure")
cat("Done!\n")
