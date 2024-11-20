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

create.exon.plot <- function(name){
  exon.aln.file <- paste0("aln/exons/", name)

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


exon.plots <- lapply(c("exon_1-6.aln", "exon_2.aln", 
                       "exon_7.aln", "exon_1.3-6.aln"), 
                     create.exon.plot)

# Create a joint figure of exons 1-6, exon 2, and exon 7
exon.joint.tree <- exon.plots[[1]] + exon.plots[[2]] + exon.plots[[3]] + 
  patchwork::plot_annotation(tag_levels = list(c("Exons 1-6", "Exon 2", "Exon 7"))) &
  theme(plot.tag = element_text(size = 6),
        plot.margin = margin(t=0, l=0, r=0, b=0))
save.double.width("figure/Figure_Sxxxx_exons_tree.png", exon.joint.tree, height=120)

#### Test selection globally in mammals ####

kaks.pairwise.plot <- plot.kaks(FILES$mammal.nt.aln, 
                                species.order = mammal.taxa.name.order,
                                kaks.limits = c(0, 1.5))
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


#### What are the hydrophobic patches in exons 2, 3 and 5 that are not 9aaTADs? ####

# Extract all nt and aa sequences from the MSA for the given region, and calculate the
# consensus sequence
extract.alignment.region <- function(nt.start=NULL, nt.end=NULL, aa.start=NULL, aa.end=NULL){
  
  if(is.null(nt.start) & is.null(nt.end)){
    nt.start <- floor(max(1, aa.start*3)-2)
    nt.end <- ceiling(max(1, aa.end*3))
  }
  
  if(is.null(aa.start) & is.null(aa.end)){
    aa.start <- floor(max(1, nt.start/3))
    aa.end <- ceiling(max(1, nt.end/3))
  }
  
  nt.aln <- as.data.frame(as.matrix(ALIGNMENTS$nt.combined.biostrings)[,nt.start:nt.end]) %>%
    dplyr::mutate(Sequence = rownames(.),
                  Sequence = factor(Sequence, levels = combined.taxa.name.order)) %>%
    dplyr::arrange(as.integer(Sequence)) %>%
    dplyr::rename_with(.cols=starts_with("V"), .fn=\(x) paste0("Site_",as.integer(gsub("V", "", x))+nt.start-1) )
  
  aa.aln <- as.data.frame(as.matrix(ALIGNMENTS$aa.combined.biostrings)[,aa.start:aa.end]) %>%
    dplyr::mutate(Sequence = rownames(.),
                  Sequence = factor(Sequence, levels = combined.taxa.name.order)) %>%
    dplyr::arrange(as.integer(Sequence)) %>%
    dplyr::rename_with(.cols=starts_with("V"), .fn=\(x) paste0("Site_",as.integer(gsub("V", "", x))+aa.start-1) )
  
  
  # find the most common nucleotide in the consensus matrix
  nt.consensus <- paste(unlist(
    apply(consensusMatrix(ALIGNMENTS$nt.combined.biostrings)[,(nt.start):(nt.end)], 
          2, \(x) names( which(x==max(x)) )[1]  ) # Take only the first consensus in case of ties
  ), 
  collapse = "")
  
  aa.consensus <- paste(unlist(
    apply(consensusMatrix(ALIGNMENTS$aa.combined.biostrings)[,(aa.start):(aa.end)], 
          2, \(x) names( which(x==max(x)) )[1] ) # Take only the first consensus in case of ties
  ), 
  collapse = "")
  
  list(
    nt.aln   = nt.aln,
    aa.aln   = aa.aln,
    aa.start = aa.start,
    aa.end   = aa.end,
    aa.length = aa.end - aa.start + 1,
    nt.start = nt.start,
    nt.end   = nt.end,
    nt.length = nt.end - nt.start + 1,
    nt.consensus = nt.consensus,
    aa.consensus = aa.consensus
  )
}

plot.hydrophobic.patch <- function(patch.data, patch.start, patch.end, patch.consensus){
  ggplot()+
    geom_tile(data=patch.data,  aes(x = position_gapped, y = sequence, fill=hydrophobicity))+
    scale_fill_paletteer_c("ggthemes::Classic Red-Blue", direction = -1, limits = c(0, 1))+
    labs(fill="Hydrophobicity (per residue)")+
    geom_text(data = patch.data, aes(x = position_gapped, y = sequence, label=character), size=2, family="mono", col="white")+
    
    # Overlay diversifying selection from MEME
    new_scale_fill()+
    labs(fill="Diversifying selection")+
    geom_tile(data = meme.cols, aes(x = site, y = 65, fill = `p-value`<0.1))+
    
    scale_y_discrete(labels = gsub("_", " ", rev(combined.taxa.name.order)))+
    coord_cartesian(xlim = c(patch.start,patch.end),  ylim = c(0, 65), clip = 'on')+
    annotate("text", label=s2c(patch.consensus), size=2, x=(patch.start):(patch.start+nchar(patch.consensus)-1), 
             y=64.2, hjust=0.5, family="mono", fontface="bold")
}


exon2.patch <- extract.alignment.region(mouse.exons$start_aa[2]+95,
                                        mouse.exons$start_aa[2]+118)

exon3.patch <- extract.alignment.region(mouse.exons$start_aa[3]+25,
                                        mouse.exons$end_aa[3])

exon5.patch <- extract.alignment.region(mouse.exons$start_aa[5]+28,
                                        mouse.exons$end_aa[5]+1)

# Export to file
merge(exon2.patch$patch.table, exon3.patch$patch.table, by=c("Sequence")) %>%
  merge(., exon5.patch$patch.table, by=c("Sequence")) %>%
  create.xlsx(., "figure/hydrophobic_patches.xlsx", cols.to.fixed.size.font = 2:4)


exon2.hydro.plot <-plot.hydrophobic.patch(exon2.patch$patch, exon2.patch$start, exon2.patch$end, exon2.patch$consensus)
exon3.hydro.plot <-plot.hydrophobic.patch(exon3.patch$patch, exon3.patch$start, exon3.patch$end, exon3.patch$consensus)
exon5.hydro.plot <-plot.hydrophobic.patch(exon5.patch$patch, exon5.patch$start, exon5.patch$end, exon5.patch$consensus)

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
  
  cat(json.file, "K=", round(hyphy.data$`test results`$`relaxation or intensification parameter`, digits = 2), 
      "p=", round(hyphy.data$`test results`$`p-value`, digits = 2), "\n")
  
  p <- ggtree(hyphy.input.tree, size=1.5) + 
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

if(file.exists("aln/hyphy/combined.rodentia.relax.json")){
  
  node.names <- c("rodentia", "eumuroida", "muridae", "murinae")
  
  make.tree <- \(name){
    relax.tree <- create.relax.k.tree(paste0("aln/hyphy/combined.", name, ".relax.json"))
    save.double.width(paste0("figure/RELAX.combined.", name, ".png"), relax.tree)
    
    relax.mammal.tree <- create.relax.k.tree(paste0("aln/hyphy/mammal.", name, ".relax.json"))
    save.double.width(paste0("figure/RELAX.mammal.", name, ".png"), relax.mammal.tree)
  } 
  
  sapply(node.names, make.tree)
}

#### Plot HyPhy MEME test for directional selection ####


# Given an alignment region produced by `extract.alignment.region`, plot it
plot.alignment.region <- function(region.data, meme.overview){
  
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
                  clade = case_when(str_detect(Sequence, "Grass_Rat_ZFY") ~ "murinae",
                                    str_detect(Sequence, "Mouse_Zfy") ~ "murinae",
                                    str_detect(Sequence, "Rat_Zfy") ~ "murinae",
                                    str_detect(Sequence, "Mongolian_gerbil_Zfx-like") ~ "muridae",
                                    str_detect(Sequence, "deer_mouse_Zfx-like") ~ "eumuroida",
                                    str_detect(Sequence, "hamster_Zfx-like") ~ "eumuroida",
                                    str_detect(Sequence, "marmot_ZFY") ~ "rodentia",
                                    str_detect(Sequence, "squirrel_Zfx-like") ~ "rodentia",
                                    str_detect(Sequence, "squirrel_Zfy") ~ "rodentia",
                                    str_detect(Sequence, "Beaver_Zfx-like") ~ "rodentia",
                                    .default = "other"
                  ),
                  isSignificantSite = case_when( (clade %in% c("rodentia")) & rodentia<0.01 ~ T,
                                                 (clade %in% c("eumuroida")) & eumuroida<0.01 ~ T,
                                                 (clade %in% c("muridae")) & muridae<0.01 ~ T,
                                                 (clade %in% c("murinae")) & murinae<0.01 ~ T,
                                                 .default = F)
                  
    )
  
  
  
  
  p <- ggplot()+
    geom_tile(data=tidy.region.nt.data,  aes(x = Site, y = Sequence, fill=Base))+
    scale_fill_manual(values = c("G"="#EB413C", "C"="#FFB340", "T"="#3C88EE", "A"="#64F73F", "-"="grey"))+ # Jalview colours
    labs(fill="Base")+
    
    new_scale_fill()+
    geom_tile(data=tidy.region.aa.data,  aes(x = Site, y = Sequence, fill=Residue, width=0.5, height=1))+
    # geom_text(data = tidy.region.aa.data, aes(x = Site, y = Sequence, label=Residue), size=2, family="mono", col="black")+
    geom_text(data = tidy.region.aa.data, aes(x = Site, y = Sequence, label=Residue, col=isSignificantSite), 
              size=2, family="mono", )+
    scale_color_manual(values = c("white", "black"))+
    guides(col = FALSE)+
    guides(fill = FALSE)+
    
    # Overlay diversifying selection from MEME
    new_scale_fill()+
    labs(fill="Diversifying selection")+
    geom_tile(data = tidy.region.aa.data, aes(x = Site, y = 64.5, width=3, height=2, fill = isSignificant))+
    scale_fill_manual(values = c(`TRUE`="lightblue", `FALSE`="grey"))+
    guides(fill = FALSE)+

    scale_x_continuous(breaks = seq(region.data$nt.start, region.data$nt.end, 3))+
    coord_cartesian(xlim = c(region.data$nt.start,region.data$nt.end),  ylim = c(0, 66), clip = 'on')+
    # NT consensus
    annotate("text", label=s2c(region.data$nt.consensus), size=2, x=(region.data$nt.start):(region.data$nt.start+nchar(region.data$nt.consensus)-1),
             y=64.2, hjust=0.5, family="mono", fontface="bold")+
    # # AA consensus
    annotate("text", label=s2c(region.data$aa.consensus), size=2,
             x=seq((region.data$nt.start)+1,
                   (region.data$nt.start+nchar(region.data$nt.consensus)),
                   3),
             y=65.2, hjust=0.5, family="mono", fontface="bold")+
    theme_bw()
  
  list(plot = p,
       region = region.data)
}


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

    # Extract the sites covered by the partition
    meme.cols$site <- unlist(meme.data$`data partitions`[[as.name(partition)]]$coverage)
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
    dplyr::filter(`p-value`<0.01) %>%
    dplyr::mutate(start = site-1,
                  end = site+1)
  
  list(meme.data = meme.data,
       meme.sites = meme.sites)
}


combined.meme.results <- lapply(node.names, read.meme.results)

mammal.meme.results <- lapply(node.names, read.meme.results, outgroup.type="mammal")

# Create plots of the MSA around the sites of interest
plot.individual.meme.sites <- function(meme.results, outgroup.type){
  
  create.meme.overview.plot <- function(meme.overview){
    
    p.value.threshold <- 0.01
    
    # Annotate MEME sites over the exons
    meme.location.plot <- ggplot()+
      # Draw the structures
      add.track(ranges.NLS.common,    1.95, 2.55, fill=NLS.COLOUR, alpha = 1)+ # +8.5
      add.track(ranges.ZF.common,     2, 2.5, fill=ZF.COLOUR)+ # 9 - 11
      add.track.labels(ranges.ZF.common, 2, 2.5, col="white", label_col = "motif_number")+   # Label the ZFs
      add.track(ranges.9aaTAD.common, 2, 2.5, fill=TAD.COLOUR,  alpha = 0.9)+  #9
      add.track.labels(ranges.9aaTAD.common, 2, 2.5, col="white")+   # Label the 9aaTADs
      
      annotate("rect", xmin = exon2.patch$start, xmax = exon2.patch$end, 
               ymin = 2, ymax = 2.5, fill = "pink")+
      annotate("rect", xmin = exon3.patch$start, xmax = exon3.patch$end, 
               ymin = 2, ymax = 2.5, fill = "pink")+
      annotate("rect", xmin = exon5.patch$start, xmax = exon5.patch$end, 
               ymin = 2, ymax = 2.5, fill = "pink")+
      
      new_scale_fill()+
      scale_fill_manual(values=c("white", "grey", "white", "grey", "white", "grey", "white"))+
      scale_pattern_color_manual(values=c("white", "white"))+
      
      scale_pattern_manual(values = c("none", "stripe")) + # which exons are patterned
      guides(fill = "none", pattern="none")+
      add.exon.track(y.start = 2.6, y.end = 3, col="black")+
      add.exon.labels(y.start = 2.6, y.end = 3)+
      
      # Show which test clades the sites show up in
      geom_point(data = meme.overview, aes(x=site, y = 3.1, alpha = rodentia<p.value.threshold), col="black",size = 0.2)+
      geom_point(data = meme.overview, aes(x=site, y = 3.2, alpha = eumuroida<p.value.threshold), col="black",size = 0.2)+
      geom_point(data = meme.overview, aes(x=site, y = 3.3, alpha = muridae<p.value.threshold), col="black",size = 0.2)+
      geom_point(data = meme.overview, aes(x=site, y = 3.4, alpha = murinae<p.value.threshold), col="black", size = 0.2)+
      scale_alpha_manual(values=c(0, 1))+
      
      geom_tile(data = meme.overview, aes(x=site, y = 2.5, height = 1, width=1), fill = "brown")+
      annotate("text", x= 750, y=3.1, label = "rodentia", hjust=0, size = 2)+
      annotate("text", x= 750, y=3.2, label = "eumuroida", hjust=0, size = 2)+
      annotate("text", x= 750, y=3.3, label = "muridae", hjust=0, size = 2)+
      annotate("text", x= 750, y=3.4, label = "murinae", hjust=0, size = 2)+
      
      
      scale_x_continuous(breaks = seq(0, 885, 50), expand = c(0, 0.05))+
      coord_cartesian(xlim = c(0, 885))+
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
    
    save.double.width(filename = paste0("figure/meme.",outgroup.type,".site.locations.png"), meme.location.plot, height = 45)
  }
  
  # Find the unique set of sites in the meme data
  meme.overview <- do.call(rbind, lapply(meme.results, \(x) x$meme.data)) %>%
    tidyr::pivot_wider(id_cols = c(site, partition), names_from = node, values_from = `p-value`) %>%
    dplyr::filter(if_any(.cols= rodentia:murinae, ~.<0.01)) # filter to sites significant in any node
  
  create.meme.overview.plot(meme.overview)
  
  site.ranges <- with(meme.overview, IRanges(start=site-2, end=site+2))
  collapsed.ranges <- as.data.frame(IRanges::reduce(site.ranges))
  unique.regions <- mapply(extract.alignment.region, aa.start = collapsed.ranges$start+1, 
                           aa.end =collapsed.ranges$end-1, SIMPLIFY = FALSE)
  
  site.plots <- lapply(unique.regions, plot.alignment.region, meme.overview=meme.overview)
  
  
  lapply(site.plots, \(x){ x$plot+theme(legend.position = "top",
                                        axis.text = element_text(size=6),
                                        axis.title = element_blank(),
                                        legend.text = element_text(size=6),
                                        legend.title = element_text(size=8))
    save.double.width(filename = paste0("figure/meme.", outgroup.type, ".", x$region$aa.start,"-",x$region$aa.end,".png"), plot = last_plot())
  } )
}

plot.individual.meme.sites(mammal.meme.results, "mammal")
plot.individual.meme.sites(combined.meme.results, "combined")

plot.meme.sites <- function(plots.to.save){
  meme.site.plots <- patchwork::wrap_plots( lapply(plots.to.save, \(x)x$plot), nrow = 1) + plot_layout(guides = "collect", axes = "collect") &
    theme_bw()+
    theme(legend.position = "none",
          axis.text = element_text(size=6),
          axis.title = element_blank(),
          legend.text = element_text(size=6),
          legend.title = element_text(size=8))
  save.double.width(filename = paste("figure/meme.site.plots.combined.png"),
                    meme.site.plots, height = 120)
}

plot.meme.sites(site.plots[1:4])
plot.meme.sites(site.plots[5:8])
plot.meme.sites(site.plots[9:11])

#### codeml site models to check for site-specific and branch-site selection ####

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

read.branch.site.codeml.output("rodentia")
read.branch.site.codeml.output("eumuroida")
read.branch.site.codeml.output("muridae")
read.branch.site.codeml.output("murinae")

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
  cat("M0 vs. M1a\n")
  cat("p=", lrts[[1]]$p.value, "lrt=", lrts[[1]]$lrt, "\n")
  
  # Compared with M1a, M2a adds a class of sites under positive selection with ω2
  # > 1 (in proportion p2). This does not improve the fit of the model
  # significantly
  # (nearly neutral vs. positive selection)
  cat("M1a, M2a\n")
  cat("p=", lrts[[2]]$p.value, "lrt=", lrts[[2]]$lrt, "\n")
  
  # Additional test for positive selection by comparing M7 (beta, null model)
  # against M8 (beta&ω, alternative model).
  # (positive selection vs null model)
  # Evidence for sites under positive selection 
  cat("M7, M8\n")
  cat("p=", lrts[[3]]$p.value, "lrt=", lrts[[3]]$lrt, "\n")
}

read.site.specific.codeml.output()





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
  
  # ZF* moves to sex chromosomes
  annotate("rect", xmin=0.02, ymin=7.8, xmax=0.04, ymax=9.2, fill="darkgreen", alpha=0.4)+
  annotate("text", x=0.02, y=8.6,label="ZF* to\nX/Y", size=2, hjust=0)+
  # Ssty box
  annotate("rect", xmin=0.12, ymin=25.8, xmax=0.19, ymax=27.5, fill="darkgreen", alpha=0.4)+
  annotate("text", x=0.13, y=27, label="Ssty appears", size=2, hjust=0)+
  # Zfy testis specific box
  annotate("text", x=0.13, y=26.2, label="Zfy testis specific", size=2, hjust=0)+
  annotate("rect", xmin=0.29, ymin=29.5, xmax=0.34, ymax=31, fill="darkgreen", alpha=0.4)+
  # Sly amplifies box
  annotate("text", x=0.295, y=30.5, label="Sly amplifies", size=2, hjust=0)+
  annotate("rect", xmin=0.24, ymin=27.4, xmax=0.27, ymax=29, fill="darkgreen", alpha=0.4)+
  # Slxl1 acquired box
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

# Make figures for presentations

# Scaled timetree view, no colours
time.plot.base <- ggtree(zfy.nt.aln.tree.time, size = 1, color="black") %<+%
  time.vals +
  geom_nodelab(size=2, nudge_x = -3, nudge_y = 0.5, hjust = 1, color = "black")+
  geom_tiplab(size=2, color = "black")+
  geom_treescale(fontsize =2, y = -1, width = 10) +
  coord_cartesian(xlim = c(-5, 210), ylim=c(-1, 32))+
  theme_tree() +
  theme(legend.position = "none",
        legend.background = element_blank(),
        legend.text = element_text(size=6),
        legend.title = element_text(size=6))
save.double.width("figure/subs.per.site.time.base.png", time.plot.base)

# Scaled timetree view, with colours
time.plot.colours <- ggtree(zfy.nt.aln.tree.time, size = 1) %<+%
  time.vals +
  aes(colour = log(subsPerMyr)) +
  scale_color_paletteer_c("ggthemes::Classic Red-Blue", 
                          direction = -1, limits =c(-9, -3))+
  labs(color = "Log substitutions per site\nper million years")+
  geom_nodelab(size=2, nudge_x = -3, nudge_y = 0.5, hjust = 1, color = "black")+
  geom_tiplab(size=2, color = "black")+
  geom_treescale(fontsize =2, y = -1, width = 10) +
  coord_cartesian(xlim = c(-5, 210), ylim=c(-1, 32))+
  theme_tree() +
  theme(legend.position = c(0.2, 0.8),
        legend.background = element_blank(),
        legend.text = element_text(size=6),
        legend.title = element_text(size=6))
save.double.width("figure/subs.per.site.time.colours.png", time.plot.colours)

# Annotate with labels
time.plot.annotated <- ggtree(zfy.nt.aln.tree.time, size = 1) %<+%
  time.vals +
  # ZF* moves to sex chromosomes
  annotate("rect", xmin=30, ymin=7.8, xmax=60, ymax=9.5, fill="darkgreen", alpha=0.4)+
  annotate("text", x=40, y=9,label="ZF* to X/Y", size=2, hjust=0)+
  # Ssty box
  annotate("rect", xmin=113, ymin=25.8, xmax=138, ymax=27.5, fill="darkgreen", alpha=0.4)+
  annotate("text", x=115, y=27, label="Ssty appears", size=2, hjust=0)+
  # Zfy testis specific
  annotate("text", x=115, y=26.2, label="Zfy testis specific", size=2, hjust=0)+
  # Sly amplifies box
  annotate("text", x=170, y=30.7, label="Sly\namplifies", size=2, hjust=0)+
  annotate("rect", xmin=170, ymin=29.8, xmax=180, ymax=31.5, fill="darkgreen", alpha=0.4)+
  # Slxl1 acquired box
  # annotate("text", x=0.242, y=28, label="Slxl1\nacquired?", size=2, hjust=0)+
  aes(colour = log(subsPerMyr)) +
  scale_color_paletteer_c("ggthemes::Classic Red-Blue", 
                          direction = -1, limits =c(-9, -3))+
  labs(color = "Log substitutions per site\nper million years")+
  geom_nodelab(size=2, nudge_x = -3, nudge_y = 0.5, hjust = 1, color = "black")+
  geom_tiplab(size=2, color = "black")+
  geom_treescale(fontsize =2, y = -1, width = 10) +
  coord_cartesian(xlim = c(-5, 210), ylim=c(-1, 32))+
  theme_tree() +
  theme(legend.position = c(0.2, 0.8),
        legend.background = element_blank(),
        legend.text = element_text(size=6),
        legend.title = element_text(size=6))
save.double.width("figure/subs.per.site.time.annotated.png", time.plot.annotated)

#### Tar the outputs ####
system2("tar", "czf figure.tar.gz figure")
cat("Done!\n")
