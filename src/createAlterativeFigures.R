# Modified figures for presentations
source("src/functions.R")
load.packages()


# dnds plot for presentations - rotate so it matches the tree
seqin.aln <- seqinr::read.alignment(FILES$mammal.nt.aln, format = "fasta")
kaks.data <- seqinr::kaks(seqin.aln)

kaks.ratio <- kaks.data$ka / kaks.data$ks

kaks.pairwise.rotated <- metagMisc::dist2list(kaks.ratio, tri = F) %>%
  dplyr::mutate(col = str_replace_all(col, "_", " "),
                row = str_replace_all(row, "_", " "),
                col = factor(col, levels = mammal.taxa.name.order),
                row = factor(row, levels = mammal.taxa.name.order),
                colnum = as.integer(col),
                rownum = as.integer(row)) %>%
  dplyr::filter(rownum < colnum)

kaks.pairwise.plot.rotated <- ggplot(kaks.pairwise.rotated, aes(x = col, y = row))+
  geom_tile(aes(fill=value))+
  scale_fill_viridis_c(limits = c(0, 1), direction = -1)+
  labs(fill="dNdS")+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        legend.background = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8))
save.double.width("figure/dnds.undecorated.png", kaks.pairwise.plot.rotated)