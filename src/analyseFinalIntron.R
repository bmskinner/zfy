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
cat("Aligning sequences\n")
run.muscle(FILES$final.intron.zfy.nt.fas, FILES$final.intron.zfy.nt.aln)
run.muscle(FILES$final.intron.zfx.nt.fas, FILES$final.intron.zfx.nt.aln)
run.muscle(FILES$final.intron.nt.fas,     FILES$final.intron.nt.aln)

#### Run divvier to improve high-confidence homologies in alignments ####

final.intron.zfy.nt.divvy.aln <- run.divvier(FILES$final.intron.zfy.nt.aln)
final.intron.zfx.nt.divvy.aln <- run.divvier(FILES$final.intron.zfx.nt.aln)

#### Make species phylogenies for Zfx and Zfy #####

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

#### Make ML tree based on raw alignments #####

cat("Creating raw alignment trees\n")
# Make the tree using the species phylogeny
final.intron.zfy.nt.aln.tree <- run.iqtree(FILES$final.intron.zfy.nt.aln, 
                        "-nt AUTO", # number of threads
                        "-keep-ident",
                        "-te aln/final.intron.zfy/zfy.nt.species.nwk" # user tree guide
                        ) %>%
        ape::read.tree(.) %>%
        reroot.tree(., c("African_bush_elephant_ZFY"), position = 0.015)

final.intron.zfx.nt.aln.tree <- run.iqtree(FILES$final.intron.zfx.nt.aln, 
                        "-nt AUTO", # number of threads
                        "-keep-ident",
                        "-te aln/final.intron.zfx/zfx.nt.species.nwk" # user tree guide
                        ) %>%
        ape::read.tree(.) %>%
        reroot.tree(., c("African_bush_elephant_ZFX"), position = 0.015)

#### Make ML tree based on divvied alignments ####
cat("Creating divvied alignment trees\n")
final.intron.zfy.nt.divvy.aln.tree <- run.iqtree(final.intron.zfy.nt.divvy.aln, 
                                           "-nt AUTO", # number of threads
                                           "-keep-ident",
                                           "-te aln/final.intron.zfy/zfy.nt.species.nwk" # user tree guide
) %>%
  ape::read.tree(.) %>%
  reroot.tree(., c("African_bush_elephant_ZFY"), position = 0.015)

final.intron.zfx.nt.divvy.aln.tree <- run.iqtree(final.intron.zfx.nt.divvy.aln, 
                                           "-nt AUTO", # number of threads
                                           "-keep-ident",
                                           "-te aln/final.intron.zfx/zfx.nt.species.nwk" # user tree guide
) %>%
  ape::read.tree(.) %>%
  reroot.tree(.,"African_bush_elephant_ZFX", position = 0.015)
#### Plot the trees ####

final.intron.zfy.nt.aln.tree.plot <- plot.tree(final.intron.zfy.nt.aln.tree)  + xlim(0, 2) + labs(title = "ZFY")
final.intron.zfx.nt.aln.tree.plot <- plot.tree(final.intron.zfx.nt.aln.tree) + xlim(0, 2)+ labs(title = "ZFX")
save.double.width("figure/final.intron.tree.png", final.intron.zfx.nt.aln.tree.plot/final.intron.zfy.nt.aln.tree.plot)

final.intron.zfy.nt.divvy.aln.tree.plot <- plot.tree(final.intron.zfy.nt.divvy.aln.tree)  + xlim(0, 2) + labs(title = "ZFY")
final.intron.zfx.nt.divvy.aln.tree.plot <- plot.tree(final.intron.zfx.nt.divvy.aln.tree) + xlim(0, 2)+ labs(title = "ZFX")
save.double.width("figure/final.intron.divvy.tree.png", final.intron.zfy.nt.divvy.aln.tree.plot/final.intron.zfx.nt.divvy.aln.tree.plot)


# final.intron.nt.aln.tree.plot <- plot.tree(final.intron.nt.aln.tree)  + xlim(0, 3) + labs(title = "Combined")
# save.double.width("figure/final.intron.combined.tree.png", final.intron.nt.aln.tree.plot)

#### End ####
cat("Done!")
