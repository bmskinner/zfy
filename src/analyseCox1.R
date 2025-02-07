# Analyse COX1
# We want to compare the evolutionary rate of an MT gene with 0.5 ploidy versus our Y genes

#### Imports #####

source("src/functions.R")
load.packages()

# Clear previous analyses
filesstrings::dir.remove("aln/cox1")
filesstrings::create_dir("aln/cox1")

#### Align sequences #####

# Read all unaligned sequence files with .fa extension
fa.files <- list.files(path = "fasta/COX1", pattern = "*.fa$", 
                       include.dirs = T, full.names = T)

# Combine to single FASTA file
fa.read  <- lapply(fa.files, ape::read.FASTA)
seq.names <- sapply(fa.read, names)
binomial.names <- gsub("_", " ", str_extract(fa.files, "fasta/COX1/(.*)\\.fa", group = 1))
names.map <- METADATA$mammal %>% 
  dplyr::select(species, Species_common_name) %>%
  dplyr::distinct()

names.map$species <- gsub(" texanus", "", names.map$species) # deer

common.names <- names.map[names.map$species==binomial.names,"Species_common_name"]

# Remove existing names from sequences
fa.read <- do.call(c, lapply(fa.read, function(x) {
  names(x)<-NA
  return(x)
}))

# Add common names to match species phylogeny
names(fa.read) <- common.names
ape::write.FASTA(fa.read, file = FILES$cox1.nt.fas)

# Use macse to align with mitochondrial code
run.macse(FILES$cox1.nt.fas, "aln/cox1/cox1",
          "-gc_def 2")  # vertebrate mito code

# system2("java", paste("-jar bin/macse_v2.07.jar -prog alignSequences",
#                       "-gc_def 2",
#                       "-seq",    FILES$cox1.nt.fas, # input
#                       "-out_NT", FILES$cox1.nt.aln,  # output nt alignment
#                       "-out_AA", FILES$cox1.aa.aln), # output aa alignment
#         stdout = paste0(FILES$cox1.nt.aln, ".macse.log"),  # logs
#         stderr = paste0(FILES$cox1.nt.aln, ".macse.log"))  # error logs

#### Make ML tree with distances #####

# Specify the true species phylogeny

cox1.phylogeny <- paste0("(Platypus, (Opossum, ", # Outgroups
                          "(", # Eutheria
                            "(Southern_two-toed_sloth, African_bush_elephant)Atlantogenata, ", # Afrotheria & Xenarthra
                            "(", # Boreoeutheria
                            "(",  # Euarchonoglires
                            "( ", # Simiiformes
                              "Common_marmoset,", # New world monkeys
                              "(", # Catarrhini (Old world monkeys & apes)
                              "(", #Cercopithecidae (Old world monkeys)
                                "Golden_snub-nosed_monkey,",   #Colobinae
                                "(Olive_baboon, Macaque)Cercopithecinae", # Cercopithecinae
                                ")Cercopithecidae,", # /Cercopithecidae (Old world monkeys)
                              "(", # Hominidae
                              "Gorilla, (Chimpanzee, Human)Hominini",
                              ")Hominidae",  # /Hominidae
                              ")Catarrhini", # /Catarrhini
                            ")Simiiformes,",  # /Simiiformes
                            "(", # Rodentia
                              "(Gray_squirrel,(Arctic_ground_squirrel, Alpine_marmot)Xerinae)Sciuridae,", # Sciuridae 
                             "(Damara_mole-rat,", 
                                "(Beaver, (", # Muroidea
                                "(North_American_deer_mouse, Desert_hamster)Cricetidae,", # Cricetidae 
                                "(Mongolian_gerbil, (Rat, (Mouse, African_Grass_Rat)Mus-Arvicanthis)Murinae)Muridae", # Muridae 
                              ")Eumuroida)Muroidea)Muroidea-Fukomys", # /Muroidea
                            ")Rodentia", # /Rodentia
                            ")Euarchonoglires,", # /Euarchonoglires
                            "(", # Laurasiatheria
                              "(",  # Carnivora
                              "Cat, (Dog, (Stoat, Polar_bear)Arctoidea)Caniformia",
                              ")Carnivora,",  # /Carnivora
                              "(",  # Euungulata 
                              "Horse, (Pig, (White_tailed_deer, (Cattle, Goat)Bovidae)Pecora)Artiodactyla",
                              ")Euungulata", # /Euungulata 
                            ")Laurasiatheria", # /Laurasiatheria
                            ")Boreoeutheria", # /Boreoeutheria
                          ")Eutheria", # /Eutheria
                        ")Theria)Mammalia;") # /Outgroups

write_file(cox1.phylogeny, "aln/cox1/cox1.nt.species.nwk")


# Make the tree using the species phylogeny
system2("iqtree", paste("-s ", "aln/cox1/cox1.nt.aln", 
                        "-nt AUTO", # number of threads
                        "-te aln/cox1/cox1.nt.species.nwk", # user tree guide
                        "-st CODON",
                        "-asr")) # ancestral sequence reconstruction

#### Plot the tree ####

cox1.nt.aln.tree <- ape::read.tree("aln/cox1/cox1.nt.aln.treefile")
# Root the trees on platypus
cox1.nt.aln.tree <- phytools::reroot(cox1.nt.aln.tree, which(cox1.nt.aln.tree$tip.label=="Platypus"), position = 0.015)

cox1.nt.aln.tree.plot <- plot.tree(cox1.nt.aln.tree) + 
  geom_nodelab(size=2, nudge_x = -0.01, nudge_y = 0.5, hjust=1,  node = "internal") + 
  xlim(0, 1.5) + 
  labs(title = "COX1 (MT)")
save.double.width("figure/cox1.tree.png", cox1.nt.aln.tree.plot)


#### Test selection on COX1 ####

cox1.taxa.name.order <- get_taxa_name(cox1.nt.aln.tree.plot) 

cox1.kaks.plot <- plot.kaks(FILES$cox1.nt.aln, 
          species.order = cox1.taxa.name.order,
          kaks.limits = c(0, 1.5))


save.double.width("figure/cox1.dnds.png", cox1.kaks.plot)
