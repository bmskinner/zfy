# Read codeml data and extract site values
source("functions.R")

#### codeml site model to check for site-specific selection ####
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
lrts[[1]]

# Compared with M1a, M2a adds a class of sites under positive selection with ω2
# > 1 (in proportion p2). This does not improve the fit of the model
# significantly
# (nearly neutral vs. positive selection)
lrts[[2]]

# Additional test for positive selection by comparing M7 (beta, null model)
# against M8 (beta&ω, alternative model).
# (positive selection vs null model)
# Evidence for sites under positive selection 
lrts[[3]]


#### codeml branch-site model to look for selection specifically in Muroidea ####

# Bayes Empirical Bayes (BEB) analysis (Yang, Wong & Nielsen 2005. Mol. Biol. Evol. 22:1107-1118)
# Extract sites under positive selection
positive.sites <- read_lines("paml/branch-site/branch.site.paml.out.txt") %>%
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

n.taxa <- 7
positive.sites.plot <-  positive.sites.plot+
  # Draw the conservation with Xenopus, chicken and opossum
  new_scale_fill()+
  scale_fill_viridis_c(limits = c(0, 1))+
  labs(fill="Conservation (5 site average)")+
  add.conservation.track(msa.aa.aln.tidy.frog.conservation,    n.taxa,   n.taxa+2)+
  add.conservation.track(msa.aa.aln.tidy.chicken.conservation, n.taxa+3, n.taxa+5)+
  add.conservation.track(msa.aa.aln.tidy.opossum.conservation, n.taxa+6, n.taxa+8)+
  
  # Draw the structures
  add.track(ranges.ZF.common,     n.taxa+9, n.taxa+11, fill="lightgrey")+
  add.track(ranges.NLS.common,    n.taxa+9, n.taxa+11, fill="green", alpha = 0.5)+
  add.track(ranges.9aaTAD.common, n.taxa+9, n.taxa+11, fill="#00366C",  alpha = 0.9)+ # fill color max from "grDevices::Blues 3"
  add.track.labels(ranges.9aaTAD.common, n.taxa+9, n.taxa+11)+   # Label the 9aaTADs
  
  new_scale_fill()+
  scale_fill_manual(values=c("white", "grey", "white", "grey", "white", "grey", "white"))+
  scale_pattern_color_manual(values=c("white", "white"))+
  
  scale_pattern_manual(values = c("none", "stripe")) + # which exons are patterned
  guides(fill = "none", pattern="none")+
  add.exon.track(n.taxa+12, n.taxa+14, col = "black")+ # color of border
  add.exon.labels(n.taxa+12, n.taxa+14)+
  
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
        panel.grid = element_blank())

save.double.width("figure/positive.sites.png", positive.sites.plot, height = 85)
