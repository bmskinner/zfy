# Calculate charge in a window surrounding each amino acid

# based on https://www.bioinformatics.nl/cgi-bin/emboss/help/charge

# aa - the amino acid sequence. May have gaps if part of a multiple sequence alignment
# sequence.name - the name of the sequence
# window.size - the size of the window (must be odd number)
# Returns - a data.frame with each character, their position in the sequence, the
# sequence name and the charge  If the input contains gaps ("-"),
# the coordinates will be with respect to the gapped alignment.
calc.charge <- function(aa, sequence.name, window.size){
  if(!require(dplyr)) stop("dplyr is required")
  if(!require(seqinr)) stop("seqinr is required")
  if(!require(slider)) stop("slider is required")
  if(window.size%%2==0) stop("Window size must be odd")
  
  # Convert ungapped coordinates back to gapped
  # site.no.gap - the integer site in an ungapped sequence to convert
  # gapped.seq - the sequence with gaps from an alignment
  convert.to.gapped.coordinate <- function(site.no.gap, gapped.seq){
    
    gapped.seq.char <- as.character(gapped.seq)
    # find gaps and stop codons
    gaps <- str_locate_all(gapped.seq.char, "-|\\*")[[1]][,1]
    n <- site.no.gap
    for(i in gaps){
      if(i<n) n <- n + 1
    }
    n
  }
  
  # Ensure this is character data in case input was from an MSA object
  aa.char <- as.character(aa)

  n.extend <- (window.size-1)/2 # number of sites ahead and behind
  
  # remove gaps and stop codons
  aa.char.nogap <- str_replace_all(aa.char, "-|\\*", "")
  
  #  'D' and 'E' a charge of -1, 'K' and 'R' a charge of +1, and the residue 'H' a charge of +0.5
  data.frame("character" = seqinr::s2c(aa.char.nogap)) %>%
    dplyr::mutate(position_ungapped = row_number(),
      charge = case_when(character=="D" ~ -1,
                     character=="K" ~ 1,
                     character=="R" ~ 1,
                     character=="H" ~ 0.5,
                     .default = 0),
                  
                  # Perform smoothing
                  charge_smoothed = slider::slide_dbl(charge,  mean, .before=n.extend, .after = n.extend)) %>%
    dplyr::rowwise() %>%
    # Now account for gaps by offsetting indexes back
    dplyr::mutate(position_gapped = convert.to.gapped.coordinate(position_ungapped, aa),
                  sequence = sequence.name) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(position_gapped)
}
