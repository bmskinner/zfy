# Calculate hydrophobicities in a window surrounding each amino acid
# Hydrophobicities via Anal. Biochem. 193:72-82(1991). 

# aa - the amino acid sequence. May have gaps if part of a multiple sequence alignment
# sequence.name - the name of the sequence
# window.size - the size of the window (must be odd number)
# Returns - a data.frame with each character, their position in the sequence, the
# sequence name and the hydrophobicity. If the input contains gaps ("-"),
# the coordinates will be with respect to the gapped alignment.
calc.hydrophobicity <- function(aa, sequence.name, window.size){
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
  
  # aa.chemistries <- data.frame(
  #   "names" = ,
  #   "hydrophobicity" = c(, , , , , , , ,
  #                         , , , , ,, , ,
  #                         , , , )
  # )
  
  aa.chemistries <- c("A"=0.616, "R"=0.000, "N"=0.236, "D"=0.028, "C"=0.680, "Q"=0.251, 
                         "E"=0.043, "G"=0.501, "H"=0.165, "I"=0.943, "L"=0.943, "K"=0.283,
                         "M"=0.738, "F"= 1.000, "P"=0.711, "S"=0.359, "T"=0.450, "W"=0.878, 
                         "Y"=0.880, "V"=0.825)
  
  # Ensure this is character data in case input was from an MSA object
  aa.char <- as.character(aa)
  
  n.extend <- (window.size-1)/2 # number of sites ahead and behind
  
  # remove gaps and stop codons
  aa.char.nogap <- str_replace_all(aa.char, "-|\\*", "")
  
  #  'D' and 'E' a charge of -1, 'K' and 'R' a charge of +1, and the residue 'H' a charge of +0.5
  data.frame("character" = seqinr::s2c(aa.char.nogap)) %>%
    dplyr::mutate(position_ungapped = row_number(),
                  hydrophobicity = aa.chemistries[character],
                  
                  # Perform smoothing
                  hydrophobicity_smoothed = slider::slide_dbl(hydrophobicity,  mean, .before=n.extend, .after = n.extend)) %>%
    dplyr::rowwise() %>%
    # Now account for gaps by offsetting indexes back
    dplyr::mutate(position_gapped = convert.to.gapped.coordinate(position_ungapped, aa),
                  sequence = sequence.name) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(position_gapped)
}
