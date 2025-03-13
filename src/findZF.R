# Find Zinc Fingers in a sequence
# Process ZF locations identified by hmmsearch

# aa - the amino acid sequence. May have gaps if part of a multiple sequence alignment
# sequence.name - the name of the sequence
# Returns - a tibble with ZFs, their coordinates in the sequence, and the
# sequence name. If the input contains gaps ("-"),
# the coordinates will be with respect to the gapped alignment.
find.zf <- function(aa, sequence.name, start, end){
  if(!require(dplyr)) stop("dplyr is required")

  aa <- aa@unmasked[[sequence.name]]
  
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
  
  # remove gaps and stop codons
  aa.char.nogap <- str_replace_all(aa.char, "-|\\*", "")
  
  # Get the location of the hits in the ungapped sequence and correct for
  # alignment gaps
  data.frame("sequence" = sequence.name, "start" = start, "end" = end) %>%
    dplyr::rowwise() %>%
    dplyr::rename(start_ungapped = start,
                  end_ungapped   = end) %>%
    # Now account for gaps by offsetting indexes back
    dplyr::mutate(start_gapped = convert.to.gapped.coordinate(start_ungapped, aa),
                  end_gapped   = convert.to.gapped.coordinate(end_ungapped, aa),
                  start_nt_ungapped = start_ungapped * 3,
                  end_nt_ungapped = end_ungapped * 3,
                  sequence = sequence.name) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(start_ungapped)
}
