# Function to find 9aaTADs

# Matches the outputs of https://www.med.muni.cz/9aaTAD/analysis.php 
#Piskacek S, Gregor M, Nemethova M, et al (2007) Nine-amino-acid transactivation
#domain: Establishment and prediction utilities. Genomics 89:756–768.
#https://doi.org/10.1016/j.ygeno.2007.02.003


# The 9aa TAD can be detected is [MDENQSTYG]{KRHCGP}[ILVFWM]{KRHCGP}{CGP}{CGP}[ILVFWM]{CGP}{CGP}
# But there are refinements to this, and the aas can vary; 
# Approved positions
# 1	2	3	4	5	6	7	8	9
# RC1	min 1 x [ILVFWMAY]	: 2, 4, 5
# RC2	min 1 x [ILVFWMAY]	: 6, 8
# RC3	min 1 x [DENQST]	 	: 4, 5, 6
# RC4	min 2 - max 4 x [DENQSTKRH]	 : 1 - 9
# RC5	min 5 - max 7 x [ILVFWMAYGPST] : 1 - 9
# RC6	max 3 x [ILVFWMAYGP] in a row	 : 1 - 9
# RC7	max 2 x [DENQST] in a row	 : 1 - 9
# RC8	max 1 x [KRH]	 : 1 - 9
# RC9	max 1 x [QNKRH]	 	 : 1 - 6
# RC10	max 1 x [QNKRH]	 	:8, 9
# RC11	any [C]	 	 	 : 2 - 8
# RC12   	any [PG]	 : 2 - 7	 	

#' Find 9aaTADs in an amino acid sequence
#'
#' @param aa the amino acid sequence (character or Biostrings alignment). May 
#' have gaps if part of a multiple sequence alignment
#' @param sequence.name  the name of the sequence
#' @param rc.threshold threshold for percent of refinement criteria matching
#'
#' @return a tibble with 9aaTADs, their coordinates in the sequence, the
# sequence name, and the rc score (as a %). If the input contains gaps ("-"),
# the coordinates will be with respect to the gapped alignment.
#' @export
#'
#' @examples
find.9aaTAD <- function(aa, sequence.name, rc.threshold){
  if(!require(dplyr)) stop("dplyr is required")
  if(!require(seqinr)) stop("seqinr is required")
  
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
  
   # make all 9aa subseqences for testing
  starts <- 1:(nchar(aa.char.nogap)-9)
  ends <- starts+8
  
  mat <- matrix(c(starts, ends), byrow=F, ncol = 2)
  colnames(mat)<-c("start", "end")
  
  # Break input into substrings
  hits <- str_sub_all(aa.char.nogap, mat)[[1]]
  
  # All 9aaTADs must meet these criteria
  essential.regex <- "[MDENQSTYG][^KRHCGP][ILVFWM][^KRHCGP][^CGP][^CGP][ILVFWM][^CGP][^CGP]"
  
  # Find the initial candidates
  candidates <- hits[grepl(essential.regex, hits)]
  
  # Define the tests for refining candidates
  test.rcs <- function(s){
    
    s <- s2c(s) # make character vector from string
    
    # Test each substring against criteria
    rc1 <- \(s) sum(s[c(2, 4, 5)] %in% s2c("ILVFWMAY"))>=1
    rc2 <- \(s) sum(s[c(6, 8)] %in% s2c("ILVFWMAY"))>=1
    rc3 <- \(s) sum(s[c(4, 5, 6)] %in% s2c("DENQST"))>=1
    rc4 <- \(s) {count <- sum(s[1:9] %in% s2c("DENQSTKRH")); count>=2 & count<=4}
    rc5 <- \(s) {count <- sum(s[1:9] %in% s2c("ILVFWMAYGPST")); count>=5 & count<=7}
    rc6 <- \(s) { 
      m <- s[1:9] %in% s2c("ILVFWMAYGP");
      !any( sapply(3:8, \(i) { all(m[i-2], m[i-1], m[i], m[i+1]) })) &
        !any( sapply(2:7, \(i) { all(m[i-1], m[i], m[i+1], m[i+2]) }))
    }
    rc7 <- \(s) { 
      m <- s[1:9] %in% s2c("DENQST");
      !any(sapply(2:8, \(i) all(m[i-1], m[i], m[i+1])))
    }
    rc8  <- \(s) sum(s[1:9] %in% s2c("KRH"))<=1
    rc9  <- \(s) sum(s[1:6] %in% s2c("QNKRH"))<=1
    rc10 <- \(s) sum(s[8:9] %in% s2c("QNKRH"))<=1
    rc11 <- \(s) sum(s[2:8] %in% s2c("C"))>=0
    rc12 <- \(s) sum(s[2:7] %in% s2c("PG"))>=0
    
    
    (rc1(s) + rc2(s) + rc3(s) + rc4(s) + rc5(s) + rc6(s) + rc7(s) + 
       rc8(s) + rc9(s) + rc10(s) + rc11(s) + rc12(s))/12*100
  }

  rc.results <- sapply(candidates, test.rcs) # score refinement criteria
  
  # Get the location of the hits in the original sequence and correct for
  # alignment gaps
  str_locate(aa.char.nogap,candidates[rc.results>=rc.threshold] )  %>%
    as.data.frame %>%
    dplyr::mutate("hit" = candidates[rc.results>=rc.threshold],
                  "rc_score" = rc.results[rc.results>=rc.threshold]) %>%
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


#' Combine overlapping 9aaTADs into 'superTAD' regions
#' 
#' For 9aaTADs, we want to combine overlapping 9aaTADs into a single wider region,
#' but only for those TADs that are in 5 or more species (otherwise we will collapse
#' separate interesting regions into one superTAD). Use overall coverage along the sequence
#' to separate the regions
#'
#' @param locations.9aaTAD tibble output from find.9aaTAD
#' @param rc.threshold 9aaTADs must have an RC score equal to or highher than 
#' this to be combined
#' @param coverage.threshold site must be in this many sequences or more to be
#' included
#' 
#' @return data frame containing the locations of the superTADs at the desired
#' refinement criteria threshold
#' @export
#'
#' @examples
find.common.aa.9aaTADs <- function(locations.9aaTAD, rc.threshold = 80, coverage.threshold = 21){
  
  # We only want to use the high confidence 9aaTADs, and ignore the site in exon 7
  tads.filtered <- locations.9aaTAD[locations.9aaTAD$rc_score>=rc.threshold & locations.9aaTAD$end_gapped < mouse.exons$start_aa[7],]
  
  # Make ranges from all of the TADs with >80% RC
  tad.ranges <- IRanges(start = tads.filtered$start_gapped, end = tads.filtered$end_gapped, names = tads.filtered$sequence)
  
  # For each base in the gapped aa alignment, count overlapping TADs
  count.overlapping.tads <- function(site){
    site.range <- IRanges(start = site, end=site)
    IRanges::countOverlaps(site.range, tad.ranges)
  }
  
  tad.coverages <- data.frame("site" = 1:max(tads.filtered$end_gapped))
  tad.coverages$coverage <- sapply(tad.coverages$site, count.overlapping.tads)
  tad.coverages <- tad.coverages[tad.coverages$coverage>=coverage.threshold,]
  
  # For high-coverage 9aaTADs, reduce to overlapping ranges
  tad.coverages.ranges <- IRanges(IRanges(start =  tad.coverages$site, end =  tad.coverages$site))
  reduced.ranges<- IRanges::reduce(tad.coverages.ranges)
  as.data.frame(reduced.ranges) %>%
    dplyr::mutate(motif_number = row_number(),
                  adj.motif.number = ifelse(motif_number>4, motif_number-1, motif_number),
                  # For labels, we want D1, D2, E ..., not D, E, F ...
                  # label = case_when(motif_number==4 ~ paste0(LETTERS[adj.motif.number], 1),
                  #                   motif_number==5 ~ paste0(LETTERS[adj.motif.number], 2),
                  #                   .default = LETTERS[adj.motif.number] )) %>%
                  label = LETTERS[motif_number]) %>%
    
    dplyr::select(-adj.motif.number) 
}

find.common.nt.9aaTADs <- function(locations.9aaTAD, rc.threshold = 80, coverage.threshold = 21){
  
  locations.9aaTAD %<>%
    dplyr::filter(!is.na(start_nt_gapped) & !is.na(end_nt_gapped))
  
  # We only want to use the high confidence 9aaTADs, and ignore the site in exon 7
  tads.filtered <- locations.9aaTAD[locations.9aaTAD$rc_score>=rc.threshold & locations.9aaTAD$end_nt_gapped < mouse.exons$start_nt[7],]
  
  # Make ranges from all of the TADs with >80% RC
  tad.ranges <- IRanges(start = tads.filtered$start_nt_gapped, end = tads.filtered$end_nt_gapped, names = tads.filtered$sequence)
  
  # For each base in the gapped aa alignment, count overlapping TADs
  count.overlapping.tads <- function(site){
    site.range <- IRanges(start = site, end=site)
    IRanges::countOverlaps(site.range, tad.ranges)
  }
  
  tad.coverages <- data.frame("site" = 1:max(tads.filtered$end_nt_gapped))
  tad.coverages$coverage <- sapply(tad.coverages$site, count.overlapping.tads)
  tad.coverages <- tad.coverages[tad.coverages$coverage>=coverage.threshold,]
  
  # For high-coverage 9aaTADs, reduce to overlapping ranges
  tad.coverages.ranges <- IRanges(IRanges(start =  tad.coverages$site, end =  tad.coverages$site))
  reduced.ranges<- IRanges::reduce(tad.coverages.ranges)
  as.data.frame(reduced.ranges) %>%
    dplyr::mutate(motif_number = row_number(),
                  adj.motif.number = ifelse(motif_number>4, motif_number-1, motif_number),
                  # For labels, we want D1, D2, E ..., not D, E, F ...
                  # label = case_when(motif_number==4 ~ paste0(LETTERS[adj.motif.number], 1),
                  #                   motif_number==5 ~ paste0(LETTERS[adj.motif.number], 2),
                  #                   .default = LETTERS[adj.motif.number] )) %>%
                  label = LETTERS[motif_number]) %>%
    dplyr::select(-adj.motif.number)  %>%
    dplyr::rename(start_nt = start, end_nt = end, width_nt = width)
}


