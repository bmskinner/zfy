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

# aa - the amino acid sequence. May have gaps if part of a multiple sequence alignment
# sequence.name - the name of the sequence
# rc.threshold - threshold for percent of refinement criteria matching
# 
# Returns - a tibble with 9aaTADs, their coordinates in the sequence, the
# sequence name, and the rc score (as a %). If the input contains gaps ("-"),
# the coordinates will be with respect to the gapped alignment.
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
    rc1 <- function(s) sum(s[c(2, 4, 5)] %in% s2c("ILVFWMAY"))>=1
    rc2 <- function(s) sum(s[c(6, 8)] %in% s2c("ILVFWMAY"))>=1
    rc3 <- function(s) sum(s[c(4, 5, 6)] %in% s2c("DENQST"))>=1
    rc4 <- function(s) {count <- sum(s[1:9] %in% s2c("DENQSTKRH")); count>=2 & count<=4}
    rc5 <- function(s) {count <- sum(s[1:9] %in% s2c("ILVFWMAYGPST")); count>=5 & count<=7}
    rc6 <- function(s) { 
      m <- s[1:9] %in% s2c("ILVFWMAYGP");
      !any( sapply(3:8, function(i) { all(m[i-2], m[i-1], m[i], m[i+1]) })) &
        !any( sapply(2:7, function(i) { all(m[i-1], m[i], m[i+1], m[i+2]) }))
    }
    rc7 <- function(s) { 
      m <- s[1:9] %in% s2c("DENQST");
      !any(sapply(2:8, function(i) all(m[i-1], m[i], m[i+1])))
    }
    rc8  <- function(s) sum(s[1:9] %in% s2c("KRH"))<=1
    rc9  <- function(s) sum(s[1:6] %in% s2c("QNKRH"))<=1
    rc10 <- function(s) sum(s[8:9] %in% s2c("QNKRH"))<=1
    rc11 <- function(s) sum(s[2:8] %in% s2c("C"))>=0
    rc12 <- function(s) sum(s[2:7] %in% s2c("PG"))>=0
    
    
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
