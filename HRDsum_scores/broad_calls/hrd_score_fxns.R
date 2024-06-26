# Annotate chromosome arms
annotate_arms <- function(cn_df, chr_info) {
  cn_df %>%
    dplyr::left_join(chr_info, by = 'chr') %>%
    dplyr::mutate(arm = dplyr::case_when(
      start <= cent_start & end <= cent_end ~ 'p',
      start >= cent_start & end >= cent_end ~ 'q',
      start <= cent_start & end >= cent_end ~ 'p_and_q')) %>%
    dplyr::filter(arm %in% c('p', 'q', 'p_and_q')) %>%
    # Split segments that span the centromere into p and q portions
    dplyr::mutate(rep = ifelse(arm == 'p_and_q', 2, 1)) %>%
    tidyr::uncount(rep, .id = 'order') %>%
    dplyr::mutate(start = ifelse(arm == 'p_and_q' & order == 2, cent_end, start),
                  end = ifelse(arm == 'p_and_q' & order == 1, cent_start, end)) %>%
    dplyr::mutate(arm = dplyr::case_when(
      arm != 'p_and_q' ~ arm,
      order == 1 ~ 'p',
      order == 2 ~ 'q')) %>%
    dplyr::select(!c(order, chr_length, dplyr::starts_with('cent')))
}

# Remove segments < 3Mb
discard_short_segments <- function(arm_df) {
  arm_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(length = end - start) %>%
    as.data.frame() %>%
    dplyr::filter(length >= 3e6) %>%
    dplyr::select(!length)
}

# Collapse neighboring segments if they have the same allele CNs
collapse_segments <- function(arm_df) {
  arm_df %>%
    dplyr::mutate(allele_pair = paste0(allele_major, ',', allele_minor)) %>%
    dplyr::mutate(segment_id = data.table::rleid(allele_pair)) %>%
    dplyr::group_by(segment_id) %>%
    dplyr::mutate(start = min(start),
                  end = max(end),
                  n_probes = sum(n_probes)) %>%
    as.data.frame() %>%
    dplyr::distinct(start, end, n_probes, .keep_all = TRUE) %>%
    dplyr::select(!c('allele_pair', 'segment_id'))
}


# LST counting functions #######################################################

# Detect CN transitions where adjacent segments are both >10Mb and <3Mb apart from each other
detect_lsts <- function(arm_df) {
  arm_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(length = end - start) %>%
    as.data.frame() %>%
    dplyr::mutate(distance = start - dplyr::lag(end, default = dplyr::first(start))) %>%
    dplyr::mutate(lst = ifelse(
      length >= 10e6 & dplyr::lag(length, default = 0) >= 10e6 & distance < 3e6, 1, 0))
}

# LST detection workflow
count_lsts <- function(arm_group_df) {
  arm_group_df %>%
    # Discard segments with <50 probes/SNPs, then collapse neighboring segments if they have the same allele CNs
    dplyr::filter(n_probes >= 50) %>%
    purrr::when(nrow(.) > 0 ~ collapse_segments(.), ~ .) %>%
    # Discard segments <3Mb, then collapse neighboring segments if they have the same allele CNs
    discard_short_segments() %>%
    purrr::when(nrow(.) > 0 ~ collapse_segments(.), ~ .) %>%
    # Detect CN transitions where adjacent segments are both >10Mb and <3Mb apart from each other
    detect_lsts()
}


# tAI counting functions #######################################################

# Reassign ploidy (the major copy number state) at the chromosome level
ploidy_by_chromosome <- function(chr_df) {
  total_cn_stats <- chr_df %>%
    # Calculate segment length
    dplyr::rowwise() %>%
    dplyr::mutate(length = end - start) %>%
    as.data.frame() %>%
    # For each unique total_cn, sum up the total segment length
    dplyr::group_by(total_cn) %>%
    dplyr::mutate(length_by_ploidy = sum(length)) %>%
    as.data.frame()
  
  # Out of all total_cn states != 0, get one with longest total segment length
  chr_ploidy <- total_cn_stats %>%
    distinct(total_cn, length_by_ploidy) %>%
    dplyr::filter(total_cn != 0) %>%
    dplyr::slice_max(length_by_ploidy, n = 1) %>%
    dplyr::pull(total_cn)
  
  # Update ploidy column to reflect chromosome-level ploidy
  total_cn_stats %>%
    dplyr::mutate(ploidy = chr_ploidy) %>%
    dplyr::select(!c(length, length_by_ploidy))
}

# Categorize different types of allelic imbalance
categorize_ai <- function(chr_df, cent_info) {
  chr_df %>%
    # Count alleles as balanced if...
    dplyr::mutate(ai = dplyr::case_when(
      # ...ploidy is 1 or is even, and either allele_major == allele_minor or both allele CNs add up to ploidy and
      # the minor allele is always lower than ploidy/2 across the chromosome
      ploidy %in% c(1,seq(2, 200, by = 2)) & 
        (allele_major == allele_minor | (allele_major + allele_minor == ploidy & all(allele_minor < ploidy/2))) ~ FALSE,
      # ...ploidy is odd (and not 1) and both allele CNs add up to ploidy. If the minor allele = 0 and 
      # major allele = ploidy, only count as balanced if the minor allele loss was chromosome-wide
      !(ploidy %in% c(1,seq(2, 200, by = 2))) & 
        ((allele_major + allele_minor == ploidy & allele_minor != 0) | (allele_major == ploidy & all(allele_minor == 0))) ~ FALSE,
      TRUE ~ TRUE)) %>%
    dplyr::left_join(cent_info, by = 'chr') %>%
    dplyr::mutate(ai_type = dplyr::case_when(
      ai == FALSE ~ 'none',
      n() == 1 ~ 'whole_chromosome',
      !(row_number() %in% c(1, n())) ~ 'interstitial',
      row_number() == 1 & end < cent_pos ~ 'telomeric',
      row_number() == n() & start > cent_pos ~ 'telomeric',
      TRUE ~ '>1_arm'
    )) %>%
    dplyr::select(!c(chr_length, dplyr::starts_with('cent'), ai))
}

# tAI detection workflow
ntai_workflow <- function(chr_df, chr_info) {
  chr_df %>%
    collapse_segments() %>%
    ploidy_by_chromosome() %>%
    categorize_ai(cent_info = chr_info)
}


# LOH counting functions #######################################################

# Count LOH events
count_lohs <- function(chr_df) {
  chr_df %>%
    # Note segments where allele_major > 0 and allele_minor == 0
    dplyr::mutate(loss = ifelse(allele_major > 0 & allele_minor == 0, TRUE, FALSE)) %>%
    # Note instances where allele_minor == 0 across the entire chromosome
    dplyr::mutate(whole_chr_loh = ifelse(all(allele_minor == 0), TRUE, FALSE)) %>%
    # Calculate segment length
    dplyr::rowwise() %>%
    dplyr::mutate(length = end - start) %>%
    as.data.frame() %>%
    # Count LOH events that are > 15MB but less than the whole chromosome
    dplyr::mutate(loh = ifelse(length > 15e6 & loss == TRUE & whole_chr_loh == FALSE, 1, 0)) %>%
    dplyr::select(!c(loss, whole_chr_loh, length))
}

# Marquard et al approach PLUS filtering/smoothing
loh_workflow <- function(chr_df) {
  chr_df %>%
    # Discard segments with < 50 probes/SNPs
    dplyr::filter(n_probes >= 50) %>%
    # If allele_major > 1, reassign CN to 1 so that segments where allele_minor == 0 and allele_major > 0 can be collapsed
    dplyr::mutate(allele_major = ifelse(allele_major > 1, 1, allele_major)) %>%
    # Collapse neighboring segments if they have the same allele CNs
    collapse_segments() %>%
    purrr::when(nrow(.) > 0 ~ count_lohs(.), ~ .)
}


