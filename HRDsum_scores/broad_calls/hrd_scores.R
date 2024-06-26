library(data.table)
library(tidyverse)

source('hrd_score_fxns.R')


# Prep data ####

# Format Broad copy number segment call data (downloaded from https://depmap.org/)
absolute_segments <- readxl::read_excel('~/data/depmap/other/CCLE_ABSOLUTE_combined_20181227.xlsx',
                                        sheet = 'ABSOLUTE_combined.segtab', na = 'NA') %>%
  dplyr::filter(sample != 'ALEXANDERCELLS_LIVER') %>%
  dplyr::transmute(depMapID,
                   chr = Chromosome,
                   start = Start,
                   end = End,
                   n_probes = Num_Probes,
                   total_cn = Modal_HSCN_1 + Modal_HSCN_2,
                   allele_major = Modal_HSCN_2,
                   allele_minor = Modal_HSCN_1)

absolute_summary <- readxl::read_excel('~/data/depmap/other/CCLE_ABSOLUTE_combined_20181227.xlsx',
                                       sheet = 'ABSOLUTE_combined.table') %>%
  dplyr::filter(CCLE_ID != 'ALEXANDERCELLS_LIVER') %>%
  dplyr::select(depMapID, ploidy, purity)

#' DepMap ID ACH-001318 was mapped to two different sample names (ALEXANDERCELLS and PLCPRF5).
#' Just kept the name that's currently assigned on depmap.org (PLCPRF5)

# Sample info (downloaded from https://depmap.org/)
sample_info <- readr::read_csv('~/data/depmap/21Q3/sample_info.csv.gz') %>%
  dplyr::transmute(depMapID = DepMap_ID,
                   cell_line = stripped_cell_line_name)

cn_profile <- absolute_segments %>%
  dplyr::left_join(absolute_summary, by = 'depMapID') %>%
  dplyr::left_join(sample_info, by = 'depMapID') %>%
  dplyr::select(!depMapID) %>%
  dplyr::relocate(cell_line)


# Assemble table with info on chromosome length and centromere positions
## Downloaded `chromInfo.txt.gz` and `gap.txt.gz` from https://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/
chr_length <- readr::read_tsv('~/data/genomes/hg19/chromInfo.txt.gz', 
                              col_names = c('chr', 'chr_length', 'build'))

centromeres <- readr::read_tsv('~/data/genomes/hg19/gap.txt.gz', 
                               col_names = c('X1', 'chr', 'start', 'end', 'X5', 'X6', 'X7', 'type', 'X9')) %>%
  dplyr::filter(type == 'centromere') %>%
  dplyr::rowwise() %>%
  dplyr::transmute(chr, cent_pos = mean(c(start, end)), cent_start = start, cent_end = end) %>%
  as.data.frame()

chrom_info <- centromeres %>%
  dplyr::left_join(chr_length, by = 'chr') %>%
  dplyr::transmute(chr = gsub('chr', '', chr), chr_length, cent_pos, cent_start, cent_end) %>%
  dplyr::filter(chr %in% c(1:22)) %>%
  dplyr::mutate(chr = as.numeric(chr)) %>%
  dplyr::arrange(chr)


# Calculate HRD scores ####

## Detect LSTs
lsts <- cn_profile %>%
  annotate_arms(chr_info = chrom_info) %>%
  dplyr::group_by(cell_line, chr, arm) %>%
  dplyr::mutate(group_id = dplyr::cur_group_id()) %>%
  split(.$group_id) %>%
  purrr::map_dfr(count_lsts)

## Detect tAIs
ntai <- cn_profile %>%
  # Remove segments with not enough SNPs, remove samples with too much contamination
  dplyr::filter(n_probes >= 500, purity > 0.36) %>%
  dplyr::group_by(cell_line, chr) %>%
  dplyr::mutate(group_id = dplyr::cur_group_id()) %>%
  split(.$group_id) %>%
  purrr::map_dfr(~ ntai_workflow(.x, chr_info = chrom_info))

## Detect LOH
lohs <- cn_profile %>%
  dplyr::group_by(cell_line, chr) %>%
  dplyr::mutate(group_id = dplyr::cur_group_id()) %>%
  split(.$group_id) %>%
  purrr::map_dfr(loh_workflow)

## Sum up events
lst_counts <- lsts %>%
  dplyr::group_by(cell_line) %>%
  dplyr::summarise(n_lst = sum(lst)) %>%
  as.data.frame()

tai_counts <- ntai %>%
  dplyr::mutate(tai = ifelse(ai_type == 'telomeric', 1, 0)) %>%
  dplyr::group_by(cell_line) %>%
  dplyr::summarise(n_tai = sum(tai)) %>%
  as.data.frame()

loh_counts <- lohs %>%
  dplyr::group_by(cell_line, ploidy) %>%
  dplyr::summarise(n_loh = sum(loh)) %>%
  as.data.frame()

all_counts <- loh_counts %>%
  dplyr::inner_join(lst_counts, by = 'cell_line') %>%
  dplyr::inner_join(tai_counts, by = 'cell_line') %>%
  dplyr::mutate(hrd_score = n_loh + n_lst + n_tai)

readr::write_csv(all_counts, 'broad_hrd_scores.csv.gz')

