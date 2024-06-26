library(data.table)
library(tidyverse)

source('hrd_score_fxns.R')


# Prep data ####

# Format Sanger copy number segment call data (downloaded from https://cellmodelpassports.sanger.ac.uk/downloads on 10/12/2022)
sanger_cn_segs <- readr::read_csv('~/data/sanger/WES_pureCN_CNV_segments_20220623.csv') %>%
  dplyr::filter(!(chr_name %in% c('chrX', 'chrY'))) %>%
  dplyr::filter(!is.na(chr_name)) %>%
  dplyr::mutate(minor_copy_number = tidyr::replace_na(minor_copy_number, 0)) %>%
  dplyr::transmute(model = paste(model_id, source, sep = '_'),
                   chr = as.numeric(gsub('chr', '', chr_name)),
                   start = chr_start,
                   end = chr_end,
                   n_snps = num_snps,
                   total_cn = round(total_copy_number),
                   allele_major = total_cn - minor_copy_number,
                   allele_minor = minor_copy_number,
                   arm)

# Downloaded BED file of centromere annotations from http://genome.ucsc.edu/cgi-bin/hgTables
centromeres <- readr::read_tsv('~/data/genomes/hg38/hg38_centromeres.bed', col_names = c('chr', 'cent_start', 'cent_end', 'X4')) %>%
  dplyr::group_by(chr) %>%
  dplyr::summarise(cent_start = min(cent_start),
                   cent_end = max(cent_end)) %>%
  as.data.frame() %>%
  dplyr::rowwise() %>%
  dplyr::mutate(cent_pos = mean(cent_start, cent_end)) %>%
  as.data.frame()

# Downloaded chromInfo.txt.gz from https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/
chrom_info <- readr::read_tsv('~/data/genomes/hg38/chromInfo.txt.gz', col_names = c('chr', 'chr_length', 'build')) %>%
  dplyr::inner_join(centromeres, by = 'chr') %>%
  dplyr::transmute(chr = gsub('chr', '', chr), chr_length, cent_pos, cent_start, cent_end) %>%
  dplyr::filter(chr %in% c(1:22)) %>%
  dplyr::mutate(chr = as.numeric(chr)) %>%
  dplyr::arrange(chr)


# Calculate HRD scores ####

## Detect LSTs
lsts <- sanger_cn_segs %>%
  annotate_arms(chr_info = chrom_info) %>%
  split_cent_segments(chr_info = chrom_info) %>%
  dplyr::group_by(model, chr, arm) %>%
  dplyr::mutate(group_id = dplyr::cur_group_id()) %>%
  split(.$group_id) %>%
  purrr::map_dfr(count_lsts)

## Detect tAIs
ntai <- sanger_cn_segs %>%
  # Remove segments <1Mb (as done in scarHRD)
  discard_short_segments(1e6) %>%
  # Remove segments with not enough SNPs (assuming ~50k SNPs total)
  dplyr::filter(n_snps >= 25) %>%
  dplyr::group_by(model, chr) %>%
  dplyr::mutate(group_id = dplyr::cur_group_id()) %>%
  split(.$group_id) %>%
  purrr::map_dfr(~ ntai_workflow(.x, chr_info = chrom_info))

## Detect LOH
lohs <- sanger_cn_segs %>%
  dplyr::group_by(model, chr) %>%
  dplyr::mutate(group_id = dplyr::cur_group_id()) %>%
  split(.$group_id) %>%
  purrr::map_dfr(loh_workflow)

## Sum up events
lst_counts <- lsts %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(n_lst = sum(lst)) %>%
  as.data.frame()

tai_counts <- ntai %>%
  dplyr::mutate(tai = ifelse(ai_type == 'telomeric', 1, 0)) %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(n_tai = sum(tai)) %>%
  as.data.frame()

loh_counts <- lohs %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(n_loh = sum(loh)) %>%
  as.data.frame()

all_counts <- loh_counts %>%
  dplyr::inner_join(lst_counts, by = 'model') %>%
  dplyr::inner_join(tai_counts, by = 'model') %>%
  dplyr::mutate(hrd_score = n_loh + n_lst + n_tai)

## Add model info (downloaded from https://cellmodelpassports.sanger.ac.uk/downloads)
model_info <- readr::read_csv('~/data/sanger/model_list_20220628.csv') %>%
  dplyr::transmute(model_id, 
                   model_name = toupper(gsub('[^[:alnum:]]', '', model_name)), 
                   tissue = gsub(' ', '_', tissue), 
                   cancer_type = gsub(' ', '_', cancer_type)) %>%
  # Differentiate duplicate model names to avoid downstream issues (will convert to CCLE names)
  dplyr::mutate(model_name = dplyr::case_when(
    model_id == 'SIDM00619' ~ 'KMHDASH2',
    model_id == 'SIDM00322' ~ 'TDOTT',
    TRUE ~ model_name
  ))

hrd_scores <- all_counts %>%
  tidyr::separate(model, c('model_id', 'wes_source'), sep = '_') %>%
  dplyr::left_join(model_info, by = 'model_id') %>%
  dplyr::select(model_id, model_name, tissue, cancer_type, wes_source, n_loh, n_lst, n_tai, hrd_score)
  

readr::write_csv(hrd_scores, 'sanger_hrd_scores.csv.gz')

