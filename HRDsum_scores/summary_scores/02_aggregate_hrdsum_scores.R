library(tidyverse)
library(splines)
library(ggpubr)
library(ggbeeswarm)


# Cell line info
sample_info <- readr::read_csv('ccle_sanger_model_info.csv')

# Get all three sources of HRD scores (generated in broad_calls and sanger_calls directories)

## Sanger line HRD scores (calculated using Sanger-generated copy number calls)
hrd_scores_sanger_sanger <- readr::read_csv('sanger_hrd_scores.csv.gz') %>%
  dplyr::select(model_id, wes_source, hrd_score) %>%
  dplyr::left_join(sample_info, by = c('model_id' = 'Sanger_ID')) %>%
  dplyr::filter(wes_source == 'Sanger') %>%
  dplyr::transmute(Sanger_ID = model_id, hrd_score_sanger_sanger = hrd_score)

## Broad line HRD scores (calculated using Sanger-generated copy number calls)
hrd_scores_sanger_broad <- readr::read_csv('sanger_hrd_scores.csv.gz') %>%
  dplyr::select(model_id, wes_source, hrd_score) %>%
  dplyr::left_join(sample_info, by = c('model_id' = 'Sanger_ID')) %>%
  dplyr::filter(wes_source == 'Broad') %>%
  dplyr::transmute(Sanger_ID = model_id, hrd_score_sanger_broad = hrd_score)

## Broad line HRD scores (calculated using CCLE-generated copy number calls)
hrd_scores_broad <- readr::read_csv('broad_hrd_scores.csv.gz') %>%
  dplyr::transmute(cell_line, hrd_score_broad = hrd_score) %>%
  dplyr::mutate(cell_line = dplyr::recode(cell_line, D341 = 'D341Med'))

## Join datasets
hrd_scores_all <- sample_info %>%
  dplyr::left_join(hrd_scores_broad, by = 'cell_line') %>%
  dplyr::left_join(hrd_scores_sanger_sanger, by = 'Sanger_ID') %>%
  dplyr::full_join(hrd_scores_sanger_broad, by = 'Sanger_ID') %>%
  dplyr::filter(!dplyr::if_all(where(is.numeric), ~ is.na(.x)))



# Integrate datasets using natural spline regression ####

dat <- hrd_scores_all

dat1 <- dat %>% 
  dplyr::select(cell_line, hrd_score_broad, hrd_score_sanger_sanger, hrd_score_sanger_broad) %>% 
  data.frame(row.names = 1)

ggplot(stack(dat1), aes(x=ind, y=values)) + theme_pubr() + 
  geom_violin() + geom_beeswarm() + ggtitle('raw all')

ggplot(stack(dat1), aes(group=ind, sample=values))  + theme_pubr() + stat_qq(aes(col=ind)) + 
  stat_qq_line(aes(col=ind)) + ggtitle('raw all') + xlab('Theoretical quantile') + ylab('Observed quantile')


comm_cell_line <- rownames(dat1)[which(rowSums(!is.na(dat1))==3)]
comm_line_dat <- dat1[comm_cell_line,]

ggplot(stack(dat1[comm_cell_line,]), aes(x=ind, y=values))  + theme_pubr() + 
  geom_violin()+ geom_beeswarm() + ggtitle('raw common')


comm_med_hrd <- apply(comm_line_dat,1,median)

plotdf <- data.frame(median=comm_med_hrd, dat1[comm_cell_line,])

plotdf %>% tidyr::pivot_longer(cols=2:4) %>% 
  ggplot(aes(median, value)) + geom_point(aes(col=name)) + 
  geom_smooth(aes(group=name, col=name), method='loess')+
  theme_pubr()


dat2 <- dat1
for(i in 1:3){
  id <- which(!is.na(dat1[,i]))
  nsmod <- lm(comm_med_hrd ~ ns(comm_line_dat[,i],3))
  yout <- round(coef(nsmod)[-1]%*%t(ns(dat1[id,i],3)))
  dat2[id,i] <- t(yout)
}


ggplot(stack(dat2[comm_cell_line,]), aes(x=ind, y=values))  + theme_pubr() + 
  geom_violin()+ geom_beeswarm() + ggtitle('normalized common')


plotdf <- data.frame(median=comm_med_hrd, dat2[comm_cell_line,])

plotdf %>% tidyr::pivot_longer(cols=2:4) %>% 
  ggplot(aes(median, value)) + geom_point(aes(col=name)) + 
  geom_smooth(aes(group=name, col=name), method='loess')+
  theme_pubr() + ggtitle('Normalized common')

ggplot(stack(dat2), aes(x=ind, y=values))  + theme_pubr() + 
  geom_violin()+ geom_beeswarm() + ggtitle('normalized all')


ggplot(stack(dat2), aes(group=ind, sample=values))  + theme_pubr() + stat_qq(aes(col=ind)) + 
  stat_qq_line(aes(col=ind)) + ggtitle('Normalized all') + xlab('Theoretical quantile') + ylab('Observed quantile')



# Calculate mean normalized score and add to table of raw HRD scores
hrd_scores_final <- dat2 %>%
  tibble::rownames_to_column('cell_line') %>%
  dplyr::rowwise() %>%
  dplyr::mutate(hrd_score_summary = round(mean(c(hrd_score_broad, hrd_score_sanger_sanger, hrd_score_sanger_broad), na.rm = TRUE))) %>%
  dplyr::ungroup() %>%
  dplyr::select(cell_line, hrd_score_summary) %>%
  dplyr::inner_join(hrd_scores_all, by = 'cell_line') %>%
  dplyr::transmute(cell_line, DepMap_ID, Sanger_ID, tissue, 
                   hrdsum_broad = hrd_score_broad, 
                   hrdsum_sanger_sangerwes = hrd_score_sanger_sanger, 
                   hrdsum_sanger_broadwes = hrd_score_sanger_broad, 
                   hrdsum_summary = hrd_score_summary)


readr::write_csv(hrd_scores_final, 'summary_hrd_scores.csv')



