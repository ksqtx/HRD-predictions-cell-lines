library(tidyverse)
library(CHORD)
library(BSgenome.Hsapiens.UCSC.hg38)

args <- commandArgs(trailingOnly = TRUE)
cell_line <- args[1]
source <- args[2]
snv_indel_vcf <- args[3]
sv_file <- args[4]

sample <- paste(cell_line, source, sep = '_')

sv_df <- readr::read_tsv(sv_file, 
                         col_types = list(sv_type = 'c', sv_len = 'd'))

# Extract mutation contexts
CHORD::extractSigsChord(
  vcf.snv = snv_indel_vcf,
  df.sv = sv_df,
  vcf.filters=list(snv = 'PASS', indel = 'PASS'),
  sample.name = sample,
  output.path = 'chord_mut_counts.txt',
  ref.genome=BSgenome.Hsapiens.UCSC.hg38,
  verbose = T
)

contexts <- read.delim('chord_mut_counts.txt', check.names = FALSE)

# Predict HRD
chord_res <- CHORD::chordPredict(contexts, do.bootstrap = T, verbose = T, min.msi.indel.rep = 10000)

# Save results
contexts %>%
  tibble::rownames_to_column('sample') %>%
  readr::write_tsv(paste(sample, 'contexts.txt', sep = '_'))

chord_res %>%
  readr::write_tsv(paste(sample, 'chord.txt', sep = '_'))

