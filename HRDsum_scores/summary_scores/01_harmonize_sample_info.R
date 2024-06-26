library(tidyverse)

# CCLE cell line info (downloaded from https://depmap.org) ####
ccle_info <- readr::read_csv('~/data/depmap/22Q4/Model.csv') %>%
  dplyr::transmute(cell_line = StrippedCellLineName, tissue = OncotreeLineage, disease = OncotreePrimaryDisease,
                   DepMap_ID = ModelID, Sanger_ID = SangerModelID)

# Sanger cell line info (downloaded from https://cellmodelpassports.sanger.ac.uk/downloads) ####
sanger_info <- readr::read_csv('~/data/sanger/model_list_20230110.csv') %>%
  dplyr::transmute(cell_line = toupper(gsub('[^[:alnum:]]', '', model_name)), 
                   tissue = gsub(' ', '_', tissue),
                   Sanger_ID = model_id,
                   DepMap_ID = BROAD_ID,
                   cancer_type) %>%
  # Differentiate duplicate model names to avoid downstream issues (will convert to CCLE names)
  dplyr::mutate(cell_line = dplyr::case_when(
    Sanger_ID == 'SIDM00408' ~ 'MSDASH1',
    Sanger_ID == 'SIDM00619' ~ 'KMHDASH2',
    Sanger_ID == 'SIDM00322' ~ 'TDOTT',
    TRUE ~ cell_line
  )) %>%
  # Remove duplicate cell line names
  dplyr::filter(!(Sanger_ID %in% c('SIDM00440', 'SIDM01656', 'SIDM01919')))


# Check model matches between Sanger and CCLE ####

# Which matches in ccle_info are absent from or don't agree with sanger_info?
ccle_match_only <- ccle_info %>%
  dplyr::rename(DepMap_ID_CCLE = DepMap_ID) %>%
  dplyr::inner_join(sanger_info, by = 'Sanger_ID') %>%
  dplyr::filter(DepMap_ID_CCLE != DepMap_ID | is.na(DepMap_ID))

# Which Sanger models match to more than one DepMap ID?
multiple_depmap_ids <- sanger_info %>%
  dplyr::filter(grepl(';', DepMap_ID))

# Are there any matches in cell line names that weren't captured by ID matches?
missed_matches <- ccle_info %>%
  dplyr::inner_join(sanger_info, by = 'cell_line', suffix = c('.ccle', '.sanger')) %>%
  dplyr::filter(is.na(DepMap_ID.sanger) & is.na(Sanger_ID.ccle))

# Fill in/correct matches in sanger_info
sanger_info$DepMap_ID[which(sanger_info$Sanger_ID == 'SIDM00096')] <- 'ACH-000338'
sanger_info$DepMap_ID[which(sanger_info$Sanger_ID == 'SIDM00969')] <- 'ACH-001063'
sanger_info$DepMap_ID[which(sanger_info$Sanger_ID == 'SIDM01095')] <- 'ACH-001189'
sanger_info$DepMap_ID[which(sanger_info$Sanger_ID == 'SIDM00399')] <- 'ACH-001227'
sanger_info$DepMap_ID[which(sanger_info$Sanger_ID == 'SIDM00603')] <- 'ACH-001543'
sanger_info$DepMap_ID[which(sanger_info$Sanger_ID == 'SIDM00069')] <- 'ACH-002024'
sanger_info$DepMap_ID[which(sanger_info$Sanger_ID == 'SIDM00246')] <- 'ACH-002084'
sanger_info$DepMap_ID[which(sanger_info$Sanger_ID == 'SIDM00104')] <- 'ACH-002123'
sanger_info$DepMap_ID[which(sanger_info$Sanger_ID == 'SIDM00518')] <- 'ACH-002136'
sanger_info$DepMap_ID[which(sanger_info$Sanger_ID == 'SIDM01219')] <- 'ACH-002158'
sanger_info$DepMap_ID[which(sanger_info$Sanger_ID == 'SIDM00842')] <- 'ACH-002219'
sanger_info$DepMap_ID[which(sanger_info$Sanger_ID == 'SIDM00517')] <- 'ACH-002239'
sanger_info$DepMap_ID[which(sanger_info$Sanger_ID == 'SIDM00619')] <- 'ACH-002397'
sanger_info$DepMap_ID[which(sanger_info$Sanger_ID == 'SIDM00400')] <- 'ACH-002392'
sanger_info$DepMap_ID[which(sanger_info$Sanger_ID == 'SIDM01748')] <- 'ACH-001712'
sanger_info$DepMap_ID[which(sanger_info$Sanger_ID == 'SIDM01845')] <- 'ACH-002014'

# Are there any matches where the cell line names are suspiciously different?
suspicious_matches <- ccle_info %>%
  dplyr::select(cell_line, tissue, DepMap_ID) %>%
  dplyr::full_join(sanger_info, by = 'DepMap_ID', suffix = c('.ccle', '.sanger')) %>%
  dplyr::filter(cell_line.ccle != cell_line.sanger)

sanger_info$DepMap_ID[which(sanger_info$Sanger_ID == 'SIDM01387')] <- NA
sanger_info$DepMap_ID[which(sanger_info$Sanger_ID == 'SIDM00911')] <- NA
sanger_info$DepMap_ID[which(sanger_info$Sanger_ID == 'SIDM01558')] <- NA


# Merge info tables ####
all_info <- ccle_info %>%
  dplyr::transmute(cell_line, tissue, disease.ccle = disease, DepMap_ID) %>%
  dplyr::full_join(sanger_info, by = 'DepMap_ID', suffix = c('.ccle', '.sanger')) %>%
  # Use CCLE cell line name unless it's missing
  dplyr::mutate(cell_line = ifelse(is.na(cell_line.ccle), cell_line.sanger, cell_line.ccle)) %>%
  dplyr::select(cell_line, DepMap_ID, Sanger_ID, tissue.ccle, tissue.sanger, disease.ccle, cancer_type)

# Harmonize tissue annotations ####
all_info_tissue <- all_info %>%
  # Recode non-cancerous lines as non-cancerous
  dplyr::mutate(tissue.ccle = ifelse(disease.ccle == 'Non-Cancerous', 'Non-Cancerous', tissue.ccle),
                tissue.sanger = ifelse(cancer_type == 'Non-Cancerous', 'Non-Cancerous', tissue.sanger)) %>%
  # Recode Sanger annotations to match CCLE
  dplyr::mutate(tissue.sanger = dplyr::recode(
    tissue.sanger, 
    Adrenal_Gland = 'Adrenal Gland',
    Bladder = 'Bladder/Urinary Tract',
    Biliary_Tract = 'Biliary Tract',
    Large_Intestine = 'Bowel',
    Small_Intestine = 'Bowel',
    Central_Nervous_System = 'CNS/Brain',
    Esophagus = 'Esophagus/Stomach',
    Stomach = 'Esophagus/Stomach',
    Head_and_Neck = 'Head and Neck',
    Ovary = 'Ovary/Fallopian Tube',
    Peripheral_Nervous_System = 'Peripheral Nervous System',
    Soft_Tissue = 'Soft Tissue',
    Endometrium = 'Uterus',
    Vulva = 'Vulva/Vagina')) %>%
  dplyr::mutate(tissue = dplyr::case_when(
    # When tissue calls agree, default to CCLE calls
    tissue.ccle == tissue.sanger ~ tissue.ccle,
    # When tissue calls don't agree, default to CCLE calls
    tissue.ccle != tissue.sanger ~ tissue.ccle,
    # When tissue.sanger is missing, use tissue.ccle
    is.na(tissue.sanger) ~ tissue.ccle,
    # If tissue.ccle is missing and tissue.sanger == 'Haematopoietic_and_Lymphoid', 
    # use cancer type to help categorize line as 'Lymphoid' or 'Myeloid'
    is.na(tissue.ccle) & tissue.sanger == 'Haematopoietic_and_Lymphoid' & grepl('Lympho', cancer_type) ~ 'Lymphoid',
    is.na(tissue.ccle) & tissue.sanger == 'Haematopoietic_and_Lymphoid' & grepl('Myelo', cancer_type) ~ 'Lymphoid'
  ))

# Of the annotations available from Sanger and not available from CCLE, are there any that don't match one of the CCLE annotations?
missing_ccle <- all_info_tissue %>%
  dplyr::filter(is.na(tissue))

unique(missing_ccle$tissue.sanger)[which(!(unique(missing_ccle$tissue.sanger) %in% unique(all_info$tissue.ccle)))]

# Categorize unresolved cell lines using help from Cellosaurus/Sanger's cancer type annotations
lymphoid <- c('MAC2B', 'NK92')
other <- c('MY', 'HS630T', 'BEWO')


# Finalize annotations
all_info_final <- all_info_tissue %>%
  dplyr::mutate(tissue = dplyr::case_when(
    cell_line %in% lymphoid ~ 'Lymphoid',
    cell_line %in% other ~ 'Other',
    is.na(tissue) & cancer_type == 'Mesothelioma' ~ 'Pleura',
    is.na(tissue) ~ tissue.sanger,
    !is.na(tissue) ~ tissue
  )) %>%
  dplyr::select(cell_line, DepMap_ID, Sanger_ID, tissue)


# Save
readr::write_csv(all_info_final, 'ccle_sanger_model_info.csv')

