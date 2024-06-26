# HRDsum_scores
Code for generating HRDsum scores in cell lines

The `broad_calls` directory contains the code used to calculate HRDsum scores for the "Broad" dataset.

The `sanger_calls` directory contains the code used to calculate HRDsum scores for the “Sanger (Broad WES)” and “Sanger (Sanger WES)” datasets.

The `summary_scores` directory contains the code used to merge the three datasets into one summary dataset, and includes the code used to harmonize tissue type annotations between the Broad and Sanger Institutes.


Scripts were executed with R version 4.1.3. Package versions include:
* tidyverse v1.3.1
* data.table v1.14.2
* ggpubr v0.4.0
* ggbeeswarm v0.6.0