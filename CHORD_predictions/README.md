# CHORD_predictions
Pipeline for CHORD pre-processing and analysis

Before running the pipeline: 
* Build the Docker image with `docker build -t ksqtx/chord .`
* Run `export NXF_DEFAULT_DSL=1`
* Configure the file paths in `nextflow.config`

## Example run

*Using Nextflow version 22.10.6*

```bash
nextflow run main.nf \
--sampleInfo /path/to/samplesheet.csv \
--outDir /path/to/chord_predictions \
-work-dir /path/to/scratch
```

The CSV samplesheet should contain the following columns (with a header), in this order:

* `sample`: sample name
* `source`: sequencing data source (e.g. CCLE)
* `sage_vcf`: path to SAGE VCF
* `sage_index`: path to index for SAGE VCF
* `gripss_vcf`: path to filtered GRIPSS VCF
* `gripss_index`: path to index for filtered GRIPSS VCF