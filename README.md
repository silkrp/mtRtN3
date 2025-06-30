# mtRtN3 - mtDNA Variant Calling with RtN

This Nextflow pipeline identifies mitochondrial SNVs and indels, estimates mtDNA copy number, and calculates sequencing depth. It removes potentially confounding reads arising from nuclear integrations of mitochondrial DNA (NUMTs) using the Remove the NUMTs (RtN) algorithm.

## Running the Pipeline

### 1. Cloning Repository

Clone the repository and cd into the directory using the following commands:

```shell
git clone https://github.com/silkrp/mtRtN3.git
cd mtRtN3
```

### 2. Uncompress Reference Files

Uncompress the NUMT reference files required by the pipeline using the following command:

```shell
tar -xvjf rtn.tar.bz2
```

### 3. Download VEP Cache

The mtRtN3 pipeline annotates mtDNA variant calls using Ensembl's Variant Effect Predictor (VEP). This requires a local VEP cache, which can be downloaded using the following commands:

```
curl -O https://ftp.ensembl.org/pub/release-114/variation/indexed_vep_cache/homo_sapiens_vep_114_GRCh38.tar.gz
tar xzf homo_sapiens_vep_114_GRCh38.tar.gz
```

### 4. Make Input Files

#### Matched-Samples
The Nextflow pipeline works from a .csv file containing the main parameters needed. The csv file should have the following format and column names:

| sample   | group  | bam                                  | index                                |
|:---------|:-------|:-------------------------------------|:-------------------------------------|
| sample1  | normal | path/to/sample1/normal/alignment.bam | path/to/sample1/normal/index.bam.bai |
| sample1  | tumour | path/to/sample1/tumour/alignment.bam | path/to/sample1/tumour/index.bam.bai |
| sample2  | normal | path/to/sample2/normal/alignment.bam | path/to/sample2/normal/index.bam.bai |
| sample2  | tumour | path/to/sample2/tumour/alignment.bam | path/to/sample2/tumour/index.bam.bai |

Alignment and index files can be for the entire genome, there is no need to subset to only the mitochondrial DNA.

#### Unmatched-Samples

The pipeline will also allow you to run samples with no matched sample. Similarly, a .csv parameters file needs to be provided containing the locations of the alignment and index files. 

| sample   | group  | bam                                  | index                                |
|:---------|:-------|:-------------------------------------|:-------------------------------------|
| sample1  | tumour | path/to/sample1/tumour/alignment.bam | path/to/sample1/tumour/index.bam.bai |
| sample2  | tumour | path/to/sample2/tumour/alignment.bam | path/to/sample2/tumour/index.bam.bai |
| sample3  | tumour | path/to/sample3/tumour/alignment.bam | path/to/sample3/tumour/index.bam.bai |
| sample4  | tumour | path/to/sample4/tumour/alignment.bam | path/to/sample4/tumour/index.bam.bai |

### 5. Run mtRtN3

Everything should now be ready to run mtRtN3. A typical command to run the pipeline on the University of Edinburgh's Eddie HPC will look something similar to the one below.

```shell
nextflow run main.nf -profile eddie,singularity \

  ## Main settings
  --input path/to/input/parameters/file.csv \         # input parameters file (see step 4)
  --outdir path/to/out/directory/ \                   # where should the results folder be created? default: run directory
  --workflows variant-calling,cn,depths \              # which workflows should be run? default: variant-calling,cn,depths
  --contig rCRS \                                     # in your bam file, how is the mitochondrial DNA contig denoted? default: rCRS 
  --fasta path/to/mtdna/reference.fa                  # path to mtDNA reference fasta file. Note that the mtDNA contig name in the reference file should match that used in your bam file. default: ./misc/rCRS.fa 
  --matched (true / false) \                          # are the samples matched? default: true
  --rnaseq (true / false)                             # is the data from RNA-sequencing? default: false

  ## VEP settings
  --vep_species = 'homo_sapiens'                      # vep species. default: homo_sapiens
  --vep_genome = 'GRCh38'                             # vep genome. default: GRCh38
  --vep_cache_version = 114                           # vep cache version. default: 114
  --vep_cache = "./"                                  # location of the vep cache. default: ./

  ## Other
  --cram_fasta path/to/cram/reference.fa              # path to reference file used to generate cram files. default: null
  --gnomad_vcf path/to/mtdna/gnomad.vcf               # path to gnomad mtDNA reference vcf. default: ./misc/gnomad.genomes.v3.1.sites.chrM.vcf.bgz
```

Make sure Nextflow and Singularity are loaded before running the pipeline. Singularity will handle containerized execution using Docker images.

### 6. Outputs

All pipeline output files will be created in a results directory, however outputs will depend on whether matched or unmatched samples were used. 

#### Matched-Samples

For matched samples, the pipeline generates both somatic and germline mtDNA variant calls (note that these are raw and unfiltered), along with mtDNA sequencing depth and copy number estimates for both tumour and normal samples. 

#### Unmatched-Samples

For unmatched samples, the pipeline generates raw, unfiltered mtDNA variant calls, as well as mtDNA sequencing depth and copy number estimates.

## Notes

This pipeline is configured to run on the University of Edinburgh’s Eddie HPC via the `nextflow.config` file. You may need to modify it to suit other systems.