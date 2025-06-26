# mtRtN - mtDNA Variant Calling with RtN

This Nextflow pipeline will call SNVs and Indels for the mitochondrial genome. It will remove potentially confounding reads introduced by germline integrations of the mitochondrial genome into the nuclear genome using the Remove the NUMTs (RtN) algorithm.

## NUMT Removal

NUMTs (pronounced new-mites) are integrations of the mitochondrial genome into the nuclear DNA. This process has occurred throughout mammalian evolution, therefore mtDNA fragments are common throughout the genome. Due to NUMTs similarity to "authentic" circular mtDNA, sequencing reads originating from NUMTs can be wrongly mapped back to the mitochondrial genome. This is an issue as these additional reads can introduce false positives into the variant calling pipeline. 

To remove these potentially confounding reads, mtRtN uses the Remove the NUMTs (RtN) algorithm to align mitochondrial reads to known NUMT sequences. A goodness of fit threshold is then set, removing reads from the alignment file that may have originated from NUMTs. See the original [manuscript](https://academic.oup.com/bioinformatics/article/36/20/5115/5876023) and [github repository](https://github.com/Ahhgust/RtN) for more details. 

## Running the Pipeline

### 1. Cloning Repository

Clone the repository and cd into the directory using the following commands:

```shell
git clone https://git.ecdf.ed.ac.uk/s2078878/mtRtN2.git
cd mtRtN2
```

### 2. Uncompress Reference Files

Uncompress the NUMT reference files for the pipeline to use using the following command:

```
tar -xvjf rtn.tar.bz2
```

### 3. Make Input Files

#### Matched-Samples
The Nextflow pipeline works from a .csv file containing the main parameters needed. The csv file should have the following format and column names:

| row | sampleID | group  | bam                                  | index                                |
|:----|:---------|:-------|:-------------------------------------|:-------------------------------------|
| 1   | sample1  | normal | path/to/sample1/normal/alignment.bam | path/to/sample1/normal/index.bam.bai |
| 2   | sample1  | tumour | path/to/sample1/tumour/alignment.bam | path/to/sample1/tumour/index.bam.bai |
| 3   | sample2  | normal | path/to/sample2/normal/alignment.bam | path/to/sample2/normal/index.bam.bai |
| 4   | sample2  | tumour | path/to/sample2/tumour/alignment.bam | path/to/sample2/tumour/index.bam.bai |

Alignment and index files can be for the entire genome, there is no need to subset to only the mitochondrial DNA.

In addition to this a sample purity file can also be supplied. This should be in the below format and should be a .csv file. sampleIDs should match those in the parameters file. A purity of 1 would indicate a completely pure sample. 

| row | sampleID | group  | purity |
|:----|:---------|:-------|:-------|
| 1   | sample1  | normal | 0.812  |
| 2   | sample1  | tumour | 0.516  |
| 3   | sample2  | normal | 0.753  |
| 4   | sample2  | tumour | 0.625  |

If purity data is not available, do not set the purity parameter when running mtRtN. 

#### Unmatched-Samples

The pipeline will also allow you to run samples with no matched sample. Similarly, a .csv parameters file needs to be provided containing the locations of the alignment and index files. Unlike the matched-samples, a group name does not need to be provided. An example of a unmatched parameters file will look something similar to this:

| row | sampleID | bam                           | index                         |
|:----|:---------|:------------------------------|:------------------------------|
| 1   | sample1  | path/to/sample1/alignment.bam | path/to/sample1/index.bam.bai |
| 2   | sample2  | path/to/sample2/alignment.bam | path/to/sample2/index.bam.bai |
| 3   | sample3  | path/to/sample3/alignment.bam | path/to/sample3/index.bam.bai |
| 4   | sample4  | path/to/sample4/alignment.bam | path/to/sample4/index.bam.bai |

However, VarScan2 will not accept purity data for unmatched samples, therefore even if purity data is provided, it will not be used in variant calling. 

### 4. Run mtRtN

Everything should now be ready to run mtRtN. A typical command to run the pipeline on the University of Edinburgh's Eddie HPC will look something similar to the one below, substituting in the path to your own parameter file. Make sure that Nextflow and Singularity are loaded beforehand. Singularity will be used to download and run the pipeline through a docker container. If desired a purity file can also be provided. Take additional note of setting the `--matched` parameter if samples do not have a corresponding matched normal.

```shell
nextflow run main.nf -profile eddie,singularity \
  --parameters path/to/parameters/file.csv \ # parameters file (see step 3)
  --purity /path/to/purity/file.csv \ # purity data, exlude if no purity data available or are unmatched samples
  --outdir path/to/out/directory/ \ # where should the results folder be created? default: run directory
  --matched (true / false) \ # are the samples matched? default: true
  --rnaseq (true / false) # is the data from RNA-sequencing? default: false
```

### 5. Outputs

All pipeline output files will be created in a results directory, however outputs will depend on whether matched or unmatched samples were used. 

#### Matched-Samples

For matched-samples, both somatic and germline variant calls will be generated. All files are the raw unfiltered variant calls.

#### Unmatched-Samples

For unmatched-samples, only one file will be generated. This will be the raw unfiltered variant calls. 

## Notes

Please be aware that this pipelines `nextflow.config` file has been configured to run on the University of Edinburgh's Eddie HPC Cluster. This may need to be reconfigured to run on other computing systems. 







