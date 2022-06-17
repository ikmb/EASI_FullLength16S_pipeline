# Usage information

Download and update this pipeline with
```
nextflow pull ikmb/EASI_FullLength16S_pipeline
```

Run the pipeline:
```
nextflow run ikmb/EASI_FullLength16S_pipeline --hifireads "/path/to/hifireads.hifi.bam" --ccga "/path/to/ccga_file.tsv" --movieid "ID123456" --outdir "results"
```
## Parameters

### Mandatory:
`--hifireads`: Location to the hifi.bam file  <br />
`--ccga`: Location of ccga_file.tsv  <br />
`--movieid`: Set a run ID  <br />

### Optional:
`--outdir`: Set a custom results directory. Default is "results".  <br />
`-work-dir`: Set a custom working directy. Default is "work".  <br />
`--keep-demul`: Will copy demultiplexed files into the results directory.  <br />
`-resume`: Will resume an earlier run.  

### Advanced:
`--scratch`: Set if scratch should be used  <br />
`--allbarcodes`: Use a different barcodes files  <br />
Set with `-profile` a profile. Default is medcluster, add new profiles or run locally with `-profile local`. Profile 'local' must be configured for your specific system!  <br />
