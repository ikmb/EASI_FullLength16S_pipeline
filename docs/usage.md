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
'--hifireads': Location to the hifi.bam file  
'--ccga': Location of ccga_file.tsv  
'--movieid': Set a run ID  

### Optional:
'--outdir': Set a custom results directory. Default is "results".  
'-work-dir': Set a custom working directy. Default is "work".  
'--keep-demul': Will copy demultiplexed files into the results directory.  
'-resume': Will resume an earlier run.  

### Advanced:
'--scratch': Set if scratch should be used  
'--allbarcodes': Use a different barcodes files  
Set with '-profile' a profile. Default is medcluster, add new profiles or run locally with '-profile local'  
