# Usage information

Download and update this pipeline with
```
nextflow pull ikmb/EASI_FullLength16S_pipleline
```
Run the pipeline as follow:
```
nextflow run main.nf --reads "/path/to/hifireads.hifi.bam" --ccga "/path/to/ccga_file.tsv" --movieid "ID123456" --outdir "results"
```
Mandatory:
'--reads': Location to the hifi.bam file
'--ccga': Location of ccga_file.tsv
'--movieid': Set a run ID

Optional:
'--outdir': Set a custom results directory. Default is "results"
'-work-dir': Set a custom working directy. Default is "work"
'--KEEP_DEMUL': Will copy demultiplexed files into the results directory.
'-resume': Will resume an earlier run.

Advanced parameters:
'--scratch': Set if scratch should be used
'--allbarcodes': Use a different barcodes files
Set with '-profile' a profile. Default is medcluster, add new profiles or run locally with '-profile local'
