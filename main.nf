#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/**
===============================
EASI_FullLength16S_pipleline
===============================

### Homepage / git
git@github.com:ikmb/EASI_FullLength16S_pipleline.git

**/

// Pipeline version

params.version = workflow.manifest.version

// Help message
helpMessage = """
===============================================================================
EASI_FullLength16S_pipleline | version ${params.version}
===============================================================================
Usage: nextflow run ikmb/EASI_FullLength16S 

Mandatory parameters:
--reads       Location to the hifi.bam file
--ccga        Location of ccga_file.tsv
--movieid     Set a run ID

Optional:
--outdir      Set a custom output directory. Default is "results"
-work-dir     Set a custom working directy. Default is "work"
--keep_demul  Will copy demultiplexed files into the results directory.
-resume       Will resume an earlier run.

Advanced parameters:
-profile      Set a profile for nextflow. Default is "medcluster", add new profiles or run locally with '-profile local'.
--scratch     Set if scratch should be used
--allbarcodes Use a different barcodes files
"""

params.help = false

// Show help when needed
if (params.help){
    log.info helpMessage
    exit 0
}

def summary = [:]

run_name = ( params.run_name == false) ? "${workflow.sessionId}" : "${params.run_name}"

include {	EASI } from './workflows/main' 
//params(params)

workflow {

	EASI()

}

workflow.onComplete {
  log.info "========================================="
  log.info "Duration:		$workflow.duration"
  log.info "========================================="
}

