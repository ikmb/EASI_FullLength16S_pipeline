include { biosample; barcode_extraction; demux_dada2 } from '../modules/pipeline'
//include { SOFTWARE_VERSIONS } from '../modules/software_versions'

workflow EASI {	
	main:
		if(params.hifireads && params.ccga && params.movieid){
			//pipeline_PacBio()

			biosample()
			barcode_extraction(biosample.out)
			demux_dada2(
				biosample.out,
				barcode_extraction.out
			)
			//dada2(demux.out)

		} else {
			exit 1, "Invalid input, check that the mandatory inputs hifireads, ccga and movieid are set!"
		}
}

	/*
	ch_software_versions = Channel.empty()

	ch_software_versions = ch_software_versions.mix(FASTP.out.version).ifEmpty(null)

	ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

	SOFTWARE_VERSIONS(		ch_software_versions	)		
	*/

