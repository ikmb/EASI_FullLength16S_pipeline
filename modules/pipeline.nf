process biosample {
    scratch params.scratch
    input:
    output:
        path("biosample.csv")
    shell:
    """
        echo "Barcode,Bio Sample" > biosample.csv

        while read line; do
        id=\$(echo \$line | awk '{print \$1}')
        fwdseq=\$(echo \$line | awk '{print \$NF}')
        revseq=\$(echo \$line | awk '{print \$(NF-1)}')
        if [ "\$id" == "subject_name" ]; then continue; fi
        fwdid=\$(grep -w \$fwdseq !{params.allbarcodes} -B 1 | head -n 1 | tr -d '>')
        revid=\$(grep -w \$revseq !{params.allbarcodes} -B 1 | head -n 1 | tr -d '>')
        echo "\${fwdid}--\${revid},\${id}" >> biosample.csv
        done < ${params.ccga}
    """
}

process barcode_extraction {
    scratch params.scratch
    input:
        path(biosample)
    output:
        path("barcodes.fasta")
    shell:
    """
        cut -d ',' -f 1 !{biosample}  | sed 's/\\-\\-/\\n/' | sort | uniq  > test.fasta
        awk '(NR>1) {print}' test.fasta | xargs -I {} grep -A 1 -w {} !{params.allbarcodes} > barcodes.fasta
        #rm test.fasta
    """
}

process demux_dada2 {
    scratch params.scratch
    //scratch false

    publishDir "${params.outdir}/${params.movieid}/", pattern: "filtered/*", mode: 'copy'
    publishDir "${params.outdir}/${params.movieid}/", pattern: "RDS/*", mode: 'copy'
    publishDir "${params.outdir}/${params.movieid}/", pattern: "tables/*", mode: 'copy'
    publishDir "${params.outdir}/${params.movieid}/", pattern: "allout*", mode: 'copy'
    publishDir "${params.outdir}/${params.movieid}/", pattern: "noprimers/*", mode: 'copy'
    if(params.keep_demul){
       publishDir "${params.outdir}/${params.movieid}/", pattern: "raw/*", mode: 'copy'
    }

    input:
        path(biosample)
        path(barcodes)
    output:
        path("*")
        path("**/*")
    shell:
    """
        mkdir -p raw
        lima -j !{task.cpus} \
                --preset HIFI-ASYMMETRIC \
                --biosample-csv !{biosample} \
                --max-input-length 2000 \
                --split-named !{params.hifireads} !{barcodes} ./raw/!{params.movieid}.fastq 

       ### gzip files
        parallel -j ${task.cpus} gzip {} ::: \$(ls raw/*.fastq)

        Rscript !{baseDir}/bin/dada2_PacBio_EASI_script.R !{params.movieid} !{task.cpus} 
    """
}