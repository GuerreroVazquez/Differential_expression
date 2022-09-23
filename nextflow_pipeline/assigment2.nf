process INDEX{

    publishDir "${params.outdir}/index", mode: 'copy'

    input:
    file(fasta) from transcriptome_created


    output:
    file("${fasta.baseName}.fa.inx") into index_created

    script:
    """    
    kalisto index  -i "${fasta.baseName}.inx" $fasta
    """

}
