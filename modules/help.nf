/*
* Help message
*/

def helpMessage() {
    println ''
    log.info """
    nucMACC   P I P E L I N E
    =============================
    Usage:

    nextflow run uschwartz/nucMACC --csvInput 'path2csvFile' --outDir 'path2outDir' --genomeIdx 'path2bowtie2_idx'

    Mandatory arguments:
      --csvInput        [string] Path to comma-separated file containing information about the samples in the experiment (see ./toyData/input.csv as example) template provided in ./input_template.csv
      --genomeIdx       [string] Path and prefix of bowtie2 index (minus .X.bt2)
      --genomeSize      [integer] effective genome size, defined as the length of the mappable genome. Used for normalisation (default: 162367812 (dm3)) 

    optional arguments:
      --outDir          [string] Name of output directory, which will be created (default: ~/nucMACC_test/)
      --blacklist       [string] A BED file containing regions that should be excluded from all nucleosome analysis (default: false)

     """.stripIndent()
     println ''
}
