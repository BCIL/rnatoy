
 
/*
 * Defines some parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 * params.reads = "$baseDir/data/10{N,T}_{1,2}.raw.fq"
 */

params.reads = "$baseDir/data/10*1000.fq"
params.annot = "$baseDir/data/chr19.gff"
params.genome = "$baseDir/data/chr19.fa"
params.outdir = 'results'


/*
 * the reference genome file
 */
genome_file = file(params.genome)
annotation_file = file(params.annot)
 

 
/*
 * A channel to simply emit file for splitting.
 */
EmitFiles = Channel.create()
EmitFiles
    .fromPath( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set {read_files}

    
/*
* A process to split the reads in smaller files.
*/
 
 process split {
 
    input:
    file reads from read_files
    
    output:
    file "${reads}.*" into split_read_files
    
    """
    split -n 4 -d ${reads} ${reads}
    """
 
 }
 
 /*
 * Create the `read_pairs` channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file 
 */
CapturePairs = Channel.create()
CapturePairs
    .from(split_read_files)
    .fromFilePairs('*{N.T}*{1,2}*')
    .set(split_read_pairs)
 
 
workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )  }