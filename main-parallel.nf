/*
 * Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'RNA-Toy'.
 *
 *   RNA-Toy is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   RNA-Toy is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with RNA-Toy.  If not, see <http://www.gnu.org/licenses/>.
 */
 
/* 
 * Proof of concept Nextflow based RNA-Seq pipeline
 * 
 * Authors
 * Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * Emilio Palumbo <emiliopalumbo@gmail.com> 
 */ 

 
/*
 * Defines some parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 */
params.reads = "$baseDir/data/10{N,T}_{1,2}.raw.fq"
params.annot = "$baseDir/data/chr19.gff"
params.genome = "$baseDir/data/chr19.fa"
params.outdir = 'results'

log.info "R N A T O Y   P I P E L I N E    "
log.info "================================="
log.info "genome             : ${params.genome}"
log.info "annotat            : ${params.annot}"
log.info "reads              : ${params.reads}"
log.info "outdir             : ${params.outdir}"

/*
 * the reference genome file
 */
genome_file = file(params.genome)
annotation_file = file(params.annot)
 
/*
 * Create the `read_pairs` channel that emits tuples containing three elements:
 * the pair ID, the first read-pair file and the second read-pair file 
 */
Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs } 

    
/*
* A channel to split the reads in smaller files.
*/
 
 process split {
 
    input:
    set pair_id, file(reads) from read_pairs
    
    output:
    set pair_id, file('10*') into split_read_pairs
    
    """
    split -n 4 ${reads} ${pair_id}
    """
 
 }
 
/*
 * Step 2. Maps each read-pair by using Tophat2 mapper tool
 */
process mapping {
    tag "$pair_id"
     
    input:
    file 'genome.index.fa' from genome_file 
    file annotation_file
    set pair_id, file(reads) from split_read_pairs
 
    output:
    set pair_id, "accepted_hits.bam" into bam
 
    """
    tophat2 -p ${task.cpus} --GTF $annotation_file "/home/bioitcore/NexFlow/rnatoy/data/chr19" ${reads}
    mv tophat_out/accepted_hits.bam .
    """
}
  
/*
 * Step 3. Assembles the transcript by using the "cufflinks" tool
 */
process makeTranscript {
    tag "$pair_id"
    publishDir params.outdir, mode: 'copy'  
       
    input:
    file 'anno.gtf' from annotation_file
    set pair_id, file(bam_file) from bam
     
    output:
    set pair_id, file('transcript_*.gtf') into transcripts
 
    """
    cufflinks --no-update-check -q -p ${task.cpus} -G anno.gtf ${bam_file}
    mv transcripts.gtf transcript_${pair_id}.gtf
    """
}
 
workflow.onComplete { 
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
