#!/usr/bin/env nextflow


// Variables and defaults
params.help = null
params.gzip = true
params.outputFolder = "./results/"
params.ext = "_paired.fastq.gz"

// Genome options
params.genomeName = "GRCh38"
params.SAindex = 14
params.index = false

//Trimming options
params.clip_r1 = 13
params.clip_r2 = 13
params.three_prime_clip_r1 = 2
params.three_prime_clip_r2 = 2
params.quality = 30
params.length = 50

log.info "------------------------------------------------------------------"
log.info "RNA Seq Pre-Processing and STAR Alignment Pipeline. Messaoudi Lab."
log.info "------------------------------------------------------------------"
if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "                     USAGE                              "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "./nextflow run rnaSeq.nf --inputFolder path/to/fastq/folder --genome /path/to/reference_genome/fasta --annot /path/to/genome/annot"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--inputFolder          FOLDER               Folder containing fastq files"
    log.info "--genome               FILE                 Reference genome fasta (.fasta or .fa) file"
    log.info "--annot                FILE                 Reference genome annotation file (.gtf)"
    log.info '--genomeDir            FOLDER               If genome indexing occurred previously,'
    log.info '                                            this folder contains the indexed genome files.'
    log.info '                                            Otherwise, this is the name of the folder where'
    log.info '                                            indexed genome files will be stored upon creation.'
    log.info ""
    log.info "Optional arguments:"
    log.info "--outputFolder         FOLDER               Output directory for result directories (default: ./results)"
    log.info "--genomeName           STRING               Name of reference genome (for labelling output files)"
    log.info "                                            (default: GRCh38)"
    log.info "--SAindex              INTEGER              Length (bases) of the SA pre-indexing string. Must be scaled"
    log.info "                                            down for smaller genomes (default: 14)"
    log.info "--gzip                 BOOLEAN              Indicates whether or not the input fastq files are gzipped (default: true)"
    log.info '--ext                  STRING               Extension of files (default: _paired.fastq.gz)'
    log.info '--clip_r1              INTEGER              Number of bp to remove from 5` end of read one (default: 13)'
    log.info '--clip_r2              INTEGER              Number of bp to remove from 5` end of read two (default: 13)'
    log.info '--three_prime_clip_R1  INTEGER              Number of bp to remove from 3` end of read one (default: 2)'
    log.info '--three_prime_clip_R2  INTEGER              Number of bp to remove from 3` end of read two (default: 2)'
    log.info '--quality              INTEGER              Phred score cutoff for trimming (default: 30)'
    log.info '--length               INTEGER              Minimum length cutoff for trimming (default: 50)'
    log.info '--index                BOOLEAN              Indicates whether or not the user requires the genome to be'
    log.info '                                            indexed (default: false)'
    log.info ""
    log.info "Flags:"
    log.info "--help                                      Display this message"
    log.info ""
    exit 0
}   

// Check for required parameters
if (!params.inputFolder || !params.genome || !params.annot || !params.genomeDir){
    exit 1, "Parameters '--inputFolder', '--genome', '--annot', and '--genomeDir' are required to run the pipeline."
}

Channel
    .fromFilePairs( params.inputFolder+'/*{1,2}'+params.ext)
    .ifEmpty { error "Cannot find any file with extension ${params.ext} in: ${params.inputFolder} " }
    .set { inFiles_trimming }

Channel
    .fromFilePairs( params.inputFolder+'/*{1,2}'+params.ext)
    .ifEmpty { error "Cannot find any file with extension ${params.ext} in: ${params.inputFolder} " }
    .set { inFiles_fastqc }

/*
 *
 * STEP 0.1 - PRE-PROCESSING: FASTQC
 *
*/

process fastqc {
  tag { name }

  publishDir "${params.outputFolder}/fastqc", mode: 'copy'

  input:
  tuple val(name), file(f) from inFiles_fastqc

  output:
  file  "*_fastqc.{zip,html}" into fastqc_results

  shell:
  '''
  fastqc -f fastq !{f.get(0)} !{f.get(1)}
  '''
}


/*
 *
 * STEP 0.2 - PRE-PROCESSING: TRIMMING
 *
*/

process trim_galore {
    tag { name}


    publishDir "${params.outputFolder}/trim_galore", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("trimming_report.txt") > 0) "trim_logs/$filename"
            else "trimmed_$filename"
        } 

    input:
    tuple val(name), file(f) from inFiles_trimming
    
    output:
    tuple val(name), file(f) into trimmed_reads
    file '*trimming_report.txt' into trimgalore_results, trimgalore_logs

    shell:
    """
    trim_galore --paired !{f.get(0)} !{f.get(1)} \
    --clip_r1 !{params.clip_r1} \
    --clip_r2 !{params.clip_r2} \
    --three_prime_clip_r1 !{params.three_prime_clip_r1} \
    --three_prime_clip_r2 !{params.three_prime_clip_r2} \
    --quality !{params.quality} \
    --length !{params.length} \
    """
}

/*
*
* STEP 0.3 - PRE-PROCESSING: Index Reference Genome (if necessary)
*
*/

process index {

    when:
    params.index == true

    shell:
    """
    STAR --runThreadN 12 \
    --runMode genomeGenerate \
    --genomeDir !{params.genomeDir} \
    --genomeFastaFiles !{params.genome} \
    --sjdbGTFfile !{params.annot} \
    --sjdbOverhang 99 \
    --genomeSAindexNbases !{params.SAindex}
    """
}


/*
*
* STEP 1.0 - ANALYSIS: Perform Alignment
*
*/

process align {
    tag { name }


    publishDir "${params.outputFolder}/STAR_alignment", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".bam") > 0) "aligned_sorted_bams/$filename"
            else "logs/$filename"
        } 

    input:
    tuple val(name), file(reads) from trimmed_reads
    
    output:
    file '*.bam' into aligned_sortedByCoord_bams
    file '*Log*' into alignment_sortedByCoord_logs

    shell:
    if( params.gzip == true )
        """
        STAR --runThreadN 12 \
        --genomeDir !{params.genomeDir} \
        --readFilesIn !{reads.get(0)} !{reads.get(1)} \
        --outFileNamePrefix !{name}_!{params.genomeName}_ \
        --outSAMtype BAM SortedByCoordinate \
        --sjdbGTFfile !{params.annot} \
        --readFilesCommand zcat

        """
    else
        """
        STAR --runThreadN 12 \
        --genomeDir !{params.genomeDir} \
        --readFilesIn !{reads.get(0)} !{reads.get(1)} \
        --outFileNamePrefix !{name}_!{params.genomeName}_ \
        --outSAMtype BAM SortedByCoordinate \
        --sjdbGTFfile !{params.annot}

        """

}

/*
*
* STEP 1.1 - ANALYSIS: Samtools sort by name
*
*/

process samtoolsSort {

    publishDir "${params.outputFolder}/samtools_sortedByName", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".bam") > 0) "$filename"
            else "logs/$filename"
        }
    input:
    file f from aligned_sortedByCoord_bams

    output:
    file f into name_sorted_bams

    shell:
    """
    samtools sort !{f} -o !{f.baseName}_sortedByName_aligned.bam
    """

}

workflow.onComplete {
    println (workflow.success ? """
        Pipeline execution summary
        --------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : $workflow.exitstatus}
        """
     )
}
