#!/usr/bin/env nextflow

/*
===============================================================
 nf-core/nanopore
===============================================================
 Nanopore RNA-Seq Analysis Pipeline. Started Nov 2018.
 #### Homepage / Documentation
 https://github.com/
 #### Authors
 Tsz Wai Pang <chris_pang429@outlook.com>
---------------------------------------------------------------
*/
 
def helpMessage() {
    log.info"""
    =======================================================
                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\'
        |\\ | |__  __ /  ` /  \\ |__) |__         }  {
        | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                              `._,._,\'
     nf-core/nanopore : v${workflow.manifest.version}
    =======================================================
    Usage:
    
    nextflow run nf-core/main.nf --query --genome --chromSize --prefix --tool --reference --refSeqBed --refflat 

    Mandatory arguments:
      --query                        Path to input data (must be surrounded with quotes)
      --genome                       Path to reference genome
      --chomSize                     Chromosome size file  
      --prefix                       Prefix for file names
      --tool                         Specify preference for GraphMap or Minimap2 
      --reference                    Specify genome or transcriptome reference 
      --refSeqBed                    Specify refseq standard bed file
      --refflat                      Specify refFlat file
      --profile                      Configuration profile to use. Can use multiple (comma separated) Available: standard, conda, docker, singularity, awsbatch, test  


    Options:
      --publishDir                  Copies the process output files to a specified folder


    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.


    AWSBatch options:
      --awsqueue                    The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

/*
 * Defines the pipeline input parameters to specify various inputs
 */
params.genome = ""
params.query = ""
params.chromSize = ""
params.tool = ""
params.prefix = ""
params.reference = ""
params.refSeqBed = ""
params.refflat = ""

/*
 * Creates objects given string parameters
 */
genome = file(params.genome)
query = file(params.query)
chromSize = file(params.chromSize)
prefix = params.prefix
tool = params.tool
reference = params.reference
refSeqBed = files(params.refflat)
refflat = file(params.refflat)

/*
 * Process 1. Maps the reads to genome/transcriptome using Minimap2/GraphMap
 */
process mapping {
    input:
    file query_file from query
    file genome_file from genome
            
    output:
    file "${prefix}.sam" into sam_files

    script:
    if(tool == 'graphmap') {
        if(reference == 'genome'){
        """
        graphmap align -r ${genome_file} -d ${query_file} -o ${prefix}.sam
        """

        } 
        else {
        """
        graphmap align -r ${genome_file} --gtf ${features} -d ${query_file} -o ${prefix}.sam
        """
        }
    }
    else {
        if(reference == 'genome'){
        """
        minimap2 -ax splice -uf -k14 -I ${task.memory.toGiga()} -t ${task.cpus} ${genome_file} ${query_file} > ${prefix}.sam
        """ 
        } 
        else {
        """
        minimap2 -ax splice -uf -k14 -I ${task.memory.toGiga()} -t ${task.cpus} ${genome_file} ${query_file} > ${prefix}.sam
        """
        }
    }
}

/*
 * Process 2. Generate sorted bam file
 */
process generateBam {
    input:
    file sam_file from sam_files

    output: 
    file "${prefix}.bam" into bam_files

    """
    samtools view -Sb ${sam_file} | samtools sort -o ${prefix}.bam
    """
} 

/*
 * Process 3. Generate bed file
 */
process generateBed {
    input:
    file bam_file from bam_files

    output:
    file "${prefix}.bed" into bed_files

    """
    bamToBed -bed12 -cigar -i ${bam_file} > ${prefix}.bed
    """
}

/*
 * Process 4. Remove non-standard chromosomes from bed file
 */
process removeNonStandardChrom {
    input:
    file bed_file from bed_files

    output:
    file "${prefix}.stdChr" into stdChr_files

    """
    Rscript '$baseDir/modifyBedFile-Rscript_R.R' -i ${bed_file} -o ${prefix}.stdChr -t removeNonStandardChromosomesNCBI
    """
}

/*
 * Process 5. Convert to UCSC chromosome names
 */
process bedToUCSCNames {
    input:
    file stdChr_file from stdChr_files

    output:
    file "${prefix}_stdChr.ucsc" into bedUCSC_files

    """
    Rscript '$baseDir/changeSeqLevelStyle.R' -i ${stdChr_file} -o ${prefix}_stdChr.ucsc -s UCSC
    """
}

/*
 * Process 6. Generate bigBed file
 */
process generateBigBed {
    input:
    file bedUCSC_file from bedUCSC_files

    output:
    file "${prefix}.bb" into bigBed_files

    """
    bedToBigBed ${bedUCSC_file} $chromSize ${prefix}.bb
    """
}

/*
 * Process 7. Generate bedgraph file
 */
process generateBedGraph {
    input:
    file bam_file from bam_files

    output:
    file "${prefix}.bedgraph" into bedgraph_files

    """
    genomeCoverageBed -split -bg -i ${bam_file} > ${prefix}.bedgraph
    """ 
}

/*
 * Process 8. Convert bedgraph to UCSC names
 */
 process bedgraphToUCSC {
    input:
    file bedgraph_file from bedgraph_files

    output:
    file "${prefix}.ucsc" into bedgraphUCSC_files

    """
    Rscript '$baseDir/changeSeqLevelStyle.R' -i ${bedgraph_file} -o ${prefix}.ucsc -s UCSC
    """
 }
  
/*
 * Process 9. Sort bedgraph.ucsc file
 */
process bedSort {
    input:
    file bedgraph_file from bedgraphUCSC_files
    
    output:
    file "${prefix}.sorted" into bedgraph_sorted_files

    """
    bedSort ${bedgraph_file} ${prefix}.sorted
    """
}

/*
 * Process 10. Generate bigWig file
 */
process generateBigWig {
    input:
    file bedgraph_sorted_file from bedgraph_sorted_files

    output:
    file "${prefix}.bw" into bigWig_files

    """
    bedGraphToBigWig ${bedgraph_sorted} $chromSize ${prefix}.bw
    """   
}

/* 
 * Process 11. Visualize fastq files with NanoPlot
 */
process visualizeFastq {
    input:
    file query_file from query

    output:
    file('*') into nanoplot_files

    """
    NanoPlot -t 2 --fastq ${query} --plots hex dot --N50
    """
}

/*
 * Process 12. Upload visualization files to S3 bucket
 */
process uploadToS3 {
    input:
    file bigBed_file from bigBed_files
    file bigWig_file from bigWig_files

    """
    aws s3 cp ${bigBed_file} s3://chrisp-publicdata.store.genome.sg/${prefix}/
    aws s3 cp ${bigWig_file} s3://chrisp-publicdata.store.genome.sg/${prefix}/
    """
}

/*
 * Process 13. Index bam file generated from Minimap2
 */
process bamIndex {
    input:
    file bam_file from bam_files

    output:
    file "${prefix}.bai" into bamIndex_files
    
    """
    samtools index ${bam_file}
    """    
 }

/*
 * Process 14. Provide statistics on a bam file using Samtools flagstat tool 
 */
process flagstat {
    input:
    file bam_file from bam_files

    output:
    file "${prefix}.flagstat" into flagstat_files

    """
    samtools flagstat ${bam_file} > ${prefix}.flagstat
    """
}

/*
 * Process 15. Compute gene body coverage
 */
process geneBodyCov {
    input:
    file bam_file from bam_files
    file refSeqBed_file from refSeqBed
    file bamIndex_file from bamIndex_files

    output:
    file "${prefix}.geneBodyCoverage*" into geneBodyCov_files

    """
    geneBody_coverage.py -r ${refSeqBed_file} -i ${bam_file} -o ${prefix}
    """
}

/* 
 * Process 16. Access transcript coverage with Picard
 */
process picard {
    input:
    file bam_file from bam_files
    file refflat_file from refflat

    output:
    file "${prefix}.metric" into picard_files

    """
    picard CollectRnaSeqMetrics REF_FLAT=${refflat} STRAND=NONE INPUT=${bam_file} OUTPUT=${prefix}.metric VALIDATION_STRINGENCY=LENIENT
    """
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
