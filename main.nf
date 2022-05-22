#!/usr/bin/env nextflow

import java.io.File
import org.apache.commons.lang.StringUtils

nextflow.enable.dsl=2

project_dir = projectDir

def helpMessage() {
    log.info"""
    Description:
        An easy to use pipeline to separate endosymbiont genomes from their host's
    Pipeline summary:
        1. Trimming using TrimGalore!
        2. Read quality control using FastQC
        3. De Novo assembly using megahit
        4. Filtering endosymbiont genome using blastn
        5. Filtering host mitogenome using blastn
        6. Read mapping for coverage estimate using bowtie2
        7. Coverage estimate
        8. Assembly quality assessment using BUSCO
    Usage:
        nextflow run main.nf --reads '*_R{1,2}\\.fastq.gz' --endosymbiont_reference '*_endosymRef\\.fna' --host_reference '*_hostRef\\.fna'
        
    Mandatory arguments:
        --reads             path to one or more sets of paired-ended reads (valid
                            file types: .fastq.gz', '.fq.gz', '.fastq', or '.fq')
        --endosymbiont_reference
                            path to one or more reference genomes for the endosymbiont
                            assembly (valid file type extensions: '.fa', '.fna', '.fasta', '.faa')
    Input/output options:
        --output            path to a directory which the results are written to
                            (default: $params.output)
    Resource allocation:
        --memory            memory limit for the assembly step in GB (default:
                            $params.max_memory)
        --threads           maximum number of threads to be used by the pipeline
                            (default: '$params.max_cpus')
    Flow control:
        --endosymbiont_only
                            skip processing of reads not belonging to the endosymbiont (default: $params.endosymbiont_only)
        --skip_coverage     skip coverage estimate step (default: $params.skip_coverage)
        --skip_trimming     skip trimming step (default: $params.skip_trimming)
        --skip_qc           skip reads quality assessment (default: $params.skip_qc)
        --skip_endosymbiont_assembly
                            skip endosymbiont assembly step (default: $params.skip_endosymbiont_assembly)
        --skip_assembly_quality
                            skip assembly quality assessment (default: $params.skip_assembly_quality)
    Miscellaneous:
        --help              display this help message and exit
        --version           display the pipeline's version number and exit
    """.stripIndent()
}

def versionNumber() {
    log.info"symbiontDivider ~ version $workflow.manifest.version"
}

// Display the version number on request
if ( params.version ) exit 0, versionNumber()

// Display a help message upon request
if ( params.help ) exit 0, helpMessage()

// Input validation reads
if ( params.reads == null) {
    exit 1, "Missing mandatory argument '--reads'\n" +
            "Launch this workflow with '--help' for more info"
}


// Input validation endosymbiont reference
if ( params.endosymbiont_reference == null) {
    exit 1, "Missing mandatory argument '--endosymbiont_reference'\n" +
            "Launch this workflow with '--help' for more info"
}

// Creation of read pair channel with file extension filter and check if empty
rawReads = Channel
    .fromFilePairs( params.reads, size: 2, type: 'file' )
    .filter { it =~/.*\.fastq\.gz|.*\.fq\.gz|.*\.fastq|.*\.fq/ }
    .ifEmpty { exit 1,
             "No FASTQ files found with pattern '${params.endosymbiont_reference}'\n" +
             "Escape dots ('.') with a backslash character ('\\')\n" +
             "Try enclosing the path in single-quotes (')\n" +
             "Valid file types: '.fastq.gz', '.fq.gz', '.fastq', or '.fq'" }

// Setting job name
extensions = ['_R{1,2}', '{1,2}', '\\.fastq\\.gz', '\\.fq\\.gz', '.fastq.gz', '.fq.gz', '\\.fastq', '\\.fq', '.fq', '.fastq' ]
String[] extensions = extensions.toArray(new String[extensions.size()])
replacement_list = ['', '', '', '', '', '', '', '', '', '']
String[] replacement_list = replacement_list.toArray(new String[replacement_list.size()])
params.job_name = StringUtils.replaceEach(new File(params.reads).getName(), extensions, replacement_list)



process RAWQC {

    /* 
        Process Description:
        Quality control of raw reads with FastQC
    */

    // Copies output file in output folder
    publishDir "${params.output}/${params.job_name}/quality_control", mode: 'copy'

    // Process label
    label 'fast'

    // If this process fails, it does not end the whole pipeline
    errorStrategy 'finish'

    // Name of read files
    tag "${params.job_name}"

    input:
    // A tuple containing the name of the raw read files and the files themself
    tuple val(name), path(reads)

    output:
    // FastQC outputs the results as .html files -> those are output
    path '*.html'

    when:
    // This process only executes when QC is not skipped
    ! skip_qc

    script:
    // FastQC command with half of input threads and in quiet mode (see FastQC documentation for details)
    """
    fastqc --threads ${task.cpus} --quiet ${reads[0]} ${reads[1]}
    """

}

process TRIMMING {

    /* 
        Process Description:
        Trimming of read files using TrimGalore!
    */
    
    // Process label
    label 'normal'

    // Name of read files
    tag "${params.job_name}"

    input:
    // A tuple containing the name of the raw read files and the files themself
    tuple val(name), path(reads)

    output:
    // TrimGalore! outputs the trimmed reads in .fg files that are output in a tuple comined with the name of the files
    tuple val("${reads[0].baseName}"), path('*.fq*'), emit: trimmed_reads

    when:
    // This process is only executed when trimming is not skipped
    ! params.skip_trimming

    script:
    // TrimGalore! command with paired read option (see TrimGalore! documentation for details)

    flagsTrimming = "--quality $params.trim_quality \
--length $params.trim_length --cores ${task.cpus}"
    if ( params.trim_phred64 )
      flagsTrimming += " --phred64"
    if ( params.trim_clip_R1 )
      flagsTrimming += " --clip_R1 $params.trim_clip_R1"
    if ( params.trim_three_prime_clip_R1 )
      flagsTrimming += " --three_prime_clip_R1 $params.trim_three_prime_clip_R1"
    if ( params.trim_clip_R2 )
      flagsTrimming += " --clip_R2 $params.trim_clip_R2"
    if ( params.trim_three_prime_clip_R2 )
      flagsTrimming += " --three_prime_clip_R2 $params.trim_three_prime_clip_R2"
    flagsTrimming += " --paired"
    commandTrimming = "trim_galore $flagsTrimming ${reads[0]} ${reads[1]}"

    """
    $commandTrimming
    """

}

process TRIMMEDQC {

    /* 
        Process Description:
        Quality control of trimmed reads with FastQC
    */

    // Copies output file in output folder
    publishDir "${params.output}/$params.job_name/quality_control", mode: 'copy'

    // Process label
    label 'fast'

    // If this process fails, it does not end the whole pipeline
    errorStrategy 'finish'

    // Name of read files
    tag "$params.job_name"

    input:
    // A tuple containing the name of the trimmed read files and the files themself
    tuple val(name), path(reads)

    output:
    // FastQC outputs the results as .html files -> those are output
    path '*.html'

    when:
    // This process is only executed if QC and trimming are not skipped
    ! skip_qc && ! params.skip_trimming

    script:
    // FastQC command with half of input threads and in quiet mode (see FastQC documentation for details)
    """
    fastqc --threads ${task.cpus} --quiet ${reads[0]} ${reads[1]}
    """

}

process DENOVOASSEMBLY {

    /* 
        Process Description:
        De novo assembly of reads using SPAdes
    */

    // Process label
    label 'big_mem'

    // Name of read files
    tag "${params.job_name}"

    input:
    // A tuple containing the name of the raw/trimmed read files and the files themself
    tuple val(name), path(reads)
    val kmers

    output:
    // SPAdes outputs the assembled genome as a .fasta file called "scaffolds.fasta" -> output of the process
    tuple val("${reads[0].baseName}"), path('scaffolds.fasta')

    script:
    // Make list of kmers SPAdes-compatible ([a, b, c] -> "a,b,c")
    kmersFormatted = kmers.toString().replaceAll("[ \\[\\]]", "")
    additionalSpadesFlags = ""
    if ( params.meta )
        additionalSpadesFlags += "--meta \\\n"

    // SPAdes command using input threads and memory (for details see SPAdes documentation)
    """
    spades.py -o . -1 ${reads[0]} -2 ${reads[1]} --threads ${task.cpus} --memory ${task.memory.toGiga()} -k "${kmersFormatted}" ${additionalSpadesFlags} --disable-gzip-output
    """
}

process ENDOSYMBIONTCONTIGFILTERING {

    /* 
        Process Description:
        Extraction of contigs belonging to the endosymbiont using blastn and grep
    */

    // Copies output file in output folder
    publishDir "${params.output}/$params.job_name/endosymbiont_assembly", mode: 'copy'
    
    // Process label
    label 'normal'

    // Name of read files
    tag "$params.job_name"

    input:
    // A tuple containing the name of the raw/trimmed read files and the contigs assembled before
    tuple val(name), path(contigs)
    // A fasta file containing the endosymbiont reference genome
    // path endosymbiont_reference

    output:
    // The process outputs a tuple with the reads name and a .fa file containing all the contigs belonging to the endosymbiont genome
    tuple val(name), path('endosymbiont_genome.fa'), emit: endosym_mapped

    script:
    /*
    Script description:
    1. A blast database from the reference genome is created
    2. A blastn search with the de novo assembled contigs is performed and found contig ids are saved in a .txt
    3. Based on the contig ids, contigs are grepped from the de novo assembled contigs using bfg
    */
    """
    cat $contigs | bfg "cov_([1-9][3-9][0-9]*|[1-9][0-9][0-9]{1,}|[2-9][0-9])\\.[0-9]+" > contigs.fa
    makeblastdb -in ${params.endosymbiont_reference} -title endosymbiont -parse_seqids -dbtype nucl -hash_index -out db
    blastn -query contigs.fa -db db -outfmt "10 qseqid" > seqid.txt
    cat contigs.fa | bfg -F -f seqid.txt > endosymbiont_genome.fa
    """
}

process HOSTMITOGENOMEFILTERING {

    /* 
        Process Description:
        Extraction of contig/s belonging to the host mitogenome using blastn and grep
    */

    // Copies output file in output folder
    publishDir "${params.output}/$params.job_name/host_assembly", mode: 'copy'
    
    // Process label
    label 'normal'

    // If this process fails, it does not end the whole pipeline
    errorStrategy 'finish'

    // Name of read files
    tag "$params.job_name"

    input:
    // A tuple containing the name of the raw/trimmed read files and the contigs assembled before
    tuple val(name), path(host_assembled)

    output:
    // The process outputs a tuple with the reads name and a .fa file containing all the contigs belonging to the host mitogenome
    tuple val(name), path('mitogenome.fa'), emit: host_filtered
    tuple val(name), path('mitogenome_candidates*')

    when:
    // This process is only executed if the endosymbont only mode is not selected
    ! params.endosymbiont_only 


    script:
    /*
    Script description:
    1. Create empty files for intermediate storage
    2. Create a blast database from a sequence for the cox1 gene (mitogenome exclusive gene)
    3. Iterate from 11 to 25
        3.1. Concatenate the previous found reads to prev_seqid.txt
        3.2. Blastn search with word size determined by iteration number, using de novo assembled reads as query
        3.3. Determined seqids are made unique and cat into unique_seqid.txt
        3.4. If the unique_seqid.txt is empty -> grep contigs based on previously found seqids by bfg -> break iteration
        3.5. If the unique_seqid.txt has 1 entry -> grep corrisponding contig by bfg -> break iteration
    */
    """
    touch mitogenome.fa
    touch prev_seqid.txt
    touch unique_seqid.txt
    touch possible_mitogenomes.fa
    cat $host_assembled | bfg "cov_[1-9][0-9][0-9]{1,}\\.[0-9]+" > possible_mitogenomes.fa
    makeblastdb -in ${params.mitogenome_bait} -title cox1 -parse_seqids -dbtype nucl -hash_index -out db
    echo "blastdb created"
    for i in {${params.min_blast_wordsize}..${params.max_blast_wordsize}..1}
      do
        echo "starting iteration with word size \$i"
        cat unique_seqid.txt > prev_seqid.txt
        blastn -query possible_mitogenomes.fa -db db -outfmt "10 qseqid" -word_size \$i > seqid.txt
        echo "blastn complete"
        cat -n seqid.txt | sort -uk2 | sort -nk1 | cut -f2- | cat > unique_seqid.txt
        echo "made seqids unique"
        cat possible_mitogenomes.fa | bfg -f unique_seqid.txt > "mitogenome_candidates_wordsize_\$i.fa"
        if [[ \$(wc -l unique_seqid.txt) = "0 unique_seqid.txt" ]];
        then
          cat possible_mitogenomes.fa | bfg -f prev_seqid.txt > mitogenome.fa
          echo "multiple possible mitogenomes found"
          break
        fi
        if [[ \$(wc -l unique_seqid.txt) = "1 unique_seqid.txt" ]];
        then
          cat possible_mitogenomes.fa | bfg -f unique_seqid.txt > mitogenome.fa
          echo "mitogenome found"
          break
        fi
      done
    echo "process successful"
    """
}

process ENDOSYMBIONTGENOMEQUALITY {

    /* 
        Process Description:
        Quality assessment of found endosymbiont genome using BUSCO
    */

    // Copies output file in output folder
    publishDir "${params.output}/$params.job_name/endosymbiont_assembly", mode: 'copy'

    // Process label
    label 'normal'

    // If this process fails, it does not end the whole pipeline
    errorStrategy 'finish'

    // Name of read files
    tag "$params.job_name"

    input:
    // A tuple containing the name of the raw/trimmed read files and the contigs belonging to the endosymbiont genome
    tuple val(name), path(endosym)

    output:
    // All files created by BUSCO are output
    path 'qc/*'

    when:
    // This process is only executed if the last quality assessment is not skipped
    ! params.skip_assembly_quality

    script:
    // BUSCO command (see BUSCO documentation for details)
    """
    busco -i $endosym -m genome -o qc -l rickettsiales_odb10 --offline --download_path $project_dir/seqs/busco_data
    """

}

process HOSTMITOGENOMEQUALITY {

    /* 
        Process Description:
        Quality assessment of found host mitogenome using quast
    */

    // Copies output file in output folder
    publishDir "${params.output}/$params.job_name/host_assembly", mode: 'copy'

    // Process label
    label 'normal'

    // If this process fails, it does not end the whole pipeline
    errorStrategy 'finish'

    // Name of read files
    tag "$params.job_name"

    input:
    // A tuple containing the name of the raw/trimmed read files and the contigs belonging to the host mitogenome
    tuple val(name), path(host)

    output:
    // All files created by quast are output
    file '*'

    when:
    // This process is only executed if the last quality assessment is not skipped and the endosymbiont only mode is not activated
    ! params.skip_assembly_quality && ! params.endosymbiont_only

    script:
    // quast command (see quast documentation for details)
    """
    quast.py $host
    """

}

process READMAPPINGFORCOVERAGE{

    /* 
        Process Description:
        Mapping of raw/trimmed reads onto the found endosymbiont genome using bowtie2
    */

    // Process label
    label 'big_mem'

    // If this process fails, it does not end the whole pipeline
    errorStrategy 'finish'

    // Name of read files
    tag "$params.job_name"

    input:
    // A tuple containing the name of the raw/trimmed read files and the files themself
    tuple val(name), path(reads)
    // A tuple containing the name of the raw/trimmed read files and the contigs belonging to the endosymbiont genome
    tuple val(name), path(assembled_endosymbiont)

    output:
    // The log file created by bowtie2 is output
    path 'log.txt', emit: alignment_stats

    script:
    /*
    Script description:
    1. Bowtie2 database is built from endosymbiont genome
    2. bowtie mapping with half the input threads (se bowtie2 documentation for details)
    3. command log is concatenated into log.txt
    */
    """
    bowtie2-build ${assembled_endosymbiont} endosymbiont
    bowtie2 -x endosymbiont -p ${task.cpus} -1 ${reads[0]} -2 ${reads[1]} -S ${name}_endosym.sam --very-sensitive
    cat .command.log > log.txt
    """

}

process COVERAGEESTIMATE {

    /* 
        Process Description:
        Estimation of sequencing coverage using alignment stats from bowtie2, length of the assembled endosymbiont and length of the input reads
    */

    // Copies output file in output folder
    publishDir "${params.output}/$params.job_name", mode: 'copy'

    // Process label
    label 'normal'

    // If this process fails, it does not end the whole pipeline
    errorStrategy 'finish'

    // Name of read files
    tag "$params.job_name"

    input:
    // log.txt file from previous process
    path stats
    // A tuple containing the name of the raw/trimmed read files and the files themself
    tuple val(name_dump), path(reads)
    // A .fa file containing the assembled endosymbiont genome
    tuple val(name), path(assembled_endosymbiont)

    output:
    // The process creates a file containign the coverage which is output
    path 'coverage.txt'

    script:
    /*
    Script description:
    1. Determining the number of bases in the endosymbiont genome
    2. Concatenating the alignment stats into a file
    3. Determining the number of bases in the read files
    4. Coverage estimate using a selfmade Python script
    */
    """
    grep -v ">" $assembled_endosymbiont | tr -d "\n" | wc -c > host_count.txt
    cat $stats > alignment_rate.txt
    cat ${reads[0]} | paste - - - - | cut -f 2 | tr -d '\n' | wc -c > base_count.txt
    python3 $project_dir/bin/coverage_estimate.py
    """
}

process CHECKENDOSYMBIONT {

    /* 
        Process Description:
        Evaluation of endosymbiont genome quality via CheckM
    */

    // Copies output file in output folder
    publishDir "${params.output}/$params.job_name/endosymbiont_assembly", mode: 'copy'

    // Process label
    label 'big_mem'

    // If this process fails, it does not end the whole pipeline
    errorStrategy 'finish'

    // Name of read files
    tag "$params.job_name"

    input:
    // A .fa file containing the endosymbiont genome
    tuple val(name), path(endosymbiont_genome)

    output:
    path '*'

    script:
    """
    checkm lineage_wf . checkm_output -t ${task.cpus} -x fa
    """

}

workflow {

    RAWQC(rawReads)
    TRIMMING(rawReads)
    TRIMMEDQC(TRIMMING.out)
    if (params.skip_trimming) {
        DENOVOASSEMBLY(rawReads, params.kmers) }
    else {
        DENOVOASSEMBLY(TRIMMING.out, params.kmers) }
    ENDOSYMBIONTCONTIGFILTERING(DENOVOASSEMBLY.out)
    endosymbiont_genome = ENDOSYMBIONTCONTIGFILTERING.out
    ENDOSYMBIONTGENOMEQUALITY(endosymbiont_genome)
    CHECKENDOSYMBIONT(endosymbiont_genome)
    if (params.skip_trimming) { 
        READMAPPINGFORCOVERAGE(rawReads, ENDOSYMBIONTCONTIGFILTERING.out)
        COVERAGEESTIMATE(READMAPPINGFORCOVERAGE.out, rawReads, ENDOSYMBIONTCONTIGFILTERING.out)
    }
    else {
        READMAPPINGFORCOVERAGE(TRIMMING.out, ENDOSYMBIONTCONTIGFILTERING.out)
        COVERAGEESTIMATE(READMAPPINGFORCOVERAGE.out, TRIMMING.out, ENDOSYMBIONTCONTIGFILTERING.out)
    }
    HOSTMITOGENOMEFILTERING(DENOVOASSEMBLY.out)
    HOSTMITOGENOMEQUALITY(HOSTMITOGENOMEFILTERING.out.host_filtered)
}

workflow.onComplete {
    // Display complete message
    log.info "Completed at: " + workflow.complete
    log.info "Duration    : " + workflow.duration
    log.info "Success     : " + workflow.success
    log.info "Exit status : " + workflow.exitStatus
}

workflow.onError {
    // Display error message
    log.info "Workflow execution stopped with the following message:"
    log.info "  " + workflow.errorMessage
}
