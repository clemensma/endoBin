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
    log.info"endoMiner ~ version $workflow.manifest.version"
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
ch_rawReads = Channel
    .fromFilePairs( params.reads, size: 2, type: 'file' )
    .filter { it =~/.*\.fastq\.gz|.*\.fq\.gz|.*\.fastq|.*\.fq/ }
    .ifEmpty { exit 1,
             "No FASTQ files found with pattern '${params.reads}'\n" +
             "Escape dots ('.') with a backslash character ('\\')\n" +
             "Try enclosing the path in single-quotes (')\n" +
             "Valid file types: '.fastq.gz', '.fq.gz', '.fastq', or '.fq'" }

if ( params.contigs ) {
  ch_contigs = Channel
      .fromPath (params.contigs)
}

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
    path('scaffolds.fasta')

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
    path contigs
    // A fasta file containing the endosymbiont reference genome
    // path endosymbiont_reference

    output:
    // The process outputs a tuple with the reads name and a .fa file containing all the contigs belonging to the endosymbiont genome
    path('endosymbiont_genome.fa'), emit: endosym_mapped

    script:
    /*
    Script description:
    1. A blast database from the reference genome is created
    2. A blastn search with the de novo assembled contigs is performed and found contig ids are saved in a .txt
    3. Based on the contig ids, contigs are grepped from the de novo assembled contigs using bfg
    */
    """
    if [[ "$params.contigs" !== 'false' ]]
    then
      cat $contigs > contigs.fa
    else
      cat $contigs | bfg "cov_([1-9][3-9][0-9]*|[1-9][0-9][0-9]{1,}|[2-9][0-9])\\.[0-9]+" > contigs.fa
    fi
      makeblastdb -in ${params.endosymbiont_reference} -title endosymbiont -parse_seqids -dbtype nucl -hash_index -out db
      blastn -query contigs.fa -db db -outfmt "10 qseqid" > seqid.txt
      cat contigs.fa | bfg -F -f seqid.txt > endosymbiont_genome.fa
    """
}


process EXTRACTMITOGENOME {
    // Copies output file in output folder
    publishDir "${params.output}/$params.job_name/mitogenome_extraction", mode: 'copy'

    // Process label
    label 'fast'

    // If this process fails, it does not end the whole pipeline
    errorStrategy 'finish'

    // Name of read files
    tag "$params.job_name"

    input:
    // Assembled contigs fasta file, reference mitogenome, forward and reverse read corresponding to contigs
    path contigs

    output:
    // Mitogenome (assembled if necessary), NOVOPlasty results, statistics
    path("mito_candidate_*"), emit: mitogenome_candidates

    script:
    """
    touch prev_seqid.txt
    touch unique_seqid.txt
    touch possible_mitogenomes.fa

    if [[ "$params.contigs" !== 'false' ]]
    then
      cat $contigs > cov_50_plus.fa
      cat $contigs > cov_100_plus.fa
    else
      cat $contigs | bfg "cov_[5-9][0-9]{1,}\\.[0-9]+" > cov_50_to_99.fa
      cat $contigs | bfg "cov_[1-9][0-9][0-9]{1,}\\.[0-9]+" > cov_100_plus.fa
      cat cov_50_to_99.fa cov_100_plus.fa > cov_50_plus.fa
    fi
    makeblastdb -in $contigs -title contig -parse_seqids -dbtype nucl -hash_index -out db
    echo "blastdb created"
    for i in {${params.min_blast_wordsize}..${params.max_blast_wordsize}..1}
      do
        echo "starting iteration with word size \$i"
        cat unique_seqid.txt > prev_seqid.txt
        blastn -query ${params.mitogenome_reference} -db db -outfmt "10 sseqid" -word_size \$i -num_threads ${task.cpus} > seqid.txt
        echo "blastn complete"
        cat -n seqid.txt | sort -uk2 | sort -nk1 | cut -f2- | cat > unique_seqid.txt
        echo "made seqids unique"
        cat cov_100_plus.fa | bfg -f unique_seqid.txt > "mg_candidate_covcut_100_wordsize_\$i.fa"
        cat cov_50_plus.fa | bfg -f unique_seqid.txt > "mg_candidate_covcut_50_wordsize_\$i.fa"        
        if [[ \$(grep -v '^>' mg_candidate_covcut_100_wordsize_\$i.fa | wc -m) ==  '0' ]] && [[ \$(grep -v '^>' mg_candidate_covcut_50_wordsize_\$i.fa | wc -m) ==  '0' ]]
        then
          break
        fi
    done

    for file in mg_candidate*
    do
      if [[ \$(grep -c  '^>' \$file) ==  '1' ]] && [[ \$(grep -v  '^>' \$file | wc -m) > '14000' ]]
      then
        cat \$file > mito_candidate_mitogenome.fa
      fi
    done
    if [[ ! -f mito_candidate_mitogenome.fa ]]
    then
      for file in mg_candidate_covcut_100_*
      do
        grep -v  '^>' \$file | wc -m
      done > nucleotide_count_covcut_100.txt
      closest_match=\$( awk -v c=1 -v t=$params.nucleotide_size 'NR==1{d=\$c-t;d=d<0?-d:d;v=\$c;next}{m=\$c-t;m=m<0?-m:m}m<d{d=m;v=\$c}END{print v}' nucleotide_count_covcut_100.txt )
      for blast_result in mg_candidate_covcut_100_*
      do
        if [[ \$(grep -v  '^>' \$blast_result | wc -m) = "\$closest_match" ]]
        then
          cat \$blast_result > mito_candidate_covcut_100_size_match.fa
          break
        fi
      done
      for file in mg_candidate_covcut_50_*
      do
        grep -v  '^>' \$file | wc -m
      done > nucleotide_count_covcut_50.txt
      closest_match=\$( awk -v c=1 -v t=$params.nucleotide_size 'NR==1{d=\$c-t;d=d<0?-d:d;v=\$c;next}{m=\$c-t;m=m<0?-m:m}m<d{d=m;v=\$c}END{print v}' nucleotide_count_covcut_50.txt )
      for blast_result in mg_candidate_covcut_50_*
      do
        if [[ \$(grep -v  '^>' \$blast_result | wc -m) = "\$closest_match" ]]
        then
          cat \$blast_result > mito_candidate_covcut_50_size_match.fa
          break
        fi
      done
    fi
    if [[ ! -f mito_candidate_mitogenome.fa ]]
    then
      for blastn_result in mg_candidate_covcut_100_*
      do
              grep '^>' "\$blastn_result" > covcut_100_header_list.txt
              while read -r header
                  do
                  bfg "\$header" "\$blastn_result" | grep -v '^>' | wc -m
              done < covcut_100_header_list.txt > "\${blastn_result%.fa}_covcut_100_nuc_per_header.txt"
              awk 'BEGIN{s=0;}{s+=\$1;}END{print s/NR;}' "\${blastn_result%.fa}_covcut_100_nuc_per_header.txt" > "\${blastn_result}_covcut_100_avg_len.txt"
      done
      cat *_covcut_100_avg_len.txt | sort -gr | head -1 | cut -d ' ' -f3 > covcut_100_highest_avg.txt
      for avg_len in *_covcut_100_avg_len.txt
      do
        if [[ \$(cat "\$avg_len") = \$(cat covcut_100_highest_avg.txt) ]]
        then
            novoplasty_seed="\${avg_len%_covcut_100_avg_len.txt}"
            cat \$novoplasty_seed > mito_candidate_covcut_100_contig_match.fa
        fi
      done
      for blastn_result in mg_candidate_covcut_50_*
      do
              grep '^>' "\$blastn_result" > covcut_50_header_list.txt
              while read -r header
                  do
                  bfg "\$header" "\$blastn_result" | grep -v '^>' | wc -m
              done < covcut_50_header_list.txt > "\${blastn_result%.fa}_covcut_50_nuc_per_header.txt"
              awk 'BEGIN{s=0;}{s+=\$1;}END{print s/NR;}' "\${blastn_result%.fa}_covcut_50_nuc_per_header.txt" > "\${blastn_result}_covcut_50_avg_len.txt"
      done
      cat *_covcut_50_avg_len.txt | sort -gr | head -1 | cut -d ' ' -f3 > covcut_50_highest_avg.txt
      for avg_len in *_covcut_50_avg_len.txt
      do
        if [[ \$(cat "\$avg_len") = \$(cat covcut_50_highest_avg.txt) ]]
        then
            novoplasty_seed="\${avg_len%_covcut_50_avg_len.txt}"
            cat \$novoplasty_seed > mito_candidate_covcut_50_contig_match.fa
        fi
      done
    fi
    if [[ \$(grep -v  '^>' mito_candidate_covcut_50_size_match.fa | wc -m) < '2000' ]] && [[ ! -f mito_candidate_mitogenome.fa ]]
    then
        cat $contigs | bfg "cov_[1-9][0-9]{1,}\\.[0-9]+" > cov_10_plus.fa
      for i in {${params.min_blast_wordsize}..${params.max_blast_wordsize}..1}
        do
          echo "starting iteration with word size \$i"
          cat unique_seqid.txt > prev_seqid.txt
          blastn -query ${params.mitogenome_reference} -db db -outfmt "10 sseqid" -word_size \$i -num_threads ${task.cpus} > seqid.txt
          echo "blastn complete"
          cat -n seqid.txt | sort -uk2 | sort -nk1 | cut -f2- | cat > unique_seqid.txt
          echo "made seqids unique"
          cat cov_10_plus.fa | bfg -f unique_seqid.txt > "mitogenome_candidates_wordsize_\$i.fa"
        done
      for file in *candidate*
      do
        grep -v  '^>' \$file | wc -m
      done > nucleotide_count.txt
      closest_match=\$( awk -v c=1 -v t=$params.nucleotide_size 'NR==1{d=\$c-t;d=d<0?-d:d;v=\$c;next}{m=\$c-t;m=m<0?-m:m}m<d{d=m;v=\$c}END{print v}' nucleotide_count.txt )
      for blast_result in *candidate*
      do
        if [[ \$(grep -v  '^>' \$blast_result | wc -m) = "\$closest_match" ]]
        then
          cat \$blast_result > identified_mitogenome.fa
          break
        fi
      done
    fi
    """
}

process REASSEMBLEMITOGENOME {
    // Copies output file in output folder
    publishDir "${params.output}/$params.job_name/mitogenome_extraction", mode: 'copy'

    // Process label
    label 'normal'

    // If this process fails, it does not end the whole pipeline
    errorStrategy 'finish'

    // Name of read files
    tag "$params.job_name"

    input:
    // Assembled contigs fasta file, reference mitogenome, forward and reverse read corresponding to contigs
    path contigs
    path mitogenomes
    tuple val(name), path(rawreads)

    output:
    path('single_contig_mitogenome.fa'), emit: mitogenome
    path("NOVOPlasty_run_*"), type: 'dir' optional true

    """
    if [[ -f mito_candidate_mitogenome.fa ]]
    then
      cat mito_candidate_mitogenome.fa > single_contig_mitogenome.fa
    elif [[ ! -f mito_candidate_mitogenome.fa ]]
    then
      echo "Project:
      -----------------------
      Project name          = Mitogenome
      Type                  = mito
      Genome Range          = $params.min_size-$params.max_size
      K-mer                 = $params.kmer_size
      Max memory            = ${task.memory.toGiga()}
      Extended log          = 0
      Save assembled reads  = no
      Seed Input            = split_mitogenome.fa
      Extend seed directly  = no
      Reference sequence    =
      Variance detection    =
      Dataset 1:
      -----------------------
      Read Length           = $params.read_length
      Insert size           = $params.insert_size
      Platform              = illumina
      Single/Paired         = PE
      Combined reads        =
      Forward reads         = ${rawreads[0]}
      Reverse reads         = ${rawreads[1]}
      Store Hash            =
      Optional:
      -----------------------
      Insert size auto      = yes
      Use Quality Scores    = no
      Output path           = " > config.txt
      candidate_list=($mitogenomes)
      for i in "\${candidate_list[@]}"
      do
        cat \$i > split_mitogenome.fa
        NOVOPlasty.pl -c config.txt
        mkdir NOVOPlasty_run_\$i
        mv contigs_tmp_Mitogenome.txt log_Mitogenome.txt NOVOPlasty_run_\$i
        if [[ -f "Merged_contigs_Mitogenome.txt" ]]
        then
          mv Merged_contigs_Mitogenome.txt NOVOPlasty_run_\$i
        fi
        if [[ -f "Circularized_assembly_1_Mitogenome.fasta" ]]
        then
          cat Circularized_assembly_1_Mitogenome.fasta > single_contig_mitogenome.fa
          mv Circularized_assembly_1_Mitogenome.fasta NOVOPlasty_run_\$i
        elif [[ -f "Uncircularized_assemblies_1_Mitogenome.fasta" ]]
        then
          cat Uncircularized_assemblies_1_Mitogenome.fasta > single_contig_mitogenome.fa
          mv Uncircularized_assemblies_1_Mitogenome.fasta NOVOPlasty_run_\$i
        elif [[ -f "Contigs_1_Mitogenome.fasta" ]]
        then
          contig=\$( head -n 1 Contigs_1_Mitogenome.fasta )
          bfg -F \$contig Contigs_1_Mitogenome.fasta > single_contig_mitogenome.fa
          mv Contigs_1_Mitogenome.fasta NOVOPlasty_run_\$i
        fi
        if [[ -f single_contig_mitogenome.fa ]] && [[ \$(grep -v  '^>' single_contig_mitogenome.fa | wc -m) > '14000' ]]
        then
          break
        fi
      done
    fi
      """
}

process STRANDCONTROL {
    // Copies output file in output folder
    publishDir "${params.output}/$params.job_name/strand_control", mode: 'copy'

    // Process label
    label 'fast'

    // If this process fails, it does not end the whole pipeline
    errorStrategy 'finish'

    // Name of read files
    tag "$params.job_name"

    input:
    // Fasta file of assembled genome
    path assembled_mitogenome

    output:
    // Mitochondrial genome
    path('single_contig_mitogenome.fa'), emit: strand_tested_mitogenome
    path('original_single_contig_mitogenome.fa') optional true
    path('blast_output.txt')
    path('blast_strands.txt')

    script:
    """
    if [[ \$( cat single_contig_mitogenome.fa | grep -v '^>' | grep -c -i -e [*] ) > '0' ]]
    then
      
      tr -d \\* < single_contig_mitogenome.fa > new_single_contig_mitogenome.fa
      cat new_single_contig_mitogenome.fa > single_contig_mitogenome.fa
      rm new_single_contig_mitogenome.fa
    fi
    makeblastdb -in ${params.mitogenome_reference} -dbtype nucl -out reference
    blastn -db reference -query $assembled_mitogenome -word_size 15 -out blast_output.txt
    cat blast_output.txt | grep 'Strand' > blast_strands.txt
    if [[ \$( head -n 1 blast_strands.txt ) == *'Strand=Plus/Minus'* ]]
    then
      cat single_contig_mitogenome.fa > original_single_contig_mitogenome.fa
      # cat original_single_contig_mitogenome.fa | while read L; do  echo \$L; read L; echo "\$L" | rev | tr "ATGC" "TACG" ; done > single_contig_mitogenome.fa

    fi
    """
}

process ANNOTATEMITOGENOME {
    // Copies output file in output folder
    publishDir "${params.output}/$params.job_name/MITOS_annotation", mode: 'copy'

    // Process label
    label 'fast'

    // If this process fails, it does not end the whole pipeline
    errorStrategy 'finish'

    // Name of read files
    tag "$params.job_name"

    input:
    // Fasta file of assembled genome
    path mitogenome

    output:
    // Mitochondrial genome
    path("mitos_output"), emit: mitos_out


    script:
    """
    mkdir -p mitos_output
    runmitos.py -i $mitogenome -o mitos_output -r $params.mitos_reference -R $baseDir -c $params.genetic_code > mitos_output.txt
    """
}

process MITOSFORMATTING {
    // Copies output file in output folder
    publishDir "${params.output}/$params.job_name/MITOS_annotation", mode: 'copy'

    // Process label
    label 'fast'

    // If this process fails, it does not end the whole pipeline
    errorStrategy 'finish'

    // Name of read files
    tag "$params.job_name"

    input:
    // Fasta file of assembled genome
    path mitos_out_dir
    tuple val(name), path(rawreads)

    output:
    // Mitochondrial genome
    path("*")

    """
    mkdir -p individual_genes_nuc
    mkdir -p individual_genes_prot
    if [[ "$params.species_id" ]]
    then
      id=\$( echo "$params.species_id" )
    else
      id=\$( echo "${rawreads[0].simpleName}" )
    fi
    sed "s/^.*\\(; \\)/>\${id}@/g" mitos_output/result.fas | sed 's/(.*//' > individual_genes_nuc/result.fas
    sed "s/^.*\\(; \\)/>\${id}@/g" mitos_output/result.faa | sed 's/(.*//' > individual_genes_prot/result.faa
    cat individual_genes_nuc/result.fas | grep '^>' | sed 's/^.*@//' > individual_genes_nuc.txt
    while read -r line; do gene=\$( echo "\$line" );  bfg "\$gene" individual_genes_nuc/result.fas > individual_genes_nuc/\$gene.fna; done < individual_genes_nuc.txt
    cat individual_genes_prot/result.faa | grep '^>' | sed 's/^.*@//' > individual_genes_prot.txt
    while read -r line; do gene=\$( echo "\$line" );  bfg "\$gene" individual_genes_prot/result.faa > individual_genes_prot/\$gene.faa; done < individual_genes_prot.txt
    if [[ -f individual_genes_nuc/nad4.fna ]]
    then
    bfg -v -F nad4l individual_genes_nuc/nad4.fna > nad4.fna
    mv nad4.fna individual_genes_nuc/nad4.fna
    fi
    if [[ -f individual_genes_prot/nad4.faa ]]
    then
    bfg -v -F nad4l individual_genes_prot/nad4.faa > nad4.faa
    mv nad4.faa individual_genes_prot/nad4.faa
    fi
    if grep -q '\\-' "mitos_output/result.geneorder"
    then
      cat mitos_output/result.geneorder > mitos_output/original_result.geneorder
      sed -i -e 's/-//g' mitos_output/result.geneorder
    fi
    if grep -q 'cox1' "mitos_output/result.geneorder"
    then
      grep -v '^>' mitos_output/result.geneorder > mitos_output/current_order.txt
      while read -r gene_order; do
          if [[ \$( cat mitos_output/current_order.txt | awk '{print \$1;}' ) == *"cox1"* ]]
          then
              cat mitos_output/current_order.txt > mitos_output/adjusted_result.geneorder
          else
              last_gene=\$( awk '{ print \$NF }' mitos_output/current_order.txt )
              new_order=\$( sed "s/\\<\$last_gene\\>//" mitos_output/current_order.txt )
              echo "\$last_gene \$new_order" > mitos_output/current_order.txt
          fi
      done < mitos_output/current_order.txt
    fi
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
    path endosym

    output:
    // All files created by BUSCO are output
    path 'qc/*'

    script:
    // BUSCO command (see BUSCO documentation for details)
    """
    busco -i $endosym -m genome -o qc -l rickettsiales_odb10 --offline --download_path $project_dir/seqs/busco_data
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
    path assembled_endosymbiont

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
    tuple val(name), path(reads)
    // A .fa file containing the assembled endosymbiont genome
    path assembled_endosymbiont 

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
    path endosymbiont_genome

    output:
    path '*'

    script:
    """
    checkm lineage_wf . checkm_output -t ${task.cpus} -x fa
    """

}

workflow {
    RAWQC(ch_rawReads)
    TRIMMING(ch_rawReads)
    TRIMMEDQC(TRIMMING.out)
    if (!params.contigs){
      DENOVOASSEMBLY(TRIMMING.out, params.kmers) 
      contigs = DENOVOASSEMBLY.out
      }
    else {
      contigs = ch_contigs
    }
    if (params.mode == "endo" || params.mode == "both" ){
    ENDOSYMBIONTCONTIGFILTERING(contigs)
    endosymbiont_genome = ENDOSYMBIONTCONTIGFILTERING.out
    ENDOSYMBIONTGENOMEQUALITY(endosymbiont_genome)
    CHECKENDOSYMBIONT(endosymbiont_genome)
    READMAPPINGFORCOVERAGE(TRIMMING.out, ENDOSYMBIONTCONTIGFILTERING.out)
    COVERAGEESTIMATE(READMAPPINGFORCOVERAGE.out, TRIMMING.out, ENDOSYMBIONTCONTIGFILTERING.out)
    }
    if (params.mode == "mito" || params.mode == "both" ){
    EXTRACTMITOGENOME(contigs)
    REASSEMBLEMITOGENOME(contigs, EXTRACTMITOGENOME.out.mitogenome_candidates, ch_rawReads)
    STRANDCONTROL(REASSEMBLEMITOGENOME.out.mitogenome)
    ANNOTATEMITOGENOME(STRANDCONTROL.out.strand_tested_mitogenome)
    MITOSFORMATTING(ANNOTATEMITOGENOME.out.mitos_out, TRIMMING.out)
    }
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
