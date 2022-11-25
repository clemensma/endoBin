#!/usr/bin/env nextflow

import java.io.File
import org.apache.commons.lang.StringUtils

nextflow.enable.dsl=2

project_dir = projectDir

def helpMessage() {
    log.info"""
    Description:
        This pipeline separates endosymbiont prokaryotic genomes from their respective eukaryotic host mitogenome based on raw Illumina reads. It performs raw- and trimmed read quality control, adapter trimming, de-novo assembly, contig sorting and assembly quality assessment on paired short read sequence data.
    Pipeline summary:
         1. Raw reads quality control using FastQC
         2. Trimming using TrimGalore!
         3. trimmed read quality control using FastQC
         4. De Novo assembly using SPAdes
         5. Contig sorting based on endosymbiont reference genome using blastn and bfg
         6. Read mapping for coverage estimate using bowtie2
         7. Endosymbiont genome coverage estimate
         8. Endosymbiont genome quality assessment using BUSCO and CheckM
         9. Mitogenome extraction using blastn and bfg
        10. Mitogenome reassembly by NOVOPlasty
        11. Mitogenome strand control
        12. Mitogenome annotation with MITOS
        13. MITOS output formatting

    Usage:
        nextflow run main.nf -profile local
        
    Input:
        nextflow run main.nf --reads '<reads>.fasta' --endosymbiont_reference <reference_endosymbiont>.fasta --mitogenome_reference <reference_mito>.fasta
        [--help] [--output <dir>] [--contigs <contigs>.fasta] [-profile <name>] [--max_memory <value>.GB] [--max_cpus <value>]
        [--max_retries <value>] [--max_time '<value>.d'] [--trim_length <value>] [--trim_quality <value>] 
        [--trim_adapter <value>] [--trim_phred64] [--trim_clip_R1 <value>] [--trim_three_prime_clip_R1 <value>]
        [--trim_clip_R2 <value>] [--trim_three_prime_clip_R2 <value>] [--kmers [<value>]] [--meta]
        [--min_blast_wordsize <value>] [--max_blast_wordsize <value>] [--nucleotide_size <value>][--min_size <value> ]
        [--max_size <value>] [--kmer_size <value>] [--read_length <value>] [--insert_size <value>] 

    Mandatory arguments:

        --reads '<reads>.fasta'                       Path to one or more sets of paired-ended reads
                                                      (valid file types: .fastq.gz', '.fq.gz', '.fastq', or '.fq')
                                                      (default: $params.reads)

        --endosymbiont_reference <endosymbiont>.fasta Path to a reference genome for the endosymbiont
                                                      (valid file type extensions: '.fa', '.fna', '.fasta', '.faa')
                                                      (default: $params.endosymbiont_reference)

        --mitogenome_reference <relatedMito>.fasta    Path to a reference genome for the host species
                                                      (valid file type extensions: '.fa', '.fna', '.fasta', '.faa')  
                                                      (default: $params.mitogenome_reference)

    Optional arguments:

      //General
        --help                             Print help menu and exit
        --version                          Display the pipeline's version number and exit
        --output <dir>                     Specify output directory (default: $params.output)  
  
      //Flow control
        --mode <name>                      Select 'endo' or 'mito' to limit endosymbiont identification to the one specified (default: both)
        --contigs <contigs>.fasta          Provide previously assembled fasta contigs to skip assembly step
        -profile <name>                    Select 'local' for Docker or 'slurm' for Singularity
  
      //Resource allocation
        --max_memory <value>.GB            Maximum memory limit in GB                   (default: $params.max_memory)
        --max_cpus <value>                 Maximum number of threads                    (default: $params.max_cpus)
        --max_retries <value>              Maximum number of process retries            (default: $params.max_retries)
        --max_time '<value>.d'             Maximum runtime per process in days on a HPC (default: $params.max_time)
  
    Process specific arguments:

      //Trimming (Trim Galore!)
        --trim_length <value>              Remove reads shorter than specified             (default: $params.trim_length)
        --trim_quality <value>             Phred score threshold for quality trimming      (default: $params.trim_quality)
        --trim_adapter <value>             Specify adapter sequence to be trimmed          (default: auto-detect)
        --trim_phred64                     Switch to Phred+64 (=Illumina 1.5) encoding for quality scores (default: Phred+33; Sanger/Illumina 1.8)
        --trim_clip_R1 <value>             Cut off bases at the start of the forward reads (default: $params.trim_clip_R1)
        --trim_three_prime_clip_R1 <value> Cut off bases at the end of the forward reads   (default: $params.trim_three_prime_clip_R1)
        --trim_clip_R2 <value>             Cut off bases at the start of the reverse reads (default: $params.trim_clip_R2)
        --trim_three_prime_clip_R2 <value> Cut off bases at the end of the reverse reads   (default: $params.trim_three_prime_clip_R2)

      //Assembly (SPAdes)
        --kmers [<value>]                  A list of K-mer sizes to use for the assembly (default: $params.kmers)
                                           Examples: one K-mer: [81], multiple K-mers: [71, 81, 91]
        --meta                             Required for metagenomic sample data (default: true)
    
      //Mitogenome extraction (BLASTn)
        --min_blast_wordsize <value>       Minimum word size word size used for blastn search (default: $params.min_blast_wordsize)
        --max_blast_wordsize <value>       Maximum word size word size used for blastn search (default: $params.max_blast_wordsize)
        --nucleotide_size <value>          Estimated nucleotide count of target mitochondrion (default: $params.mito_size)
  
      //Mitogenome reassembly (NOVOPlasty)
        --min_size <value>                 Estimated minimum nucleotide count of target mitochondrion (default: $params.min_size)
        --max_size <value>                 Estimated maximum nucleotide count of target mitochondrion (default: $params.max_size)
        --kmer_size <value>                K-mer size used for reassembly (default: $params.kmer_size)
        --read_length <value>              Read length of Illumina short reads (default: $params.read_length)
        --insert_size <value>              Insert size of paired-end reads (default: $params.insert_size)
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

    // Copies output file in output folder
    publishDir "${params.output}/$params.job_name/host_assembly", mode: 'copy'

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
    if [[ "$params.contigs" != 'false' ]]
    then
      cat $contigs > contigs.fa
    else
      cat $contigs | better_fasta_grep "cov_([1-9][3-9][0-9]*|[1-9][0-9][0-9]{1,}|[2-9][0-9])\\.[0-9]+" > contigs.fa
    fi
      makeblastdb -in ${params.endosymbiont_reference} -title endosymbiont -parse_seqids -dbtype nucl -hash_index -out db
      blastn -query contigs.fa -db db -outfmt "10 qseqid" > seqid.txt
      cat contigs.fa | better_fasta_grep -F -f seqid.txt > endosymbiont_genome.fa
    """
}


process EXTRACTMITOGENOME {
    // Copies output file in output folder
    publishDir "${params.output}/$params.job_name/mitogenome/extraction", mode: 'copy'

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
    path('stats.txt') optional true

    script:
    """
    touch prev_seqid.txt
    touch unique_seqid.txt
    touch possible_mitogenomes.fa

    if [[ "$params.mito_min_size" = 'false' ]] || [[ "$params.mito_min_size" -gt "$params.mito_size" ]]
    then
      calc_threshold_085=\$(( "$params.mito_size*17/20" ))
      threshold_085=\$( echo \$calc_threshold_085 | awk '{printf("%d\\n",\$1 + 0.5)}' )
    else
      threshold_085="$params.mito_min_size"
    fi

    threshold_100=\$( echo "$params.mito_size" | awk '{printf("%d\\n",\$1 + 0.5)}' )
    calc_threshold_135=\$(( "$params.mito_size*27/20" ))
    threshold_135=\$( echo \$calc_threshold_135 | awk '{printf("%d\\n",\$1 + 0.5)}' )
    calc_threshold_200=\$(( "$params.mito_size*2" ))
    threshold_200=\$( echo \$calc_threshold_200 | awk '{printf("%d\\n",\$1 + 0.5)}' )

    if [[ "$params.contigs" = 'false' ]]
    then
      cat $contigs | better_fasta_grep "cov_[5-9][0-9]{1,}\\.[0-9]+" > cov_50_to_99.fa
      cat $contigs | better_fasta_grep "cov_[1-9][0-9][0-9]{1,}\\.[0-9]+" > cov_100_plus.fa
      cat cov_50_to_99.fa cov_100_plus.fa > cov_50_plus.fa
    fi
    cat $contigs > cov_0_plus.fa

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
        if [[ "$params.contigs" = 'false' ]]
        then
          cat cov_100_plus.fa | better_fasta_grep -f unique_seqid.txt > "blastn_covcut_100_wordsize_\$i.fa"
          cat cov_50_plus.fa | better_fasta_grep -f unique_seqid.txt > "blastn_covcut_50_wordsize_\$i.fa"
        fi
        cat cov_0_plus.fa | better_fasta_grep -f unique_seqid.txt > "blastn_covcut_0_wordsize_\$i.fa"
    done

    for file in blastn_*
    do
      if [[ \$(grep -c  '^>' \$file) -eq '1' ]] && [[ \$(grep -v  '^>' \$file | wc -m) -gt "\$threshold_085" ]]
      then
        cat \$file > mito_candidate_mitogenome.fa
        echo "Found the mitogenome on a single contig."
      fi
    done

    size_match () {
      echo "Starting search for closest mitogenome size match."
      for file in blastn_covcut_\${covcut}_*
      do
        grep -v  '^>' \$file | wc -m
      done > nucleotide_count_covcut_\${covcut}.txt
      closest_match=\$( awk -v c=1 -v t=\$threshold 'NR==1{d=\$c-t;d=d<0?-d:d;v=\$c;next}{m=\$c-t;m=m<0?-m:m}m<d{d=m;v=\$c}END{print v}' nucleotide_count_covcut_\${covcut}.txt )
      for blast_result in blastn_covcut_\${covcut}_*
      do
        if [[ \$(grep -v  '^>' \$blast_result | wc -m) = "\$closest_match" ]]
        then
          cat \$blast_result > mito_candidate_\${counter}_covcut_\${covcut}_size_match.fa
          break
        fi
      done
      echo "Finished search for closest mitogenome size match (cov \${covcut})."
    }
    
    if [[ ! -f mito_candidate_mitogenome.fa ]]
    then
      echo "Start size script"
      if [[ "$params.contigs" = 'false' ]]
      then
        covcut='100'; threshold=\$( echo "\$threshold_100" ); counter='1'; size_match
        covcut='100'; threshold=\$( echo "\$threshold_135" ); counter='7'; size_match
        covcut='50'; threshold=\$( echo "\$threshold_100" ); counter='3'; size_match
        covcut='50'; threshold=\$( echo "\$threshold_135" ); counter='8'; size_match
      fi
      covcut='0'; threshold=\$( echo "\$threshold_100" ); counter='5'; size_match
      covcut='0'; threshold=\$( echo "\$threshold_200" ); counter='9'; size_match
      echo "End size script"

      contig_match () {
      for blastn_result in blastn_covcut_\${covcut}_*
        do
                if [[ \$(cat \$blastn_result | wc -m) != '0' ]]
                then
                grep '^>' "\$blastn_result" > covcut_\${covcut}_header_list.txt
                while read -r header
                    do
                    better_fasta_grep "\$header" "\$blastn_result" | grep -v '^>' | wc -m
                done < covcut_\${covcut}_header_list.txt > "\${blastn_result%.fa}_covcut_\${covcut}_nuc_per_header.txt"
                awk 'BEGIN{s=0;}{s+=\$1;}END{print s/NR;}' "\${blastn_result%.fa}_covcut_\${covcut}_nuc_per_header.txt" > "\${blastn_result}_covcut_\${covcut}_avg_len.txt"
                fi
        done
        echo "Determined the average nucleotide size per contig for each blast result (cov \${covcut})."
        cat *_covcut_\${covcut}_avg_len.txt | sort -gr | head -1 | cut -d ' ' -f3 > covcut_\${covcut}_highest_avg.txt
        echo "Saved the highest average to the file covcut_\${covcut}_highest_avg.txt (cov 100)."
        for avg_len in *_covcut_\${covcut}_avg_len.txt
        do
          if [[ \$(cat "\$avg_len") = \$(cat covcut_\${covcut}_highest_avg.txt) ]]
          then
              novoplasty_seed="\${avg_len%_covcut_\${covcut}_avg_len.txt}"
              cat \$novoplasty_seed > mito_candidate_\${counter}_covcut_\${covcut}_contig_match.fa
          fi
        done

      }

        echo "Start contig script"
        if [[ "$params.contigs" = 'false' ]]
        then
          covcut='100'; threshold=\$( echo "\$threshold_100" ); counter='2'; contig_match
          covcut='50'; threshold=\$( echo "\$threshold_100" ); counter='4'; contig_match
        fi
        covcut='0'; threshold=\$( echo "\$threshold_100" ); counter='6'; contig_match
        echo "End contig script"

      seqkit stats *.fa > stats.txt

      rm blastn_*

      echo '0' > candidate_size_list.txt
      for candidate in mito_candidate_*
      do
        nucleotide_count=\$( grep -v '^>' \$candidate | wc -m)
        if grep -Fxq "\$nucleotide_count" candidate_size_list.txt
        then
          rm \$candidate
          echo "Removed candidate \$candidate. A file with \$nucleotide_count nucleotides is already included."
        else
          echo "\$nucleotide_count" >> candidate_size_list.txt
        fi
      done
    fi
    """
}

process REASSEMBLEMITOGENOME {
    // Copies output file in output folder
    publishDir "${params.output}/${params.job_name}/mitogenome/reassembly", mode: 'copy'

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
    separate_contigs () {
    COUNT="0"
    grep "^>" \$input | sort | uniq | while read -r header || [ -n "\$header" ]
    do
      COUNT=\$((\$COUNT + 1))
      echo \$header > contig_name_\${COUNT}.txt
    done
    COUNT="0"
    for header in contig_name_*.txt
    do
        COUNT=\$((\$COUNT + 1))
        search=\$( cat "\$header" )
        PRINT="0"
        while read line || [ -n "\$line" ]
        do
            if [[ "\$line" = "\$search" ]]
            then
                PRINT="1"
                echo \$line
                search='re_set_variable'
                continue
            fi
            if [[ \$PRINT = "1" ]] && [[ \${line:0:1} != ">" ]]
            then
                echo \$line
            else
                PRINT='0'
            fi
        done < \$input > \${step}_NOVOPlasty_contig_\${COUNT}.fa
    done
    rm contig_name_*.txt
    }
    
    select_largest_contig () {
    for contig in *_NOVOPlasty_contig_*.fa
    do
      grep -v "^>" \$contig | wc -m
    done > contig_sizes.txt
    largest_contig=\$( cat contig_sizes.txt | sort -gr | uniq | head -n 1 )
    rm contig_sizes.txt
    for contig in *_NOVOPlasty_contig_*.fa
    do
      if [[ \$(grep -v "^>" \$contig | wc -m) = "\$largest_contig" ]]
      then
        cat \$contig > largest_single_contig.fa
      fi
    done
    }
    
    create_stats () {
    for file in *.fa
    do
      contig_value=\$( grep '^>' \$file | wc -l )
      nucleotide_value=\$( grep -v '^>' \$file | wc -m )
      echo "\$contig_value \$nucleotide_value \$file" > \${file}_stats.txt
    done
    stats=\$( cat *_stats.txt )
    echo "Contigs | Nucleotides | Filename
    \$stats" > stats.txt
    rm *_stats.txt
    }

    if [[ -f mito_candidate_mitogenome.fa ]]
    then
      cat mito_candidate_mitogenome.fa > single_contig_mitogenome.fa
    elif [[ ! -f mito_candidate_mitogenome.fa ]]
    then
      if [[ "$params.mito_min_size" = 'false' ]] || [[ "$params.mito_min_size" -gt "$params.mito_size" ]]
      then
        calc_threshold_085=\$(( "$params.mito_size*17/20" ))
        threshold_085=\$( echo \$calc_threshold_085 | awk '{printf("%d\\n",\$1 + 0.5)}' )
      else
        threshold_085="$params.mito_min_size"
      fi

      calc_threshold_300=\$(( "$params.mito_size*3" ))
      threshold_300=\$( echo \$calc_threshold_300 | awk '{printf("%d\\n",\$1 + 0.5)}' )

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

      counter='0'
      candidate_list=($mitogenomes)    
      for i in "\${candidate_list[@]}"
      do
        counter=\$((\$counter + 1))
        if [[ \$(grep -c '^>' \$i) -eq '1' ]]
        then 
          mkdir -p NOVOPlasty_run_\$counter
          cat \$i > largest_single_contig.fa
          mv \$i largest_single_contig.fa NOVOPlasty_run_\$counter
        else
          if [[ \$(grep -v '^>' \$i | wc -m) -eq '0' ]] || [[ \$(grep -v '^>' \$i | wc -m) -gt "\$threshold_300" ]]
          then
            rm \$i
            continue
          fi
          cat \$i > split_mitogenome.fa
          perl /usr/bin/NOVOPlasty3.7.pl -c config.txt
          mkdir NOVOPlasty_run_\$counter
          mv contigs_tmp_Mitogenome.txt log_Mitogenome.txt NOVOPlasty_run_\$counter
          if [[ -f "Merged_contigs_Mitogenome.txt" ]]
          then
            mv Merged_contigs_Mitogenome.txt NOVOPlasty_run_\$counter
          fi

          if [[ -f "Circularized_assembly_1_Mitogenome.fasta" ]]
          then
            input="\$i"; step='pre'; separate_contigs
            input='Circularized_assembly_1_Mitogenome.fasta'; step='post'; separate_contigs
            select_largest_contig
            create_stats
            mv \$i *_NOVOPlasty_contig_*.fa largest_single_contig.fa Circularized_assembly_1_Mitogenome.fasta stats.txt NOVOPlasty_run_\$counter

          elif [[ -f "Uncircularized_assemblies_1_Mitogenome.fasta" ]]
          then
            input="\$i"; step='pre'; separate_contigs
            input='Uncircularized_assemblies_1_Mitogenome.fasta'; step='post'; separate_contigs
            select_largest_contig
            create_stats
            mv \$i *_NOVOPlasty_contig_*.fa largest_single_contig.fa Uncircularized_assemblies_1_Mitogenome.fasta stats.txt NOVOPlasty_run_\$counter

          elif [[ -f "Contigs_1_Mitogenome.fasta" ]]
          then
              echo "Mitogenome was not circularized."
              input="\$i"; step='pre'; separate_contigs
              input='Contigs_1_Mitogenome.fasta'; step='post'; separate_contigs
              select_largest_contig
              create_stats
              mv \$i *_NOVOPlasty_contig_*.fa largest_single_contig.fa Contigs_1_Mitogenome.fasta stats.txt NOVOPlasty_run_\$counter

          fi
          if [[ -f NOVOPlasty_run_\${counter}/largest_single_contig.fa ]] && [[ \$(grep -v '^>' NOVOPlasty_run_\${counter}/largest_single_contig.fa | wc -m) -gt "\$threshold_085" ]]
          then
            cat NOVOPlasty_run_\${counter}/largest_single_contig.fa > single_contig_mitogenome.fa
            cp single_contig_mitogenome.fa NOVOPlasty_run_\${counter}
            break
          fi
        fi
      done

      check_contigs=( NOVOPlasty_run_*/largest_single_contig.fa )
      if [[ ! -f single_contig_mitogenome.fa ]] && [[ -f "\$check_contigs" ]]
      then
        touch largest_contigs_list.txt
        for dir in NOVOPlasty_run_*/
        do
          if [[ -f "\${dir}"largest_single_contig.fa ]]
          then
            grep -v "^>" \${dir}largest_single_contig.fa | wc -m
          fi
        done > contig_sizes.txt
        largest_contig=\$( cat contig_sizes.txt | sort -gr | uniq | head -n 1 )
        rm contig_sizes.txt
        for dir in NOVOPlasty_run_*/
        do
          if [[ \$(grep -v "^>" "\${dir}"largest_single_contig.fa | wc -m) = "\$largest_contig" ]]
          then
            cat "\${dir}"largest_single_contig.fa > single_contig_mitogenome.fa
          fi
        done
      fi

    fi

    create_stats
    """

}

process STRANDCONTROL {
    // Copies output file in output folder
    publishDir "${params.output}/$params.job_name/mitogenome", mode: 'copy'

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
    path('mitogenome.fa'), emit: strand_tested_mitogenome

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
      grep '^>' -v single_contig_mitogenome.fa | tr -d '\n' | rev | tr "ATGC" "TACG" > mitogenome_seq.fa
      head -n 1 single_contig_mitogenome.fa > mitogenome_header.fa
      cat mitogenome_header.fa mitogenome_seq.fa > mitogenome.fa
      rm mitogenome_header.fa mitogenome_seq.fa
    else
      cat single_contig_mitogenome.fa > mitogenome.fa
    fi
    """
}

process ANNOTATEMITOGENOME {
    // Copies output file in output folder
    publishDir "${params.output}/$params.job_name/mitogenome/MITOS_annotation", mode: 'copy'

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
    mv mitos_output.txt mitos_output
    """
}

process MITOSFORMATTING {
    // Copies output file in output folder
    publishDir "${params.output}/$params.job_name/mitogenome/MITOS_annotation", mode: 'copy'

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
    while read -r line; do gene=\$( echo "\$line" );  better_fasta_grep "\$gene" individual_genes_nuc/result.fas > individual_genes_nuc/\$gene.fna; done < individual_genes_nuc.txt
    cat individual_genes_prot/result.faa | grep '^>' | sed 's/^.*@//' > individual_genes_prot.txt
    while read -r line; do gene=\$( echo "\$line" );  better_fasta_grep "\$gene" individual_genes_prot/result.faa > individual_genes_prot/\$gene.faa; done < individual_genes_prot.txt
    if [[ -f individual_genes_nuc/nad4.fna ]]
    then
    better_fasta_grep -v -F nad4l individual_genes_nuc/nad4.fna > nad4.fna
    mv nad4.fna individual_genes_nuc/nad4.fna
    fi
    if [[ -f individual_genes_prot/nad4.faa ]]
    then
    better_fasta_grep -v -F nad4l individual_genes_prot/nad4.faa > nad4.faa
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
    busco -i $endosym -m genome -o qc -l bacteria_odb10 --offline --download_path $project_dir/seqs/busco_data
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
