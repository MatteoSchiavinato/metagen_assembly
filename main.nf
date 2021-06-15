#!/usr/bin/env nextflow

// help section
if (params.help) {

  println """

  ### misc ###

  --threads               Number of parallel threads

  ### input ###

  --fasta-ccs-dir         Path to a directory containing the CCS reads in FASTA format
                          The directory will be scanned for files ending in {fa,fasta}
                          Other files will not be considered
                          The sample ID will be inferred from the file name, so name it
                          wisely
                          (mandatory)

  ### remove host/human reads ###

  --host_genome_seq       FASTA genome sequence of the host genome
  --human_genome_seq      FASTA genome sequence of the human genome

  ### output ###

  --output_dir            Base directory where all the sub-directories with the outputs
                          will be placed.
                          Ideally, this is the directory containing /raw_data and /scripts

  """
  exit 0
}


// ccs reads channel

Channel
  .fromPath("${params.fasta_ccs_dir}/*.{fa,fasta}")
  .map{ it -> [it.baseName, "ccs", it] }
  .set{ Minimap }

// remove host reads

process remove_host_reads {

  executor = "local"
  cpus = params.threads

  publishDir "${params.output_dir}/filt_reads", mode: "copy"

  input:
    tuple val(sample_id), val(read_type), file(fasta) from Minimap

  output:
    tuple val(sample_id), val(read_type), file("${sample_id}.${read_type}.no_host.fasta") into NoHost
    file "${sample_id}.${read_type}.no_host.fasta"  into NoHost_count

  script:
    """
    ${MINIMAP2} \
    -H \
    -k ${params.minimap_kmer_size} \
    -d host_genome.mm2_index \
    -a \
    -t ${params.threads} \
    -x map-pb \
    ${params.host_genome_seq} \
    ${fasta} | \
    ${SAMTOOLS} fasta \
    -@ ${params.threads} \
    -f 0x4 -F 0x0100 \
    - \
    > ${sample_id}.${read_type}.no_host.fasta \

    """
}

// count removed reads

process count_reads_wo_host {

  executor = "local"
  cpus = 1

  publishDir "${params.output_dir}/statistics", mode: "copy"

  input:
    file all_fastas from NoHost_count.collect()

  output:
    path "no_host.read_counts.tsv"

  script:
    """
    for k in ${all_fastas}
    do
      NAME=\$(echo \${k} | cut -d "." -f 1) &&
      ${BIOAWK} \
      -c fastx \
      -v x=\${NAME} \
      '{totseq+=1; totlen+=length(\$seq); } END {print x"\\t"totseq"\\t"totlen"\\t"totlen/totseq}' \${k}
    done \
    > no_host.read_counts.tsv
    """
}


// remove human reads

process remove_human_reads {

  executor = "local"
  cpus = params.threads

  publishDir "${params.output_dir}/filt_reads", mode: "copy"

  input:
    tuple val(sample_id), val(read_type), file(fasta) from NoHost

  output:
    tuple val(sample_id), val(read_type), file("${sample_id}.${read_type}.no_host.no_human.fasta") into NoHuman_canu
    file "${sample_id}.${read_type}.no_host.no_human.fasta" into NoHuman_count

  script:
  """
  ${MINIMAP2} \
  -H \
  -k ${params.minimap_kmer_size} \
  -d human_genome.mm2_index \
  -a \
  -t ${params.threads} \
  -x map-pb \
  ${params.human_genome_seq} \
  ${fasta} | \
  ${SAMTOOLS} fasta \
  -@ ${params.threads} \
  -f 0x4 -F 0x0100 \
  - \
  > ${sample_id}.${read_type}.no_host.no_human.fasta

  """
}


// count reads after removal of human genome

process count_reads_wo_human {

  executor = "local"
  cpus = 1

  publishDir "${params.output_dir}/statistics", mode: "copy"

  input:
    file all_fastas from NoHuman_count.collect()

  output:
    path "no_host.no_human.read_counts.tsv"

  script:
    """
    for k in ${all_fastas}
    do
      NAME=\$(echo \${k} | cut -d "." -f 1) &&
      ${BIOAWK} \
      -c fastx \
      -v x=\${NAME} \
      '{totseq+=1; totlen+=length(\$seq); } END {print x"\\t"totseq"\\t"totlen"\\t"totlen/totseq}' \${k}
    done \
    > no_host.no_human.read_counts.tsv
    """
}


// run HiCanu

process canu_assembly {

  executor = 'slurm'
  module = ['gnuplot/5.2.8-gcc-9.1.0-2dqwfce']
  clusterOptions = "-N 1 --ntasks-per-node ${params.threads} --partition mem_0096 --qos mem_0096 --account p71579"

  publishDir "${params.output_dir}/assembly/canu", mode: "copy"

  input:
    tuple val(sample_id), val(read_type), file(fasta) from NoHuman_canu

  output:
    tuple val(sample_id), \
    val(read_type), \
    file("${sample_id}.${read_type}.contigs.fasta"), \
    file("${sample_id}.${read_type}.unassembled.fasta"), \
    file("${sample_id}.${read_type}.report") \
    into Canu_out
    file("${sample_id}.${read_type}.contigs.fasta") into Canu_quast
    tuple val(sample_id), val(read_type), file("${sample_id}.${read_type}.contigs.fasta") into Contigs_plot
    tuple val(sample_id), val(read_type), file("${sample_id}.${read_type}.contigs.fasta") into Contigs_prokka

  script:
    """
    ${CANU} \
    -d . \
    -p ${sample_id}.${read_type} \
    useGrid=false \
    genomeSize=${params.avg_genome_size} \
    maxThreads=${params.threads} \
    maxMemory=${params.max_mem} \
    merylThreads=${params.meryl_threads} \
    merylMemory=${params.meryl_memory} \
    corMinCoverage=${params.cor_min_coverage} \
    corOutCoverage=${params.cor_out_coverage} \
    minReadLength=${params.min_read_length} \
    minOverlapLength=${params.min_overlap_length} \
    corOvlErrorRate=${params.cor_ovl_error_rate} \
    obtOvlErrorRate=${params.obt_ovl_error_rate} \
    utgOvlErrorRate=${params.utg_ovl_error_rate} \
    corErrorRate=${params.cor_error_rate} \
    obtErrorRate=${params.obt_error_rate} \
    utgErrorRate=${params.utg_error_rate} \
    cnsErrorRate=${params.cns_error_rate} \
    MMapMerSize=${params.minimap_kmer_size} \
    ovlMerSize=${params.kmer_size} \
    mhapMerSize=${params.kmer_size} \
    maxInputCoverage=${params.max_input_coverage} \
    utgReAlign=true \
    MhapSensitivity=high \
    -pacbio-hifi ${fasta}
    """
}


// check assembly quality

process metaquast {

    executor = 'local'
    cpus = params.threads

    publishDir "${params.output_dir}/assembly", mode: "copy"

  input:
    file all_fastas from Canu_quast.collect()

  output:
    path 'metaquast', hidden:false, type:"dir" into Metaquast

  script:
  """
    ${METAQUAST} \
    --output-dir metaquast \
    --contig-thresholds ${params.contig_thresholds} \
    --min-alignment ${params.quast_min_alignment} \
    --min-contig ${params.min_read_length} \
    --threads ${params.threads} \
    --min-identity ${params.min_identity} \
    ${all_fastas}
  """
}


// calculate contig lengths
process get_ctg_lengths {

  executor = "local"
  cpus = 1

  publishDir "${params.output_dir}/assembly/canu", mode: "copy"

  input:
    tuple val(sample_id), val(read_type), file(contigs_fasta) from Contigs_plot

  output:
    tuple val(sample_id), val(read_type), file("${sample_id}.${read_type}.contig_lengths.tsv") into Contig_lengths

  script:
    """
    ${BIOAWK} \
    -c fastx \
    '{print \$name"\\t"length(\$seq)}' \
    ${contigs_fasta} \
    > ${sample_id}.${read_type}.contig_lengths.tsv \

    """
}


// plot contig length

process plot_ctg_lengths {

  executor = "local"
  cpus = params.threads

  publishDir "${params.output_dir}/assembly/canu", mode: "copy"

  input:
    tuple val(sample_id), val(read_type), file(contig_lengths) from Contig_lengths

  output:
    file "${sample_id}.${read_type}.contig_lengths.png" into Ctg_len_pngs
    file "${sample_id}.${read_type}.contig_lengths.svg" into Ctg_len_svgs

  script:
    """
    ${RSCRIPT} \
    ${params.source_dir}/plot-contig-lengths.Rscript \
    ${contig_lengths} \
    ${sample_id} \
    ${read_type} \
    100000 \

    """
}



// run prokka to annotate the whole contigs at once

process run_prokka {

  executor = "local"
  cpus = params.threads

  publishDir "${params.output_dir}/pangenome/prokka", mode: "copy"

  input:
    tuple val(sample_id), val(read_type), file(fasta) from Contigs_prokka

  output:
    tuple val(sample_id), path("${sample_id}", type:"dir") into Prokka_out

  script:
    """
    ${PROKKA} \
    --metagenome \
    --cpus ${params.threads} \
    --kingdom ${params.kingdom} \
    --outdir ${sample_id} \
    --prefix ${sample_id} \
    --addgenes \
    --addmrna \
    --gffver 3 \
    --evalue ${params.evalue} \
    --mincontiglen ${params.min_contig_len} \
    ${fasta} \

    """
}
