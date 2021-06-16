# Workflow for the assembly of one or more metagenomic samples

This workflow will remove any host or human genome circular consensus sequences (ccs) from the read files provided, and then will assemble each read set into a draft genome assembly using [Canu](https://canu.readthedocs.io/en/latest/). At the end, it will run [Metaquast](http://quast.sourceforge.net/metaquast.html) to assess the metagenome's quality parameters and composition. It will also produce a plot with the contig lengths for each assembly.

### generating circular consensus sequences

Most sequencing facilities will provide pacbio reads in the form of raw subreads. These shall not be used as such, as their error rate is usually pretty high (10-15%). Instead, they shall be processed to obtain circular consensus sequences (ccs), which are the standard pacbio read data. Ccs are obtained by overlapping all the subreads from the same DNA stretch, to obtain a consensus that reduces the error rate from 10-15% to almost 1%.

##### Using Lima

The first step, if not done by your sequencing core facility, is to demultiplex the subreads into individual files. This is an optional step and depends on your experiment layout. To do so, I used [lima](https://github.com/PacificBiosciences/barcoding) which is the official pacbio tool. Lima wants subreads in **bam** format and a file containing barcodes in **fasta** format that are necessary information to demultiplex the reads. Note: lima's documentation has moved [here](https://lima.how/).

```
lima \
--split-named \
--ccs \
--num-threads 16 \
--min-score 45 \
/path/to/subreads.bam \
/path/to/barcodes.fasta \
/path/to/output.bam
```

This command takes the subreads in `subreads.bam`, it takes the barcodes in `barcodes.fasta`, and uses this information to demultiplex the reads and place them in different files depending on the sample name (`--split-named`). The minimum score (`--min-score 45`) is the one suggested from PacBio guidelines for PacBio reads.

##### Using ccs

```
ccs \
--min-passes 3 \
--top-passes 60 \
--num-threads 16 \
--report-file /path/to/output.report \
--log-file /path/to/output.log \
/path/to/subreads.bam \
/path/to/ccs.bam
```

```
bam2fasta \
--output /path/to/output/folder/and_prefix \
-u \
${PREFIX}.bam
```


### Install dependencies

This pipeline depends on the following programs:

| Software    | Version | Type        | Link                                                          |
|-------------|---------|-------------|---------------------------------------------------------------|
| Minimap2    | 2.18    | Program     | https://github.com/lh3/minimap2                               |
| Samtools    | 1.12    | Program     | http://www.htslib.org/download/                               |
| Canu        | 2.1.1   | Program     | https://canu.readthedocs.io/en/latest/                        |
| Metaquast   | 5.0.2   | Program     | http://quast.sourceforge.net/docs/manual.html                 |
| Bioawk      | 1.0     | Program     | https://anaconda.org/bioconda/bioawk                          |
| Prokka      | 1.14.5  | Program     | https://github.com/tseemann/prokka                            |
| Roary       | 3.13.0  | Program     | https://sanger-pathogens.github.io/Roary/                     |
| Mcl         | 14-137  | Program     | https://micans.org/mcl/                                       |
| Mcxdeblast  | 12-068  | Program     | https://micans.org/mcl/man/mcxdeblast.html                    |
| Blastp      | 2.11.0  | Program     | https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ |
| Makeblastdb | 2.11.0  | Program     | https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ |
| Rscript     | 3.6.2   | Interpreter | https://cran.r-project.org/doc/manuals/r-release/R-admin.html |
| Python      | 3.8     | Interpreter | https://www.python.org/downloads/                             |

Make sure you have them installed before using this pipeline. The installation path of each executable must be modified accordingly in the `nextflow.config` file contained in this repository. Main executables don't need to be in the environment `$PATH`, but you must be sure that they work before using this pipeline. Most of them have dependencies of their own, so take your time to install each of them in your system (or have your sysadmin do it, probably a better choice).

### Run the pipeline

To see the pipeline help section, run:

```
nextflow run main.nf --help
```

Then, run the pipeline with a command like this one, which will also generate a report, a timeline and a direct acyclic graph for the proceeding of the pipeline itself:

```
nextflow \
run \
main.nf \
--output_dir /path/to/output \
--threads 48 \
--fasta_ccs_dir ${WD}/raw_data/demux_ccs_fasta \
--output_dir ${WD} \
--host_genome_seq ${WD}/raw_data/host_genomes/Gallus_gallus.GRCg6a.dna.toplevel.fa \
--human_genome_seq ${WD}/raw_data/host_genomes/hg38.fa \
--max_mem 90 \

```

##### Nextflow tips

- You can resume a crashed run by re-running the same command and adding `-resume` as an option after `run`. More info [here](https://www.nextflow.io/docs/latest/getstarted.html).
- You can specify where to save the temporary files of the pipeline by specifying a `-work-dir` directory right after `run`. More info on that and on other options available in `nextflow run` can be found [here](https://www.nextflow.io/docs/latest/cli.html#clean).
- The `work` directory tends to become quite crowded and full thousands of internal nextflow files, so every now and then make sure you clean it with `nextflow clean`. More info [here](https://www.nextflow.io/docs/latest/cli.html#clean)

##### Editing the pipeline

If you don't like the plots, if you want to change something in the code to accustom it to your own taste, you can edit directly the `Rscript` or `py` files inside the `/src` directory of the git repository that you cloned.


### Pipeline steps and output

Below are the steps performed by the pipeline, briefly described, and their output.

##### Removing reads from unwanted sources

Reads from the host genome and the human genome are removed from each ccs pool. These two genomes are passed as arguments (see `--help`). Reads that do not align to these two genomes are retained and used for the assembly. These operations are performed with **minimap2** and **samtools**. The output of this step is placed in the `/filt_reads` directory inside the specified `--output_dir`.

##### Metagenomic assembly

Metagenomic reads are assembled into a metagenomic assembly using the celera assembler **canu**. The pipeline is currently optimized for **pacbio hifi** reads only, but upon request I can modify it and make it able to use also earlier pacbio reads, illumina reads, and nanopore reads. It's a matter of time, and I didn't have any when I made it. The canu assembler has many options, which can be controlled from the nextflow command line parameters and from the `nextflow.config` file. Make sure you set all the error rates the way you want, as they can make or break an assembly in terms of quality. Hint: little is known about what paramters to use for metagenomic samples. Test them out, see what works for you.

The results of the assembly are assessed with **metaquast**, which is a metagenomic version of the popular assembly quality control tool quast. The contig lengths obtained from each sample are then extracted and a plot is generated for each of them.

The output of this step is placed in the `/assembly` directory inside the specified `--output_dir`. The `canu` subdirectory contains the assembly results and the contig lengths plots. The `metaquast` subdirectory contains the output of metaquast, which is famously a lot of files. You'll have to look through them yourself, but you may be interested in opening the `report.html` file with your internet browser.

##### Detection of genes

In the final step of the pipeline, the prokariotic gene predictor **prokka** is used to extract all candidate coding sequences from the assembled contigs. This tool requires a lot of `$PATH` dependencies, so make sure to have it installed properly (see chapter above). Then, the core genome calculator **roary** is used to obtain the number and identity of genes that can make up the core genome of the provided samples. It's important that you set the `--perc_core_isolates` option properly. This option controls the % of samples that have to feature a certain gene, for that gene to be part of the core genome. By default, the value is 99. However, if you have only 5 samples, there is no way you can have a gene in 99% of the isolates: it's either 80% or 100%. Make your own choice but keep this thing in mind.

The results of this step are placed inside the `/pangenome` directory inside the specified `--output_dir`. The `prokka` subdirectory will contain one folder per sample, named as the sample. Inside it, there will be all the prokka output files, which are a lot in many formats. Look up in their official documentation and what is what, but you may be interested into the `*.gff` file (containing the gene annotations), the `*.fna` file (containing gene sequences) and the `*.faa` file (containing protein sequences). The `*.txt` file contains a few numbers describing how many genes were predicted.

The **gff** files produced by prokka are then used to calculate the core genome with **roary**. The results are placed inside `/pangenome/roary`. 
