#!/usr/bin/env sh
#SBATCH -N 1
#SBATCH -J as-MetaPac
#SBATCH --mail-type FAIL,END
#SBATCH --mail-user matteo.schiavinato@boku.ac.at
#SBATCH --account p71579
#SBATCH --qos mem_0096
#SBATCH --partition mem_0096

cd /gpfs/data/fs71579/schmat90/CF/projects/meta/M00001_PacBio_chicken_gut_metagenomics/scripts/wf-assembly

nextflow run main.nf -resume -work-dir /gpfs/data/fs71579/schmat90/CF/projects/meta/M00001_PacBio_chicken_gut_metagenomics/work -with-report cmd.sbatch.report.html -with-timeline cmd.sbatch.timeline.html -with-dag cmd.sbatch.dag.png --fasta_ccs_dir /gpfs/data/fs71579/schmat90/CF/projects/meta/M00001_PacBio_chicken_gut_metagenomics/raw_data/demux_ccs_fasta --output_dir /gpfs/data/fs71579/schmat90/CF/projects/meta/M00001_PacBio_chicken_gut_metagenomics --host_genome_seq /gpfs/data/fs71579/schmat90/CF/projects/meta/M00001_PacBio_chicken_gut_metagenomics/raw_data/host_genomes/Gallus_gallus.GRCg6a.dna.toplevel.fa --human_genome_seq /gpfs/data/fs71579/schmat90/CF/projects/meta/M00001_PacBio_chicken_gut_metagenomics/raw_data/host_genomes/hg38.fa --threads 48 --max_mem 90 

