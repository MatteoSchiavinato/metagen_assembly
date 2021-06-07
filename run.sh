#!/usr/bin/env sh

PROJECT="M00001_PacBio_chicken_gut_metagenomics"
WD="/gpfs/data/fs71579/schmat90/CF/projects/meta/${PROJECT}"

echo """\
#!/usr/bin/env sh
#SBATCH -N 1
#SBATCH -J as-MetaPac
#SBATCH --mail-type FAIL,END
#SBATCH --mail-user matteo.schiavinato@boku.ac.at
#SBATCH --account p71579
#SBATCH --qos mem_0096
#SBATCH --partition mem_0096

cd ${WD}/scripts/wf-assembly

nextflow \
run \
main.nf \
-resume \
-work-dir ${WD}/work \
-with-report cmd.sbatch.report.html \
-with-timeline cmd.sbatch.timeline.html \
-with-dag cmd.sbatch.dag.png \
--fasta_ccs_dir ${WD}/raw_data/demux_ccs_fasta \
--output_dir ${WD} \
--host_genome_seq ${WD}/raw_data/host_genomes/Gallus_gallus.GRCg6a.dna.toplevel.fa \
--human_genome_seq ${WD}/raw_data/host_genomes/hg38.fa \
--threads 48 \
--max_mem 90 \

""" \
> ${WD}/scripts/wf-assembly/cmd.sbatch

sbatch \
--output ${WD}/scripts/wf-assembly/cmd.sbatch.out \
${WD}/scripts/wf-assembly/cmd.sbatch
