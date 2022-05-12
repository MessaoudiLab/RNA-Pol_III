# SINE-paper
Human-VZV RNASeq Analysis
This analysis includes (1) a Nextflow pipeline used for reference genome preparation, RNA-Seq data quality control and processing, RNA-Seq data alignment using STAR, and sorting, (2) a bash script used for calculating raw gene counts with FeatureCounts, (3) an Rscript used for performing the differential expression analysis with DESeq2.

Nextflow Installation
To install Nextflow, make sure Java 8 or later is installed. Then enter this command into your terminal:

curl -s https://get.nextflow.io | bash
Singularity Container
The singularity container, rnaSeq.sif, will be employed by each of the scripts. The only installation necessary is for Java and Nextflow.

Nextflow Usage
For more information on how to use the Nextflow pipeline run the following command in your terminal:

./nextflow run rnaSeq.nf --help
Run the Nextflow pipeline with the following mandatory parameters: - --inputFolder Path to folder containing fastq files - --genome Path to reference genome fasta file - --annot Path to reference genome annotation (.gtf) file - --genomeDir Path to folder containing indexed genome files (if genome indexing occurred previously) OR Path to folder where indexed genome files will be stored upon creation

./nextflow run rnaSeq.nf --inputFolder <pathToFastqFilesDir> --genome <pathToReferenceFasta>
--annot <pathToGTF> --genomeDir <pathToIndexedGenomeFiles> [options]
The following optional parameters can also be used: - --outputFolder Path to output directory for result directories (default: ./results) - --genomeName Name of reference genome used for labelling output files (default: GRCh38) - --SAindex Length(bases) of the SA pre-indexing string. Must be scaled down for smaller genomes (default: 14) - --gzip Boolean indicating whether the input fastq files are gzipped (default: true) - --ext Extension of fastq files (default: _paired.fastq.gz) - --clip_r1 Number of bp to remove from 5' end of read one (default: 13) - --clip_r2 Number of bp to remove from 5' end of read two (default: 13) - --three_prime_clip_R1 Number of bp to remove from 3' end of read one (default: 2) - --three_prime_clip_R2 Number of bp to remove from 3' end of read two (default: 2) - --quality Phred score cutoff for trimming (default: 30) - --length Minimum length cutoff for trimming (default: 50) - --index Indicates whether the user requires the genome to indexed (default: false)

If needed, you can run the the above bash command within an sbatch script in order to send the entire Nextflow pipeline to the cluster.

Each individual Nextflow process is managed as a separate job that is submitted to the cluster by using the sbatch command. Resource requests and job characteristics can be maniupulated in the "nextflow.config" file.

FeatureCounts Usage
The aligned and sorted bam output files from the Nextflow pipeline will be used to produce raw gene counts. This part of the analysis consists of an sbatch script, runFeatureCounts.sh, that is used to run a single bash script, featureCounts.sh, within the singularity container.

The following parameters are required in this order to run the script: - Path to where the aligned and sorted bam files are located - Path to the reference annotation GTF file - Desired output file name

sbatch runFeatureCounts.sh <pathToBamFiles> <pathToGTF> <outputFile>
Resource requests and job configurations can be manipulated at the top of the 'runFeatureCounts.sh' file.

DESeq2 Differential Expresion Analysis
This part of the analysis consists of an Rscript that can be sent to the cluster using the sbatch command with the runRscript.sh script. An R Markdown html file has also been provided that provides further details.

This analysis requires the following input: 1. Raw count data (where each row contains the GeneID and the corresponding counts for each sample), 'data/rawCounts.tsv' 2. Sample data (where each row contains the sample name and the corresponding sample characteristics--Condition, Transfection, Cell line, etc.), 'data/sampleData.csv' 3. Ensembl ID file used for labelling GeneIDs, 'data/gns.tsv'

sbatch runRscript_DESeq2.sh <rawCounts> <sampleData> <gnsData>
