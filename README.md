# RNAseq project, Osmotic tolerance in multiple species of killifish

[![Binder](http://mybinder.org/badge.svg)](http://mybinder.org:/repo/ljcohen/osmotic)

References
===========

Whitehead, A. (2010).  The evolutionary radiation of diverse osmotolerant physiologies in killifish (Fundulus sp.).  Evolution, 64(7): 2070-2085.  https://whiteheadresearch.files.wordpress.com/2012/05/whitehead-2010-evolution.pdf

Experimental Background
==========
- Goal: What are the gene expression differences between 3 treatments, by species?
- 3 Genera: Fundulus, Adinia, Lucania
- 12 species
- Sequencing done on 3 flowcells, 3 sample replicates per treatment
- 3 treatments/physiologies: freshwater (FW), brackish water (BW), transfer (brackish->fresh, sacrificed after 1 day)
- Data stored on farm.cse.ucdavis.edu here:

	/home/nreid/rnaseq/rawdata/

in three separate flowcells:

	141203_HS4B, 
	flowcell C5NU8ACXX

	141212_HS2A,
	flowcell C5VPYACXX

	141126_HS3A,
	flowcell C5NGRACXX
	
Basecall stats here:

	/home/nreid/rnaseq/rawdata/141212_HS2A/Basecall_Stats_C5VPYACXX.zip
	/home/nreid/rnaseq/rawdata/141203_HS4B/Basecall_Stats_C5NU8ACXX.tar
	/home/nreid/rnaseq/rawdata/141126_HS3A/??
	
Want to make a summary chart for each genus, species, treatment, sample:
	
	* flow cell
	* lane
	* read number
	* file numbers

Need to combine files from same lane (because default demultiplexing setting was used, threshold max # reads splitting into seperate files)
	
- info on using slurm: http://wiki.cse.ucdavis.edu/support:hpc:software:slurm
- files named, e.g.

	AWJRDD003_L_parva_transfer_3_TCTCGCGC-TATAGCCT_L006_R1_002.fastq.gz
	AWJRDD002_L_goodei_transfer_3_CGCTCATT-CCTATCCT_L004_R1_002.fastq.gz
	AWJRDD001_F_similis_BW_1_TCCGCGAA-ATAGAGGC_L002_R1_001.fastq.gz
	
	flowcell#_Genus_species_treatment_sample#_barcode_lane_read_file#.fastq.gz
	
- 2 populations of heteroclitus: MDPL, MDPL

	AWJRDD002_F_heteroclitus_MDPL_transfer_2_ATTCAGAA-ATAGAGGC_L004_R1_002.fastq.gz
	AWJRDD002_F_heteroclitus_MDPP_BW_2_ATTCAGAA-CCTATCCT_L003_R2_001.fastq.gz
	
	flowcell#_Genus_species_population_treatment_sample#_barcode_lane_read_file#.fastq.gz

- Subgoal: Since there is an annotated reference for Fundulus heteroclitus, it is not clear with a multispecies dataset like this whether a reference-guided or de novo transcriptome assembly approach would be better for analyzing differential expression of transcripts for this experiment. We will look at both and compare.

De novo transcriptome assembly
===================

Adapting the Eel pond khmer protocols by C. Titus Brown et al:

https://khmer-protocols.readthedocs.org/en/ctb/mrnaseq/

1. trim: https://khmer-protocols.readthedocs.org/en/ctb/mrnaseq/1-quality.html

	- Trimmomatic
	- combine lanes

2. diginorm: https://khmer-protocols.readthedocs.org/en/ctb/mrnaseq/2-diginorm.html
https://github.com/dib-lab/khmer-protocols/blob/jem-streaming/mrnaseq/1-quality.rst


3. Trinity: 
https://khmer-protocols.readthedocs.org/en/ctb/mrnaseq/3-big-assembly.html

	- by species
	- memory allocation

4. annotation/evaluation with transrate: http://hibberdlab.com/transrate/
	
	- use annotated translated aa from reference transcriptome?

5. differential expression

	- salmon
	- edgeR, DESEq2, limma
	
6. comparative transcriptomics

	
	
Reference-guided
===============

from Dr. Noah Reid, has already run these:

- Fundulus heteroclitus reference:
http://arthropods.eugenes.org/EvidentialGene/killifish/
- run .sh script in Noah's directory to expand and correct reference files: 
/home/nreid/popgen/kfish3/15.09.22.modify_gff.sh
- Send Don Gilbert email if questions
- gff reference version 20130322 we use slightly diffrent than public ncbi version
- avoid using ncbi set, a lot of effort already put in to use version 2rae5g before format changed for submitting to ncbi

1: run trimmomatic

	15.02.17.trimmomatic_array.sh 

2: run bwa mem, merge and index bams

	bwa_mem_paired_array.sh
	bwa_mem_unpaired_array.sh
	mergebams.sh
	indexbam.sh

3: run featurecounts

	15.11.13.featurecounts.sh

4: analyze counts in edgeR

	osmotic_edgeR_script.R

5: make some plots

	osmotic_edgeR_plots.R


