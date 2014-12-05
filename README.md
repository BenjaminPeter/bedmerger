run the following command for an up-to-date help message and an overview of all
commands:
    
    python bedmerger.py -h

requires plink 1.9, numpy and pandas, on spudhead/spudlings, run using python
2.7.6

#### known issue:
logging is very messy and not consistent

#### example usages:

- from the human origins data set, get all snp on chromosome 1 that have an ancestral allele typed 

    python bedmerger.py --bed /data/external_public/human_origins_affy/EuropeFullyPublic/vdata
			--check-reference 
			--reference_path ancestral_hg19 
			--pwd test 
			--keep_snp_id false 
			--out vdata_polarized
			--chromosome 1


    - merge the data sets from bhakar ,the human origins data set and 1000g, check the reference
	 and keep only snp present in all three data sets 

	python bedmerger.py --bed /data/external_public/human_origins_affy/EuropeFullyPublic/vdata

The help message output is pasted below for convenience:

    usage: bedmerger.py [-h] [--bed [BED [BED ...]]] [--vcf [VCF [VCF ...]]]
                        [--out OUT] [--output_type {bed,vcf}]
                        [--merge_type {inner,outer,left,right}]
                        [--retained_snp RETAINED_SNP]
                        [--keep_snp_id {left,false,merge}] [--pwd PWD] [--twd TWD]
                        [--check_reference] [--reference_path REFERENCE_PATH]
                        [--dropped_snp DROPPED_SNP] [--set_missing_to_reference]
                        [--chromosomes [CHROMOSOMES [CHROMOSOMES ...]]]
                        [--plink PLINK] [--logfile LOGFILE] [-d] [-v]

    A python wrapper around plink 1.9 for merging data sets in various file
    formats. Assumptions (currently NOT checked): - same positions implies same
    variant - indels are removed - alleles given on same strand in all files

    optional arguments:
      -h, --help            show this help message and exit
      --bed [BED [BED ...]], --bedfiles [BED [BED ...]]
                            The bedfiles to merge. The bim and fam files are
                            assumed to have the same name
      --vcf [VCF [VCF ...]], --vcffiles [VCF [VCF ...]]
                            The files in vcf format to merge
      --out OUT             the output file name prefix. Will have ending
                            bed/bim/fam (for bed)or vcf.gz (for vcf output)
      --output_type {bed,vcf}
                            the format of the output file
      --merge_type {inner,outer,left,right}
                            how the SNP are merged. Default is an outer join,
                            keeping all SNP. `inner` only keeps SNP that are
                            shared between populations. `left` and `right` keep
                            all SNP from the first/second data set, respectively
                            in each merge.
      --retained_snp RETAINED_SNP
                            file name of the file with all snp included in the
                            analyssis
      --keep_snp_id {left,false,merge}
                            if SNP id's should be retained. `false` will set it to
                            chr_pos_a1_a2, left keeps the id from the leftmost
                            file, as long as it is not na. `merge` concatenates
                            the SNP ids.
      --pwd PWD, --working-directory PWD
                            The working directory for relative paths
      --twd TWD, --temp-directory TWD
                            The directory where temporary files are stored.
                            created if it does not exist
      --check_reference, --check-reference
                            should SNP alleles be checked against reference
                            sequence?
      --reference_path REFERENCE_PATH, --reference-path REFERENCE_PATH
                            A folder with the reference sequence files in fasta
                            format, only required when --check_reference is used.
                            hg18, hg19, ancestral_hg19 and ancestral_hg38 lead to
                            reference sequences and ancestral sequences for these
                            genomes, respectively.
      --dropped_snp DROPPED_SNP
                            File with all SNP that were dropped from the analysis
      --set_missing_to_reference
                            should SNP absent from a data set be called as
                            reference?
      --chromosomes [CHROMOSOMES [CHROMOSOMES ...]], --chromosome [CHROMOSOMES [CHROMOSOMES ...]]
                            list of chromosomes to be included in data
      --plink PLINK         path to the plink executable to be used
      --logfile LOGFILE     file name for log file, default is stdout
      -d, --debug           Print lots of debugging statements
      -v, --verbose         Be verbose
