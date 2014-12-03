import numpy as np

import utils


def write_plink_reference_file(reference_alleles, file_name="dump.txt"):
    """writes file with reference information for plink.
    Run in plink with the --a2-allele flag
    
    Parameters
    ----------
    reference_alleles : refseq
        reference sequence data to be written
    file_name : string
        file to write
    
    """
    snps = reference_alleles[reference_alleles['consistency']]
    s = ""
    for line in (snps[['chrom', 'pos', 'ref']]).drop_duplicates().values:
        s += "%s_%s\t%s\n" % tuple(line)

    with open(file_name, 'w') as handle:
        handle.write(s)


def write_exclusion_file(reference_alleles, file_name="exclude.txt"):
    """writes file with alleles to exclude in plink
    in chrom_pos format
    
    Parameters
    ----------
    reference_alleles : refseq
        reference sequence data to be written
    file_name : string
        file to write
    
    """
    snps = reference_alleles[np.logical_not(reference_alleles['consistency'])]
    with open(file_name, 'w') as handle:
        for line in snps[['chrom', 'pos']].values:
            handle.write("%s_%s\n" % tuple(line))


def write_inclusion_file(reference_alleles, file_name="include.txt"):
    """writes file with snp to be included in the analysis
    Run in plink with the --extract flag
    
    Parameters
    ----------
    reference_alleles : refseq
        reference sequence data to be written
    file_name : string
        file to write
    
    """
    cols = ['dataid']
    reference_alleles.to_csv(file_name, sep="\t",
                             columns=cols,
                             header=False,
                             index=False,
                             na_rep='NA')


def write_original_id_file(reference_alleles, file_name):
    """writes file to restore original sample id
    
    Parameters
    ----------
    reference_alleles : refseq
        reference sequence data to be written
    file_name : string
        file to write
    
    """

    cols = ['dataid', 'snpid']
    reference_alleles.to_csv(file_name, sep="\t",
                             columns=cols,
                             header=False,
                             index=False,
                             na_rep='NA')


def write_all_plink_files(reference_alleles, params):
    """writes all the files to interact with plink
    
    Parameters
    ----------
    reference_alleles : refseq
        reference sequence data to be written
    params : Parameters
        Parameters object with file names
    
    """
    write_exclusion_file(reference_alleles, params.dropped_snp)
    write_inclusion_file(reference_alleles, params.retained_snp)
    write_plink_reference_file(reference_alleles,
                               "dump.txt")
    write_original_id_file(reference_alleles,
                           utils.path(params.twd,"original_ids.txt"))
