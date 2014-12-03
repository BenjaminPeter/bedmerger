import logging
import os

"""module to handle everything with saving and restoring SNP names
"""

def restore_snp_names( *args ):
    pass

def unify_snp_names( params ):
    for bim_file in params.bed:
        unify_bim_name( "%s.bim"%bim_file )

    for vcf_file in params.vcf:
        unify_vcf_name( vcf_file )

def unify_bim_name( bim_file ):
    logging.info("replacing snp ids in %s", bim_file )
    os.system("bash replace_bim_snp_id.sh %s"%bim_file)

def unify_vcf_name( vcf_file, twd ):
    """
    handles snp names.

    does the following:
    1. creates 
    """

    logging.info("replacing snp ids in %s", vcf_file )
    os.system("bash replace_vcf_snp_id.sh %s"%vcf_file)



############################################################
### functions not used in current version of program
############################################################
"""
these functions were used in the first, very basic version
and were replaced by the pandas pipeline
"""

def get_intersecting_snp( files, outfile="intersection.txt" ):
    """
    returns the SNP intersection from a bunch of bed files
    """
    files = [ add_extension( f, '.bim') for f in files ]
    fstreams = [ open( f, 'r' ) for f in files ]
    snp_sets = [ read_snp_id_from_bim( f ) for f in fstreams]

    snp_is = list( set.intersection( *snp_sets ) )

    snp_is.sort( key=sort_function )

    with open(params.out, "w") as f:
        for snp in snp_is:
            f.write( "%s\n"%snp )

def read_snp_id_from_bim( fstream ):
    s = set()
    for line in fstream:
        line = line.split()
        s.add( line[1] ) 


    fstream.close()
    return s

def sort_function( snp ):
    c, pos = snp.split("_")
    try:
        return int(c), int(pos)
    except ValueError:
        if c.upper() == "X":
            return 23, int(pos)
        if c.upper() == "Y":
            return 24, int(pos)
        if c.lower() == "mt":
            return 90, int(pos)

def add_extension( fname, extension=".bim"):
    base, ext = os.path.splitext( fname )
    if ext == extension:
        return fname

    return "%s.%s"%(fname, extension)
