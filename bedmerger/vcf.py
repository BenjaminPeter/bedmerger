import os
from os.path import exists
import utils


def is_vcf(vcf_file):
    """ check if file ending is corresponding to vcf """
    if exists( vcf_file ) and vcf_file.endswith("vcf"):
        return True
    if exists( vcf_file ) and vcf_file.endswith("vcf.gz"):
        return True

    return False

def is_zipped( vcf_file ):
    """check if file is bgzipped"""
    if exists( vcf_file ) and vcf_file.endswith("gz"):
        return True
    return False
    
def is_indexed( vcf_file ):
    """ checks if tabix index exists """
    if exists( vcf_file + ".tbi") and vcf_file.endswith("gz"):
        return True
    return False

def zip_and_index(vcf_file, tabix="tabix", bgzip="bgzip"):
    """ bgzips and indexes a vcf file """
    os.system("%s %s"%( bgzip, vcf_file ))
    vcf_file += ".gz"
    index( vcf_file, tabix )

def index( vcf_file, tabix="tabix"):
    """ creates vcf index using tabix """
    os.system("%s -p vcf %s"%(tabix, vcf_file ))


def to_bed( vcf_file, bed_name, exclude=None, plink="plink" ):
    """ transforms a vcf into a bed file """

    flags=dict()
    if exclude is not None:
        flags['exclude'] = exclude
    flags['vcf'] =  vcf_file
    flags['make-bed'] = ''
    flags['out'] = bed_name

    utils.run_plink( plink, flags )

