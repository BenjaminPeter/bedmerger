import os
import logging


def standardize_bed_file_name(bedfile):
    """ removes the bed extension from bed name if present"""
    fname, ext = os.path.splitext(bedfile)
    assert ext == ".bed" or ext == ""

    return fname


def assert_bed_file_exists(bedfile):
    pass


def assert_vcf_file_exists(bedfile):
    pass


def run_plink(plink, flags):
    """ constructs a plink command line with flags from dict"""
    
    s = "%s " % plink
    for item in flags.iteritems():
        s += "--%s %s " % item

    # s += " >>plink.log"

    os.system(s)
    logging.debug("ran plink with the following command\n%s" % s)


def path(*args):
    return os.path.sep.join(args)
