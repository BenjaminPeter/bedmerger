import os
import argparse
import logging

import utils


def reference_shortcuts(fname):
    """
    loads some shortcuts for reference_path
    """
    if fname == "hg18":
        return "/data/external_public/reference_genomes/"\
            "hg18/chr%s.fa.gz"
    elif fname == "hg19":
        return "/data/external_public/reference_genomes/"\
            "hg19/chr%s.fa.gz"
    elif fname == "ancestral_hg38":
        return "/data/external_public/ancestral_alleles/"\
            "hg19/homo_sapiens_ancestor_%s.fa.gz"
    elif fname == "ancestral_hg19":
        return "/data/external_public/ancestral_alleles/"\
            "hg19/homo_sapiens_ancestor_%s.fa.gz"
    return fname


class Parameters(object):
    """ class that handles all parameters

    the main idea is that this class unifies the different ways a project can
    be loaded. Current options are:
        - a dict like object, passed to init
        - from an input file
        - from the command line
        
    the main input format is from the command line
    """

    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)

        self.create_defaults()

    @staticmethod
    def create_defaults():
        """create_default values of dict"""
        if hasattr(Parameters, "defaults"):
            return
    
        Parameters.create_parser()
        Parameters.default = Parameters.args.parse_args(args=[])

    @staticmethod
    def create_parser():
        """
        generates ArgumentParser and reads options from CLI

        use -h flag for details
        """

        if hasattr(Parameters, "args"):
            return
        args = argparse.ArgumentParser(description="""
        A python wrapper around plink 1.9 for merging data sets in various file
        formats. Assumptions (currently NOT checked):
            - same positions implies same variant
            - indels are removed
            - alleles given on same strand in all files
        """)

        # --------------------------------------------------
        #  formatting & input files
        # --------------------------------------------------
        args.add_argument("--bed", '--bedfiles', nargs="*",
                          default=[],
                          action='append',
                          type=utils.standardize_bed_file_name,
                          help="""The bedfiles to merge. The bim and fam files
                              are assumed to have the same name""")

        args.add_argument("--vcf", '--vcffiles', nargs="*",
                          action='append',
                          default=[],
                          help="The files in vcf format to merge")

        args.add_argument("--out", default=None,
                          help="""the output file name prefix. Will have ending
                          bed/bim/fam (for
                          bed)or vcf.gz (for vcf output)"""
                          )

        args.add_argument("--output_type", default="bed",
                          choices=['bed', 'vcf'],
                          help="the format of the output file")

        args.add_argument("--merge_type", default="outer",
                          choices=['inner', 'outer', 'left', 'right'],
                          help="""how the SNP are merged. Default is an outer
                          join, keeping all SNP. `inner` only keeps SNP that
                          are shared between populations. `left` and `right`
                          keep all SNP from the first/second
                          data set, respectively in each merge.
                          """)

        args.add_argument("--retained_snp", default="retained.txt",
                          help=""" file name of the file with all snp included
                          in the analyssis
                          """)

        args.add_argument("--keep_snp_id", default='false',
                          choices=['left', 'false', 'merge'],
                          help="""if SNP id's should be retained. `false` will set
                          it to chr_pos_a1_a2, left keeps the id from the
                          leftmost file, as long as it is not na.
                          `merge` concatenates
                          the SNP ids. 
                          WARNING: everything except false is very
                          slow at the moment!""")

#        args.add_argument("--pwd", '--working-directory', default="./",
#                          help="""The working directory for relative paths
#                          """)

        args.add_argument("--twd", '--temp-directory', default="/tmp",
                          help="""The directory where temporary files are
                          stored.
                          created if it does not exist
                          """)

        # --------------------------------------------------
        # reference checking
        # --------------------------------------------------
        args.add_argument("--check_reference", "--check-reference",
                          default=False, action='store_true',
                          help="""should SNP alleles be checked against 
                          reference sequence?
                          """)

        args.add_argument("--reference_path", "--reference-path",
                          default='/data/external_public/reference_genomes/'
                          'hg19/chr%s.fa.gz',
                          help="""A folder with the reference sequence files in
                          fasta format, only required when --check_reference is
                          used. hg18, hg19, ancestral_hg19 and ancestral_hg38
                          lead to reference sequences and ancestral sequences
                          for these genomes, respectively.
                          """,
                          type=reference_shortcuts)

        args.add_argument("--dropped_snp", default="dropped.txt",
                          help="""File with all SNP that were dropped from the
                          analysis""")

        args.add_argument("--set_missing_to_reference", default=False,
                          action='store_true',
                          help="""should SNP absent from a data set be called
                          as reference?
                          """)

        #--------------------------------------------------
        # subsetting
        #--------------------------------------------------
        args.add_argument("--chromosomes", "--chromosome", nargs="*",
                          default=None,
                          help="""list of chromosomes to be included in data
                          """)

        #--------------------------------------------------
        # executable paths
        #--------------------------------------------------
        args.add_argument("--plink", default="plink",
                          help="""path to the plink executable to be used""")


        #--------------------------------------------------
        # logging info
        #--------------------------------------------------
        args.add_argument("--logfile", default=None,
                          help="""file name for log file, default is stdout""")

        args.add_argument('-d', '--debug',
                          help='Print lots of debugging statements',
                          action="store_const", dest="loglevel",
                          const=logging.DEBUG, default=logging.WARNING
                          )
        args.add_argument('-v', '--verbose', help='Be verbose',
                          action="store_const", dest="loglevel",
                          const=logging.INFO
                          )
    
        Parameters.args = args
    
        return args

    @staticmethod
    def from_command_line():
        """loads arguments from command line
        
        Returns
        -------
        p : Parameters
            the parameters object read from the command line
        """

        Parameters.create_defaults()
    
        args = Parameters.args
        params = args.parse_args()

        params.bed = [item for sublist in params.bed for item in sublist]
        params.vcf = [item for sublist in params.vcf for item in sublist]

        if params.out is None:
            vcf_names = [os.path.splitext(v)[0] for v in params.vcf]
            vcf_names = [v.split(os.path.sep)[-1] for v in vcf_names]
            bed_names = [b.split(os.path.sep)[-1] for b in params.bed]
            files = bed_names + vcf_names
            params.out = "merged_" + "+".join(files)

        params.ref_allele = "dump.txt"
            
        p = Parameters(**params.__dict__)
        return p
    
    def setup(params):
        """setup operations, currently just creates directories
        """
#        if not os.path.exists(params.pwd):
#                os.makedirs(params.pwd)
#        os.chdir(params.pwd)
        if not os.path.exists(params.twd):
                os.makedirs(params.twd)
        out_dir = os.path.dirname(params.out)
        if not os.path.exists(out_dir) and not out_dir == '':
                os.makedirs(out_dir)
    
    def sanity_checks(params):
        """perform some checks that the input arguments make sense"""

        for bedfile in params.bed:
            utils.assert_bed_file_exists(bedfile)
        for vcffile in params.vcf:
            utils.assert_vcf_file_exists(vcffile)
        
        params.n_files = len(params.bed)+ len(params.vcf)
        
        if params.n_files > 2 and \
            params.merge_type in ("left", "right"):
            s=("""merge mode left/right 
            doesn't make sense with
            more than two files""")
            logging.error(s)
            raise ValueError(s)
        
        if params.set_missing_to_reference and not \
                params.check_reference:
            error_msg = """cannot set missing to reference when
            check reference is not enabled"""
            logging.error(error_msg)
            raise ValueError(error_msg)


        if params.n_files < 1:
            error_msg = """no input files specified"""
            logging.error(error_msg)
            raise ValueError(error_msg)
