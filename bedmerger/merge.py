import socket
import logging

import vcf
import utils
import refseq
import grid

"""
This module handles the main merging pipeline
"""


def merge(params):
    """
    merges all input files into a single file

    merges all input files. The basic approach here is to
    1) construct reference seq table
    2) transform vcf file to plink format
    3) rename all snp names in all plink files
    4) filter plink files to only keep desired snp
    5) merge the files using plink
    6) restore snpids

    Parameters
    ----------

    params : Parameters
        a Parameters object with all the options, see the -h output of the
        argparser
    """

    if socket.gethostname() == "spudhead":
        return merge_grid(params)

    # 1
    construct_reference_sequence(params)
    
    # 2
    transformed_vcf_bed = transform_vcf_to_bed(params.vcf,
                                               params.twd,
                                               params.plink)

    # 3
    tmp_bim_files = refseq.rename_snpids_from_data(transformed_vcf_bed)
    tmp_bim_files += refseq.rename_snpids_from_data(params.bed,
                                                    params.twd)

    # 4
    filtered_beds = filter_bed_files(bed_files=transformed_vcf_bed + params.bed,
                                     bim_files=tmp_bim_files,
                                     params=params)

    # 5
    merge_beds(filtered_beds, params)

    # 6
    refseq.restore_snpids(params.out,
                          utils.path(params.twd, "original_ids.txt"))


def construct_reference_sequence(params):
    """
    gets the reference allele for each SNP position

    queries the files in params.reference_path to get the
    reference allele for all sites, then writes plink files
    to be used subsequently

    Parameters
    ----------
    params : Parameters
        parameters object

    Returns
    -------
    ref : pd.DataFrame
        reference sequence
    """
    if params.check_reference:
        logging.info("starting reference construction")
        ref = refseq.make_reference_allele_table(
                vcf_files = params.vcf, 
                bim_files = params.bed,
                reference_path = params.reference_path,
                merge_type = params.merge_type,
                chromosomes=params.chromosomes
               )
        refseq.write_all_plink_files(ref, params)
        logging.info("finished reference construction")

        return ref


def transform_vcf_to_bed(vcf_files, twd, plink="plink"):
    """
    transforms the input vcf files to bed files for merging

    as plink only merges at most one vcf file, we convert
    the vcf files to bed first. The new files are written in
    tmp_bed_files
    
    Parameters
    ----------
    vcf_files : list of strings
        names of vcf files
    twd : string
        temporary working directory
    plink : string
        name of plink exe
    
    Returns
    -------

    tmp_bed_files : list of strings
        bed file id of new bed files

    """
    logging.info("starting vcf to bed transformation")
    tmp_bed_files = []
    for i, vcf_file in enumerate(vcf_files):
        tmpname = utils.path(twd, "tmp_vcf%s"%i)


        vcf.to_bed (vcf_file, bed_name = tmpname,
                    plink = plink)
        
        tmp_bed_files.append(tmpname)


    logging.info("finished vcf to bed transformation")
    return tmp_bed_files


def filter_bed_files(bed_files, bim_files, params):
    """ filters SNP in all bed files
    
    Parameters
    ----------
    bed_files : string
        identifier for the bed and fam files, without extension
    bim_files : string
        name of the bim file, with extension
    params : Parameters object
        Parameters object
    
    Returns
    -------

    new_bed_ids : list of strings
        the file ids for the filtered bed files
    """

    logging.info("starting bed filtering")

    new_bed_ids = []
    for i, bed_file in enumerate(bed_files):


        filtered_bed_name = utils.path(params.twd, 
                "filtered_bed%s"%i)

        filter_bed(bed_file, bim_files[i],
                filtered_bed_name, 
                params.retained_snp,
                params.plink)
        new_bed_ids.append(filtered_bed_name)




    logging.info("finished bed filtering")
    return new_bed_ids


def filter_single_bed(bed_file, bim_file, new_name, extract,
                      plink="plink"):
    """filters a single bed file by removing the snp in exclude
    
    Parameters
    ----------
    bed_file : string
        the name of the bed and fam files, without extension
    bim_file : string
        name of the bim file, with extension
    new_name : string
        output name for the bim, bam and fam files to be created
    extract : string
        name of file with snp to be extracted/kept
    plink : string
        name of plink exe
    
    """

    flags = dict()
    flags['bed'] = bed_file + ".bed"
    flags['bim'] = bim_file
    flags['fam'] = bed_file + ".fam"
    flags['out'] = new_name
    flags['extract'] = extract
    flags['make-bed'] = ''

    utils.run_plink(plink, flags)


def merge_beds(bed_files, params, merge_file="merge.txt"):
    """
    merges all the bed files
    
    sets up a plink command that merges all input bed files.
    in order to do so, we have to create a merge file containing
    all the bed files
    
    Parameters
    ----------
    bed_files : string
        all the bed files to merge
    params : Parameters
        parameters object
    merge_file : string
        name of the plink merge file to be created in the temporary
        working directory
    
    Returns
    -------
    """

    logging.info("starting merging")

    final_filter = False
    if params.output_type is not 'bed':
        final_filter = True
    elif params.subset_snp is not None:
        final_filter = True
    elif params.subset_individuals is not None:
        final_filter = True
    elif params.set_missing_to_reference:
        final_filter = True

    merge_file = utils.path(params.twd, "merge.txt")
    with open(merge_file, 'w') as handle:
        for bed_file in bed_files[1:]:
            handle.write("%s\n" % bed_file)

    flags = dict()
    flags['merge-list'] = merge_file
    flags['bfile'] = bed_files[0]
    flags['merge-mode'] = 1

    if final_filter:
        flags['out'] = utils.path(params.twd, "tmp_out")
    else:
        flags['out'] = params.out

    if params.check_reference:
        flags['a2-allele'] = params.ref_allele

    utils.run_plink(params.plink, flags)

    if final_filter:
        logging.info(" applying final filter ")
        filter_flags = dict()
        filter_flags['bfile'] = flags['out']
        
        if params.output_type is 'bed':
            filter_flags['make-bed'] = ''
        
        if params.output_type is 'vcf':
            filter_flags['recode'] = 'vcf-iid'
            
        if params.output_type is 'ped':
            filter_flags['recode'] = ''

        if params.subset_snp is not None:
            filter_flags['extract'] = params.subset_snp

        if params.subset_individuals is not None:
            raise NotImplementedError()
            filter_flags['extract'] = params.subset_individuals

        filter_flags['out'] = params.out

        if params.check_reference:
            filter_flags['a2_allele'] = params.ref_allele

        if params.set_missing_to_reference:
            filter_flags['fill-missing-a2'] = ''
    
        utils.run_plink(params.plink, filter_flags)
    
    logging.info("finished merging")


def merge_grid(params):
    """merge_grid
    starts an SGE job and launches itself
    
    
    Parameters
    ----------
    params : type
        Desc
    
    """
    grid.merge(params)
