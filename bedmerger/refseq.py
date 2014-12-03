import gzip
import numpy as np
import pandas as pd
import logging


from datatypes import DataFile, BimFile
import utils


def make_reference_allele_table(vcf_files=[], bim_files=[],
                                reference_path='/data/reference/hg19/',
                                merge_type="outer", id_mode="pos",
                                chromosomes=None):
    """make_reference_allele_table

    function to create a file with all reference alleles

    takes all vcf files and bim files and creates a pandas table with 5 columns
    - chromosome
    - position along chromosome
    - snp_name (final)
    - allele 1
    - allele 2
    - reference allele
    
    Parameters
    ----------
    vcf_files : str
        file names of vcf files, may be gzipped
    bim_files : str
        file id of bim files, without extension
    reference_path : str
        path to where reference alleles are found
    merge_type : str
        how data sets are merged, either inner or outer
    id_mode : str
        mode how to unify snp ids
    chromosomes : list of str
        chromosomes to be retained
    
    Returns
    -------
    reference_alleles : pd.DataFrame
        refseq table generated
    """

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    logging.info('creating file with reference sequence')

    # loads data sets
    data_sets = load_all_data_sets(vcf_files + bim_files,
                                   chromosomes)

    # makes the giant table of variants from all files
    variants = merge_data_sets(data_sets, merge_type)

    # makes a single id for each snp
    unify_ids(variants, id_mode)
    
    # gets the reference sequence allele for all snp pos
    reference_alleles = get_reference_alleles(variants,
                                              reference_path)

    # checks if ref allele is present as one of the snp alleles
    check_reference_alleles(reference_alleles)

    # creates for each snp an id of form chrom_pos_a1_a2
    reference_alleles['dataid'] = create_data_id(reference_alleles)

    logging.info('found reference_alleles ')
    return reference_alleles


def load_all_data_sets(files, chromosomes=None):
    """loads all data sets
    
    Parameters
    ----------
    files : list of str
        the names of files to be loaded
    chromosomes : list of str
        chromosomes to be included
    
    Returns
    -------
    data_sets : list of pd.DataFrame
        all the data sets loaded as data frame
    """
    
    data_sets = []
    for f in files:
        data_set = DataFile(f)
        if not data_set.exists():
            raise IOError("file %s not found" % data_set.fname)
        new_data_set = data_set.load_variants()
        new_data_set = filter_data_set(new_data_set, chromosomes)
        
        data_sets.append(new_data_set)

    return data_sets


def filter_data_set(data, chromosomes=None):
    """filter_data_set
    keeps only snp in chromosomes
    
    Parameters
    ----------
    data : pd.DataFrame
        pandas data frame to be filtered
    chromosomes : list of string
        chromosomes to be filtered
    
    Returns
    -------
    filtered_data_set : pd.DataFrame
        input data with only snp on chromosomes in chromosomes
    """
    if chromosomes is None:
        return data
    else:
        mask = np.zeros(len(data), dtype="bool")
        for c in chromosomes:
            mask[np.array(data.chrom == c)] = True
        return data[mask]


def merge_data_sets(data_sets, merge_type="outer"):
    """merges data sets into one big table of snps
    
    Parameters
    ----------
    data_sets : list of pd.DataTable
        a list of the data sets to be merged, return value from
        load_all_data_sets
    merge_type : str
        merge type, see pandas.merge; one of `inner`, `outer`, `left`, `right`
    
    Returns
    -------
    variants: pd.DataTable
        table with all the data sets merged

    """
    
    variants = data_sets[0]
    for i, data_set in enumerate(data_sets[1:]):

        suffixes = '', '_%s' % (i + 1)

        variants = pd.merge(variants, data_set, merge_type,
                            on=('chrom', 'pos', 'a1', 'a2'),
                            suffixes = suffixes)
        variants['pos'] = variants['pos'].astype(np.int64)

    logging.info('read coordinates from all files')

    return variants


def unify_ids(variants, id_mode="left", na_strings=[".", "NA"]):
    """ gets a single id for the merged snp, works inplace
    
    Parameters
    ----------
    variants : pd.DataTable
        table with snp on rows
    id_mode : str
        filter mode, options are `left`, `merge`, and `pos`
        if left, the left most snpid is kept
        if merge, then the snp ids are merged
        if pos, snp ids are generated from chrom_pos_a1_a2 format
    na_strings : list of str
        strings that signify no data
    
    """
    if id_mode not in ["pos", "left", "merge"]:
        raise NotImplementedError("id merge type not supported")

    id_cols2 = [s.startswith('snpid') for s in variants.columns.values]

    new_ids = np.empty(len(variants), dtype="S50")
    for i in xrange(len(new_ids)):
        u = np.unique(variants.ix[i, id_cols2])

        not_na_ids = np.ones_like(u, dtype="bool")
        not_na_ids[u == 'nan'] = False
        for s in na_strings:
            not_na_ids[s == u] = False

        u = u[not_na_ids]
            
        if len(u) == 1:
            new_ids[i] = u[0]

        elif id_mode == "merge":
            new_ids[i] = "+".join(u)
    
        elif id_mode == "left":
            if variants.snpid[i] not in na_strings:
                new_ids[i] = variants.snpid[i]
            else:
                for j in xrange(1, np.sum(id_cols2)):
                    if variants['snpid_%s' % j] not in na_strings:
                        new_ids[i] = variants['snpid_%s' % j][i]
                        break

                else:
                    new_ids[i] = variants.snpid[i]
    
    variants.snpid = new_ids
    for j in xrange(1, np.sum(id_cols2)):
        del variants['snpid_%s' % j]


def get_reference_alleles(pos_alleles, reference_path,
                          sort=False):
    """
    gets the reference allele for a set of positions

    pos_alleles  : structured np.array with 'chrom' and 'pos' denoting
                   chromosome and position, respectively.
    reference_path : path to folder with reference file
    sort : does the array need to be sorted first?

    """
    if sort:
        pos_alleles.sort(order=('chrom', 'pos'))

    ref_call = np.zeros(len(pos_alleles), dtype="S1")

    for chrom in np.unique(pos_alleles['chrom']):
        ref_chr = read_reference_chromosome(chrom, reference_path)

        logging.debug('loaded chromosome %s, %d reference alleles,'
                      ' %d remaining', chrom, len(pos_alleles),
                      np.sum(ref_chr is None is np.nan))
                
        chrom_loc = pos_alleles['chrom'] == chrom
        cur_chrom_pos = pos_alleles['pos'][chrom_loc] - 1
        ref_alleles = [(ref_chr[i]).upper() for i in cur_chrom_pos]
        ref_call[np.array(chrom_loc)] = ref_alleles
    
        logging.debug('added %d new entries from chromosome %s',
                      np.sum(chrom_loc), chrom)

    # end of loop

    pos_alleles['ref'] = ref_call
    
    return pos_alleles


def read_reference_chromosome(chromosome_id, reference_path):
    """read_reference_chromosome
    reads a reference genome fasta file and returns the resulting
    sequence
    
    Parameters
    ----------
    chromosome_id : str
        name/number of the chromosome
    reference_path : str
        path where chromosome files are loaded, we assume one fa.gz file per
        chromosome at location reference_path % chromosome_id
    
    Returns
    -------
    refseq : str
        string with reference sequence
    """

    chromosome_id = str(chromosome_id)
    chromosome_id = chromosome_id.upper()

    if chromosome_id == "MT" or chromosome_id == "90" or chromosome_id == "26":
        chromosome_id = "M"

    if chromosome_id == "23" or chromosome_id == "23":
        chromosome_id = "X"

    if chromosome_id == "24":
        chromosome_id = "Y"

    ref_chr_name = reference_path % chromosome_id
    
    with gzip.open(ref_chr_name) as f:
        f.readline()
        ref1 = f.readlines()

    ref1 = [l.strip() for l in ref1]
    ref1 = "".join(ref1)
    
    return ref1


def check_reference_alleles(reference_alleles):
    """check_reference_alleles

    checks whether the reference allele is one of hte alleles at the
    position
    
    Parameters
    ----------
    reference_alleles : refseq
        data table with columns `a1`, `a2`, and `ref`
    
    Returns
    -------
    None, but adds a col `consistency` to the table
    """
    a1 = reference_alleles['ref'] == reference_alleles['a1']
    a2 = reference_alleles['ref'] == reference_alleles['a2']
        
    con2 = np.logical_or(a1, a2)
    reference_alleles['consistency'] = con2


def create_data_id(ref_seq):
    """create_data_id
    
    Parameters
    ----------
    ref_seq : refseq
        reference allele table
    
    Returns
    -------
    for each snp in the table, an id of form chrom_pos_a1_a2
    """
    return ref_seq.chrom + "_" + \
        ref_seq.pos.astype(str) + "_" + \
        ref_seq.a1 + "_" +\
        ref_seq.a2


def rename_snpids_from_data(tmp_bed_files, twd=""):
    """ renames snpids to chrom_pos_a1_a2 """
    new_tmp_files = []

    for fname in tmp_bed_files:
        if twd == "":
            out_path = "%s_tmp.bim" % fname
        else:
            out_path = utils.path(twd, "%s_tmp.bim" % fname)

        bim = BimFile(fname)
        data = bim.load_variants()
        data['snpid'] = create_data_id(data)
        BimFile.write_file(out_path, data)
        new_tmp_files.append(out_path)
    
    return new_tmp_files


def restore_snpids(bed_file, ref_file):
    """ restores the  in the final merged file """
    snpid_file = load_snpid_file(ref_file)
    bim = BimFile(bed_file)
    data = bim.load_variants()
    
    data2 = data.merge(snpid_file, on='snpid', how="left")
    data2.snpid = data2.oldid
    BimFile.write_file(bed_file + ".bim", data2)


def load_snpid_file(snpid_file):
    data = pd.read_csv(snpid_file, sep="\t", header=0,
                       names=["snpid", "oldid"])
    return data
