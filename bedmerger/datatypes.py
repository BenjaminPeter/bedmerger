import os
import numpy as np
import pandas as pd

from vcf import is_vcf


class DataFile(object):
    """class representing generic data file type to read allele data from,
    subclass this for other data types """
    def __new__(cls, fname):
        if cls is DataFile:
            if is_vcf(fname):
                return super(DataFile, cls).__new__(VCFFile)
            else:
                return super(DataFile, cls).__new__(BimFile)
        else:
            print "not creating"
            return super(DataFile, cls).__new__(cls, fname)
    
    def __init__(self, fname):
        print "data init"
        self.fname = fname

    def load_variants(self):
        raise NotImplemetedError()

    def exists(self):
        """checks if file exists"""
        return os.path.isfile(self.fname)


class VCFFile(DataFile):
    """class representing loader/writer for a vcf file"""
    def __init__(self, fname):
        super(VCFFile, self).__init__(fname)
    
    def load_variants(self):
        """ loads the allele info from a vcf file into memory and gets variant info

        Returns
        -------
        alleles : pd.DataFrame
            a pandas data frame with the alleles
        """
        dt = np.dtype([('chrom', 'S5'), ('pos', np.int64),
                      ('snpid', 'S20'), ('a1', 'S10'), ('a2', 'S10')])

        alleles = np.loadtxt(self.fname, dtype=dt,
                             comments="#", usecols=(0, 1, 2, 3, 4))

        alleles = pd.DataFrame(alleles)
        alleles['pos'] = alleles['pos'].astype(np.int64)

        sort_alleles(alleles)
        return alleles


class BimFile(DataFile):
    """class representing loader/writer for a bim file"""
    def __init__(self, fname):
        super(BimFile, self).__init__(fname)

    def exists(self):
        """checks if file exists"""
        if not os.path.isfile(self.fname + ".bed"):
            return False
        if not os.path.isfile(self.fname + ".bim"):
            return False
        if not os.path.isfile(self.fname + ".fam"):
            return False
        return True

    def load_variants(self):
        """ loads a bim file into memory and gets variant info

        Returns
        -------
        alleles : pd.DataFrame
            a pandas data frame with the alleles
        """
        dt = np.dtype([('chrom', 'S5'), ('pos', np.int64),
                      ('snpid', 'S20'), ('a1', 'S10'), ('a2', 'S10')])
        alleles = np.loadtxt(self.fname + ".bim", dtype=dt,
                             usecols=(0, 3, 1, 5, 4))
        alleles = pd.DataFrame(alleles)
        alleles['pos'] = alleles['pos'].astype(np.int64)

        sort_alleles(alleles)

        return alleles
        
    @staticmethod
    def write_file(out_path, data):
        """writes file to out_path
        
        Parameters
        ----------
        out_path : str
            location to write to
        data : pd.DataFrame
            data to write, typically loaded from load_variants
        
        """
        data['X'] = 0
        cols = ['chrom', 'snpid', 'X', 'pos', 'a1', 'a2']

        data.to_csv(out_path, sep="\t", columns=cols,
                    header=False,
                    index=False,
                    na_rep='NA')

        del data['X']


def sort_alleles( a ):
    """ sorts allele names alphanumerically.

    """
    to_swap = a.a1 > a.a2
    a.a1[to_swap], a.a2[to_swap] = a.a2[to_swap], a.a1[to_swap]

    return a
