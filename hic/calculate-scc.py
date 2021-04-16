import sys, cooler, joblib
from typing import Optional, Iterable, List, Tuple
import hicreppy.utils.mat_process as cu
import numpy as np
import scipy.stats as ss
from scipy.sparse import SparseEfficiencyWarning
import warnings

def get_scc(
    mat1: 'scipy.sparse.csr_matrix', mat2: 'scipy.sparse.csr_matrix', max_bins: int
) -> float:
    """
    Compute the stratum-adjusted correlation coefficient (SCC) between two
    Hi-C matrices up to max_dist. A Pearson correlation coefficient is computed
    for each diagonal in the range of 0 to max_dist and a weighted sum of those
    coefficients is returned.
    Parameters
    ----------
    mat1 : scipy.sparse.csr_matrix
        First matrix to compare.
    mat2 : scipy.sparse.csr_matrix
        Second matrix to compare.
    max_bins : int
        Maximum distance at which to consider, in bins.
    Returns
    -------
    scc : float
        Stratum adjusted correlation coefficient.
    """
    corr_diag = np.zeros(len(range(max_bins)))
    weight_diag = corr_diag.copy()
    for d in range(max_bins):
        d1 = mat1.diagonal(d)
        d2 = mat2.diagonal(d)
        # Silence NaN warnings: this happens for empty diagonals and will
        # not be used in the end.
        with warnings.catch_warnings():
            # Warning does not exist in older scipy versions (<=1.2)
            try:
                warnings.filterwarnings(
                    "ignore", category=ss.PearsonRConstantInputWarning
                )
            except AttributeError:
                pass
            # Compute raw pearson coeff for this diag
            corr_diag[d] = ss.pearsonr(d1, d2)[0]
        # Compute weight for this diag
        r2k = cu.vstrans(d1, d2)
        weight_diag[d] = len(d1) * r2k
    # Normalize weights
    weight_diag /= sum(weight_diag)

    # Weighted sum of coefficients to get SCCs
    scc = np.sum(corr_diag * weight_diag)

    return scc

def make_chromlist(
    c: 'cooler.Cooler',
    whitelist: Optional[Iterable[str]] = None,
    blacklist: Optional[Iterable[str]] = None,
    min_size: Optional[int] = None,
) -> Tuple[List[str], List[int]]:
    """Given a cool object, a blacklist and whitelist of chromosomes, return
    the list of chromosomes to include in the analysis.
    Parameters
    ----------
    c : cooler.Cooler
        First matrix to compare.
    whitelist : None or list of strs
        If given, only compare those chromosomes.
    blacklist : None or list of strs
        If given, do not compare those chromosomes.
    min_size : int
        Chromosomes smaller than this value will be removed.
    Returns
    -------
    chromlist : list of strs
        Names of chromosomes to process
    chromlist : list of strs
        Lengths of chromosomes to process
    """
    chroms = c.chroms()[:].set_index("name")
    chromlist = chroms.index.values.tolist()
    if whitelist is not None:
        chromlist = whitelist
    if blacklist is not None:
        for black in blacklist:
            chromlist.remove(black)
    # Remove small chromosomes
    orig_chromlist = np.array(chromlist).copy()
    if min_size is not None:
        for chrom in orig_chromlist:
            if chroms.loc[chrom].length < min_size:
                chromlist.remove(chrom)
    n_removed = len(orig_chromlist) - len(chromlist)
    if n_removed > 0:
        print(
            f"Removed {n_removed} chromosomes shorter than max_dist",
            file=sys.stderr,
        )
    chrom_lengths = chroms.loc[chromlist].length.values.tolist()
    return chromlist, chrom_lengths

def array2sparse(arr):

    from scipy.sparse import coo_matrix

    x, y = np.nonzero(arr)
    value = arr[x, y]
    shape = arr.shape

    M = coo_matrix((value, (x, y)), shape=shape)

    return M

def genome_scc(
    mat1: 'cooler.Cooler',
    mat2: 'cooler.Cooler',
    max_dist: int,
    h: int,
    correct1: str,
    correct2: str,
    whitelist: Optional[Iterable[str]] = None,
    blacklist: Optional[Iterable[str]] = None,
) -> list:
    
    if mat1.binsize != mat2.binsize:
        raise ValueError("Both matrices must be binned at the same resolution")

    max_bins = max_dist // mat1.binsize

    # Define chromosomes to scan
    # NOTE: chromosomes smaller than the kernel used for smoothing or the
    # max_dist must be removed.
    min_size = max((2 * h + 1) * mat1.binsize, max_dist)
    chromlist, chroms_lengths = make_chromlist(
        mat1, whitelist, blacklist, min_size=min_size
    )
    if not len(chromlist):
        raise KeyError(
            "All chromosomes were too short and have been discarded. "
            "Try reducing max_dist."
        )

    # Compute SCC values separately for each chromosome
    chroms_scc = np.zeros(len(chromlist))
    for c, chrom in enumerate(chromlist):
        if correct1=='False':
            chrom_1 = mat1.matrix(sparse=True, balance=False).fetch(chrom)
        else:
            chrom_1 = mat1.matrix(sparse=False, balance=correct1).fetch(chrom)
            chrom_1[np.isnan(chrom_1)] = 0
            chrom_1 = array2sparse(chrom_1)
        if correct2=='False':
            chrom_2 = mat2.matrix(sparse=True, balance=False).fetch(chrom)
        else:
            chrom_2 = mat2.matrix(sparse=False, balance=correct2).fetch(chrom)
            chrom_2[np.isnan(chrom_2)] = 0
            chrom_2 = array2sparse(chrom_2)
        
        # Trim diagonals which are too far to be scanned to reduce
        # compute time and memory usage
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=SparseEfficiencyWarning)
            chrom_1 = cu.diag_trim(chrom_1.todia(), max_bins + h).tocoo()
            chrom_2 = cu.diag_trim(chrom_2.todia(), max_bins + h).tocoo()
        smooth_1 = cu.smooth(chrom_1, h)
        smooth_2 = cu.smooth(chrom_2, h)

        chroms_scc[c] = get_scc(smooth_1, smooth_2, max_bins=max_bins)
    # Compute the genome SCC using the weighted averge of chromosomes
    # SCC by their lengths. NaN values of SCC are not considered
    # This happens when comparing empty diagonals
    nan_scc_mask = ~np.isnan(chroms_scc)
    trunc_scc = chroms_scc[nan_scc_mask]
    trunc_lengths = np.array(chroms_lengths)[nan_scc_mask]
    scc = np.average(trunc_scc, weights=trunc_lengths)
    
    return scc

cools = [
    'K562-R2.bwa.mcool',
    'K562-R1.bwa.mcool',
    'K562-R2.chromap.mcool'
]

res = 25000
h = 5
max_dist = 2000000

scc = np.zeros((len(cools), len(cools)))
for i in range(len(cools)):
    for j in range(len(cools)):
        if i < j:
            q1 = cools[i]
            q2 = cools[j]
            uri1 = '{0}::resolutions/{1}'.format(q1, res)
            uri2 = '{0}::resolutions/{1}'.format(q2, res)
            clr1 = cooler.Cooler(uri1)
            clr2 = cooler.Cooler(uri2)
            tmp = genome_scc(clr1, clr2, max_dist=max_dist, h=h, correct1='weight', correct2='weight', blacklist=['chrY'])
            scc[i, j] = tmp
            scc[j, i] = tmp
    scc[i, i] = 1

joblib.dump(scc, 'scc-matrix.pkl')

