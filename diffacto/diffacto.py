#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
from __future__ import division, print_function

"""diffacto.diffacto: provides entry point main()."""

__version__ = "1.0.3"

import csv
import re
import warnings
from collections import defaultdict
from multiprocessing import Pool

from scipy import optimize, stats
import networkx as nx
import numpy as np
import pandas
from numpy import array, isfinite, nanmean, nansum
from pyteomics import fasta

# from numba import jit  # # Enable just-in-time compiler for speeding up
#@jit
def fast_farms(probes, weight=0.5, mu=0, max_iter=1000,
               force_iter=False, min_noise=1e-4, fill_nan=0.0):
    '''Bayesian Factor Analysis for Proteomics Summarization
       A python translation of function "generateExprVal.method.farms" from
       Bioconductor FARMS.
       [http://www.bioconductor.org/packages/release/bioc/html/farms.html]
       [http://www.bioinf.jku.at/publications/papers/farms/supplementary.ps]

    Reference:
       Hochreiter S, Clevert D and Obermayer K (2006). A new summarization
       method for affymetrix probe level data. Bioinformatics, 22(8),
       http://bioinformatics.oxfordjournals.org/cgi/content/abstract/22/8/943.

    Inputs:
       probes: peptide abundance array (N peptides, M samples) in log scale.
       weight: Hyperparameter (backscale factor) value in the range of [0,1]
                which determines the influence of the prior.
       mu:     Hyperparameter value which allows to quantify different aspects
                of potential prior knowledge. A value near zero assumes that
                most genes do not contain a signal, and introduces a bias for
                loading matrix elements near zero. '''

    readouts = np.array(probes)
    if fill_nan != 0:
        readouts[np.isnan(readouts)] = fill_nan

    # normalize and transform X
    X = np.nan_to_num(readouts).T
    X = X - np.nanmean(X, axis=0)
    xsd = np.nanstd(X, axis=0)
    xsd[xsd < min_noise] = 1.0
    X /= xsd
    X[~isfinite(X)] = 0

    n_samples, n_features = X.shape
    C = np.cov(X.T, ddof=0)

    # positive definite
    C = 0.5 * (C + C.T)
    C[np.where(C < 0)] = 0

    # robustness
    U, s, V = np.linalg.svd(C)
    s[s < min_noise] = min_noise
    C = U.dot(np.diag(s)).dot(V)

    # initiation
    λ = np.sqrt(np.diag(C) * 0.75)
    ψ = np.diag(C) - λ ** 2
    old_psi = ψ
    old_lambda = λ
    alpha = weight * n_features
    E = 1.0
    min_noise_square = min_noise ** 2
    C_diag = np.diag(C)
    for i in range(max_iter):
        # E step
        φ = λ / ψ
        a = 1 + np.matrix(λ) * np.matrix(φ).T
        η = φ / a
        ζ = C.dot(η.T)
        E = 1 - η.dot(λ) + η.dot(ζ)
        # M step
        λ = ζ.T / (E + ψ * alpha)
        λ = np.asarray(λ)[0]
        ψ = C_diag - np.asarray(ζ)[0] * λ + ψ * alpha * λ * (mu - λ)
        ψ = np.maximum(ψ, min_noise_square)
        if ψ[-1] == old_psi[-1] and ψ[0] == old_psi[0] and np.array_equal(ψ, old_psi) and np.array_equal(λ, old_lambda):
            break
        if not force_iter:
            if abs(ψ - old_psi).max() / old_psi.max() < min_noise / 10:
                break
        old_psi = ψ
        old_lambda = λ
    loading = np.sqrt(E[0, 0]) * λ
    φ = loading / ψ
    weights = loading / loading.max()  # rescale loadings to the range of [0,1]
    noise = 1 / (1 + np.matrix(loading) * np.matrix(φ).T)
    return weights, noise[0, 0]


#@jit(nogil=True)
def fast_gmean_nomissing(weights, pep_abd, group_ix):
    abd_w = pep_abd * weights[..., None]
    one_w = abd_w / abd_w * weights[..., None]
    a_sums = np.nansum(abd_w, axis=0)
    w_sums = np.nansum(one_w, axis=0)
    expr = a_sums[group_ix].sum(axis=1) / w_sums[group_ix].sum(axis=1)
    return expr


#@jit(nogil=True)
def sum_squares(pep_abd, group_ix, estimates):
    global nGroups
    residual = 0.0
    for i in range(nGroups):
        res = pep_abd[:, group_ix[i]] - estimates[i]
        residual += np.nansum(res * res)
    return residual


#@jit(nogil=True)
def f_ANOVA(pep_abd, group_ix, estimates, null_ave, dof_loss=0):
    global nGroups
    ss_total = sum_squares(pep_abd, group_ix, null_ave)
    ss_resid = sum_squares(pep_abd, group_ix, estimates)
    dof1 = nGroups - 1
    dof2 = isfinite(pep_abd).sum() - nGroups - dof_loss
    f = ((ss_total - ss_resid) / dof1) / (ss_resid / dof2)
    return f, dof1, dof2


def mv_impute(pep_abd, group_ix, least_missing=0.99, impute_as=0.001):
    ''' Impute missing values when having a large proportion in a sample group.
    Inputs:
        pep_abd:     n peptides, m samples, in linear scale
        group_ix:    grouping index for each of the m samples
        least_missing: set the minimum threshold of missng rate to trigger the
                    imputation (Default: 99%).
        impute_as: set missng values in the sample to this value '''
    aT = np.array(pep_abd).T
    for ix in group_ix:
        if np.isnan(aT[ix]).sum() > least_missing * len(aT[ix].flatten()):
            val = aT[ix]
            val[np.where(np.isnan(val))] = impute_as
            aT[ix] = val
    return aT.T


#@jit(nogil=True)
def weighted_average(weights, pep_abd, group_ix):
    '''
    Calculate weighted geometric means for sample groups
    Inputs:
        weights:    weights of peptides after filtering by loading threshold
        pep_abd:    peptide abundances after filtering by loading threshold
        group_ix:   array indexes of sample groups
    '''
    global nGroups
    abd_w = pep_abd * weights[..., None]
    one_w = abd_w / abd_w * weights[..., None]
    a_sums = np.nansum(abd_w, axis=0)
    w_sums = np.nansum(one_w, axis=0)
    expr = np.empty(nGroups)
    for i in range(expr.shape[0]):
        expr[i] = a_sums[group_ix[i]].sum() / w_sums[group_ix[i]].sum()
    return expr

def _load_fasta(db, id_regex):
    global prot_dict
    prot_dict = dict()
    for header, seq, in fasta.read(db):
        seq = seq.replace('I', 'L').upper()  # convert DB sequence I -> L
        prot_id = header.split()[0]
        if id_regex is not None:
            find_id = re.findall(id_regex, header)
            if len(find_id) > 0:
                prot_id = find_id[0]
        prot_dict[prot_id] = seq

def _map_seq(p):
    pairs = []
    for prot_id, seq, in prot_dict.items():
        if p in seq:
            pairs.append([p, prot_id])
    return pairs


def peptide_db_graph(peps, db, id_regex=None):
    ''' search a set of peptides against a FASTA database  '''
    g = nx.Graph()

    with Pool(initializer = _load_fasta(db, id_regex)) as pool:
        mapped_ppps = pool.map(_map_seq, peps)

    for ppps in mapped_ppps:
        if len(ppps):
            g.add_edges_from(ppps)
    return g


def parsimony_grouping(g, peps):
    ''' Group peptides to proteins using the rule of parsimony
    Inputs:
        g:  an undirected graph with peptide <-> protein as edges
        peps: the set of peptide sequences, nodes not listed in the peptide set
              are protein IDs.
    '''
    not_peps = set(g.nodes()) - set(peps)
    prot_groups = dict()
    for cc in nx.connected_component_subgraphs(g):
        in_group_peptides = set(cc.nodes()) - not_peps
        in_group_proteins = not_peps.intersection(cc.nodes())

        if len(in_group_proteins) == 1:
            prot_groups[in_group_proteins.pop()] = in_group_peptides
        elif len(in_group_proteins) > 1:
            reported = set()
            while len(in_group_proteins - reported) > 0:
                candidate_proteins = sorted(in_group_proteins - reported,
                                            key=lambda p: (
                                                len(set(cc[p].keys()) - reported), p),
                                            reverse=True)
                p = candidate_proteins[0]
                current_peps = set(cc[p].keys())
                plabel = [p]
                for i in range(1, len(candidate_proteins)):
                    _p = candidate_proteins[i]
                    _peps = set(cc[_p].keys())
                    if _peps == current_peps:
                        plabel.append(_p)
                    if len(_peps - current_peps) == 0:
                        reported.add(_p)

                plabel = ';'.join(sorted(plabel))
                if len(current_peps - reported) > 0:
                    prot_groups[plabel] = current_peps
                    reported = reported.union(current_peps)
                reported.add(p)
    return prot_groups


def protein_grouping(df, proteinDb):
    '''
    Grouping peptide sequences in the given dataframe (df)
        by mapping to a protein database (FASTA);
        or by the first column of dataframe when the database is absent
    '''
    peptides = sorted(set(df.index))
    if not proteinDb:
        g = nx.Graph()
        for i, x in df.iterrows():
            for prot in x.values.astype('str')[0].split(';'):
                if len(prot) > 0:
                    g.add_edge(i, prot)
    else:
        g = peptide_db_graph(peptides, proteinDb)
    pg = parsimony_grouping(g, peptides)
    return pg


def zero_center_normalize(df, samples, logInput=False, method='median'):
    '''
    Transforming input peptide abundance table into log2-scale and centralize to zero.
    Inputs:
        df :        dataframe of peptide abundaces
        samples:    column names of selected samples
        logInput:   input abundances are already in log scale
        method:     method for estimating zero point
    '''
    assert method in ('median', 'average', 'GMM'), \
        'Zero centering method has to be among median, average or GMM!'
    if not logInput:
        # convert abundances to log2 scale
        df[samples] = df[samples].apply(np.log2)
    if method == 'average':
        norm_scale = np.nanmean(df[samples], axis=0)
    elif method == 'median':
        norm_scale = np.nanmedian(df[samples], axis=0)
    elif method == 'GMM':
        ''' two-component Gaussian mixture model '''
        from sklearn.mixture import GMM
        gmm = GMM(2)
        norm_scale = []
        for sp in samples:
            v = df[sp].values
            v = v[np.logical_not(np.isnan(v))]
            v = v[np.logical_not(np.isinf(v))]
            try:
                gmm.fit(np.matrix(v.values).T)
                vmean = gmm.means_[np.argmin(gmm.covars_)]
                norm_scale.append(vmean)
            except:
                norm_scale.append(np.nanmean(v))
        norm_scale = np.array(norm_scale)
    df[samples] = df[samples] - norm_scale
    return df


def pqpq(peptide_abundances, metric='correlation', method='complete', t=0.4):
    ''' The essential PQPQ2 process from @yafeng
        [https://github.com/yafeng/pqpq_python/blob/master/pqpq2.py]
    '''
    from scipy import cluster
    d = cluster.hierarchy.distance.pdist(peptide_abundances, metric)
    if metric == "correlation":
        D = np.clip(d, 0, 2)
    else:
        D = d
    L = cluster.hierarchy.linkage(D, method, metric)
    ind = cluster.hierarchy.fcluster(L, t, 'distance')
    return ind


# =====================
# Monte Carlo permutation tests
def monte_carlo_permutation(samp_index, n):
    '''
    Generating a batch of random permutations of sample indexes
    Inputs:
        samp_index: array indexes of sample groups
        n:          size of the batch of permutations
    '''
    flat = np.hstack(samp_index)
    ix = [0]
    [ix.append(ix[-1] + len(i)) for i in samp_index]
    for i in range(n):
        permute = np.random.permutation(flat)
        new_ix = [permute[ix[i - 1]:ix[i]] for i in range(1, len(ix))]
        yield np.array(new_ix)


def calc_q(pvals):
    '''
    Calculationg q-values based on a list of p-values, with a conservative estimate
        of the proportion of true null hypotheses (pi0_hat) based on the given p-values.
    '''
    pv = np.array(pvals)
    pv = pv[isfinite(pv)]
    pi0_hat = min(1, np.sum(pv) * 2 / len(pv))
    ranking = pv.argsort().argsort() + 1
    qlist = pv * pi0_hat * len(pv) / ranking
    for i, rank in enumerate(ranking):
        qlist[i] = min(qlist[ranking >= rank])
    qlist = list(qlist)
    qvals = np.ones_like(pvals)
    for i, e in enumerate(pvals):
        if isfinite(e):
            qvals[i] = qlist.pop(0)
    return qvals


def perform_mcfdr(diffacto_res, sampIx, max_mc=1e5, batch_size=100,
                  terminate_t=50, target_fdr=0.05, sn_threshold=-20):
    '''
    Sequential Monte Carlo permutation test
    Inputs:
        diffacto_res:   a dictionary of Diffacto statistics for each protein
        sampIx:         array indexes of sample groups
        max_mc:         maximun number of random permutations
        batch_size:     number of permutations for every iteration
        terminate_t:    target number of permutation tests with better statistics to terminate the simulation for one protein
        target_fdr:     target level of FDR to stop simulation for the remaining proteins.
        sn_threshold:   sinal-to-noise threshold for exclusion of non-informative proteins.
    '''
    proteins = sorted(diffacto_res.keys())
    preTermination = set()
    for batch in range(1, int(max_mc / batch_size) + 2):
        mc_pvals = []
        for prot in proteins:
            grand_ave, weight, abd_qc, sn, f, T, N = diffacto_res[prot]
            if sn <= sn_threshold:
                mc_pvals.append(np.nan)
                preTermination.add(prot)
                continue
            if prot in preTermination:
                mc_pvals.append(T / N)
                continue
            for ix in monte_carlo_permutation(sampIx, batch_size):
                N += 1
                try:
                    yhat = weighted_average(weight, abd_qc, ix)
                    f_mc, _, _ = f_ANOVA(abd_qc, ix, yhat, grand_ave)
                except:
                    f_mc = f
                if f_mc >= f:
                    T += 1
            diffacto_res[prot][-1] = N  # 1 + Total MC simulations performed
            diffacto_res[prot][-2] = T  # 1 + MC simulations with better stats
            mc_pvals.append(T / N)
            if T >= terminate_t:
                preTermination.add(prot)

        mc_fdr = calc_q(mc_pvals)
        curr_prot = [proteins.index(p) for p in proteins
                     if p not in preTermination]

        if len(curr_prot) == 0 \
                or max(mc_fdr[curr_prot]) < target_fdr \
                or batch * batch_size >= max_mc:
            print("Monte Carlo permutation test finished.")
            return zip(proteins, mc_pvals, mc_fdr)
        else:
            print("%d times simulation, %d proteins remaining (FDR %.3f)" %
                  (batch * batch_size, len(curr_prot), max(mc_fdr[curr_prot])))




# =================================================
#  Main
# =================================================
def main():
    import argparse
    import sys

    DEBUG = False
    SUMMARIZE_EACH_RUN = False
    TOPN = 3
    T_PQPQ = 0.4
    EXAMPLE = 'HUMAN'

    MC_SIMULATION = True
    MC_MAX_N = 200000
    MC_BATCH_SIZE = 100
    MC_MAX_HIT = MC_MAX_N / 1000

    apars = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    apars.add_argument('-i', required=True, nargs=1,
                       help='''Peptides abundances in CSV format.
                               The first row should contain names for all samples.
                               The first column should contain unique peptide sequences.
                               Missing values should be empty instead of zeros.
                       ''')
    # The first column contains unique peptide sequences
    # Missing values should be empty instead of zeros

    apars.add_argument('-db', nargs='?',
                       help='''Protein database in FASTA format.
        If None, the peptide file must have protein ID(s) in the second column.
        ''')

    apars.add_argument('-samples', nargs='?',
                       help='''File of the sample list.
        One run and its sample group per line, separated by tab.
        If None, read from peptide file headings,
           then each run will be summarized as a group.
        ''')

    apars.add_argument('-log2', default='False',
                       help='Input abundances are in log scale (True) or linear scale (False)')

    apars.add_argument('-normalize',
                       choices=['average', 'median', 'GMM', 'None'], default='None',
                       help='Method for sample-wise normalization.')
    # Normalize input abundances (per sample) to zero-centered in log-scale
    # Valid methods include: 'average', 'median' or 'GMM' (two-component
    # Gaussian mixture model).  If None (default), do not normalize.

    apars.add_argument('-farms_mu', type=float, default=0.1,
                       help='Hyperparameter mu')
    # Hyperparameter mu of the FARMS algorithm: prior knowledge of the
    # expected loading.

    apars.add_argument('-farms_alpha', type=float, default=0.1,
                       help='Hyperparameter weight of prior probability')
    # Hyperparameter weight of the FARMS algorithm: weight of prior
    # probability in EM calculation.

    apars.add_argument('-reference', default='average',
                       help='Names of reference sample groups (separated by semicolon)')
    # If average (default) calculate average abundance as the reference.
    # Otherwise, keep peptide abundance values as is.

    apars.add_argument('-min_samples', type=int, default=1,
                       help='Minimum number of samples peptides needed to be quantified in')
    # Peptides quantified in less than the minimum number will be discarded

    apars.add_argument('-use_unique', default='False',
                       help='Use unique peptides only')

    apars.add_argument('-impute_threshold', type=float,
                       default=0.99,
                       help='''Minimum fraction of missing values in the group. Impute missing values if missing fraction is larger than the threshold. ''')

    apars.add_argument('-cutoff_weight', type=float, default=0.5,
                       help='Peptides weighted lower than the cutoff will be excluded')

    apars.add_argument('-fast', default='False',
                       help='Allow early termination in EM calculation when noise is sufficiently small.')

    apars.add_argument('-out', type=argparse.FileType('w'),
                       default=sys.stdout,
                       help='Path to output file (writing in TSV format).')

    apars.add_argument('-mc_out', default=None,
                       help='Path to MCFDR output (writing in TSV format).')

    # ------------------------------------------------
    args = apars.parse_args()

    def boolparam(p):
        ''' convert a string parameter to boolean value'''
        if str(p).lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        else:
            return False

    args.log2 = boolparam(args.log2)
    args.fast = boolparam(args.fast)
    args.use_unique = boolparam(args.use_unique)
    print(args)
    diffacto_res = dict()
    df = pandas.read_csv(args.i[0], index_col=0)
    df.index = [i.upper().replace('I', 'L') for i in df.index]
    print("Abundance matrix loaded: %d peptides" % len(df.index))

    if not args.samples:  # read sample names from header
        samples = df.columns.tolist()
        if args.db is None:
            samples.pop(0)
        groups = samples
    else:  # read sample labels
        samples, groups = ([], [])
        with open(args.samples) as fh:
            for line in fh.readlines():
                try:
                    _s, _g = line.rstrip().split('\t')
                    samples.append(_s)
                    groups.append(_g)
                except ValueError:
                    pass

    # per sample normalization of peptide abundances
    logInput = args.log2
    if not args.normalize == 'None':
        df = zero_center_normalize(
            df, samples, logInput=logInput, method=args.normalize)
        args.log2 = True

    # select reference runs if specified
    ref_samples = []
    if args.reference:
        for r in args.reference.split(';'):
            for i in range(len(groups)):
                if groups[i] == r:
                    ref_samples.append(i)
    ref_samples = [samples[i] for i in ref_samples]

    print("Number of runs: %d" % len(samples))

    # sample grouping
    group_names = [i for i in sorted(set(groups),
                                     key=lambda k: "{0:0>50}".format(k))
                   if i not in args.reference.split(';')]
    if len(group_names) == len(samples):
        group_names = samples

    sampIx = np.array([[j for j in range(len(groups)) if groups[j] == i]
                       for i in group_names])
    global nGroups
    nGroups = len(group_names)
    print("Number of sample groups: %d" % nGroups)
    print("Reference runs (%d): " % len(ref_samples), *ref_samples, sep='\t')

    # protein grouping
    pg = protein_grouping(df, args.db)
    print("Number of protein groups: %d" % len(pg.keys()))

    # coverage filtering
    df = df[[np.count_nonzero(np.nan_to_num(v)) >=
             args.min_samples for v in df[samples].values]]

    # reversed mapping (peptide to protein group) for checking peptide
    # uniqueness.
    pep2prot = defaultdict(list)
    for prot_ids, bseqs in pg.items():
        for s in bseqs:
            pep2prot[s] += prot_ids.split()

    # use unique peptides
    if args.use_unique:
        df = df[[len(pep2prot[p]) == 1 for p in df.index]]

    #Check that we don't have any peptides with a single non-missing value. These tend to break diffacto, because in fast_farms we end up with a covariance matrix of less than full rank. Which the algorithm is not set up to handle.
    nonZeroNonMissing = np.vectorize(lambda x : ~np.isnan(x) and x > 0, otypes = [np.bool])
    if df.shape[0] > 0:
        for prot in sorted(pg.keys()):
            if prot == 'nan':
                continue
            if DEBUG and EXAMPLE not in prot:
                continue
            # =====----=====-----=====-----=====
            peps = pg[prot]  # constituent peptides
            dx = df.ix[[p for p in sorted(peps) if p in df.index]]  # dataframe
            pep_count = len(dx)  # number of peptides
            pep_abd = dx[samples].values
            counts = np.sum(nonZeroNonMissing(pep_abd), axis = 1)
            if any(counts < 2):
                print("Protein {} contained peptides with fewer than two non-missing or non-zero values. Please remove these peptides".format(prot))
                return

    # -------------------------------------------------------------------------
    # perform differential analysis
    output_header = ['Protein', 'N.Pept', 'Q.Pept', 'S/N', 'P(PECA)']
    output_header += group_names
    if SUMMARIZE_EACH_RUN:
        output_header += ['P(Top-%d)' % TOPN, 'P(Median)', 'P(PQPQ)']
        output_header += ["Top-%d_%s" % (TOPN, s) for s in samples]
        output_header += ["Median_%s" % s for s in samples]
        output_header += ["PQPQ_%s" % s for s in samples]

    print(*output_header, sep='\t', file=args.out)
    for prot in sorted(pg.keys()):
        if prot == 'nan':
            continue
        if DEBUG and EXAMPLE not in prot:
            continue
        # =====----=====-----=====-----=====
        peps = pg[prot]  # constituent peptides
        dx = df.ix[[p for p in sorted(peps) if p in df.index]]  # dataframe
        pep_count = len(dx)  # number of peptides
        pep_abd = dx[samples].values

        if len(ref_samples):  # rescale peptide abundances by reference runs
            reference_abundance = dx[ref_samples].mean(
                axis=1).fillna(np.nanmean(dx[samples])).values
        elif args.reference.lower() == 'average':  # rescale by average values
            reference_abundance = dx[samples].mean(axis=1).values
        else:
            if not args.log2:
                reference_abundance = 1.0
            else:
                reference_abundance = 0

        if not args.log2:
            pep_abd = np.log2(pep_abd)
            reference_abundance = np.log2(reference_abundance)

        pep_abd = (pep_abd.T - reference_abundance).T

        if pep_count == 1:  # single peptide group
            loading = array([1 for _ in dx.index])
            noise = 1.0
            continue  # do not report
        elif pep_count > 1:
            loading, noise = fast_farms(pep_abd,
                                        mu=args.farms_mu,
                                        weight=args.farms_alpha,
                                        max_iter=1000,
                                        force_iter=not args.fast)
        else:
            continue

        sn = 10 * np.log10((1 - noise) / noise)
        qc = loading > args.cutoff_weight
        abd_qc = mv_impute(pep_abd[qc], sampIx,
                           least_missing=args.impute_threshold,
                           impute_as=np.nanmin(pep_abd) - 1)
        protein_summary_group = weighted_average(loading[qc], abd_qc, sampIx)

        if SUMMARIZE_EACH_RUN:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                # Top-N averaging
                v = dx[samples].values
                if logInput:
                    v = 2 ** v
                protein_summary_topn = np.array([np.mean(np.sort(v[:, i][isfinite(v[:, i])])[-TOPN:])
                                                 for i in range(len(samples))])
                p_ave = stats.f_oneway(*[protein_summary_topn[s][isfinite(protein_summary_topn[s])]
                                         for s in sampIx])[1]

                # Median
                v = dx[samples].values
                protein_summary_median = np.nanmedian(v, axis=0)
                p_med = stats.f_oneway(
                    *[protein_summary_median[s][isfinite(protein_summary_median[s])] for s in sampIx])[1]

                # PQPQ clustering and averaging
                v = np.nan_to_num(pep_abd)
                clusters = pqpq(v, t=T_PQPQ)
                major = sorted([(len(clusters[clusters == i]), i)
                                for i in set(clusters.tolist())])[-1]
                if major[0] >= 2:
                    clusters[clusters != major[1]] = 0
                    clusters[clusters != 0] = 1
                else:
                    clusters = np.ones(*clusters.shape)
                protein_summary_pqpq = np.nanmean(
                    dx[samples].values[clusters > 0], axis=0)
                p_pqpq = stats.f_oneway(*[protein_summary_pqpq[s][isfinite(protein_summary_pqpq[s])]
                                          for s in sampIx])[1]

        # ================================================================
        # PECA: grouping peptide-level p-values based on beta distribution
        # https://www.bioconductor.org/packages/release/bioc/html/PECA.html
        '''
        Calculates Probe-level Expression Change Averages (PECA)
        to identify differential expression in Affymetrix gene expression
        microarray studies or in proteomic studies using peptide-level
        mesurements respectively.
        '''
        pep_pvals = []
        for pep_v in abd_qc:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                ave_0 = np.broadcast_to(np.nanmean(pep_v), (nGroups, 1))
                ave_1 = np.array([np.nanmean(pep_v[i]) for i in sampIx])
            try:
                f, d1, d2 = f_ANOVA(pep_v[None, ...], sampIx, ave_1, ave_0)
                pv = stats.f.sf(f, d1, d2)
                pep_pvals.append(pv)
            except:
                pass

        pep_pvals = np.array(pep_pvals)
        pep_pvals = pep_pvals[isfinite(pep_pvals)]
        beta_ab = len(pep_pvals) / 2 + 0.5
        p_peca = stats.beta.cdf(np.median(pep_pvals), beta_ab, beta_ab)

        if MC_SIMULATION:
            grand_ave = np.broadcast_to(np.nanmean(abd_qc), (nGroups, 1))
            f, _, _ = f_ANOVA(abd_qc, sampIx, protein_summary_group, grand_ave)
            diffacto_res[prot] = [grand_ave, loading[qc], abd_qc, sn, f, 1, 1]

        # =============================
        if not logInput:
            protein_summary_group = 2 ** protein_summary_group

        output_row = [prot, pep_count, sum(qc), sn, p_peca] \
            + list(protein_summary_group)

        if SUMMARIZE_EACH_RUN:
            output_row += [p_ave, p_med, p_pqpq] \
                + list(protein_summary_topn) \
                + list(protein_summary_median) \
                + list(protein_summary_pqpq)

        print(*output_row, sep='\t', file=args.out)

    if MC_SIMULATION and args.mc_out:
        try:
            mc_out = open(args.mc_out, 'w')
        except:
            print('Cannot open file: ', args.mc_out, '. Use stdout instead.')
            mc_out = sys.stdout

        print('Protein', 'P(MC)', 'MCFDR', sep='\t', file=mc_out)
        mc_result = perform_mcfdr(diffacto_res, sampIx,
                                  max_mc=MC_MAX_N,
                                  batch_size=MC_BATCH_SIZE,
                                  terminate_t=MC_MAX_HIT,
                                  target_fdr=0.05)

        for prot, p, q in mc_result:
            print(prot, p, q, sep='\t', file=mc_out)

if __name__ == '__main__':
    main()
