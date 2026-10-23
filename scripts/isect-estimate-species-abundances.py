#! /usr/bin/env python
"""
@CTB
"""
import sys
import argparse
import polars as pl
import sourmash
import os.path
import numpy as np
from sklearn.mixture import GaussianMixture
import scipy
import collections

### GMM Code.

GMM_Peak = collections.namedtuple('GMM_Peak', ['idx', 'mean', 'weight', 'covariance'])

def find_best_bic(values, *, verbose=False, max_n_components=5):
    # calculate best num components based on Bayesian Information Criterion score:
    bics = dict()
    for n_components in range(1, max_n_components+1):
        clf = GaussianMixture(n_components=n_components, covariance_type="spherical")
        labels = clf.fit(values).predict(values)
        bics[n_components] = round(float(clf.bic(values)), 2)

    bics = list(bics.items())
    bics.sort(key=lambda x: x[1])
    best_n = bics[0][0]
    if verbose:
        print('BICS:', best_n, bics)

    return best_n


def gmm_fit(values, num_components, *, verbose=False):
    clf = GaussianMixture(n_components=num_components, covariance_type="spherical")
    clf.fit(values)

    orig_peaks = []
    for i in range(num_components):
        mean = float(clf.means_[i][0])
        weight = float(clf.weights_[i])
        cov = float(clf.covariances_[i])
        peak = GMM_Peak(i, mean, weight, cov)
        orig_peaks.append(peak)

    # choose largest mean with weight above 0.1
    peaks = [ p for p in orig_peaks if p.weight >= 0.1]
    peaks.sort(key=lambda x: -x.mean)
    best_peak = peaks[0]

    if verbose:
        for peak in orig_peaks:
            ch = ''
            if peak == best_peak:
                ch = '*'
            print(f'{ch}{peak.idx}: {peak.mean:.1f} {peak.weight:.3f} {peak.covariance:.3f}')

    return best_peak


def fit_isect(ss, *, min_mean=5, min_values=100):
    # extract hash abundances into numpy array from sig.
    values = np.array(list(ss.minhash.hashes.values())).reshape(-1, 1)
    mean = np.mean(values)

    best_peak = None
    if mean >= min_mean and len(values) >= min_values:
        try:
            best_n = find_best_bic(values)
            best_peak = gmm_fit(values, best_n)
        except:
            pass

    if best_peak:
        return best_peak.mean, best_peak.weight

    return None, None


## kmer coverage curve code

def coverage_curve(x, G, C):
    return G*(1 - np.exp(-C*x))

def do_covcurve_fit(contain_df, *, fn=coverage_curve, p0=None):
    "Fit to given function, which takes two parameters, G and C. Return fit values & params."

    if p0 is None:
        p0 = (max(contain_df['weighted_bp']), .001)

    popt, pcov = scipy.optimize.curve_fit(fn, contain_df['weighted_bp'], contain_df['intersect_bp'], p0)
    G, C = popt

    y = np.asarray(contain_df['weighted_bp'])
    yfit = coverage_curve(y, G, C)

    # residual sum of squares
    ss_res = np.sum((y - yfit) ** 2)
    # total sum of squares
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    # r-squared
    r2 = 1 - (ss_res / ss_tot)

    return yfit, (G, C), r2


def covcurve_estimate_values(intersect_bp, weighted_bp, gsize, cexp):
    fit_breadth = intersect_bp / gsize
    fit_depth = weighted_bp / gsize
    curve_breadth = coverage_curve(weighted_bp, gsize, cexp) / gsize

    return fit_breadth, fit_depth, curve_breadth


def main():
    p = argparse.ArgumentParser()
    p.add_argument('isect_zips', nargs='+')
    p.add_argument('--output-kmer-stats-csv', required=False)
    p.add_argument('--output-abund-estimates', required=True)
    p.add_argument('--no-unmangle-names', action='store_true')
    args = p.parse_args()

    # load all the isects and calculate stuff -
    # * record basic k-mer stats
    # * use gmm stuff to estimate k-mer abundances

    # we calculate k-mer coverage curve fit in a second loop

    dflist = []
    for filename in args.isect_zips:
        species_kmers_xx = []
        species_abund_est = []

        for n, ss in enumerate(sourmash.load_file_as_signatures(filename)):
            if n % 500 == 0:
                print('...', n, os.path.basename(filename))
#                if n > 0:
#                    break

            # basic k-mer stats: intersect_bp and weighted_bp
            assert ss.name.endswith(' isect')
            # periods not allowed in species name :)
            species_with_acc, _, metag = ss.name[:-6].split('.', 2)
            if args.no_unmangle_names:
                species = species_with_acc
            else:
                species_with_singlehash  = species_with_acc.split(' ', 1)[1]
                assert species_with_singlehash.endswith(' singlehash')
                species = species_with_singlehash[:-len(' singlehash')]

            mh = ss.minhash
            species_kmers_xx.append(dict(species=species, metag=metag,
                                         intersect_bp=len(mh) * mh.scaled,
                                         weighted_bp=mh.sum_abundances*mh.scaled))
            if n == 0:
                print(species_kmers_xx)

            # also run GMM fit & record
            gmm_mean, gmm_weight = fit_isect(ss)
            if gmm_mean is not None:
                dd = dict(species=species,
                          metag=metag,
                          abund=gmm_mean,
                          method='gmm')
                species_abund_est.append(dd)

        df = pl.DataFrame(species_kmers_xx)
        dflist.append(df)

    all_df = pl.concat(dflist)

    if args.output_kmer_stats_csv:
        all_df.write_csv(args.output_kmer_stats_csv)

    # second loop: 2 rounds of k-mer coverage fit
    for species_name in all_df['species'].unique():
        xx_df = all_df.filter(pl.col('species') == species_name)

        yfit, params, r2 = do_covcurve_fit(xx_df)
        (GSIZE, CEXP) = params

        # use the new information to estimate breadth and depth for every row,
        # and then make a new dataframe

        round2_xx = []
        for row in xx_df.iter_rows(named=True):
            species = row["species"]
            metag = row["metag"]
            intersect_bp = row["intersect_bp"]
            weighted_bp = row["weighted_bp"]

            fit_breadth, fit_depth, curve_breadth = covcurve_estimate_values(
                intersect_bp, weighted_bp, GSIZE, CEXP
            )

            # record for second round fit
            dd = dict(species=species,
                      metag=metag,
                      intersect_bp=intersect_bp,
                      weighted_bp=weighted_bp,
                      fit_breadth=fit_breadth,
                      fit_depth=fit_depth,
                      curve_breadth=curve_breadth)
            round2_xx.append(dd)

            # record in species_abund_est, too
            dd = dict(species=species,
                      metag=metag,
                      abund=fit_depth,
                      method='covcurve_round1')
            species_abund_est.append(dd)

        round2_df = pl.DataFrame(round2_xx)

        # now, do a second round fit for those with low abundances est < 10.
        round2_df = round2_df.filter(pl.col('fit_depth') < 10)
        
        yfit2, params2, r2_2 = do_covcurve_fit(round2_df)
        (GSIZE2, CEXP2) = params2

        # use to estimate depth etc for ALL of the original for this species
        for row in xx_df.iter_rows(named=True):
            species = row["species"]
            metag = row["metag"]
            intersect_bp = row["intersect_bp"]
            weighted_bp = row["weighted_bp"]

            fit_breadth, fit_depth, curve_breadth = covcurve_estimate_values(
                intersect_bp, weighted_bp, GSIZE2, CEXP2
            )

            # record in species_abund_est
            dd = dict(species=species,
                      metag=metag,
                      abund=fit_depth,
                      method='covcurve_round2')
            species_abund_est.append(dd)

    # now, build a composite estimate! where GMM fits exist, use those;
    # otherwise, use round2. Let's see what happens!

    est_df = pl.DataFrame(species_abund_est)

    unique_pairs = est_df.select(['species', 'metag']).unique()

    with pl.Config(tbl_rows=-1):
        print(unique_pairs.sort('metag'))
    
    for row in unique_pairs.iter_rows(named=True):
        species = row["species"]
        metag = row["metag"]
                      
        ff = (pl.col('species') == species) & (pl.col('metag') == metag)
        gmm_est = est_df.filter((pl.col('method') == 'gmm') & ff)
        curve_est = est_df.filter((pl.col('method') == 'covcurve_round2') & ff)
        assert len(curve_est) or len(gmm_est), (species, metag)

        if len(gmm_est):
            composite_abund = gmm_est.select('abund').item()
        else:
            composite_abund = curve_est.select('abund').item()

        species_abund_est.append(dict(species=species, metag=metag,
                                      abund=composite_abund, method='composite'))

    est_df = pl.DataFrame(species_abund_est)
    est_df.write_csv(args.output_abund_estimates)


if __name__ == '__main__':
    sys.exit(main())
