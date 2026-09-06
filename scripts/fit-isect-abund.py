#! /usr/bin/env python
import sys
import argparse
import polars as pl
import os.path
import numpy as np
from sklearn.mixture import GaussianMixture
import sourmash
import collections

GMM_Peak = collections.namedtuple('GMM_Peak', ['idx', 'mean', 'weight', 'covariance'])

def find_best_bic(values, *, verbose=False):
    # calculate best num components based on Bayesian Information Criterion score:
    bics = dict()
    for n_components in range(1, 6):
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


def fit_isect(ss):
    species, _, metag = ss.name.split('.', 3)
    # get rid of accession & 'singlehash'
    species = " ".join(species.split(' ')[1:3])
    metag = metag.split(' ')[0]

    values = np.array(list(ss.minhash.hashes.values())).reshape(-1, 1)
    mean = np.mean(values)

    if mean >= 2 and len(values) > 100:
        try:
            best_n = find_best_bic(values)
            best_peak = gmm_fit(values, best_n)
        except:
            return metag, species, -1, -1
    else:
        return metag, species, mean, 1

    # automatically sort out best.
    return metag, species, best_peak.mean, best_peak.weight


def main():
    p = argparse.ArgumentParser()
    p.add_argument('filenames', nargs='+')
    p.add_argument('-o', '--output-csv', required=True)
    args = p.parse_args()

    xx = []
    n = 0
    for filename in args.filenames:
        for ss in sourmash.load_file_as_signatures(filename):
            if n % 100 == 0:
                print('...', n, args.output_csv)
            metag, species, abund, _ = fit_isect(ss)
            xx.append(dict(metag=metag, species=species, abund=abund))
            n += 1

    abund_df = pl.DataFrame(xx)
    abund_df.write_csv(args.output_csv)


if __name__ == '__main__':
    sys.exit(main())
