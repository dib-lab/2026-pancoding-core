#! /usr/bin/env python
import sys
import argparse
import polars as pl
import os.path
import numpy as np
from sklearn.mixture import GaussianMixture
import sourmash


def fit_isect(ss):
    species, _, metag = ss.name.split('.', 3)
    species = " ".join(species.split(' ')[1:3])
    metag = metag.split(' ')[0]

    values = np.array(list(ss.minhash.hashes.values())).reshape(-1, 1)

    try:
        best_n = 2
        clf = GaussianMixture(n_components=best_n, covariance_type="spherical")
        clf.fit(values)
    except:
        return metag, species, -1, -1

    # automatically sort out best.
    xx = sorted(zip(range(2), clf.weights_, clf.means_), key=lambda x: -x[1])
    xx = list(xx)
    #print(xx)

    right_idx = xx[0][0]
    right_weight = xx[0][1]
    right_mean = xx[0][2][0]

    return metag, species, float(right_mean), float(right_weight)


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
