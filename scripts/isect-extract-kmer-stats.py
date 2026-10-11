#! /usr/bin/env python
"""
Use the collected intersections to calculate intersect_bp and weighted_bp for
collections of species-x-metag intersections.
"""
import sys
import argparse
import polars as pl
import sourmash
import os.path


def main():
    p = argparse.ArgumentParser()
    p.add_argument('isect_zips', nargs='+')
    p.add_argument('-o', '--output-csv', required=True)
    args = p.parse_args()

    dflist = []
    for filename in args.isect_zips:
        species_kmers_xx = []
        for n, ss in enumerate(sourmash.load_file_as_signatures(filename)):
            if n % 500 == 0:
                print('...', n, os.path.basename(filename))

            assert ss.name.endswith(' isect')
            species_with_acc, _, metag = ss.name[:-6].split('.')
            species_with_singlehash  = species_with_acc.split(' ', 1)[1]
            assert species_with_singlehash.endswith(' singlehash')
            species = species_with_singlehash[:-len(' singlehash')]
            mh = ss.minhash
            species_kmers_xx.append(dict(species=species, metag=metag,
                                         intersect_bp=len(mh) * mh.scaled,
                                         weighted_bp=mh.sum_abundances*mh.scaled))
            if n == 0:
                print(species_kmers_xx)
            

        df = pl.DataFrame(species_kmers_xx)
        dflist.append(df)

    all_df = pl.concat(dflist)

    all_df.write_csv(args.output_csv)


if __name__ == '__main__':
    sys.exit(main())
