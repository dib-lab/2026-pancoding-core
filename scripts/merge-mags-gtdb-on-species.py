#! /usr/bin/env python
"""
Merge MAGs (individual sig zip) and GTDB species (pangenome) into a single
set of signatures.
"""
import sourmash
import sys
import argparse
import polars as pl
from sourmash import sourmash_args


def main():
    p = argparse.ArgumentParser()
    p.add_argument('gtdb_sigzip')
    p.add_argument('mags_sigzip')
    p.add_argument('mags_lincsv')
    p.add_argument('-k', '--ksize', type=int, default=21)
    p.add_argument('-o', '--output', required=True)
    args = p.parse_args()

    gtdb_db = sourmash.load_file_as_index(args.gtdb_sigzip)
    gtdb_db = gtdb_db.select(ksize=args.ksize)
    print(f'GTDB: found {len(gtdb_db)} sketches')

    mags_d = {}
    mags_sigs = sourmash.load_file_as_signatures(args.mags_sigzip,
                                                 ksize=args.ksize)

    print('loading mags...')
    for n, ss in enumerate(mags_sigs):
        if n % 1000 == 0:
            print('...', n)
        mags_d[ss.name] = ss
    print(f'MAGs: loaded {len(mags_d)}')

    lin_df = pl.read_csv(args.mags_lincsv)

    with sourmash_args.SaveSignaturesToLocation(args.output) as save_ss:
        # iterate over all GTDB signatures
        for n, gtdb_ss in enumerate(gtdb_db.signatures()):
            species_mh = gtdb_ss.minhash.to_mutable()

            species = gtdb_ss.name.split(' ', 1)[1]

            # find matching species in MAGs and pull ident
            xx_df = lin_df.filter(pl.col('species') == species)

            # merge all species-matching sketches
            for ident in xx_df['ident']:
                mag_ss = mags_d[ident]
                mag_mh = mag_ss.minhash
                species_mh += mag_mh

            # save updated sketch
            updated_ss = sourmash.SourmashSignature(species_mh, name=gtdb_ss.name)
            if len(gtdb_ss.minhash) != len(species_mh):
                print(n, species, len(gtdb_ss.minhash), len(updated_ss.minhash))
            save_ss.add(updated_ss)

    # done!


if __name__ == '__main__':
    sys.exit(main())
