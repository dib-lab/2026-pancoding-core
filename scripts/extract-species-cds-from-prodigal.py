#! /usr/bin/env python
"""
Take taxonomy information and a list of species, load in a pile of
ffn files, and concatenate all the ffn files for the given list of species.
"""
import sys
import argparse
import pathlib
from collections import defaultdict
from sourmash.tax.tax_utils import MultiLineageDB
import polars as pl
import os
import screed
import gzip
import polars as pl


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--prodigal-dirs', nargs='+', required=True)
    p.add_argument('--taxonomy-csvs', nargs='+', required=True)
    p.add_argument('-o', '--output-dir', required=True,
                   help="output concatenated CDS gzip files here")
    p.add_argument('--species-list', required=True,
                   help="list of species")
    p.add_argument('--csv-summary', required=True, type=str,
                   help="output a CSV summary here")
    args = p.parse_args()

    SPECIES=[ x.strip() for x in open(args.species_list) ]

    dflist = []
    for csv in args.taxonomy_csvs:
        dflist.append(pl.read_csv(csv))

    tax_df = pl.concat(dflist)

    print(f'loaded {len(tax_df)} taxonomy entries from lineage sheets.')

    ident_to_ffn = {}

    for dirname in args.prodigal_dirs:
        p = pathlib.Path(dirname)
        ffn_files = p.glob('**/*.ffn')
        n = None
        for n, pp in enumerate(ffn_files):
            if n % 10000 == 0:
                print(f'reading: {dirname} (at {n})')
            ident = pp.name.rsplit('.', 1)[0]
            assert ident not in ident_to_ffn, ident
            ident_to_ffn[ident] = pp.absolute()

            if n == 0:
                print('example found:', ident, ident_to_ffn[ident])

        if n is None:
            print(f"ERROR: '{dirname}' was empty of prodigal files",
                  file=sys.stderr)
            sys.exit(-1)

    print(f'loaded {len(ident_to_ffn)} idents -> ffn files.')

    if len(ident_to_ffn) == 0:
        print("ERROR: loaded no ident-to-ffn mappings!?")
        sys.exit(-1)

    seen = set()
    species_to_idents = defaultdict(list)
    for n, row in enumerate(tax_df.iter_rows(named=True)):
        ident = row["ident"]
        species = row["species"]

        if ident in ident_to_ffn:
            species_to_idents[species].append(ident)
            seen.add(ident)

    if len(seen) != len(ident_to_ffn):
        print(f'ERROR: saw {len(seen)} idents, but {len(ident_to_ffn)} FFN idents.')

    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)

    summary_xx = []
    for species in SPECIES:
        ident_xx = species_to_idents[species]
        print(f'{species} - {len(ident_xx)} genomes.')
        outpath = f'{args.output_dir}/{species}.cds.fa.gz'

        with gzip.open(outpath, 'wt') as fp:
            for ident in ident_xx:
                print('...', ident)
                ffn_file = ident_to_ffn[ident]

                with open(ffn_file, 'r') as infp:
                    data = infp.read()
                    fp.write(data)

                summary_xx.append(dict(species=species,
                                       ident=ident,
                                       ffn_file=str(ffn_file)))

    summary_df = pl.DataFrame(summary_xx)
    summary_df.write_csv(args.csv_summary)


if __name__ == '__main__':
    sys.exit(main())
