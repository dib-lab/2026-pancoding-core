#! /usr/bin/env python
"""
Take dedup output from cd-hit, attribute sequences back to original split
species-level FASTAs, output CSV linking the dedup sequences to species.
"""
import sys
import screed
import hashlib
import argparse
import polars as pl
import os


def calc_md5(record):
    hh = hashlib.md5()
    hh.update(record.sequence.encode('utf-8'))
    return hh.hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('fasta', nargs='+')
    p.add_argument('--species-fa', nargs='+', required=True)
    p.add_argument('-o', '--output-csv', required=True)
    args = p.parse_args()

    d = {}
    for filename in args.fasta:
        for record in screed.open(filename):
            md5 = calc_md5(record)
            assert md5 not in d
            d[md5] = record.name

    md5_to_species = {}
    for species_fa in args.species_fa:
        for record in screed.open(species_fa):
            md5 = calc_md5(record)
            md5_to_species[md5] = species_fa

    xx = []
    found = set()
    for md5, record_name in d.items():
        species_fa=md5_to_species[md5]
        species = os.path.basename(species_fa).split('.')[0]
        ident = record_name.split(' ')[0]
        dd = dict(record_name=record_name,
                  ident=ident,
                  md5=md5,
                  species_fa=species_fa,
                  species=species)
        xx.append(dd)

        found.add(md5)

    xx_df = pl.DataFrame(xx)
    xx_df.write_csv(args.output_csv)

    orig = set(d)
    remaining = orig - found
    assert not remaining, remaining


if __name__ == '__main__':
    sys.exit(main())
