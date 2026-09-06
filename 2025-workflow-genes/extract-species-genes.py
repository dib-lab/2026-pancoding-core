#! /usr/bin/env python
"""
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
    p.add_argument('csv')
    p.add_argument('fasta')
    p.add_argument('-o', '--output', required=True)
    p.add_argument('--save-md5-csv', required=True)
    args = p.parse_args()

    df = pl.read_csv(args.csv).filter(pl.col('good') == 1)
    gene_names = set(df['gene_name'])
    found = set()

    xx = []

    with open(args.output, 'wt') as fp:
        for record in screed.open(args.fasta):
            ident = record.name.split(' ')[0]
            if ident in gene_names:
                found.add(ident)
                fp.write(f'>{ident}\n{record.sequence}\n')

                xx.append(dict(ident=ident,
                               seq=record.sequence,
                               md5=calc_md5(record)))

    assert found == gene_names
    print(f'found {len(found)}')

    md5_df = pl.DataFrame(xx)
    md5_df.write_csv(args.save_md5_csv)


if __name__ == '__main__':
    sys.exit(main())
