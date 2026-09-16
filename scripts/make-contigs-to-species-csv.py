#! /usr/bin/env python
"Build contig <=> species mapping CSV for GTDB representatives."
import sys
import argparse
import polars as pl
import screed


def main():
    p = argparse.ArgumentParser()
    p.add_argument('gtdb_metadata_reps')
    p.add_argument('--genome-dir', required=True)
    p.add_argument('-o', '--output-csv', required=True)
    args = p.parse_args()

    df = pl.read_csv(args.gtdb_metadata_reps)

    contig_xx = []
    for row in df.iter_rows(named=True):
        acc = row['acc']
        species = row['gtdb_taxonomy'].split(';')[-1]

        filename = f'{args.genome_dir}/{acc}_genomic.fna.gz'
        for record in screed.open(filename):
            contig_xx.append(dict(contig=record.name.split(' ')[0],
                                  species=species,
                                  acc=acc))

    contig_df = pl.DataFrame(contig_xx)
    contig_df.write_csv(args.output_csv)


if __name__ == '__main__':
    sys.exit(main())
