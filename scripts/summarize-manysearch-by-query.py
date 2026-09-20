#! /usr/bin/env python
import argparse
import polars as pl
import sys


def main():
    p = argparse.ArgumentParser()
    p.add_argument('manysearch_pq')
    p.add_argument('--min-intersect', default=3)
    p.add_argument('--min-frequency', default=0.8)
    p.add_argument('-o', '--output-csv', required=True)
    args = p.parse_args()

    # load and do base filtering of intersect_hashes
    base_df = (
        pl.scan_parquet(args.manysearch_pq)
        .filter(pl.col('intersect_hashes') >= args.min_intersect)
    ).collect()

    # summarize by frequency of match
    n_metag = base_df['match_name'].n_unique()
    by_species = base_df.group_by('query_name').agg(
        freq=pl.len() / n_metag
    )

    # filter by frequency
    by_species = by_species.filter(pl.col('freq') >= args.min_frequency)

    # filter entire df
    df = base_df.join(by_species, how='semi', on='query_name')

    # save!
    df.write_csv(args.output_csv)


if __name__ == '__main__':
    sys.exit(main())
