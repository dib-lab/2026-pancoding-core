#! /usr/bin/env python
import argparse
import polars as pl
import sys
from scipy.stats import binom
import sourmash


def p_of_match(query_size, metag_size, n_matches, ksize, scaled):
    r = 2*(4.0 ** (-ksize)) * metag_size * scaled

    prob = binom.sf(n_matches - 1, query_size, r)
    return float(prob)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('manysearch_csv')
    p.add_argument('--queries', nargs='+', help='query sketches',
                   required=True)
    p.add_argument('--against', nargs='+', help='against sketches',
                   required=True)
    p.add_argument('-o', '--output-csv')
    p.add_argument('-P', '--output-parquet')
    args = p.parse_args()

    if not (args.output_csv or args.output_parquet):
        print("ERROR: no output specified!? Give either -o or -P",
              file=sys.stderr)
        sys.exit(-1)

    df = pl.read_csv(args.manysearch_csv)
    print(f"loaded {len(df)} rows from '{args.manysearch_csv}'")
    assert 'match_name' in df.columns
    assert 'query_name' in df.columns
    assert 'intersect_hashes' in df.columns
    assert 'ksize' in df.columns
    assert 'scaled' in df.columns

    first_row = next(iter(df.iter_rows(named=True)))
    scaled = first_row["scaled"]
    ksize = first_row["ksize"]
    print(f"scaled={scaled} ksize={ksize}")

    # load query sizes
    query_sizes = {}
    for qfile in args.queries:
        db = sourmash.load_file_as_index(qfile)
        db.select(ksize=ksize)
        mf = db.manifest
        for row in mf.rows:
            name = row["name"]
            size = row["n_hashes"]
            query_sizes[name] = size

    # load against sizes
    match_sizes = {}
    for afile in args.against:
        db = sourmash.load_file_as_index(afile)
        db.select(ksize=ksize)
        mf = db.manifest
        for row in mf.rows:
            name = row["name"]
            size = row["n_hashes"]
            match_sizes[name] = size

    xx = []
    for row in df.iter_rows(named=True):
        query_name = row["query_name"]
        match_name = row["match_name"]
        query_size = query_sizes[query_name]
        match_size = match_sizes[match_name]

        intersect_hashes = row["intersect_hashes"]
        pval = p_of_match(query_size, match_size, intersect_hashes,
                          ksize, scaled)
        xx.append(dict(query_name=query_name, match_name=match_name, pval=pval))

    pval_df = pl.DataFrame(xx)

    a = len(df)
    df = df.join(pval_df, on=['query_name', 'match_name'], how='inner',
                 validate="1:1")
    b = len(df)
    assert a == b, (a, b)

    if args.output_csv:
        print(f"wrote CSV to '{args.output_csv}'")
        df.write_csv(args.output_csv)

    if args.output_parquet:
        print(f"wrote parquet to '{args.output_parquet}'")
        df.write_parquet(args.output_parquet)


if __name__ == '__main__':
    sys.exit(main())
