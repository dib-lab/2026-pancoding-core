#! /usr/bin/env python
"""
CTB todo: output a single parquet file? OR, separate script to
combine parquet files into one??
"""
import polars as pl
import sys
import argparse
import os.path
import traceback


def main():
    p = argparse.ArgumentParser()
    p.add_argument('csvs', nargs='+')
    p.add_argument('-o', '--outdir', help='put parquet files in this directory')
    p.add_argument('--redo-existing', action='store_true')
    p.add_argument('-f', '--force', action='store_true')
    args = p.parse_args()

    failed = []
    for n, csvfile in enumerate(args.csvs):
        assert csvfile.endswith('.csv')
        pq_file = csvfile[:-4] + '.parquet'
        if args.outdir:
            pq_file = os.path.join(args.outdir, os.path.basename(pq_file))

        if os.path.exists(pq_file) and not args.redo_existing:
            print(f"{n+1}: parquet file already exists: '{pq_file}'; skipping.")
            continue

        try:
            df = pl.scan_csv(csvfile, ignore_errors=args.force)
            df.sink_parquet(pq_file)
            print(f"{n+1}: saved to '{pq_file}'")
        except pl.exceptions.NoDataError:
            print(f"{n}: failed '{csvfile}': empty!'")
            failed.append(csvfile)
            continue
        except:
            print(f"{n+1}: failed '{csvfile}': ???!'")
            traceback.print_exc()
            failed.append(csvfile)
            continue

    if failed:
        with open('failed.txt', 'wt') as fp:
            fp.write("\n".join(failed))
        print(f"wrote {len(failed)} failed CSV file names to failed.txt")


if __name__ == '__main__':
    sys.exit(main())
