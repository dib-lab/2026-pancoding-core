#! /usr/bin/env python
"""
Summarize samtools depth files (from short-read mapping) into depth estimates
for each/all of the reference sequences.
"""
import sys
import argparse
import polars as pl
import glob
import os
import screed
import csv
import screed


def read_depth_txt(depth_file, *, exclude_ends):
    """
    Read an individual depth file in with polars, calculating:

    - len (sliced length, excluding ends)
    - cov (depth numbers for sliced portion)
    - hits (number of bases 'hit' by a read in the sliced portion)

    and then convert them into:
    - depth_all (average across slice)
    - depth_cov (average across hits)
    - breadth (fraction hit)
    """
    basename = os.path.basename(depth_file)
    x_loc = basename.find('.x.')
    assert x_loc > 0
    metag = basename[:x_loc]
    
    df = (pl.scan_csv(depth_file,
                      separator='\t', has_header=False,
                      new_columns=('gene', 'pos', 'cov', 'foo'))
          .select(['gene', 'pos', 'cov'])
          ).collect()

    sum_df = df.group_by('gene').all().with_columns(
        # select slice [75:-75]
        len=(pl.col("pos").list.slice(exclude_ends, -exclude_ends).list.len()),
        cov=(pl.col("cov").list.slice(exclude_ends, -exclude_ends)),
        hits=((pl.col("cov").list.slice(exclude_ends, -exclude_ends).list.filter(pl.element() > 0)).list.len()),
    ).with_columns(
        (pl.lit(metag).alias("metag")),
        # summarize: average depth across contig,
        depth_all=(pl.col("cov").list.sum() / pl.col("len")),
        # average depth across covered bases,
        depth_cov=(pl.col("cov").list.sum() / pl.col("hits")),
        # fraction of bases covered
        breadth=(pl.col("hits") / pl.col("cov").list.len()),
    ).select(["metag", "gene", "len", "hits",
              "breadth", "depth_all", "depth_cov"])
    return sum_df

def main():
    p = argparse.ArgumentParser()
    p.add_argument('depth_files', nargs='+')
    p.add_argument('--exclude-ends', type=int, default=75)
    p.add_argument('-o', '--output-csv', required=True)
    args = p.parse_args()

    depth_files = args.depth_files

    dflist = []
    for n, filename in enumerate(depth_files):
        if n % 10 == 0:
            print(f'reading depth file {n} of {len(depth_files)}')
        sum_df = read_depth_txt(filename, exclude_ends=args.exclude_ends)
        dflist.append(sum_df)
    df = pl.concat(dflist)

    df.write_csv(args.output_csv)


if __name__ == '__main__':
    sys.exit(main())
