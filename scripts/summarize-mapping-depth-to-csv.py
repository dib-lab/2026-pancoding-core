#! /usr/bin/env python
import sys
import argparse
import polars as pl
import glob
import os
import screed
import csv
import screed


def read_depth_txt(depth_file, *, exclude_ends):
    basename = os.path.basename(depth_file)
    metag = basename.split('.')[0]
    
    df = (pl.scan_csv(depth_file,
                      separator='\t', has_header=False,
                      new_columns=('gene', 'pos', 'cov', 'foo'))
          .select(['gene', 'pos', 'cov'])
          ).collect()

    sum_df = df.group_by('gene').all().with_columns(
        # select slice [75:-75]
        (pl.col("pos").list.slice(exclude_ends, -exclude_ends).list.len()).alias("len"),
        (pl.col("cov").list.slice(exclude_ends, -exclude_ends)),
        (pl.col("cov").list.slice(exclude_ends, -exclude_ends).list.filter(pl.element() > 0)).list.len().alias("hits"),
    ).with_columns(
        (pl.lit(metag).alias("metag")),
        # summarize: average depth across contig,
        (pl.col("cov").list.sum() / pl.col("len")).alias("depth_all"),
        # average depth across covered bases,
        (pl.col("cov").list.sum() / pl.col("hits")).alias("depth_cov"),
        # fraction of bases covered
        (pl.col("hits") / pl.col("cov").list.len()).alias("breadth"),
    ).select(["metag", "gene", "len", "hits",
              "breadth", "depth_all", "depth_cov"])
    #xx_df = xx_df.join(seq_to_species_df, on='gene', how='inner')
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
