#! /usr/bin/env python
"""
Produce a manysketch CSV from a set of prodigal files.
"""
import sys
import argparse
from pathlib import Path
import polars as pl


def main():
    p = argparse.ArgumentParser()
    p.add_argument('ffn_dir')
    p.add_argument('-o', '--output-csv', required=True)
    args = p.parse_args()

    ffn_dir = args.ffn_dir.rstrip('/')

    ffn_files = Path(ffn_dir).glob('*.ffn')
    ffn_files = list(ffn_files)

    print(f'found {len(ffn_files)} files at path {ffn_dir}')

    xx = []
    for filename in ffn_files:
        ident = filename.name[:-4] # strip off .ffn
        print(filename.absolute())
        xx.append(dict(name=ident,
                       genome_filename=str(filename.absolute()),
                       protein_filename='')
                  )
                  

    xx_df = pl.DataFrame(xx)
    xx_df.write_csv(args.output_csv)


if __name__ == '__main__':
    sys.exit(main())
