#! /usr/bin/env python
"""
"""
import sys
import argparse
import os, os.path
import glob
from collections import defaultdict


import sourmash
from sourmash.save_load import SaveSignaturesToLocation


def main():
    p = argparse.ArgumentParser()
    p.add_argument('input_dir')
    p.add_argument('-o', '--output-dir', required=True)
    p.add_argument('-n', '--dry-run', action='store_true')
    args = p.parse_args()

    if args.dry_run:
        print(f'DRY-RUN MODE. Nothing will be written.')

    args.input_dir = args.input_dir.rstrip('/')
    args.output_dir = args.output_dir.rstrip('/')

    print(f'{args.input_dir}/*.x.*.sig.zip')
    xx = glob.glob(f'{args.input_dir}/*.x.*.sig.zip')
    print(f'found {len(xx)} individual intersection sigs')
    assert xx

    metags_set = set()
    species_d = defaultdict(list)

    for filename in xx:
        x = os.path.basename(filename)
        x = x[:-len('.sig.zip')]
        midpos = x.find('.x.')
        assert midpos > 0

        species = x[:midpos]
        metag = x[midpos + 3:]
        metags_set.add(metag)

        # remove accession + 'singlehash'
        species = " ".join(species.split(' ')[1:3])

        assert species.startswith('s__'), species

        # track files by species
        species_d[species].append(filename)

    print(f"found {len(species_d)} species x {len(metags_set)} metags")
    assert len(species_d) * len(metags_set) == len(xx)

    if not args.dry_run:
        try:
            os.mkdir(args.output_dir)
        except FileExistsError:
            pass

    print(f"saving to species-specific collections...")
    for species in species_d:
        print(f'working on {species}')
        output_name = f"{args.output_dir}/{species}.collected.sig.zip"

        if not args.dry_run:
            with SaveSignaturesToLocation(output_name) as save_ss:
                for isect_file in species_d[species]:
                    for ss in sourmash.load_file_as_signatures(isect_file):
                        save_ss.add(ss)
            

if __name__ == '__main__':
    sys.exit(main())
