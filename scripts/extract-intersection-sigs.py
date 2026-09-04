#! /usr/bin/env python
"""
Compute many intersections at once between species sketches and
metagenomes, efficiently (but still in Python).

Loads all species sketches once, then loads each metagenome, calculates
intersection, and inflates the intersection with the metagenome abundances.

Does approximately the same thing as sourmash sig intersect ... -A ...,
but much more efficiently for large operations.

Outputs to either a directory full of nicely-named files (-o), or a zip file
(with -O) with nicely named signatures.

@CTB TODO:
- resume? if -O
- optimize? do intersection once, inflation once, then sub?
"""
import sys
import argparse
import os


import sourmash
from sourmash.save_load import SaveSignaturesToLocation

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--query-sigs', nargs='+', required=True,
                   help='list of genome/species signatures')
    p.add_argument('--metagenomes', required=True,
                   help='list of paths to metagenome signatures')
    p.add_argument('-k', '--ksize', default=21, type=int)
    p.add_argument('--scaled', default=1000, type=float)
    p.add_argument('-o', '--output-dir',
                   help='save intersection signatures to this directory as individual files')
    p.add_argument('-O', '--output-file',
                   help='save intersection signatures to this file (e.g. a zip file)')
    args = p.parse_args()

    if args.output_dir and args.output_file:
        print("ERROR: can only specify one of --output-dir and --output-file",
              file=sys.stderr)
        sys.exit(-1)

    if not (args.output_dir or args.output_file):
        print("ERROR: must specify one of --output-dir or --output-file",
              file=sys.stderr)
        sys.exit(-1)

    if args.output_dir:
        try:
            os.mkdir(args.output_dir)
        except FileExistsError:
            pass

    bulk_save = None
    if args.output_file:
        bulk_save = SaveSignaturesToLocation(args.output_file)
        bulk_save.open()

    scaled = int(args.scaled)
    ksize = args.ksize

    query_mh_list = []
    merged_query_mh = None      # build a merged set of k-mers for optimize
    for q in args.query_sigs:
        db = sourmash.load_file_as_index(q)
        db = db.select(ksize=ksize)
        for ss in db.signatures():
            mh = ss.minhash
            name = ss.name
            query_mh = mh.downsample(scaled=scaled).flatten()
            query_mh_list.append((name, query_mh))

            if merged_query_mh is None:
                merged_query_mh = query_mh.to_mutable()
            else:
                merged_query_mh += query_mh
                

    print(f'loaded {len(query_mh_list)} query signatures')

    lines = [ x.strip() for x in open(args.metagenomes) ]

    for n, metagenome_file in enumerate(lines):
        print(f'{n+1}/{len(lines)} {os.path.basename(metagenome_file)}')
        sigs = sourmash.load_file_as_signatures(metagenome_file,
                                                ksize=ksize)

        # for every metagenome,
        for ss in sigs:
            metag_name = ss.name
            metag_mh = ss.minhash.downsample(scaled=scaled)
            metag_flat_mh = metag_mh.flatten()

            # intersect with merged queries
            metag_flat_mh = metag_flat_mh.intersection(merged_query_mh)

            # inflate that intersection from the original
            metag_mh = metag_flat_mh.inflate(metag_mh)

            # loop across all query sketches,
            for (query_name, query_mh) in query_mh_list:
                # find intersection & borrow abundances from metagenome,
                isect_mh = metag_flat_mh.intersection(query_mh)
                isect_mh = isect_mh.inflate(metag_mh)

                # then create new sig and save.
                isect_name = f'{query_name}.x.{metag_name} isect'
                isect_ss = sourmash.SourmashSignature(isect_mh, isect_name)

                # save to either file-in-directory, or add to large file.
                if args.output_dir:
                    filename = f'{args.output_dir}/{query_name}.x.{metag_name}.sig.zip'
                    with SaveSignaturesToLocation(filename) as save_ss:
                        save_ss.add(isect_ss)
                else:
                    bulk_save.add(isect_ss)

    if bulk_save:
        bulk_save.close()


if __name__ == '__main__':
    sys.exit(main())
