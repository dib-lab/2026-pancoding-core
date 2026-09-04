#! /usr/bin/env python
"""
Keep only hashes that occur more than m times in a collection of other sigs.
"""
import sys
import argparse
import sourmash
from collections import Counter, defaultdict
from sourmash.sourmash_args import SaveSignaturesToLocation
import os


def main():
    p = argparse.ArgumentParser()
    p.add_argument('query')
    p.add_argument('dbs', nargs='+')
    p.add_argument('-k', '--ksize', default=31, type=int)
    p.add_argument('-o', '--output-sig', required=True)
    p.add_argument('-m', '--min-count', type=int, default=2)
    args = p.parse_args()

    query_sig = sourmash.load_file_as_signatures(args.query, ksize=args.ksize)
    query_sig = list(query_sig)
    assert len(query_sig) == 1
    query_sig = query_sig[0]

    query_mh = query_sig.minhash
    query_counts = {}

    for hashval in query_mh.hashes:
        query_counts[hashval] = 0

    # load and count.
    n = 0
    for db in args.dbs:
        db = sourmash.load_file_as_index(db)
        db = db.select(ksize=args.ksize)
        for ss in db.signatures():
            n += 1
            print('...', n)

            mh = ss.minhash.flatten()
            db_hashes = mh.hashes
            for hashval in query_counts:
                if hashval in db_hashes:
                    query_counts[hashval] += 1
        break

    # output query above threshold
    new_query_mh = query_mh.copy_and_clear()
    for hashval, count in query_counts.items():
        if count >= args.min_count:
            new_query_mh.add_hash(hashval)

    print(f'kept {len(new_query_mh)} of {len(query_mh)}')

    ss = sourmash.SourmashSignature(new_query_mh, name=query_sig.name)

    with SaveSignaturesToLocation(args.output_sig) as save_sig:
        save_sig.add(ss)


if __name__ == '__main__':
    sys.exit(main())
