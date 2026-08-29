#! /usr/bin/env python
"""
Remove any hash that occurs more than once across all signatures.
"""
import sys
import argparse
import sourmash
from collections import Counter
from sourmash.sourmash_args import SaveSignaturesToLocation
import os


def main():
    p = argparse.ArgumentParser()
    p.add_argument('pangenome_db', nargs='+')
    p.add_argument('-k', '--ksize', default=31, type=int)
    p.add_argument('-o', '--output-directory', required=True)
    args = p.parse_args()

    try:
        os.mkdir(args.output_directory)
    except:
        pass

    # load and count.
    cnt = Counter()
    for filename in args.pangenome_db:
        db = sourmash.load_file_as_index(filename)
        db = db.select(ksize=args.ksize)
        for n, ss in enumerate(db.signatures()):
            if n % 1000 == 0:
                print('...', n)
            mh = ss.minhash.flatten()
            cnt.update(mh.hashes)

    # calculate all non-singletons
    nonsingletons = set()
    for hashval, count in cnt.most_common():
        if count < 2:
            break
        nonsingletons.add(hashval)

    # print out counts
    print(len(nonsingletons), len(cnt))

    # remove non-singletons
    for filename in args.pangenome_db:
        db = sourmash.load_file_as_index(filename)
        db = db.select(ksize=args.ksize)
        outfile = os.path.join(args.output_directory, os.path.basename(filename))
        with SaveSignaturesToLocation(outfile) as save_sig:
            for n, ss in enumerate(db.signatures()):
                if n % 1000 == 0:
                    print('writing ...', n)

                mh = ss.minhash
                hashes = set(mh.hashes)
                hashes -= nonsingletons
                mh = mh.copy_and_clear()
                mh.add_many(hashes)

                ss = sourmash.SourmashSignature(mh, name=ss.name + ' singlehash')
                save_sig.add(ss)


if __name__ == '__main__':
    sys.exit(main())
