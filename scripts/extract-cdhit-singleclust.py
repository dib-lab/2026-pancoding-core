#! /usr/bin/env python
import sys
import argparse
from collections import defaultdict
import os

import screed


def main():
    p = argparse.ArgumentParser()
    p.add_argument('cdhit_clstr', help='cd-hit cluster output file')
    p.add_argument('fasta_files', nargs='+', help='files containing sequences')
    p.add_argument('-o', '--outdir', '--output-directory',
                   required=True)
    args = p.parse_args()

    names_to_clusters = {}
    cluster_sizes = defaultdict(int)

    cluster = None
    # for each cluster,
    for line in open(args.cdhit_clstr):
        line = line.strip()

        # get name of cluster in 'cluster',
        if line[0] == '>':      # start of cluster
            assert line[1:9] == 'Cluster ', line[1:9]
            cluster = int(line[9:])
        else:
            # for every sequence, save name of sequence => cluster in
            # 'names_to_clusters'.
            assert cluster is not None
            cluster_sizes[cluster] += 1
            line = line.split('\t')[1]
            seq = line.split(' ')[1].rstrip('...').lstrip('>')
            assert seq not in names_to_clusters
            assert len(seq) == 14, (len(seq), seq)
            names_to_clusters[seq] = cluster

    # now choose all cluster names that have precisely one sequence.
    n_single = 0
    singleclust = set()
    for cluster, size in cluster_sizes.items():
        if size == 1:
            n_single += 1
            singleclust.add(cluster)

    print(n_single)
    print(len(cluster_sizes))

    try:
        os.mkdir(args.outdir)
    except:
        pass

    n_written = 0
    for filename in args.fasta_files:
        basename = os.path.basename(filename)
        outfile = os.path.join(args.outdir, basename)
        fp = open(outfile, 'wt')
        for record in screed.open(filename):
            thisname = record.name[:14]
            thisclust = names_to_clusters[thisname]
            if thisclust in singleclust:
                fp.write(f'>{record.name}\n{record.sequence}\n')
                n_written += 1

    print(n_written, n_single, n_written == n_single)

if __name__ == '__main__':
    sys.exit(main())
