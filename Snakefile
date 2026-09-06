# note: --keep-incomplete

CORE_NAMES=[ x.strip() for x in open('inputs.branchwater/names.list') ]

GTDB_SIG_ZIP='/home/ctbrown/scratch3/2026-gtdb-dl/gtdb-cds-rs226.species.sig.zip'
GTDB_SINGLETON_SIG_ZIP='/home/ctbrown/scratch3/2026-gtdb-dl/gtdb-cds-rs226.species.singleton.sig.zip'
GTDB_LIN_CSV='/group/ctbrowngrp5/sourmash-db.new/gtdb-rs226/gtdb-rs226.lineages.csv'
MAG_NAMES_FASTA='/home/ctbrown/scratch3/sourmash-midgie-raker/outputs.ath/rename/manysketch-renamed.csv'
MAG_LIN_CSV='/home/ctbrown/scratch3/sourmash-midgie-raker/outputs.ath/rename/bin-sketches.lineages.csv'

BW_ROCKSDB='/group/ctbrowngrp5/sra-metagenomes/20241128-k21-s1000'

GRIST_RAND100='/home/ctbrown/scratch3/2025-grist-annie/outputs.rand100/'
GRIST_HIGHCOV='/home/ctbrown/scratch3/2025-grist-annie/outputs.highcov-core/'

RAND_METAG = [ x.strip() for x in open('inputs.mapping/rand_subset.3216.100.txt') ]
print(f'loaded {len(RAND_METAG)} metagenome names for rand mapping')
print(RAND_METAG[:3])

HIGHCOV_METAG = [ x.strip() for x in open('inputs.mapping/highcov-metags.txt') ]
print(f'loaded {len(HIGHCOV_METAG)} metagenome names for highcov mapping')
print(HIGHCOV_METAG[:3])


include: "workflows/cds3.smk"
include: "workflows/branchwater.smk"
include: "workflows/intersections.smk"
include: "workflows/mapping.smk"
