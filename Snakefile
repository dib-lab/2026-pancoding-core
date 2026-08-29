# note: --keep-incomplete

CORE_NAMES=[ x.strip() for x in open('inputs.branchwater/names.list') ]

GTDB_SIG_ZIP='/home/ctbrown/scratch3/2026-gtdb-dl/gtdb-cds-rs226.species.sig.zip'
GTDB_SINGLETON_SIG_ZIP='/home/ctbrown/scratch3/2026-gtdb-dl/gtdb-cds-rs226.species.singleton.sig.zip'
MAG_NAMES_FASTA='/home/ctbrown/scratch3/sourmash-midgie-raker/outputs.ath/rename/manysketch-renamed.csv'
MAG_LIN_CSV='/home/ctbrown/scratch3/sourmash-midgie-raker/outputs.ath/rename/bin-sketches.lineages.csv'

BW_ROCKSDB='/group/ctbrowngrp5/sra-metagenomes/20241128-k21-s1000'

include: "workflows/cds3.smk"

include: "workflows/branchwater.smk"
