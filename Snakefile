# note: --keep-incomplete

GTDB_SIG_ZIP='/home/ctbrown/scratch3/2026-gtdb-dl/gtdb-cds-rs226.species.sig.zip'
GTDB_SINGLETON_SIG_ZIP='/home/ctbrown/scratch3/2026-gtdb-dl/gtdb-cds-rs226.species.singleton.sig.zip'
MAG_NAMES_FASTA='/home/ctbrown/scratch3/sourmash-midgie-raker/outputs.ath/rename/manysketch-renamed.csv'
MAG_LIN_CSV='/home/ctbrown/scratch3/sourmash-midgie-raker/outputs.ath/rename/bin-sketches.lineages.csv'

include: "workflows/cds3.smk"
