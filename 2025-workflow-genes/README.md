# NOTES:

'all-dedup-95.fa' in the top-level directory is the old set of manually
curated genes from 2025-workflow-core99.

'extract-species-genes.py' takes the (old) manually curated species-genes.csv
file and 'all-dedup-95.fa' and extracts the sequences into species-genes.fa.

This is then used in 'new/' to BLAST against the _new_ all-dedup-95.fa
from 2026-pancoding-core.
