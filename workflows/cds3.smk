"""
Build pancoding singlehash sketches and search with them.
"""
import polars as pl

mag_fasta_df = (pl.read_csv(MAG_NAMES_FASTA)
    .with_columns(ident=pl.col('name').str.split(' ').list.get(0))
)
mag_lin_df = pl.read_csv(MAG_LIN_CSV)

# AtH_MAG_N_ etc.
MAG_NAMES = (mag_fasta_df
    ['ident'].to_list()
    )


rule do_cds3:
    input:
        'outputs.cds3/mag-cds.sig.zip',
        'outputs.cds3/mag+gtdb.cds.sig.zip',
        'outputs.cds3/mag+gtdb.cds.singleton.sig.zip',
        'outputs.cds3/gtdb.cds3.x.3216.manysearch.csv',
        'outputs.cds3/mag+gtdb.cds3.x.3216.manysearch.csv',
        'outputs.cds3/mag+gtdb.cds.singleton.describe.csv',

# retrieve the fasta for wildcards.g, where g is a MAG ident
def _get_mag_fasta_path(w):
    name = w.g
    xx_df = mag_fasta_df.filter(pl.col('ident') == name)
    assert len(xx_df) == 1
    filename = xx_df['genome_filename'].item()
    return filename
    
rule run_prodigal:
    input:
        g=_get_mag_fasta_path
    output:
        ffn="outputs.cds3/prodigal/{g}.ffn",
        faa="outputs.cds3/prodigal/{g}.faa",
    threads: 1
    shell: """
        prodigal -i {input.g} -d {output.ffn} -a {output.faa} -o /dev/null
    """

rule make_manysketch_csv:
    input:
        expand('outputs.cds3/prodigal/{n}.ffn', n=MAG_NAMES)
    output:
        touch("outputs.cds3/.prodigal.done")
    shell: """
        scripts/make-prodigal-manysketch.py outputs.cds3/prodigal/ -o \
            outputs.cds3/mag-manysketch.csv
    """

rule run_manysketch:
    input:
        'outputs.cds3/mag-manysketch.csv'
    output:
        'outputs.cds3/mag-cds.sig.zip'
    threads: 32
    shell: """
        sourmash scripts manysketch -p k=21,dna,scaled=1000 {input} -o {output}
    """

rule update_species:
    input:
        gtdb_sigzip=GTDB_SIG_ZIP,
        mag_sigzip='outputs.cds3/mag-cds.sig.zip',
        mag_lincsv=MAG_LIN_CSV,
    output:
        'outputs.cds3/mag+gtdb.cds.sig.zip'
    params:
        ksize=21,
    shell: """
         scripts/merge-mags-gtdb-on-species.py {input.gtdb_sigzip} {input.mag_sigzip} {input.mag_lincsv} -o {output} -k {params.ksize}
    """

rule remove_multihash_d:
    input:
        'outputs.cds3/mag+gtdb.cds.sig.zip'
    output:
        protected(directory('outputs.cds3/mag+gtdb.cds.singleton.sig.d'))
    shell: """
        scripts/remove-multihash-by-sig.py -k 21 {input} -o {output}
    """

rule remove_multihash_cp:
    input:
        'outputs.cds3/mag+gtdb.cds.singleton.sig.d'
    output:
        protected('outputs.cds3/mag+gtdb.cds.singleton.sig.zip',)
    shell: """
        cp outputs.cds3/mag+gtdb.cds.singleton.sig.d/mag+gtdb.cds.sig.zip {output}
    """

rule describe_db:
    input:
        'outputs.cds3/mag+gtdb.cds.singleton.sig.zip'
    output:
        'outputs.cds3/mag+gtdb.cds.singleton.describe.csv'
    shell: """
        sourmash sig describe {input} --csv {output} > /dev/null
    """
        
rule manysearch_mag_gtdb:
    input:
        db='outputs.cds3/mag+gtdb.cds.singleton.sig.zip',
        manifest='3216.manifest.csv',
    output:
        csv=protected('outputs.cds3/mag+gtdb.cds3.x.3216.manysearch.csv')
    threads: 32
    shell: """
        sourmash scripts manysearch -k 21 --scaled=1000 --threshold=0 \
           {input.db} {input.manifest} -o {output.csv} -c {threads}
    """

rule manysearch_gtdb:
    input:
        db=GTDB_SINGLETON_SIG_ZIP,
        manifest='3216.manifest.csv',
    output:
        csv=protected('outputs.cds3/gtdb.cds3.x.3216.manysearch.csv')
    threads: 32
    shell: """
        sourmash scripts manysearch -k 21 --scaled=1000 --threshold=0 \
           {input.db} {input.manifest} -o {output.csv} -c {threads}
    """
